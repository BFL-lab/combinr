//! Per-base scoring vectors and the intron-score map for one region.
//!
//! [`RegionVectors`] holds the EVM `@CODING_SCORES`, `@INTERGENIC_SCORES` and repeat
//! `@MASK` arrays, but region-local (indexed by offset from the region origin) rather
//! than whole-contig. All three are built fresh for each strand pass — the engine
//! reverse-complements the region and rebuilds the candidates for the minus strand.

use fixedbitset::FixedBitSet;
use std::collections::HashMap;

/// Region-local per-base vectors. Index `k` corresponds to genomic position
/// `origin + k`, for `k in 0..len`.
pub struct RegionVectors {
    origin: i64,
    len: usize,
    coding: Vec<f64>,
    intergenic: Vec<f64>,
    mask: FixedBitSet,
}

impl RegionVectors {
    pub fn new(origin: i64, len: usize) -> Self {
        RegionVectors {
            origin,
            len,
            coding: vec![0.0; len],
            intergenic: vec![0.0; len],
            mask: FixedBitSet::with_capacity(len),
        }
    }

    fn idx(&self, g: i64) -> Option<usize> {
        if g < self.origin {
            return None;
        }
        let k = (g - self.origin) as usize;
        (k < self.len).then_some(k)
    }

    /// The half-open offset range covering the genomic span `[lend, rend]`, clipped to the
    /// region. Out-of-range bases simply fall outside the range (they contribute nothing),
    /// matching the per-base `idx` bounds check the mask-free scoring loops used to do.
    fn clamp(&self, lend: i64, rend: i64) -> std::ops::Range<usize> {
        let lo = (lend - self.origin).clamp(0, self.len as i64) as usize;
        let hi = (rend - self.origin + 1).clamp(0, self.len as i64) as usize;
        lo..hi.max(lo)
    }

    /// Mark `[lend, rend]` as repeat-masked (excluded from scoring).
    /// Bases outside the region are ignored.
    pub fn set_masked(&mut self, lend: i64, rend: i64) {
        let r = self.clamp(lend, rend);
        if !r.is_empty() {
            self.mask.insert_range(r);
        }
    }

    pub fn is_masked(&self, g: i64) -> bool {
        self.idx(g).is_some_and(|k| self.mask.contains(k))
    }

    /// `add_match_coverage` (~line 2620): add `weight` to each unmasked base in
    /// `[lend, rend]`. A negative weight clamps the running coding score at 0 so it
    /// cannot defeat the coding potential.
    pub fn add_coverage(&mut self, lend: i64, rend: i64, weight: f64) {
        for g in lend..=rend {
            if let Some(k) = self.idx(g) {
                if self.mask.contains(k) {
                    continue;
                }
                if weight < 0.0 {
                    self.coding[k] = (self.coding[k] + weight).max(0.0);
                } else {
                    self.coding[k] += weight;
                }
            }
        }
    }

    /// `score_exons` coding contribution (~line 2230): `Σ max(0, coding[i])` over
    /// `[lend, rend]` (negative coding scores do not subtract). The mask is not consulted
    /// here — masked bases were never painted, so they are already 0.
    pub fn coding_sum(&self, lend: i64, rend: i64) -> f64 {
        self.coding[self.clamp(lend, rend)]
            .iter()
            .map(|&v| v.max(0.0))
            .sum()
    }

    /// Add `weight` to the intergenic score of each unmasked base in `[lend, rend]`.
    pub fn add_intergenic(&mut self, lend: i64, rend: i64, weight: f64) {
        for g in lend..=rend {
            if let Some(k) = self.idx(g)
                && !self.mask.contains(k)
            {
                self.intergenic[k] += weight;
            }
        }
    }

    /// Raise the intergenic score of each unmasked base in `[lend, rend]` to at least
    /// `level` (peak augmentation — `augment_intergenic_from_start_stop_peaks`).
    pub fn set_intergenic_max(&mut self, lend: i64, rend: i64, level: f64) {
        for g in lend..=rend {
            if let Some(k) = self.idx(g)
                && !self.mask.contains(k)
            {
                self.intergenic[k] = self.intergenic[k].max(level);
            }
        }
    }

    /// `calc_intergenic_score` (~line 3294): `Σ intergenic[i]` over `[lend, rend]`,
    /// scaled by `adjust` (EVM's `INTERGENIC_SCORE_ADJUST_FACTOR`).
    pub fn intergenic_score(&self, lend: i64, rend: i64, adjust: f64) -> f64 {
        self.intergenic[self.clamp(lend, rend)].iter().sum::<f64>() * adjust
    }

    /// Number of unmasked bases in `[lend, rend]` (for length-aware scoring).
    /// Bases outside the region are not counted.
    pub fn unmasked_len(&self, lend: i64, rend: i64) -> i64 {
        let r = self.clamp(lend, rend);
        (r.len() - self.mask.count_ones(r)) as i64
    }
}

/// Evidence-supported introns, keyed by `(intron_lend, intron_rend)` in the strand pass's
/// working coordinates (forward, or RC-local for the minus pass): the gap between two
/// exons. Ports EVM's `%INTRONS_TO_SCORE`, plus the per-intron `(accession, ev_type)`
/// attribution the evidence report lists.
#[derive(Default)]
pub struct IntronScores {
    score: HashMap<(i64, i64), f64>,
    evidence: HashMap<(i64, i64), Vec<(String, String)>>,
}

impl IntronScores {
    /// Accumulate one evidence chain's support for the intron `(lend, rend)`:
    /// `weight * unmasked_len` added to the intron's score.
    pub fn add(&mut self, intron: (i64, i64), weight: f64, unmasked_len: i64) {
        *self.score.entry(intron).or_default() += weight * unmasked_len as f64;
    }

    /// Record that the chain `(accession, ev_type)` supports `intron` (deduped, in
    /// insertion order).
    pub fn attribute(&mut self, intron: (i64, i64), accession: &str, ev_type: &str) {
        let ev = self.evidence.entry(intron).or_default();
        if !ev.iter().any(|(a, t)| a == accession && t == ev_type) {
            ev.push((accession.to_string(), ev_type.to_string()));
        }
    }

    /// The `(accession, ev_type)` pairs supporting `intron` (empty if none recorded).
    pub fn evidence(&self, intron: (i64, i64)) -> &[(String, String)] {
        self.evidence.get(&intron).map_or(&[], Vec::as_slice)
    }

    pub fn score(&self, intron: (i64, i64)) -> Option<f64> {
        self.score.get(&intron).copied()
    }

    pub fn contains(&self, intron: (i64, i64)) -> bool {
        self.score.contains_key(&intron)
    }

    pub fn len(&self) -> usize {
        self.score.len()
    }

    pub fn is_empty(&self) -> bool {
        self.score.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coverage_accumulates_and_clamps_negative() {
        let mut v = RegionVectors::new(100, 50); // covers genomic 100..149
        v.add_coverage(100, 109, 5.0);
        v.add_coverage(105, 114, 3.0);
        // base 107 has 5+3=8; base 102 has 5; base 112 has 3
        assert_eq!(v.coding_sum(107, 107), 8.0);
        assert_eq!(v.coding_sum(102, 102), 5.0);
        assert_eq!(v.coding_sum(112, 112), 3.0);
        // negative weight clamps at 0, never below
        v.add_coverage(102, 102, -100.0);
        assert_eq!(v.coding_sum(102, 102), 0.0);
    }

    #[test]
    fn mask_excludes_bases_from_coverage_and_intergenic() {
        let mut v = RegionVectors::new(1, 20);
        v.set_masked(5, 10);
        v.add_coverage(1, 20, 1.0);
        // bases 5..10 masked -> not painted; 14 unmasked bases contribute
        assert_eq!(v.coding_sum(1, 20), 14.0);
        assert!(v.is_masked(7) && !v.is_masked(4));
        assert_eq!(v.unmasked_len(1, 20), 14);
    }

    #[test]
    fn set_masked_clips_to_region() {
        let mut v = RegionVectors::new(100, 50); // covers genomic 100..149
        v.set_masked(90, 104); // overhangs the left end
        v.set_masked(145, 200); // overhangs the right end
        v.set_masked(1, 99); // entirely before: no-op
        v.set_masked(150, 300); // entirely after: no-op
        let masked: Vec<i64> = (90..160).filter(|&g| v.is_masked(g)).collect();
        assert_eq!(
            masked,
            vec![100, 101, 102, 103, 104, 145, 146, 147, 148, 149]
        );
        let mut all = RegionVectors::new(100, 50);
        all.set_masked(0, 1000); // spans both ends
        assert_eq!(all.unmasked_len(0, 1000), 0);
    }

    #[test]
    fn unmasked_len_matches_per_base_count() {
        let mut v = RegionVectors::new(100, 50);
        v.set_masked(103, 107);
        v.set_masked(120, 120);
        v.set_masked(140, 160);
        for (lend, rend) in [
            (90, 200),
            (100, 149),
            (104, 121),
            (50, 99),
            (150, 170),
            (120, 120),
        ] {
            let expect = (lend..=rend)
                .filter(|&g| (100..150).contains(&g) && !v.is_masked(g))
                .count() as i64;
            assert_eq!(v.unmasked_len(lend, rend), expect, "[{lend},{rend}]");
        }
    }

    #[test]
    fn coding_sum_ignores_negative_bases() {
        let mut v = RegionVectors::new(1, 10);
        v.add_coverage(1, 5, 2.0);
        // a lone negative paint can't push the *sum* below the positive region
        assert_eq!(v.coding_sum(1, 10), 10.0);
    }

    #[test]
    fn intergenic_score_scales_by_adjust() {
        let mut v = RegionVectors::new(1, 10);
        v.add_intergenic(1, 5, 2.0);
        assert_eq!(v.intergenic_score(1, 10, 1.0), 10.0);
        assert_eq!(v.intergenic_score(1, 10, 0.5), 5.0);
    }

    #[test]
    fn set_intergenic_max_raises_floor_only() {
        let mut v = RegionVectors::new(1, 100);
        v.add_intergenic(1, 50, 2.0);
        v.set_intergenic_max(1, 100, 5.0);
        assert_eq!(v.intergenic_score(1, 1, 1.0), 5.0); // 2 -> 5
        assert_eq!(v.intergenic_score(60, 60, 1.0), 5.0); // 0 -> 5
        v.set_intergenic_max(1, 1, 3.0); // lower level does not reduce
        assert_eq!(v.intergenic_score(1, 1, 1.0), 5.0);
    }

    #[test]
    fn intron_scores_accumulate() {
        let mut introns = IntronScores::default();
        introns.add((200, 299), 1.0, 100);
        introns.add((200, 299), 5.0, 100);
        assert!(introns.contains((200, 299)));
        assert_eq!(introns.score((200, 299)), Some(600.0)); // 1*100 + 5*100
        assert_eq!(introns.score((1, 2)), None);
    }
}
