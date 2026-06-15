//! Per-base scoring vectors and the intron-score map for one region.
//!
//! [`RegionVectors`] holds the EVM `@CODING_SCORES`, `@INTERGENIC_SCORES` and repeat
//! `@MASK` arrays, but region-local (indexed by offset from the region origin) rather
//! than whole-contig. The coding vector is rebuilt per strand pass; the intergenic
//! vector and [`IntronScores`] are forward-coordinate and shared across both strands.

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

    /// Mark `[lend, rend]` as repeat-masked (excluded from scoring).
    pub fn set_masked(&mut self, lend: i64, rend: i64) {
        for g in lend..=rend {
            if let Some(k) = self.idx(g) {
                self.mask.insert(k);
            }
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
    /// `[lend, rend]` (negative coding scores do not subtract).
    pub fn coding_sum(&self, lend: i64, rend: i64) -> f64 {
        let mut s = 0.0;
        for g in lend..=rend {
            if let Some(k) = self.idx(g) {
                s += self.coding[k].max(0.0);
            }
        }
        s
    }

    /// Reset the coding vector before a new strand pass.
    pub fn clear_coding(&mut self) {
        self.coding.iter_mut().for_each(|v| *v = 0.0);
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
        let mut s = 0.0;
        for g in lend..=rend {
            if let Some(k) = self.idx(g) {
                s += self.intergenic[k];
            }
        }
        s * adjust
    }

    /// Number of unmasked bases in `[lend, rend]` (for length-aware scoring).
    pub fn unmasked_len(&self, lend: i64, rend: i64) -> i64 {
        let mut n = 0;
        for g in lend..=rend {
            if let Some(k) = self.idx(g)
                && !self.mask.contains(k)
            {
                n += 1;
            }
        }
        n
    }
}

/// Evidence-supported introns, keyed by `(intron_lend, intron_rend)` in forward genomic
/// coordinates (the gap between two exons). Ports `%INTRONS_TO_SCORE`,
/// `%INTRONS_TO_EVIDENCE` and `%PREDICTED_INTRONS` (abinitio-only).
#[derive(Default)]
pub struct IntronScores {
    score: HashMap<(i64, i64), f64>,
    evidence: HashMap<(i64, i64), Vec<(String, String)>>,
    predicted: HashMap<(i64, i64), f64>,
}

impl IntronScores {
    /// Accumulate one evidence chain's support for the intron `(lend, rend)`:
    /// `weight * unmasked_len` added to the intron's score (and to the abinitio-only
    /// predicted-intron pool when `is_abinitio`).
    pub fn add(
        &mut self,
        intron: (i64, i64),
        weight: f64,
        unmasked_len: i64,
        acc: &str,
        ev_type: &str,
        is_abinitio: bool,
    ) {
        let contrib = weight * unmasked_len as f64;
        *self.score.entry(intron).or_default() += contrib;
        self.evidence
            .entry(intron)
            .or_default()
            .push((acc.to_string(), ev_type.to_string()));
        if is_abinitio {
            *self.predicted.entry(intron).or_default() += contrib;
        }
    }

    pub fn score(&self, intron: (i64, i64)) -> Option<f64> {
        self.score.get(&intron).copied()
    }

    pub fn predicted(&self, intron: (i64, i64)) -> f64 {
        self.predicted.get(&intron).copied().unwrap_or(0.0)
    }

    pub fn contains(&self, intron: (i64, i64)) -> bool {
        self.score.contains_key(&intron)
    }

    pub fn evidence(&self, intron: (i64, i64)) -> &[(String, String)] {
        self.evidence.get(&intron).map(Vec::as_slice).unwrap_or(&[])
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
    fn coding_sum_ignores_negative_bases() {
        let mut v = RegionVectors::new(1, 10);
        v.add_coverage(1, 5, 2.0);
        // a lone negative paint can't push the *sum* below the positive region
        assert_eq!(v.coding_sum(1, 10), 10.0);
    }

    #[test]
    fn clear_coding_resets() {
        let mut v = RegionVectors::new(1, 10);
        v.add_coverage(1, 10, 4.0);
        v.clear_coding();
        assert_eq!(v.coding_sum(1, 10), 0.0);
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
    fn intron_scores_accumulate_and_track_abinitio() {
        let mut introns = IntronScores::default();
        introns.add((200, 299), 1.0, 100, "a1", "fgenesh", true);
        introns.add((200, 299), 5.0, 100, "p1", "genewise", false);
        assert!(introns.contains((200, 299)));
        assert_eq!(introns.score((200, 299)), Some(600.0)); // 1*100 + 5*100
        assert_eq!(introns.predicted((200, 299)), 100.0); // abinitio only
        assert_eq!(introns.evidence((200, 299)).len(), 2);
        assert_eq!(introns.score((1, 2)), None);
    }
}
