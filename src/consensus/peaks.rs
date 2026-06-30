//! Optional start/stop peak augmentation (EVM `augment_intergenic_from_start_stop_peaks`).
//!
//! Opt-in via `--peak-augment` (off by default). When enabled, protein chain termini feed
//! a per-base start/stop signal; clusters of that signal that clear the gene-prediction
//! weight threshold become peaks, and the intergenic score between each peak and the
//! nearest matching exon edge is maxed out, promoting a gene boundary there. When the
//! feature is off, none of this is allocated or computed (the builder holds `None`).

use crate::consensus::exon::{ExonCandidate, ExonType};
use crate::consensus::vectors::RegionVectors;

/// Window (bp) for clustering start/stop termini signals into peaks.
const CHAIN_TERMINI_WINDOW: i64 = 250;
/// Max bp between a start/stop peak and the exon it augments the intergenic score toward.
const START_STOP_RANGE: i64 = 500;

/// Per-base start/stop termini signal (from protein chain ends), accumulated only when
/// peak augmentation is enabled, then resolved into intergenic-score peaks by [`augment`].
///
/// [`augment`]: PeakSignal::augment
pub struct PeakSignal {
    origin: i64,
    begins: Vec<f64>,
    ends: Vec<f64>,
}

impl PeakSignal {
    pub fn new(origin: i64, len: usize) -> Self {
        PeakSignal {
            origin,
            begins: vec![0.0; len],
            ends: vec![0.0; len],
        }
    }

    fn idx(&self, g: i64) -> Option<usize> {
        let k = g - self.origin;
        (k >= 0 && (k as usize) < self.begins.len()).then_some(k as usize)
    }

    /// Record a protein chain's 5' terminus signal at genomic `g`.
    pub fn bump_begin(&mut self, g: i64, w: f64) {
        if let Some(k) = self.idx(g) {
            self.begins[k] += w;
        }
    }

    /// Record a protein chain's 3' terminus signal at genomic `g`.
    pub fn bump_end(&mut self, g: i64, w: f64) {
        if let Some(k) = self.idx(g) {
            self.ends[k] += w;
        }
    }

    /// `augment_intergenic_from_start_stop_peaks` (~4262): find start/stop termini peaks
    /// and max out the intergenic score between each peak and the nearest matching exon,
    /// promoting a gene boundary there. `level` is `SUM_GENEPRED_WEIGHTS` (the peak
    /// threshold and the augmented intergenic level).
    pub fn augment(&self, exons: &[ExonCandidate], vectors: &mut RegionVectors, level: f64) {
        let starts = find_peaks(&self.begins, self.origin, CHAIN_TERMINI_WINDOW, level);
        let stops = find_peaks(&self.ends, self.origin, CHAIN_TERMINI_WINDOW, level);
        for p in starts {
            if let Some(edge) = nearest_exon_edge(exons, p, true) {
                vectors.set_intergenic_max(edge.min(p), edge.max(p), level);
            }
        }
        for p in stops {
            if let Some(edge) = nearest_exon_edge(exons, p, false) {
                vectors.set_intergenic_max(edge.min(p), edge.max(p), level);
            }
        }
    }
}

/// Nearest initial/single exon `lend` to a start peak (`start=true`) or terminal/single
/// exon `rend` to a stop peak, within `START_STOP_RANGE`.
fn nearest_exon_edge(exons: &[ExonCandidate], peak: i64, start: bool) -> Option<i64> {
    let mut best: Option<(i64, i64)> = None; // (distance, edge)
    for e in exons {
        let edge = match (start, e.exon_type) {
            (true, ExonType::Initial | ExonType::Single) => e.coords.lend,
            (false, ExonType::Terminal | ExonType::Single) => e.coords.rend,
            _ => continue,
        };
        let dist = (edge - peak).abs();
        if dist <= START_STOP_RANGE && best.map(|(d, _)| dist < d).unwrap_or(true) {
            best = Some((dist, edge));
        }
    }
    best.map(|(_, edge)| edge)
}

/// Cluster non-zero per-base signal within `window` bp and emit the max-signal genomic
/// position of each cluster whose total exceeds `threshold` (a simplified `analyze_peaks`).
fn find_peaks(v: &[f64], origin: i64, window: i64, threshold: f64) -> Vec<i64> {
    let sig: Vec<(i64, f64)> = v
        .iter()
        .enumerate()
        .filter(|&(_, x)| *x > 0.0)
        .map(|(k, &x)| (origin + k as i64, x))
        .collect();
    let mut peaks = Vec::new();
    let mut i = 0;
    while i < sig.len() {
        let start = sig[i].0;
        let (mut sum, mut best_pos, mut best_val) = (0.0, sig[i].0, sig[i].1);
        let mut j = i;
        while j < sig.len() && sig[j].0 - start <= window {
            sum += sig[j].1;
            if sig[j].1 > best_val {
                best_val = sig[j].1;
                best_pos = sig[j].0;
            }
            j += 1;
        }
        if sum > threshold {
            peaks.push(best_pos);
        }
        i = j;
    }
    peaks
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn find_peaks_clusters_and_thresholds() {
        let mut v = vec![0.0; 600];
        v[0] = 2.0; // genomic 100
        v[5] = 3.0; // genomic 105 (same cluster within window 250); sum 5, max at 105
        v[500] = 10.0; // genomic 600 (separate cluster); sum 10
        assert_eq!(find_peaks(&v, 100, 250, 1.0), vec![105, 600]);
        assert_eq!(find_peaks(&v, 100, 250, 6.0), vec![600]); // only the sum-10 cluster passes
    }
}
