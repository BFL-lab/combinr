//! Low-support gene assessment — **flag, don't drop** by default.
//!
//! EVM's `filter_predictions_low_support` (~3436) eliminates a gene when its
//! coding/noncoding score ratio < `MIN_CODING_NONCODING_SCORE_RATIO` (0.75) or its
//! coding length is too short, which can leave a locus blank. Per
//! `prefer-flag-not-drop-genes`, we compute the same metrics but DEFAULT to keeping the
//! gene and tagging it `low_support`; `--strict` restores EVM's dropping. The recursion
//! that consumes this only treats *retained* genes' spans as covered, fixing EVM's bug
//! where an eliminated gene still blocks re-search.

use crate::consensus::exon::ExonCandidate;
use crate::consensus::trellis::ConsensusGene;
use crate::consensus::vectors::RegionVectors;

/// Filtering knobs.
#[derive(Clone, Copy, Debug)]
pub struct FilterParams {
    pub min_score_ratio: f64,
    pub min_coding_length: i64,
    pub intergenic_adjust: f64,
    /// When true, low-support genes are dropped (EVM behaviour) rather than flagged.
    pub strict: bool,
}

/// Support metrics for one gene.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SupportFlags {
    pub score_ratio: f64,
    pub coding_length: i64,
    pub low_support: bool,
}

/// Assess a gene's support. The noncoding-equivalent is the intergenic score over the
/// gene span (a simplified `filter_predictions_low_support`; the predicted-intron offset
/// refinement rides along with the fuller intergenic model). Floors at `0.0001 × score`
/// to avoid divide-by-zero, exactly as EVM does.
pub fn assess(
    gene: &ConsensusGene,
    exons: &[ExonCandidate],
    vectors: &RegionVectors,
    p: &FilterParams,
) -> SupportFlags {
    let lend = gene
        .exon_indices
        .iter()
        .map(|&i| exons[i].coords.lend)
        .min()
        .expect("gene has exons");
    let rend = gene
        .exon_indices
        .iter()
        .map(|&i| exons[i].coords.rend)
        .max()
        .expect("gene has exons");
    let coding_length: i64 = gene
        .exon_indices
        .iter()
        .map(|&i| exons[i].coords.len())
        .sum();

    let raw_noncoding = vectors.intergenic_score(lend, rend, p.intergenic_adjust);
    let noncoding_equiv = if raw_noncoding <= 0.0 {
        0.0001 * gene.score.max(1.0)
    } else {
        raw_noncoding
    };
    let score_ratio = gene.score / noncoding_equiv;

    let low_support = score_ratio < p.min_score_ratio || coding_length < p.min_coding_length;
    SupportFlags {
        score_ratio,
        coding_length,
        low_support,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consensus::exon::{ExonType, end_frame_for};
    use crate::model::{Coordset, Strand};

    fn exon(lend: i64, rend: i64) -> ExonCandidate {
        ExonCandidate {
            coords: Coordset::new(lend, rend),
            orient: Strand::Plus,
            exon_type: ExonType::Single,
            start_frame: 1,
            end_frame: end_frame_for(1, rend - lend + 1),
            left_seq_boundary: [b'A', b'A'],
            right_seq_boundary: [b'A', b'A'],
            evidence: vec![],
        }
    }

    fn gene(score: f64, lend: i64, rend: i64) -> (ConsensusGene, Vec<ExonCandidate>) {
        (
            ConsensusGene {
                exon_indices: vec![0],
                orient: Strand::Plus,
                score,
            },
            vec![exon(lend, rend)],
        )
    }

    fn params() -> FilterParams {
        FilterParams {
            min_score_ratio: 0.75,
            min_coding_length: 150,
            intergenic_adjust: 1.0,
            strict: false,
        }
    }

    #[test]
    fn short_gene_is_flagged_low_support() {
        let (g, exons) = gene(1000.0, 100, 200); // 101 bp < 150
        let v = RegionVectors::new(1, 1000);
        let s = assess(&g, &exons, &v, &params());
        assert!(s.low_support);
        assert_eq!(s.coding_length, 101);
    }

    #[test]
    fn long_high_scoring_gene_not_flagged() {
        let (g, exons) = gene(1000.0, 100, 400); // 301 bp, no intergenic noncoding
        let v = RegionVectors::new(1, 1000); // intergenic all zero -> tiny noncoding -> huge ratio
        let s = assess(&g, &exons, &v, &params());
        assert!(!s.low_support);
    }

    #[test]
    fn low_ratio_is_flagged() {
        let (g, exons) = gene(10.0, 100, 400); // 301 bp coding, fine length
        let mut v = RegionVectors::new(1, 1000);
        v.add_intergenic(100, 400, 1.0); // noncoding 301 >> score 10 -> ratio ~0.03
        let s = assess(&g, &exons, &v, &params());
        assert!(s.score_ratio < 0.75);
        assert!(s.low_support);
    }
}
