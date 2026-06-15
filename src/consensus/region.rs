//! Per-locus regions for the consensus DP.
//!
//! Unlike the PASA path (which clusters by `(contig, strand)`), the consensus engine
//! runs a single trellis spanning **both strands** of a locus, so evidence is piled
//! into independent regions by **contig only** (see `cluster_spans_by_contig`). Each
//! region is the evidence span padded by a flank. EVM's sliding-window partitioning +
//! recombine DP are intentionally dropped (see `consensus-evm-port`).

use crate::cluster::cluster_spans_by_contig;
use crate::consensus::evidence::EvidenceChain;
use crate::model::Coordset;

/// One independent locus: a contig span (evidence + flank) and the evidence chains
/// (indices into the input slice) that fall in it.
#[derive(Debug, Clone)]
pub struct ConsensusRegion {
    pub contig: String,
    pub span: Coordset,
    pub chain_indices: Vec<usize>,
}

/// Group `chains` into per-contig regions by overlap (single linkage), padding each
/// locus by `flank` bp. The lower bound is clamped at 1; the upper bound is clamped to
/// the contig length in a later milestone (needs the genome).
pub fn build_regions(chains: &[EvidenceChain], flank: i64) -> Vec<ConsensusRegion> {
    let spans: Vec<(String, Coordset)> =
        chains.iter().map(|c| (c.contig.clone(), c.span)).collect();
    cluster_spans_by_contig(&spans)
        .into_iter()
        .map(|cl| {
            let lend = cl
                .member_indices
                .iter()
                .map(|&i| chains[i].span.lend)
                .min()
                .expect("non-empty cluster");
            let rend = cl
                .member_indices
                .iter()
                .map(|&i| chains[i].span.rend)
                .max()
                .expect("non-empty cluster");
            let lo = (lend - flank).max(1);
            let hi = rend + flank;
            ConsensusRegion {
                contig: cl.contig,
                span: Coordset::new(lo, hi),
                chain_indices: cl.member_indices,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consensus::weights::EvClass;
    use crate::model::Strand;

    fn chain(contig: &str, lend: i64, rend: i64) -> EvidenceChain {
        EvidenceChain {
            accession: "x".into(),
            ev_type: "t".into(),
            ev_class: EvClass::AbinitioPrediction,
            weight: 1.0,
            contig: contig.into(),
            orient: Strand::Plus,
            span: Coordset::new(lend, rend),
            links: vec![Coordset::new(lend, rend)],
        }
    }

    #[test]
    fn merges_overlapping_into_region_with_flank() {
        let chains = vec![
            chain("chr1", 1000, 2000),
            chain("chr1", 1900, 2500),
            chain("chr1", 9000, 9500),
        ];
        let regions = build_regions(&chains, 100);
        assert_eq!(regions.len(), 2);
        assert_eq!(regions[0].span, Coordset::new(900, 2600));
        assert_eq!(regions[0].chain_indices.len(), 2);
        assert_eq!(regions[1].chain_indices.len(), 1);
    }

    #[test]
    fn flank_clamped_at_one() {
        let regions = build_regions(&[chain("chr1", 50, 100)], 1000);
        assert_eq!(regions[0].span.lend, 1);
    }

    #[test]
    fn separate_contigs_split() {
        let chains = vec![chain("chr1", 1, 100), chain("chr2", 1, 100)];
        assert_eq!(build_regions(&chains, 10).len(), 2);
    }
}
