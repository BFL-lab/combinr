//! Single-linkage clustering of alignments by genomic overlap.
//!
//! Replaces PASA's external `slclust` (`SingleLinkageClusterer.pm`). Alignments
//! are bucketed by `(contig, strand)` and, within each bucket, grouped into
//! maximal runs of transitively-overlapping spans. For intervals the overlap
//! relation's connected components are exactly the runs produced by sorting on
//! `lend` and extending a run while the next span starts at or before the run's
//! running maximum `rend` — so a linear sweep computes single linkage without a
//! separate union-find.
//!
//! Cross-strand merging never happens here, and it does not need to: the
//! assembler's [`can_merge`](crate::assemble::compat::can_merge) already rejects
//! opposite orientations, and the orientation wrapper resolves strand per
//! assembly. Strandless (`Unknown`) alignments form their own bucket.

use crate::model::{Alignment, Strand};
use std::collections::BTreeMap;

/// A set of alignments to be assembled together (indices into the input slice).
#[derive(Debug, Clone)]
pub struct Cluster {
    pub contig: String,
    pub strand: Strand,
    pub member_indices: Vec<usize>,
}

/// The strand used for bucketing: the transcribed orientation if known, else the
/// aligned orientation, else `Unknown`.
fn cluster_strand(a: &Alignment) -> Strand {
    if a.spliced_orient.is_known() {
        a.spliced_orient
    } else if a.aligned_orient.is_known() {
        a.aligned_orient
    } else {
        Strand::Unknown
    }
}

/// Cluster `alignments` by `(contig, strand)` single linkage on overlapping
/// spans. Clusters are returned in a deterministic order (by contig, strand,
/// then leftmost coordinate).
pub fn cluster_alignments(alignments: &[Alignment]) -> Vec<Cluster> {
    // Bucket indices by (contig, strand), preserving deterministic key order.
    let mut buckets: BTreeMap<(String, u8), Vec<usize>> = BTreeMap::new();
    for (i, a) in alignments.iter().enumerate() {
        let key = (a.contig.clone(), strand_key(cluster_strand(a)));
        buckets.entry(key).or_default().push(i);
    }

    let mut clusters = Vec::new();
    for ((contig, skey), mut idxs) in buckets {
        let strand = strand_from_key(skey);
        // Sort by (lend, rend) for the sweep; ties stable.
        idxs.sort_by_key(|&i| (alignments[i].coords.lend, alignments[i].coords.rend));

        let mut current: Vec<usize> = Vec::new();
        let mut run_max_rend = i64::MIN;
        for idx in idxs {
            let span = alignments[idx].coords;
            if current.is_empty() || span.lend <= run_max_rend {
                // Inclusive overlap with the running merged interval.
                current.push(idx);
                run_max_rend = run_max_rend.max(span.rend);
            } else {
                clusters.push(Cluster {
                    contig: contig.clone(),
                    strand,
                    member_indices: std::mem::take(&mut current),
                });
                current.push(idx);
                run_max_rend = span.rend;
            }
        }
        if !current.is_empty() {
            clusters.push(Cluster {
                contig: contig.clone(),
                strand,
                member_indices: current,
            });
        }
    }
    clusters
}

fn strand_key(s: Strand) -> u8 {
    match s {
        Strand::Plus => 0,
        Strand::Minus => 1,
        Strand::Unknown => 2,
    }
}
fn strand_from_key(k: u8) -> Strand {
    match k {
        0 => Strand::Plus,
        1 => Strand::Minus,
        _ => Strand::Unknown,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Segment;

    fn al(acc: &str, contig: &str, orient: Strand, segs: &[(i64, i64)]) -> Alignment {
        let mut a = Alignment::new(
            acc,
            segs.iter().map(|&(l, r)| Segment::new(l, r)).collect(),
            orient,
        );
        a.contig = contig.to_string();
        a.spliced_orient = orient;
        a
    }

    fn cluster_accs(c: &Cluster, aligns: &[Alignment]) -> Vec<String> {
        let mut v: Vec<String> = c
            .member_indices
            .iter()
            .map(|&i| aligns[i].acc.clone())
            .collect();
        v.sort();
        v
    }

    #[test]
    fn transitive_overlap_chains_into_one_cluster() {
        // a—b overlap, b—c overlap, a—c do NOT: single linkage joins all three.
        let aligns = vec![
            al("a", "chr1", Strand::Plus, &[(100, 200)]),
            al("b", "chr1", Strand::Plus, &[(180, 300)]),
            al("c", "chr1", Strand::Plus, &[(280, 400)]),
        ];
        let clusters = cluster_alignments(&aligns);
        assert_eq!(clusters.len(), 1);
        assert_eq!(cluster_accs(&clusters[0], &aligns), vec!["a", "b", "c"]);
    }

    #[test]
    fn non_overlapping_genes_are_separate_clusters() {
        let aligns = vec![
            al("a", "chr1", Strand::Plus, &[(100, 200)]),
            al("b", "chr1", Strand::Plus, &[(500, 600)]),
        ];
        let clusters = cluster_alignments(&aligns);
        assert_eq!(clusters.len(), 2);
    }

    #[test]
    fn strand_and_contig_separate_clusters() {
        let aligns = vec![
            al("p", "chr1", Strand::Plus, &[(100, 300)]),
            al("m", "chr1", Strand::Minus, &[(100, 300)]), // antisense, same span
            al("o", "chr2", Strand::Plus, &[(100, 300)]),  // other contig
        ];
        let clusters = cluster_alignments(&aligns);
        assert_eq!(clusters.len(), 3, "split by both strand and contig");
    }

    #[test]
    fn adjacent_non_overlapping_intervals_do_not_cluster() {
        // touching at a single base counts as overlap; a gap of one does not.
        let overlap = vec![
            al("a", "chr1", Strand::Plus, &[(100, 200)]),
            al("b", "chr1", Strand::Plus, &[(200, 300)]), // shares base 200
        ];
        assert_eq!(cluster_alignments(&overlap).len(), 1);

        let gap = vec![
            al("a", "chr1", Strand::Plus, &[(100, 200)]),
            al("b", "chr1", Strand::Plus, &[(201, 300)]), // adjacent, no shared base
        ];
        assert_eq!(cluster_alignments(&gap).len(), 2);
    }
}
