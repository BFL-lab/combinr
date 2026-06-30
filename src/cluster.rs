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

use crate::model::{Alignment, Coordset, Strand};
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
///
/// `min_overlap_frac` is PASA's `--stringent_alignment_overlap`: two spans only
/// link when their overlap is at least that percent of the **shorter** span. At
/// `0.0` (the default) any single shared base links them — identical to the
/// historical any-overlap sweep.
pub fn cluster_alignments(alignments: &[Alignment], min_overlap_frac: f64) -> Vec<Cluster> {
    // Bucket indices by (contig, strand). `Strand` derives `Ord` (Plus < Minus <
    // Unknown), so the `BTreeMap` key order is deterministic without a manual table.
    let mut buckets: BTreeMap<(String, Strand), Vec<usize>> = BTreeMap::new();
    for (i, a) in alignments.iter().enumerate() {
        let key = (a.contig.clone(), cluster_strand(a));
        buckets.entry(key).or_default().push(i);
    }

    let mut clusters = Vec::new();
    for ((contig, strand), idxs) in buckets {
        let spans: Vec<Coordset> = idxs.iter().map(|&i| alignments[i].coords).collect();
        for group in single_linkage_groups(&spans, min_overlap_frac) {
            clusters.push(Cluster {
                contig: contig.clone(),
                strand,
                member_indices: group.into_iter().map(|l| idxs[l]).collect(),
            });
        }
    }
    clusters
}

/// A set of spans grouped into one region (indices into the input slice), keyed by
/// contig only. Used by the consensus path, whose trellis spans both strands of a locus.
#[derive(Debug, Clone)]
pub struct ContigCluster {
    pub contig: String,
    pub member_indices: Vec<usize>,
}

/// Single-linkage clustering of `(contig, span)` items by inclusive genomic overlap,
/// bucketed by **contig only** (strand ignored) — the consensus DP integrates both
/// strands of a locus in one pass. Always any-overlap (no stringency knob): the EVM
/// region partitioner intentionally piles all overlapping evidence into one region and
/// lets the trellis separate genes.
pub fn cluster_spans_by_contig(spans: &[(String, Coordset)]) -> Vec<ContigCluster> {
    let mut buckets: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    for (i, (contig, _)) in spans.iter().enumerate() {
        buckets.entry(contig.as_str()).or_default().push(i);
    }

    let mut clusters = Vec::new();
    for (contig, idxs) in buckets {
        let bucket_spans: Vec<Coordset> = idxs.iter().map(|&i| spans[i].1).collect();
        for group in single_linkage_groups(&bucket_spans, 0.0) {
            clusters.push(ContigCluster {
                contig: contig.to_string(),
                member_indices: group.into_iter().map(|l| idxs[l]).collect(),
            });
        }
    }
    clusters
}

// ---------------------------------------------------------------------------
// Shared single-linkage core (overlap-fraction aware)
// ---------------------------------------------------------------------------

/// Edge test for single linkage: two spans link iff they overlap by at least
/// `min_frac` percent of the **shorter** span's genomic length. `min_frac <= 0.0`
/// links on any single shared base (the historical any-overlap behavior).
fn span_overlap_edge(a: Coordset, b: Coordset, min_frac: f64) -> bool {
    let overlap_bp = a.overlap_len(&b);
    if overlap_bp <= 0 {
        return false; // no shared base → never linked
    }
    if min_frac <= 0.0 {
        return true; // any-overlap (default): exactly Coordset::overlaps_inclusive
    }
    let shorter = a.len().min(b.len());
    (overlap_bp as f64) / (shorter as f64) * 100.0 >= min_frac
}

/// Disjoint-set forest (path compression + union by rank) for single linkage.
struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }
    fn find(&mut self, x: usize) -> usize {
        let mut root = x;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        let mut cur = x;
        while self.parent[cur] != root {
            let next = self.parent[cur];
            self.parent[cur] = root;
            cur = next;
        }
        root
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return;
        }
        match self.rank[ra].cmp(&self.rank[rb]) {
            std::cmp::Ordering::Less => self.parent[ra] = rb,
            std::cmp::Ordering::Greater => self.parent[rb] = ra,
            std::cmp::Ordering::Equal => {
                self.parent[rb] = ra;
                self.rank[ra] += 1;
            }
        }
    }
}

/// Single-linkage grouping of `spans` by genomic-span overlap, requiring overlap
/// `>= min_frac` percent of the shorter span (PASA `--stringent_alignment_overlap`).
/// At `min_frac <= 0.0` this reduces EXACTLY to single linkage on any 1-bp overlap, so
/// the connected components — and their ordering — match the historical lend-sorted
/// sweep byte-for-byte.
///
/// Returns groups of indices into `spans`. Members within a group are ordered by
/// `(lend, rend, index)`; groups are ordered by their leftmost member. A fraction
/// threshold needs true single linkage (union-find), not a running-max sweep, because
/// two spans both overlapping a third need not meet the fraction with each other. The
/// `(lend, rend)`-sorted active window keeps it ~`O(n · pile-up depth)`.
pub(crate) fn single_linkage_groups(spans: &[Coordset], min_frac: f64) -> Vec<Vec<usize>> {
    let n = spans.len();
    if n == 0 {
        return Vec::new();
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&i| (spans[i].lend, spans[i].rend, i));

    let mut uf = UnionFind::new(n);
    let mut active: Vec<usize> = Vec::new(); // spans whose rend can still reach forward
    for &i in &order {
        let lend_i = spans[i].lend;
        active.retain(|&j| spans[j].rend >= lend_i);
        for &j in &active {
            if span_overlap_edge(spans[i], spans[j], min_frac) {
                uf.union(i, j);
            }
        }
        active.push(i);
    }

    let mut by_root: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..n {
        let r = uf.find(i);
        by_root.entry(r).or_default().push(i);
    }
    let mut groups: Vec<Vec<usize>> = by_root.into_values().collect();
    for g in &mut groups {
        g.sort_by_key(|&i| (spans[i].lend, spans[i].rend, i));
    }
    groups.sort_by_key(|g| (spans[g[0]].lend, spans[g[0]].rend, g[0]));
    groups
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
        let clusters = cluster_alignments(&aligns, 0.0);
        assert_eq!(clusters.len(), 1);
        assert_eq!(cluster_accs(&clusters[0], &aligns), vec!["a", "b", "c"]);
    }

    #[test]
    fn non_overlapping_genes_are_separate_clusters() {
        let aligns = vec![
            al("a", "chr1", Strand::Plus, &[(100, 200)]),
            al("b", "chr1", Strand::Plus, &[(500, 600)]),
        ];
        let clusters = cluster_alignments(&aligns, 0.0);
        assert_eq!(clusters.len(), 2);
    }

    #[test]
    fn strand_and_contig_separate_clusters() {
        let aligns = vec![
            al("p", "chr1", Strand::Plus, &[(100, 300)]),
            al("m", "chr1", Strand::Minus, &[(100, 300)]), // antisense, same span
            al("o", "chr2", Strand::Plus, &[(100, 300)]),  // other contig
        ];
        let clusters = cluster_alignments(&aligns, 0.0);
        assert_eq!(clusters.len(), 3, "split by both strand and contig");
    }

    #[test]
    fn adjacent_non_overlapping_intervals_do_not_cluster() {
        // touching at a single base counts as overlap; a gap of one does not.
        let overlap = vec![
            al("a", "chr1", Strand::Plus, &[(100, 200)]),
            al("b", "chr1", Strand::Plus, &[(200, 300)]), // shares base 200
        ];
        assert_eq!(cluster_alignments(&overlap, 0.0).len(), 1);

        let gap = vec![
            al("a", "chr1", Strand::Plus, &[(100, 200)]),
            al("b", "chr1", Strand::Plus, &[(201, 300)]), // adjacent, no shared base
        ];
        assert_eq!(cluster_alignments(&gap, 0.0).len(), 2);
    }

    #[test]
    fn stringent_overlap_splits_collinear_tip_overlap() {
        // Two long single-exon spans that overlap only at their tips:
        // a=(100,1100) len 1001, b=(1000,2000) len 1001, overlap = 1100-1000+1 = 101
        // → 101/1001 ≈ 10% of the shorter. Below 30% they must NOT cluster; at the
        // default 0.0 (any overlap) they do.
        let aligns = vec![
            al("a", "chr1", Strand::Plus, &[(100, 1100)]),
            al("b", "chr1", Strand::Plus, &[(1000, 2000)]),
        ];
        assert_eq!(cluster_alignments(&aligns, 0.0).len(), 1);
        assert_eq!(cluster_alignments(&aligns, 30.0).len(), 2);
    }

    #[test]
    fn stringent_overlap_keeps_contained_short_transcript() {
        // A short span fully inside a long one: overlap = 100% of the shorter span,
        // so it stays clustered even at a high threshold (legitimate evidence).
        let aligns = vec![
            al("long", "chr1", Strand::Plus, &[(100, 2000)]),
            al("short", "chr1", Strand::Plus, &[(500, 600)]),
        ];
        assert_eq!(cluster_alignments(&aligns, 90.0).len(), 1);
    }

    #[test]
    fn stringent_overlap_threshold_is_inclusive() {
        // shorter span len 100; overlap exactly 30 bp = 30% → clusters at min_frac=30
        // (predicate uses >=). a=(1,100) len 100, b=(71,500): overlap = 100-71+1 = 30.
        let aligns = vec![
            al("a", "chr1", Strand::Plus, &[(1, 100)]),
            al("b", "chr1", Strand::Plus, &[(71, 500)]),
        ];
        assert_eq!(cluster_alignments(&aligns, 30.0).len(), 1);
        // a hair stricter and the 30% edge no longer qualifies.
        assert_eq!(cluster_alignments(&aligns, 30.5).len(), 2);
    }

    #[test]
    fn stringent_overlap_preserves_transitive_single_linkage() {
        // a–b and b–c each overlap >=50% of their shorter span, but a–c not at all:
        // single linkage (union-find) still joins all three into one cluster.
        let aligns = vec![
            al("a", "chr1", Strand::Plus, &[(100, 200)]),
            al("b", "chr1", Strand::Plus, &[(150, 250)]), // shares 51 bp with a
            al("c", "chr1", Strand::Plus, &[(200, 300)]), // shares 51 bp with b, 1 bp with a
        ];
        assert_eq!(cluster_alignments(&aligns, 50.0).len(), 1);
    }

    #[test]
    fn span_clustering_by_contig_ignores_strand() {
        let spans = vec![
            ("chr1".to_string(), Coordset::new(100, 200)),
            ("chr1".to_string(), Coordset::new(180, 300)), // overlaps prev
            ("chr1".to_string(), Coordset::new(500, 600)), // separate
            ("chr2".to_string(), Coordset::new(1, 50)),
        ];
        let cs = cluster_spans_by_contig(&spans);
        // chr1: {0,1} merged + {2}; chr2: {3}
        assert_eq!(cs.len(), 3);
        assert_eq!(cs[0].contig, "chr1");
        assert_eq!(cs[0].member_indices.len(), 2);
        assert_eq!(cs[1].member_indices, vec![2]);
        assert_eq!(cs[2].contig, "chr2");
    }
}
