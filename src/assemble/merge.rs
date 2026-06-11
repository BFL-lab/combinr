//! Folding compatible alignments into a consensus structure.
//!
//! Port of `pasa_cpp/cdna_alignment_assembler.cpp`: `mergeAlignments`
//! (lines 432–551). `create_assembly` lives on the assembler (it needs the
//! alignment table) but folds left-to-right using this function.

use crate::model::{Alignment, Coordset, Segment};
use std::collections::HashSet;

/// Port of `mergeAlignments`: merge two compatible alignments `a`/`b` into one
/// consensus alignment. Interior boundaries that are splice junctions in either
/// input are preserved exactly; non-junction termini extend to the min/max; and
/// `b` segments overlapping nothing in `a` are appended.
pub fn merge_alignments(a: &Alignment, b: &Alignment) -> Alignment {
    let orient = a.aligned_orient; // a and b share orientation

    // Collect splice-junction-anchored coordinates from BOTH alignments.
    let mut left_splice: HashSet<i64> = HashSet::new();
    let mut right_splice: HashSet<i64> = HashSet::new();
    for seg in a.segments.iter().chain(b.segments.iter()) {
        if seg.left_splice {
            left_splice.insert(seg.coords.lend);
        }
        if seg.right_splice {
            right_splice.insert(seg.coords.rend);
        }
    }

    let mut merged: Vec<Coordset> = Vec::new();

    // For each a-segment, find the first overlapping b-segment and merge bounds.
    for a_seg in &a.segments {
        let a1l = a_seg.coords.lend;
        let a1r = a_seg.coords.rend;
        let mut merged_coord: Option<Coordset> = None;

        for b_seg in &b.segments {
            let a2l = b_seg.coords.lend;
            let a2r = b_seg.coords.rend;
            if a_seg.coords.overlaps_inclusive(&b_seg.coords) {
                let merged_lend = if left_splice.contains(&a1l) {
                    a1l
                } else if left_splice.contains(&a2l) {
                    a2l
                } else {
                    a1l.min(a2l)
                };
                let merged_rend = if right_splice.contains(&a1r) {
                    a1r
                } else if right_splice.contains(&a2r) {
                    a2r
                } else {
                    a1r.max(a2r)
                };
                merged_coord = Some(Coordset {
                    lend: merged_lend,
                    rend: merged_rend,
                });
                break;
            }
        }

        // No overlap found → keep the a-segment coords unchanged.
        merged.push(merged_coord.unwrap_or(Coordset {
            lend: a1l,
            rend: a1r,
        }));
    }

    // Append b-segments that overlap none of the merged coords.
    for b_seg in &b.segments {
        let overlaps_any = merged.iter().any(|m| b_seg.coords.overlaps_inclusive(m));
        if !overlaps_any {
            merged.push(b_seg.coords);
        }
    }

    let segs: Vec<Segment> = merged
        .into_iter()
        .map(|c| Segment::new(c.lend, c.rend))
        .collect();

    // Alignment::new runs refine() → sort by lend + reclassify junctions.
    Alignment::new(String::new(), segs, orient)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Alignment, Segment, Strand};

    fn al(segs: &[(i64, i64)]) -> Alignment {
        Alignment::new(
            "x",
            segs.iter().map(|&(l, r)| Segment::new(l, r)).collect(),
            Strand::Plus,
        )
    }

    fn coords(a: &Alignment) -> Vec<(i64, i64)> {
        a.segments
            .iter()
            .map(|s| (s.coords.lend, s.coords.rend))
            .collect()
    }

    #[test]
    fn merge_extends_termini_and_keeps_shared_intron() {
        // shared intron 200..300; left terminus extends to min, right to max.
        let a = al(&[(120, 200), (300, 400)]);
        let b = al(&[(100, 200), (300, 450)]);
        let m = merge_alignments(&a, &b);
        assert_eq!(coords(&m), vec![(100, 200), (300, 450)]);
    }

    #[test]
    fn merge_appends_non_overlapping_tail_segment() {
        // b contributes a downstream exon a doesn't have.
        let a = al(&[(100, 200), (300, 400)]);
        let b = al(&[(100, 200), (300, 400), (500, 600)]);
        let m = merge_alignments(&a, &b);
        assert_eq!(coords(&m), vec![(100, 200), (300, 400), (500, 600)]);
    }

    #[test]
    fn merge_self_is_idempotent() {
        let a = al(&[(100, 200), (300, 400)]);
        let m = merge_alignments(&a, &a);
        assert_eq!(coords(&m), coords(&a));
    }
}
