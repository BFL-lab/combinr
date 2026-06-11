//! Pairwise compatibility + encapsulation tests.
//!
//! Direct port of `pasa_cpp/cdna_alignment_assembler.cpp`:
//! `canMerge` (lines 319–428) and `encapsulates` (640–658).

use crate::model::Alignment;

/// Port of `canMerge`. Two alignments are compatible (can be assembled together)
/// iff their spans overlap, they share orientation, and every overlapping
/// internal splice junction agrees (terminal exons tolerated within `fuzz`).
///
/// Junction flags are positional (set in `Alignment::refine`); no genome
/// sequence is consulted (see `avoid-canonical-splice-bias`).
pub fn can_merge(a1: &Alignment, a2: &Alignment, fuzz: i64) -> bool {
    // (a) spans must overlap (inclusive).
    if !a1.coords.overlaps_inclusive(&a2.coords) {
        return false;
    }
    // (b) same orientation.
    if a1.aligned_orient != a2.aligned_orient {
        return false;
    }

    let s1 = &a1.segments;
    let s2 = &a2.segments;

    // (c) phase: find the FIRST overlapping segment pair (a1 outer, a2 inner).
    let mut start1: Option<usize> = None;
    let mut start2: Option<usize> = None;
    'outer: for (i, seg1) in s1.iter().enumerate() {
        for (j, seg2) in s2.iter().enumerate() {
            if seg1.coords.overlaps_inclusive(&seg2.coords) {
                start1 = Some(i);
                start2 = Some(j);
                break 'outer;
            }
        }
    }
    let (mut i, mut j) = match (start1, start2) {
        (Some(i), Some(j)) => (i, j),
        _ => return false, // couldn't align two segments
    };
    // The first overlapping pair must include segment 0 of at least one side.
    if !(i == 0 || j == 0) {
        return false;
    }

    // (d) walk the overlapping segments in lockstep.
    while i < s1.len() && j < s2.len() {
        let seg1 = &s1[i];
        let seg2 = &s2[j];
        let (a1l, a1r) = (seg1.coords.lend, seg1.coords.rend);
        let (a2l, a2r) = (seg2.coords.lend, seg2.coords.rend);

        if seg1.coords.overlaps_inclusive(&seg2.coords) {
            // Left splice junction. The C++ uses an `if / else if / else if`
            // chain that all `return false`; collapsed here into one
            // short-circuiting condition (logically identical):
            //   - both sides are junctions but their left coords differ; or
            //   - only one side is a junction and the other terminus is more
            //     than `fuzz` short of it.
            if seg1.left_splice || seg2.left_splice {
                let conflict = (seg1.left_splice && seg2.left_splice && a1l != a2l)
                    || (seg1.left_splice && (a2l + fuzz < a1l))
                    || (seg2.left_splice && (a1l + fuzz < a2l));
                if conflict {
                    return false;
                }
            }
            // Right splice junction (mirror of the left).
            if seg1.right_splice || seg2.right_splice {
                let conflict = (seg1.right_splice && seg2.right_splice && a1r != a2r)
                    || (seg1.right_splice && (a2r - fuzz > a1r))
                    || (seg2.right_splice && (a1r - fuzz > a2r));
                if conflict {
                    return false;
                }
            }
        } else {
            // two ordered segments do not overlap each other
            return false;
        }
        i += 1;
        j += 1;
    }

    true
}

/// Port of `encapsulates`: `true` if `b`'s span is fully within `a`'s span.
pub fn encapsulates(a: &Alignment, b: &Alignment) -> bool {
    a.coords.encapsulates(&b.coords)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Alignment, Segment, Strand};

    fn al(segs: &[(i64, i64)], orient: Strand) -> Alignment {
        Alignment::new(
            "x",
            segs.iter().map(|&(l, r)| Segment::new(l, r)).collect(),
            orient,
        )
    }

    #[test]
    fn identical_internal_junctions_merge() {
        // share the internal intron exactly; terminal exons differ in length.
        let a = al(&[(100, 200), (300, 400)], Strand::Plus);
        let b = al(&[(120, 200), (300, 450)], Strand::Plus);
        assert!(can_merge(&a, &b, 20));
    }

    #[test]
    fn differing_internal_junction_blocks_merge() {
        // a's first exon ends at 200 (right junction), b's at 205 → conflict.
        let a = al(&[(100, 200), (300, 400)], Strand::Plus);
        let b = al(&[(100, 205), (300, 400)], Strand::Plus);
        assert!(!can_merge(&a, &b, 20));
    }

    #[test]
    fn opposite_orientation_blocks_merge() {
        let a = al(&[(100, 200), (300, 400)], Strand::Plus);
        let b = al(&[(100, 200), (300, 400)], Strand::Minus);
        assert!(!can_merge(&a, &b, 20));
    }

    #[test]
    fn non_overlapping_spans_block_merge() {
        let a = al(&[(100, 200)], Strand::Plus);
        let b = al(&[(500, 600)], Strand::Plus);
        assert!(!can_merge(&a, &b, 20));
    }

    #[test]
    fn one_sided_junction_fuzz_boundary() {
        // a: single exon ending 220 (no junction, terminal).
        // b: first exon ends at 200 (right junction). a's rend extends past b's.
        // Rule (a1 here = a, has NO right junction; a2 = b, HAS right junction):
        //   `a2.right && (a1_rend - fuzz > a2_rend)` → blocks when 220 - fuzz > 200.
        // fuzz=20 → 200 > 200 false → compatible (exactly fuzz).
        // fuzz=19 → 201 > 200 true  → blocked.
        let a = al(&[(150, 220)], Strand::Plus);
        let b = al(&[(150, 200), (300, 400)], Strand::Plus);
        assert!(can_merge(&a, &b, 20), "exactly fuzz distance is compatible");
        assert!(!can_merge(&a, &b, 19), "one bp beyond fuzz blocks");
        assert!(can_merge(&a, &b, 21));
    }

    #[test]
    fn encapsulation_span() {
        let a = al(&[(100, 500)], Strand::Plus);
        let b = al(&[(200, 300)], Strand::Plus);
        assert!(encapsulates(&a, &b));
        assert!(!encapsulates(&b, &a));
    }
}
