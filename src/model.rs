//! Core data types shared across all stages.
//!
//! These mirror the C++ structures in `PASApipeline/pasa_cpp`
//! (`common_structs.h`, `alignment_segment.h`, `cdna_alignment.h`) and the Perl
//! alignment objects, but only the fields the two ported algorithms actually use.
//!
//! Coordinates are **1-based, inclusive** everywhere (as in PASA), with
//! `lend <= rend` after construction.

use std::sync::Arc;

/// A genomic coordinate interval, 1-based inclusive, `lend <= rend`.
///
/// The C++ analog is `struct coordset` (`common_structs.h`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Coordset {
    pub lend: i64,
    pub rend: i64,
}

impl Coordset {
    /// Construct, normalizing so `lend <= rend`.
    pub fn new(a: i64, b: i64) -> Self {
        if a <= b {
            Coordset { lend: a, rend: b }
        } else {
            Coordset { lend: b, rend: a }
        }
    }

    /// Inclusive overlap, exactly `common_subs.cpp` `overlap()`:
    /// `a.lend <= b.rend && a.rend >= b.lend`.
    ///
    /// This is the assembler's overlap test. A single shared base counts as
    /// overlap. **Do not** substitute [`Coordset::overlaps_strict`] here.
    pub fn overlaps_inclusive(&self, other: &Coordset) -> bool {
        self.lend <= other.rend && self.rend >= other.lend
    }

    /// Strict overlap, exactly the Perl alt-splice classifiers'
    /// `i_lend < j_rend && i_rend > j_lend`.
    ///
    /// A single shared base does **not** count. This is the alt-splice test.
    pub fn overlaps_strict(&self, other: &Coordset) -> bool {
        self.lend < other.rend && self.rend > other.lend
    }

    /// `true` if `other`'s span is fully within `self`'s span
    /// (`pasa_cpp/cdna_alignment_assembler.cpp::encapsulates`).
    pub fn encapsulates(&self, other: &Coordset) -> bool {
        other.lend >= self.lend && other.rend <= self.rend
    }

    /// Length in bp (inclusive): `rend - lend + 1`.
    pub fn len(&self) -> i64 {
        self.rend - self.lend + 1
    }

    /// Always `false` for a normalized interval (kept for clippy's `len`/`is_empty`
    /// pairing; a `Coordset` always spans at least one base).
    pub fn is_empty(&self) -> bool {
        self.rend < self.lend
    }

    /// The gap span to the **next** interval in a `lend`-sorted list:
    /// `(self.rend + 1, next.lend - 1)`. This is the project's intron / inter-segment
    /// convention (CLAUDE.md "Deliberate divergences") and the single source for it.
    /// It is intentionally **not** normalized — abutting/overlapping intervals yield an
    /// "empty" gap (`lend > rend`), exactly as the open-coded form did.
    pub fn gap_to(&self, next: &Coordset) -> Coordset {
        Coordset {
            lend: self.rend + 1,
            rend: next.lend - 1,
        }
    }

    /// The intersection of two intervals, or `None` when they are disjoint (share no base).
    pub fn intersect(&self, other: &Coordset) -> Option<Coordset> {
        let lend = self.lend.max(other.lend);
        let rend = self.rend.min(other.rend);
        (lend <= rend).then_some(Coordset { lend, rend })
    }

    /// Number of shared bases with `other` (`0` when disjoint).
    pub fn overlap_len(&self, other: &Coordset) -> i64 {
        (self.rend.min(other.rend) - self.lend.max(other.lend) + 1).max(0)
    }

    /// Reflect the interval about `axis`, mapping each genomic position `g` to
    /// `axis - g + 1` — used to transpose between forward and reverse-complement
    /// coordinates. The ends swap so the result stays `lend <= rend` when `axis >= rend`.
    pub fn reflect(&self, axis: i64) -> Coordset {
        Coordset {
            lend: axis - self.rend + 1,
            rend: axis - self.lend + 1,
        }
    }
}

/// The introns (gaps) of a `lend`-sorted interval list: `seg[i].gap_to(seg[i + 1])` for
/// each adjacent pair. The single source for the intron convention shared by the
/// assembler, alt-splice, ORF, and consensus paths.
pub fn introns_between(segs: &[Coordset]) -> impl Iterator<Item = Coordset> + '_ {
    segs.windows(2).map(|w| w[0].gap_to(&w[1]))
}

/// Transcribed / aligned orientation.
///
/// `Unknown` corresponds to PASA's `?` (and GFF/GTF `.`): an alignment whose
/// transcribed strand is not determined by the input. Only legal for
/// single-segment "flex" alignments handled by the orientation wrapper.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Strand {
    Plus,
    Minus,
    Unknown,
}

impl Strand {
    pub fn from_char(c: char) -> Strand {
        match c {
            '+' => Strand::Plus,
            '-' => Strand::Minus,
            _ => Strand::Unknown, // '?', '.', anything else
        }
    }

    pub fn to_char(self) -> char {
        match self {
            Strand::Plus => '+',
            Strand::Minus => '-',
            Strand::Unknown => '?',
        }
    }

    pub fn is_known(self) -> bool {
        !matches!(self, Strand::Unknown)
    }
}

/// Segment position within its alignment. Assigned purely by position in
/// [`Alignment::refine`] — the genome sequence is never consulted (see
/// `avoid-canonical-splice-bias`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SegType {
    Single,
    First,
    Internal,
    Last,
}

/// One aligned exon segment.
///
/// `left_splice` / `right_splice` mark whether that boundary is an internal
/// splice junction. Set **positionally** in [`Alignment::refine`]
/// (`pasa_cpp/cdna_alignment.cpp::refineAlignment`, lines 56–74):
/// single → neither, first → right only, last → left only, internal → both.
#[derive(Debug, Clone, PartialEq)]
pub struct Segment {
    pub coords: Coordset,
    pub seg_type: SegType,
    pub left_splice: bool,
    pub right_splice: bool,
    /// cDNA/transcript-space coords (forward axis), when known from GFF3/GTF.
    pub mcoords: Option<Coordset>,
    /// Per-segment percent identity, when present in the input.
    pub per_id: Option<f64>,
}

impl Segment {
    /// A bare segment with only genomic coords (the C++ token-format case).
    pub fn new(lend: i64, rend: i64) -> Self {
        Segment {
            coords: Coordset::new(lend, rend),
            seg_type: SegType::Single,
            left_splice: false,
            right_splice: false,
            mcoords: None,
            per_id: None,
        }
    }
}

/// Provenance: which input file a transcript came from. Assemblies carry the
/// union of their members' sources.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Provenance {
    pub source_file: Arc<str>,
}

impl Provenance {
    pub fn new(source_file: impl Into<Arc<str>>) -> Self {
        Provenance {
            source_file: source_file.into(),
        }
    }
}

/// A spliced cDNA/transcript alignment: an ordered list of exon segments on one
/// contig with an orientation. Mirrors `CDNA::CDNA_alignment` /
/// `pasa_cpp/cdna_alignment.cpp`.
#[derive(Debug, Clone)]
pub struct Alignment {
    pub acc: String,
    pub contig: String,
    pub segments: Vec<Segment>,
    /// Span: `(min lend, last-segment rend)` — see [`Alignment::refine`].
    pub coords: Coordset,
    /// Orientation used by the assembler (the "forced"/aligned orient).
    pub aligned_orient: Strand,
    /// Transcribed orientation as given by the input (`Unknown` if absent).
    pub spliced_orient: Strand,
    pub is_fl: bool,
    pub per_id: Option<f64>,
    pub percent_aligned: Option<f64>,
    pub provenance: Option<Provenance>,
    pub gene_id: Option<String>,
}

impl Alignment {
    /// Build from segments + orientation and immediately [`refine`](Self::refine).
    pub fn new(acc: impl Into<String>, segments: Vec<Segment>, aligned_orient: Strand) -> Self {
        let mut a = Alignment {
            acc: acc.into(),
            contig: String::new(),
            segments,
            coords: Coordset { lend: 0, rend: 0 },
            aligned_orient,
            spliced_orient: Strand::Unknown,
            is_fl: false,
            per_id: None,
            percent_aligned: None,
            provenance: None,
            gene_id: None,
        };
        a.refine();
        a
    }

    /// Port of `CDNA_alignment::refineAlignment` (`pasa_cpp/cdna_alignment.cpp`
    /// 40–75): sort segments by `lend`, set the span, and classify each
    /// segment's type + splice-junction flags **by position only**.
    ///
    /// Note: the span's `rend` is the **last segment's** `rend` (after sorting
    /// by `lend`), matching the C++ `seglist[size-1].get_coords().rend` — not
    /// the maximum `rend` across all segments.
    pub fn refine(&mut self) {
        let n = self.segments.len();
        debug_assert!(n > 0, "alignment with no segments");

        if n > 1 {
            // Stable sort by lend (Rust's sort is stable; the C++ std::sort is
            // not, but its comparator `a.lend <= b.lend` keeps already-sorted
            // input in place, so stable matches for well-formed exon lists).
            self.segments.sort_by_key(|s| s.coords.lend);
        }

        let min_lend = self.segments[0].coords.lend;
        let max_rend = self.segments[n - 1].coords.rend;
        self.coords = Coordset {
            lend: min_lend,
            rend: max_rend,
        };

        if n == 1 {
            let s = &mut self.segments[0];
            s.seg_type = SegType::Single;
            s.left_splice = false;
            s.right_splice = false;
            return;
        }

        for (i, seg) in self.segments.iter_mut().enumerate() {
            if i == 0 {
                seg.seg_type = SegType::First;
                seg.left_splice = false;
                seg.right_splice = true;
            } else if i == n - 1 {
                seg.seg_type = SegType::Last;
                seg.left_splice = true;
                seg.right_splice = false;
            } else {
                seg.seg_type = SegType::Internal;
                seg.left_splice = true;
                seg.right_splice = true;
            }
        }
    }

    pub fn num_segments(&self) -> usize {
        self.segments.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seg(lend: i64, rend: i64) -> Segment {
        Segment::new(lend, rend)
    }

    #[test]
    fn overlap_inclusive_vs_strict_single_base_touch() {
        let a = Coordset::new(100, 200);
        let b = Coordset::new(200, 300); // share exactly base 200
        assert!(a.overlaps_inclusive(&b), "inclusive counts a shared base");
        assert!(
            !a.overlaps_strict(&b),
            "strict does not count a shared base"
        );

        let c = Coordset::new(201, 300); // adjacent, no shared base
        assert!(!a.overlaps_inclusive(&c));
        assert!(!a.overlaps_strict(&c));

        let d = Coordset::new(150, 250); // genuine overlap
        assert!(a.overlaps_inclusive(&d));
        assert!(a.overlaps_strict(&d));
    }

    #[test]
    fn gap_to_is_the_intron_convention() {
        // adjacent exons → the gap between them, not normalized
        let a = Coordset::new(100, 200);
        let b = Coordset::new(301, 400);
        assert_eq!(
            a.gap_to(&b),
            Coordset {
                lend: 201,
                rend: 300
            }
        );
        // abutting exons → an "empty" (lend > rend) gap, preserved verbatim
        let touching = Coordset::new(201, 300);
        let g = a.gap_to(&touching);
        assert_eq!((g.lend, g.rend), (201, 200));
    }

    #[test]
    fn introns_between_walks_adjacent_pairs() {
        let exons = vec![
            Coordset::new(100, 200),
            Coordset::new(301, 400),
            Coordset::new(450, 500),
        ];
        let introns: Vec<_> = introns_between(&exons).collect();
        assert_eq!(
            introns,
            vec![Coordset::new(201, 300), Coordset::new(401, 449)]
        );
    }

    #[test]
    fn intersect_and_overlap_len() {
        let a = Coordset::new(100, 200);
        assert_eq!(
            a.intersect(&Coordset::new(150, 300)),
            Some(Coordset::new(150, 200))
        );
        assert_eq!(a.overlap_len(&Coordset::new(150, 300)), 51);
        // single shared base
        assert_eq!(
            a.intersect(&Coordset::new(200, 300)),
            Some(Coordset::new(200, 200))
        );
        assert_eq!(a.overlap_len(&Coordset::new(200, 300)), 1);
        // disjoint (adjacent, no shared base)
        assert_eq!(a.intersect(&Coordset::new(201, 300)), None);
        assert_eq!(a.overlap_len(&Coordset::new(201, 300)), 0);
    }

    #[test]
    fn reflect_is_its_own_inverse() {
        let c = Coordset::new(50, 120);
        let axis = 200;
        let r = c.reflect(axis); // hi - rend + 1 .. hi - lend + 1 = 81 .. 151
        assert_eq!(r, Coordset::new(81, 151));
        assert_eq!(
            r.reflect(axis),
            c,
            "reflecting twice about the same axis restores it"
        );
    }

    #[test]
    fn encapsulation_is_span_containment() {
        let a = Coordset::new(100, 500);
        assert!(a.encapsulates(&Coordset::new(100, 500))); // equal spans encapsulate
        assert!(a.encapsulates(&Coordset::new(200, 300)));
        assert!(!a.encapsulates(&Coordset::new(50, 300)));
        assert!(!a.encapsulates(&Coordset::new(300, 600)));
    }

    #[test]
    fn refine_single_segment() {
        let a = Alignment::new("x", vec![seg(100, 200)], Strand::Plus);
        assert_eq!(a.coords, Coordset::new(100, 200));
        assert_eq!(a.segments[0].seg_type, SegType::Single);
        assert!(!a.segments[0].left_splice);
        assert!(!a.segments[0].right_splice);
    }

    #[test]
    fn refine_classifies_positions_and_junctions() {
        // three exons, supplied out of order to exercise the sort
        let a = Alignment::new(
            "x",
            vec![seg(300, 400), seg(100, 150), seg(220, 260)],
            Strand::Plus,
        );
        assert_eq!(a.coords, Coordset::new(100, 400));
        let types: Vec<_> = a.segments.iter().map(|s| s.seg_type).collect();
        assert_eq!(
            types,
            vec![SegType::First, SegType::Internal, SegType::Last]
        );
        // first: right junction only
        assert!(!a.segments[0].left_splice && a.segments[0].right_splice);
        // internal: both
        assert!(a.segments[1].left_splice && a.segments[1].right_splice);
        // last: left junction only
        assert!(a.segments[2].left_splice && !a.segments[2].right_splice);
    }

    #[test]
    fn refine_span_rend_is_last_segment_rend() {
        // Last-by-lend segment is also last-by-rend here; matches C++ exactly.
        let a = Alignment::new("x", vec![seg(100, 150), seg(200, 999)], Strand::Plus);
        assert_eq!(a.coords.rend, 999);
    }

    #[test]
    fn strand_roundtrip() {
        assert_eq!(Strand::from_char('+'), Strand::Plus);
        assert_eq!(Strand::from_char('-'), Strand::Minus);
        assert_eq!(Strand::from_char('?'), Strand::Unknown);
        assert_eq!(Strand::from_char('.'), Strand::Unknown);
        assert_eq!(Strand::Plus.to_char(), '+');
        assert_eq!(Strand::Unknown.to_char(), '?');
    }
}
