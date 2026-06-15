//! The gene-structure grammar: which exon→exon transitions are legal, and which are
//! phased (intron) vs unphased (intergenic).
//!
//! Direct port of EVM's `%ACCEPTABLE_EXON_LINKAGES`, `%PHASED_CONNECTIONS` and
//! `%INTERGENIC_CONNECTIONS` tables (`initialization`, evidence_modeler.pl ~2467) plus
//! the internal-exon frame-transition pairs (~2511). A single trellis spans both
//! strands, so the tables include forward→reverse and reverse→forward gene transitions.

use crate::consensus::exon::ExonType;
use crate::model::Strand;

/// How two consecutive exons may connect.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Linkage {
    /// Within a gene, across an intron (frames + intron evidence checked).
    Phased,
    /// Between genes (or strands), across an intergenic gap.
    Intergenic,
    /// Not a legal transition.
    Incompatible,
}

/// Classify the transition from exon `a` to exon `b` (a precedes b), each given as
/// `(type, orientation)`. Mirrors the `%ACCEPTABLE_EXON_LINKAGES` entries: those with
/// the phased flag are [`Linkage::Phased`], the rest are [`Linkage::Intergenic`].
pub fn classify(a: (ExonType, Strand), b: (ExonType, Strand)) -> Linkage {
    use ExonType::{Initial, Internal, Single, Terminal};
    use Linkage::{Incompatible, Intergenic, Phased};
    use Strand::{Minus as M, Plus as P};
    match (a, b) {
        // forward strand, within-gene (phased / intron)
        ((Initial, P), (Terminal, P)) => Phased,
        ((Initial, P), (Internal, P)) => Phased,
        ((Internal, P), (Internal, P)) => Phased,
        ((Internal, P), (Terminal, P)) => Phased,
        // reverse strand, within-gene (phased / intron)
        ((Terminal, M), (Initial, M)) => Phased,
        ((Internal, M), (Initial, M)) => Phased,
        ((Internal, M), (Internal, M)) => Phased,
        ((Terminal, M), (Internal, M)) => Phased,
        // forward strand, between-gene (intergenic)
        ((Terminal, P), (Initial, P)) => Intergenic,
        ((Terminal, P), (Single, P)) => Intergenic,
        ((Single, P), (Single, P)) => Intergenic,
        ((Single, P), (Initial, P)) => Intergenic,
        // reverse strand, between-gene (intergenic)
        ((Initial, M), (Terminal, M)) => Intergenic,
        ((Single, M), (Terminal, M)) => Intergenic,
        ((Single, M), (Single, M)) => Intergenic,
        ((Initial, M), (Single, M)) => Intergenic,
        // forward -> reverse transitions
        ((Single, P), (Single, M)) => Intergenic,
        ((Single, P), (Terminal, M)) => Intergenic,
        ((Terminal, P), (Terminal, M)) => Intergenic,
        ((Terminal, P), (Single, M)) => Intergenic,
        // reverse -> forward transitions
        ((Single, M), (Single, P)) => Intergenic,
        ((Single, M), (Initial, P)) => Intergenic,
        ((Initial, M), (Initial, P)) => Intergenic,
        ((Initial, M), (Single, P)) => Intergenic,
        _ => Incompatible,
    }
}

/// Frame transition across an intron is legal: forward `1→2→3→1`, reverse `4→5→6→4`
/// (the `[1,2],[2,3],[3,1],[4,5],[5,6],[6,4]` pairs in `%ACCEPTABLE_EXON_LINKAGES`).
pub fn frame_transition_ok(end_a: u8, start_b: u8) -> bool {
    matches!(
        (end_a, start_b),
        (1, 2) | (2, 3) | (3, 1) | (4, 5) | (5, 6) | (6, 4)
    )
}

/// Whether an exon ends a gene in `traverse_path`'s split rule
/// (`terminal+ | single | initial-`, evidence_modeler.pl ~2803).
pub fn is_gene_boundary(exon_type: ExonType, orient: Strand) -> bool {
    matches!(
        (exon_type, orient),
        (ExonType::Terminal, Strand::Plus)
            | (ExonType::Single, _)
            | (ExonType::Initial, Strand::Minus)
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use ExonType::{Initial, Internal, Single, Terminal};
    use Strand::{Minus as M, Plus as P};

    #[test]
    fn forward_within_gene_is_phased() {
        assert_eq!(classify((Initial, P), (Internal, P)), Linkage::Phased);
        assert_eq!(classify((Internal, P), (Internal, P)), Linkage::Phased);
        assert_eq!(classify((Internal, P), (Terminal, P)), Linkage::Phased);
        assert_eq!(classify((Initial, P), (Terminal, P)), Linkage::Phased);
    }

    #[test]
    fn forward_between_gene_is_intergenic() {
        assert_eq!(classify((Terminal, P), (Initial, P)), Linkage::Intergenic);
        assert_eq!(classify((Single, P), (Single, P)), Linkage::Intergenic);
    }

    #[test]
    fn strand_switch_is_intergenic() {
        assert_eq!(classify((Terminal, P), (Terminal, M)), Linkage::Intergenic);
        assert_eq!(classify((Initial, M), (Initial, P)), Linkage::Intergenic);
    }

    #[test]
    fn illegal_transitions_are_incompatible() {
        // can't go initial -> initial within a gene, or terminal -> internal forward
        assert_eq!(classify((Initial, P), (Initial, P)), Linkage::Incompatible);
        assert_eq!(
            classify((Terminal, P), (Internal, P)),
            Linkage::Incompatible
        );
        // wrong-strand within-gene
        assert_eq!(classify((Initial, P), (Internal, M)), Linkage::Incompatible);
    }

    #[test]
    fn frame_cycle() {
        assert!(
            frame_transition_ok(1, 2) && frame_transition_ok(2, 3) && frame_transition_ok(3, 1)
        );
        assert!(
            frame_transition_ok(4, 5) && frame_transition_ok(5, 6) && frame_transition_ok(6, 4)
        );
        assert!(!frame_transition_ok(1, 3));
        assert!(!frame_transition_ok(1, 1));
        assert!(!frame_transition_ok(3, 4)); // no cross-strand frame transition
    }

    #[test]
    fn gene_boundaries() {
        assert!(is_gene_boundary(Terminal, P));
        assert!(is_gene_boundary(Single, P));
        assert!(is_gene_boundary(Single, M));
        assert!(is_gene_boundary(Initial, M));
        assert!(!is_gene_boundary(Initial, P));
        assert!(!is_gene_boundary(Internal, P));
        assert!(!is_gene_boundary(Terminal, M));
    }
}
