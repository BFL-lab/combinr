//! Candidate exons and the EVM frame model.
//!
//! EVM types each candidate exon as initial/internal/terminal/single, on a strand,
//! with a reading frame encoded 1-6 (forward 1/2/3, reverse 4/5/6 — see
//! `evidence_modeler.pl` `Exon` package ~4657 and the grammar in `initialization`
//! ~2445). During a strand pass the engine works in forward frames 1-3; the minus
//! pass reverse-complements the region, runs the same machinery, then transposes
//! (frames 1->4/2->5/3->6) — that transpose lands in M3.
//!
//! Non-canonical (see `avoid-canonical-splice-bias`): [`determine_good_phases`]
//! scans the supplied subsequence for in-frame stops under the *configurable*
//! [`GeneticCode`] instead of consulting a pre-scanned canonical STOP vector, and no
//! GT-AG/ATG gating is applied to exon boundaries.

use crate::model::{Coordset, Strand};
use crate::orf::GeneticCode;

/// Position/role of a candidate exon within a gene.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum ExonType {
    Initial,
    Internal,
    Terminal,
    Single,
    /// A trellis boundary sentinel (added in M2).
    Bound,
}

impl ExonType {
    pub fn as_str(self) -> &'static str {
        match self {
            ExonType::Initial => "initial",
            ExonType::Internal => "internal",
            ExonType::Terminal => "terminal",
            ExonType::Single => "single",
            ExonType::Bound => "bound",
        }
    }
}

/// One candidate exon: a typed, framed genomic interval with the evidence that
/// supports it and the 2-bp sequence boundaries used for the junction stop-check.
#[derive(Clone, Debug, PartialEq)]
pub struct ExonCandidate {
    pub coords: Coordset,
    pub orient: Strand,
    pub exon_type: ExonType,
    /// Reading frame at the 5' end, encoded 1-6 (forward 1-3, reverse 4-6).
    pub start_frame: u8,
    /// Reading frame state at the 3' end, encoded 1-6.
    pub end_frame: u8,
    /// First two bases of the exon (genomic 5'->3' on `orient`), for the across-junction
    /// stop-codon test in `are_compatible_exons` (M2).
    pub left_seq_boundary: [u8; 2],
    /// Last two bases of the exon.
    pub right_seq_boundary: [u8; 2],
    /// `(accession, ev_type)` of every evidence chain supporting this exon.
    pub evidence: Vec<(String, String)>,
}

/// EVM `Exon::setStartFrame` (~line 4740): `end_frame = ((len + start_frame - 1) % 3)`,
/// with a result of 0 mapped to 3. `start_frame` is a forward frame in `1..=3`.
pub fn end_frame_for(start_frame: u8, len: i64) -> u8 {
    let m = (len + start_frame as i64 - 1).rem_euclid(3);
    if m == 0 { 3 } else { m as u8 }
}

/// Map a forward frame (`1..=3`) to its reverse-strand encoding (`4..=6`).
pub fn to_reverse_frame(f: u8) -> u8 {
    f + 3
}

/// Forward reading frames (1,2,3) in which `seq` contains no in-frame stop codon
/// under `code`. Port of `determine_good_phases` (~line 1616): a stop starting at
/// offset `off` kills the phase given by `off % 3` (0->phase 1, 1->phase 3, 2->phase 2).
///
/// Callers pass the exon subsequence to check; for terminal/single exons they trim the
/// trailing stop codon (3 bp) first so the gene's own stop isn't counted.
pub fn determine_good_phases(seq: &[u8], code: &GeneticCode) -> Vec<u8> {
    let mut ok = [true; 3]; // index = phase - 1
    if seq.len() >= 3 {
        for off in 0..=(seq.len() - 3) {
            if code.is_stop(&seq[off..off + 3]) {
                let phase = match off % 3 {
                    0 => 1usize,
                    1 => 3,
                    _ => 2,
                };
                ok[phase - 1] = false;
            }
        }
    }
    (1u8..=3).filter(|&p| ok[(p - 1) as usize]).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn end_frame_matches_evm_formula() {
        // len 3 (one full codon), frame 1 -> (3+1-1)%3 = 0 -> 3 (connects to start 1)
        assert_eq!(end_frame_for(1, 3), 3);
        // len 1, frame 1 -> (1+1-1)%3 = 1
        assert_eq!(end_frame_for(1, 1), 1);
        // len 2, frame 1 -> (2)%3 = 2
        assert_eq!(end_frame_for(1, 2), 2);
        // len 100, frame 2 -> (100+2-1)%3 = 101%3 = 2
        assert_eq!(end_frame_for(2, 100), 2);
    }

    #[test]
    fn reverse_frame_remap() {
        assert_eq!(to_reverse_frame(1), 4);
        assert_eq!(to_reverse_frame(3), 6);
    }

    #[test]
    fn good_phases_no_stop() {
        let code = GeneticCode::default();
        // free of TAA/TAG/TGA at every offset and frame
        assert_eq!(determine_good_phases(b"AAAGGGCCC", &code), vec![1, 2, 3]);
    }

    #[test]
    fn good_phases_stop_in_phase1() {
        let code = GeneticCode::default();
        // TAA at offset 3 (off%3==0 -> kills phase 1)
        assert_eq!(determine_good_phases(b"ATGTAAGGG", &code), vec![2, 3]);
    }

    #[test]
    fn good_phases_respects_genetic_code() {
        // Ciliate code 6: TAA/TAG are NOT stops, so a TAA no longer kills a phase.
        let c6 = GeneticCode::from_ncbi_id(6).unwrap();
        assert_eq!(determine_good_phases(b"ATGTAAGGG", &c6), vec![1, 2, 3]);
        // ...but TGA still stops under code 6 (TGA at offset 3 -> kills phase 1).
        assert_eq!(determine_good_phases(b"ATGTGAGGG", &c6), vec![2, 3]);
    }

    #[test]
    fn good_phases_short_seq_all_ok() {
        let code = GeneticCode::default();
        assert_eq!(determine_good_phases(b"AT", &code), vec![1, 2, 3]);
    }
}
