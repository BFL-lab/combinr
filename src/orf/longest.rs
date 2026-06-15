//! De-novo longest-ORF finding on a spliced cDNA sequence.
//!
//! PASA/EVidenceModeler-style (reference: `EVidenceModeler/PerlLib/Longest_orf.pm`
//! `capture_all_ORFs` / `identify_putative_starts`). The transcript is already stranded
//! (the caller passes the spliced 5'→3' cDNA), so only the three forward frames are
//! scanned. Stops use the configurable [`GeneticCode`]; starts are `ATG` (PASA
//! convention). 5'-partial (no upstream `ATG`) and 3'-partial (no stop) ORFs are allowed,
//! so a coding transcript whose ends weren't fully captured still yields an ORF.

use crate::orf::translate::GeneticCode;

/// The coding span of an ORF, in transcript offsets (0-based, half-open). `t_end`
/// includes the stop codon when `has_stop`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OrfSpan {
    pub t_start: usize,
    pub t_end: usize,
    pub has_start: bool,
    pub has_stop: bool,
}

impl OrfSpan {
    pub fn len(&self) -> usize {
        self.t_end - self.t_start
    }

    pub fn is_empty(&self) -> bool {
        self.t_end <= self.t_start
    }
}

/// The single longest ORF in `seq` (uppercase nucleotides), or `None` if shorter than a
/// codon. Ties break toward a complete (start+stop) ORF for determinism.
pub fn find_longest_orf(seq: &[u8], code: &GeneticCode) -> Option<OrfSpan> {
    let n = seq.len();
    if n < 3 {
        return None;
    }
    let mut cands: Vec<OrfSpan> = Vec::new();

    for frame in 0..3 {
        if frame + 3 > n {
            continue;
        }
        let last_codon_end = frame + ((n - frame) / 3) * 3; // exclusive end of last full codon
        let mut region_start = frame; // start of the current inter-stop region
        let mut first_atg: Option<usize> = None;
        let mut pos = frame;

        while pos + 3 <= n {
            let codon = &seq[pos..pos + 3];
            if first_atg.is_none() && codon == b"ATG" {
                first_atg = Some(pos);
            }
            if code.is_stop(codon) {
                let stop_end = pos + 3;
                // complete ORF: first ATG of this region -> this stop
                if let Some(a) = first_atg {
                    push(&mut cands, a, stop_end, true, true);
                }
                // 5'-partial: only the first region (runs off the transcript's 5' end)
                if region_start == frame && first_atg.is_none() {
                    push(&mut cands, frame, stop_end, false, true);
                }
                region_start = stop_end;
                first_atg = None;
            }
            pos += 3;
        }

        // trailing region (no stop): 3'-partial
        if let Some(a) = first_atg {
            push(&mut cands, a, last_codon_end, true, false);
        } else if region_start == frame {
            // whole frame: no stop and no ATG -> 5'+3' partial
            push(&mut cands, frame, last_codon_end, false, false);
        }
    }

    cands.into_iter().max_by(|a, b| {
        a.len()
            .cmp(&b.len())
            .then((a.has_start && a.has_stop).cmp(&(b.has_start && b.has_stop)))
    })
}

fn push(cands: &mut Vec<OrfSpan>, t_start: usize, t_end: usize, has_start: bool, has_stop: bool) {
    if t_end > t_start {
        cands.push(OrfSpan {
            t_start,
            t_end,
            has_start,
            has_stop,
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn code() -> GeneticCode {
        GeneticCode::default()
    }

    #[test]
    fn complete_orf_atg_to_stop() {
        // ATG AAA AAA TAA  (CDS 0..12, includes stop)
        let orf = find_longest_orf(b"ATGAAAAAATAA", &code()).unwrap();
        assert_eq!((orf.t_start, orf.t_end), (0, 12));
        assert!(orf.has_start && orf.has_stop);
    }

    #[test]
    fn picks_longest_across_frames() {
        // short ORF in frame 0, longer one in frame 1
        // frame0: ATG TAA (2 codons -> CDS 0..6); insert a long frame-1 ORF after
        let seq = b"ATGTAA" // frame0 ATG..stop (len 6)
            .iter()
            .chain(b"AATGAAAAAAAAAAAATAAA") // frame? craft a longer ORF
            .copied()
            .collect::<Vec<u8>>();
        let orf = find_longest_orf(&seq, &code()).unwrap();
        assert!(orf.len() > 6, "longest ORF wins, got len {}", orf.len());
    }

    #[test]
    fn three_prime_partial_when_no_stop() {
        // ATG AAA AAA  (no stop) -> 3'-partial, runs to last full codon
        let orf = find_longest_orf(b"ATGAAAAAA", &code()).unwrap();
        assert_eq!((orf.t_start, orf.t_end), (0, 9));
        assert!(orf.has_start && !orf.has_stop);
    }

    #[test]
    fn five_prime_partial_when_no_atg() {
        // AAA AAA TAA — no ATG, but a stop in frame 0 -> 5'-partial 0..9
        let orf = find_longest_orf(b"AAAAAATAA", &code()).unwrap();
        assert_eq!((orf.t_start, orf.t_end), (0, 9));
        assert!(!orf.has_start && orf.has_stop);
    }

    #[test]
    fn respects_genetic_code() {
        // ATG AAA TAA AAA: under code 1 the TAA stops the ORF at 9; under code 6
        // (TAA = Gln) translation reads through to the end (3'-partial, len 12).
        let seq = b"ATGAAATAAAAA";
        let c1 = find_longest_orf(seq, &GeneticCode::from_ncbi_id(1).unwrap()).unwrap();
        assert_eq!((c1.t_start, c1.t_end, c1.has_stop), (0, 9, true));
        let c6 = find_longest_orf(seq, &GeneticCode::from_ncbi_id(6).unwrap()).unwrap();
        assert_eq!((c6.t_start, c6.t_end, c6.has_stop), (0, 12, false));
    }

    #[test]
    fn too_short_is_none() {
        assert!(find_longest_orf(b"AT", &code()).is_none());
    }
}
