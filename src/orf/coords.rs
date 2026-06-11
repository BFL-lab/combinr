//! Transcript ↔ genome coordinate projection for a spliced isoform, and ORF
//! projection from a start codon. Concept ported from
//! `PerlLib/Exons_to_geneobj.pm` (mapping a CDS range onto exon structure).
//!
//! Transcript coordinates are 0-based offsets along the spliced mRNA, 5'→3'.
//! Genomic coordinates are 1-based inclusive.

use crate::io::fasta::Fasta;
use crate::model::{Coordset, Strand};
use crate::orf::translate::{is_stop, reverse_complement};

struct ExonPiece {
    genomic: Coordset, // lend <= rend
    t_start: usize,    // transcript offset of this piece's 5' base
    len: usize,
}

/// A spliced transcript: its exon pieces in transcription order plus the
/// genome↔transcript coordinate maps.
pub struct SplicedTranscript {
    plus: bool,
    pieces: Vec<ExonPiece>,
    total_len: usize,
}

/// Result of projecting an ORF from a start codon along a transcript.
#[derive(Debug, Clone, Copy)]
pub struct OrfProjection {
    pub cds_t_start: usize,
    pub cds_t_end: usize, // exclusive; includes the stop codon when one is hit
    pub hit_stop: bool,
}

impl SplicedTranscript {
    /// Build from genomic exons (any order) and a strand. `Unknown` is treated
    /// as `+`.
    pub fn new(exons: &[Coordset], strand: Strand) -> SplicedTranscript {
        let mut ex: Vec<Coordset> = exons.to_vec();
        ex.sort_by_key(|c| c.lend);
        let plus = strand != Strand::Minus;
        if !plus {
            ex.reverse(); // transcription order is descending genomic on '-'
        }
        let mut pieces = Vec::with_capacity(ex.len());
        let mut t = 0usize;
        for c in ex {
            let len = c.len() as usize;
            pieces.push(ExonPiece {
                genomic: c,
                t_start: t,
                len,
            });
            t += len;
        }
        SplicedTranscript {
            plus,
            pieces,
            total_len: t,
        }
    }

    pub fn len(&self) -> usize {
        self.total_len
    }

    pub fn is_empty(&self) -> bool {
        self.total_len == 0
    }

    /// Transcript offset of a genomic coordinate, or `None` if it lies in no exon.
    pub fn genomic_to_tpos(&self, g: i64) -> Option<usize> {
        for p in &self.pieces {
            if g >= p.genomic.lend && g <= p.genomic.rend {
                let off = if self.plus {
                    (g - p.genomic.lend) as usize
                } else {
                    (p.genomic.rend - g) as usize
                };
                return Some(p.t_start + off);
            }
        }
        None
    }

    /// Genomic coordinate of a transcript offset (must be in range).
    pub fn tpos_to_genomic(&self, t: usize) -> Option<i64> {
        for p in &self.pieces {
            if t >= p.t_start && t < p.t_start + p.len {
                let off = (t - p.t_start) as i64;
                return Some(if self.plus {
                    p.genomic.lend + off
                } else {
                    p.genomic.rend - off
                });
            }
        }
        None
    }

    /// The spliced cDNA sequence (5'→3', uppercased) read from the genome, or
    /// `None` if any exon is out of bounds.
    pub fn sequence(&self, fa: &Fasta, contig: &str) -> Option<Vec<u8>> {
        let mut out = Vec::with_capacity(self.total_len);
        for p in &self.pieces {
            let chunk = fa.subseq(contig, p.genomic.lend, p.genomic.rend)?;
            if self.plus {
                out.extend_from_slice(chunk);
            } else {
                out.extend(reverse_complement(chunk));
            }
        }
        Some(out)
    }

    /// Genomic exon pieces covering the transcript half-open span `[t_lo, t_hi)`,
    /// returned as `lend`-sorted genomic coordsets.
    pub fn genomic_segments_for_tspan(&self, t_lo: usize, t_hi: usize) -> Vec<Coordset> {
        let mut out = Vec::new();
        if t_hi <= t_lo {
            return out;
        }
        for p in &self.pieces {
            let p_lo = p.t_start;
            let p_hi = p.t_start + p.len;
            let a = t_lo.max(p_lo);
            let b = t_hi.min(p_hi);
            if a >= b {
                continue;
            }
            // map transcript sub-interval [a, b) to genomic
            let (g1, g2) = if self.plus {
                (
                    p.genomic.lend + (a - p.t_start) as i64,
                    p.genomic.lend + (b - 1 - p.t_start) as i64,
                )
            } else {
                (
                    p.genomic.rend - (b - 1 - p.t_start) as i64,
                    p.genomic.rend - (a - p.t_start) as i64,
                )
            };
            out.push(Coordset::new(g1, g2));
        }
        out.sort_by_key(|c| c.lend);
        out
    }

    /// Project an ORF starting at transcript offset `cds_t_start` (the 5' base of
    /// the start codon), translating `seq` until the first in-frame stop. With no
    /// stop, the CDS runs to the last complete codon (3'-partial).
    pub fn project_orf(&self, seq: &[u8], cds_t_start: usize) -> OrfProjection {
        let mut t = cds_t_start;
        let mut hit_stop = false;
        let mut cds_t_end = self.total_len;
        while t + 3 <= seq.len() {
            if is_stop(&seq[t..t + 3]) {
                cds_t_end = t + 3; // include the stop codon
                hit_stop = true;
                break;
            }
            t += 3;
        }
        if !hit_stop {
            let n_codons = (seq.len().saturating_sub(cds_t_start)) / 3;
            cds_t_end = cds_t_start + n_codons * 3;
        }
        OrfProjection {
            cds_t_start,
            cds_t_end,
            hit_stop,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cs(l: i64, r: i64) -> Coordset {
        Coordset::new(l, r)
    }

    #[test]
    fn plus_strand_mapping_and_segments() {
        // two exons: 100-104 (len 5), 200-204 (len 5); transcript 0..10
        let t = SplicedTranscript::new(&[cs(100, 104), cs(200, 204)], Strand::Plus);
        assert_eq!(t.len(), 10);
        assert_eq!(t.genomic_to_tpos(100), Some(0));
        assert_eq!(t.genomic_to_tpos(104), Some(4));
        assert_eq!(t.genomic_to_tpos(200), Some(5));
        assert_eq!(t.genomic_to_tpos(204), Some(9));
        assert_eq!(t.genomic_to_tpos(150), None); // intronic
        assert_eq!(t.tpos_to_genomic(0), Some(100));
        assert_eq!(t.tpos_to_genomic(5), Some(200));
        // transcript [3, 7) spans the junction: 103-104 and 200-201
        let segs = t.genomic_segments_for_tspan(3, 7);
        assert_eq!(segs, vec![cs(103, 104), cs(200, 201)]);
    }

    #[test]
    fn minus_strand_mapping_is_reversed() {
        // two exons; on '-' transcription starts at the rightmost base.
        let t = SplicedTranscript::new(&[cs(100, 104), cs(200, 204)], Strand::Minus);
        assert_eq!(t.len(), 10);
        // 5' base is 204
        assert_eq!(t.genomic_to_tpos(204), Some(0));
        assert_eq!(t.genomic_to_tpos(200), Some(4));
        assert_eq!(t.genomic_to_tpos(104), Some(5));
        assert_eq!(t.genomic_to_tpos(100), Some(9));
        assert_eq!(t.tpos_to_genomic(0), Some(204));
        assert_eq!(t.tpos_to_genomic(9), Some(100));
        // transcript [3, 7): 201-204 (3 bases of first exon) wait -> compute
        let segs = t.genomic_segments_for_tspan(3, 7);
        // t3..t6 inclusive: t3=201,t4=200 (exon2), t5=104,t6=103 (exon1)
        assert_eq!(segs, vec![cs(103, 104), cs(200, 201)]);
    }

    #[test]
    fn project_orf_hits_stop() {
        // single exon, sequence ATG AAA TAA ...
        let t = SplicedTranscript::new(&[cs(1, 12)], Strand::Plus);
        let seq = b"ATGAAATAACC";
        let p = t.project_orf(seq, 0);
        assert!(p.hit_stop);
        assert_eq!(p.cds_t_start, 0);
        assert_eq!(p.cds_t_end, 9); // includes TAA
    }

    #[test]
    fn project_orf_runs_off_end_when_no_stop() {
        let t = SplicedTranscript::new(&[cs(1, 11)], Strand::Plus);
        let seq = b"ATGAAACCCGG"; // 11 nt, no stop
        let p = t.project_orf(seq, 0);
        assert!(!p.hit_stop);
        // last complete codon: 11/3 = 3 codons → 9
        assert_eq!(p.cds_t_end, 9);
    }
}
