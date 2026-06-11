//! BAM input: spliced cDNA/transcript-to-genome alignments.
//!
//! Parses a BAM produced by a splice-aware mapper (e.g. segemehl,
//! `minimap2 -ax splice`, STAR, HISAT2) into the same [`Alignment`] structures
//! the GTF/GFF3 readers produce, so everything downstream of `load_sources` is
//! format-agnostic. The model is cDNA-to-genome (few alignments per locus), so
//! a single linear pass over the records is sufficient.
//!
//! Conventions:
//! - One primary mapped record → one [`Alignment`]. Unmapped, secondary, and
//!   supplementary records are skipped.
//! - Exon blocks come from the CIGAR: a run of reference-consuming non-skip ops
//!   (`M`/`=`/`X`/`D`) is one exon; `N` (skip) splits exons; `I` stays within an
//!   exon; soft/hard clips and pads do not affect the genomic span. Only `N`
//!   introduces an intron.
//! - `aligned_orient` is the FLAG strand (`0x10`). `spliced_orient` — the
//!   transcribed strand that drives output — is read from the `XS:A` tag
//!   (segemehl/HISAT2/STAR convention: absolute genomic strand), or `Unknown`
//!   when absent. The genome is never consulted and splice motifs are never
//!   re-derived (see `avoid-canonical-splice-bias`).
//! - `per_id` is left `None`: BAM carries `NM` (edit distance), not PASA's
//!   percent identity, and there is no GTF analog to match against, so the
//!   optional `--min-avg-per-id` filter simply never evaluates BAM input.
//! - CRAM and plain-text SAM are out of scope.

use std::path::Path;
use std::sync::Arc;

use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Value;

use super::{RawSegment, build_alignment};
use crate::error::{CombinrError, Result};
use crate::model::{Alignment, Coordset, Strand};

/// Parse every primary mapped record in `path` into an [`Alignment`], stamped
/// with `source` as provenance.
pub fn parse(path: &Path, source: &str) -> Result<Vec<Alignment>> {
    let src: Arc<str> = Arc::from(source);
    let mut reader = bam::io::reader::Builder.build_from_path(path)?;
    let header = reader.read_header()?;

    let mut out: Vec<Alignment> = Vec::new();
    for result in reader.records() {
        let record = result?;
        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }

        let acc = record
            .name()
            .map(|n| n.to_string())
            .ok_or_else(|| parse_err(source, "mapped record has no QNAME".to_string()))?;

        let ref_id = match record.reference_sequence_id() {
            Some(r) => r?,
            None => continue, // mapped but no reference id; skip defensively
        };
        let contig = header
            .reference_sequences()
            .get_index(ref_id)
            .map(|(name, _)| name.to_string())
            .ok_or_else(|| {
                parse_err(
                    source,
                    format!("{acc}: reference id {ref_id} not in header"),
                )
            })?;

        let pos = record
            .alignment_start()
            .transpose()?
            .ok_or_else(|| parse_err(source, format!("{acc}: mapped record has no POS")))?;
        let pos_1based = usize::from(pos) as i64;

        let mut ops: Vec<(Kind, i64)> = Vec::new();
        for op in record.cigar().iter() {
            let op = op?;
            ops.push((op.kind(), op.len() as i64));
        }
        let exons = cigar_to_exons(pos_1based, &ops);
        if exons.is_empty() {
            continue; // no aligned block (e.g. all clips); nothing to model
        }

        let aligned = if flags.is_reverse_complemented() {
            Strand::Minus
        } else {
            Strand::Plus
        };
        let spliced = xs_strand(&record);

        let raw: Vec<RawSegment> = exons
            .into_iter()
            .map(|g| RawSegment {
                genomic: g,
                mcoords: None,
                per_id: None,
            })
            .collect();

        // `build_alignment` collapses aligned/spliced into one strand and uses
        // it to order the synthesized mcoords. Drive that order by the spliced
        // (transcribed) strand when known, else by the FLAG strand; then restore
        // the two distinct orientations. Downstream only consults aligned_orient
        // when spliced_orient is Unknown, so this faithfully separates them.
        let order_strand = if spliced.is_known() { spliced } else { aligned };
        let mut a = build_alignment(acc, contig, order_strand, None, raw, &src);
        a.spliced_orient = spliced;
        a.aligned_orient = aligned;
        out.push(a);
    }

    // Deterministic order, matching the text parsers: (contig, leftmost, acc).
    out.sort_by(|a, b| {
        (a.contig.as_str(), a.coords.lend, a.acc.as_str()).cmp(&(
            b.contig.as_str(),
            b.coords.lend,
            b.acc.as_str(),
        ))
    });
    Ok(out)
}

fn parse_err(file: &str, msg: String) -> CombinrError {
    CombinrError::Parse {
        file: file.to_string(),
        line: 0,
        msg,
    }
}

/// Read the `XS:A` strand tag (`+`/`-`), returning [`Strand::Unknown`] when the
/// tag is absent or not a recognized strand character.
fn xs_strand(record: &bam::Record) -> Strand {
    match record.data().get(b"XS") {
        Some(Ok(Value::Character(c))) => Strand::from_char(c as char),
        _ => Strand::Unknown,
    }
}

/// Convert a CIGAR (as `(kind, len)` ops) plus the 1-based leftmost mapped
/// position into genomic exon spans (1-based inclusive). A run of
/// reference-consuming non-skip ops (`M`/`=`/`X`/`D`) is one exon; `N` (skip)
/// closes the current exon and starts the next; insertions, clips, and pads
/// never split an exon and never advance the genomic cursor (except that `D`,
/// being a reference deletion, does advance it but stays within the exon).
fn cigar_to_exons(pos_1based: i64, ops: &[(Kind, i64)]) -> Vec<Coordset> {
    let mut exons = Vec::new();
    let mut cursor = pos_1based; // next genomic base to be consumed (1-based)
    let mut block_start: Option<i64> = None;
    for &(kind, len) in ops {
        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                if block_start.is_none() {
                    block_start = Some(cursor);
                }
                cursor += len;
            }
            Kind::Skip => {
                if let Some(start) = block_start.take() {
                    exons.push(Coordset::new(start, cursor - 1));
                }
                cursor += len;
            }
            // Insertion / SoftClip / HardClip / Pad: no genomic advance, no split.
            Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
        }
    }
    if let Some(start) = block_start.take() {
        exons.push(Coordset::new(start, cursor - 1));
    }
    exons
}

#[cfg(test)]
mod tests {
    use super::*;

    // Oracle 1 — TRINITY_DN4001_c0_g2_i1, POS=331213, FLAG=16, XS:+
    // CIGAR 16S58M188N130M190N27M153N61M435N72M98N28M210N214M → 7 exons.
    #[test]
    fn cigar_to_exons_seven_exon_oracle() {
        let ops = [
            (Kind::SoftClip, 16),
            (Kind::Match, 58),
            (Kind::Skip, 188),
            (Kind::Match, 130),
            (Kind::Skip, 190),
            (Kind::Match, 27),
            (Kind::Skip, 153),
            (Kind::Match, 61),
            (Kind::Skip, 435),
            (Kind::Match, 72),
            (Kind::Skip, 98),
            (Kind::Match, 28),
            (Kind::Skip, 210),
            (Kind::Match, 214),
        ];
        let got = cigar_to_exons(331213, &ops);
        let want = vec![
            Coordset::new(331213, 331270),
            Coordset::new(331459, 331588),
            Coordset::new(331779, 331805),
            Coordset::new(331959, 332019),
            Coordset::new(332455, 332526),
            Coordset::new(332625, 332652),
            Coordset::new(332863, 333076),
        ];
        assert_eq!(got, want);
    }

    // Oracle 2 — TRINITY_DN4001_c0_g4_i1, POS=342777, FLAG=0, XS:+
    // CIGAR 4S23M1D14M1D172M112N26M1I217N72M → 3 exons.
    // The two 1D (deletions) stay inside exon 1; the 1I (insertion) stays inside
    // exon 2; only the two N ops split. Confirms D/I absorption.
    #[test]
    fn cigar_to_exons_deletion_insertion_absorbed() {
        let ops = [
            (Kind::SoftClip, 4),
            (Kind::Match, 23),
            (Kind::Deletion, 1),
            (Kind::Match, 14),
            (Kind::Deletion, 1),
            (Kind::Match, 172),
            (Kind::Skip, 112),
            (Kind::Match, 26),
            (Kind::Insertion, 1),
            (Kind::Skip, 217),
            (Kind::Match, 72),
        ];
        let got = cigar_to_exons(342777, &ops);
        let want = vec![
            Coordset::new(342777, 342987),
            Coordset::new(343100, 343125),
            Coordset::new(343343, 343414),
        ];
        assert_eq!(got, want);
    }

    // A single ungapped block: one exon, soft clips ignored.
    #[test]
    fn cigar_to_exons_single_block() {
        let ops = [(Kind::SoftClip, 5), (Kind::Match, 100), (Kind::SoftClip, 3)];
        assert_eq!(cigar_to_exons(1000, &ops), vec![Coordset::new(1000, 1099)]);
    }

    // No aligned bases (all clips) → no exons.
    #[test]
    fn cigar_to_exons_all_clips_empty() {
        let ops = [(Kind::SoftClip, 10), (Kind::HardClip, 4)];
        assert!(cigar_to_exons(500, &ops).is_empty());
    }
}
