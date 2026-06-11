//! Input/output: multi-source GTF/GFF3 parsing (with provenance) and GFF3
//! output. Parsers are hand-rolled for byte-faithful parity with PASA's
//! idiosyncratic readers (1-based inclusive coords, `.` percent-id → 100,
//! Target-token accessions) and to keep the binary self-contained.

pub mod bam;
pub mod fasta;
pub mod gff3;
pub mod gtf;
pub mod out_model;
pub mod writer_events;
pub mod writer_gff3;
pub mod writer_gtf;

use crate::error::{CombinrError, Result};
use crate::model::{Alignment, Coordset, Provenance, Segment, Strand};
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// Supported input formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Gff3,
    Gtf,
    /// Binary BAM. Detected and dispatched before any text read (see
    /// [`is_bam`]); never produced by [`detect_format`], which only sees text.
    Bam,
}

/// Load and concatenate alignments from several source files, stamping each with
/// its file-basename provenance. Format is auto-detected per file.
pub fn load_sources(paths: &[PathBuf]) -> Result<Vec<Alignment>> {
    let mut out = Vec::new();
    for path in paths {
        let source = basename(path);
        // BAM is binary: dispatch it before any attempt to read the file as text.
        if is_bam(path)? {
            let mut aligns = bam::parse(path, &source)?;
            out.append(&mut aligns);
            continue;
        }
        let text = std::fs::read_to_string(path).map_err(|e| CombinrError::Parse {
            file: path.display().to_string(),
            line: 0,
            msg: format!("cannot read: {e}"),
        })?;
        let fmt = detect_format(path, &text);
        let mut aligns = match fmt {
            Format::Gff3 => gff3::parse(&text, &source)?,
            Format::Gtf => gtf::parse(&text, &source)?,
            Format::Bam => unreachable!("BAM is handled before read_to_string"),
        };
        out.append(&mut aligns);
    }
    Ok(out)
}

fn basename(path: &Path) -> String {
    path.file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("<input>")
        .to_string()
}

/// `true` if `path` is a BAM file: a `.bam` extension (case-insensitive) backed
/// by the BGZF/gzip magic prefix `1f 8b`. A `.bam` file lacking the magic is a
/// parse error (so a mislabeled text file fails clearly instead of being read as
/// text). CRAM and plain-text SAM are not recognized here.
fn is_bam(path: &Path) -> Result<bool> {
    let has_bam_ext = path
        .extension()
        .and_then(|e| e.to_str())
        .is_some_and(|s| s.eq_ignore_ascii_case("bam"));
    if !has_bam_ext {
        return Ok(false);
    }
    use std::io::Read;
    let mut magic = [0u8; 2];
    match std::fs::File::open(path)?.read_exact(&mut magic) {
        Ok(()) if magic == [0x1f, 0x8b] => Ok(true),
        Ok(()) => Err(CombinrError::Parse {
            file: path.display().to_string(),
            line: 0,
            msg: "has a .bam extension but is not BGZF/BAM (missing 1f 8b magic)".to_string(),
        }),
        Err(e) => Err(CombinrError::Io(e)),
    }
}

/// Detect format by extension, falling back to a content sniff.
pub fn detect_format(path: &Path, text: &str) -> Format {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|s| s.to_ascii_lowercase());
    match ext.as_deref() {
        Some("gtf") => return Format::Gtf,
        Some("gff") | Some("gff3") => return Format::Gff3,
        _ => {}
    }
    sniff_format(text)
}

/// Content sniff: GTF attributes look like `key "value";`; GFF3 like `key=value`.
fn sniff_format(text: &str) -> Format {
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }
        let attrs = cols[8];
        if attrs.contains('=')
            && (attrs.contains("ID=") || attrs.contains("Parent=") || attrs.contains("Target="))
        {
            return Format::Gff3;
        }
        if attrs.contains("transcript_id ") || attrs.contains("gene_id ") {
            return Format::Gtf;
        }
    }
    Format::Gff3
}

/// A raw exon segment accumulated during parsing: genomic span, optional
/// cDNA/transcript span, optional percent identity.
pub(crate) struct RawSegment {
    pub genomic: Coordset,
    pub mcoords: Option<Coordset>,
    pub per_id: Option<f64>,
}

/// Build a finished [`Alignment`] from accumulated raw segments. Synthesizes
/// transcript (mcoords) coordinates when absent (gene-model GFF3 / GTF) and
/// computes the length-weighted average percent identity.
pub(crate) fn build_alignment(
    acc: String,
    contig: String,
    strand: Strand,
    gene_id: Option<String>,
    mut raw: Vec<RawSegment>,
    source: &Arc<str>,
) -> Alignment {
    // genomic order
    raw.sort_by_key(|s| (s.genomic.lend, s.genomic.rend));

    // Synthesize mcoords if any are missing: walk in transcription order
    // (reverse on minus strand) accumulating cDNA length.
    if raw.iter().any(|s| s.mcoords.is_none()) {
        let mut order: Vec<usize> = (0..raw.len()).collect();
        if strand == Strand::Minus {
            order.reverse();
        }
        let mut cum = 0i64;
        for &i in &order {
            let len = raw[i].genomic.len();
            let m = Coordset::new(cum + 1, cum + len);
            cum += len;
            raw[i].mcoords = Some(m);
        }
    }

    let segments: Vec<Segment> = raw
        .iter()
        .map(|s| Segment {
            coords: s.genomic,
            seg_type: crate::model::SegType::Single, // reset by refine()
            left_splice: false,
            right_splice: false,
            mcoords: s.mcoords,
            per_id: s.per_id,
        })
        .collect();

    // length-weighted average percent identity, when present on all segments.
    let per_id = weighted_avg_per_id(&raw);

    let mut a = Alignment::new(acc, segments, strand);
    a.contig = contig;
    a.spliced_orient = strand;
    a.gene_id = gene_id;
    a.per_id = per_id;
    a.provenance = Some(Provenance {
        source_file: source.clone(),
    });
    a
}

fn weighted_avg_per_id(raw: &[RawSegment]) -> Option<f64> {
    let mut num = 0.0;
    let mut den = 0.0;
    for s in raw {
        let pid = s.per_id?;
        let w = s.genomic.len() as f64;
        num += pid * w;
        den += w;
    }
    if den > 0.0 { Some(num / den) } else { None }
}
