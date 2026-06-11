//! GTF alignment parser. Port of `PerlLib/GTF_alignment_utils.pm`.
//!
//! Only `exon` rows are consumed, grouped by `transcript_id`. Strand is taken
//! verbatim and becomes the transcribed orientation directly (no splice-site
//! validation); a `.` score column becomes 100. cDNA coordinates are synthesized
//! from exon lengths in transcription order.

use super::{RawSegment, build_alignment};
use crate::error::{CombinrError, Result};
use crate::model::{Alignment, Coordset, Strand};
use std::collections::HashMap;
use std::sync::Arc;

struct Group {
    contig: String,
    strand: Strand,
    gene_id: Option<String>,
    segs: Vec<RawSegment>,
}

pub fn parse(text: &str, source: &str) -> Result<Vec<Alignment>> {
    let source: Arc<str> = Arc::from(source);
    let mut order: Vec<String> = Vec::new();
    let mut groups: HashMap<String, Group> = HashMap::new();

    for (lineno, raw) in text.lines().enumerate() {
        let line = raw.trim_end();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }
        if cols[2] != "exon" {
            continue;
        }
        let attrs = parse_attrs(cols[8]);
        let Some(tid) = attrs.get("transcript_id") else {
            continue;
        };

        let err = |msg: &str| CombinrError::Parse {
            file: source.to_string(),
            line: lineno + 1,
            msg: msg.to_string(),
        };
        let lend: i64 = cols[3].parse().map_err(|_| err("bad start coordinate"))?;
        let rend: i64 = cols[4].parse().map_err(|_| err("bad end coordinate"))?;
        let per_id = parse_score(cols[5]);
        let strand = Strand::from_char(cols[6].chars().next().unwrap_or('.'));

        let key = tid.clone();
        let entry = groups.entry(key.clone()).or_insert_with(|| {
            order.push(key.clone());
            Group {
                contig: cols[0].to_string(),
                strand,
                gene_id: attrs.get("gene_id").cloned(),
                segs: Vec::new(),
            }
        });
        entry.segs.push(RawSegment {
            genomic: Coordset::new(lend, rend),
            mcoords: None,
            per_id,
        });
    }

    let mut out = Vec::with_capacity(order.len());
    for key in order {
        let g = groups.remove(&key).expect("group present");
        out.push(build_alignment(
            key, g.contig, g.strand, g.gene_id, g.segs, &source,
        ));
    }
    out.sort_by(|a, b| {
        (a.contig.as_str(), a.coords.lend, a.acc.as_str()).cmp(&(
            b.contig.as_str(),
            b.coords.lend,
            b.acc.as_str(),
        ))
    });
    Ok(out)
}

/// Parse `key "value"; key "value";` GTF attributes.
fn parse_attrs(s: &str) -> HashMap<String, String> {
    let mut m = HashMap::new();
    for part in s.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        if let Some((k, v)) = part.split_once(char::is_whitespace) {
            m.insert(k.trim().to_string(), v.trim().trim_matches('"').to_string());
        }
    }
    m
}

fn parse_score(s: &str) -> Option<f64> {
    if s == "." {
        Some(100.0)
    } else {
        s.parse().ok()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_stringtie_like_gtf() {
        let gtf = "\
chr2\tStringTie\ttranscript\t1000\t2000\t.\t+\t.\tgene_id \"G.1\"; transcript_id \"G.1.1\";
chr2\tStringTie\texon\t1000\t1200\t.\t+\t.\tgene_id \"G.1\"; transcript_id \"G.1.1\";
chr2\tStringTie\texon\t1500\t2000\t.\t+\t.\tgene_id \"G.1\"; transcript_id \"G.1.1\";";
        let aligns = parse(gtf, "st.gtf").unwrap();
        assert_eq!(aligns.len(), 1, "transcript row ignored, exons grouped");
        let a = &aligns[0];
        assert_eq!(a.acc, "G.1.1");
        assert_eq!(a.gene_id.as_deref(), Some("G.1"));
        assert_eq!(a.contig, "chr2");
        assert_eq!(a.spliced_orient, Strand::Plus);
        assert_eq!(a.num_segments(), 2);
        assert_eq!(a.coords, Coordset::new(1000, 2000));
        assert!(a.segments.iter().all(|s| s.mcoords.is_some()));
    }

    #[test]
    fn minus_strand_sets_orientation() {
        let gtf = "\
chrX\ts\texon\t10\t20\t.\t-\t.\tgene_id \"g\"; transcript_id \"t\";
chrX\ts\texon\t40\t60\t.\t-\t.\tgene_id \"g\"; transcript_id \"t\";";
        let aligns = parse(gtf, "m.gtf").unwrap();
        assert_eq!(aligns[0].spliced_orient, Strand::Minus);
    }
}
