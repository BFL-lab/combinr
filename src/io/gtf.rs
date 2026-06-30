//! GTF alignment parser. Port of `PerlLib/GTF_alignment_utils.pm`.
//!
//! Only `exon` rows are consumed, grouped by `transcript_id`. Strand is taken
//! verbatim and becomes the transcribed orientation directly (no splice-site
//! validation); a `.` score column becomes 100. cDNA coordinates are synthesized
//! from exon lengths in transcription order.

use super::gff::{self, AttrSep};
use super::{RawSegment, build_alignment, sort_alignments_canonical};
use crate::error::Result;
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
        let Some(rec) = gff::record(raw, &source, lineno + 1)? else {
            continue;
        };
        if rec.ftype != "exon" {
            continue;
        }
        let attrs = gff::parse_attrs(rec.attrs, AttrSep::Space);
        let Some(tid) = attrs.get("transcript_id") else {
            continue;
        };

        let key = tid.clone();
        let entry = groups.entry(key.clone()).or_insert_with(|| {
            order.push(key.clone());
            Group {
                contig: rec.contig.to_string(),
                strand: rec.strand,
                gene_id: attrs.get("gene_id").cloned(),
                segs: Vec::new(),
            }
        });
        entry.segs.push(RawSegment {
            genomic: Coordset::new(rec.lend, rec.rend),
            mcoords: None,
            per_id: rec.score,
        });
    }

    let mut out = Vec::with_capacity(order.len());
    for key in order {
        let g = groups.remove(&key).expect("group present");
        out.push(build_alignment(
            key, g.contig, g.strand, g.gene_id, g.segs, &source,
        ));
    }
    sort_alignments_canonical(&mut out);
    Ok(out)
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
