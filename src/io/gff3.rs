//! GFF3 alignment parser. Port of `PerlLib/GFF3_alignment_utils.pm`, extended to
//! also read plain gene-model GFF3 (`exon` rows with `Parent=`).
//!
//! Two row shapes are accepted, decided per row:
//! - **cDNA_match**: any feature row carrying a `Target=acc lend rend` attribute.
//!   Segments group by their `ID`; the cDNA span comes from `Target`.
//! - **gene model**: `exon` rows carrying `Parent=`. Segments group by `Parent`
//!   (the mRNA id); cDNA coordinates are synthesized from exon lengths.
//!
//! Coordinates are 1-based inclusive; a `.` in the score column becomes 100
//! (matching PASA's assumption). Strand is taken verbatim — no splice-site
//! validation (see `avoid-canonical-splice-bias`).

use super::{RawSegment, build_alignment};
use crate::error::{CombinrError, Result};
use crate::model::{Alignment, Coordset, Strand};
use std::collections::HashMap;
use std::sync::Arc;

struct Group {
    contig: String,
    strand: Strand,
    segs: Vec<RawSegment>,
}

pub fn parse(text: &str, source: &str) -> Result<Vec<Alignment>> {
    let source: Arc<str> = Arc::from(source);

    // mRNA/transcript id -> gene (Parent), for gene-model gene_id resolution.
    let mut mrna_to_gene: HashMap<String, String> = HashMap::new();
    // ordered group keys + data.
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
        let attrs = parse_attrs(cols[8]);
        let ftype = cols[2];

        // Track mRNA -> gene for later gene_id assignment.
        if matches!(ftype, "mRNA" | "transcript") {
            if let (Some(id), Some(parent)) = (attrs.get("ID"), attrs.get("Parent")) {
                mrna_to_gene.insert(id.clone(), parent.clone());
            }
            continue;
        }

        let err = |msg: &str| CombinrError::Parse {
            file: source.to_string(),
            line: lineno + 1,
            msg: msg.to_string(),
        };

        let (key, mcoords): (String, Option<Coordset>) = if let Some(target) = attrs.get("Target") {
            // cDNA_match form: group by ID (unique alignment id), cDNA from Target.
            let key = attrs
                .get("ID")
                .cloned()
                .unwrap_or_else(|| target_acc(target));
            (key, parse_target_coords(target))
        } else if ftype == "exon" {
            // gene-model form: group by Parent (the mRNA id).
            let Some(parent) = attrs.get("Parent") else {
                continue;
            };
            (parent.clone(), None)
        } else {
            continue; // gene/CDS/region/etc.
        };

        let lend: i64 = cols[3].parse().map_err(|_| err("bad start coordinate"))?;
        let rend: i64 = cols[4].parse().map_err(|_| err("bad end coordinate"))?;
        let per_id = parse_score(cols[5]);
        let strand = Strand::from_char(cols[6].chars().next().unwrap_or('.'));

        let entry = groups.entry(key.clone()).or_insert_with(|| {
            order.push(key.clone());
            Group {
                contig: cols[0].to_string(),
                strand,
                segs: Vec::new(),
            }
        });
        entry.segs.push(RawSegment {
            genomic: Coordset::new(lend, rend),
            mcoords,
            per_id,
        });
    }

    // Resolve gene ids for gene-model groups (key == mRNA id).
    let mut out = Vec::with_capacity(order.len());
    for key in order {
        let g = groups.remove(&key).expect("group present");
        let gene_id = mrna_to_gene.get(&key).cloned();
        out.push(build_alignment(
            key, g.contig, g.strand, gene_id, g.segs, &source,
        ));
    }
    // deterministic order
    out.sort_by(|a, b| {
        (a.contig.as_str(), a.coords.lend, a.acc.as_str()).cmp(&(
            b.contig.as_str(),
            b.coords.lend,
            b.acc.as_str(),
        ))
    });
    Ok(out)
}

/// Parse `key=value;key=value` GFF3 attributes.
pub(crate) fn parse_attrs(s: &str) -> HashMap<String, String> {
    let mut m = HashMap::new();
    for part in s.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        if let Some((k, v)) = part.split_once('=') {
            m.insert(k.trim().to_string(), v.trim().trim_matches('"').to_string());
        }
    }
    m
}

/// First whitespace token of a `Target` value (the cDNA accession).
fn target_acc(target: &str) -> String {
    target.split_whitespace().next().unwrap_or("").to_string()
}

/// Extract the two cDNA coordinates from `Target=acc lend rend [strand]`.
fn parse_target_coords(target: &str) -> Option<Coordset> {
    let mut it = target.split_whitespace();
    let _acc = it.next()?;
    let l: i64 = it.next()?.parse().ok()?;
    let r: i64 = it.next()?.parse().ok()?;
    Some(Coordset::new(l, r))
}

/// Score column → percent identity; `.` becomes 100 (PASA assumption).
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
    fn parses_cdna_match_form() {
        let gff = "\
chr1\tgmap\tcDNA_match\t100\t200\t98.5\t+\t.\tID=m1;Target=cdnaA 1 101
chr1\tgmap\tcDNA_match\t300\t400\t99.0\t+\t.\tID=m1;Target=cdnaA 102 202";
        let aligns = parse(gff, "src.gff3").unwrap();
        assert_eq!(aligns.len(), 1);
        let a = &aligns[0];
        assert_eq!(a.acc, "m1");
        assert_eq!(a.contig, "chr1");
        assert_eq!(a.spliced_orient, Strand::Plus);
        assert_eq!(a.num_segments(), 2);
        assert_eq!(a.coords, Coordset::new(100, 400));
        assert_eq!(
            a.provenance.as_ref().unwrap().source_file.as_ref(),
            "src.gff3"
        );
        // length-weighted avg per id
        let pid = a.per_id.unwrap();
        assert!((pid - 98.75).abs() < 1e-6, "got {pid}");
    }

    #[test]
    fn parses_gene_model_form_and_resolves_gene_id() {
        let gff = "\
chr1\tsrc\tgene\t100\t400\t.\t-\t.\tID=g1
chr1\tsrc\tmRNA\t100\t400\t.\t-\t.\tID=t1;Parent=g1
chr1\tsrc\texon\t100\t200\t.\t-\t.\tID=e1;Parent=t1
chr1\tsrc\texon\t300\t400\t.\t-\t.\tID=e2;Parent=t1";
        let aligns = parse(gff, "g.gff3").unwrap();
        assert_eq!(aligns.len(), 1);
        let a = &aligns[0];
        assert_eq!(a.acc, "t1");
        assert_eq!(a.gene_id.as_deref(), Some("g1"));
        assert_eq!(a.spliced_orient, Strand::Minus);
        assert_eq!(a.num_segments(), 2);
        // mcoords synthesized in transcription order (minus strand).
        assert!(a.segments.iter().all(|s| s.mcoords.is_some()));
    }

    #[test]
    fn dot_score_becomes_100() {
        let gff = "chr1\ts\tcDNA_match\t1\t10\t.\t+\t.\tID=x;Target=c 1 10";
        let aligns = parse(gff, "s").unwrap();
        assert_eq!(aligns[0].per_id, Some(100.0));
    }
}
