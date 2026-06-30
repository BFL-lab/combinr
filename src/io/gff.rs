//! Shared GFF/GTF record tokenizer used by every parser in the crate (the GFF3 and GTF
//! alignment readers, the consensus evidence reader, the gene-prediction CDS reader, and
//! the format sniffer). It single-sources the row skeleton — comment/blank/short-line
//! skipping, the tab split, the 1-based coordinate parse and its [`CombinrError::Parse`],
//! the strand and score columns — so the per-format readers carry only their distinct
//! grouping-key logic. Values are borrowed from the input line (no allocation here).

use crate::error::{CombinrError, Result};
use crate::model::Strand;
use std::collections::HashMap;

/// Which delimiter separates an attribute key from its value: GFF3 uses `=`, GTF uses
/// whitespace.
#[derive(Clone, Copy, Debug)]
pub(crate) enum AttrSep {
    /// GFF3 `key=value`.
    Eq,
    /// GTF `key "value"`.
    Space,
}

/// One tab-delimited GFF/GTF record, borrowed from the source line.
pub(crate) struct GffRecord<'a> {
    pub contig: &'a str,
    /// Column 2 (the source / program), which keys the consensus weights file.
    pub source: &'a str,
    /// Column 3 (the feature type).
    pub ftype: &'a str,
    pub lend: i64,
    pub rend: i64,
    /// Column 6 score → percent identity (`.` becomes 100, the PASA assumption).
    pub score: Option<f64>,
    pub strand: Strand,
    /// Column 9, the raw attribute string (parse with [`parse_attrs`]).
    pub attrs: &'a str,
}

/// Split a record line into its tab columns, or `None` for a blank/comment line or one
/// with fewer than the 9 mandatory GFF columns (silently skipped, as all readers do).
pub(crate) fn columns(line: &str) -> Option<Vec<&str>> {
    let line = line.trim_end();
    if line.is_empty() || line.starts_with('#') {
        return None;
    }
    let cols: Vec<&str> = line.split('\t').collect();
    (cols.len() >= 9).then_some(cols)
}

/// Tokenize one GFF/GTF record. Returns `Ok(None)` for a blank/comment/short line and
/// `Err` (with `file` + 1-based `line_number`) for malformed coordinates.
pub(crate) fn record<'a>(
    line: &'a str,
    file: &str,
    line_number: usize,
) -> Result<Option<GffRecord<'a>>> {
    let Some(cols) = columns(line) else {
        return Ok(None);
    };
    let bad = |msg: &str| CombinrError::Parse {
        file: file.to_string(),
        line: line_number,
        msg: msg.to_string(),
    };
    let lend: i64 = cols[3].parse().map_err(|_| bad("bad start coordinate"))?;
    let rend: i64 = cols[4].parse().map_err(|_| bad("bad end coordinate"))?;
    Ok(Some(GffRecord {
        contig: cols[0],
        source: cols[1],
        ftype: cols[2],
        lend,
        rend,
        score: parse_score(cols[5]),
        strand: Strand::from_char(cols[6].chars().next().unwrap_or('.')),
        attrs: cols[8],
    }))
}

/// Parse a `key<sep>value;...` attribute column into a map. Values are unquoted.
pub(crate) fn parse_attrs(s: &str, sep: AttrSep) -> HashMap<String, String> {
    let mut m = HashMap::new();
    for part in s.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let pair = match sep {
            AttrSep::Eq => part.split_once('='),
            AttrSep::Space => part.split_once(char::is_whitespace),
        };
        if let Some((k, v)) = pair {
            m.insert(k.trim().to_string(), v.trim().trim_matches('"').to_string());
        }
    }
    m
}

/// Score column → percent identity; `.` becomes 100 (the PASA assumption).
pub(crate) fn parse_score(s: &str) -> Option<f64> {
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
    fn skips_blank_comment_and_short_lines() {
        assert!(columns("").is_none());
        assert!(columns("# a comment").is_none());
        assert!(columns("chr1\tsrc\texon\t1\t10").is_none()); // < 9 columns
        assert!(columns("chr1\ts\texon\t1\t10\t.\t+\t.\tID=x").is_some());
    }

    #[test]
    fn record_parses_columns_and_errors_on_bad_coords() {
        let r = record(
            "chr1\tgmap\tcDNA_match\t100\t200\t98.5\t+\t.\tID=m1",
            "f",
            1,
        )
        .unwrap()
        .unwrap();
        assert_eq!(
            (r.contig, r.source, r.ftype),
            ("chr1", "gmap", "cDNA_match")
        );
        assert_eq!((r.lend, r.rend), (100, 200));
        assert_eq!(r.score, Some(98.5));
        assert_eq!(r.strand, Strand::Plus);
        assert!(record("chr1\ts\texon\tNOPE\t200\t.\t+\t.\tID=x", "f", 7).is_err());
    }

    #[test]
    fn attrs_honour_the_separator() {
        let g = parse_attrs("ID=m1;Target=cdnaA 1 101", AttrSep::Eq);
        assert_eq!(g.get("ID").map(String::as_str), Some("m1"));
        assert_eq!(g.get("Target").map(String::as_str), Some("cdnaA 1 101"));
        let t = parse_attrs("gene_id \"G.1\"; transcript_id \"G.1.1\";", AttrSep::Space);
        assert_eq!(t.get("gene_id").map(String::as_str), Some("G.1"));
        assert_eq!(t.get("transcript_id").map(String::as_str), Some("G.1.1"));
    }

    #[test]
    fn dot_score_becomes_100() {
        assert_eq!(parse_score("."), Some(100.0));
        assert_eq!(parse_score("90.5"), Some(90.5));
        assert_eq!(parse_score("junk"), None);
    }
}
