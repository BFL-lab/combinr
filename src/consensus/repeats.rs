//! Repeat-mask GFF3 parsing.
//!
//! Masked bases are excluded from coding/intron/intergenic scoring (EVM `repeatMask`).
//! Any feature row contributes its `[start, end]` span on its sequence (columns 1/4/5);
//! the feature type is not inspected.
//!
//! Invariant: each contig's interval list is sorted by `lend` and non-overlapping
//! (overlapping and adjacent intervals are merged; the union of masked bases is unchanged),
//! so [`overlapping`] can binary-search the intervals intersecting a region.

use crate::error::{CombinrError, Result};
use crate::model::Coordset;
use std::collections::HashMap;
use std::path::Path;

/// Parse a repeats GFF3 into per-contig masked intervals: sorted, non-overlapping, per
/// contig (overlapping/adjacent rows are merged).
pub fn parse_repeats(path: &Path) -> Result<HashMap<String, Vec<Coordset>>> {
    let text = std::fs::read_to_string(path).map_err(|e| CombinrError::Parse {
        file: path.display().to_string(),
        line: 0,
        msg: format!("cannot read repeats GFF3: {e}"),
    })?;
    let mut map: HashMap<String, Vec<Coordset>> = HashMap::new();
    for (lineno, raw) in text.lines().enumerate() {
        let line = raw.trim_end();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 5 {
            continue;
        }
        let err = |m: &str| CombinrError::Parse {
            file: path.display().to_string(),
            line: lineno + 1,
            msg: m.to_string(),
        };
        let lend: i64 = cols[3].parse().map_err(|_| err("bad repeat start"))?;
        let rend: i64 = cols[4].parse().map_err(|_| err("bad repeat end"))?;
        map.entry(cols[0].to_string())
            .or_default()
            .push(Coordset::new(lend, rend));
    }
    for ivs in map.values_mut() {
        ivs.sort_unstable_by_key(|iv| (iv.lend, iv.rend));
        let mut merged: Vec<Coordset> = Vec::with_capacity(ivs.len());
        for iv in ivs.drain(..) {
            match merged.last_mut() {
                Some(last) if iv.lend <= last.rend.saturating_add(1) => {
                    last.rend = last.rend.max(iv.rend);
                }
                _ => merged.push(iv),
            }
        }
        *ivs = merged;
    }
    Ok(map)
}

/// The sub-slice of `sorted` intersecting `span` (inclusive). Requires the
/// [`parse_repeats`] invariant: sorted by `lend`, non-overlapping.
pub fn overlapping(sorted: &[Coordset], span: Coordset) -> &[Coordset] {
    let start = sorted.partition_point(|iv| iv.rend < span.lend);
    let end = start + sorted[start..].partition_point(|iv| iv.lend <= span.rend);
    &sorted[start..end]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_repeat_intervals_by_contig() {
        let gff = "\
chr1\tRepeatMasker\tmatch\t100\t200\t.\t+\t.\tName=ALU
chr1\tRepeatMasker\tmatch\t500\t600\t.\t-\t.\tName=LINE
chr2\tRepeatMasker\tmatch\t10\t20\t.\t+\t.\tName=ALU";
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut tmp, gff.as_bytes()).unwrap();
        let m = parse_repeats(tmp.path()).unwrap();
        assert_eq!(
            m["chr1"],
            vec![Coordset::new(100, 200), Coordset::new(500, 600)]
        );
        assert_eq!(m["chr2"], vec![Coordset::new(10, 20)]);
    }

    #[test]
    fn sorts_and_merges_overlapping_and_adjacent() {
        let gff = "\
chr1\tRM\tmatch\t500\t600\t.\t+\t.\t.
chr1\tRM\tmatch\t100\t200\t.\t+\t.\t.
chr1\tRM\tmatch\t150\t180\t.\t+\t.\t.
chr1\tRM\tmatch\t201\t250\t.\t+\t.\t.
chr1\tRM\tmatch\t590\t700\t.\t+\t.\t.
chr1\tRM\tmatch\t252\t260\t.\t+\t.\t.";
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut tmp, gff.as_bytes()).unwrap();
        let m = parse_repeats(tmp.path()).unwrap();
        assert_eq!(
            m["chr1"],
            vec![
                Coordset::new(100, 250),
                Coordset::new(252, 260),
                Coordset::new(500, 700)
            ]
        );
    }

    #[test]
    fn overlapping_selects_intersecting_intervals() {
        let ivs = [
            Coordset::new(10, 20),
            Coordset::new(30, 40),
            Coordset::new(50, 60),
        ];
        let q = |a, b| overlapping(&ivs, Coordset::new(a, b)).to_vec();
        // no overlap (gap between intervals), before all, after all
        assert!(q(21, 29).is_empty());
        assert!(q(1, 9).is_empty());
        assert!(q(61, 100).is_empty());
        // touching boundaries are inclusive
        assert_eq!(
            q(20, 30),
            vec![Coordset::new(10, 20), Coordset::new(30, 40)]
        );
        assert_eq!(q(60, 60), vec![Coordset::new(50, 60)]);
        // an interval spanning the whole query
        assert_eq!(q(33, 37), vec![Coordset::new(30, 40)]);
        // query spanning everything
        assert_eq!(q(1, 100), ivs.to_vec());
        assert!(overlapping(&[], Coordset::new(1, 10)).is_empty());
    }
}
