//! Repeat-mask GFF3 parsing.
//!
//! Masked bases are excluded from coding/intron/intergenic scoring (EVM `repeatMask`).
//! Any feature row contributes its `[start, end]` span on its sequence (columns 1/4/5);
//! the feature type is not inspected.

use crate::error::{CombinrError, Result};
use crate::model::Coordset;
use std::collections::HashMap;
use std::path::Path;

/// Parse a repeats GFF3 into per-contig masked intervals.
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
    Ok(map)
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
}
