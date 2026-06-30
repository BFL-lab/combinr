//! Parser for the C++ `pasa` token format:
//! `accession,orientation,lend1-rend1,lend2-rend2,...`
//!
//! This is exactly the format `pasa.cpp::main` reads and the format
//! `PASA_alignment_assembler.pm` writes, so it doubles as the golden-test
//! harness input and the orientation wrapper's interchange format.

use crate::error::{CombinrError, Result};
use crate::model::{Alignment, Segment, Strand};

/// Line prefixes that mark non-data lines — `//` separators and the `input:`/`assembly`/
/// `Individual` headers `pasa` emits — so its output can be fed straight back in.
const SKIP_PREFIXES: &[&str] = &["//", "input:", "assembly", "Individual"];

/// Parse all token-format alignments from `input`. Blank lines, `//` separators,
/// and `input:`/`assembly:` lines (so `pasa` output can be fed back) are
/// skipped. Like `pasa`, only lines containing a comma are treated as data.
pub fn parse_tokens(input: &str, source: &str) -> Result<Vec<Alignment>> {
    let mut out = Vec::new();
    for (i, raw) in input.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() || SKIP_PREFIXES.iter().any(|p| line.starts_with(p)) {
            continue;
        }
        if !line.contains(',') {
            continue;
        }
        out.push(parse_line(line, source, i + 1)?);
    }
    Ok(out)
}

fn parse_line(line: &str, source: &str, lineno: usize) -> Result<Alignment> {
    let err = |msg: &str| CombinrError::Parse {
        file: source.to_string(),
        line: lineno,
        msg: msg.to_string(),
    };

    let mut it = line.split(',');
    let acc = it
        .next()
        .ok_or_else(|| err("empty line"))?
        .trim()
        .to_string();
    let orient_str = it.next().ok_or_else(|| err("missing orientation"))?;
    let orient = orient_str.trim().chars().next().unwrap_or('?');

    let mut segs = Vec::new();
    for cp in it {
        let cp = cp.trim();
        if cp.is_empty() {
            continue;
        }
        let (l, r) = cp
            .split_once('-')
            .ok_or_else(|| err(&format!("bad coordinate pair {cp:?}")))?;
        let lend: i64 = l
            .trim()
            .parse()
            .map_err(|_| err(&format!("bad lend in {cp:?}")))?;
        let rend: i64 = r
            .trim()
            .parse()
            .map_err(|_| err(&format!("bad rend in {cp:?}")))?;
        segs.push(Segment::new(lend, rend));
    }

    if segs.is_empty() {
        return Err(err("no coordinate segments"));
    }

    let strand = Strand::from_char(orient);
    let mut a = Alignment::new(acc, segs, strand);
    a.spliced_orient = strand;
    Ok(a)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_a_simple_record() {
        let aligns = parse_tokens("acc1,+,100-200,300-500", "test").unwrap();
        assert_eq!(aligns.len(), 1);
        let a = &aligns[0];
        assert_eq!(a.acc, "acc1");
        assert_eq!(a.aligned_orient, Strand::Plus);
        assert_eq!(a.num_segments(), 2);
        assert_eq!(a.coords, crate::model::Coordset::new(100, 500));
    }

    #[test]
    fn skips_separators_and_output_lines() {
        let input = "//\ninput: acc9,+,1-2\nacc1,-,100-200\n\nacc2,+,50-60";
        let aligns = parse_tokens(input, "test").unwrap();
        let accs: Vec<&str> = aligns.iter().map(|a| a.acc.as_str()).collect();
        assert_eq!(accs, vec!["acc1", "acc2"]);
    }

    #[test]
    fn reports_line_number_on_bad_coords() {
        let input = "acc1,+,100-200\nacc2,+,oops";
        let e = parse_tokens(input, "f").unwrap_err();
        match e {
            CombinrError::Parse { line, .. } => assert_eq!(line, 2),
            other => panic!("unexpected error: {other:?}"),
        }
    }
}
