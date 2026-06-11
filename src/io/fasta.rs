//! Minimal in-memory FASTA reader with 1-based inclusive subsequence access.
//! Used only by the ORF/UTR step.

use crate::error::{CombinrError, Result};
use std::collections::HashMap;
use std::path::Path;

/// A loaded set of contig sequences (uppercased).
pub struct Fasta {
    seqs: HashMap<String, Vec<u8>>,
}

impl Fasta {
    /// Load all records from a FASTA file. The contig key is the first
    /// whitespace token of each header line.
    pub fn load(path: &Path) -> Result<Fasta> {
        let text = std::fs::read_to_string(path).map_err(|e| CombinrError::Parse {
            file: path.display().to_string(),
            line: 0,
            msg: format!("cannot read genome FASTA: {e}"),
        })?;
        let mut seqs: HashMap<String, Vec<u8>> = HashMap::new();
        let mut current: Option<String> = None;
        for line in text.lines() {
            if let Some(header) = line.strip_prefix('>') {
                let name = header.split_whitespace().next().unwrap_or("").to_string();
                current = Some(name.clone());
                seqs.entry(name).or_default();
            } else if let Some(name) = &current {
                let seq = seqs.get_mut(name).unwrap();
                for b in line.trim().bytes() {
                    if !b.is_ascii_whitespace() {
                        seq.push(b.to_ascii_uppercase());
                    }
                }
            }
        }
        Ok(Fasta { seqs })
    }

    /// 1-based inclusive subsequence `[lend, rend]` on `contig`, or `None` if the
    /// contig is absent or the range is out of bounds.
    pub fn subseq(&self, contig: &str, lend: i64, rend: i64) -> Option<&[u8]> {
        if lend < 1 || rend < lend {
            return None;
        }
        let seq = self.seqs.get(contig)?;
        let (lo, hi) = ((lend - 1) as usize, rend as usize);
        if hi > seq.len() {
            return None;
        }
        Some(&seq[lo..hi])
    }

    pub fn contig_len(&self, contig: &str) -> Option<usize> {
        self.seqs.get(contig).map(|s| s.len())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn loads_and_fetches_subseq_1based() {
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        writeln!(tmp, ">chr1 description here").unwrap();
        writeln!(tmp, "ACGTACGT").unwrap();
        writeln!(tmp, "TTTTGGGG").unwrap();
        let fa = Fasta::load(tmp.path()).unwrap();
        assert_eq!(fa.contig_len("chr1"), Some(16));
        // 1-based inclusive: bases 1..4 = ACGT
        assert_eq!(fa.subseq("chr1", 1, 4).unwrap(), b"ACGT");
        // spanning the line break
        assert_eq!(fa.subseq("chr1", 7, 10).unwrap(), b"GTTT");
        assert!(fa.subseq("chr1", 10, 100).is_none(), "out of bounds");
        assert!(fa.subseq("chrX", 1, 2).is_none(), "missing contig");
    }
}
