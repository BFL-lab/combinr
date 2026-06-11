//! Codon translation + reverse complement. Port of the relevant parts of
//! `PerlLib/Nuc_translator.pm` (standard genetic code, table 1).

/// Reverse-complement a nucleotide sequence, handling IUPAC ambiguity codes
/// (mirrors the `tr/ACGTacgtyrkmYRKM/TGCAtgcarymkRYMK/` mapping).
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'a' => b't',
        b't' => b'a',
        b'g' => b'c',
        b'c' => b'g',
        b'Y' => b'R',
        b'R' => b'Y',
        b'K' => b'M',
        b'M' => b'K',
        b'y' => b'r',
        b'r' => b'y',
        b'k' => b'm',
        b'm' => b'k',
        b'N' => b'N',
        b'n' => b'n',
        other => other,
    }
}

/// `true` if `codon` (3 uppercase bytes) is a stop codon under the standard code.
pub fn is_stop(codon: &[u8]) -> bool {
    matches!(codon, b"TAA" | b"TAG" | b"TGA")
}

/// Translate a single codon to a one-letter amino acid (`*` = stop, `X` =
/// unknown/contains ambiguity). Standard genetic code.
pub fn translate_codon(codon: &[u8]) -> u8 {
    if codon.len() != 3 {
        return b'X';
    }
    match codon {
        b"TTT" | b"TTC" => b'F',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => b'L',
        b"ATT" | b"ATC" | b"ATA" => b'I',
        b"ATG" => b'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => b'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => b'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => b'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => b'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => b'A',
        b"TAT" | b"TAC" => b'Y',
        b"TAA" | b"TAG" | b"TGA" => b'*',
        b"CAT" | b"CAC" => b'H',
        b"CAA" | b"CAG" => b'Q',
        b"AAT" | b"AAC" => b'N',
        b"AAA" | b"AAG" => b'K',
        b"GAT" | b"GAC" => b'D',
        b"GAA" | b"GAG" => b'E',
        b"TGT" | b"TGC" => b'C',
        b"TGG" => b'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => b'R',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => b'G',
        _ => b'X',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn revcomp_with_ambiguity() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"ATGC"), b"GCAT");
        // ambiguity code R(A/G) <-> Y(C/T)
        assert_eq!(reverse_complement(b"R"), b"Y");
    }

    #[test]
    fn stops_and_translation() {
        assert!(is_stop(b"TAA") && is_stop(b"TAG") && is_stop(b"TGA"));
        assert!(!is_stop(b"ATG"));
        assert_eq!(translate_codon(b"ATG"), b'M');
        assert_eq!(translate_codon(b"TAA"), b'*');
        assert_eq!(translate_codon(b"GGG"), b'G');
        assert_eq!(translate_codon(b"NNN"), b'X');
    }
}
