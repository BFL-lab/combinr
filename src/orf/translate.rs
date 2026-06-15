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

/// A genetic code, reduced to what ORF reconciliation needs: which codons are
/// stops. Only stop assignments differ between the tables we support; sense
/// reassignments (e.g. CTG→Ser/Ala in tables 12/26) don't affect CDS/UTR bounds,
/// which depend solely on where translation hits the first in-frame stop. Stop
/// sets follow the NCBI translation tables (BioPython's
/// `unambiguous_dna_by_id[id].stop_codons`).
#[derive(Clone, Copy, Debug)]
pub struct GeneticCode {
    ncbi_id: u32,
    stops: &'static [&'static [u8]],
}

impl GeneticCode {
    /// Build from an NCBI translation-table id, or `Err` with a human-readable
    /// message for a table whose stop set we have not vetted.
    ///
    /// Stop-codon sets (differences from the standard code noted):
    /// - `1`  Standard ............ TAA TAG TGA
    /// - `4`  Mold/Protozoan Mito . TAA TAG      (TGA = Trp)
    /// - `6`  Ciliate Nuclear ..... TGA          (TAA TAG = Gln)
    /// - `10` Euplotid Nuclear .... TAA TAG      (TGA = Cys)
    /// - `12` Alt. Yeast Nuclear .. TAA TAG TGA  (CTG = Ser, a sense change)
    /// - `26` Pachysolen Nuclear .. TAA TAG TGA  (CTG = Ala, a sense change)
    pub fn from_ncbi_id(id: u32) -> std::result::Result<GeneticCode, String> {
        let stops: &'static [&'static [u8]] = match id {
            1 | 12 | 26 => &[b"TAA", b"TAG", b"TGA"],
            4 | 10 => &[b"TAA", b"TAG"],
            6 => &[b"TGA"],
            _ => {
                return Err(format!(
                    "unsupported NCBI genetic code {id}; supported tables: 1, 4, 6, 10, 12, 26"
                ));
            }
        };
        Ok(GeneticCode { ncbi_id: id, stops })
    }

    /// The NCBI translation-table id.
    pub fn ncbi_id(&self) -> u32 {
        self.ncbi_id
    }

    /// `true` if `codon` (3 uppercase bytes) is a stop under this code.
    pub fn is_stop(&self, codon: &[u8]) -> bool {
        self.stops.iter().any(|&s| s == codon)
    }
}

impl Default for GeneticCode {
    /// The standard code (NCBI table 1).
    fn default() -> GeneticCode {
        GeneticCode {
            ncbi_id: 1,
            stops: &[b"TAA", b"TAG", b"TGA"],
        }
    }
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
        let std = GeneticCode::default();
        assert!(std.is_stop(b"TAA") && std.is_stop(b"TAG") && std.is_stop(b"TGA"));
        assert!(!std.is_stop(b"ATG"));
        assert_eq!(translate_codon(b"ATG"), b'M');
        assert_eq!(translate_codon(b"TAA"), b'*');
        assert_eq!(translate_codon(b"GGG"), b'G');
        assert_eq!(translate_codon(b"NNN"), b'X');
    }

    #[test]
    fn genetic_code_stop_sets() {
        let c1 = GeneticCode::from_ncbi_id(1).unwrap();
        assert!(c1.is_stop(b"TAA") && c1.is_stop(b"TAG") && c1.is_stop(b"TGA"));

        // Ciliate nuclear (6): TAA/TAG are Gln; only TGA stops.
        let c6 = GeneticCode::from_ncbi_id(6).unwrap();
        assert!(c6.is_stop(b"TGA"));
        assert!(!c6.is_stop(b"TAA") && !c6.is_stop(b"TAG"));

        // Pachysolen (26): CTG→Ala is a sense change; stops match the standard code.
        let c26 = GeneticCode::from_ncbi_id(26).unwrap();
        assert!(c26.is_stop(b"TAA") && c26.is_stop(b"TAG") && c26.is_stop(b"TGA"));

        // TGA-as-sense tables (4 Mold Mito, 10 Euplotid): only TAA/TAG stop.
        let c4 = GeneticCode::from_ncbi_id(4).unwrap();
        assert!(c4.is_stop(b"TAA") && c4.is_stop(b"TAG") && !c4.is_stop(b"TGA"));

        assert!(GeneticCode::from_ncbi_id(99).is_err());
    }
}
