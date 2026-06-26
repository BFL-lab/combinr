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
/// stops. Only stop assignments differ between the tables; sense reassignments
/// (e.g. CTG→Ser/Ala in tables 12/26) don't affect CDS/UTR bounds, which depend
/// solely on where translation hits the first in-frame stop. Stop sets follow the
/// NCBI translation tables (as declared in BioPython's `CodonTable`,
/// `register_ncbi_table(... stop_codons=[...])`).
#[derive(Clone, Copy, Debug)]
pub struct GeneticCode {
    ncbi_id: u32,
    stops: &'static [&'static [u8]],
}

impl GeneticCode {
    /// Build from an NCBI translation-table id, or `Err` with a human-readable
    /// message for an id NCBI has never assigned.
    ///
    /// Only stop assignments matter here — sense reassignments (e.g. CTG→Ser/Ala
    /// in tables 12/26) don't move CDS/UTR bounds, which depend solely on the first
    /// in-frame stop — so tables are grouped by their stop-codon set, taken from the
    /// NCBI translation tables (BioPython `register_ncbi_table(... stop_codons=[...])`):
    ///
    /// - TAA TAG TGA (standard) ...... 1, 11, 12, 26, 28
    /// - TAA TAG ..................... 3, 4, 5, 9, 10, 13, 21, 24, 25, 31
    /// - TGA ......................... 6, 27, 29, 30   (ciliate/karyorelict nuclear; TAA TAG = Gln/Tyr/Glu)
    /// - TAA TAG AGA AGG ............. 2   (vertebrate mito)
    /// - TAG ......................... 14, 33   (alt. flatworm mito; Cephalodiscidae "UAA-Tyr" mito)
    /// - TAA TGA ..................... 15, 16, 32   (TAG = Gln/Leu/Trp)
    /// - TCA TAA TGA ................. 22  (Scenedesmus obliquus mito; TAG = Leu)
    /// - TTA TAA TAG TGA ............. 23  (Thraustochytrium mito)
    ///
    /// Tables 27/28/31 have *dual-function* codons (a codon serves as both sense and
    /// stop depending on context); we use NCBI's declared stop set, i.e. such a codon
    /// terminates translation. That is the conservative choice for CDS bounds —
    /// translation halts at the first occurrence rather than reading through a real
    /// stop — and still honours each table's reassignments (e.g. table 27 reads
    /// TAA/TAG as Gln and stops only at TGA). Ids 7, 8, and 17–20 were deleted or
    /// never assigned by NCBI and are rejected.
    pub fn from_ncbi_id(id: u32) -> std::result::Result<GeneticCode, String> {
        let stops: &'static [&'static [u8]] = match id {
            1 | 11 | 12 | 26 | 28 => &[b"TAA", b"TAG", b"TGA"],
            3 | 4 | 5 | 9 | 10 | 13 | 21 | 24 | 25 | 31 => &[b"TAA", b"TAG"],
            6 | 27 | 29 | 30 => &[b"TGA"],
            2 => &[b"TAA", b"TAG", b"AGA", b"AGG"],
            14 | 33 => &[b"TAG"],
            15 | 16 | 32 => &[b"TAA", b"TGA"],
            22 => &[b"TCA", b"TAA", b"TGA"],
            23 => &[b"TTA", b"TAA", b"TAG", b"TGA"],
            _ => {
                return Err(format!(
                    "unknown NCBI genetic code {id}; valid tables are 1-6, 9-16, \
                     and 21-33 (ids 7, 8, and 17-20 are unassigned)"
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

        // Bacterial/plastid (11): stops match the standard code (only starts differ).
        let c11 = GeneticCode::from_ncbi_id(11).unwrap();
        assert!(c11.is_stop(b"TAA") && c11.is_stop(b"TAG") && c11.is_stop(b"TGA"));

        // Vertebrate mito (2): AGA/AGG become stops alongside TAA/TAG.
        let c2 = GeneticCode::from_ncbi_id(2).unwrap();
        assert!(c2.is_stop(b"AGA") && c2.is_stop(b"AGG") && c2.is_stop(b"TAA"));

        // Chlorophycean mito (16): TAG → Leu, so only TAA/TGA stop.
        let c16 = GeneticCode::from_ncbi_id(16).unwrap();
        assert!(c16.is_stop(b"TAA") && c16.is_stop(b"TGA") && !c16.is_stop(b"TAG"));
    }

    #[test]
    fn context_dependent_and_late_tables() {
        // Karyorelict (27): TAA/TAG → Gln, only TGA terminates.
        let c27 = GeneticCode::from_ncbi_id(27).unwrap();
        assert!(c27.is_stop(b"TGA") && !c27.is_stop(b"TAA") && !c27.is_stop(b"TAG"));

        // Condylostoma (28): TAA/TAG/TGA are all stop-capable (dual function) →
        // conservative stop set is the standard one.
        let c28 = GeneticCode::from_ncbi_id(28).unwrap();
        assert!(c28.is_stop(b"TAA") && c28.is_stop(b"TAG") && c28.is_stop(b"TGA"));

        // Mesodinium (29) & Peritrich (30): TAA/TAG reassigned (Tyr/Glu), only TGA stops.
        for id in [29, 30] {
            let c = GeneticCode::from_ncbi_id(id).unwrap();
            assert!(c.is_stop(b"TGA") && !c.is_stop(b"TAA") && !c.is_stop(b"TAG"));
        }

        // Blastocrithidia (31): TGA → Trp; TAA/TAG terminate.
        let c31 = GeneticCode::from_ncbi_id(31).unwrap();
        assert!(c31.is_stop(b"TAA") && c31.is_stop(b"TAG") && !c31.is_stop(b"TGA"));

        // Balanophoraceae plastid (32): TAG → Trp; TAA/TGA stop.
        let c32 = GeneticCode::from_ncbi_id(32).unwrap();
        assert!(c32.is_stop(b"TAA") && c32.is_stop(b"TGA") && !c32.is_stop(b"TAG"));

        // Cephalodiscidae mito "UAA-Tyr" (33): TAA → Tyr, only TAG stops.
        let c33 = GeneticCode::from_ncbi_id(33).unwrap();
        assert!(c33.is_stop(b"TAG") && !c33.is_stop(b"TAA") && !c33.is_stop(b"TGA"));

        // Deleted / never-assigned ids and out-of-range ids are rejected.
        for id in [0, 7, 8, 17, 18, 19, 20, 34, 99] {
            assert!(
                GeneticCode::from_ncbi_id(id).is_err(),
                "id {id} should be rejected"
            );
        }
    }
}
