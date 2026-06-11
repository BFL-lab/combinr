//! Assembler golden parity (Algorithm 1) against committed PASA references.
//!
//! For every vendored `pasa_cpp_sample_input*`, the original C++ `pasa` binary's
//! assembly output is committed under `tests/data/assembler/*.golden`. Three
//! combinr code paths — the raw assembler, the orientation wrapper, and the full
//! GFF3 pipeline — must each reproduce that golden assembly set (member sets +
//! merged structures). Set comparison, not line order: tied `lend` values make
//! the C++ `std::sort` and Rust's stable sort number members differently without
//! ever changing which assemblies are produced. No PASA code runs here.

mod common;

use combinr::assemble::{Assembler, assemble_cluster};
use combinr::filter::Filters;
use combinr::model::{Alignment, Strand};
use combinr::pipeline::assemble_sources;
use combinr::token::parse_tokens;
use std::collections::BTreeSet;
use std::io::Write;

fn read(p: &std::path::Path) -> String {
    std::fs::read_to_string(p).unwrap_or_else(|e| panic!("reading {}: {e}", p.display()))
}

/// The raw assembler (`assemble-tokens` path) vs golden.
#[test]
fn raw_assembler_matches_golden() {
    let fixtures = common::assembler_fixtures();
    assert!(!fixtures.is_empty(), "no assembler fixtures vendored");
    for (input, golden) in fixtures {
        let aligns = parse_tokens(&read(&input), input.to_str().unwrap()).unwrap();
        let mut asm = Assembler::new(aligns);
        asm.assemble().unwrap();
        let got = common::parse_assemblies(&asm.format_pasa_assemblies());
        let expected = common::parse_assemblies(&read(&golden));
        assert_eq!(
            got,
            expected,
            "raw assembler mismatch on {}",
            input.display()
        );
    }
}

/// The orientation wrapper (`assemble_cluster`) vs golden.
#[test]
fn orientation_wrapper_matches_golden() {
    for (input, golden) in common::assembler_fixtures() {
        let aligns = parse_tokens(&read(&input), input.to_str().unwrap()).unwrap();
        let got: BTreeSet<common::CanonAssembly> = assemble_cluster(&aligns, 20)
            .unwrap()
            .iter()
            .map(common::canon_from_cluster)
            .collect();
        let expected = common::parse_assemblies(&read(&golden));
        assert_eq!(got, expected, "wrapper mismatch on {}", input.display());
    }
}

/// The full GFF3 pipeline (`assemble_sources`) vs golden: round-trip each token
/// input through a cDNA_match GFF3 to exercise parsing + clustering too.
#[test]
fn gff3_pipeline_matches_golden() {
    for (input, golden) in common::assembler_fixtures() {
        let aligns = parse_tokens(&read(&input), input.to_str().unwrap()).unwrap();
        let gff = token_to_gff3(&aligns, "chr1");
        let mut tmp = tempfile::Builder::new().suffix(".gff3").tempfile().unwrap();
        tmp.write_all(gff.as_bytes()).unwrap();

        let asms = assemble_sources(&[tmp.path().to_path_buf()], 20, &Filters::none()).unwrap();
        let got: BTreeSet<common::CanonAssembly> =
            asms.iter().map(common::canon_from_cluster).collect();
        let expected = common::parse_assemblies(&read(&golden));
        assert_eq!(got, expected, "pipeline mismatch on {}", input.display());
    }
}

fn token_to_gff3(aligns: &[Alignment], contig: &str) -> String {
    let mut s = String::from("##gff-version 3\n");
    for a in aligns {
        let orient = match a.aligned_orient {
            Strand::Plus => '+',
            Strand::Minus => '-',
            Strand::Unknown => '.',
        };
        let mut cdna = 0i64;
        for seg in &a.segments {
            let len = seg.coords.rend - seg.coords.lend + 1;
            let (cl, cr) = (cdna + 1, cdna + len);
            cdna += len;
            s.push_str(&format!(
                "{contig}\tt\tcDNA_match\t{}\t{}\t.\t{orient}\t.\tID={};Target={} {cl} {cr}\n",
                seg.coords.lend, seg.coords.rend, a.acc, a.acc
            ));
        }
    }
    s
}
