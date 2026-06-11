//! End-to-end parity (M3): convert each `pasa_cpp_sample_input*` into a
//! cDNA_match GFF3, run the real pipeline (`load → cluster → assemble`), and
//! assert the resulting assembly set still matches the reference `pasa` binary.
//! Exercises GFF3 parsing, single-linkage clustering, and the orientation
//! wrapper together — each sample is one contig/one orientation, one cluster.

mod common;

use combinr::filter::Filters;
use combinr::model::{Alignment, Strand};
use combinr::pipeline::assemble_sources;
use combinr::token::parse_tokens;
use std::collections::BTreeSet;
use std::io::Write;
use std::process::Command;

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

#[test]
fn pipeline_matches_golden_pasa_via_gff3() {
    let Some(pasa) = common::pasa_binary() else {
        eprintln!("SKIP: reference `pasa` binary not found");
        return;
    };
    let dir = pasa.parent().unwrap().to_path_buf();

    for input in common::sample_inputs(&dir) {
        let text = std::fs::read_to_string(&input).unwrap();
        let aligns = parse_tokens(&text, input.to_str().unwrap()).unwrap();

        let gff = token_to_gff3(&aligns, "chr1");
        let mut tmp = tempfile::Builder::new().suffix(".gff3").tempfile().unwrap();
        tmp.write_all(gff.as_bytes()).unwrap();

        let asms = assemble_sources(&[tmp.path().to_path_buf()], 20, &Filters::none()).unwrap();
        let got: BTreeSet<common::CanonAssembly> = asms
            .iter()
            .map(|a| {
                let members: BTreeSet<String> = a.contained_accs.iter().cloned().collect();
                let mut fields = vec![match a.orient {
                    Strand::Plus => "+".to_string(),
                    Strand::Minus => "-".to_string(),
                    Strand::Unknown => "?".to_string(),
                }];
                for seg in &a.structure.segments {
                    fields.push(format!("{}-{}", seg.coords.lend, seg.coords.rend));
                }
                (members, fields)
            })
            .collect();

        let golden = Command::new(&pasa).arg(&input).output().expect("run pasa");
        let expected = common::parse_assemblies(&String::from_utf8_lossy(&golden.stdout));

        assert_eq!(
            got,
            expected,
            "pipeline assembly set mismatch on {}\n  only in combinr: {:?}\n  only in pasa: {:?}",
            input.display(),
            got.difference(&expected).collect::<Vec<_>>(),
            expected.difference(&got).collect::<Vec<_>>(),
        );
    }
}
