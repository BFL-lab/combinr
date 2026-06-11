//! End-to-end parity (M3): convert each `pasa_cpp_sample_input*` into a
//! cDNA_match GFF3, run the real pipeline (`load → cluster → assemble`), and
//! assert the resulting assembly set still matches the golden C++ `pasa` binary.
//! This exercises GFF3 parsing, single-linkage clustering, and the orientation
//! wrapper together — each sample is one contig/one orientation, so it forms a
//! single cluster that must reproduce the reference assemblies.

use combinr::model::{Alignment, Strand};
use combinr::pipeline::assemble_sources;
use combinr::token::parse_tokens;
use std::collections::BTreeSet;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

type CanonAssembly = (BTreeSet<String>, Vec<String>);

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

fn parse_golden(text: &str) -> BTreeSet<CanonAssembly> {
    let mut out = BTreeSet::new();
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("assembly:") else {
            continue;
        };
        let (Some(members), Some(structure)) = (
            between(rest, "contains alignments: [", "]"),
            between(rest, "with structure [", "]"),
        ) else {
            continue;
        };
        let member_set: BTreeSet<String> = members.split(',').map(String::from).collect();
        let fields: Vec<String> = structure.split(',').skip(1).map(String::from).collect();
        out.insert((member_set, fields));
    }
    out
}

fn between<'a>(s: &'a str, start: &str, end: &str) -> Option<&'a str> {
    let i = s.find(start)? + start.len();
    let rest = &s[i..];
    let j = rest.find(end)?;
    Some(&rest[..j])
}

fn pasa_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("PASApipeline/pasa_cpp")
}

fn sample_inputs(dir: &Path) -> Vec<PathBuf> {
    let mut v: Vec<PathBuf> = std::fs::read_dir(dir)
        .unwrap()
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .map(|n| n.starts_with("pasa_cpp_sample_input"))
                .unwrap_or(false)
                && !matches!(
                    p.extension().and_then(|e| e.to_str()),
                    Some("o") | Some("output")
                )
        })
        .collect();
    v.sort();
    v
}

#[test]
fn pipeline_matches_golden_pasa_via_gff3() {
    let dir = pasa_dir();
    let pasa = dir.join("pasa");
    if !pasa.exists() {
        eprintln!(
            "SKIP: golden `pasa` binary not built — run `make` in {}",
            dir.display()
        );
        return;
    }

    for input in sample_inputs(&dir) {
        let text = std::fs::read_to_string(&input).unwrap();
        let aligns = parse_tokens(&text, input.to_str().unwrap()).unwrap();

        // Write a cDNA_match GFF3 to a temp file and run the real pipeline.
        let gff = token_to_gff3(&aligns, "chr1");
        let mut tmp = tempfile::Builder::new().suffix(".gff3").tempfile().unwrap();
        tmp.write_all(gff.as_bytes()).unwrap();
        let path = tmp.path().to_path_buf();

        let asms = assemble_sources(&[path], 20, &combinr::filter::Filters::none()).unwrap();
        let got: BTreeSet<CanonAssembly> = asms
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
        let expected = parse_golden(&String::from_utf8_lossy(&golden.stdout));

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
