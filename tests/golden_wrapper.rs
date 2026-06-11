//! Golden parity for the orientation wrapper (M2), via the library API.
//!
//! The sample inputs are single-orientation clusters, so the wrapper's two
//! forced passes are identical and its greedy cover must reproduce exactly the
//! reference `pasa` assembly set (members + merged structure). This is the
//! degenerate case that proves the wrapper does not regress Algorithm 1; mixed-
//! orientation/flex behavior is covered by unit tests in `assemble::wrapper`.

use combinr::assemble::assemble_cluster;
use combinr::model::Strand;
use combinr::token::parse_tokens;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::Command;

type CanonAssembly = (BTreeSet<String>, Vec<String>);

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
fn wrapper_matches_golden_pasa_on_all_samples() {
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
        let asms = assemble_cluster(&aligns, 20).unwrap();

        let got: BTreeSet<CanonAssembly> = asms
            .iter()
            .map(|a| {
                let members: BTreeSet<String> = a.contained_accs.iter().cloned().collect();
                let mut fields = vec![match a.orient {
                    Strand::Plus => "+".to_string(),
                    Strand::Minus => "-".to_string(),
                    Strand::Unknown => "?".to_string(),
                }];
                for s in &a.structure.segments {
                    fields.push(format!("{}-{}", s.coords.lend, s.coords.rend));
                }
                (members, fields)
            })
            .collect();

        let golden = Command::new(&pasa).arg(&input).output().expect("run pasa");
        let expected = parse_golden(&String::from_utf8_lossy(&golden.stdout));

        assert_eq!(
            got,
            expected,
            "wrapper assembly set mismatch on {}\n  only in combinr: {:?}\n  only in pasa: {:?}",
            input.display(),
            got.difference(&expected).collect::<Vec<_>>(),
            expected.difference(&got).collect::<Vec<_>>(),
        );
    }
}
