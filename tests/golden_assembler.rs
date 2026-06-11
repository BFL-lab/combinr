//! Golden parity test for the assembler core (M1).
//!
//! For every `pasa_cpp_sample_input*`, run both the compiled C++ `pasa` binary
//! and `combinr assemble-tokens`, then assert the produced **assembly sets** are
//! identical. Sets — not line order — because tied `lend` values make the C++
//! `std::sort` (unstable) and Rust's stable sort assign different indices, which
//! only reorders how members are listed, never which assemblies are produced.
//!
//! If the `pasa` binary has not been built (`cd PASApipeline/pasa_cpp && make`),
//! the test skips with a notice rather than failing.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::Command;

/// One assembly canonicalized for order-independent comparison:
/// (set of member accessions, (orient, segment strings...)).
type CanonAssembly = (BTreeSet<String>, Vec<String>);

fn parse_assemblies(text: &str) -> BTreeSet<CanonAssembly> {
    let mut out = BTreeSet::new();
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("assembly:") else {
            continue;
        };
        // ... contains alignments: [A,B,C] with structure [title,orient,seg,seg]
        let members = between(rest, "contains alignments: [", "]");
        let structure = between(rest, "with structure [", "]");
        let (Some(members), Some(structure)) = (members, structure) else {
            continue;
        };
        let member_set: BTreeSet<String> = members.split(',').map(|s| s.to_string()).collect();
        // Drop the title (field 0); keep orient + segments.
        let struct_fields: Vec<String> = structure
            .split(',')
            .skip(1)
            .map(|s| s.to_string())
            .collect();
        out.insert((member_set, struct_fields));
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
        })
        .filter(|p| {
            // only the raw text inputs, not object/output files
            !matches!(
                p.extension().and_then(|e| e.to_str()),
                Some("o") | Some("output")
            )
        })
        .collect();
    v.sort();
    v
}

#[test]
fn assembler_matches_golden_pasa_on_all_samples() {
    let dir = pasa_dir();
    let pasa = dir.join("pasa");
    if !pasa.exists() {
        eprintln!(
            "SKIP: golden `pasa` binary not built at {} — run `cd PASApipeline/pasa_cpp && make`",
            pasa.display()
        );
        return;
    }
    let combinr = env!("CARGO_BIN_EXE_combinr");

    let inputs = sample_inputs(&dir);
    assert!(
        !inputs.is_empty(),
        "no sample inputs found in {}",
        dir.display()
    );

    let mut checked = 0;
    for input in &inputs {
        let golden = Command::new(&pasa).arg(input).output().expect("run pasa");
        let expected = parse_assemblies(&String::from_utf8_lossy(&golden.stdout));

        let mine = Command::new(combinr)
            .args(["assemble-tokens", input.to_str().unwrap()])
            .output()
            .expect("run combinr");
        assert!(
            mine.status.success(),
            "combinr failed on {}: {}",
            input.display(),
            String::from_utf8_lossy(&mine.stderr)
        );
        let got = parse_assemblies(&String::from_utf8_lossy(&mine.stdout));

        assert_eq!(
            got,
            expected,
            "assembly set mismatch on {}\n  only in combinr: {:?}\n  only in pasa: {:?}",
            input.display(),
            got.difference(&expected).collect::<Vec<_>>(),
            expected.difference(&got).collect::<Vec<_>>(),
        );
        checked += 1;
    }
    eprintln!("golden parity: {checked} sample inputs matched");
}
