//! Shared helpers for the golden-parity tests.
//!
//! Locates the reference PASA `pasa_cpp` directory (and its compiled `pasa`
//! binary) so the golden tests run wherever the reference lives:
//!   1. `$COMBINR_PASA_DIR` (the PASApipeline root), if set;
//!   2. `PASApipeline/` inside the crate (vendored);
//!   3. `../PASApipeline/` beside the crate (sibling checkout);
//!   4. `/home/matt/PASApipeline` (this workstation's location).

#![allow(dead_code)]

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

/// One assembly canonicalized for order-independent comparison:
/// (member accessions, (orient, segment strings...)).
pub type CanonAssembly = (BTreeSet<String>, Vec<String>);

/// The `pasa_cpp` directory containing a built `pasa` binary, if findable.
pub fn pasa_cpp_dir() -> Option<PathBuf> {
    let mut candidates: Vec<PathBuf> = Vec::new();
    if let Ok(d) = std::env::var("COMBINR_PASA_DIR") {
        candidates.push(PathBuf::from(d).join("pasa_cpp"));
    }
    let manifest = Path::new(env!("CARGO_MANIFEST_DIR"));
    candidates.push(manifest.join("PASApipeline/pasa_cpp"));
    candidates.push(manifest.join("../PASApipeline/pasa_cpp"));
    candidates.push(PathBuf::from("/home/matt/PASApipeline/pasa_cpp"));

    candidates.into_iter().find(|d| d.join("pasa").exists())
}

/// The compiled `pasa` binary, if findable.
pub fn pasa_binary() -> Option<PathBuf> {
    pasa_cpp_dir().map(|d| d.join("pasa"))
}

/// All `pasa_cpp_sample_input*` text inputs in `dir`, sorted.
pub fn sample_inputs(dir: &Path) -> Vec<PathBuf> {
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

/// Parse the `pasa`/`assemble-tokens` `assembly: ...` lines into a canonical set
/// (member set + structure orient/segments, ignoring the assembly title).
pub fn parse_assemblies(text: &str) -> BTreeSet<CanonAssembly> {
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
