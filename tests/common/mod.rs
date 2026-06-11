//! Shared helpers for the golden tests.
//!
//! Golden references are committed under `tests/data/` (generated once from the
//! original PASA C++ `pasa` binary and the PASA `Alternative_splice_comparer`
//! Perl module). The tests compare combinr's output against these fixtures, so
//! **no PASA code runs at test time**. They are a frozen snapshot; the one-off
//! generation tooling is preserved in git history.

#![allow(dead_code)]

use combinr::altsplice::EventRecord;
use combinr::assemble::ClusterAssembly;
use combinr::model::Strand;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

/// One assembly canonicalized for order-independent comparison:
/// (member accessions, (orient, segment strings...)).
pub type CanonAssembly = (BTreeSet<String>, Vec<String>);

pub fn data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

/// `(input, golden)` pairs for the assembler golden (vendored
/// `pasa_cpp_sample_input*` + their C++ `pasa` reference output).
pub fn assembler_fixtures() -> Vec<(PathBuf, PathBuf)> {
    let dir = data_dir().join("assembler");
    let mut v = Vec::new();
    for e in std::fs::read_dir(&dir).unwrap() {
        let p = e.unwrap().path();
        let name = p.file_name().unwrap().to_str().unwrap().to_string();
        if name.starts_with("pasa_cpp_sample_input") && !name.ends_with(".golden") {
            v.push((dir.join(&name), dir.join(format!("{name}.golden"))));
        }
    }
    v.sort();
    v
}

/// `(input, golden)` pairs for the alt-splice golden (input GTF/GFF3 +
/// PASA-derived event reference).
pub fn altsplice_fixtures() -> Vec<(PathBuf, PathBuf)> {
    let dir = data_dir().join("altsplice");
    let mut v = Vec::new();
    for e in std::fs::read_dir(&dir).unwrap() {
        let p = e.unwrap().path();
        let name = p.file_name().unwrap().to_str().unwrap().to_string();
        if (name.ends_with(".gtf") || name.ends_with(".gff3")) && !name.contains(".events.golden") {
            let stem = name.rsplit_once('.').unwrap().0;
            v.push((dir.join(&name), dir.join(format!("{stem}.events.golden"))));
        }
    }
    v.sort();
    v
}

fn strand_str(s: Strand) -> String {
    match s {
        Strand::Plus => "+",
        Strand::Minus => "-",
        Strand::Unknown => "?",
    }
    .to_string()
}

/// Parse `pasa`/`assemble-tokens` `assembly: ...` lines into a canonical set
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

/// Canonicalize a combinr [`ClusterAssembly`] the same way [`parse_assemblies`]
/// canonicalizes a golden line.
pub fn canon_from_cluster(a: &ClusterAssembly) -> CanonAssembly {
    let members: BTreeSet<String> = a.contained_accs.iter().cloned().collect();
    let mut fields = vec![strand_str(a.orient)];
    for s in &a.structure.segments {
        fields.push(format!("{}-{}", s.coords.lend, s.coords.rend));
    }
    (members, fields)
}

/// Canonicalize combinr events to the golden event line form
/// `event_type \t coords \t isoform_a \t isoform_b` (coords lend-sorted).
pub fn canon_events(events: &[EventRecord]) -> BTreeSet<String> {
    events
        .iter()
        .map(|e| {
            let mut segs = e.coords.clone();
            segs.sort_by_key(|c| c.lend);
            let coords = segs
                .iter()
                .map(|c| format!("{}-{}", c.lend, c.rend))
                .collect::<Vec<_>>()
                .join(",");
            format!(
                "{}\t{}\t{}\t{}",
                e.kind.as_str(),
                coords,
                e.isoform_a,
                e.isoform_b
            )
        })
        .collect()
}

/// Load a golden text file as a set of non-empty trimmed lines.
pub fn load_lines(path: &Path) -> BTreeSet<String> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("reading {}: {e}", path.display()))
        .lines()
        .map(|l| l.trim_end().to_string())
        .filter(|l| !l.is_empty())
        .collect()
}

fn between<'a>(s: &'a str, start: &str, end: &str) -> Option<&'a str> {
    let i = s.find(start)? + start.len();
    let rest = &s[i..];
    let j = rest.find(end)?;
    Some(&rest[..j])
}
