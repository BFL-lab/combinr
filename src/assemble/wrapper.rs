//! Orientation wrapper — port of `PerlLib/CDNA/PASA_alignment_assembler.pm`.
//!
//! The bare [`Assembler`](super::Assembler) assembles alignments of a single
//! fixed orientation. This wrapper runs it once forcing `+` and once forcing
//! `-` (so strandless single-exon "flex" alignments can join either), pools the
//! results, then greedily keeps assemblies in decreasing size as long as each
//! introduces a not-yet-covered member — yielding the minimal non-redundant set
//! that covers every input. Each kept assembly's transcribed orientation, FL
//! status, and source provenance are resolved from its members.

use super::Assembler;
use crate::error::{CombinrError, Result};
use crate::model::{Alignment, Strand};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::sync::Arc;

/// One non-redundant assembly produced from a cluster.
#[derive(Debug, Clone)]
pub struct ClusterAssembly {
    /// Merged exon structure. `acc` is the `/`-joined member accessions
    /// (PASA's "unity" accession).
    pub structure: Alignment,
    /// Member transcript accessions (ascending by source index).
    pub contained_accs: Vec<String>,
    pub num_contained: usize,
    /// Transcribed orientation resolved from members (`Unknown` only if no
    /// member had a known spliced orientation).
    pub spliced_orient: Strand,
    /// Orientation for output: the spliced orientation if known, else the
    /// majority aligned orientation of the members.
    pub orient: Strand,
    pub is_fl: bool,
    /// Union of the source files contributing members.
    pub member_provenance: BTreeSet<Arc<str>>,
}

/// Per-accession metadata captured before assembly.
struct AccMeta {
    spliced: Strand,
    aligned: Strand,
    is_fl: bool,
    provenance: Option<Arc<str>>,
}

/// A pooled raw assembly from one forced-orientation pass.
struct RawAssembly {
    structure: Alignment,
    members: Vec<String>,
}

/// Assemble one cluster of alignments into a non-redundant set.
///
/// `alignments` should already be a single cluster (same contig, overlapping).
pub fn assemble_cluster(alignments: &[Alignment], fuzzlength: i64) -> Result<Vec<ClusterAssembly>> {
    if alignments.is_empty() {
        return Ok(Vec::new());
    }

    // All members of a cluster share a contig; the merged structures inherit it
    // (merge_alignments builds contig-less structures).
    let contig = alignments[0].contig.clone();

    // Capture per-acc metadata (spliced/aligned orient, FL, provenance).
    let mut acc_meta: HashMap<String, AccMeta> = HashMap::new();
    for a in alignments {
        acc_meta.insert(
            a.acc.clone(),
            AccMeta {
                spliced: a.spliced_orient,
                aligned: a.aligned_orient,
                is_fl: a.is_fl,
                provenance: a.provenance.as_ref().map(|p| p.source_file.clone()),
            },
        );
    }
    let all_accs: HashSet<&String> = acc_meta.keys().collect();

    // Two passes: force '+' then '-'. Flex (unknown-spliced) alignments take the
    // forced orient; alignments with a known spliced orient keep it in both.
    let mut pooled: Vec<RawAssembly> = Vec::new();
    for forced in [Strand::Plus, Strand::Minus] {
        let forced_aligns = force_flexorient(alignments, forced);
        let mut asm = Assembler::new(forced_aligns);
        asm.set_fuzzlength(fuzzlength);
        asm.assemble()?;
        for (i, structure) in asm.assemblies.iter().enumerate() {
            let members: Vec<String> = asm.containment[i]
                .iter()
                .map(|&idx| asm.alignments[idx].acc.clone())
                .collect();
            pooled.push(RawAssembly {
                structure: structure.clone(),
                members,
            });
        }
    }

    // Decreasing size; stable so ties keep pass order (+ before -).
    pooled.sort_by(|a, b| b.members.len().cmp(&a.members.len()));

    // Greedy unseen-member cover.
    let mut seen: HashSet<String> = HashSet::new();
    let mut report: Vec<ClusterAssembly> = Vec::new();

    for raw in pooled {
        let mut have_unseen = false;
        let mut spliced = Strand::Unknown;
        let mut is_fl = false;
        let mut aligned_counts: HashMap<Strand, usize> = HashMap::new();
        let mut provenance: BTreeSet<Arc<str>> = BTreeSet::new();

        for acc in &raw.members {
            if !seen.contains(acc) {
                have_unseen = true;
            }
            seen.insert(acc.clone());

            let meta = &acc_meta[acc];
            if meta.spliced.is_known() {
                if spliced.is_known() && spliced != meta.spliced {
                    return Err(CombinrError::ConflictingSplicedOrientation {
                        have: spliced.to_char(),
                        acc: acc.clone(),
                        found: meta.spliced.to_char(),
                    });
                }
                spliced = meta.spliced;
            }
            if meta.is_fl {
                is_fl = true;
            }
            *aligned_counts.entry(meta.aligned).or_insert(0) += 1;
            if let Some(p) = &meta.provenance {
                provenance.insert(p.clone());
            }
        }

        let orient = if spliced.is_known() {
            spliced
        } else {
            majority_orient(&aligned_counts)
        };

        if have_unseen {
            let mut structure = raw.structure;
            structure.acc = raw.members.join("/");
            structure.contig = contig.clone();
            structure.aligned_orient = orient;
            structure.spliced_orient = spliced;
            structure.is_fl = is_fl;
            report.push(ClusterAssembly {
                structure,
                num_contained: raw.members.len(),
                contained_accs: raw.members,
                spliced_orient: spliced,
                orient,
                is_fl,
                member_provenance: provenance,
            });
        }

        if seen.len() == all_accs.len() {
            break; // all members covered
        }
    }

    if seen.len() != all_accs.len() {
        return Err(CombinrError::UncoveredAlignments(all_accs.len()));
    }

    Ok(report)
}

/// Port of `force_flexorient`: clone the alignments, setting each one's working
/// (aligned) orientation to its known spliced orientation, or to `forced` if it
/// has none.
fn force_flexorient(alignments: &[Alignment], forced: Strand) -> Vec<Alignment> {
    alignments
        .iter()
        .map(|a| {
            let mut c = a.clone();
            c.aligned_orient = if a.spliced_orient.is_known() {
                a.spliced_orient
            } else {
                forced
            };
            c
        })
        .collect()
}

/// Highest-count aligned orientation; deterministic tie-break (Plus < Minus <
/// Unknown, via `Strand`'s derived `Ord`). Returns `Unknown` only if there are no
/// members.
fn majority_orient(counts: &HashMap<Strand, usize>) -> Strand {
    counts
        .iter()
        .max_by(|a, b| a.1.cmp(b.1).then_with(|| b.0.cmp(a.0)))
        .map(|(s, _)| *s)
        .unwrap_or(Strand::Unknown)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Coordset, Segment};

    fn spliced(acc: &str, orient: Strand, segs: &[(i64, i64)]) -> Alignment {
        let mut a = Alignment::new(
            acc,
            segs.iter().map(|&(l, r)| Segment::new(l, r)).collect(),
            orient,
        );
        a.spliced_orient = orient;
        a
    }

    #[test]
    fn single_orientation_cluster_assembles_like_raw() {
        // Two compatible '-' alignments sharing an intron should yield one
        // assembly covering both, even across the +/- two-pass wrapper.
        let aligns = vec![
            spliced("a", Strand::Minus, &[(100, 200), (300, 400)]),
            spliced("b", Strand::Minus, &[(120, 200), (300, 450)]),
        ];
        let asms = assemble_cluster(&aligns, 20).unwrap();
        assert_eq!(asms.len(), 1);
        let a = &asms[0];
        assert_eq!(a.num_contained, 2);
        let mut members = a.contained_accs.clone();
        members.sort();
        assert_eq!(members, vec!["a", "b"]);
        assert_eq!(a.orient, Strand::Minus);
        assert_eq!(
            a.structure.segments[0].coords,
            Coordset::new(100, 200),
            "termini extended to the union"
        );
        assert_eq!(a.structure.segments[1].coords, Coordset::new(300, 450));
    }

    #[test]
    fn flex_singleton_joins_a_known_orientation_assembly() {
        // A strandless single-exon read overlapping a '+' spliced assembly's
        // first exon should be covered (joining the '+' pass) and not surface
        // as its own assembly.
        let mut flex = Alignment::new("flex", vec![Segment::new(130, 180)], Strand::Unknown);
        flex.spliced_orient = Strand::Unknown;
        let aligns = vec![spliced("p", Strand::Plus, &[(100, 200), (300, 400)]), flex];
        let asms = assemble_cluster(&aligns, 20).unwrap();
        // Every input acc is covered.
        let covered: HashSet<&str> = asms
            .iter()
            .flat_map(|a| a.contained_accs.iter().map(|s| s.as_str()))
            .collect();
        assert!(covered.contains("p"));
        assert!(covered.contains("flex"));
        // The assembly carrying both is '+' (the only known spliced orient).
        let joint = asms
            .iter()
            .find(|a| a.contained_accs.iter().any(|m| m == "flex"))
            .unwrap();
        assert_eq!(joint.orient, Strand::Plus);
    }

    #[test]
    fn fl_status_propagates_from_any_member() {
        let mut a = spliced("a", Strand::Plus, &[(100, 200), (300, 400)]);
        let mut b = spliced("b", Strand::Plus, &[(120, 200), (300, 420)]);
        a.is_fl = false;
        b.is_fl = true;
        let asms = assemble_cluster(&[a, b], 20).unwrap();
        assert_eq!(asms.len(), 1);
        assert!(asms[0].is_fl, "FL propagates if any member is FL");
    }
}
