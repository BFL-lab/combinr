//! Promote assemblies to isoforms (de-duplicating identical structures) and
//! group them into gene loci by transitive exon overlap on the same strand.

use super::{Isoform, Locus};
use crate::assemble::ClusterAssembly;
use crate::model::{Coordset, Strand};
use std::collections::{BTreeMap, BTreeSet, HashMap};

/// Build isoforms from assemblies, collapsing any with identical
/// `(strand, exon vector)` and merging their provenance + members.
pub fn build_isoforms(assemblies: &[ClusterAssembly]) -> Vec<Isoform> {
    let mut by_key: HashMap<(Strand, Vec<Coordset>), usize> = HashMap::new();
    let mut isoforms: Vec<Isoform> = Vec::new();

    for asm in assemblies {
        let exons: Vec<Coordset> = asm.structure.segments.iter().map(|s| s.coords).collect();
        let key = (asm.orient, exons.clone());
        if let Some(&idx) = by_key.get(&key) {
            // merge into the existing identical isoform
            let iso = &mut isoforms[idx];
            for acc in &asm.contained_accs {
                if !iso.contained_accs.contains(acc) {
                    iso.contained_accs.push(acc.clone());
                }
            }
            for src in &asm.member_provenance {
                iso.source_set.insert(src.clone());
            }
            iso.is_fl |= asm.is_fl;
            continue;
        }
        let mut source_set = BTreeSet::new();
        for src in &asm.member_provenance {
            source_set.insert(src.clone());
        }
        by_key.insert(key, isoforms.len());
        isoforms.push(Isoform {
            id: String::new(), // assigned during locus grouping
            contig: asm.structure.contig.clone(),
            strand: asm.orient,
            exons,
            contained_accs: asm.contained_accs.clone(),
            source_set,
            is_fl: asm.is_fl,
        });
    }
    isoforms
}

/// Group isoforms into loci by `(contig, strand)` single-linkage span overlap,
/// assigning stable `locus_N` / `locus_N.isoM` ids in place.
///
/// `min_overlap_frac` is PASA's `--stringent_alignment_overlap`: two isoforms only
/// share a locus when their genomic spans overlap by at least that percent of the
/// **shorter** span. At `0.0` any overlap groups them (the historical behavior).
pub fn group_into_loci(isoforms: &mut [Isoform], min_overlap_frac: f64) -> Vec<Locus> {
    // bucket indices by (contig, strand). `Strand` derives `Ord`, so the key order is
    // deterministic (Plus < Minus < Unknown) without a manual u8 table.
    let mut buckets: BTreeMap<(String, Strand), Vec<usize>> = BTreeMap::new();
    for (i, iso) in isoforms.iter().enumerate() {
        buckets
            .entry((iso.contig.clone(), iso.strand))
            .or_default()
            .push(i);
    }

    // ordered list of runs (each becomes a locus), sorted deterministically.
    let mut runs: Vec<(String, Strand, i64, Vec<usize>)> = Vec::new();
    for ((contig, strand), idxs) in buckets {
        let spans: Vec<Coordset> = idxs
            .iter()
            .map(|&i| Coordset::new(isoforms[i].exons[0].lend, span_rend(&isoforms[i])))
            .collect();
        for group in crate::cluster::single_linkage_groups(&spans, min_overlap_frac) {
            let members: Vec<usize> = group.into_iter().map(|l| idxs[l]).collect();
            let run_lend = members
                .iter()
                .map(|&i| isoforms[i].exons[0].lend)
                .min()
                .expect("non-empty locus group");
            runs.push((contig.clone(), strand, run_lend, members));
        }
    }

    // assign deterministic locus + isoform ids.
    runs.sort_by(|a, b| (a.0.as_str(), a.1, a.2).cmp(&(b.0.as_str(), b.1, b.2)));

    let mut loci = Vec::with_capacity(runs.len());
    for (n, (contig, strand, _lend, idxs)) in runs.into_iter().enumerate() {
        let locus_id = format!("locus_{}", n + 1);
        for (m, &idx) in idxs.iter().enumerate() {
            isoforms[idx].id = format!("{locus_id}.iso{}", m + 1);
        }
        loci.push(Locus {
            id: locus_id,
            contig,
            strand,
            isoform_indices: idxs,
        });
    }
    loci
}

fn span_rend(iso: &Isoform) -> i64 {
    iso.exons.last().map(|e| e.rend).unwrap_or(i64::MIN)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assemble::assemble_cluster;
    use crate::model::{Alignment, Segment};

    fn spliced(acc: &str, contig: &str, orient: Strand, segs: &[(i64, i64)]) -> Alignment {
        let mut a = Alignment::new(
            acc,
            segs.iter().map(|&(l, r)| Segment::new(l, r)).collect(),
            orient,
        );
        a.contig = contig.to_string();
        a.spliced_orient = orient;
        a
    }

    #[test]
    fn overlapping_isoforms_share_a_locus() {
        // Two distinct isoforms (different middle exon) that overlap.
        let asm1 = assemble_cluster(
            &[spliced(
                "a",
                "chr1",
                Strand::Plus,
                &[(100, 200), (300, 400), (700, 800)],
            )],
            20,
        )
        .unwrap();
        let asm2 = assemble_cluster(
            &[spliced(
                "b",
                "chr1",
                Strand::Plus,
                &[(100, 200), (500, 600), (700, 800)],
            )],
            20,
        )
        .unwrap();
        let all: Vec<_> = asm1.into_iter().chain(asm2).collect();

        let mut isoforms = build_isoforms(&all);
        let loci = group_into_loci(&mut isoforms, 0.0);
        assert_eq!(loci.len(), 1, "overlapping isoforms group together");
        assert_eq!(loci[0].isoform_indices.len(), 2);
        assert!(isoforms.iter().all(|i| i.id.starts_with("locus_1.iso")));
    }

    #[test]
    fn identical_structures_are_deduplicated() {
        let a = assemble_cluster(
            &[spliced(
                "a",
                "chr1",
                Strand::Plus,
                &[(100, 200), (300, 400)],
            )],
            20,
        )
        .unwrap();
        let b = assemble_cluster(
            &[spliced(
                "b",
                "chr1",
                Strand::Plus,
                &[(100, 200), (300, 400)],
            )],
            20,
        )
        .unwrap();
        let all: Vec<_> = a.into_iter().chain(b).collect();
        let isoforms = build_isoforms(&all);
        assert_eq!(isoforms.len(), 1, "identical exon structures collapse");
        // both members preserved
        let mut accs = isoforms[0].contained_accs.clone();
        accs.sort();
        assert_eq!(accs, vec!["a", "b"]);
    }

    #[test]
    fn distant_isoforms_are_separate_loci() {
        let a = assemble_cluster(
            &[spliced(
                "a",
                "chr1",
                Strand::Plus,
                &[(100, 200), (300, 400)],
            )],
            20,
        )
        .unwrap();
        let b = assemble_cluster(
            &[spliced(
                "b",
                "chr1",
                Strand::Plus,
                &[(5000, 5100), (5300, 5400)],
            )],
            20,
        )
        .unwrap();
        let all: Vec<_> = a.into_iter().chain(b).collect();
        let mut isoforms = build_isoforms(&all);
        let loci = group_into_loci(&mut isoforms, 0.0);
        assert_eq!(loci.len(), 2);
    }

    #[test]
    fn stringent_overlap_splits_collinear_isoforms() {
        // Two distinct isoforms whose genomic spans (100..1100 and 1000..2000, each
        // len 1001) overlap only at the tips (101 bp ≈ 10% of the shorter). At the
        // default 0.0 they share a locus; at 30% they split into two genes.
        let a = assemble_cluster(
            &[spliced(
                "a",
                "chr1",
                Strand::Plus,
                &[(100, 200), (300, 1100)],
            )],
            20,
        )
        .unwrap();
        let b = assemble_cluster(
            &[spliced(
                "b",
                "chr1",
                Strand::Plus,
                &[(1000, 1100), (1300, 2000)],
            )],
            20,
        )
        .unwrap();
        let all: Vec<_> = a.into_iter().chain(b).collect();

        let mut iso0 = build_isoforms(&all);
        assert_eq!(group_into_loci(&mut iso0, 0.0).len(), 1);

        let mut iso30 = build_isoforms(&all);
        assert_eq!(group_into_loci(&mut iso30, 30.0).len(), 2);
    }

    #[test]
    fn stringent_overlap_keeps_contained_isoform_in_locus() {
        // A short isoform fully inside a long one (overlap = 100% of the shorter span)
        // still shares the locus even at a high threshold.
        let long = assemble_cluster(
            &[spliced(
                "long",
                "chr1",
                Strand::Plus,
                &[(100, 200), (2000, 2100)],
            )],
            20,
        )
        .unwrap();
        let short = assemble_cluster(
            &[spliced(
                "short",
                "chr1",
                Strand::Plus,
                &[(500, 600), (700, 800)],
            )],
            20,
        )
        .unwrap();
        let all: Vec<_> = long.into_iter().chain(short).collect();
        let mut isoforms = build_isoforms(&all);
        assert_eq!(group_into_loci(&mut isoforms, 90.0).len(), 1);
    }
}
