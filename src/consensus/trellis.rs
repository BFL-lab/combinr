//! The gene-structure dynamic program: score candidate exons, build the trellis, and
//! traverse the best path into consensus genes. Forward-strand only for now (the
//! minus-strand reverse-complement + transpose lands in M3).
//!
//! Port of `score_exons` (~2199), `build_trellis` (~2265), `are_compatible_exons`
//! (~2646), `score_boundary_condition` (~3248) and `traverse_path` (~2765).
//!
//! Intron keying deviates from EVM deliberately: an intron is keyed on its *actual gap
//! span* `(A.rend+1, B.lend-1)` in both `add_introns` (M1) and here, rather than EVM's
//! canonical `lend-2` dinucleotide convention — self-consistent for the non-canonical
//! model.

use crate::consensus::exon::ExonCandidate;
use crate::consensus::grammar::{Linkage, classify, frame_transition_ok, is_gene_boundary};
use crate::consensus::vectors::{IntronScores, RegionVectors};
use crate::consensus::weights::Weights;
use crate::model::Strand;
use crate::orf::GeneticCode;

/// A consensus gene: ordered exon indices (into the region's exon slice), overall
/// orientation, and the gene's score (Σ exon base scores + Σ intron scores).
#[derive(Debug, Clone, PartialEq)]
pub struct ConsensusGene {
    pub exon_indices: Vec<usize>,
    pub orient: Strand,
    pub score: f64,
}

/// `score_exons` (~2230): exon base score = `Σ max(0, coding[i])` over the exon, plus,
/// for each TRANSCRIPT/OTHER_PREDICTION evidence on the exon, `weight × unmasked_len`.
pub fn score_exon(exon: &ExonCandidate, vectors: &RegionVectors, weights: &Weights) -> f64 {
    let (lend, rend) = (exon.coords.lend, exon.coords.rend);
    let mut score = vectors.coding_sum(lend, rend);
    let mut specific = 0.0;
    for (_acc, ev_type) in &exon.evidence {
        if let Some((class, weight)) = weights.lookup(ev_type)
            && class.exon_specific()
        {
            specific += weight;
        }
    }
    if specific != 0.0 {
        score += specific * vectors.unmasked_len(lend, rend) as f64;
    }
    score
}

/// Base scores for every exon, in the same order as `exons`.
pub fn score_all_exons(
    exons: &[ExonCandidate],
    vectors: &RegionVectors,
    weights: &Weights,
) -> Vec<f64> {
    exons
        .iter()
        .map(|e| score_exon(e, vectors, weights))
        .collect()
}

/// `are_compatible_exons` (~2646): returns the join score (intron score for phased
/// links, intergenic score for unphased), or a negative value when incompatible.
/// `a` precedes `b` (sorted by 5' end); forward strand only for now.
fn are_compatible_exons(
    a: &ExonCandidate,
    b: &ExonCandidate,
    introns: &IntronScores,
    vectors: &RegionVectors,
    code: &GeneticCode,
    intergenic_adjust: f64,
) -> f64 {
    let link = classify((a.exon_type, a.orient), (b.exon_type, b.orient));
    if link == Linkage::Incompatible {
        return -1.0;
    }

    // no overlap allowed
    if a.coords.lend <= b.coords.rend && a.coords.rend >= b.coords.lend {
        return -1.0;
    }

    match link {
        Linkage::Phased => {
            // intron = the gap between the two exons (forward genomic coords)
            let intron = (a.coords.rend + 1, b.coords.lend - 1);
            let Some(intron_score) = introns.score(intron) else {
                return -1.0; // no evidence-supported intron here
            };
            if !frame_transition_ok(a.end_frame, b.start_frame) {
                return -1.0;
            }
            if creates_stop_across_junction(a, b, code) {
                return -1.0;
            }
            intron_score
        }
        Linkage::Intergenic => {
            vectors.intergenic_score(a.coords.rend + 1, b.coords.lend - 1, intergenic_adjust)
        }
        Linkage::Incompatible => -1.0,
    }
}

/// The across-junction stop test (~2719): concatenate `a`'s last two bases with `b`'s
/// first two, then, depending on `a`'s end frame, check whether the codon spanning the
/// splice is a stop under `code`. Forward strand.
fn creates_stop_across_junction(a: &ExonCandidate, b: &ExonCandidate, code: &GeneticCode) -> bool {
    let junction = [
        a.right_seq_boundary[0],
        a.right_seq_boundary[1],
        b.left_seq_boundary[0],
        b.left_seq_boundary[1],
    ];
    let end_frame = match a.end_frame % 3 {
        0 => 3,
        m => m,
    };
    let codon: Option<&[u8]> = match end_frame {
        1 => Some(&junction[1..4]),
        2 => Some(&junction[0..3]),
        _ => None, // frame 3: the codon does not span the junction
    };
    codon.is_some_and(|c| code.is_stop(c))
}

/// Build the trellis over the exons within `range`, traverse the best path, and split
/// it into consensus genes. `base_scores` is parallel to `exons` (from
/// [`score_all_exons`]). Indices in the returned genes refer back into `exons`.
#[allow(clippy::too_many_arguments)]
pub fn run_trellis(
    exons: &[ExonCandidate],
    base_scores: &[f64],
    range: (i64, i64),
    introns: &IntronScores,
    vectors: &RegionVectors,
    code: &GeneticCode,
    intergenic_adjust: f64,
    max_prev: usize,
) -> Vec<ConsensusGene> {
    // exons within range, sorted by 5' end (lend on the forward strand)
    let mut idxs: Vec<usize> = (0..exons.len())
        .filter(|&i| exons[i].coords.lend >= range.0 && exons[i].coords.rend <= range.1)
        .collect();
    if idxs.is_empty() {
        return Vec::new();
    }
    idxs.sort_by_key(|&i| (exons[i].coords.lend, exons[i].coords.rend));

    let k = idxs.len();
    let last = k + 1; // right-bound node index; node 0 = left bound
    let num_nodes = k + 2;

    let is_bound = |n: usize| n == 0 || n == last;
    let exon_at = |n: usize| &exons[idxs[n - 1]]; // valid for 1..=k
    let base_of = |n: usize| {
        if is_bound(n) {
            0.0
        } else {
            base_scores[idxs[n - 1]]
        }
    };

    // join score from node `j` (precedes) to node `i`
    let join = |j: usize, i: usize| -> f64 {
        let j_bound = is_bound(j);
        let i_bound = is_bound(i);
        if j_bound && i_bound {
            return -1.0;
        }
        if j_bound {
            // left bound -> exon i: intergenic over [range.0, exon.lend - 1]
            return vectors.intergenic_score(
                range.0,
                exon_at(i).coords.lend - 1,
                intergenic_adjust,
            );
        }
        if i_bound {
            // exon j -> right bound: intergenic over [exon.rend + 1, range.1]
            return vectors.intergenic_score(
                exon_at(j).coords.rend + 1,
                range.1,
                intergenic_adjust,
            );
        }
        are_compatible_exons(
            exon_at(j),
            exon_at(i),
            introns,
            vectors,
            code,
            intergenic_adjust,
        )
    };

    let mut sum: Vec<f64> = (0..num_nodes).map(base_of).collect();
    let mut link: Vec<Option<usize>> = vec![None; num_nodes];

    let mut highest = f64::NEG_INFINITY;
    let mut highest_node = 0usize;

    for i in 1..num_nodes {
        let base_i = base_of(i);
        let mut best = sum[i]; // initialized to base_i (standalone)
        let mut compare = 0usize;
        let mut found = false;
        let mut jj = i;
        while jj > 0 {
            jj -= 1; // candidate predecessor j = i-1, i-2, ..., 0
            // EVM look-back guard: stop once we've compared the cap AND found a link
            if compare >= max_prev && found {
                break;
            }
            compare += 1;
            let j_score = join(jj, i);
            if j_score >= 0.0 {
                found = true;
                let score = base_i + sum[jj] + j_score;
                if score > best {
                    best = score;
                    link[i] = Some(jj);
                }
            }
        }
        sum[i] = best;
        if best >= highest {
            highest = best;
            highest_node = i;
        }
    }

    // backtrace from the highest-scoring node
    let mut path = Vec::new();
    let mut cur = Some(highest_node);
    while let Some(n) = cur {
        path.push(n);
        cur = link[n];
    }
    path.reverse();

    // split into genes at gene-boundary exons; drop bound sentinels
    let mut genes = Vec::new();
    let mut current: Vec<usize> = Vec::new();
    for &n in &path {
        if is_bound(n) {
            continue;
        }
        let orig = idxs[n - 1];
        current.push(orig);
        let e = &exons[orig];
        if is_gene_boundary(e.exon_type, e.orient) {
            genes.push(finish_gene(
                std::mem::take(&mut current),
                exons,
                base_scores,
                introns,
            ));
        }
    }
    if !current.is_empty() {
        genes.push(finish_gene(current, exons, base_scores, introns));
    }
    genes
}

/// Assemble a [`ConsensusGene`] from its ordered exon indices, scoring it as
/// `Σ exon base scores + Σ intron scores` (EVM `EVM_prediction::_init`).
fn finish_gene(
    exon_indices: Vec<usize>,
    exons: &[ExonCandidate],
    base_scores: &[f64],
    introns: &IntronScores,
) -> ConsensusGene {
    let orient = exons[exon_indices[0]].orient;
    let mut score: f64 = exon_indices.iter().map(|&i| base_scores[i]).sum();
    for w in exon_indices.windows(2) {
        let a = &exons[w[0]];
        let b = &exons[w[1]];
        let intron = (a.coords.rend + 1, b.coords.lend - 1);
        score += introns.score(intron).unwrap_or(0.0);
    }
    ConsensusGene {
        exon_indices,
        orient,
        score,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consensus::exon::{ExonType, end_frame_for};
    use crate::model::Coordset;

    fn exon(lend: i64, rend: i64, et: ExonType, start_frame: u8) -> ExonCandidate {
        ExonCandidate {
            coords: Coordset::new(lend, rend),
            orient: Strand::Plus,
            exon_type: et,
            start_frame,
            end_frame: end_frame_for(start_frame, rend - lend + 1),
            // AA / AA -> any spanning codon is AAA (Lys), never a stop
            left_seq_boundary: [b'A', b'A'],
            right_seq_boundary: [b'A', b'A'],
            evidence: vec![("acc".into(), "src".into())],
        }
    }

    fn vectors() -> RegionVectors {
        RegionVectors::new(1, 1000)
    }

    #[test]
    fn score_exon_adds_coding_and_exon_specific() {
        let mut v = RegionVectors::new(1, 100);
        v.add_coverage(10, 19, 2.0); // 10 bases * 2 = 20 coding
        let w = Weights::parse_str("TRANSCRIPT src 3\n", "w").unwrap();
        let mut e = exon(10, 19, ExonType::Internal, 1);
        e.evidence = vec![("a".into(), "src".into())];
        // 20 coding + (transcript weight 3 * 10 unmasked bases) = 50
        assert_eq!(score_exon(&e, &v, &w), 50.0);
    }

    #[test]
    fn two_exon_gene_links_across_evidence_intron() {
        // initial 100..120 (len 21, end_frame 3) -> terminal 200..230 (start_frame 1)
        let exons = vec![
            exon(100, 120, ExonType::Initial, 1),
            exon(200, 230, ExonType::Terminal, 1),
        ];
        assert_eq!(exons[0].end_frame, 3); // (21)%3 == 0 -> 3, connects to start 1
        let base = vec![10.0, 15.0];
        let mut introns = IntronScores::default();
        introns.add((121, 199), 1.0, 5, "p", "fgenesh", true); // score 5
        let v = vectors();
        let code = GeneticCode::default();

        let genes = run_trellis(&exons, &base, (1, 300), &introns, &v, &code, 1.0, 500);
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].exon_indices, vec![0, 1]);
        assert_eq!(genes[0].orient, Strand::Plus);
        assert_eq!(genes[0].score, 30.0); // 10 + 15 + intron 5
    }

    #[test]
    fn missing_intron_blocks_the_link() {
        // same exons but no evidence intron -> they cannot be joined into one gene
        let exons = vec![
            exon(100, 120, ExonType::Initial, 1),
            exon(200, 230, ExonType::Terminal, 1),
        ];
        let base = vec![10.0, 15.0];
        let introns = IntronScores::default(); // empty
        let v = vectors();
        let code = GeneticCode::default();
        let genes = run_trellis(&exons, &base, (1, 300), &introns, &v, &code, 1.0, 500);
        // initial+ and terminal+ can't connect intergenically (not in the table), so the
        // best path is a single exon (the higher-scoring terminal)
        assert!(genes.iter().all(|g| g.exon_indices.len() == 1));
    }

    #[test]
    fn wrong_frame_blocks_the_link() {
        // terminal start_frame 2 violates the (3 -> 1) transition from the initial exon
        let exons = vec![
            exon(100, 120, ExonType::Initial, 1), // end_frame 3
            exon(200, 230, ExonType::Terminal, 2),
        ];
        let base = vec![10.0, 15.0];
        let mut introns = IntronScores::default();
        introns.add((121, 199), 1.0, 5, "p", "fgenesh", true);
        let v = vectors();
        let code = GeneticCode::default();
        let genes = run_trellis(&exons, &base, (1, 300), &introns, &v, &code, 1.0, 500);
        assert!(
            genes.iter().all(|g| g.exon_indices.len() == 1),
            "frame mismatch must not link"
        );
    }

    #[test]
    fn junction_stop_blocks_the_link() {
        // engineer a stop across the junction: end_frame 1 checks junction[1..4].
        // make a 22-bp initial exon so end_frame = (22)%3 = 1.
        let mut a = exon(100, 121, ExonType::Initial, 1);
        assert_eq!(a.end_frame, 1);
        // a.right = "xT", b.left = "AA" -> junction = x T A A; codon[1..4] = "TAA" = stop
        a.right_seq_boundary = [b'C', b'T'];
        let mut b = exon(200, 230, ExonType::Terminal, 2); // start_frame to satisfy (1 -> 2)
        b.left_seq_boundary = [b'A', b'A'];
        let exons = vec![a, b];
        let base = vec![10.0, 15.0];
        let mut introns = IntronScores::default();
        introns.add((122, 199), 1.0, 5, "p", "fgenesh", true);
        let v = vectors();
        let code = GeneticCode::default();
        let genes = run_trellis(&exons, &base, (1, 300), &introns, &v, &code, 1.0, 500);
        assert!(
            genes.iter().all(|g| g.exon_indices.len() == 1),
            "junction stop must not link"
        );
    }

    #[test]
    fn trellis_invariants_under_random_inputs() {
        // Deterministic pseudo-random, non-overlapping exon sets with random types,
        // frames and introns. The trellis must never panic, never reuse an exon across
        // genes, and only chain within-gene exons across a real (evidence-supported) intron.
        let mut seed: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut rng = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            (seed >> 33) as i64
        };
        for _ in 0..50 {
            let n = 3 + (rng() % 8) as usize;
            let mut exons = Vec::new();
            let mut pos = 1i64;
            for _ in 0..n {
                let lend = pos + 5 + rng() % 50;
                let rend = lend + 6 + rng() % 60;
                pos = rend;
                let et = match rng() % 4 {
                    0 => ExonType::Initial,
                    1 => ExonType::Internal,
                    2 => ExonType::Terminal,
                    _ => ExonType::Single,
                };
                exons.push(exon(lend, rend, et, (1 + rng() % 3) as u8));
            }
            let base: Vec<f64> = (0..exons.len()).map(|i| (i + 1) as f64).collect();
            let mut introns = IntronScores::default();
            for w in 0..exons.len().saturating_sub(1) {
                if rng() % 2 == 0 {
                    let (a, b) = (&exons[w], &exons[w + 1]);
                    introns.add(
                        (a.coords.rend + 1, b.coords.lend - 1),
                        1.0,
                        5,
                        "x",
                        "s",
                        true,
                    );
                }
            }
            let v = RegionVectors::new(1, (pos + 100) as usize);
            let genes = run_trellis(
                &exons,
                &base,
                (1, pos + 50),
                &introns,
                &v,
                &GeneticCode::default(),
                1.0,
                500,
            );

            let mut seen = std::collections::HashSet::new();
            for g in &genes {
                for &i in &g.exon_indices {
                    assert!(seen.insert(i), "an exon was reused across genes");
                }
                for w in g.exon_indices.windows(2) {
                    let (a, b) = (exons[w[0]].coords, exons[w[1]].coords);
                    assert!(
                        a.rend < b.lend,
                        "within-gene exons overlap or are out of order"
                    );
                    assert!(
                        introns.contains((a.rend + 1, b.lend - 1)),
                        "within-gene exons linked without a supporting intron"
                    );
                }
            }
        }
    }

    #[test]
    fn two_single_genes_split_via_intergenic() {
        // two single-exon genes; intergenic connection puts both on one path, then
        // traverse splits them at the single+ gene boundaries.
        let exons = vec![
            exon(100, 150, ExonType::Single, 1),
            exon(300, 350, ExonType::Single, 1),
        ];
        let base = vec![10.0, 10.0];
        let introns = IntronScores::default();
        let v = vectors();
        let code = GeneticCode::default();
        let genes = run_trellis(&exons, &base, (1, 400), &introns, &v, &code, 1.0, 500);
        assert_eq!(genes.len(), 2);
        assert_eq!(genes[0].exon_indices, vec![0]);
        assert_eq!(genes[1].exon_indices, vec![1]);
    }
}
