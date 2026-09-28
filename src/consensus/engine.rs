//! Per-region consensus: run the trellis on both strands and recurse into large tails.
//!
//! The minus strand is handled by **reverse-complementing the region and running the
//! exact same forward machinery**, then transposing the resulting genes back to forward
//! coordinates — so the trellis and compatibility logic never need per-strand branches.
//! Forward and minus are explored independently across the whole region (leaning
//! anti-false-negative; see `prefer-flag-not-drop-genes`), rather than EVM's single
//! combined trellis that picks one strand per antisense locus.
//!
//! Recursion mirrors `generate_consensus_gene_predictions` (~483): after the best path,
//! re-search the left/right tails beyond the *retained* genes' span (≥ `research_size`),
//! which is also the fix for EVM's bug where an eliminated gene's span blocked re-search.

use crate::consensus::candidates::{CandidateParams, RegionData, build_candidates};
use crate::consensus::evidence::EvidenceChain;
use crate::consensus::exon::{ExonType, frame_offset};
use crate::consensus::filter::{FilterParams, SupportFlags, assess};
use crate::consensus::region::ConsensusRegion;
use crate::consensus::trellis::{ConsensusGene, run_trellis, score_all_exons};
use crate::consensus::weights::Weights;
use crate::io::fasta::Fasta;
use crate::model::{Coordset, Strand};
use crate::orf::GeneticCode;
use crate::orf::coords::SplicedTranscript;
use crate::orf::translate::reverse_complement;

/// Trellis/recursion knobs.
#[derive(Clone, Copy, Debug)]
pub struct EngineParams {
    pub max_prev_exons: usize,
    /// Minimum tail size (bp) to re-search (EVM `terminal_intergenic_re_search`, 10000).
    pub research_size: i64,
    /// Minimum intergenic-gap size (bp) between called genes to re-search for additional
    /// genes (EVM `--re_search_intergenic`). 0 = off (the EVM default).
    pub research_intergenic: i64,
    /// Minimum intron length (bp) to re-search for nested genes (EVM
    /// `--search_long_introns`). 0 = off (the EVM default).
    pub search_long_introns: i64,
    pub intergenic_adjust: f64,
}

/// A consensus gene resolved to forward genomic coordinates, with its coding region
/// projected (CDS + flanking UTR).
#[derive(Debug, Clone, PartialEq)]
pub struct CalledGene {
    pub contig: String,
    pub orient: Strand,
    pub exons: Vec<Coordset>,
    pub cds: Vec<Coordset>,
    pub five_utr: Vec<Coordset>,
    pub three_utr: Vec<Coordset>,
    pub partial5: bool,
    pub partial3: bool,
    /// GFF3 phase of the 5'-most CDS base: the number of bases before the first complete
    /// codon. Nonzero only for a 5'-partial gene, whose CDS starts at its first exon base
    /// (the partial codon is coding, not UTR); 0 for a gene with a start codon.
    pub cds_start_phase: u8,
    pub score: f64,
    pub support: SupportFlags,
    /// True for a transcript-ORF gene recovered at a locus with no consensus CDS
    /// (`--promote-transcript-orfs`), vs. a trellis-derived consensus gene.
    pub promoted: bool,
    /// Per-exon / per-intron evidence attribution for the evidence report, in ascending
    /// forward `lend` order.
    pub features: Vec<FeatureSupport>,
}

/// One exon or intron of a called gene with the evidence supporting it (forward genomic
/// coordinates). Rendered by [`crate::consensus::report`].
#[derive(Debug, Clone, PartialEq)]
pub struct FeatureSupport {
    pub coords: Coordset,
    pub kind: FeatureKind,
    /// `(accession, ev_type)` of each supporting evidence chain, in insertion order.
    pub evidence: Vec<(String, String)>,
}

/// An exon (typed, with EVM frames: 1-3 forward, 4-6 reverse) or an intron.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FeatureKind {
    Exon {
        exon_type: ExonType,
        start_frame: u8,
        end_frame: u8,
    },
    Intron,
}

impl CalledGene {
    /// Genomic span `(min lend, max rend)` over the gene's exons. Robust to ordering
    /// (the exons are lend-sorted in practice, but this does not rely on it).
    pub fn span(&self) -> (i64, i64) {
        let lend = self.exons.iter().map(|c| c.lend).min().unwrap_or(0);
        let rend = self.exons.iter().map(|c| c.rend).max().unwrap_or(0);
        (lend, rend)
    }
}

const MAX_RECURSION_DEPTH: usize = 32;

/// Find consensus genes in `region` on both strands.
#[allow(clippy::too_many_arguments)]
pub fn consensus_region(
    region: &ConsensusRegion,
    chains: &[EvidenceChain],
    genome: &Fasta,
    cand: &CandidateParams,
    filt: &FilterParams,
    weights: &Weights,
    eng: &EngineParams,
    mask: &[Coordset],
) -> Vec<CalledGene> {
    // clamp the region span to the contig
    let contig_len = genome
        .contig_len(&region.contig)
        .map(|n| n as i64)
        .unwrap_or(region.span.rend);
    let lo = region.span.lend.max(1);
    let hi = region.span.rend.min(contig_len);
    if hi < lo {
        return Vec::new();
    }
    let region = ConsensusRegion {
        contig: region.contig.clone(),
        span: Coordset::new(lo, hi),
        chain_indices: region.chain_indices.clone(),
    };

    let mut out = Vec::new();

    // ---- forward strand ----
    let fwd = build_candidates(&region, chains, genome, cand, mask);
    let fwd_base = score_all_exons(&fwd.exons, &fwd.vectors, weights);
    let mut fwd_genes = Vec::new();
    collect_genes(
        &fwd,
        &fwd_base,
        (lo, hi),
        &cand.code,
        filt,
        eng,
        &mut fwd_genes,
        0,
    );
    resolve_all(
        &fwd_genes,
        &fwd,
        &region.contig,
        Strand::Plus,
        None,
        genome,
        &cand.code,
        &mut out,
    );

    // ---- minus strand: reverse-complement the region, run forward, transpose back ----
    if let Some((rc_fasta, rc_region, rc_chains, rc_mask)) =
        build_rc_region(&region, chains, mask, genome, lo, hi)
    {
        let len = hi - lo + 1;
        let rc = build_candidates(&rc_region, &rc_chains, &rc_fasta, cand, &rc_mask);
        let rc_base = score_all_exons(&rc.exons, &rc.vectors, weights);
        let mut rc_genes = Vec::new();
        collect_genes(
            &rc,
            &rc_base,
            (1, len),
            &cand.code,
            filt,
            eng,
            &mut rc_genes,
            0,
        );
        resolve_all(
            &rc_genes,
            &rc,
            &region.contig,
            Strand::Minus,
            Some(hi),
            genome,
            &cand.code,
            &mut out,
        );
    }

    dedup_genes(out)
}

/// Reverse-complement transform for the minus strand: build the RC genome window over
/// `[lo, hi]`, transpose every minus-strand chain and the repeat mask into RC-local
/// coordinates (via [`Coordset::reflect`]), and return the RC region ready to run through
/// the same forward machinery. `None` when the genome window is unavailable.
fn build_rc_region(
    region: &ConsensusRegion,
    chains: &[EvidenceChain],
    mask: &[Coordset],
    genome: &Fasta,
    lo: i64,
    hi: i64,
) -> Option<(Fasta, ConsensusRegion, Vec<EvidenceChain>, Vec<Coordset>)> {
    let fwd_bytes = genome.subseq(&region.contig, lo, hi)?;
    let len = hi - lo + 1;
    let rc_fasta = Fasta::from_seq("rc", reverse_complement(fwd_bytes));

    let mut rc_chains = Vec::new();
    for &ci in &region.chain_indices {
        let c = &chains[ci];
        if c.orient != Strand::Minus {
            continue;
        }
        let mut links: Vec<Coordset> = c.links.iter().map(|s| s.reflect(hi)).collect();
        links.sort_by_key(|c| c.lend);
        let span = Coordset::new(links.first().unwrap().lend, links.last().unwrap().rend);
        rc_chains.push(EvidenceChain {
            accession: c.accession.clone(),
            ev_type: c.ev_type.clone(),
            ev_class: c.ev_class,
            weight: c.weight,
            contig: "rc".into(),
            orient: Strand::Plus,
            span,
            links,
        });
    }

    let rc_region = ConsensusRegion {
        contig: "rc".into(),
        span: Coordset::new(1, len),
        chain_indices: (0..rc_chains.len()).collect(),
    };
    let rc_mask: Vec<Coordset> = mask.iter().map(|iv| iv.reflect(hi)).collect();
    Some((rc_fasta, rc_region, rc_chains, rc_mask))
}

/// Resolve a strand's trellis genes to forward coordinates and append them to `out`.
#[allow(clippy::too_many_arguments)]
fn resolve_all(
    genes: &[(ConsensusGene, SupportFlags)],
    rd: &RegionData,
    contig: &str,
    orient: Strand,
    transpose_hi: Option<i64>,
    genome: &Fasta,
    code: &GeneticCode,
    out: &mut Vec<CalledGene>,
) {
    for (g, support) in genes {
        out.push(resolve(
            g,
            rd,
            contig,
            orient,
            transpose_hi,
            *support,
            genome,
            code,
        ));
    }
}

/// Drop genes with an identical (orientation, exon structure) — re-search can rediscover
/// the same gene via overlapping tail/gap/intron sub-ranges.
fn dedup_genes(mut genes: Vec<CalledGene>) -> Vec<CalledGene> {
    let mut seen = std::collections::HashSet::new();
    genes.retain(|g| seen.insert((g.orient, g.exons.clone())));
    genes
}

/// Resolve a trellis gene to forward genomic coordinates (transposing RC-local minus
/// coordinates when `transpose_hi` is set) and project its CDS/UTR.
///
/// The consensus exons ARE the coding structure: build a [`SplicedTranscript`] over them
/// and project the ORF from the gene's 5' start, in its reading frame, to the first
/// in-frame stop (reusing `orf::project_orf` — decision 6). For a clean model the CDS
/// equals the exon structure; a non-canonical terminus without a real stop is 3'-partial.
#[allow(clippy::too_many_arguments)]
fn resolve(
    g: &ConsensusGene,
    rd: &RegionData,
    contig: &str,
    orient: Strand,
    transpose_hi: Option<i64>,
    support: SupportFlags,
    genome: &Fasta,
    code: &GeneticCode,
) -> CalledGene {
    let mut exons: Vec<Coordset> = g
        .exon_indices
        .iter()
        .map(|&i| {
            let c = rd.exons[i].coords;
            match transpose_hi {
                Some(hi) => c.reflect(hi),
                None => c,
            }
        })
        .collect();
    exons.sort_by_key(|c| c.lend);

    // 5'/3'-most exons in the gene's reading order (path order is 5'→3' in each pass).
    let first = &rd.exons[g.exon_indices[0]];
    let last = &rd.exons[*g.exon_indices.last().unwrap()];
    let has_start = matches!(first.exon_type, ExonType::Initial | ExonType::Single);
    let has_stop = matches!(last.exon_type, ExonType::Terminal | ExonType::Single);
    // first complete codon within the spliced transcript: a 5'-partial path begins with
    // a partial codon of `p` bases (0 for an initial/single exon, which is always frame 1)
    let p = frame_offset(first.start_frame);

    let st = SplicedTranscript::new(&exons, orient);
    let (cds, five_utr, three_utr, hit_stop, cds_start_phase) = match st.sequence(genome, contig) {
        Some(seq) => {
            // the ORF (and its stop) is read from the first complete codon ...
            let proj = st.project_orf(&seq, p.min(seq.len()), code);
            // ... but a 5'-partial CDS starts at the transcript's first base, the leading
            // partial codon carried as the first CDS row's phase (GFF3), not as a 5'UTR
            let (cds_t_start, phase) = if !has_start && proj.cds_t_end > proj.cds_t_start {
                (0, p as u8)
            } else {
                (proj.cds_t_start, 0)
            };
            let (cds, five_utr, three_utr) = st.cds_and_utrs(cds_t_start, proj.cds_t_end);
            (cds, five_utr, three_utr, proj.hit_stop, phase)
        }
        None => {
            let phase = if has_start { 0 } else { p as u8 };
            (exons.clone(), Vec::new(), Vec::new(), has_stop, phase)
        }
    };

    CalledGene {
        contig: contig.into(),
        orient,
        exons,
        cds,
        five_utr,
        three_utr,
        partial5: !has_start,
        partial3: !has_stop || !hit_stop,
        cds_start_phase,
        score: g.score,
        support,
        promoted: false,
        features: path_features(g, rd, transpose_hi),
    }
}

/// The evidence attribution of a trellis gene's path: each path exon (its candidate's
/// type, frames and evidence) and each intron between consecutive path exons. Intron
/// evidence is looked up in the pass's WORKING coordinates (where `rd.introns` is keyed),
/// and only then are the coordinates transposed to forward space for the minus pass
/// (whose frames become EVM's reverse 4-6).
fn path_features(
    g: &ConsensusGene,
    rd: &RegionData,
    transpose_hi: Option<i64>,
) -> Vec<FeatureSupport> {
    let to_fwd = |c: Coordset| transpose_hi.map_or(c, |hi| c.reflect(hi));
    let frame_shift = if transpose_hi.is_some() { 3 } else { 0 };
    let mut features = Vec::with_capacity(2 * g.exon_indices.len());
    for (k, &i) in g.exon_indices.iter().enumerate() {
        let ex = &rd.exons[i];
        if k > 0 {
            let gap = rd.exons[g.exon_indices[k - 1]].coords.gap_to(&ex.coords);
            features.push(FeatureSupport {
                coords: to_fwd(gap),
                kind: FeatureKind::Intron,
                evidence: rd.introns.evidence((gap.lend, gap.rend)).to_vec(),
            });
        }
        features.push(FeatureSupport {
            coords: to_fwd(ex.coords),
            kind: FeatureKind::Exon {
                exon_type: ex.exon_type,
                start_frame: ex.start_frame + frame_shift,
                end_frame: ex.end_frame + frame_shift,
            },
            evidence: ex.evidence.clone(),
        });
    }
    features.sort_by_key(|f| f.coords.lend);
    features
}

/// Build the best path over `range`, flag low support, keep the retained genes, then
/// re-search the tails beyond their span. Flag-not-drop keeps everything by default;
/// `--strict` drops low-support genes, and only retained spans block re-search.
#[allow(clippy::too_many_arguments)]
fn collect_genes(
    rd: &RegionData,
    base: &[f64],
    range: (i64, i64),
    code: &GeneticCode,
    filt: &FilterParams,
    eng: &EngineParams,
    out: &mut Vec<(ConsensusGene, SupportFlags)>,
    depth: usize,
) {
    if depth >= MAX_RECURSION_DEPTH {
        return;
    }
    let genes = run_trellis(
        &rd.exons,
        base,
        range,
        &rd.introns,
        &rd.vectors,
        code,
        eng.intergenic_adjust,
        eng.max_prev_exons,
    );
    if genes.is_empty() {
        return;
    }

    let mut retained: Vec<(ConsensusGene, SupportFlags)> = Vec::new();
    for g in genes {
        let support = assess(&g, &rd.exons, &rd.vectors, filt);
        if filt.strict && support.low_support {
            continue; // EVM-style drop
        }
        retained.push((g, support));
    }
    if retained.is_empty() {
        return;
    }

    // Per-gene spans and (optionally) long introns, captured before `retained` moves.
    let mut gene_spans: Vec<(i64, i64)> = retained
        .iter()
        .map(|(g, _)| {
            let l = g
                .exon_indices
                .iter()
                .map(|&i| rd.exons[i].coords.lend)
                .min()
                .unwrap();
            let r = g
                .exon_indices
                .iter()
                .map(|&i| rd.exons[i].coords.rend)
                .max()
                .unwrap();
            (l, r)
        })
        .collect();
    let long_introns: Vec<(i64, i64)> = if eng.search_long_introns > 0 {
        retained
            .iter()
            .flat_map(|(g, _)| {
                g.exon_indices.windows(2).filter_map(|w| {
                    let intron = rd.exons[w[0]].coords.gap_to(&rd.exons[w[1]].coords);
                    (intron.len() >= eng.search_long_introns).then_some((intron.lend, intron.rend))
                })
            })
            .collect()
    } else {
        Vec::new()
    };

    let span_lend = gene_spans.iter().map(|s| s.0).min().unwrap();
    let span_rend = gene_spans.iter().map(|s| s.1).max().unwrap();
    out.extend(retained);

    // left / right tail re-search (EVM default behaviour, gated on `research_size`)
    if span_lend - range.0 >= eng.research_size {
        collect_genes(
            rd,
            base,
            (range.0, span_lend - 1),
            code,
            filt,
            eng,
            out,
            depth + 1,
        );
    }
    if range.1 - span_rend >= eng.research_size {
        collect_genes(
            rd,
            base,
            (span_rend + 1, range.1),
            code,
            filt,
            eng,
            out,
            depth + 1,
        );
    }

    // intergenic-gap re-search between called genes (EVM `--re_search_intergenic`, opt-in)
    if eng.research_intergenic > 0 {
        gene_spans.sort_unstable();
        for w in gene_spans.windows(2) {
            let (gap_lend, gap_rend) = (w[0].1 + 1, w[1].0 - 1);
            if gap_rend - gap_lend + 1 >= eng.research_intergenic {
                collect_genes(
                    rd,
                    base,
                    (gap_lend, gap_rend),
                    code,
                    filt,
                    eng,
                    out,
                    depth + 1,
                );
            }
        }
    }

    // long-intron re-search for nested genes (EVM `--search_long_introns`, opt-in)
    for &(il, ir) in &long_introns {
        collect_genes(rd, base, (il, ir), code, filt, eng, out, depth + 1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consensus::weights::EvClass;
    use std::io::Write;

    fn write_fasta(seq: &[u8]) -> Fasta {
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        writeln!(tmp, ">chr1").unwrap();
        tmp.write_all(seq).unwrap();
        writeln!(tmp).unwrap();
        Fasta::load(tmp.path()).unwrap()
    }

    fn chain(class: EvClass, weight: f64, orient: Strand, segs: &[(i64, i64)]) -> EvidenceChain {
        let links: Vec<Coordset> = segs.iter().map(|&(l, r)| Coordset::new(l, r)).collect();
        let span = Coordset::new(links.first().unwrap().lend, links.last().unwrap().rend);
        EvidenceChain {
            accession: "g".into(),
            ev_type: "fgenesh".into(),
            ev_class: class,
            weight,
            contig: "chr1".into(),
            orient,
            span,
            links,
        }
    }

    fn config() -> (CandidateParams, FilterParams, Weights, EngineParams) {
        (
            CandidateParams {
                min_intron_length: 20,
                code: GeneticCode::default(),
                extend_terminal_stop: false,
                peak_augment: false,
                sum_genepred_weights: 0.0,
            },
            FilterParams {
                min_score_ratio: 0.75,
                min_coding_length: 150,
                intergenic_adjust: 1.0,
                strict: false,
            },
            Weights::parse_str("ABINITIO_PREDICTION fgenesh 1\n", "w").unwrap(),
            EngineParams {
                max_prev_exons: 500,
                research_size: 10_000,
                research_intergenic: 0,
                search_long_introns: 0,
                intergenic_adjust: 1.0,
            },
        )
    }

    #[test]
    fn forward_minus_strand_gene_is_transposed_back() {
        // A 200-bp genome of 'C' with a minus-strand single-exon prediction 50..120.
        // On the minus strand the coding reads revcomp(C-run)=G-run -> no stop, so the
        // single exon survives and must come back at its forward coordinates with '-'.
        let g = write_fasta(&[b'C'; 200]);
        let chains = vec![chain(
            EvClass::AbinitioPrediction,
            1.0,
            Strand::Minus,
            &[(50, 120)],
        )];
        let region = ConsensusRegion {
            contig: "chr1".into(),
            span: Coordset::new(1, 200),
            chain_indices: vec![0],
        };
        let (cand, filt, w, eng) = config();
        let genes = consensus_region(&region, &chains, &g, &cand, &filt, &w, &eng, &[]);
        let minus: Vec<_> = genes.iter().filter(|g| g.orient == Strand::Minus).collect();
        assert_eq!(minus.len(), 1, "the minus-strand gene is called");
        assert_eq!(
            minus[0].exons,
            vec![Coordset::new(50, 120)],
            "transposed back to forward coords"
        );
    }

    #[test]
    fn forward_strand_gene_called_in_place() {
        let g = write_fasta(&[b'C'; 200]);
        // forward single-exon prediction; revcomp not needed
        let chains = vec![chain(
            EvClass::AbinitioPrediction,
            1.0,
            Strand::Plus,
            &[(50, 120)],
        )];
        let region = ConsensusRegion {
            contig: "chr1".into(),
            span: Coordset::new(1, 200),
            chain_indices: vec![0],
        };
        let (cand, filt, w, eng) = config();
        let genes = consensus_region(&region, &chains, &g, &cand, &filt, &w, &eng, &[]);
        let plus: Vec<_> = genes.iter().filter(|g| g.orient == Strand::Plus).collect();
        assert_eq!(plus.len(), 1);
        assert_eq!(plus[0].exons, vec![Coordset::new(50, 120)]);
    }

    fn genome_with(len: usize, edits: &[(i64, &[u8])]) -> Fasta {
        let mut seq = vec![b'C'; len];
        for &(pos1, bytes) in edits {
            for (k, &b) in bytes.iter().enumerate() {
                seq[(pos1 - 1) as usize + k] = b;
            }
        }
        write_fasta(&seq)
    }

    fn call_single(g: &Fasta, segs: &[(i64, i64)], orient: Strand) -> Option<CalledGene> {
        let chains = vec![chain(EvClass::AbinitioPrediction, 1.0, orient, segs)];
        let region = ConsensusRegion {
            contig: "chr1".into(),
            span: Coordset::new(1, 200),
            chain_indices: vec![0],
        };
        let (cand, filt, w, eng) = config();
        consensus_region(&region, &chains, g, &cand, &filt, &w, &eng, &[])
            .into_iter()
            .find(|x| x.orient == orient)
    }

    fn ev_g() -> Vec<(String, String)> {
        vec![("g".to_string(), "fgenesh".to_string())]
    }

    fn exon_feat(l: i64, r: i64, t: ExonType, sf: u8, ef: u8) -> FeatureSupport {
        FeatureSupport {
            coords: Coordset::new(l, r),
            kind: FeatureKind::Exon {
                exon_type: t,
                start_frame: sf,
                end_frame: ef,
            },
            evidence: ev_g(),
        }
    }

    fn intron_feat(l: i64, r: i64) -> FeatureSupport {
        FeatureSupport {
            coords: Coordset::new(l, r),
            kind: FeatureKind::Intron,
            evidence: ev_g(),
        }
    }

    #[test]
    fn plus_gene_features_carry_exon_and_intron_evidence() {
        let g = genome_with(200, &[]);
        let gene = call_single(&g, &[(10, 18), (50, 60)], Strand::Plus).unwrap();
        assert_eq!(
            gene.features,
            vec![
                exon_feat(10, 18, ExonType::Initial, 1, 3),
                intron_feat(19, 49),
                exon_feat(50, 60, ExonType::Terminal, 1, 2),
            ]
        );
    }

    #[test]
    fn minus_gene_features_transposed_with_reverse_frames() {
        // Forward exons A=50..58 (3', terminal) and B=89..98 (5', initial), intron 59..88.
        // In the RC-local pass (hi = 200) the intron is keyed 113..142, so a lookup with
        // the forward coordinates would find no evidence. Frames: B initial 1->1 (len 10),
        // A terminal starts at cum 10 -> 2, ends 1; reverse encoding adds 3.
        let g = genome_with(200, &[]);
        let gene = call_single(&g, &[(50, 58), (89, 98)], Strand::Minus).unwrap();
        assert_eq!(
            gene.exons,
            vec![Coordset::new(50, 58), Coordset::new(89, 98)]
        );
        assert_eq!(
            gene.features,
            vec![
                exon_feat(50, 58, ExonType::Terminal, 5, 4),
                intron_feat(59, 88),
                exon_feat(89, 98, ExonType::Initial, 4, 4),
            ]
        );
    }

    #[test]
    fn clean_forward_orf_yields_cds_no_utr() {
        // 10..21 = ATG AAA CCC TAA (stop at 19..21)
        let g = genome_with(40, &[(10, b"ATGAAACCCTAA")]);
        let gene = call_single(&g, &[(10, 21)], Strand::Plus).unwrap();
        assert_eq!(gene.cds, vec![Coordset::new(10, 21)]);
        assert!(gene.five_utr.is_empty() && gene.three_utr.is_empty());
        assert!(!gene.partial5 && !gene.partial3);
        assert_eq!(
            gene.cds_start_phase, 0,
            "a gene with a start codon has phase 0"
        );
    }

    #[test]
    fn no_stop_codon_is_three_prime_partial() {
        // all-C exon: translation never hits a stop
        let g = genome_with(40, &[]);
        let gene = call_single(&g, &[(10, 21)], Strand::Plus).unwrap();
        assert!(gene.partial3, "no in-frame stop -> 3' partial");
    }

    #[test]
    fn minus_gene_cds_projected_in_forward_coords() {
        // minus single-exon, length a multiple of 3 (50..118 = 69 bp) so no remainder
        let g = genome_with(200, &[]);
        let gene = call_single(&g, &[(50, 118)], Strand::Minus).unwrap();
        assert_eq!(gene.orient, Strand::Minus);
        assert_eq!(gene.cds, vec![Coordset::new(50, 118)]);
    }

    #[test]
    fn long_intron_research_finds_nested_gene() {
        // Gene A: two CDS exons across a ~34 kb intron (its intron score dominates, so
        // the global best path is A alone). Gene B is a single exon nested inside that
        // intron, unreachable from A's exons by the grammar.
        let g = genome_with(40_000, &[]);
        let a = chain(
            EvClass::AbinitioPrediction,
            1.0,
            Strand::Plus,
            &[(1000, 1100), (35000, 35100)],
        );
        let b = chain(
            EvClass::AbinitioPrediction,
            1.0,
            Strand::Plus,
            &[(10000, 10300)],
        );
        let chains = vec![a, b];
        let region = ConsensusRegion {
            contig: "chr1".into(),
            span: Coordset::new(1, 40_000),
            chain_indices: vec![0, 1],
        };
        let (cand, filt, w, mut eng) = config();

        let has_b = |genes: &[CalledGene]| {
            genes
                .iter()
                .any(|x| x.exons.iter().any(|c| c.lend == 10000))
        };

        // off (default): B is hidden inside A's long intron
        let off = consensus_region(&region, &chains, &g, &cand, &filt, &w, &eng, &[]);
        assert!(
            !has_b(&off),
            "nested gene hidden when long-intron re-search is off"
        );

        // on: re-search introns >= 20 kb finds B
        eng.search_long_introns = 20_000;
        let on = consensus_region(&region, &chains, &g, &cand, &filt, &w, &eng, &[]);
        assert!(has_b(&on), "long-intron re-search recovers the nested gene");
    }

    /// A hand-built region (working coordinates) holding the given path exons, and the
    /// trellis gene that walks them in order.
    fn partial_path(len: usize, path: &[(i64, i64, ExonType, u8)]) -> (RegionData, ConsensusGene) {
        use crate::consensus::exon::{ExonCandidate, end_frame_for};
        use crate::consensus::vectors::{IntronScores, RegionVectors};
        let exons = path
            .iter()
            .map(|&(l, r, t, sf)| ExonCandidate {
                coords: Coordset::new(l, r),
                orient: Strand::Plus,
                exon_type: t,
                start_frame: sf,
                end_frame: end_frame_for(sf, r - l + 1),
                left_seq_boundary: *b"CC",
                right_seq_boundary: *b"CC",
                evidence: ev_g(),
            })
            .collect();
        let rd = RegionData {
            contig: "chr1".into(),
            span: Coordset::new(1, len as i64),
            exons,
            vectors: RegionVectors::new(1, len),
            introns: IntronScores::default(),
        };
        let g = ConsensusGene {
            exon_indices: (0..path.len()).collect(),
            orient: Strand::Plus,
            score: 1.0,
        };
        (rd, g)
    }

    fn resolve_partial(
        genome: &Fasta,
        len: usize,
        path: &[(i64, i64, ExonType, u8)],
        orient: Strand,
        transpose_hi: Option<i64>,
    ) -> CalledGene {
        let (rd, g) = partial_path(len, path);
        let support = SupportFlags {
            raw_noncoding: 0.0,
            noncoding_equivalent: 0.0,
            score_ratio: 0.0,
            coding_length: 0,
            low_support: false,
        };
        resolve(
            &g,
            &rd,
            "chr1",
            orient,
            transpose_hi,
            support,
            genome,
            &GeneticCode::default(),
        )
    }

    #[test]
    fn five_prime_partial_frame2_starts_cds_at_exon_with_phase_2() {
        // Internal 10..30 (start frame 2: its first base is codon position 2) -> Terminal
        // 50..60. The next codon starts 2 bases in, at 12: AGC CCC CCC TAA (stop 21..23).
        // Reading 1 base in (from 11) would hit TAG at 11..13 at once. The CDS starts at
        // the exon's first base (10) with phase 2; its 3' end (the stop) is unchanged.
        let g = genome_with(80, &[(11, b"TAG"), (21, b"TAA")]);
        let path = [
            (10, 30, ExonType::Internal, 2),
            (50, 60, ExonType::Terminal, 1),
        ];
        let gene = resolve_partial(&g, 80, &path, Strand::Plus, None);
        assert_eq!(
            gene.cds,
            vec![Coordset::new(10, 23)],
            "CDS from exon.lend, ORF read from exon.lend + 2 (stop at 21..23)"
        );
        assert!(
            gene.five_utr.is_empty(),
            "the partial codon is coding, not UTR"
        );
        assert_eq!(gene.cds_start_phase, 2);
        assert!(gene.partial5, "path starts at an internal exon");
        assert!(!gene.partial3, "terminal exon + in-frame stop");
    }

    #[test]
    fn five_prime_partial_frame3_starts_cds_at_exon_with_phase_1() {
        // start frame 3: the next codon starts 1 base in, at 11: CTA GCC CCC TAA (stop
        // 20..22). Reading 2 bases in (from 12) would hit TAG at 12..14 at once.
        let g = genome_with(80, &[(12, b"TAG"), (20, b"TAA")]);
        let path = [
            (10, 30, ExonType::Internal, 3),
            (50, 60, ExonType::Terminal, 1),
        ];
        let gene = resolve_partial(&g, 80, &path, Strand::Plus, None);
        assert_eq!(
            gene.cds,
            vec![Coordset::new(10, 22)],
            "CDS from exon.lend, ORF read from exon.lend + 1 (stop at 20..22)"
        );
        assert!(gene.five_utr.is_empty());
        assert_eq!(gene.cds_start_phase, 1);
        assert!(gene.partial5);
        assert!(!gene.partial3);
    }

    #[test]
    fn five_prime_partial_frame2_minus_strand() {
        // The frame-2 case in the RC-local pass (hi = 200): Internal 11..31 (sf 2) ->
        // Terminal 51..61, i.e. forward 170..190 (5') and 140..150. RC-local codons from
        // 13: AGC CCC CCC TAA (stop 22..24); from 12 it would hit TAG at 12..14 at once.
        // The forward genome is the reverse complement of that RC-local sequence.
        let mut rc = vec![b'C'; 200];
        rc[11..14].copy_from_slice(b"TAG"); // RC-local 12..14
        rc[21..24].copy_from_slice(b"TAA"); // RC-local 22..24
        let g = write_fasta(&reverse_complement(&rc));
        let path = [
            (11, 31, ExonType::Internal, 2),
            (51, 61, ExonType::Terminal, 1),
        ];
        let gene = resolve_partial(&g, 200, &path, Strand::Minus, Some(200));
        assert_eq!(
            gene.exons,
            vec![Coordset::new(140, 150), Coordset::new(170, 190)]
        );
        // ORF RC-local 13..24 -> forward 177..188 (first complete codon at forward
        // rend - 2 = 188); the CDS starts at the exon's 5' base (forward 190), phase 2
        assert_eq!(gene.cds, vec![Coordset::new(177, 190)]);
        assert!(gene.five_utr.is_empty());
        assert_eq!(gene.cds_start_phase, 2);
        assert!(gene.partial5);
        assert!(!gene.partial3);
    }
}
