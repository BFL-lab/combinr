//! End-to-end orchestration: load sources → (filter) → cluster → assemble each
//! cluster in parallel, and (for alt-splice) group into loci and classify events.

use crate::altsplice::{AltSpliceResult, EventRecord, Isoform, Locus, analyze};
use crate::assemble::{ClusterAssembly, assemble_cluster};
use crate::cluster::cluster_alignments;
use crate::consensus::evidence::load_evidence;
use crate::consensus::region::build_regions;
use crate::consensus::{
    CalledGene, CandidateParams, EngineParams, EvClass, EvidenceChain, FilterParams, Weights,
    consensus_region,
};
use crate::error::{CombinrError, Result};
use crate::filter::{self, Filters};
use crate::io::fasta::Fasta;
use crate::io::load_sources;
use crate::io::out_model::OutGene;
use crate::orf::{GeneticCode, ReconcileResult, parse_cds_models, reconcile};
use rayon::prelude::*;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

/// PASA pairwise-compatibility fuzz default (bp) used when the consensus path
/// internally assembles its transcript alignments. Mirrors the CLI default.
pub const DEFAULT_FUZZLENGTH: i64 = 20;
/// EVM `terminal_intergenic_re_search` default (bp): minimum tail size to re-search.
pub const DEFAULT_RESEARCH_SIZE: i64 = 10_000;
/// Minimum consensus coding length (bp) below which a gene is flagged low-support.
pub const DEFAULT_MIN_CODING_LENGTH: i64 = 150;
/// Scaling applied to intergenic scores (EVM `INTERGENIC_SCORE_ADJUST_FACTOR`).
pub const DEFAULT_INTERGENIC_ADJUST: f64 = 1.0;

/// Resolve an NCBI genetic-code id, mapping the table-validation message into a
/// [`CombinrError::Parse`]. Single source for the consensus entry points.
pub fn parse_genetic_code(id: u32) -> Result<GeneticCode> {
    GeneticCode::from_ncbi_id(id).map_err(|msg| CombinrError::Parse {
        file: "genetic-code".into(),
        line: 0,
        msg,
    })
}

/// Load alignments from `paths`, apply optional quality filters, cluster
/// genome-wide, and assemble each cluster (in parallel) into a non-redundant
/// set. Provenance from all sources is preserved.
pub fn assemble_sources(
    paths: &[PathBuf],
    fuzzlength: i64,
    min_overlap_frac: f64,
    filters: &Filters,
) -> Result<Vec<ClusterAssembly>> {
    let alignments = load_sources(paths)?;
    let (alignments, _dropped) = filter::apply(alignments, filters);
    let clusters = cluster_alignments(&alignments, min_overlap_frac);

    // Clusters are independent → assemble in parallel.
    let per_cluster: Result<Vec<Vec<ClusterAssembly>>> = clusters
        .par_iter()
        .map(|cluster| {
            let members: Vec<_> = cluster
                .member_indices
                .iter()
                .map(|&i| alignments[i].clone())
                .collect();
            assemble_cluster(&members, fuzzlength)
        })
        .collect();

    Ok(per_cluster?.into_iter().flatten().collect())
}

/// Full Algorithm 2: assemble sources, then group into loci and classify
/// alternative-splicing events.
pub fn analyze_sources(
    paths: &[PathBuf],
    fuzzlength: i64,
    min_overlap_frac: f64,
    filters: &Filters,
) -> Result<AltSpliceResult> {
    let assemblies = assemble_sources(paths, fuzzlength, min_overlap_frac, filters)?;
    Ok(analyze(&assemblies, fuzzlength, min_overlap_frac))
}

/// Algorithm 2 + the optional ORF/UTR step: analyze, then graft an external
/// gene-prediction CDS onto each isoform and tag events by region. `code`
/// selects the genetic code for divergent-isoform stop-codon detection.
pub fn reconcile_sources(
    paths: &[PathBuf],
    gene_pred: &Path,
    genome: &Path,
    fuzzlength: i64,
    min_overlap_frac: f64,
    filters: &Filters,
    code: GeneticCode,
) -> Result<(Vec<Isoform>, Vec<Locus>, ReconcileResult)> {
    let asr = analyze_sources(paths, fuzzlength, min_overlap_frac, filters)?;
    let models = parse_cds_models(gene_pred)?;
    let genome = Fasta::load(genome)?;
    let recon = reconcile(
        &asr.isoforms,
        &asr.loci,
        asr.events,
        &models,
        &genome,
        &code,
    );
    Ok((asr.isoforms, asr.loci, recon))
}

/// Counts reported by [`augment_sources`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AugmentStats {
    /// Genes in the input gene-model GFF3 (every one is emitted).
    pub input_genes: usize,
    /// mRNAs in the input gene-model GFF3 (every one is emitted verbatim).
    pub input_mrnas: usize,
    /// Isoform mRNAs appended: output mRNAs minus input mRNAs.
    pub appended: usize,
}

/// `assemble --models`: analyze the transcript sources (assembly + alt-splice), load an
/// existing gene-model GFF3, and append each gene's genuine alternative isoforms as extra
/// mRNAs, the input genes kept verbatim (see [`crate::consensus::altsplice::augment`]).
/// Returns the output genes, the region-tagged alt-splice events and the counts.
pub fn augment_sources(
    paths: &[PathBuf],
    models: &Path,
    genome: &Path,
    fuzzlength: i64,
    min_overlap_frac: f64,
    filters: &Filters,
    code: GeneticCode,
) -> Result<(Vec<OutGene>, Vec<EventRecord>, AugmentStats)> {
    let asr = analyze_sources(paths, fuzzlength, min_overlap_frac, filters)?;
    let models = crate::io::gene_models::load_gene_models(models)?;
    let mrnas = |genes: &[OutGene]| genes.iter().map(|g| g.transcripts.len()).sum::<usize>();
    let (input_genes, input_mrnas) = (models.len(), mrnas(&models));
    let genome = Fasta::load(genome)?;
    let (genes, events) = crate::consensus::altsplice::augment(models, asr, &genome, &code);
    let stats = AugmentStats {
        input_genes,
        input_mrnas,
        appended: mrnas(&genes) - input_mrnas,
    };
    Ok((genes, events, stats))
}

/// Configuration for the EVM-style `consensus` path.
pub struct ConsensusConfig {
    pub weights: PathBuf,
    pub gene_predictions: Vec<PathBuf>,
    pub protein_alignments: Vec<PathBuf>,
    pub transcript_alignments: Vec<PathBuf>,
    pub genome: PathBuf,
    pub repeats: Option<PathBuf>,
    pub genetic_code: u32,
    pub flank: i64,
    pub strict: bool,
    pub max_prev_exons: usize,
    pub min_score_ratio: f64,
    pub min_intron_length: i64,
    pub research_size: i64,
    pub research_intergenic: i64,
    pub search_long_introns: i64,
    pub extend_terminal_stop: bool,
    pub peak_augment: bool,
    pub promote_transcript_orfs: bool,
    pub alt_splice: bool,
    pub min_coding_length: i64,
    /// PASA `--stringent_alignment_overlap` for the `--alt-splice` isoform grouping:
    /// isoforms share a gene only when their spans overlap `>=` this percent of the
    /// shorter span. 0.0 = any overlap (off). Does not affect the EVM region partitioner.
    pub stringent_overlap: f64,
}

/// Everything `consensus_inner` builds once, shared by both consensus entry points so
/// the `--alt-splice` path never re-loads the genome, re-parses the code, or re-assembles
/// the transcripts that the core run already produced.
struct ConsensusInner {
    genes: Vec<CalledGene>,
    genome: Fasta,
    code: GeneticCode,
    /// The transcript-assembly result, present iff `--promote-transcript-orfs` ran it.
    asr: Option<AltSpliceResult>,
}

/// Core consensus run: parse weights, ingest weighted evidence, cluster into per-contig
/// regions, run the both-strand trellis on each region in parallel, and (opt-in) promote
/// transcript-only loci via de-novo ORF — returning the genome/code/assembly it loaded.
fn consensus_inner(cfg: &ConsensusConfig) -> Result<ConsensusInner> {
    let weights = Weights::parse_file(&cfg.weights)?;
    let chains = load_evidence(
        &cfg.gene_predictions,
        &cfg.protein_alignments,
        &cfg.transcript_alignments,
        &weights,
    )?;
    let genome = Fasta::load(&cfg.genome)?;
    let code = parse_genetic_code(cfg.genetic_code)?;
    let cand = CandidateParams {
        min_intron_length: cfg.min_intron_length,
        code,
        extend_terminal_stop: cfg.extend_terminal_stop,
        peak_augment: cfg.peak_augment,
        sum_genepred_weights: weights.sum_genepred_weights(),
    };
    let filt = FilterParams {
        min_score_ratio: cfg.min_score_ratio,
        min_coding_length: cfg.min_coding_length,
        intergenic_adjust: DEFAULT_INTERGENIC_ADJUST,
        strict: cfg.strict,
    };
    let eng = EngineParams {
        max_prev_exons: cfg.max_prev_exons,
        research_size: cfg.research_size,
        research_intergenic: cfg.research_intergenic,
        search_long_introns: cfg.search_long_introns,
        intergenic_adjust: DEFAULT_INTERGENIC_ADJUST,
    };

    let mask = match &cfg.repeats {
        Some(path) => crate::consensus::repeats::parse_repeats(path)?,
        None => std::collections::HashMap::new(),
    };

    let regions = build_regions(&chains, cfg.flank);
    // Regions are independent → process in parallel (as `assemble_sources` does).
    let per_region: Vec<Vec<CalledGene>> = regions
        .par_iter()
        .map(|r| {
            // Only the repeats intersecting this region (the per-contig lists are sorted
            // and merged by `parse_repeats`; `set_masked` clips to the engine's window).
            let region_mask = mask.get(&r.contig).map_or(&[][..], |m| {
                crate::consensus::repeats::overlapping(m, r.span)
            });
            consensus_region(
                r,
                &chains,
                &genome,
                &cand,
                &filt,
                &weights,
                &eng,
                region_mask,
            )
        })
        .collect();
    let mut genes: Vec<CalledGene> = per_region.into_iter().flatten().collect();

    // Recover transcript-only loci via de-novo ORF (opt-in). Assemble the transcript
    // alignments (the PASA path) and promote loci with no overlapping consensus gene.
    let asr = if cfg.promote_transcript_orfs && !cfg.transcript_alignments.is_empty() {
        let asr = analyze_sources(
            &cfg.transcript_alignments,
            DEFAULT_FUZZLENGTH,
            cfg.stringent_overlap,
            &Filters::none(),
        )?;
        genes = crate::consensus::promote::promote_and_merge(
            genes,
            &asr,
            &genome,
            &code,
            cfg.min_coding_length,
            &transcript_sources(&chains),
        );
        Some(asr)
    } else {
        None
    };

    Ok(ConsensusInner {
        genes,
        genome,
        code,
        asr,
    })
}

/// Transcript accession → GFF column-2 source, from the loaded transcript chains: the
/// evidence attribution for promoted transcript-ORF genes.
fn transcript_sources(chains: &[EvidenceChain]) -> HashMap<String, String> {
    chains
        .iter()
        .filter(|c| c.ev_class == EvClass::Transcript)
        .map(|c| (c.accession.clone(), c.ev_type.clone()))
        .collect()
}

/// Build consensus gene models: parse weights, ingest weighted evidence, cluster into
/// per-contig regions, and run the both-strand trellis on each region in parallel.
pub fn consensus_sources(cfg: &ConsensusConfig) -> Result<Vec<CalledGene>> {
    Ok(consensus_inner(cfg)?.genes)
}

/// Output of [`consensus_with_isoforms`].
pub struct ConsensusIsoforms {
    /// The called consensus genes (for the evidence report; IDs via
    /// [`crate::consensus::output::ordered_with_ids`], identical to `out_genes`').
    pub genes: Vec<CalledGene>,
    /// Consensus mRNA + alternative-isoform mRNAs per gene.
    pub out_genes: Vec<OutGene>,
    /// Region-tagged alt-splice events.
    pub events: Vec<EventRecord>,
}

/// `consensus --alt-splice`: build the consensus genes, then attach each locus's
/// alternative transcript isoforms as extra mRNAs (CDS grafted from the consensus). Returns
/// the called genes, the output genes and the region-tagged alt-splice events.
pub fn consensus_with_isoforms(cfg: &ConsensusConfig) -> Result<ConsensusIsoforms> {
    let ConsensusInner {
        genes,
        genome,
        code,
        asr,
    } = consensus_inner(cfg)?;
    // Reuse the transcript assembly when `--promote` already produced it; otherwise build
    // it now (identical args, so the output is unchanged either way).
    let asr = match asr {
        Some(asr) => asr,
        None => analyze_sources(
            &cfg.transcript_alignments,
            DEFAULT_FUZZLENGTH,
            cfg.stringent_overlap,
            &Filters::none(),
        )?,
    };
    let (out_genes, events) = crate::consensus::altsplice::annotate(&genes, asr, &genome, &code);
    Ok(ConsensusIsoforms {
        genes,
        out_genes,
        events,
    })
}
