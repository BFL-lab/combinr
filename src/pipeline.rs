//! End-to-end orchestration: load sources → (filter) → cluster → assemble each
//! cluster in parallel, and (for alt-splice) group into loci and classify events.

use crate::altsplice::{AltSpliceResult, EventRecord, Isoform, Locus, analyze};
use crate::assemble::{ClusterAssembly, assemble_cluster};
use crate::cluster::cluster_alignments;
use crate::consensus::evidence::load_evidence;
use crate::consensus::region::build_regions;
use crate::consensus::{
    CalledGene, CandidateParams, EngineParams, FilterParams, Weights, consensus_region,
};
use crate::error::{CombinrError, Result};
use crate::filter::{self, Filters};
use crate::io::fasta::Fasta;
use crate::io::load_sources;
use crate::io::out_model::OutGene;
use crate::orf::{GeneticCode, ReconcileResult, parse_cds_models, reconcile};
use rayon::prelude::*;
use std::path::{Path, PathBuf};

/// Load alignments from `paths`, apply optional quality filters, cluster
/// genome-wide, and assemble each cluster (in parallel) into a non-redundant
/// set. Provenance from all sources is preserved.
pub fn assemble_sources(
    paths: &[PathBuf],
    fuzzlength: i64,
    filters: &Filters,
) -> Result<Vec<ClusterAssembly>> {
    let alignments = load_sources(paths)?;
    let (alignments, _dropped) = filter::apply(alignments, filters);
    let clusters = cluster_alignments(&alignments);

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
    filters: &Filters,
) -> Result<AltSpliceResult> {
    let assemblies = assemble_sources(paths, fuzzlength, filters)?;
    Ok(analyze(&assemblies, fuzzlength))
}

/// Algorithm 2 + the optional ORF/UTR step: analyze, then graft an external
/// gene-prediction CDS onto each isoform and tag events by region. `code`
/// selects the genetic code for divergent-isoform stop-codon detection.
pub fn reconcile_sources(
    paths: &[PathBuf],
    gene_pred: &Path,
    genome: &Path,
    fuzzlength: i64,
    filters: &Filters,
    code: GeneticCode,
) -> Result<(Vec<Isoform>, Vec<Locus>, ReconcileResult)> {
    let asr = analyze_sources(paths, fuzzlength, filters)?;
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
}

/// Build consensus gene models: parse weights, ingest weighted evidence, cluster into
/// per-contig regions, and run the both-strand trellis on each region in parallel.
pub fn consensus_sources(cfg: &ConsensusConfig) -> Result<Vec<CalledGene>> {
    let weights = Weights::parse_file(&cfg.weights)?;
    let chains = load_evidence(
        &cfg.gene_predictions,
        &cfg.protein_alignments,
        &cfg.transcript_alignments,
        &weights,
    )?;
    let genome = Fasta::load(&cfg.genome)?;
    let code = GeneticCode::from_ncbi_id(cfg.genetic_code).map_err(|msg| CombinrError::Parse {
        file: "genetic-code".into(),
        line: 0,
        msg,
    })?;
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
        intergenic_adjust: 1.0,
        strict: cfg.strict,
    };
    let eng = EngineParams {
        max_prev_exons: cfg.max_prev_exons,
        research_size: cfg.research_size,
        research_intergenic: cfg.research_intergenic,
        search_long_introns: cfg.search_long_introns,
        intergenic_adjust: 1.0,
    };

    let mask = match &cfg.repeats {
        Some(path) => crate::consensus::repeats::parse_repeats(path)?,
        None => std::collections::HashMap::new(),
    };
    let empty: Vec<crate::model::Coordset> = Vec::new();

    let regions = build_regions(&chains, cfg.flank);
    // Regions are independent → process in parallel (as `assemble_sources` does).
    let per_region: Vec<Vec<CalledGene>> = regions
        .par_iter()
        .map(|r| {
            let region_mask = mask.get(&r.contig).unwrap_or(&empty);
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
    if cfg.promote_transcript_orfs && !cfg.transcript_alignments.is_empty() {
        let asr = analyze_sources(&cfg.transcript_alignments, 20, &Filters::none())?; // PASA fuzz default
        genes = crate::consensus::promote::promote_and_merge(
            genes,
            &asr,
            &genome,
            &code,
            cfg.min_coding_length,
        );
    }
    Ok(genes)
}

/// `consensus --alt-splice`: build the consensus genes, then attach each locus's
/// alternative transcript isoforms as extra mRNAs (CDS grafted from the consensus). Returns
/// the output genes and the region-tagged alt-splice events.
pub fn consensus_with_isoforms(cfg: &ConsensusConfig) -> Result<(Vec<OutGene>, Vec<EventRecord>)> {
    let genes = consensus_sources(cfg)?;
    let genome = Fasta::load(&cfg.genome)?;
    let code = GeneticCode::from_ncbi_id(cfg.genetic_code).map_err(|msg| CombinrError::Parse {
        file: "genetic-code".into(),
        line: 0,
        msg,
    })?;
    let asr = analyze_sources(&cfg.transcript_alignments, 20, &Filters::none())?; // PASA fuzz default
    Ok(crate::consensus::altsplice::annotate(
        &genes, asr, &genome, &code,
    ))
}
