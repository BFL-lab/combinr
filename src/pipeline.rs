//! End-to-end orchestration: load sources → (filter) → cluster → assemble each
//! cluster in parallel, and (for alt-splice) group into loci and classify events.

use crate::altsplice::{AltSpliceResult, Isoform, Locus, analyze};
use crate::assemble::{ClusterAssembly, assemble_cluster};
use crate::cluster::cluster_alignments;
use crate::error::Result;
use crate::filter::{self, Filters};
use crate::io::fasta::Fasta;
use crate::io::load_sources;
use crate::orf::{ReconcileResult, parse_cds_models, reconcile};
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
/// gene-prediction CDS onto each isoform and tag events by region.
pub fn reconcile_sources(
    paths: &[PathBuf],
    gene_pred: &Path,
    genome: &Path,
    fuzzlength: i64,
    filters: &Filters,
) -> Result<(Vec<Isoform>, Vec<Locus>, ReconcileResult)> {
    let asr = analyze_sources(paths, fuzzlength, filters)?;
    let models = parse_cds_models(gene_pred)?;
    let genome = Fasta::load(genome)?;
    let recon = reconcile(&asr.isoforms, &asr.loci, asr.events, &models, &genome);
    Ok((asr.isoforms, asr.loci, recon))
}
