use anyhow::{Context, Result};
use clap::Parser;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

mod cli;

use cli::{AssembleArgs, AssembleTokensArgs, Cli, Command, ConsensusArgs, OutputFormat};
use combinr::altsplice::EventRecord;
use combinr::assemble::Assembler;
use combinr::consensus::{CalledGene, report};
use combinr::filter::Filters;
use combinr::io::out_model::{OutGene, from_annotated_loci, from_assemblies, from_loci};
use combinr::io::{writer_events, writer_gff3, writer_gtf};
use combinr::pipeline::{
    ConsensusConfig, DEFAULT_MIN_CODING_LENGTH, DEFAULT_RESEARCH_SIZE, analyze_sources,
    assemble_sources, consensus_sources, consensus_with_isoforms, parse_genetic_code,
    reconcile_sources,
};
use combinr::token::parse_tokens;

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::Assemble(a) => run_assemble(a),
        Command::Consensus(a) => run_consensus(a),
        Command::AssembleTokens(a) => run_assemble_tokens(a),
    }
}

/// Configure the global rayon pool. A positive `threads` caps the pool at that
/// many workers (default 4); `0` leaves rayon's own default (all available cores).
fn init_threads(threads: usize) {
    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .ok();
    }
}

/// `assemble`: the PASA driver. Bare emits a non-redundant assembly set; `--alt-splice`
/// also groups into loci and classifies alternative-splicing events; `--gene-pred`
/// with `--genome` reconciles an external CDS onto the isoforms for CDS + 5'/3' UTRs
/// (the former `orf`/`run` step).
fn run_assemble(a: AssembleArgs) -> Result<()> {
    init_threads(a.common.threads);
    let fmt = a.common.format;
    let output = a.common.output.clone();
    let fuzz = a.tuning.fuzzlength;
    let overlap = a.tuning.stringent_overlap;
    let filters = Filters {
        min_avg_per_id: a.tuning.min_avg_per_id,
        min_intron: a.tuning.min_intron,
        // --max-intron defaults to 100000; 0 or negative disables the cap.
        max_intron: a.tuning.max_intron.filter(|&n| n > 0),
    };

    match (a.reconcile.gene_pred, a.reconcile.genome) {
        // CDS/UTR reconcile: graft an external prediction's CDS onto the isoforms.
        // Supersedes --alt-splice; the reconcile path emits region-tagged events too.
        (Some(gene_pred), Some(genome)) => {
            let code = parse_genetic_code(a.pipeline.genetic_code)?;
            let (isoforms, loci, recon) = reconcile_sources(
                &a.inputs.input,
                &gene_pred,
                &genome,
                fuzz,
                overlap,
                &filters,
                code,
            )
            .with_context(|| "reconciling ORF/UTR")?;
            let genes = from_annotated_loci(&isoforms, &loci, &recon.isoform_codings);
            write_models(&genes, fmt, output.as_deref())?;
            write_events_file(&recon.events, &a.pipeline.events)?;
            let coding = recon
                .isoform_codings
                .iter()
                .filter(|c| !c.is_empty())
                .count();
            eprintln!(
                "combinr: {} isoform(s) in {} loci, {coding} coding, {} event(s) -> {}",
                isoforms.len(),
                loci.len(),
                recon.events.len(),
                a.pipeline.events.display()
            );
            Ok(())
        }
        // Alt-splice classification only.
        (None, None) if a.pipeline.alt_splice => {
            let r = analyze_sources(&a.inputs.input, fuzz, overlap, &filters)
                .with_context(|| "analyzing alt-splicing")?;
            let genes = from_loci(&r.isoforms, &r.loci);
            write_models(&genes, fmt, output.as_deref())?;
            write_events_file(&r.events, &a.pipeline.events)?;
            eprintln!(
                "combinr: {} isoform(s) in {} loci, {} alt-splice event(s) -> {}",
                r.isoforms.len(),
                r.loci.len(),
                r.events.len(),
                a.pipeline.events.display()
            );
            Ok(())
        }
        // Bare assembly: non-redundant set, no events.
        (None, None) => {
            let assemblies = assemble_sources(&a.inputs.input, fuzz, overlap, &filters)
                .with_context(|| "assembling input sources")?;
            let genes = from_assemblies(&assemblies);
            write_models(&genes, fmt, output.as_deref())?;
            eprintln!(
                "combinr: {} non-redundant assemblies from {} source file(s)",
                assemblies.len(),
                a.inputs.input.len()
            );
            Ok(())
        }
        // --gene-pred and --genome are paired by clap `requires`, so a lone one
        // never reaches here.
        _ => unreachable!("--gene-pred and --genome are paired by clap `requires`"),
    }
}

impl ConsensusArgs {
    /// Lower the parsed CLI args into the library's [`ConsensusConfig`]. The two
    /// non-CLI knobs come from named defaults rather than inline literals.
    fn into_config(self) -> ConsensusConfig {
        ConsensusConfig {
            weights: self.inputs.weights,
            gene_predictions: self.inputs.gene_predictions,
            protein_alignments: self.inputs.protein_alignments,
            transcript_alignments: self.inputs.transcript_alignments,
            genome: self.inputs.genome,
            repeats: self.inputs.repeats,
            genetic_code: self.tuning.genetic_code,
            flank: self.behavior.flank,
            strict: self.behavior.strict,
            max_prev_exons: self.tuning.max_prev_exons,
            min_score_ratio: self.tuning.min_score_ratio,
            min_intron_length: self.tuning.min_intron_length,
            research_size: DEFAULT_RESEARCH_SIZE,
            research_intergenic: self.behavior.research_intergenic,
            search_long_introns: self.behavior.search_long_introns,
            extend_terminal_stop: self.tuning.extend_terminal_stop,
            peak_augment: self.tuning.peak_augment,
            promote_transcript_orfs: self.behavior.promote_transcript_orfs,
            alt_splice: self.behavior.alt_splice,
            min_coding_length: DEFAULT_MIN_CODING_LENGTH,
            stringent_overlap: self.tuning.stringent_overlap,
        }
    }
}

/// `consensus`: build EVM-style consensus gene models by integrating weighted evidence
/// across both strands, then emit them as GFF3 (or GTF). Low-support genes are flagged,
/// not dropped, unless `--strict` is given.
fn run_consensus(a: ConsensusArgs) -> Result<()> {
    use combinr::consensus::{ordered_with_ids, to_out_genes};
    use combinr::model::Strand;

    init_threads(a.common.threads);
    let fmt = a.common.format;
    let output = a.common.output.clone();
    let strict = a.behavior.strict;
    let alt_splice = a.behavior.alt_splice;
    let events_path = a.behavior.events.clone();
    let report_path = a.behavior.evidence_report.clone();
    let cfg = a.into_config();
    let report_note = |p: &Option<std::path::PathBuf>| {
        p.as_ref()
            .map(|p| format!("; evidence report -> {}", p.display()))
            .unwrap_or_default()
    };

    // --alt-splice: emit consensus + transcript-isoform mRNAs and a region-tagged events TSV.
    if alt_splice {
        let r = consensus_with_isoforms(&cfg)
            .with_context(|| "building consensus alt-splice models")?;
        write_models(&r.out_genes, fmt, output.as_deref())?;
        write_events_file(&r.events, &events_path)?;
        if let Some(path) = &report_path {
            write_report_file(&ordered_with_ids(&r.genes), path)?;
        }
        let mrnas: usize = r.out_genes.iter().map(|g| g.transcripts.len()).sum();
        eprintln!(
            "combinr consensus (--alt-splice): {} gene(s), {mrnas} mRNA(s), {} event(s) -> {}{}",
            r.out_genes.len(),
            r.events.len(),
            events_path.display(),
            report_note(&report_path)
        );
        return Ok(());
    }

    let genes = consensus_sources(&cfg).with_context(|| "building consensus gene models")?;
    let plus = genes.iter().filter(|g| g.orient == Strand::Plus).count();
    let low = genes.iter().filter(|g| g.support.low_support).count();
    let promoted = genes.iter().filter(|g| g.promoted).count();

    let out_genes = to_out_genes(&genes);
    write_models(&out_genes, fmt, output.as_deref())?;
    if let Some(path) = &report_path {
        write_report_file(&ordered_with_ids(&genes), path)?;
    }

    eprintln!(
        "combinr consensus: {} gene(s) ({} +, {} -); {} flagged low_support{}; {promoted} promoted transcript-ORF{}",
        genes.len(),
        plus,
        genes.len() - plus,
        low,
        if strict {
            " (dropped, --strict)"
        } else {
            " (kept)"
        },
        report_note(&report_path)
    );
    Ok(())
}

/// Write transcript models in the selected format to `output` (a file) or, when
/// `None`, to stdout.
fn write_models(genes: &[OutGene], fmt: OutputFormat, output: Option<&Path>) -> Result<()> {
    match output {
        Some(path) => {
            let file =
                File::create(path).with_context(|| format!("creating {}", path.display()))?;
            let mut out = BufWriter::new(file);
            write_models_to(&mut out, genes, fmt)?;
            out.flush()?;
        }
        None => {
            let stdout = std::io::stdout();
            let mut out = BufWriter::new(stdout.lock());
            write_models_to(&mut out, genes, fmt)?;
            out.flush()?;
        }
    }
    Ok(())
}

/// Serialize `genes` to `out` in the selected format.
fn write_models_to<W: Write>(out: &mut W, genes: &[OutGene], fmt: OutputFormat) -> Result<()> {
    match fmt {
        OutputFormat::Gff3 => writer_gff3::write(out, genes)?,
        OutputFormat::Gtf => writer_gtf::write(out, genes)?,
    }
    Ok(())
}

/// Write the alt-splice event TSV to `path`.
fn write_events_file(events: &[EventRecord], path: &std::path::Path) -> Result<()> {
    let mut ev =
        BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
    writer_events::write_events(&mut ev, events)?;
    ev.flush()?;
    Ok(())
}

/// Write the consensus evidence report (`(gene_id, gene)` in output order) to `path`.
fn write_report_file(genes: &[(String, &CalledGene)], path: &Path) -> Result<()> {
    let mut w =
        BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
    report::write_report(&mut w, genes)?;
    w.flush()?;
    Ok(())
}

/// Hidden golden-diff harness: read the C++ `pasa` token format, assemble, and
/// print the `pasa`-compatible assembly lines.
fn run_assemble_tokens(a: AssembleTokensArgs) -> Result<()> {
    let (text, source) = match a.input {
        Some(path) => {
            let s = std::fs::read_to_string(&path)
                .with_context(|| format!("reading {}", path.display()))?;
            (s, path.display().to_string())
        }
        None => {
            let mut s = String::new();
            std::io::stdin()
                .read_to_string(&mut s)
                .context("reading stdin")?;
            (s, "<stdin>".to_string())
        }
    };

    let alignments = parse_tokens(&text, &source)?;
    let mut assembler = Assembler::new(alignments);
    assembler.set_fuzzlength(a.fuzzlength);
    assembler.assemble()?;
    print!("{}", assembler.format_pasa_assemblies());
    Ok(())
}
