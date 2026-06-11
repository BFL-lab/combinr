use anyhow::{Context, Result};
use clap::Parser;
use std::fs::File;
use std::io::{BufWriter, Read, Write};

mod cli;

use cli::{
    AltspliceArgs, AssembleArgs, AssembleTokensArgs, Cli, Command, OrfArgs, OutputFormat, RunArgs,
};
use combinr::altsplice::EventRecord;
use combinr::assemble::Assembler;
use combinr::filter::Filters;
use combinr::io::out_model::{OutGene, from_annotated_loci, from_assemblies, from_loci};
use combinr::io::{writer_events, writer_gff3, writer_gtf};
use combinr::pipeline::{analyze_sources, assemble_sources, reconcile_sources};
use combinr::token::parse_tokens;

fn main() -> Result<()> {
    let args = Cli::parse();

    if let Some(t) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .ok();
    }
    let filters = Filters {
        min_avg_per_id: args.min_avg_per_id,
        min_intron: args.min_intron,
        // --max-intron defaults to 100000; 0 or negative disables the cap.
        max_intron: args.max_intron.filter(|&n| n > 0),
    };
    let fuzz = args.fuzzlength;
    let fmt = args.format;

    match args.command {
        Command::Assemble(a) => run_assemble(a, fuzz, &filters, fmt),
        Command::Altsplice(a) => run_altsplice(a, fuzz, &filters, fmt),
        Command::Orf(a) => run_orf(a, fuzz, &filters, fmt),
        Command::Run(a) => run_run(a, fuzz, &filters, fmt),
        Command::AssembleTokens(a) => run_assemble_tokens(a, fuzz),
    }
}

/// `assemble`: load GTF/GFF3 sources, cluster, assemble, write models.
fn run_assemble(a: AssembleArgs, fuzz: i64, filters: &Filters, fmt: OutputFormat) -> Result<()> {
    let assemblies =
        assemble_sources(&a.input, fuzz, filters).with_context(|| "assembling input sources")?;
    let genes = from_assemblies(&assemblies);
    write_models(&genes, fmt)?;
    eprintln!(
        "combinr: {} non-redundant assemblies from {} source file(s)",
        assemblies.len(),
        a.input.len()
    );
    Ok(())
}

/// `altsplice`: assemble + group into loci + classify events.
fn run_altsplice(a: AltspliceArgs, fuzz: i64, filters: &Filters, fmt: OutputFormat) -> Result<()> {
    let r = analyze_sources(&a.input, fuzz, filters).with_context(|| "analyzing alt-splicing")?;
    let genes = from_loci(&r.isoforms, &r.loci);
    write_models(&genes, fmt)?;
    write_events_file(&r.events, &a.events)?;
    eprintln!(
        "combinr: {} isoform(s) in {} loci, {} alt-splice event(s) -> {}",
        r.isoforms.len(),
        r.loci.len(),
        r.events.len(),
        a.events.display()
    );
    Ok(())
}

/// `orf`: alt-splice + reconcile an external CDS into CDS/UTR.
fn run_orf(a: OrfArgs, fuzz: i64, filters: &Filters, fmt: OutputFormat) -> Result<()> {
    let (isoforms, loci, recon) =
        reconcile_sources(&a.input, &a.gene_pred, &a.genome, fuzz, filters)
            .with_context(|| "reconciling ORF/UTR")?;
    let genes = from_annotated_loci(&isoforms, &loci, &recon.isoform_codings);
    write_models(&genes, fmt)?;
    write_events_file(&recon.events, &a.events)?;
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
        a.events.display()
    );
    Ok(())
}

/// `run`: full pipeline; reconciles ORF/UTR when both --gene-pred and --genome
/// are given, otherwise stops after alt-splice.
fn run_run(a: RunArgs, fuzz: i64, filters: &Filters, fmt: OutputFormat) -> Result<()> {
    match (a.gene_pred, a.genome) {
        (Some(gp), Some(g)) => run_orf(
            OrfArgs {
                input: a.input,
                gene_pred: gp,
                genome: g,
                events: a.events,
            },
            fuzz,
            filters,
            fmt,
        ),
        (None, None) => run_altsplice(
            AltspliceArgs {
                input: a.input,
                events: a.events,
            },
            fuzz,
            filters,
            fmt,
        ),
        _ => anyhow::bail!("--gene-pred and --genome must be provided together"),
    }
}

/// Write transcript models to stdout in the selected format.
fn write_models(genes: &[OutGene], fmt: OutputFormat) -> Result<()> {
    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    match fmt {
        OutputFormat::Gff3 => writer_gff3::write(&mut out, genes)?,
        OutputFormat::Gtf => writer_gtf::write(&mut out, genes)?,
    }
    out.flush()?;
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

/// Hidden golden-diff harness: read the C++ `pasa` token format, assemble, and
/// print the `pasa`-compatible assembly lines.
fn run_assemble_tokens(a: AssembleTokensArgs, fuzz: i64) -> Result<()> {
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
    assembler.set_fuzzlength(fuzz);
    assembler.assemble()?;
    print!("{}", assembler.format_pasa_assemblies());
    Ok(())
}
