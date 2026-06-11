//! Command-line interface (clap).

use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

/// Default pairwise-compatibility fuzz distance (bp) at alignment termini.
/// Matches `CDNA_alignment_assembler` `fuzzlength = 20`.
pub const DEFAULT_FUZZLENGTH: i64 = 20;

#[derive(Parser, Debug)]
#[command(name = "combinr", version, about = "PASA core algorithms in Rust")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,

    /// Pairwise-compatibility fuzz distance (bp) at alignment termini.
    #[arg(long, global = true, default_value_t = DEFAULT_FUZZLENGTH)]
    pub fuzzlength: i64,

    /// Output format for transcript models.
    #[arg(long, global = true, value_enum, default_value_t = OutputFormat::Gff3)]
    pub format: OutputFormat,

    /// Worker threads for per-cluster assembly (default: all available cores).
    #[arg(short = 't', long, global = true)]
    pub threads: Option<usize>,

    // ---- optional, off-by-default quality filters ----
    /// Drop alignments whose average percent identity is below this.
    #[arg(long, global = true)]
    pub min_avg_per_id: Option<f64>,
    /// Drop alignments containing an intron shorter than this (bp).
    #[arg(long, global = true)]
    pub min_intron: Option<i64>,
    /// Drop alignments containing an intron longer than this (bp).
    #[arg(long, global = true)]
    pub max_intron: Option<i64>,

    #[arg(short, long, global = true, action = clap::ArgAction::Count)]
    pub verbose: u8,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum OutputFormat {
    Gff3,
    Gtf,
}

#[derive(Subcommand, Debug)]
pub enum Command {
    /// Combine transcript sources into a non-redundant assembly set.
    Assemble(AssembleArgs),

    /// Assemble, then classify alternative-splicing events between isoforms.
    Altsplice(AltspliceArgs),

    /// Add CDS/UTR by reconciling an external gene-prediction GFF3.
    Orf(OrfArgs),

    /// Full pipeline: assemble -> altsplice -> (optional) orf.
    Run(RunArgs),

    /// Hidden: assemble from the C++ `pasa` token format on stdin/file and emit
    /// the `pasa`-compatible illustration lines (golden-diff harness).
    #[command(hide = true)]
    AssembleTokens(AssembleTokensArgs),
}

#[derive(Parser, Debug)]
pub struct AssembleArgs {
    /// One or more transcript files (GTF and/or GFF3); repeatable.
    #[arg(short, long, required = true)]
    pub input: Vec<PathBuf>,
}

#[derive(Parser, Debug)]
pub struct AltspliceArgs {
    #[arg(short, long, required = true)]
    pub input: Vec<PathBuf>,
    /// Path for the tab-separated alt-splice event report.
    #[arg(long, default_value = "combinr.alt_splice_events.tsv")]
    pub events: PathBuf,
}

#[derive(Parser, Debug)]
pub struct OrfArgs {
    #[arg(short, long, required = true)]
    pub input: Vec<PathBuf>,
    /// Gene-prediction GFF3 carrying CDS features.
    #[arg(long, required = true)]
    pub gene_pred: PathBuf,
    /// Genome FASTA (required for the ORF/UTR step).
    #[arg(long, required = true)]
    pub genome: PathBuf,
    /// Path for the tab-separated alt-splice event report.
    #[arg(long, default_value = "combinr.alt_splice_events.tsv")]
    pub events: PathBuf,
}

#[derive(Parser, Debug)]
pub struct RunArgs {
    #[arg(short, long, required = true)]
    pub input: Vec<PathBuf>,
    /// Gene-prediction GFF3 with CDS (enables the ORF/UTR step when given).
    #[arg(long)]
    pub gene_pred: Option<PathBuf>,
    /// Genome FASTA (required together with --gene-pred).
    #[arg(long)]
    pub genome: Option<PathBuf>,
    /// Path for the tab-separated alt-splice event report.
    #[arg(long, default_value = "combinr.alt_splice_events.tsv")]
    pub events: PathBuf,
}

#[derive(Parser, Debug)]
pub struct AssembleTokensArgs {
    /// Token-format input file (`acc,orient,lend-rend,...`); stdin if omitted.
    pub input: Option<PathBuf>,
}
