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

    // ---- pre-assembly quality filters ----
    /// Drop alignments whose average percent identity is below this (off by default).
    #[arg(long, global = true)]
    pub min_avg_per_id: Option<f64>,
    /// Drop alignments containing an intron shorter than this (bp) (off by default).
    #[arg(long, global = true)]
    pub min_intron: Option<i64>,
    /// Drop alignments containing an intron longer than this (bp). Defaults to
    /// 100000 (matching PASA's MAX_INTRON_LENGTH) to discard spurious long-range
    /// junctions; pass 0 (or a negative value) to disable the cap.
    #[arg(long, global = true, default_value = "100000")]
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

    /// Build consensus gene models by integrating weighted evidence (EVM-style).
    Consensus(ConsensusArgs),

    /// Full pipeline: assemble -> altsplice -> (optional) orf.
    Run(RunArgs),

    /// Hidden: assemble from the C++ `pasa` token format on stdin/file and emit
    /// the `pasa`-compatible illustration lines (golden-diff harness).
    #[command(hide = true)]
    AssembleTokens(AssembleTokensArgs),
}

#[derive(Parser, Debug)]
pub struct AssembleArgs {
    /// One or more transcript files (GTF, GFF3, and/or BAM); repeatable.
    /// BAM is auto-detected by the `.bam` extension; CRAM and text SAM are not
    /// supported.
    #[arg(short, long, required = true)]
    pub input: Vec<PathBuf>,
}

#[derive(Parser, Debug)]
pub struct AltspliceArgs {
    /// One or more transcript files (GTF, GFF3, and/or BAM); repeatable.
    #[arg(short, long, required = true)]
    pub input: Vec<PathBuf>,
    /// Path for the tab-separated alt-splice event report.
    #[arg(long, default_value = "combinr.alt_splice_events.tsv")]
    pub events: PathBuf,
}

#[derive(Parser, Debug)]
pub struct OrfArgs {
    /// One or more transcript files (GTF, GFF3, and/or BAM); repeatable.
    #[arg(short, long, required = true)]
    pub input: Vec<PathBuf>,
    /// Gene-prediction GFF3 carrying CDS features.
    #[arg(long, required = true)]
    pub gene_pred: PathBuf,
    /// Genome FASTA (required for the ORF/UTR step).
    #[arg(long, required = true)]
    pub genome: PathBuf,
    /// NCBI genetic code (translation table) for the ORF/UTR step's stop-codon
    /// detection. Default 1 (standard). Supported: 1, 4, 6, 10, 12, 26.
    #[arg(short = 'g', long = "genetic-code", default_value_t = 1)]
    pub genetic_code: u32,
    /// Path for the tab-separated alt-splice event report.
    #[arg(long, default_value = "combinr.alt_splice_events.tsv")]
    pub events: PathBuf,
}

#[derive(Parser, Debug)]
pub struct ConsensusArgs {
    /// Evidence weights file: 3 whitespace columns `CLASS TYPE WEIGHT`, where TYPE is
    /// the GFF column-2 source of the corresponding evidence rows.
    #[arg(long, required = true)]
    pub weights: PathBuf,
    /// Gene-prediction GFF3 (consensus reads CDS features); repeatable.
    #[arg(long)]
    pub gene_predictions: Vec<PathBuf>,
    /// Protein-alignment GFF3 (match chains carrying `Target=`); repeatable.
    #[arg(long)]
    pub protein_alignments: Vec<PathBuf>,
    /// Transcript-alignment GFF3 (match chains carrying `Target=`); repeatable.
    #[arg(long)]
    pub transcript_alignments: Vec<PathBuf>,
    /// Genome FASTA (used by the trellis junction stop-check and CDS/UTR graft).
    #[arg(long, required = true)]
    pub genome: PathBuf,
    /// Repeat-mask GFF3 (optional; masked bases are excluded from scoring).
    #[arg(long)]
    pub repeats: Option<PathBuf>,
    /// NCBI genetic code for stop-codon detection. Supported: 1, 4, 6, 10, 12, 26.
    #[arg(short = 'g', long = "genetic-code", default_value_t = 1)]
    pub genetic_code: u32,
    /// Flank (bp) added to each evidence locus to form a region.
    #[arg(long, default_value_t = 10_000)]
    pub flank: i64,
    /// Drop low-support genes (EVM behaviour). Default keeps them and tags
    /// `low_support=true` so a locus is never left blank.
    #[arg(long)]
    pub strict: bool,
    /// DP look-back limit: max previous exons compared per trellis node.
    #[arg(long, default_value_t = 500)]
    pub max_prev_exons: usize,
    /// Coding/noncoding score ratio below which a gene is flagged low-support.
    #[arg(long, default_value_t = 0.75)]
    pub min_score_ratio: f64,
    /// Minimum intron length (bp).
    #[arg(long, default_value_t = 20)]
    pub min_intron_length: i64,
    /// Re-search intergenic gaps of at least this size (bp) between called genes for
    /// additional genes (EVM `--re_search_intergenic`). 0 = off.
    #[arg(long, default_value_t = 0)]
    pub research_intergenic: i64,
    /// Re-search introns of at least this length (bp) for nested genes (EVM
    /// `--search_long_introns`). 0 = off.
    #[arg(long, default_value_t = 0)]
    pub search_long_introns: i64,
    /// Extend protein/transcript 3' termini to a genetic-code stop to form terminal exons
    /// (lets evidence-only loci gain a 3' end). Off by default.
    #[arg(long)]
    pub extend_terminal_stop: bool,
    /// Boost intergenic scores near start/stop evidence-density peaks. Off by default.
    #[arg(long)]
    pub peak_augment: bool,
}

#[derive(Parser, Debug)]
pub struct RunArgs {
    /// One or more transcript files (GTF, GFF3, and/or BAM); repeatable.
    #[arg(short, long, required = true)]
    pub input: Vec<PathBuf>,
    /// Gene-prediction GFF3 with CDS (enables the ORF/UTR step when given).
    #[arg(long)]
    pub gene_pred: Option<PathBuf>,
    /// Genome FASTA (required together with --gene-pred).
    #[arg(long)]
    pub genome: Option<PathBuf>,
    /// NCBI genetic code (translation table) for the ORF/UTR step's stop-codon
    /// detection. Default 1 (standard). Supported: 1, 4, 6, 10, 12, 26.
    #[arg(short = 'g', long = "genetic-code", default_value_t = 1)]
    pub genetic_code: u32,
    /// Path for the tab-separated alt-splice event report.
    #[arg(long, default_value = "combinr.alt_splice_events.tsv")]
    pub events: PathBuf,
}

#[derive(Parser, Debug)]
pub struct AssembleTokensArgs {
    /// Token-format input file (`acc,orient,lend-rend,...`); stdin if omitted.
    pub input: Option<PathBuf>,
}
