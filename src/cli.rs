//! Command-line interface (clap).

use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

/// Default pairwise-compatibility fuzz distance (bp) at alignment termini.
/// Matches `CDNA_alignment_assembler` `fuzzlength = 20`.
pub const DEFAULT_FUZZLENGTH: i64 = 20;

/// Default stringent-overlap fraction (percent of the shorter span). 0 = off
/// (any 1-bp overlap clusters), i.e. PASA `--stringent_alignment_overlap` disabled.
pub const DEFAULT_STRINGENT_OVERLAP: f64 = 0.0;

#[derive(Parser, Debug)]
#[command(
    name = "combinr",
    version,
    about = "PASA + EVidenceModeler core algorithms in Rust"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum OutputFormat {
    Gff3,
    Gtf,
}

/// Output/execution options shared by every subcommand (the "Execution and
/// logging" help group). Per-field `help_heading` keeps the heading correct
/// regardless of where this struct is flattened.
#[derive(Parser, Debug)]
pub struct CommonOpts {
    /// Write transcript models here instead of stdout (the default).
    #[arg(short = 'o', long, help_heading = "Execution and logging")]
    pub output: Option<PathBuf>,
    /// Output format for transcript models.
    #[arg(
        long,
        value_enum,
        default_value_t = OutputFormat::Gff3,
        help_heading = "Execution and logging"
    )]
    pub format: OutputFormat,
    /// Worker threads for per-cluster/-region work. Pass 0 to use all available cores.
    #[arg(
        short = 't',
        long,
        default_value_t = 4,
        help_heading = "Execution and logging"
    )]
    pub threads: usize,
}

#[derive(Subcommand, Debug)]
pub enum Command {
    /// Combine transcript sources into a non-redundant assembly set. Add
    /// `--alt-splice` to also classify alternative-splicing events between the
    /// isoforms, `--gene-pred` + `--genome` to reconcile an external CDS onto
    /// the isoforms for CDS + 5'/3' UTRs, or `--models` + `--genome` to append the
    /// alternative isoforms to an existing gene set kept unchanged.
    Assemble(AssembleArgs),

    /// Build consensus gene models by integrating weighted evidence (EVM-style).
    Consensus(ConsensusArgs),

    /// Hidden: assemble from the C++ `pasa` token format on stdin/file and emit
    /// the `pasa`-compatible illustration lines (golden-diff harness).
    #[command(hide = true)]
    AssembleTokens(AssembleTokensArgs),
}

// ---------------------------------------------------------------------------
// assemble — the PASA driver (assembly / alt-splice / CDS-UTR reconcile)
// ---------------------------------------------------------------------------

#[derive(Parser, Debug)]
#[command(next_help_heading = "Inputs")]
pub struct AssembleInputs {
    /// One or more transcript files (GTF, GFF3, and/or BAM); repeatable.
    /// BAM is auto-detected by the `.bam` extension; CRAM and text SAM are not
    /// supported.
    #[arg(short, long, required = true)]
    pub input: Vec<PathBuf>,
}

/// Optional inputs that switch `assemble` from plain non-redundant assembly into
/// the CDS/UTR reconcile step (--gene-pred) or into augmenting an existing gene set
/// (--models). Either needs --genome, and --genome needs one of them; they are not part
/// of the default transcript-combining behavior.
#[derive(Parser, Debug)]
#[command(next_help_heading = "CDS/UTR reconcile & model augmentation (optional)")]
#[command(group = clap::ArgGroup::new("anchor").args(["gene_pred", "models"]))]
pub struct AssembleReconcile {
    /// Gene-prediction GFF3 with mapped CDS. Supplying this (together with
    /// --genome) grafts the prediction's CDS onto each assembled isoform to add
    /// CDS + 5'/3' UTRs, instead of emitting coordinate-only transcripts. Requires
    /// --genome.
    #[arg(long, requires = "genome")]
    pub gene_pred: Option<PathBuf>,
    /// Gene-model GFF3 (gene/mRNA/exon/CDS[/UTR]) to keep unchanged and augment. Every
    /// input gene and mRNA is emitted as-is (structure, IDs, attributes); each assembled
    /// transcript isoform that overlaps a gene on its strand, hosts the gene's CDS start
    /// and is a genuine alternative (a novel splice junction — not a CDS-identical
    /// structure or a mere terminal truncation of any existing mRNA) is appended as an
    /// extra mRNA of that gene, with the gene's CDS grafted on (inherited when the
    /// structure matches, re-projected from the same start codon when it diverges).
    /// Isoforms overlapping no model are dropped. Requires --genome; mutually exclusive
    /// with --gene-pred.
    #[arg(long, requires = "genome", conflicts_with = "gene_pred")]
    pub models: Option<PathBuf>,
    /// Genome FASTA, read for codon/UTR sequence during the reconcile / augmentation
    /// step. Required with, and only used by, --gene-pred or --models.
    #[arg(long, requires = "anchor")]
    pub genome: Option<PathBuf>,
}

#[derive(Parser, Debug)]
#[command(next_help_heading = "Pipeline behavior")]
pub struct AssemblePipeline {
    /// Classify alternative-splicing events between isoforms and write them to
    /// `--events`. Superseded by the reconcile / augmentation step when --gene-pred or
    /// --models is given (which emits region-tagged events anyway). Off by default.
    #[arg(long)]
    pub alt_splice: bool,
    /// Path for the tab-separated alt-splice event report (written by
    /// `--alt-splice` or the --gene-pred / --models step).
    #[arg(long, default_value = "combinr.alt_splice_events.tsv")]
    pub events: PathBuf,
    /// NCBI genetic code for the reconcile / augmentation step's stop-codon detection
    /// (only used with --gene-pred or --models). Supported: NCBI tables 1-6, 9-16, 21-33.
    #[arg(short = 'g', long = "genetic-code", default_value_t = 1)]
    pub genetic_code: u32,
}

#[derive(Parser, Debug)]
#[command(next_help_heading = "Assembly tuning")]
pub struct AssembleTuning {
    /// Pairwise-compatibility fuzz distance (bp) at alignment termini.
    #[arg(long, default_value_t = DEFAULT_FUZZLENGTH)]
    pub fuzzlength: i64,
    /// Require two transcripts' genomic-span overlap to be >= this percent of the
    /// SHORTER span before they cluster into one gene (PASA
    /// --stringent_alignment_overlap). 0 = off (any 1-bp overlap clusters, the default).
    #[arg(long, default_value_t = DEFAULT_STRINGENT_OVERLAP)]
    pub stringent_overlap: f64,
    /// Drop alignments whose average percent identity is below this (off by default).
    #[arg(long)]
    pub min_avg_per_id: Option<f64>,
    /// Drop alignments containing an intron shorter than this (bp) (off by default).
    #[arg(long)]
    pub min_intron: Option<i64>,
    /// Drop alignments containing an intron longer than this (bp). Defaults to
    /// 100000 (matching PASA's MAX_INTRON_LENGTH) to discard spurious long-range
    /// junctions; pass 0 (or a negative value) to disable the cap.
    #[arg(long, default_value = "100000")]
    pub max_intron: Option<i64>,
}

#[derive(Parser, Debug)]
pub struct AssembleArgs {
    #[command(flatten)]
    pub inputs: AssembleInputs,
    #[command(flatten)]
    pub reconcile: AssembleReconcile,
    #[command(flatten)]
    pub pipeline: AssemblePipeline,
    #[command(flatten)]
    pub tuning: AssembleTuning,
    #[command(flatten)]
    pub common: CommonOpts,
}

// ---------------------------------------------------------------------------
// consensus — the EVidenceModeler driver
// ---------------------------------------------------------------------------

#[derive(Parser, Debug)]
#[command(next_help_heading = "Inputs")]
pub struct ConsensusInputs {
    /// Evidence weights file: 3 whitespace columns `CLASS TYPE WEIGHT`, where TYPE is
    /// the GFF column-2 source of the corresponding evidence rows.
    #[arg(long, required = true)]
    pub weights: PathBuf,
    /// Genome FASTA (used by the trellis junction stop-check and CDS/UTR graft).
    #[arg(long, required = true)]
    pub genome: PathBuf,
    /// Gene-prediction GFF3 (consensus reads CDS features); repeatable.
    #[arg(long)]
    pub gene_predictions: Vec<PathBuf>,
    /// Protein-alignment GFF3 (match chains carrying `Target=`); repeatable.
    #[arg(long)]
    pub protein_alignments: Vec<PathBuf>,
    /// Transcript-alignment GFF3 (match chains carrying `Target=`); repeatable.
    #[arg(long)]
    pub transcript_alignments: Vec<PathBuf>,
    /// Repeat-mask GFF3 (optional; masked bases are excluded from scoring).
    #[arg(long)]
    pub repeats: Option<PathBuf>,
}

#[derive(Parser, Debug)]
#[command(next_help_heading = "Consensus tuning")]
pub struct ConsensusTuning {
    /// NCBI genetic code for stop-codon detection. Supported: NCBI tables 1-6, 9-16, 21-33.
    #[arg(short = 'g', long = "genetic-code", default_value_t = 1)]
    pub genetic_code: u32,
    /// Minimum intron length (bp).
    #[arg(long, default_value_t = 20)]
    pub min_intron_length: i64,
    /// DP look-back limit: max previous exons compared per trellis node.
    #[arg(long, default_value_t = 500)]
    pub max_prev_exons: usize,
    /// Coding/noncoding score ratio below which a gene is flagged low-support.
    #[arg(long, default_value_t = 0.75)]
    pub min_score_ratio: f64,
    /// (--alt-splice only) Require two transcript isoforms' genomic-span overlap to be
    /// >= this percent of the SHORTER span before they attach to one gene (PASA
    /// --stringent_alignment_overlap). 0 = off (the default). The EVM consensus region
    /// partitioner is unaffected.
    #[arg(long, default_value_t = DEFAULT_STRINGENT_OVERLAP)]
    pub stringent_overlap: f64,
    /// Extend protein/transcript 3' termini to a genetic-code stop to form terminal exons
    /// (lets evidence-only loci gain a 3' end). Off by default.
    #[arg(long)]
    pub extend_terminal_stop: bool,
    /// Boost intergenic scores near start/stop evidence-density peaks. Off by default.
    #[arg(long)]
    pub peak_augment: bool,
}

#[derive(Parser, Debug)]
#[command(next_help_heading = "Pipeline behavior")]
pub struct ConsensusBehavior {
    /// Padding (bp) added to each side of an evidence locus to form its DP region.
    /// This only widens the per-locus window — giving the trellis room to place
    /// UTRs and start/stop ends beyond the evidence span. It does NOT affect
    /// clustering: loci are grouped by raw evidence overlap *before* the flank is
    /// applied, so a larger flank never merges separate loci into one gene.
    #[arg(long, default_value_t = 10_000)]
    pub flank: i64,
    /// Drop low-support genes (EVM behaviour). Default keeps them and tags
    /// `low_support=true` so a locus is never left blank.
    #[arg(long)]
    pub strict: bool,
    /// Re-search intergenic gaps of at least this size (bp) between called genes for
    /// additional genes (EVM `--re_search_intergenic`). 0 = off.
    #[arg(long, default_value_t = 0)]
    pub research_intergenic: i64,
    /// Re-search introns of at least this length (bp) for nested genes (EVM
    /// `--search_long_introns`). 0 = off.
    #[arg(long, default_value_t = 0)]
    pub search_long_introns: i64,
    /// At loci with only transcript evidence (no consensus CDS), find the longest ORF in
    /// the transcripts and promote it into the output. Off by default.
    #[arg(long)]
    pub promote_transcript_orfs: bool,
    /// Emit each consensus locus's alternative transcript isoforms as additional mRNAs,
    /// with CDS derived from the consensus; also writes a region-tagged events TSV
    /// (`--events`). Off by default.
    #[arg(long)]
    pub alt_splice: bool,
    /// Path for the alt-splice events TSV (used with `--alt-splice`).
    #[arg(long, default_value = "combinr.alt_splice_events.tsv")]
    pub events: PathBuf,
    /// Also write a per-gene evidence report (EVidenceModeler `.evm.out` style): for
    /// every consensus gene, a header line with its GFF3 ID and scores, then one line per
    /// exon and per intron listing the evidence features `{accession;source}` that
    /// support it.
    #[arg(long, value_name = "FILE")]
    pub evidence_report: Option<PathBuf>,
}

#[derive(Parser, Debug)]
pub struct ConsensusArgs {
    #[command(flatten)]
    pub inputs: ConsensusInputs,
    #[command(flatten)]
    pub tuning: ConsensusTuning,
    #[command(flatten)]
    pub behavior: ConsensusBehavior,
    #[command(flatten)]
    pub common: CommonOpts,
}

// ---------------------------------------------------------------------------
// assemble-tokens (hidden golden-diff harness)
// ---------------------------------------------------------------------------

#[derive(Parser, Debug)]
pub struct AssembleTokensArgs {
    /// Token-format input file (`acc,orient,lend-rend,...`); stdin if omitted.
    pub input: Option<PathBuf>,
    /// Pairwise-compatibility fuzz distance (bp) at alignment termini.
    #[arg(long, default_value_t = DEFAULT_FUZZLENGTH)]
    pub fuzzlength: i64,
}
