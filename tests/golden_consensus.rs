//! Consensus regression golden.
//!
//! Runs `combinr consensus` (via the library) on the EvidenceModeler `testing/` data set
//! and compares its GFF3 against a committed frozen snapshot. Because our consensus
//! deliberately diverges from EVM (non-canonical sites + two independent per-strand
//! trellises), this is a SELF-golden: it guards against regressions in OUR algorithm, not
//! parity with EVM. Input data is vendored from `EVidenceModeler/testing/`; the golden
//! was generated once from this binary. No EVM/Perl runs at test time.

mod common;

use combinr::consensus::to_out_genes;
use combinr::io::writer_gff3;
use combinr::pipeline::{ConsensusConfig, consensus_sources};
use std::collections::BTreeSet;
use std::path::PathBuf;

fn fixture(name: &str) -> PathBuf {
    common::data_dir().join("consensus").join(name)
}

/// Config mirroring the CLI defaults (all opt-in heuristics off).
fn default_config() -> ConsensusConfig {
    ConsensusConfig {
        weights: fixture("weights.txt"),
        gene_predictions: vec![fixture("gene_predictions.gff3")],
        protein_alignments: vec![fixture("protein_alignments.gff3")],
        transcript_alignments: vec![fixture("transcript_alignments.gff3")],
        genome: fixture("genome.fasta"),
        repeats: None,
        genetic_code: 1,
        flank: 10_000,
        strict: false,
        max_prev_exons: 500,
        min_score_ratio: 0.75,
        min_intron_length: 20,
        research_size: 10_000,
        research_intergenic: 0,
        search_long_introns: 0,
        extend_terminal_stop: false,
        peak_augment: false,
        promote_transcript_orfs: false,
        alt_splice: false,
        min_coding_length: 150,
        stringent_overlap: 0.0,
    }
}

fn render(cfg: &ConsensusConfig) -> String {
    let genes = consensus_sources(cfg).expect("consensus_sources");
    let out = to_out_genes(&genes);
    let mut buf = Vec::new();
    writer_gff3::write(&mut buf, &out).expect("write gff3");
    String::from_utf8(buf).expect("utf8")
}

fn line_set(s: &str) -> BTreeSet<String> {
    s.lines()
        .map(|l| l.trim_end().to_string())
        .filter(|l| !l.is_empty())
        .collect()
}

#[test]
fn consensus_matches_committed_golden() {
    let got = line_set(&render(&default_config()));
    let expected = common::load_lines(&fixture("smalltest.consensus.gff3.golden"));
    assert_eq!(
        got, expected,
        "consensus output drifted from the committed golden \
         (regenerate tests/data/consensus/smalltest.consensus.gff3.golden if intended)"
    );
}

#[test]
fn consensus_output_is_deterministic() {
    let cfg = default_config();
    assert_eq!(
        render(&cfg),
        render(&cfg),
        "consensus output is non-deterministic"
    );
}

/// With only transcript evidence and `--promote-transcript-orfs`, every locus is
/// transcript-only, so the output is all promoted (de-novo ORF) genes, each tagged.
#[test]
fn promotion_recovers_transcript_only_loci() {
    let mut cfg = default_config();
    cfg.gene_predictions = Vec::new();
    cfg.protein_alignments = Vec::new();
    cfg.promote_transcript_orfs = true;

    let gff = render(&cfg);
    let genes = gff.lines().filter(|l| l.contains("\tgene\t")).count();
    let promoted = gff.matches("support=transcript_orf").count();
    assert!(genes > 0, "transcript-only input + promotion yields genes");
    assert_eq!(
        genes, promoted,
        "every gene here is a promoted transcript-ORF"
    );
}

/// `--alt-splice` attaches each consensus locus's transcript isoforms as extra mRNAs and
/// emits region-tagged events.
#[test]
fn alt_splice_emits_isoform_mrnas_and_events() {
    use combinr::pipeline::consensus_with_isoforms;

    let mut cfg = default_config();
    cfg.alt_splice = true;
    let (genes, events) = consensus_with_isoforms(&cfg).expect("consensus_with_isoforms");

    let mrnas: usize = genes.iter().map(|g| g.transcripts.len()).sum();
    assert!(
        mrnas > genes.len(),
        "alt-splice adds isoform mRNAs beyond one per gene"
    );
    assert!(!events.is_empty(), "alt-splice emits region-tagged events");
    assert!(
        genes.iter().all(|g| {
            g.transcripts
                .iter()
                .any(|t| t.transcript_id.ends_with(".consensus"))
        }),
        "every gene keeps its consensus mRNA"
    );
}
