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

fn repeats_file(gff: &str) -> tempfile::NamedTempFile {
    let mut tmp = tempfile::NamedTempFile::new().unwrap();
    std::io::Write::write_all(&mut tmp, gff.as_bytes()).unwrap();
    tmp
}

/// Repeats that intersect no region (an unknown contig, and coordinates past the end of
/// the 63,304 bp `Contig1`, beyond the flank) leave the output unchanged.
#[test]
fn repeats_outside_every_region_leave_golden_unchanged() {
    let tmp = repeats_file(
        "\
NoSuchContig\tRM\tmatch\t1\t100000\t.\t+\t.\t.
Contig1\tRM\tmatch\t90000\t95000\t.\t+\t.\t.
Contig1\tRM\tmatch\t80000\t85000\t.\t+\t.\t.
",
    );
    let mut cfg = default_config();
    cfg.repeats = Some(tmp.path().to_path_buf());
    let got = line_set(&render(&cfg));
    let expected = common::load_lines(&fixture("smalltest.consensus.gff3.golden"));
    assert_eq!(got, expected, "out-of-region repeats changed the output");
}

/// Masking the CDS of the first golden gene changes the output (the mask is applied).
#[test]
fn repeats_masking_first_gene_cds_change_output() {
    let golden = std::fs::read_to_string(fixture("smalltest.consensus.gff3.golden")).unwrap();
    let cds: Vec<Vec<&str>> = golden
        .lines()
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .filter(|c| c.len() == 9 && c[2] == "CDS")
        .collect();
    let first_parent = cds[0][8].rsplit("Parent=").next().unwrap();
    let rows: String = cds
        .iter()
        .filter(|c| c[8].ends_with(&format!("Parent={first_parent}")))
        .map(|c| format!("{}\tRM\tmatch\t{}\t{}\t.\t+\t.\t.\n", c[0], c[3], c[4]))
        .collect();
    assert!(!rows.is_empty(), "golden has a CDS");
    let tmp = repeats_file(&rows);
    let mut cfg = default_config();
    cfg.repeats = Some(tmp.path().to_path_buf());
    let got = line_set(&render(&cfg));
    let expected = common::load_lines(&fixture("smalltest.consensus.gff3.golden"));
    assert_ne!(got, expected, "masking the first gene's CDS had no effect");
}
