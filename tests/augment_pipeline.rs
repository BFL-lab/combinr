//! `assemble --models` end to end: a consensus gene set written to GFF3, re-read and
//! augmented with the transcript isoforms of the vendored EVM test data.
//!
//! (1) every input gene/mRNA is re-emitted byte-identically (bar a widened gene span) and
//!     the appended mRNAs are all `support=transcript_isoform`;
//! (2) the appended isoforms are exactly the `.iso*` mRNAs `consensus --alt-splice`
//!     produces for the same genes (same assembly parameters, same predicate).

mod common;

use combinr::consensus::to_out_genes;
use combinr::filter::Filters;
use combinr::io::out_model::OutGene;
use combinr::io::writer_gff3;
use combinr::model::Coordset;
use combinr::orf::GeneticCode;
use combinr::pipeline::{
    AugmentStats, ConsensusConfig, DEFAULT_FUZZLENGTH, augment_sources, consensus_sources,
    consensus_with_isoforms,
};
use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;
use std::path::PathBuf;

fn fixture(name: &str) -> PathBuf {
    common::data_dir().join("consensus").join(name)
}

/// Same as `tests/golden_consensus.rs::default_config`.
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

fn render(genes: &[OutGene]) -> String {
    let mut buf = Vec::new();
    writer_gff3::write(&mut buf, genes).unwrap();
    String::from_utf8(buf).unwrap()
}

/// GFF3 text split into per-gene blocks of lines (gene row first).
fn blocks(gff3: &str) -> Vec<Vec<&str>> {
    let mut out: Vec<Vec<&str>> = Vec::new();
    for line in gff3.lines().filter(|l| !l.starts_with('#')) {
        if line.split('\t').nth(2) == Some("gene") {
            out.push(Vec::new());
        }
        out.last_mut().unwrap().push(line);
    }
    out
}

fn col(line: &str, i: usize) -> &str {
    line.split('\t').nth(i).unwrap()
}

/// Consensus genes → GFF3 tempfile → `augment_sources` with the consensus path's
/// assembly parameters. Returns (input GFF3 text, augmented genes, counts).
fn augmented() -> (String, Vec<OutGene>, AugmentStats) {
    let cfg = default_config();
    let input = render(&to_out_genes(&consensus_sources(&cfg).unwrap()));
    let mut tmp = tempfile::Builder::new().suffix(".gff3").tempfile().unwrap();
    tmp.write_all(input.as_bytes()).unwrap();
    let (genes, _events, stats) = augment_sources(
        &cfg.transcript_alignments,
        tmp.path(),
        &cfg.genome,
        DEFAULT_FUZZLENGTH,
        cfg.stringent_overlap,
        &Filters::none(),
        GeneticCode::from_ncbi_id(1).unwrap(),
    )
    .unwrap();
    (input, genes, stats)
}

fn is_appended(line: &str) -> bool {
    col(line, 2) == "mRNA" && col(line, 8).contains(";support=transcript_isoform")
}

#[test]
fn input_genes_are_re_emitted_verbatim_with_isoforms_appended() {
    let (input, genes, stats) = augmented();
    let output = render(&genes);
    let (inb, outb) = (blocks(&input), blocks(&output));
    assert!(!inb.is_empty());
    assert_eq!(inb.len(), outb.len(), "one output gene per input gene");

    let mut appended = 0;
    for (i, o) in inb.iter().zip(&outb) {
        // gene row: same id/attrs/contig/strand; span only ever widened
        assert_eq!(col(i[0], 8), col(o[0], 8));
        assert_eq!((col(i[0], 0), col(i[0], 6)), (col(o[0], 0), col(o[0], 6)));
        let n = |l: &str, c| col(l, c).parse::<i64>().unwrap();
        assert!(n(o[0], 3) <= n(i[0], 3) && n(o[0], 4) >= n(i[0], 4));
        // every original mRNA line and child row, byte-identical, in order
        assert_eq!(&o[1..i.len()], &i[1..], "gene {}", col(i[0], 8));
        // everything after is appended isoforms
        let mrnas: Vec<&&str> = o[i.len()..]
            .iter()
            .filter(|l| col(l, 2) == "mRNA")
            .collect();
        assert!(mrnas.iter().all(|l| is_appended(l)), "{mrnas:?}");
        appended += mrnas.len();
    }
    assert!(appended > 0, "the test data yields alternative isoforms");
    // the reported counts are exact: input genes/mRNAs, and every mRNA beyond the input's
    let n_mrna = |s: &str| {
        blocks(s)
            .iter()
            .flatten()
            .filter(|l| col(l, 2) == "mRNA")
            .count()
    };
    let input_mrnas = n_mrna(&input);
    let output_mrnas = n_mrna(&output);
    assert_eq!(stats.input_genes, inb.len());
    assert_eq!(stats.input_mrnas, input_mrnas);
    assert_eq!(stats.appended, output_mrnas - input_mrnas);
    assert_eq!(stats.appended, appended);
}

type Structures = BTreeMap<String, BTreeSet<Vec<Coordset>>>;

/// gene id → the exon structures of its `transcript_isoform` mRNAs.
fn iso_structures(genes: &[OutGene]) -> Structures {
    genes
        .iter()
        .map(|g| {
            let set = g
                .transcripts
                .iter()
                .filter(|t| {
                    t.attrs
                        .iter()
                        .any(|(k, v)| k == "support" && v[0] == "transcript_isoform")
                })
                .map(|t| t.exons.clone())
                .collect();
            (g.gene_id.clone(), set)
        })
        .collect()
}

#[test]
fn appended_isoforms_match_consensus_alt_splice() {
    let (_, genes, _) = augmented();
    let mut cfg = default_config();
    cfg.alt_splice = true;
    let alt = consensus_with_isoforms(&cfg).unwrap();
    let (a, b) = (iso_structures(&genes), iso_structures(&alt.out_genes));
    assert!(a.values().any(|s| !s.is_empty()));
    assert_eq!(a, b);
}
