//! Consensus regression golden.
//!
//! Runs `combinr consensus` (via the library) on the EvidenceModeler `testing/` data set
//! and compares its GFF3 against a committed frozen snapshot. Because our consensus
//! deliberately diverges from EVM (non-canonical sites + two independent per-strand
//! trellises), this is a SELF-golden: it guards against regressions in OUR algorithm, not
//! parity with EVM. Input data is vendored from `EVidenceModeler/testing/`; the golden
//! was generated once from this binary. No EVM/Perl runs at test time.

mod common;

use combinr::consensus::report::write_report;
use combinr::consensus::{CalledGene, ordered_with_ids, to_out_genes};
use combinr::io::out_model::OutGene;
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
    let r = consensus_with_isoforms(&cfg).expect("consensus_with_isoforms");
    let (genes, events) = (r.out_genes, r.events);

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

fn render_report(genes: &[CalledGene]) -> String {
    let mut buf = Vec::new();
    write_report(&mut buf, &ordered_with_ids(genes)).expect("write report");
    String::from_utf8(buf).expect("utf8")
}

/// Check an evidence report against the output genes it describes: one header per gene
/// with identical IDs (same order), every exon line's `(min, max)` equals an exon or CDS
/// segment of that gene, and every exon/intron line carries >= 1 `{acc;source}` token.
fn check_report(report: &str, out: &[OutGene]) {
    let mut blocks: Vec<(String, Vec<&str>)> = Vec::new();
    for line in report.lines() {
        if line.starts_with("##") || line.is_empty() {
            continue;
        }
        if let Some(h) = line.strip_prefix("# ") {
            blocks.push((h.split(' ').next().unwrap().to_string(), Vec::new()));
        } else {
            blocks
                .last_mut()
                .expect("feature before header")
                .1
                .push(line);
        }
    }
    let report_ids: Vec<&str> = blocks.iter().map(|(id, _)| id.as_str()).collect();
    let gff_ids: Vec<&str> = out.iter().map(|g| g.gene_id.as_str()).collect();
    assert_eq!(report_ids, gff_ids, "report gene IDs == GFF3 gene IDs");

    let mut exon_lines = 0;
    for ((id, lines), g) in blocks.iter().zip(out) {
        assert!(!lines.is_empty(), "{id}: gene has feature lines");
        let segs: BTreeSet<(i64, i64)> = g
            .transcripts
            .iter()
            .flat_map(|t| t.exons.iter().chain(t.cds.iter()))
            .map(|c| (c.lend, c.rend))
            .collect();
        for l in lines {
            let c: Vec<&str> = l.split('\t').collect();
            assert_eq!(c.len(), 6, "{id}: six columns: {l:?}");
            let tokens: Vec<&str> = c[5].split(',').collect();
            assert!(
                !c[5].is_empty()
                    && tokens
                        .iter()
                        .all(|t| t.starts_with('{') && t.ends_with('}') && t.contains(';')),
                "{id}: evidence tokens present: {l:?}"
            );
            if c[2] == "INTRON" {
                continue;
            }
            exon_lines += 1;
            let (a, b): (i64, i64) = (c[0].parse().unwrap(), c[1].parse().unwrap());
            assert!(
                segs.contains(&(a.min(b), a.max(b))),
                "{id}: exon line {l:?} matches an exon/CDS of the gene"
            );
        }
    }
    assert!(exon_lines > 0);
}

/// The evidence report lists every GFF3 gene under the same ID, with evidence-backed
/// exon/intron lines matching the gene structure.
#[test]
fn evidence_report_matches_gff3_genes() {
    let genes = consensus_sources(&default_config()).expect("consensus_sources");
    let out = to_out_genes(&genes);
    let report = render_report(&genes);
    assert_eq!(
        report.lines().filter(|l| l.starts_with("# ")).count(),
        out.len()
    );
    check_report(&report, &out);
}

/// Same guarantees on the `--alt-splice` path (IDs are the `.g{n}` gene IDs there too).
#[test]
fn evidence_report_matches_alt_splice_genes() {
    use combinr::pipeline::consensus_with_isoforms;

    let mut cfg = default_config();
    cfg.alt_splice = true;
    let r = consensus_with_isoforms(&cfg).expect("consensus_with_isoforms");
    let report = render_report(&r.genes);
    assert_eq!(
        report.lines().filter(|l| l.starts_with("# ")).count(),
        r.out_genes.len()
    );
    check_report(&report, &r.out_genes);
}

/// Promoted transcript-ORF genes: CDS-segment exon rows attributed to their transcripts.
#[test]
fn evidence_report_covers_promoted_genes() {
    let mut cfg = default_config();
    cfg.gene_predictions = Vec::new();
    cfg.protein_alignments = Vec::new();
    cfg.promote_transcript_orfs = true;
    let genes = consensus_sources(&cfg).expect("consensus_sources");
    let out = to_out_genes(&genes);
    let report = render_report(&genes);
    assert!(report.contains("support(transcript_orf)"));
    check_report(&report, &out);
}
