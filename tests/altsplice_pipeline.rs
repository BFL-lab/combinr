//! Integration test for Algorithm 2 (M4): the full analyze path over a locus
//! with three isoforms exhibiting an exon skip and an alternate acceptor.

use combinr::altsplice::EventKind;
use combinr::pipeline::analyze_sources;
use std::io::Write;

const GTF: &str = "\
chr1\td\texon\t1000\t1200\t.\t+\t.\tgene_id \"g\"; transcript_id \"iso1\";
chr1\td\texon\t1400\t1600\t.\t+\t.\tgene_id \"g\"; transcript_id \"iso1\";
chr1\td\texon\t1800\t2000\t.\t+\t.\tgene_id \"g\"; transcript_id \"iso1\";
chr1\td\texon\t1000\t1200\t.\t+\t.\tgene_id \"g\"; transcript_id \"iso2\";
chr1\td\texon\t1800\t2000\t.\t+\t.\tgene_id \"g\"; transcript_id \"iso2\";
chr1\td\texon\t1000\t1200\t.\t+\t.\tgene_id \"g\"; transcript_id \"iso3\";
chr1\td\texon\t1420\t1600\t.\t+\t.\tgene_id \"g\"; transcript_id \"iso3\";
chr1\td\texon\t1800\t2000\t.\t+\t.\tgene_id \"g\"; transcript_id \"iso3\";
";

#[test]
fn analyzes_locus_with_skip_and_alt_acceptor() {
    let mut tmp = tempfile::Builder::new().suffix(".gtf").tempfile().unwrap();
    tmp.write_all(GTF.as_bytes()).unwrap();

    let result = analyze_sources(
        &[tmp.path().to_path_buf()],
        20,
        &combinr::filter::Filters::none(),
    )
    .unwrap();

    // Three distinct isoforms collapse into a single locus.
    assert_eq!(result.isoforms.len(), 3);
    assert_eq!(result.loci.len(), 1);
    assert_eq!(result.loci[0].isoform_indices.len(), 3);

    let kinds: Vec<EventKind> = result.events.iter().map(|e| e.kind).collect();
    // The skipping isoform vs each spliced isoform → two exon-skip events.
    assert_eq!(
        kinds.iter().filter(|k| **k == EventKind::ExonSkip).count(),
        2,
        "expected two exon-skip events, got {kinds:?}"
    );
    // The two middle-exon variants differ at the acceptor (left boundary, +).
    assert!(
        kinds.contains(&EventKind::AltAcceptor),
        "expected an alt-acceptor event, got {kinds:?}"
    );
    // No spurious donor call in this construction.
    assert!(!kinds.contains(&EventKind::AltDonor));
}
