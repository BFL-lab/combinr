//! Integration test for the optional ORF/UTR reconciliation (M5): a clean
//! isoform inherits the predicted CDS verbatim, while a retained-intron isoform
//! is projected from the same start codon to a premature stop, and the event is
//! tagged as CDS-affecting.

use combinr::altsplice::{EventKind, RegionClass};
use combinr::model::Coordset;
use combinr::pipeline::reconcile_sources;
use std::io::Write;

// genome: 120 bp of 'C' (Pro, no stops) with ATG@11, premature TAA@41 (intron),
// clean TAA@67 (exon2).
fn genome_fa() -> String {
    let mut g = vec![b'C'; 120];
    let put = |g: &mut Vec<u8>, pos1: usize, s: &[u8]| {
        for (k, &b) in s.iter().enumerate() {
            g[pos1 - 1 + k] = b;
        }
    };
    put(&mut g, 11, b"ATG");
    put(&mut g, 41, b"TAA");
    put(&mut g, 67, b"TAA");
    format!(">chr1\n{}\n", String::from_utf8(g).unwrap())
}

const PRED: &str = "\
##gff-version 3
chr1\tpred\tmRNA\t11\t90\t.\t+\t.\tID=pred1;Parent=g1
chr1\tpred\tCDS\t11\t40\t.\t+\t0\tID=pred1.cds;Parent=pred1
chr1\tpred\tCDS\t61\t69\t.\t+\t0\tID=pred1.cds;Parent=pred1
";

const ISOS: &str = "\
chr1\td\texon\t11\t40\t.\t+\t.\tgene_id \"g\"; transcript_id \"clean\";
chr1\td\texon\t61\t90\t.\t+\t.\tgene_id \"g\"; transcript_id \"clean\";
chr1\td\texon\t11\t90\t.\t+\t.\tgene_id \"g\"; transcript_id \"retained\";
";

fn write(suffix: &str, content: &str) -> tempfile::NamedTempFile {
    let mut f = tempfile::Builder::new().suffix(suffix).tempfile().unwrap();
    f.write_all(content.as_bytes()).unwrap();
    f
}

#[test]
fn clean_inherits_and_retained_projects_premature_stop() {
    let genome = write(".fa", &genome_fa());
    let pred = write(".gff3", PRED);
    let isos = write(".gtf", ISOS);

    let (isoforms, loci, recon) = reconcile_sources(
        &[isos.path().to_path_buf()],
        pred.path(),
        genome.path(),
        20,
        &combinr::filter::Filters::none(),
        combinr::orf::GeneticCode::default(),
    )
    .unwrap();

    assert_eq!(isoforms.len(), 2);
    assert_eq!(loci.len(), 1);

    // locate each isoform by its contained transcript accession
    let find = |acc: &str| {
        isoforms
            .iter()
            .position(|i| i.contained_accs.iter().any(|a| a == acc))
            .unwrap()
    };
    let clean = find("clean");
    let retained = find("retained");

    // clean isoform: inherits the predicted CDS verbatim, not altered.
    let cc = &recon.isoform_codings[clean];
    assert_eq!(cc.len(), 1);
    assert!(!cc[0].coding_altered);
    assert_eq!(
        cc[0].cds_segments,
        vec![Coordset::new(11, 40), Coordset::new(61, 69)]
    );
    assert_eq!(cc[0].three_utr, vec![Coordset::new(70, 90)]);

    // retained-intron isoform: diverges, projects to the premature stop at 41-43.
    let rc = &recon.isoform_codings[retained];
    assert_eq!(rc.len(), 1);
    assert!(rc[0].coding_altered);
    assert_eq!(rc[0].cds_segments, vec![Coordset::new(11, 43)]);

    // the retained-intron event is CDS-affecting.
    let ri = recon
        .events
        .iter()
        .find(|e| e.kind == EventKind::RetainedIntron)
        .expect("a retained-intron event");
    assert_eq!(ri.region, Some(RegionClass::Cds));
}
