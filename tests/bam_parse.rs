//! BAM input parsing: the segemehl-mapped `locus.bam` fixture is decoded into
//! the same `Alignment` structures the GTF/GFF3 readers produce. Exon
//! coordinates and strands are checked against values hand-derived from
//! `samtools view` (see the CIGAR oracles in `src/io/bam.rs`).
//!
//! `load_sources` performs no filtering, so all six primary mapped records are
//! returned here; the default `--max-intron` cap is applied later in the
//! pipeline, not at parse time.

mod common;

use combinr::io::load_sources;
use combinr::model::Strand;
use std::path::PathBuf;

fn locus_bam() -> PathBuf {
    common::data_dir().join("locus.bam")
}

#[test]
fn parses_all_primary_records() {
    let aligns = load_sources(&[locus_bam()]).unwrap();
    assert_eq!(
        aligns.len(),
        6,
        "expected 6 alignments, got {}",
        aligns.len()
    );
    assert!(
        aligns.iter().all(|a| a.contig == "contig_6362"),
        "all fixture reads map to contig_6362"
    );
}

#[test]
fn seven_exon_read_coords_and_strands() {
    let aligns = load_sources(&[locus_bam()]).unwrap();
    let a = aligns
        .iter()
        .find(|a| a.acc == "TRINITY_DN4001_c0_g2_i1")
        .expect("g2_i1 present");
    assert_eq!(a.num_segments(), 7);
    // Segments are in ascending genomic order (Alignment::refine sorts by lend).
    let first = a.segments[0].coords;
    let last = a.segments[6].coords;
    assert_eq!((first.lend, first.rend), (331213, 331270));
    assert_eq!((last.lend, last.rend), (332863, 333076));
    // FLAG=16 → aligned Minus; XS:+ → spliced Plus.
    assert_eq!(a.aligned_orient, Strand::Minus);
    assert_eq!(a.spliced_orient, Strand::Plus);
}

#[test]
fn three_exon_read_absorbs_indels() {
    let aligns = load_sources(&[locus_bam()]).unwrap();
    let a = aligns
        .iter()
        .find(|a| a.acc == "TRINITY_DN4001_c0_g4_i1")
        .expect("g4_i1 present");
    // CIGAR 4S23M1D14M1D172M112N26M1I217N72M: the 1D/1I are absorbed; only N splits.
    assert_eq!(a.num_segments(), 3);
    let first = a.segments[0].coords;
    let last = a.segments[2].coords;
    assert_eq!((first.lend, first.rend), (342777, 342987));
    assert_eq!((last.lend, last.rend), (343343, 343414));
    assert_eq!(a.aligned_orient, Strand::Plus); // FLAG=0
    assert_eq!(a.spliced_orient, Strand::Plus); // XS:+
}

#[test]
fn aligned_and_spliced_orient_are_independent() {
    // DN6617...i3 aligns to the + strand (FLAG=0) but its gene is transcribed on
    // the - strand (XS:-): aligned_orient and spliced_orient must differ.
    let aligns = load_sources(&[locus_bam()]).unwrap();
    let a = aligns
        .iter()
        .find(|a| a.acc == "TRINITY_DN6617_c1_g1_i3")
        .expect("DN6617 present");
    assert_eq!(a.aligned_orient, Strand::Plus);
    assert_eq!(a.spliced_orient, Strand::Minus);
}
