//! Alternative-splicing modeling — port of
//! `PerlLib/CDNA/Alternative_splice_comparer.pm`.
//!
//! Assemblies from Algorithm 1 are promoted to [`Isoform`]s, de-duplicated, and
//! grouped into gene loci by transitive exon overlap on the same strand. Within
//! each locus every isoform pair is compared, classifying five event types
//! (retained intron, alternate donor/acceptor, transcription-start/polyA within
//! an intron, exon skipping, alternate terminal exons). All overlap tests here
//! use the **strict** rule (`overlaps_strict`), matching the Perl classifiers.

pub mod events;
pub mod locus;

use crate::assemble::ClusterAssembly;
use crate::model::{Coordset, Strand};
use std::collections::BTreeSet;
use std::sync::Arc;

/// A distinct splice isoform (one assembled transcript structure).
#[derive(Debug, Clone)]
pub struct Isoform {
    pub id: String,
    pub contig: String,
    pub strand: Strand,
    /// Genomic exons, sorted ascending by `lend`.
    pub exons: Vec<Coordset>,
    pub contained_accs: Vec<String>,
    pub source_set: BTreeSet<Arc<str>>,
    pub is_fl: bool,
}

impl Isoform {
    /// Genomic introns: the gaps between consecutive (lend-sorted) exons.
    pub fn introns(&self) -> Vec<Coordset> {
        self.exons
            .windows(2)
            .map(|w| Coordset {
                lend: w[0].rend + 1,
                rend: w[1].lend - 1,
            })
            .collect()
    }

    fn plus_like(&self) -> bool {
        // Unknown is treated as '+' for terminal/labeling purposes.
        self.strand != Strand::Minus
    }

    /// 5'-terminal genomic coordinate (transcription start site).
    pub fn tss_coord(&self) -> i64 {
        if self.plus_like() {
            self.exons.first().unwrap().lend
        } else {
            self.exons.last().unwrap().rend
        }
    }

    /// 3'-terminal genomic coordinate (polyadenylation site).
    pub fn tes_coord(&self) -> i64 {
        if self.plus_like() {
            self.exons.last().unwrap().rend
        } else {
            self.exons.first().unwrap().lend
        }
    }

    /// The exon containing the 5' terminus (in transcript order).
    pub fn first_transcript_exon(&self) -> Coordset {
        if self.plus_like() {
            *self.exons.first().unwrap()
        } else {
            *self.exons.last().unwrap()
        }
    }

    /// The exon containing the 3' terminus (in transcript order).
    pub fn last_transcript_exon(&self) -> Coordset {
        if self.plus_like() {
            *self.exons.last().unwrap()
        } else {
            *self.exons.first().unwrap()
        }
    }
}

/// A gene locus: a set of overlapping, same-strand isoforms (indices into the
/// analysis's isoform table).
#[derive(Debug, Clone)]
pub struct Locus {
    pub id: String,
    pub contig: String,
    pub strand: Strand,
    pub isoform_indices: Vec<usize>,
}

/// The kind of alternative-splicing event detected between two isoforms.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum EventKind {
    RetainedIntron,
    AltAcceptor,
    AltDonor,
    StartWithinIntron,
    EndWithinIntron,
    ExonSkip,
    AlternateExon,
}

impl EventKind {
    pub fn as_str(self) -> &'static str {
        match self {
            EventKind::RetainedIntron => "retained_intron",
            EventKind::AltAcceptor => "alt_acceptor",
            EventKind::AltDonor => "alt_donor",
            EventKind::StartWithinIntron => "start_within_intron",
            EventKind::EndWithinIntron => "end_within_intron",
            EventKind::ExonSkip => "exon_skip",
            EventKind::AlternateExon => "alternate_exon",
        }
    }
}

/// Whether an event falls in the 5' UTR, CDS, or 3' UTR of the reference
/// isoform. Filled only when the ORF/UTR step runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum RegionClass {
    FivePrimeUtr,
    Cds,
    ThreePrimeUtr,
}

impl RegionClass {
    pub fn as_str(self) -> &'static str {
        match self {
            RegionClass::FivePrimeUtr => "5UTR",
            RegionClass::Cds => "CDS",
            RegionClass::ThreePrimeUtr => "3UTR",
        }
    }
}

/// One detected alternative-splicing event. Directional: `isoform_a` is the
/// isoform the feature is found *in* relative to `isoform_b`.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct EventRecord {
    pub kind: EventKind,
    pub contig: String,
    pub strand: Strand,
    pub coords: Vec<Coordset>,
    pub isoform_a: String,
    pub isoform_b: String,
    pub detail: Option<String>,
    /// Region tag (5'UTR / CDS / 3'UTR), set by the ORF/UTR step.
    pub region: Option<RegionClass>,
}

/// Result of an alternative-splicing analysis.
pub struct AltSpliceResult {
    pub isoforms: Vec<Isoform>,
    pub loci: Vec<Locus>,
    pub events: Vec<EventRecord>,
}

/// Analyze a set of assemblies: build isoforms, group into loci, and classify
/// alternative-splicing events between isoforms of each locus.
pub fn analyze(assemblies: &[ClusterAssembly], fuzzlength: i64) -> AltSpliceResult {
    let mut isoforms = locus::build_isoforms(assemblies);
    let loci = locus::group_into_loci(&mut isoforms);

    let mut events = Vec::new();
    for locus in &loci {
        let idxs = &locus.isoform_indices;
        for a in 0..idxs.len() {
            for b in (a + 1)..idxs.len() {
                events.extend(events::compare_isoforms(
                    &isoforms[idxs[a]],
                    &isoforms[idxs[b]],
                    fuzzlength,
                ));
            }
        }
    }

    // Deterministic order + de-duplication of identical records.
    events.sort();
    events.dedup();

    AltSpliceResult {
        isoforms,
        loci,
        events,
    }
}
