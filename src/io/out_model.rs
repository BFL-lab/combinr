//! A format-neutral output model. The GFF3 and GTF writers both consume
//! [`OutGene`], so the gene/mRNA/exon/CDS/UTR structure is built once and
//! rendered in either format.

use crate::altsplice::{Isoform, Locus};
use crate::assemble::ClusterAssembly;
use crate::model::{Coordset, Strand};
use crate::orf::CodingAnnotation;

/// One transcript to emit, with optional CDS/UTR features.
pub struct OutTranscript {
    pub transcript_id: String,
    pub contig: String,
    pub strand: Strand,
    pub exons: Vec<Coordset>,
    pub cds: Vec<Coordset>,
    pub five_utr: Vec<Coordset>,
    pub three_utr: Vec<Coordset>,
    /// Extra attributes (key → values; multiple values render as a comma list).
    pub attrs: Vec<(String, Vec<String>)>,
    /// Phase (0-2) of the 5'-most CDS segment: the bases to skip before the first
    /// complete codon. 0 for every CDS combinr builds; nonzero only for a 5'-partial
    /// model read verbatim from an input gene set (`--models`).
    pub cds_start_phase: u8,
}

/// One gene (locus) with its transcripts.
pub struct OutGene {
    pub gene_id: String,
    pub contig: String,
    pub strand: Strand,
    pub lend: i64,
    pub rend: i64,
    /// Extra gene attributes (key → values), written after `ID=` by the GFF3 writer
    /// (ignored by GTF). Empty for every gene combinr builds; carries an input gene
    /// set's attributes through `--models`.
    pub attrs: Vec<(String, Vec<String>)>,
    pub transcripts: Vec<OutTranscript>,
}

fn span(exons: &[Coordset]) -> (i64, i64) {
    (
        exons.iter().map(|c| c.lend).min().unwrap_or(0),
        exons.iter().map(|c| c.rend).max().unwrap_or(0),
    )
}

/// Build genes from non-redundant assemblies (one gene per assembly).
pub fn from_assemblies(assemblies: &[ClusterAssembly]) -> Vec<OutGene> {
    let mut order: Vec<&ClusterAssembly> = assemblies.iter().collect();
    order.sort_by(|a, b| {
        (
            a.structure.contig.as_str(),
            a.structure.coords.lend,
            a.structure.coords.rend,
            a.structure.acc.as_str(),
        )
            .cmp(&(
                b.structure.contig.as_str(),
                b.structure.coords.lend,
                b.structure.coords.rend,
                b.structure.acc.as_str(),
            ))
    });

    order
        .into_iter()
        .enumerate()
        .map(|(i, asm)| {
            let gene_id = format!("asm{}", i + 1);
            let exons: Vec<Coordset> = asm.structure.segments.iter().map(|s| s.coords).collect();
            let (lend, rend) = span(&exons);
            let mut attrs = Vec::new();
            push_sources(&mut attrs, asm.member_provenance.iter().map(|p| p.as_ref()));
            attrs.push(("contains".into(), asm.contained_accs.clone()));
            attrs.push(("num_contained".into(), vec![asm.num_contained.to_string()]));
            attrs.push(("full_length".into(), vec![asm.is_fl.to_string()]));
            OutGene {
                gene_id: gene_id.clone(),
                contig: asm.structure.contig.clone(),
                strand: asm.orient,
                lend,
                rend,
                attrs: Vec::new(),
                transcripts: vec![OutTranscript {
                    transcript_id: format!("{gene_id}.t1"),
                    contig: asm.structure.contig.clone(),
                    strand: asm.orient,
                    exons,
                    cds: Vec::new(),
                    five_utr: Vec::new(),
                    three_utr: Vec::new(),
                    attrs,
                    cds_start_phase: 0,
                }],
            }
        })
        .collect()
}

/// Build genes from isoform loci (no CDS).
pub fn from_loci(isoforms: &[Isoform], loci: &[Locus]) -> Vec<OutGene> {
    build_loci(isoforms, loci, None)
}

/// Build genes from isoform loci with CDS/UTR annotations.
pub fn from_annotated_loci(
    isoforms: &[Isoform],
    loci: &[Locus],
    codings: &[Vec<CodingAnnotation>],
) -> Vec<OutGene> {
    build_loci(isoforms, loci, Some(codings))
}

fn build_loci(
    isoforms: &[Isoform],
    loci: &[Locus],
    codings: Option<&[Vec<CodingAnnotation>]>,
) -> Vec<OutGene> {
    let mut genes = Vec::new();
    for locus in loci {
        if locus.isoform_indices.is_empty() {
            continue;
        }
        let mut lend = i64::MAX;
        let mut rend = i64::MIN;
        let mut transcripts = Vec::new();

        for &i in &locus.isoform_indices {
            let iso = &isoforms[i];
            let (l, r) = span(&iso.exons);
            lend = lend.min(l);
            rend = rend.max(r);
            let anns: &[CodingAnnotation] = codings.map(|c| c[i].as_slice()).unwrap_or(&[]);

            if anns.is_empty() {
                transcripts.push(make_transcript(iso, &iso.id, None));
            } else {
                let multi = anns.len() > 1;
                for ann in anns {
                    let tid = if multi {
                        format!("{}.{}", iso.id, ann.model_id)
                    } else {
                        iso.id.clone()
                    };
                    transcripts.push(make_transcript(iso, &tid, Some(ann)));
                }
            }
        }

        genes.push(OutGene {
            gene_id: locus.id.clone(),
            contig: locus.contig.clone(),
            strand: locus.strand,
            lend,
            rend,
            attrs: Vec::new(),
            transcripts,
        });
    }
    genes
}

fn make_transcript(iso: &Isoform, tid: &str, ann: Option<&CodingAnnotation>) -> OutTranscript {
    let mut attrs = Vec::new();
    push_sources(&mut attrs, iso.source_set.iter().map(|s| s.as_ref()));
    attrs.push(("contains".into(), iso.contained_accs.clone()));
    attrs.push(("full_length".into(), vec![iso.is_fl.to_string()]));
    if let Some(a) = ann {
        attrs.push(("cds_model".into(), vec![a.model_id.clone()]));
        attrs.push(("coding_altered".into(), vec![a.coding_altered.to_string()]));
        attrs.push(("partial3".into(), vec![a.partial3.to_string()]));
    }
    OutTranscript {
        transcript_id: tid.to_string(),
        contig: iso.contig.clone(),
        strand: iso.strand,
        exons: iso.exons.clone(),
        cds: ann.map(|a| a.cds_segments.clone()).unwrap_or_default(),
        five_utr: ann.map(|a| a.five_utr.clone()).unwrap_or_default(),
        three_utr: ann.map(|a| a.three_utr.clone()).unwrap_or_default(),
        attrs,
        cds_start_phase: 0,
    }
}

fn push_sources<'a>(
    attrs: &mut Vec<(String, Vec<String>)>,
    sources: impl Iterator<Item = &'a str>,
) {
    let v: Vec<String> = sources.map(String::from).collect();
    if !v.is_empty() {
        attrs.push(("sources".into(), v));
    }
}

/// CDS phase per segment, in transcription order, given the phase of the 5'-most
/// segment (`start_phase`, 0 for a CDS that begins on a complete codon). Returns a map
/// keyed by `(lend, rend)`.
pub fn cds_phases(
    cds: &[Coordset],
    strand: Strand,
    start_phase: u8,
) -> std::collections::HashMap<(i64, i64), u8> {
    let mut order: Vec<&Coordset> = cds.iter().collect();
    if strand == Strand::Minus {
        order.sort_by(|a, b| b.lend.cmp(&a.lend));
    } else {
        order.sort_by_key(|c| c.lend);
    }
    let mut phases = std::collections::HashMap::new();
    // Coding bases preceding each segment, counted from the first complete codon: the
    // `start_phase` bases skipped at the 5' end count negatively.
    let mut cumulative = -i64::from(start_phase);
    for seg in order {
        let phase = ((3 - cumulative.rem_euclid(3)) % 3) as u8;
        phases.insert((seg.lend, seg.rend), phase);
        cumulative += seg.len();
    }
    phases
}

/// Strand character for output (`.` for unknown).
pub fn strand_char(s: Strand) -> char {
    match s {
        Strand::Plus => '+',
        Strand::Minus => '-',
        Strand::Unknown => '.',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cs(l: i64, r: i64) -> Coordset {
        Coordset::new(l, r)
    }

    #[test]
    fn cds_phases_start_phase_zero_is_unchanged() {
        // segments of 10, 11, 9 bases, plus strand: 0, then (3 - 10%3)%3 = 2, then
        // (3 - 21%3)%3 = 0.
        let cds = [cs(1, 10), cs(21, 31), cs(41, 49)];
        let p = cds_phases(&cds, Strand::Plus, 0);
        assert_eq!((p[&(1, 10)], p[&(21, 31)], p[&(41, 49)]), (0, 2, 0));
    }

    #[test]
    fn cds_phases_honours_start_phase_on_both_strands() {
        let cds = [cs(1, 10), cs(21, 31), cs(41, 49)];
        // plus, start phase 1: 10-1 = 9 coding bases -> next 0; 9+11 = 20 -> next 1.
        let p = cds_phases(&cds, Strand::Plus, 1);
        assert_eq!((p[&(1, 10)], p[&(21, 31)], p[&(41, 49)]), (1, 0, 1));
        // plus, start phase 2: 10-2 = 8 -> next 1; 8+11 = 19 -> next 2.
        let p = cds_phases(&cds, Strand::Plus, 2);
        assert_eq!((p[&(1, 10)], p[&(21, 31)], p[&(41, 49)]), (2, 1, 2));
        // minus: transcription order is (41,49), (21,31), (1,10).
        // start phase 1: 9-1 = 8 -> next 1; 8+11 = 19 -> next 2.
        let p = cds_phases(&cds, Strand::Minus, 1);
        assert_eq!((p[&(41, 49)], p[&(21, 31)], p[&(1, 10)]), (1, 1, 2));
        // start phase 2: 9-2 = 7 -> next 2; 7+11 = 18 -> next 0.
        let p = cds_phases(&cds, Strand::Minus, 2);
        assert_eq!((p[&(41, 49)], p[&(21, 31)], p[&(1, 10)]), (2, 2, 0));
    }
}
