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
}

/// One gene (locus) with its transcripts.
pub struct OutGene {
    pub gene_id: String,
    pub contig: String,
    pub strand: Strand,
    pub lend: i64,
    pub rend: i64,
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
                transcripts: vec![OutTranscript {
                    transcript_id: format!("{gene_id}.t1"),
                    contig: asm.structure.contig.clone(),
                    strand: asm.orient,
                    exons,
                    cds: Vec::new(),
                    five_utr: Vec::new(),
                    three_utr: Vec::new(),
                    attrs,
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

/// CDS phase per segment, in transcription order. Returns a map keyed by
/// `(lend, rend)`.
pub fn cds_phases(cds: &[Coordset], strand: Strand) -> std::collections::HashMap<(i64, i64), u8> {
    let mut order: Vec<&Coordset> = cds.iter().collect();
    if strand == Strand::Minus {
        order.sort_by(|a, b| b.lend.cmp(&a.lend));
    } else {
        order.sort_by_key(|c| c.lend);
    }
    let mut phases = std::collections::HashMap::new();
    let mut cumulative = 0i64;
    for seg in order {
        let phase = ((3 - (cumulative % 3)) % 3) as u8;
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
