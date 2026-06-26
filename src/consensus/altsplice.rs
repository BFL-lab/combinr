//! Phase 2: emit a consensus gene's alternative transcript isoforms as extra mRNAs,
//! with each isoform's CDS derived from the consensus ("the alt-splice set is based off
//! the gene-prediction consensus").
//!
//! Each consensus gene becomes a `CdsModel`; the existing `orf::reconcile` grafts that CDS
//! onto the transcript isoforms (inherit if the structure matches, re-project from the
//! consensus start if divergent) and region-tags the alt-splice events. We then assemble
//! one `OutGene` per consensus gene whose mRNAs are the consensus model PLUS every isoform
//! that hosts its CDS (deduping an isoform identical to the consensus). Opt-in via
//! `--alt-splice`; the default one-model-per-locus output is untouched.

use crate::altsplice::{AltSpliceResult, EventRecord, Isoform};
use crate::consensus::engine::CalledGene;
use crate::io::fasta::Fasta;
use crate::io::out_model::{OutGene, OutTranscript};
use crate::model::Strand;
use crate::orf::{CdsModel, CodingAnnotation, GeneticCode, reconcile};

/// Annotate consensus genes with their alternative transcript isoforms. Returns the
/// output genes (consensus mRNA + hosting isoform mRNAs each) and the region-tagged
/// alt-splice events.
pub fn annotate(
    genes: &[CalledGene],
    asr: AltSpliceResult,
    genome: &Fasta,
    code: &GeneticCode,
) -> (Vec<OutGene>, Vec<EventRecord>) {
    // deterministic output order + gene ids
    let mut order: Vec<usize> = (0..genes.len()).collect();
    order.sort_by(|&a, &b| gkey(&genes[a]).cmp(&gkey(&genes[b])));

    // each consensus gene -> a CdsModel id'd by its output rank
    let mut models = Vec::new();
    let mut model_id_for: Vec<Option<String>> = vec![None; genes.len()];
    for (rank, &gi) in order.iter().enumerate() {
        if let Some(m) = to_cds_model(&genes[gi], format!("cons{rank}")) {
            model_id_for[gi] = Some(m.id.clone());
            models.push(m);
        }
    }

    let recon = reconcile(&asr.isoforms, &asr.loci, asr.events, &models, genome, code);

    let mut out_genes = Vec::with_capacity(order.len());
    for (rank, &gi) in order.iter().enumerate() {
        let g = &genes[gi];
        let gene_id = format!("consensus.{}.g{}", g.contig, rank + 1);

        // Partition the consensus gene's hosting isoforms. An isoform that leaves the
        // CDS unchanged (`coding_altered == false`) is the SAME coding model as the
        // bare consensus — it only clothes it in UTRs — so emitting both produces two
        // mRNAs with identical CDS, one with UTRs and one without. Instead we fold the
        // best such isoform's UTRs onto the single consensus mRNA and drop the rest.
        // Only `coding_altered == true` isoforms are genuine alternatives that earn
        // their own mRNA. (The old `iso.exons == g.exons` test only deduped the no-UTR
        // case, since the consensus exons are CDS-only and a UTR-bearing isoform's
        // exons never match them.)
        let mut utr_donor: Option<(&Isoform, &CodingAnnotation)> = None;
        let mut alt_isoforms: Vec<(&Isoform, &CodingAnnotation)> = Vec::new();
        if let Some(mid) = &model_id_for[gi] {
            for (idx, iso) in asr.isoforms.iter().enumerate() {
                let Some(ann) = recon.isoform_codings[idx]
                    .iter()
                    .find(|c| &c.model_id == mid)
                else {
                    continue;
                };
                if ann.coding_altered {
                    alt_isoforms.push((iso, ann));
                } else if iso.exons != g.exons && more_utr(utr_donor, iso) {
                    // skips an isoform whose exons exactly equal the CDS-only
                    // consensus (no UTRs to contribute)
                    utr_donor = Some((iso, ann));
                }
            }
        }

        // The consensus mRNA — dressed in UTRs from the best CDS-identical isoform when
        // the base consensus has none of its own.
        let mut transcripts = vec![consensus_mrna(g, &gene_id, utr_donor)];
        for (iso, ann) in alt_isoforms {
            let n = transcripts.len();
            transcripts.push(iso_mrna(iso, ann, &gene_id, n));
        }

        let (lend, rend) = span_of(&transcripts);
        out_genes.push(OutGene {
            gene_id,
            contig: g.contig.clone(),
            strand: g.orient,
            lend,
            rend,
            transcripts,
        });
    }

    (out_genes, recon.events)
}

/// A consensus gene's CDS as a `CdsModel` for grafting, or `None` if it has no CDS.
fn to_cds_model(g: &CalledGene, id: String) -> Option<CdsModel> {
    let first = g.cds.first()?;
    let last = g.cds.last()?;
    let (start_genomic, stop_genomic) = if g.orient == Strand::Minus {
        (last.rend, first.lend)
    } else {
        (first.lend, last.rend)
    };
    Some(CdsModel {
        id,
        contig: g.contig.clone(),
        strand: g.orient,
        cds_segments: g.cds.clone(),
        start_genomic,
        stop_genomic,
    })
}

/// The consensus model itself as an mRNA (mirrors `consensus::output::to_out_genes`).
///
/// When the base consensus carries no UTRs of its own but a CDS-identical transcript
/// isoform does (`utr_donor`), the consensus adopts that isoform's exon structure and
/// UTRs — so the single consensus mRNA carries them rather than being duplicated by a
/// separate UTR-bearing isoform. The donor's transcript provenance is recorded so the
/// UTR support is traceable.
fn consensus_mrna(
    g: &CalledGene,
    gene_id: &str,
    utr_donor: Option<(&Isoform, &CodingAnnotation)>,
) -> OutTranscript {
    let ratio = if g.promoted {
        "NA".to_string()
    } else if g.support.score_ratio.is_finite() {
        format!("{:.2}", g.support.score_ratio)
    } else {
        "inf".to_string()
    };
    let support = if g.promoted {
        "transcript_orf"
    } else {
        "consensus"
    };
    let mut attrs = vec![
        ("support".into(), vec![support.to_string()]),
        ("score".into(), vec![format!("{:.1}", g.score)]),
        ("score_ratio".into(), vec![ratio]),
        (
            "coding_length".into(),
            vec![g.support.coding_length.to_string()],
        ),
        (
            "low_support".into(),
            vec![g.support.low_support.to_string()],
        ),
        ("partial5".into(), vec![g.partial5.to_string()]),
        ("partial3".into(), vec![g.partial3.to_string()]),
    ];

    // Adopt the donor's UTRs only when the base consensus has none of its own.
    let graft = utr_donor.filter(|_| g.five_utr.is_empty() && g.three_utr.is_empty());
    let (exons, five_utr, three_utr) = match graft {
        Some((iso, ann)) => {
            let sources: Vec<String> = iso.source_set.iter().map(|s| s.to_string()).collect();
            if !sources.is_empty() {
                attrs.push(("sources".into(), sources));
            }
            if !iso.contained_accs.is_empty() {
                attrs.push(("contains".into(), iso.contained_accs.clone()));
            }
            (iso.exons.clone(), ann.five_utr.clone(), ann.three_utr.clone())
        }
        None => (g.exons.clone(), g.five_utr.clone(), g.three_utr.clone()),
    };

    OutTranscript {
        transcript_id: format!("{gene_id}.consensus"),
        contig: g.contig.clone(),
        strand: g.orient,
        exons,
        cds: g.cds.clone(),
        five_utr,
        three_utr,
        attrs,
    }
}

/// Whether `cand` contributes more UTR coverage than the current donor (larger total
/// exon span; the consensus CDS is identical across all candidates, so more exon span
/// means more UTR). Deterministic: ties keep the incumbent.
fn more_utr(current: Option<(&Isoform, &CodingAnnotation)>, cand: &Isoform) -> bool {
    let span = |ex: &[crate::model::Coordset]| -> i64 { ex.iter().map(|c| c.rend - c.lend).sum() };
    match current {
        None => true,
        Some((iso, _)) => span(&cand.exons) > span(&iso.exons),
    }
}

/// A transcript isoform as an alternative mRNA, with the consensus CDS grafted on.
fn iso_mrna(iso: &Isoform, ann: &CodingAnnotation, gene_id: &str, n: usize) -> OutTranscript {
    let mut attrs = vec![
        ("support".into(), vec!["transcript_isoform".to_string()]),
        (
            "coding_altered".into(),
            vec![ann.coding_altered.to_string()],
        ),
        ("partial3".into(), vec![ann.partial3.to_string()]),
    ];
    let sources: Vec<String> = iso.source_set.iter().map(|s| s.to_string()).collect();
    if !sources.is_empty() {
        attrs.push(("sources".into(), sources));
    }
    if !iso.contained_accs.is_empty() {
        attrs.push(("contains".into(), iso.contained_accs.clone()));
    }
    OutTranscript {
        transcript_id: format!("{gene_id}.iso{n}"),
        contig: iso.contig.clone(),
        strand: iso.strand,
        exons: iso.exons.clone(),
        cds: ann.cds_segments.clone(),
        five_utr: ann.five_utr.clone(),
        three_utr: ann.three_utr.clone(),
        attrs,
    }
}

fn gkey(g: &CalledGene) -> (String, i64, i64, char) {
    (
        g.contig.clone(),
        g.exons.first().map(|c| c.lend).unwrap_or(0),
        g.exons.last().map(|c| c.rend).unwrap_or(0),
        g.orient.to_char(),
    )
}

fn span_of(transcripts: &[OutTranscript]) -> (i64, i64) {
    let lend = transcripts
        .iter()
        .flat_map(|t| t.exons.iter())
        .map(|c| c.lend)
        .min()
        .unwrap_or(0);
    let rend = transcripts
        .iter()
        .flat_map(|t| t.exons.iter())
        .map(|c| c.rend)
        .max()
        .unwrap_or(0);
    (lend, rend)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::altsplice::{AltSpliceResult, Locus};
    use crate::consensus::filter::SupportFlags;
    use crate::model::Coordset;
    use std::collections::BTreeSet;
    use std::io::Write;

    fn write_fasta(seq: &[u8]) -> Fasta {
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        writeln!(tmp, ">chr1").unwrap();
        tmp.write_all(seq).unwrap();
        writeln!(tmp).unwrap();
        Fasta::load(tmp.path()).unwrap()
    }

    fn cs(l: i64, r: i64) -> Coordset {
        Coordset::new(l, r)
    }

    fn consensus_gene(exons: &[(i64, i64)]) -> CalledGene {
        let ex: Vec<Coordset> = exons.iter().map(|&(l, r)| cs(l, r)).collect();
        CalledGene {
            contig: "chr1".into(),
            orient: Strand::Plus,
            exons: ex.clone(),
            cds: ex,
            five_utr: vec![],
            three_utr: vec![],
            partial5: false,
            partial3: false,
            score: 100.0,
            support: SupportFlags {
                score_ratio: 5.0,
                coding_length: 60,
                low_support: false,
            },
            promoted: false,
        }
    }

    fn iso(id: &str, segs: &[(i64, i64)]) -> Isoform {
        Isoform {
            id: id.into(),
            contig: "chr1".into(),
            strand: Strand::Plus,
            exons: segs.iter().map(|&(l, r)| cs(l, r)).collect(),
            contained_accs: vec![id.into()],
            source_set: BTreeSet::new(),
            is_fl: false,
        }
    }

    #[test]
    fn to_cds_model_uses_cds_bounds() {
        let g = consensus_gene(&[(10, 40), (60, 90)]);
        let m = to_cds_model(&g, "cons0".into()).unwrap();
        assert_eq!(m.start_genomic, 10);
        assert_eq!(m.stop_genomic, 90);
        assert_eq!(m.cds_segments, vec![cs(10, 40), cs(60, 90)]);
    }

    #[test]
    fn gene_gets_consensus_plus_divergent_isoform_mrna() {
        // all-C genome with an ATG at the consensus start (10); no stop -> divergent
        // isoform projects 3'-partial. iso1 == consensus structure (deduped); iso2 is a
        // retained-intron variant (kept, coding_altered).
        let mut g = vec![b'C'; 100];
        g[9] = b'A';
        g[10] = b'T';
        g[11] = b'G';
        let genome = write_fasta(&g);

        let gene = consensus_gene(&[(10, 40), (60, 90)]);
        let isoforms = vec![iso("t1", &[(10, 40), (60, 90)]), iso("t2", &[(10, 90)])];
        let loci = vec![Locus {
            id: "L1".into(),
            contig: "chr1".into(),
            strand: Strand::Plus,
            isoform_indices: vec![0, 1],
        }];
        let asr = AltSpliceResult {
            isoforms,
            loci,
            events: vec![],
        };

        let (out, _events) = annotate(&[gene], asr, &genome, &GeneticCode::default());
        assert_eq!(out.len(), 1);
        let mrnas = &out[0].transcripts;
        // consensus mRNA + iso2 (iso1 deduped as identical to the consensus)
        assert_eq!(mrnas.len(), 2);
        assert!(mrnas[0].transcript_id.ends_with(".consensus"));
        assert!(
            mrnas[0]
                .attrs
                .iter()
                .any(|(k, v)| k == "support" && v == &vec!["consensus".to_string()])
        );
        assert!(
            mrnas[1]
                .attrs
                .iter()
                .any(|(k, v)| k == "support" && v == &vec!["transcript_isoform".to_string()])
        );
        // the retained-intron isoform diverges from the consensus CDS
        assert!(
            mrnas[1]
                .attrs
                .iter()
                .any(|(k, v)| k == "coding_altered" && v == &vec!["true".to_string()])
        );
    }

    #[test]
    fn gene_without_transcripts_is_single_mrna() {
        let genome = write_fasta(&[b'C'; 100]);
        let gene = consensus_gene(&[(10, 40), (60, 90)]);
        let asr = AltSpliceResult {
            isoforms: vec![],
            loci: vec![],
            events: vec![],
        };
        let (out, _) = annotate(&[gene], asr, &genome, &GeneticCode::default());
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].transcripts.len(), 1);
    }

    #[test]
    fn utr_only_isoform_folds_into_single_consensus_mrna() {
        // The consensus CDS is [(10,40),(60,90)] with no UTRs of its own. A transcript
        // isoform hosts the IDENTICAL CDS (same intron 41..59, contains the model start
        // 10 and stop 90) but extends the flanking exons to add UTRs. It must NOT spawn
        // a second mRNA — the single consensus mRNA adopts its UTRs instead.
        let genome = write_fasta(&[b'C'; 110]);
        let gene = consensus_gene(&[(10, 40), (60, 90)]);
        let isoforms = vec![iso("t1", &[(1, 40), (60, 100)])];
        let loci = vec![Locus {
            id: "L1".into(),
            contig: "chr1".into(),
            strand: Strand::Plus,
            isoform_indices: vec![0],
        }];
        let asr = AltSpliceResult {
            isoforms,
            loci,
            events: vec![],
        };

        let (out, _events) = annotate(&[gene], asr, &genome, &GeneticCode::default());
        assert_eq!(out.len(), 1);
        let mrnas = &out[0].transcripts;
        // Exactly one mRNA: the consensus, now carrying the isoform's UTRs.
        assert_eq!(mrnas.len(), 1, "CDS-identical UTR isoform must not duplicate the consensus");
        let m = &mrnas[0];
        assert!(m.transcript_id.ends_with(".consensus"));
        assert!(
            m.attrs
                .iter()
                .any(|(k, v)| k == "support" && v == &vec!["consensus".to_string()])
        );
        assert!(!m.five_utr.is_empty(), "5' UTR grafted from the isoform");
        assert!(!m.three_utr.is_empty(), "3' UTR grafted from the isoform");
        // The consensus mRNA adopts the isoform's UTR-extended exon structure...
        assert_eq!(m.exons.first().unwrap().lend, 1);
        assert_eq!(m.exons.last().unwrap().rend, 100);
        // ...while keeping the consensus CDS unchanged...
        assert_eq!(m.cds, vec![cs(10, 40), cs(60, 90)]);
        // ...and records the UTR donor's provenance.
        assert!(m.attrs.iter().any(|(k, _)| k == "contains"));
    }
}
