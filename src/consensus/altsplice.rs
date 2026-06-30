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
use crate::model::{Coordset, Strand, introns_between};
use crate::orf::{CdsModel, CodingAnnotation, GeneticCode, reconcile};
use std::collections::BTreeSet;

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

    // Precomputed once: each isoform's introns, reused by the per-gene UTR salvage.
    let iso_introns: Vec<Vec<Coordset>> = asr.isoforms.iter().map(|iso| iso.introns()).collect();

    let mut out_genes = Vec::with_capacity(order.len());
    for (rank, &gi) in order.iter().enumerate() {
        let g = &genes[gi];
        let gene_id = format!("consensus.{}.g{}", g.contig, rank + 1);
        // The consensus CDS introns, computed once per gene and reused by the
        // truncation test and UTR salvage below.
        let cds_introns = seg_introns(&g.cds);

        // An isoform earns its OWN mRNA only when it is a genuine alternative: it
        // introduces a novel splice junction (`coding_altered` AND not a mere terminal
        // truncation; see `is_terminal_truncation`). Everything else — the CDS-identical
        // isoforms and the no-novel-junction terminal truncations — is the same coding
        // model the consensus already represents, so it spawns no second mRNA; its only
        // contribution is UTR, folded onto the single consensus mRNA by `salvage_utrs`.
        // (Dropping the spurious truncation mRNAs mirrors PASA subsuming a contained
        // alignment and EVM never minting a terminal exon at a truncation point.)
        let mut alt_isoforms: Vec<(&Isoform, &CodingAnnotation)> = Vec::new();
        if let Some(mid) = &model_id_for[gi] {
            for (idx, iso) in asr.isoforms.iter().enumerate() {
                let Some(ann) = recon.isoform_codings[idx]
                    .iter()
                    .find(|c| &c.model_id == mid)
                else {
                    continue;
                };
                if ann.coding_altered
                    && !is_terminal_truncation(
                        &iso.exons,
                        &iso_introns[idx],
                        &g.cds,
                        &cds_introns,
                        ann,
                    )
                {
                    alt_isoforms.push((iso, ann));
                }
            }
        }

        // Assemble the consensus mRNA's 5'/3' UTRs from every transcript compatible with
        // its CDS — taking each end from whichever transcript extends it furthest, even a
        // transcript truncated at the other end (a 3'-truncated fragment still yields a
        // valid 5'UTR; a 5'-truncated one a valid 3'UTR).
        let utrs = salvage_utrs(g, &asr.isoforms, &iso_introns, &cds_introns);
        let mut transcripts = vec![consensus_mrna(g, &gene_id, &utrs)];
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
/// When the base consensus carries no UTRs of its own, it adopts the UTRs assembled
/// from compatible transcript evidence (`utrs`, from [`salvage_utrs`]): its CDS-only
/// exon structure is extended with the salvaged 5'/3' UTR segments, and the donor
/// transcripts' provenance is recorded so the UTR support is traceable.
fn consensus_mrna(g: &CalledGene, gene_id: &str, utrs: &Utrs) -> OutTranscript {
    // `support` first (always present here), then the six shared consensus attributes.
    let support = if g.promoted {
        "transcript_orf"
    } else {
        "consensus"
    };
    let mut attrs = vec![("support".into(), vec![support.to_string()])];
    attrs.extend(crate::consensus::output::consensus_core_attrs(g));

    // Adopt the salvaged UTRs only when the base consensus has none of its own; extend
    // the CDS-only exon structure with the UTR segments (merging the UTR that abuts a
    // terminal CDS exon into it, keeping any spliced UTR exon separate).
    let use_salvage = g.five_utr.is_empty() && g.three_utr.is_empty() && !utrs.is_empty();
    let (exons, five_utr, three_utr) = if use_salvage {
        if !utrs.sources.is_empty() {
            attrs.push(("sources".into(), utrs.sources.clone()));
        }
        if !utrs.contains.is_empty() {
            attrs.push(("contains".into(), utrs.contains.clone()));
        }
        let mut all = g.cds.clone();
        all.extend_from_slice(&utrs.five);
        all.extend_from_slice(&utrs.three);
        (merge_coords(all), utrs.five.clone(), utrs.three.clone())
    } else {
        (g.exons.clone(), g.five_utr.clone(), g.three_utr.clone())
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

/// UTRs assembled for a consensus mRNA from compatible transcript evidence.
#[derive(Default)]
struct Utrs {
    five: Vec<Coordset>,
    three: Vec<Coordset>,
    sources: Vec<String>,
    contains: Vec<String>,
}

impl Utrs {
    fn is_empty(&self) -> bool {
        self.five.is_empty() && self.three.is_empty()
    }
}

/// Introns (gaps) between consecutive lend-sorted segments.
fn seg_introns(segs: &[Coordset]) -> Vec<Coordset> {
    introns_between(segs).collect()
}

/// The introns lying wholly within `[lo, hi]`.
fn within(introns: &[Coordset], lo: i64, hi: i64) -> Vec<Coordset> {
    introns
        .iter()
        .copied()
        .filter(|i| i.lend >= lo && i.rend <= hi)
        .collect()
}

/// Whether `iso_exons` shares EVERY splice junction with the consensus CDS over the
/// genomic span the two overlap — i.e. it introduces no novel junction (no alternate
/// donor/acceptor, retained intron or exon skip) there. Differences are confined to the
/// terminal exons, so it is the same spliced structure, just possibly truncated and/or
/// UTR-extended. `iso_introns`/`cds_introns` are precomputed [`seg_introns`].
fn shares_junctions(
    iso_exons: &[Coordset],
    iso_introns: &[Coordset],
    cds: &[Coordset],
    cds_introns: &[Coordset],
) -> bool {
    let (Some(i0), Some(il)) = (iso_exons.first(), iso_exons.last()) else {
        return false;
    };
    let (Some(c0), Some(cl)) = (cds.first(), cds.last()) else {
        return false;
    };
    // The genomic span the iso and the CDS share; no shared span → no shared structure.
    // (first.lend <= last.rend for lend-sorted exons, so `new` never swaps.)
    let iso_span = Coordset::new(i0.lend, il.rend);
    let cds_span = Coordset::new(c0.lend, cl.rend);
    let Some(ov) = iso_span.intersect(&cds_span) else {
        return false;
    };
    within(iso_introns, ov.lend, ov.rend) == within(cds_introns, ov.lend, ov.rend)
}

/// `exons` material strictly beyond `bound`: with `above`, the parts right of `bound`
/// (clipped to start at `bound+1`); otherwise the parts left of `bound` (clipped to end
/// at `bound-1`). Carves a UTR out of a transcript relative to the start/stop codon.
fn clip_beyond(exons: &[Coordset], bound: i64, above: bool) -> Vec<Coordset> {
    let mut out = Vec::new();
    for e in exons {
        if above {
            if e.rend > bound {
                out.push(Coordset {
                    lend: e.lend.max(bound + 1),
                    rend: e.rend,
                });
            }
        } else if e.lend < bound {
            out.push(Coordset {
                lend: e.lend,
                rend: e.rend.min(bound - 1),
            });
        }
    }
    out
}

/// Whether any exon spans the genomic position `pos`.
fn covers(exons: &[Coordset], pos: i64) -> bool {
    exons.iter().any(|e| e.lend <= pos && pos <= e.rend)
}

/// Merge overlapping or directly adjacent segments (sorted, coalesced).
fn merge_coords(mut segs: Vec<Coordset>) -> Vec<Coordset> {
    segs.sort_by_key(|c| (c.lend, c.rend));
    let mut out: Vec<Coordset> = Vec::new();
    for s in segs {
        match out.last_mut() {
            Some(last) if s.lend <= last.rend + 1 => {
                if s.rend > last.rend {
                    last.rend = s.rend;
                }
            }
            _ => out.push(s),
        }
    }
    out
}

fn total_len(segs: &[Coordset]) -> i64 {
    segs.iter().map(|c| c.rend - c.lend + 1).sum()
}

/// Assemble the consensus gene's UTRs from the transcript isoforms compatible with its
/// CDS. A transcript donates a 5'UTR when it spans the start codon (its material beyond
/// the start) and a 3'UTR when it spans the stop codon (its material beyond the stop);
/// each end takes the longest such donor independently, so a 3'-truncated fragment can
/// still supply the 5'UTR and a 5'-truncated one the 3'UTR (the latter never reaches the
/// start codon, so the start-anchored [`graft`] cannot see it — hence this separate
/// scan over all isoforms). Only transcripts sharing every junction with the CDS
/// contribute; genuine alternatives are excluded.
fn salvage_utrs(
    g: &CalledGene,
    isoforms: &[Isoform],
    iso_introns: &[Vec<Coordset>],
    cds_introns: &[Coordset],
) -> Utrs {
    let mut out = Utrs::default();
    let (Some(cds0), Some(cdsl)) = (g.cds.first(), g.cds.last()) else {
        return out;
    };
    let minus = g.orient == Strand::Minus;
    // start/stop genomic = the 5'/3'-most CDS base for the strand.
    let (start_g, stop_g) = if minus {
        (cdsl.rend, cds0.lend)
    } else {
        (cds0.lend, cdsl.rend)
    };
    let (mut best_five, mut best_three) = (0i64, 0i64);
    let mut five_donor: Option<&Isoform> = None;
    let mut three_donor: Option<&Isoform> = None;
    for (i, iso) in isoforms.iter().enumerate() {
        if iso.strand != g.orient {
            continue;
        }
        let (Some(i0), Some(il)) = (iso.exons.first(), iso.exons.last()) else {
            continue;
        };
        if il.rend < cds0.lend || i0.lend > cdsl.rend {
            continue; // no genomic overlap with the CDS
        }
        if !shares_junctions(&iso.exons, &iso_introns[i], &g.cds, cds_introns) {
            continue;
        }
        if covers(&iso.exons, start_g) {
            let five = clip_beyond(&iso.exons, start_g, minus);
            let l = total_len(&five);
            if l > best_five {
                best_five = l;
                out.five = five;
                five_donor = Some(iso);
            }
        }
        if covers(&iso.exons, stop_g) {
            let three = clip_beyond(&iso.exons, stop_g, !minus);
            let l = total_len(&three);
            if l > best_three {
                best_three = l;
                out.three = three;
                three_donor = Some(iso);
            }
        }
    }
    // Provenance from the (up to two) winning donors.
    let mut sources: BTreeSet<String> = BTreeSet::new();
    for d in [five_donor, three_donor].into_iter().flatten() {
        for s in &d.source_set {
            sources.insert(s.to_string());
        }
        for a in &d.contained_accs {
            if !out.contains.contains(a) {
                out.contains.push(a.clone());
            }
        }
    }
    out.sources = sources.into_iter().collect();
    out
}

/// Whether `iso` is a pure terminal truncation of consensus gene `g` rather than a
/// genuine alternative isoform: the grafted ORF ran off an incomplete end
/// (`ann.partial3`) and the isoform introduces no novel splice junction
/// ([`shares_junctions`]), so it is the same coding model, merely incomplete —
/// `coding_altered` was set only because it lost terminal exons and the stop codon, so
/// it must not earn its own mRNA. (Its UTR, if any, is still salvaged onto the consensus
/// by [`salvage_utrs`].) Mirrors PASA subsuming a contained alignment and EVM never
/// minting a terminal exon at a truncation point; a retained intron / alternate
/// donor-acceptor / exon skip perturbs the shared-span junctions and is kept.
fn is_terminal_truncation(
    iso_exons: &[Coordset],
    iso_introns: &[Coordset],
    cds: &[Coordset],
    cds_introns: &[Coordset],
    ann: &CodingAnnotation,
) -> bool {
    ann.partial3 && shares_junctions(iso_exons, iso_introns, cds, cds_introns)
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
    let (lend, rend) = g.span();
    (g.contig.clone(), lend, rend, g.orient.to_char())
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

    #[test]
    fn truncated_isoform_is_dropped_not_emitted() {
        // A 3'-truncated transcript FRAGMENT: it shares the consensus's 5' exons and
        // introns but is MISSING the terminal exon, so its ORF runs off the end with no
        // stop (partial3) and `graft` marks it coding_altered. It introduces NO novel
        // junction, so it is an incomplete copy of the consensus, not a genuine
        // alternative — it must be dropped, leaving a single consensus mRNA. (PASA
        // subsumes such a contained fragment; EVM never creates a terminal exon at the
        // truncation point. Regression for the GEEHB_01841 spurious-isoform bug.)
        let mut g = vec![b'C'; 150];
        g[9] = b'A';
        g[10] = b'T';
        g[11] = b'G'; // start codon at the consensus start (genomic 10)
        let genome = write_fasta(&g);

        let gene = consensus_gene(&[(10, 40), (60, 90), (110, 140)]);
        // fragment: top two exons only; the terminal (110,140) exon is lost.
        let isoforms = vec![iso("frag", &[(10, 40), (60, 90)])];
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
        assert_eq!(
            mrnas.len(),
            1,
            "a no-novel-junction 3'-truncated fragment must not earn its own mRNA"
        );
        assert!(mrnas[0].transcript_id.ends_with(".consensus"));
    }

    #[test]
    fn single_exon_truncated_isoform_is_dropped() {
        // The single-exon analogue (147 of the 246 GEEHB run cases): a single-exon
        // consensus and a single-exon transcript that covers only part of it and runs
        // off the end (partial3). Empty intron chains are trivially equal, so the
        // truncation test fires and the fragment is dropped.
        let mut g = vec![b'C'; 160];
        g[9] = b'A';
        g[10] = b'T';
        g[11] = b'G';
        let genome = write_fasta(&g);

        let gene = consensus_gene(&[(10, 140)]);
        let isoforms = vec![iso("frag", &[(10, 90)])];
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
        assert_eq!(
            out[0].transcripts.len(),
            1,
            "a single-exon 3'-truncated fragment must not earn its own mRNA"
        );
    }

    #[test]
    fn truncated_fragment_donates_its_five_utr_to_the_consensus() {
        // A 3'-truncated fragment (dropped as its own mRNA) still carries a valid 5'UTR.
        // That UTR must be folded onto the single consensus mRNA — not lost.
        let mut g = vec![b'C'; 120];
        g[19] = b'A';
        g[20] = b'T';
        g[21] = b'G'; // start codon at genomic 20
        let genome = write_fasta(&g);

        let gene = consensus_gene(&[(20, 40), (60, 90)]);
        // fragment shares the 5' CDS exon, extends it upstream to 10 (5'UTR 10..19), and
        // is truncated at the 3' end (missing the (60,90) CDS exon) -> partial3, dropped.
        let isoforms = vec![iso("frag", &[(10, 40)])];
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

        let (out, _e) = annotate(&[gene], asr, &genome, &GeneticCode::default());
        assert_eq!(out.len(), 1);
        let mrnas = &out[0].transcripts;
        assert_eq!(mrnas.len(), 1, "the fragment must not spawn a second mRNA");
        let m = &mrnas[0];
        assert!(m.transcript_id.ends_with(".consensus"));
        assert_eq!(m.cds, vec![cs(20, 40), cs(60, 90)], "CDS unchanged");
        assert_eq!(m.five_utr, vec![cs(10, 19)], "fragment's 5'UTR salvaged");
        assert!(m.three_utr.is_empty(), "no bogus 3'UTR from the truncated end");
        assert_eq!(m.exons.first().unwrap().lend, 10, "5'-terminal exon extended");
    }

    #[test]
    fn five_prime_truncated_transcript_donates_its_three_utr() {
        // A 5'-truncated transcript never reaches the start codon, so `graft` cannot host
        // it; the standalone UTR scan must still recover its 3'UTR onto the consensus.
        // This is the case the start-anchored machinery alone would miss.
        let genome = write_fasta(&[b'C'; 130]);

        let gene = consensus_gene(&[(20, 40), (60, 90)]); // stop codon at genomic 90 (+)
        // covers the 3' CDS exon and extends past the stop to 110 (3'UTR 91..110), but
        // never reaches the start codon at 20.
        let isoforms = vec![iso("threep", &[(60, 110)])];
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

        let (out, _e) = annotate(&[gene], asr, &genome, &GeneticCode::default());
        assert_eq!(out.len(), 1);
        let mrnas = &out[0].transcripts;
        assert_eq!(mrnas.len(), 1, "a 5'-truncated UTR donor is not its own mRNA");
        let m = &mrnas[0];
        assert_eq!(m.cds, vec![cs(20, 40), cs(60, 90)], "CDS unchanged");
        assert_eq!(
            m.three_utr,
            vec![cs(91, 110)],
            "3'UTR recovered from a non-grafting transcript"
        );
        assert!(m.five_utr.is_empty());
        assert_eq!(m.exons.last().unwrap().rend, 110, "3'-terminal exon extended");
    }
}
