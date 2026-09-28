//! Phase 2: emit a consensus gene's — or an externally supplied gene set's — alternative
//! transcript isoforms as extra mRNAs, with each isoform's CDS derived from the gene's
//! own CDS ("the alt-splice set is based off the gene-prediction consensus").
//!
//! Each gene's CDS becomes a `CdsModel`; the existing `orf::reconcile` grafts that CDS
//! onto the transcript isoforms (inherit if the structure matches, re-project from the
//! gene's start codon if divergent) and region-tags the alt-splice events. We then
//! assemble one `OutGene` per gene whose mRNAs are its own model(s) PLUS every isoform
//! that hosts its CDS and is a genuine alternative ([`is_genuine_alternative`]).
//!
//! - [`annotate`] (`consensus --alt-splice`): the consensus genes; the consensus mRNA
//!   also adopts UTRs salvaged from compatible transcripts. Opt-in; the default
//!   one-model-per-locus output is untouched.
//! - [`augment`] (`assemble --models`): an input gene set, emitted verbatim, with the
//!   genuine alternatives appended (no UTR salvage, input transcripts untouched).

use crate::altsplice::{AltSpliceResult, EventRecord, Isoform};
use crate::consensus::engine::CalledGene;
use crate::consensus::output::ordered_with_ids;
use crate::io::fasta::Fasta;
use crate::io::out_model::{OutGene, OutTranscript};
use crate::model::{Coordset, Strand, introns_between, merge_coords};
use crate::orf::{CdsModel, CodingAnnotation, GeneticCode, reconcile};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

/// Annotate consensus genes with their alternative transcript isoforms. Returns the
/// output genes (consensus mRNA + hosting isoform mRNAs each) and the region-tagged
/// alt-splice events.
pub fn annotate(
    genes: &[CalledGene],
    asr: AltSpliceResult,
    genome: &Fasta,
    code: &GeneticCode,
) -> (Vec<OutGene>, Vec<EventRecord>) {
    // deterministic output order + gene ids (shared with `to_out_genes` and the report)
    let order = ordered_with_ids(genes);

    // each consensus gene -> a CdsModel id'd by its output rank
    let mut models = Vec::new();
    let mut model_id_for: Vec<Option<String>> = vec![None; order.len()];
    for (rank, (_, g)) in order.iter().enumerate() {
        if let Some(m) = to_cds_model(g, format!("cons{rank}")) {
            model_id_for[rank] = Some(m.id.clone());
            models.push(m);
        }
    }

    let recon = reconcile(&asr.isoforms, &asr.loci, asr.events, &models, genome, code);

    // Precomputed once: each isoform's introns, reused by the per-gene UTR salvage.
    let iso_introns: Vec<Vec<Coordset>> = asr.isoforms.iter().map(|iso| iso.introns()).collect();

    let mut out_genes = Vec::with_capacity(order.len());
    for (rank, (gene_id, g)) in order.into_iter().enumerate() {
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
        if let Some(mid) = &model_id_for[rank] {
            for (idx, iso) in asr.isoforms.iter().enumerate() {
                let Some(ann) = recon.isoform_codings[idx]
                    .iter()
                    .find(|c| &c.model_id == mid)
                else {
                    continue;
                };
                if is_genuine_alternative(&iso.exons, &iso_introns[idx], &g.cds, &cds_introns, ann)
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
            attrs: Vec::new(),
            transcripts,
        });
    }

    (out_genes, recon.events)
}

/// Augment an input gene set (`assemble --models`) with its alternative transcript
/// isoforms. Every input gene and transcript is kept verbatim, in input order; each
/// isoform that hosts the gene's primary CDS (the transcript with the longest total CDS,
/// ties → first) and is a genuine alternative ([`is_genuine_alternative`]) against EVERY
/// model of the gene it hosts is appended as `{gene_id}.iso{n}` (the first `n` whose id
/// is not already taken in that gene) carrying the primary model's grafted CDS. The gene
/// span is only ever widened. As on the consensus path, an isoform may attach to several
/// overlapping same-strand genes. Returns the genes and the region-tagged events.
pub fn augment(
    mut models: Vec<OutGene>,
    asr: AltSpliceResult,
    genome: &Fasta,
    code: &GeneticCode,
) -> (Vec<OutGene>, Vec<EventRecord>) {
    // one CdsModel per input transcript with a CDS; `owner` maps its id to (gene, tx)
    let mut cds_models = Vec::new();
    let mut owner: HashMap<String, (usize, usize)> = HashMap::new();
    for (gi, g) in models.iter().enumerate() {
        if g.strand == Strand::Unknown {
            continue;
        }
        for (ti, t) in g.transcripts.iter().enumerate() {
            let id = format!("m{gi}.{ti}");
            if let Some(m) = cds_model(&id, &t.contig, t.strand, &t.cds, t.cds_start_phase) {
                owner.insert(id, (gi, ti));
                cds_models.push(m);
            }
        }
    }

    let recon = reconcile(
        &asr.isoforms,
        &asr.loci,
        asr.events,
        &cds_models,
        genome,
        code,
    );

    // per gene: isoform index -> the (transcript, annotation) pairs it hosts
    type Hits<'a> = BTreeMap<usize, Vec<(usize, &'a CodingAnnotation)>>;
    let mut hits: Vec<Hits> = (0..models.len()).map(|_| BTreeMap::new()).collect();
    for (idx, anns) in recon.isoform_codings.iter().enumerate() {
        for ann in anns {
            let (gi, ti) = owner[&ann.model_id];
            hits[gi].entry(idx).or_default().push((ti, ann));
        }
    }

    for (g, gene_hits) in models.iter_mut().zip(&hits) {
        let Some(primary) = (0..g.transcripts.len())
            .filter(|&ti| !g.transcripts[ti].cds.is_empty())
            .rev()
            .max_by_key(|&ti| total_len(&g.transcripts[ti].cds))
        else {
            continue; // no CDS: emitted verbatim
        };
        let cds_introns: Vec<Vec<Coordset>> =
            g.transcripts.iter().map(|t| seg_introns(&t.cds)).collect();
        let mut taken: HashSet<String> = g
            .transcripts
            .iter()
            .map(|t| t.transcript_id.clone())
            .collect();
        let mut n = 1;
        let mut appended = Vec::new();
        for (&idx, anns) in gene_hits {
            let Some(&(_, primary_ann)) = anns.iter().find(|(ti, _)| *ti == primary) else {
                continue; // cannot host the primary CDS start
            };
            let iso = &asr.isoforms[idx];
            let iso_introns = iso.introns();
            let genuine = anns.iter().all(|&(ti, ann)| {
                is_genuine_alternative(
                    &iso.exons,
                    &iso_introns,
                    &g.transcripts[ti].cds,
                    &cds_introns[ti],
                    ann,
                )
            });
            if !genuine {
                continue;
            }
            while taken.contains(&format!("{}.iso{n}", g.gene_id)) {
                n += 1;
            }
            let t = iso_mrna(iso, primary_ann, &g.gene_id, n);
            taken.insert(t.transcript_id.clone());
            appended.push(t);
        }
        let (lend, rend) = span_of(&appended);
        if !appended.is_empty() {
            g.lend = g.lend.min(lend);
            g.rend = g.rend.max(rend);
        }
        g.transcripts.extend(appended);
    }

    (models, recon.events)
}

/// A consensus gene's CDS as a `CdsModel` for grafting, or `None` if it has no CDS.
fn to_cds_model(g: &CalledGene, id: String) -> Option<CdsModel> {
    cds_model(&id, &g.contig, g.orient, &g.cds, g.cds_start_phase)
}

/// A lend-sorted CDS as a `CdsModel`, or `None` if empty. The start codon is the 5'-most
/// CDS base shifted INTO frame by `start_phase` (the 5' phase of a 5'-partial model), so
/// a divergent isoform is projected in the model's reading frame.
fn cds_model(
    id: &str,
    contig: &str,
    strand: Strand,
    cds: &[Coordset],
    start_phase: u8,
) -> Option<CdsModel> {
    let first = cds.first()?;
    let last = cds.last()?;
    let phase = i64::from(start_phase);
    let (start_genomic, stop_genomic) = if strand == Strand::Minus {
        (last.rend - phase, first.lend)
    } else {
        (first.lend + phase, last.rend)
    };
    Some(CdsModel {
        id: id.to_string(),
        contig: contig.to_string(),
        strand,
        cds_segments: cds.to_vec(),
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
    // A CDS that begins mid-codon (`cds_start_phase > 0`) has no UTR room at its 5' end
    // — this keeps such genes unsalvaged, as when their partial codon was a 5'UTR stub.
    let use_salvage = g.five_utr.is_empty()
        && g.three_utr.is_empty()
        && g.cds_start_phase == 0
        && !utrs.is_empty();
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
        cds_start_phase: g.cds_start_phase,
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

/// Whether an isoform hosting model CDS `cds` (annotation `ann`) is a genuine
/// alternative to it — one that earns its own mRNA: it diverges within the CDS
/// (`coding_altered`, so not CDS-identical) and is not a mere terminal truncation
/// ([`is_terminal_truncation`]), i.e. it introduces a novel splice junction. The single
/// predicate shared by [`annotate`] and [`augment`].
fn is_genuine_alternative(
    iso_exons: &[Coordset],
    iso_introns: &[Coordset],
    cds: &[Coordset],
    cds_introns: &[Coordset],
    ann: &CodingAnnotation,
) -> bool {
    ann.coding_altered && !is_terminal_truncation(iso_exons, iso_introns, cds, cds_introns, ann)
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
        cds_start_phase: 0,
    }
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
            cds_start_phase: 0,
            score: 100.0,
            support: SupportFlags {
                raw_noncoding: 0.0,
                noncoding_equivalent: 0.0,
                score_ratio: 5.0,
                coding_length: 60,
                low_support: false,
            },
            promoted: false,
            features: vec![],
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
    fn mid_codon_cds_start_takes_no_salvaged_utr() {
        let mut g = consensus_gene(&[(10, 40), (60, 90)]);
        let utrs = Utrs {
            five: vec![cs(1, 5)],
            three: vec![],
            sources: vec![],
            contains: vec![],
        };
        // phase 0: the salvaged 5'UTR is adopted (unchanged behavior)
        assert_eq!(consensus_mrna(&g, "g1", &utrs).five_utr, vec![cs(1, 5)]);
        // a 5'-partial CDS starting mid-codon keeps its own (empty) UTRs and phase
        g.partial5 = true;
        g.cds_start_phase = 2;
        let m = consensus_mrna(&g, "g1", &utrs);
        assert!(m.five_utr.is_empty());
        assert_eq!(m.exons, vec![cs(10, 40), cs(60, 90)]);
        assert_eq!(m.cds_start_phase, 2);
    }

    #[test]
    fn to_cds_model_of_a_five_prime_partial_gene_starts_in_frame() {
        // CDS from the exon's first base, phase 2: the first complete codon is 2 bases in
        let mut g = consensus_gene(&[(10, 40), (60, 90)]);
        g.partial5 = true;
        g.cds_start_phase = 2;
        let m = to_cds_model(&g, "cons0".into()).unwrap();
        assert_eq!((m.start_genomic, m.stop_genomic), (12, 90));
        assert_eq!(
            m.cds_segments,
            vec![cs(10, 40), cs(60, 90)],
            "segments verbatim"
        );
        // minus strand: the 5' end is the last segment's rend
        g.orient = Strand::Minus;
        g.cds_start_phase = 1;
        let m = to_cds_model(&g, "cons0".into()).unwrap();
        assert_eq!((m.start_genomic, m.stop_genomic), (89, 10));
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
        assert_eq!(
            mrnas.len(),
            1,
            "CDS-identical UTR isoform must not duplicate the consensus"
        );
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
        assert!(
            m.three_utr.is_empty(),
            "no bogus 3'UTR from the truncated end"
        );
        assert_eq!(
            m.exons.first().unwrap().lend,
            10,
            "5'-terminal exon extended"
        );
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
        assert_eq!(
            mrnas.len(),
            1,
            "a 5'-truncated UTR donor is not its own mRNA"
        );
        let m = &mrnas[0];
        assert_eq!(m.cds, vec![cs(20, 40), cs(60, 90)], "CDS unchanged");
        assert_eq!(
            m.three_utr,
            vec![cs(91, 110)],
            "3'UTR recovered from a non-grafting transcript"
        );
        assert!(m.five_utr.is_empty());
        assert_eq!(
            m.exons.last().unwrap().rend,
            110,
            "3'-terminal exon extended"
        );
    }

    // ---- augment (assemble --models) ----

    fn mrna(tid: &str, exons: &[(i64, i64)], cds: &[(i64, i64)], phase: u8) -> OutTranscript {
        OutTranscript {
            transcript_id: tid.into(),
            contig: "chr1".into(),
            strand: Strand::Plus,
            exons: exons.iter().map(|&(l, r)| cs(l, r)).collect(),
            cds: cds.iter().map(|&(l, r)| cs(l, r)).collect(),
            five_utr: vec![],
            three_utr: vec![],
            attrs: vec![("Name".into(), vec![format!("{tid}-name")])],
            cds_start_phase: phase,
        }
    }

    fn model(id: &str, transcripts: Vec<OutTranscript>) -> OutGene {
        OutGene {
            gene_id: id.into(),
            contig: "chr1".into(),
            strand: Strand::Plus,
            lend: span_of(&transcripts).0,
            rend: span_of(&transcripts).1,
            attrs: vec![("Note".into(), vec!["kept".into()])],
            transcripts,
        }
    }

    /// One locus per isoform, on the isoform's strand.
    fn asr_of(isoforms: Vec<Isoform>) -> AltSpliceResult {
        let loci = isoforms
            .iter()
            .enumerate()
            .map(|(i, iso)| Locus {
                id: format!("L{i}"),
                contig: iso.contig.clone(),
                strand: iso.strand,
                isoform_indices: vec![i],
            })
            .collect();
        AltSpliceResult {
            isoforms,
            loci,
            events: vec![],
        }
    }

    fn atg_genome(len: usize, at: usize) -> Fasta {
        let mut g = vec![b'C'; len];
        g[at - 1..at + 2].copy_from_slice(b"ATG");
        write_fasta(&g)
    }

    fn tids(g: &OutGene) -> Vec<&str> {
        g.transcripts
            .iter()
            .map(|t| t.transcript_id.as_str())
            .collect()
    }

    fn attr<'a>(t: &'a OutTranscript, k: &str) -> Option<&'a str> {
        t.attrs
            .iter()
            .find(|(key, _)| key == k)
            .map(|(_, v)| v[0].as_str())
    }

    fn render(genes: &[OutGene]) -> String {
        let mut buf = Vec::new();
        crate::io::writer_gff3::write(&mut buf, genes).unwrap();
        String::from_utf8(buf).unwrap()
    }

    #[test]
    fn augment_keeps_base_verbatim_and_skips_cds_identical_utr_isoform() {
        let genome = atg_genome(110, 10);
        let input = || {
            vec![model(
                "g",
                vec![mrna(
                    "g.t1",
                    &[(10, 40), (60, 90)],
                    &[(10, 40), (60, 90)],
                    0,
                )],
            )]
        };
        // CDS-identical, UTR-extended isoform: not appended, and NO UTR salvage
        let asr = asr_of(vec![iso("t1", &[(1, 40), (60, 100)])]);
        let (out, _) = augment(input(), asr, &genome, &GeneticCode::default());
        assert_eq!(render(&out), render(&input()), "gene set emitted verbatim");
    }

    #[test]
    fn augment_appends_novel_junction_isoform_and_widens_span() {
        let genome = atg_genome(110, 10);
        let input = vec![model(
            "g",
            vec![mrna(
                "g.t1",
                &[(10, 40), (60, 90)],
                &[(10, 40), (60, 90)],
                0,
            )],
        )];
        let base = render(&input);
        let asr = asr_of(vec![iso("ri", &[(5, 95)])]); // retained intron, extended ends
        let (out, _) = augment(input, asr, &genome, &GeneticCode::default());
        assert_eq!(tids(&out[0]), vec!["g.t1", "g.iso1"]);
        let t = &out[0].transcripts[1];
        assert_eq!(attr(t, "support"), Some("transcript_isoform"));
        assert_eq!(attr(t, "coding_altered"), Some("true"));
        assert_eq!(
            t.cds.first().unwrap().lend,
            10,
            "CDS grafted from the model start"
        );
        assert_eq!((out[0].lend, out[0].rend), (5, 95), "span widened");
        // the input block is a verbatim prefix apart from the widened gene row
        let rendered = render(&out);
        let body = |s: &str| {
            s.lines()
                .skip(2)
                .take(4)
                .map(String::from)
                .collect::<Vec<_>>()
        };
        assert_eq!(body(&rendered), body(&base));
    }

    #[test]
    fn augment_skips_terminal_truncation() {
        let genome = atg_genome(150, 10);
        let ex = [(10, 40), (60, 90), (110, 140)];
        let input = vec![model("g", vec![mrna("g.t1", &ex, &ex, 0)])];
        let asr = asr_of(vec![iso("frag", &[(10, 40), (60, 90)])]);
        let (out, _) = augment(input, asr, &genome, &GeneticCode::default());
        assert_eq!(tids(&out[0]), vec!["g.t1"]);
    }

    #[test]
    fn augment_skips_opposite_strand_and_unattached_isoforms() {
        let genome = atg_genome(400, 10);
        let input = vec![model(
            "g",
            vec![mrna(
                "g.t1",
                &[(10, 40), (60, 90)],
                &[(10, 40), (60, 90)],
                0,
            )],
        )];
        let mut minus = iso("minus", &[(10, 90)]);
        minus.strand = Strand::Minus;
        let far = iso("far", &[(200, 250), (300, 350)]);
        let (out, _) = augment(
            input,
            asr_of(vec![minus, far]),
            &genome,
            &GeneticCode::default(),
        );
        assert_eq!(out.len(), 1, "an isoform overlapping no model is dropped");
        assert_eq!(tids(&out[0]), vec!["g.t1"]);
    }

    #[test]
    fn augment_dedups_against_every_mrna_of_the_gene() {
        let genome = atg_genome(150, 10);
        let full = [(10, 40), (60, 90), (110, 140)];
        let skip = [(10, 40), (110, 140)];
        let input = vec![model(
            "g",
            vec![mrna("g.t1", &full, &full, 0), mrna("g.t2", &skip, &skip, 0)],
        )];
        // `same2` is the second mRNA's structure (novel vs the primary, identical to t2):
        // not appended. `ri` retains the first intron: appended once.
        let asr = asr_of(vec![
            iso("same2", &skip),
            iso("ri", &[(10, 90), (110, 140)]),
        ]);
        let (out, _) = augment(input, asr, &genome, &GeneticCode::default());
        assert_eq!(tids(&out[0]), vec!["g.t1", "g.t2", "g.iso1"]);
        assert_eq!(out[0].transcripts[2].exons, vec![cs(10, 90), cs(110, 140)]);
    }

    #[test]
    fn augment_iso_id_avoids_existing_ids() {
        let genome = atg_genome(110, 10);
        let input = vec![model(
            "g",
            vec![
                mrna("g.t1", &[(10, 40), (60, 90)], &[(10, 40), (60, 90)], 0),
                mrna("g.iso1", &[(10, 40), (60, 90)], &[], 0), // non-coding, id taken
            ],
        )];
        let asr = asr_of(vec![iso("ri", &[(10, 90)])]);
        let (out, _) = augment(input, asr, &genome, &GeneticCode::default());
        assert_eq!(tids(&out[0]), vec!["g.t1", "g.iso1", "g.iso2"]);
    }

    #[test]
    fn augment_projects_a_five_prime_partial_model_in_frame() {
        // 5'-partial model: CDS rows start at 10 with phase 1, so the first complete
        // codon is at 11. A stop TAA at 41..43 is in frame from 11 (30 bases on) but not
        // from 10, so the retained-intron isoform's projected CDS must end exactly there.
        let mut g = vec![b'C'; 110];
        g[40..43].copy_from_slice(b"TAA");
        let genome = write_fasta(&g);
        let input = vec![model(
            "g",
            vec![mrna(
                "g.t1",
                &[(10, 40), (60, 90)],
                &[(10, 40), (60, 90)],
                1,
            )],
        )];
        let asr = asr_of(vec![iso("ri", &[(10, 90)])]);
        let (out, _) = augment(input, asr, &genome, &GeneticCode::default());
        assert_eq!(tids(&out[0]), vec!["g.t1", "g.iso1"]);
        let t = &out[0].transcripts[1];
        assert_eq!(t.cds, vec![cs(11, 43)], "projected from the in-frame start");
        assert_eq!(attr(t, "partial3"), Some("false"), "hit the in-frame stop");
        assert_eq!(out[0].transcripts[0].cds_start_phase, 1, "input phase kept");
    }

    #[test]
    fn cds_model_shifts_the_start_into_frame_on_both_strands() {
        let cds = [cs(10, 40), cs(60, 90)];
        let p = cds_model("m", "chr1", Strand::Plus, &cds, 2).unwrap();
        assert_eq!((p.start_genomic, p.stop_genomic), (12, 90));
        let m = cds_model("m", "chr1", Strand::Minus, &cds, 1).unwrap();
        assert_eq!((m.start_genomic, m.stop_genomic), (89, 10));
        let z = cds_model("m", "chr1", Strand::Minus, &cds, 0).unwrap();
        assert_eq!((z.start_genomic, z.stop_genomic), (90, 10));
    }
}
