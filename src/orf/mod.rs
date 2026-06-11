//! Optional ORF/UTR reconciliation.
//!
//! Given an external gene-prediction GFF3 (already-mapped CDS, possibly several
//! coding models per locus) and the genome FASTA, graft each predicted CDS onto
//! the isoforms of its locus:
//!
//! - an isoform whose intron structure across the CDS region matches the model
//!   **inherits** the CDS verbatim;
//! - a **divergent** isoform is projected from the *same* predicted start codon
//!   and translated along its own exon structure to the first in-frame stop (a
//!   premature stop simply becomes that isoform's ORF stop).
//!
//! 5'/3' UTRs follow from the resulting CDS bounds, and each alt-splice event is
//! tagged by the region (5'UTR / CDS / 3'UTR) it falls in. See
//! `orf-from-external-gff3`.

pub mod coords;
pub mod translate;

use crate::altsplice::{EventRecord, Isoform, Locus, RegionClass};
use crate::error::{CombinrError, Result};
use crate::io::fasta::Fasta;
use crate::io::gff3::parse_attrs;
use crate::model::{Coordset, Strand};
use coords::SplicedTranscript;
use std::collections::HashMap;
use std::path::Path;

/// A predicted coding model: CDS already mapped to the genome.
#[derive(Debug, Clone)]
pub struct CdsModel {
    pub id: String,
    pub contig: String,
    pub strand: Strand,
    /// CDS segments, lend-sorted.
    pub cds_segments: Vec<Coordset>,
    /// 5' base of the start codon (genomic).
    pub start_genomic: i64,
    /// 3' base of the CDS (genomic).
    pub stop_genomic: i64,
}

impl CdsModel {
    fn span(&self) -> Coordset {
        Coordset {
            lend: self.cds_segments.first().unwrap().lend,
            rend: self.cds_segments.last().unwrap().rend,
        }
    }
    /// Introns of the CDS (gaps between segments) — robust to terminal
    /// stop-codon conventions.
    fn introns(&self) -> Vec<Coordset> {
        cds_introns(&self.cds_segments)
    }
}

fn cds_introns(segs: &[Coordset]) -> Vec<Coordset> {
    segs.windows(2)
        .map(|w| Coordset {
            lend: w[0].rend + 1,
            rend: w[1].lend - 1,
        })
        .collect()
}

/// One coding annotation grafted onto an isoform.
#[derive(Debug, Clone)]
pub struct CodingAnnotation {
    pub model_id: String,
    pub cds_segments: Vec<Coordset>,
    pub five_utr: Vec<Coordset>,
    pub three_utr: Vec<Coordset>,
    /// `true` when the isoform diverges from the model within the CDS.
    pub coding_altered: bool,
    pub partial3: bool,
}

/// Result of reconciliation: per-isoform coding annotations (parallel to the
/// isoform table) and the alt-splice events with region tags filled in.
pub struct ReconcileResult {
    pub isoform_codings: Vec<Vec<CodingAnnotation>>,
    pub events: Vec<EventRecord>,
}

/// Parse a gene-prediction GFF3 into coding models (CDS rows grouped by mRNA
/// `Parent`).
pub fn parse_cds_models(path: &Path) -> Result<Vec<CdsModel>> {
    let text = std::fs::read_to_string(path).map_err(|e| CombinrError::Parse {
        file: path.display().to_string(),
        line: 0,
        msg: format!("cannot read gene-prediction GFF3: {e}"),
    })?;

    struct Acc {
        contig: String,
        strand: Strand,
        segs: Vec<Coordset>,
    }
    let mut order: Vec<String> = Vec::new();
    let mut groups: HashMap<String, Acc> = HashMap::new();

    for (lineno, raw) in text.lines().enumerate() {
        let line = raw.trim_end();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 || cols[2] != "CDS" {
            continue;
        }
        let attrs = parse_attrs(cols[8]);
        let parent = match attrs.get("Parent").or_else(|| attrs.get("ID")) {
            Some(p) => p.clone(),
            None => continue,
        };
        let err = |msg: &str| CombinrError::Parse {
            file: path.display().to_string(),
            line: lineno + 1,
            msg: msg.to_string(),
        };
        let lend: i64 = cols[3].parse().map_err(|_| err("bad CDS start"))?;
        let rend: i64 = cols[4].parse().map_err(|_| err("bad CDS end"))?;
        let strand = Strand::from_char(cols[6].chars().next().unwrap_or('.'));
        let acc = groups.entry(parent.clone()).or_insert_with(|| {
            order.push(parent.clone());
            Acc {
                contig: cols[0].to_string(),
                strand,
                segs: Vec::new(),
            }
        });
        acc.segs.push(Coordset::new(lend, rend));
    }

    let mut models = Vec::with_capacity(order.len());
    for id in order {
        let mut acc = groups.remove(&id).unwrap();
        acc.segs.sort_by_key(|c| c.lend);
        let (first, last) = (acc.segs.first().unwrap(), acc.segs.last().unwrap());
        let (start_genomic, stop_genomic) = if acc.strand == Strand::Minus {
            (last.rend, first.lend)
        } else {
            (first.lend, last.rend)
        };
        models.push(CdsModel {
            id,
            contig: acc.contig,
            strand: acc.strand,
            cds_segments: acc.segs,
            start_genomic,
            stop_genomic,
        });
    }
    Ok(models)
}

/// Reconcile predicted CDS models onto isoforms and tag events by region.
pub fn reconcile(
    isoforms: &[Isoform],
    loci: &[Locus],
    mut events: Vec<EventRecord>,
    models: &[CdsModel],
    genome: &Fasta,
) -> ReconcileResult {
    let mut isoform_codings: Vec<Vec<CodingAnnotation>> = vec![Vec::new(); isoforms.len()];

    for locus in loci {
        let isos = &locus.isoform_indices;
        if isos.is_empty() {
            continue;
        }
        let locus_lend = isos
            .iter()
            .map(|&i| isoforms[i].exons.first().unwrap().lend)
            .min()
            .unwrap();
        let locus_rend = isos
            .iter()
            .map(|&i| isoforms[i].exons.last().unwrap().rend)
            .max()
            .unwrap();
        let locus_span = Coordset {
            lend: locus_lend,
            rend: locus_rend,
        };

        let locus_models: Vec<&CdsModel> = models
            .iter()
            .filter(|m| {
                m.contig == locus.contig
                    && m.strand == locus.strand
                    && m.span().overlaps_inclusive(&locus_span)
            })
            .collect();

        for &iso_idx in isos {
            for model in &locus_models {
                if let Some(ann) = graft(&isoforms[iso_idx], model, genome) {
                    isoform_codings[iso_idx].push(ann);
                }
            }
        }
    }

    // tag each event by region relative to isoform_a's (first) coding annotation
    let id_to_idx: HashMap<&str, usize> = isoforms
        .iter()
        .enumerate()
        .map(|(i, iso)| (iso.id.as_str(), i))
        .collect();
    for e in &mut events {
        if let Some(&idx) = id_to_idx.get(e.isoform_a.as_str())
            && let Some(ann) = isoform_codings[idx].first()
        {
            e.region = Some(region_of_event(e, ann, isoforms[idx].strand));
        }
    }

    ReconcileResult {
        isoform_codings,
        events,
    }
}

/// Graft one model onto one isoform, returning the coding annotation, or `None`
/// if the isoform cannot host the model (start codon not in any exon).
fn graft(iso: &Isoform, model: &CdsModel, genome: &Fasta) -> Option<CodingAnnotation> {
    let st = SplicedTranscript::new(&iso.exons, iso.strand);
    let t_start = st.genomic_to_tpos(model.start_genomic)?;

    // structural match across the CDS: same introns within the CDS span.
    let span = model.span();
    let iso_cds_introns: Vec<Coordset> = iso
        .introns()
        .into_iter()
        .filter(|i| i.lend >= span.lend && i.rend <= span.rend)
        .collect();
    let stop_tpos = st.genomic_to_tpos(model.stop_genomic);
    let clean = iso_cds_introns == model.introns() && stop_tpos.is_some();

    let (cds_t_start, cds_t_end, partial3, coding_altered) = if clean {
        let t_stop = stop_tpos.unwrap();
        (t_start, t_stop.max(t_start) + 1, false, false)
    } else {
        let seq = st.sequence(genome, &iso.contig)?;
        let proj = st.project_orf(&seq, t_start);
        (proj.cds_t_start, proj.cds_t_end, !proj.hit_stop, true)
    };

    Some(CodingAnnotation {
        model_id: model.id.clone(),
        cds_segments: st.genomic_segments_for_tspan(cds_t_start, cds_t_end),
        five_utr: st.genomic_segments_for_tspan(0, cds_t_start),
        three_utr: st.genomic_segments_for_tspan(cds_t_end, st.len()),
        coding_altered,
        partial3,
    })
}

/// Classify an event's location relative to a CDS (strand-aware).
fn region_of_event(e: &EventRecord, ann: &CodingAnnotation, strand: Strand) -> RegionClass {
    let cds_lend = ann.cds_segments.first().map(|c| c.lend).unwrap_or(i64::MAX);
    let cds_rend = ann.cds_segments.last().map(|c| c.rend).unwrap_or(i64::MIN);
    let ev_lend = e.coords.iter().map(|c| c.lend).min().unwrap_or(0);
    let ev_rend = e.coords.iter().map(|c| c.rend).max().unwrap_or(0);
    let plus = strand != Strand::Minus;

    if ev_rend >= cds_lend && ev_lend <= cds_rend {
        RegionClass::Cds
    } else if ev_rend < cds_lend {
        if plus {
            RegionClass::FivePrimeUtr
        } else {
            RegionClass::ThreePrimeUtr
        }
    } else if plus {
        RegionClass::ThreePrimeUtr
    } else {
        RegionClass::FivePrimeUtr
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::altsplice::locus::build_isoforms;
    use crate::assemble::assemble_cluster;
    use crate::model::{Alignment, Segment};
    use std::io::Write;

    // Build a tiny genome where a clean 2-exon ORF and a retained-intron variant
    // can be projected.
    fn genome() -> Fasta {
        // chr1: design so the spliced ORF starts with ATG and the retained-intron
        // isoform hits a premature stop in the intron.
        // We construct exons: e1=1..6, e2=16..30 (intron 7..15).
        // Spliced (+) ORF for clean isoform e1+e2:
        //   e1: ATG AAA (1..6)
        //   e2: ... designed so a stop appears only after several codons.
        // Retained-intron isoform e1..e2 contiguous (1..30): the intron seq
        //   contains an in-frame stop → premature.
        let mut seq = vec![b'A'; 40];
        let set = |s: &mut Vec<u8>, pos1: i64, bytes: &[u8]| {
            for (k, &b) in bytes.iter().enumerate() {
                s[(pos1 - 1) as usize + k] = b;
            }
        };
        // exon1 1..6: ATGAAA
        set(&mut seq, 1, b"ATGAAA");
        // intron 7..15: put TAA at 7..9 (in-frame for the retained variant: after
        // ATGAAA = 2 codons, next codon starts at 7) → premature stop.
        set(&mut seq, 7, b"TAACCCGGG"); // 7..15
        // exon2 16..30: GGG CCC ... TAA somewhere later for the clean ORF
        //   clean spliced seq: ATGAAA + (16..30). Put a stop at 16..18? that would
        //   be right after the 2 codons too. Put first stop at 25..27.
        set(&mut seq, 16, b"GGGCCCAAATAAGGG"); // 16..30 (15 nt); TAA at 25..27
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        writeln!(tmp, ">chr1").unwrap();
        tmp.write_all(&seq).unwrap();
        writeln!(tmp).unwrap();
        Fasta::load(tmp.path()).unwrap()
    }

    fn iso_from(acc: &str, segs: &[(i64, i64)]) -> Isoform {
        let mut a = Alignment::new(
            acc,
            segs.iter().map(|&(l, r)| Segment::new(l, r)).collect(),
            Strand::Plus,
        );
        a.contig = "chr1".to_string();
        a.spliced_orient = Strand::Plus;
        let asms = assemble_cluster(&[a], 20).unwrap();
        build_isoforms(&asms).pop().unwrap()
    }

    #[test]
    fn clean_isoform_inherits_and_divergent_projects_premature_stop() {
        let g = genome();
        // model CDS = clean spliced ORF over exons e1 (1..6) and e2 (16..30),
        // start at 1, stop at 27 (the TAA at 25..27).
        let model = CdsModel {
            id: "m1".into(),
            contig: "chr1".into(),
            strand: Strand::Plus,
            cds_segments: vec![Coordset::new(1, 6), Coordset::new(16, 27)],
            start_genomic: 1,
            stop_genomic: 27,
        };

        // clean isoform: same 2-exon structure → inherits, not altered.
        let clean = iso_from("clean", &[(1, 6), (16, 30)]);
        let ann_clean = graft(&clean, &model, &g).unwrap();
        assert!(!ann_clean.coding_altered, "matching structure inherits CDS");
        assert_eq!(
            ann_clean.cds_segments,
            vec![Coordset::new(1, 6), Coordset::new(16, 27)]
        );
        assert_eq!(ann_clean.three_utr, vec![Coordset::new(28, 30)]);

        // retained-intron isoform: single exon 1..30 → diverges, projects from the
        // same start (pos 1) and hits the premature TAA at 7..9.
        let retained = iso_from("retained", &[(1, 30)]);
        let ann_div = graft(&retained, &model, &g).unwrap();
        assert!(ann_div.coding_altered, "retained intron diverges");
        // CDS = 1..9 (ATG AAA TAA), premature stop.
        assert_eq!(ann_div.cds_segments, vec![Coordset::new(1, 9)]);
    }
}
