//! Evidence ingestion for the consensus path.
//!
//! Parses the three evidence kinds straight from GFF3 into [`EvidenceChain`]s,
//! capturing the GFF column-2 source as the `ev_type` that keys the weights file.
//! This is intentionally self-contained (it does not route through the PASA
//! `Alignment`/`load_sources` path) because the consensus model is different:
//! gene predictions are read from **CDS** rows, while protein/transcript evidence is
//! read from spliced **match** chains (`Target=`), and every chain carries its source.
//!
//! Non-canonical (see `avoid-canonical-splice-bias`): no genome scanning happens here.
//! Candidate splice/start/stop sites are derived from these chains' structure in later
//! milestones; this layer only groups rows into chains and attaches `(class, weight)`.

use crate::consensus::weights::{EvClass, Weights};
use crate::error::{CombinrError, Result};
use crate::io::first_token;
use crate::io::gff::{self, AttrSep};
use crate::model::{Coordset, Strand};
use std::collections::HashMap;
use std::path::PathBuf;

/// One ingested evidence chain (a gene prediction's CDS, or a protein/transcript
/// match), with its source weight resolved.
#[derive(Debug, Clone, PartialEq)]
pub struct EvidenceChain {
    pub accession: String,
    pub ev_type: String,
    pub ev_class: EvClass,
    pub weight: f64,
    pub contig: String,
    /// Input strand, verbatim.
    pub orient: Strand,
    /// Chain span (forward genomic): `links` first lend .. last rend.
    pub span: Coordset,
    /// Exon/CDS segments, lend-sorted (forward genomic).
    pub links: Vec<Coordset>,
}

/// How a GFF3 file's rows become chains.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EvidenceKind {
    /// Gene predictions: group `CDS` rows by `Parent` (mRNA id).
    Prediction,
    /// Protein/transcript alignments: group rows carrying `Target=` by `ID`.
    Alignment,
}

/// Load all three evidence kinds, attaching `(class, weight)` from `weights` and
/// dropping (with a count) any chain whose source has no weight (EVM behaviour).
pub fn load_evidence(
    predictions: &[PathBuf],
    proteins: &[PathBuf],
    transcripts: &[PathBuf],
    weights: &Weights,
) -> Result<Vec<EvidenceChain>> {
    let mut chains = Vec::new();
    let mut skipped = 0usize;
    load_kind(
        predictions,
        EvidenceKind::Prediction,
        weights,
        &mut chains,
        &mut skipped,
    )?;
    load_kind(
        proteins,
        EvidenceKind::Alignment,
        weights,
        &mut chains,
        &mut skipped,
    )?;
    load_kind(
        transcripts,
        EvidenceKind::Alignment,
        weights,
        &mut chains,
        &mut skipped,
    )?;
    if skipped > 0 {
        eprintln!(
            "combinr consensus: skipped {skipped} evidence chain(s) whose source has no weight"
        );
    }
    Ok(chains)
}

fn load_kind(
    paths: &[PathBuf],
    kind: EvidenceKind,
    weights: &Weights,
    out: &mut Vec<EvidenceChain>,
    skipped: &mut usize,
) -> Result<()> {
    for p in paths {
        let text = std::fs::read_to_string(p).map_err(|e| CombinrError::Parse {
            file: p.display().to_string(),
            line: 0,
            msg: format!("cannot read evidence file: {e}"),
        })?;
        let (mut chains, sk) = ingest_text(&text, &p.display().to_string(), kind, weights)?;
        out.append(&mut chains);
        *skipped += sk;
    }
    Ok(())
}

/// Parse one GFF3 text into weighted chains, returning `(chains, num_skipped)`.
fn ingest_text(
    text: &str,
    file: &str,
    kind: EvidenceKind,
    weights: &Weights,
) -> Result<(Vec<EvidenceChain>, usize)> {
    let mut out = Vec::new();
    let mut skipped = 0usize;
    for rc in parse_chains(text, file, kind)? {
        match weights.lookup(&rc.ev_type) {
            Some((ev_class, weight)) => {
                let lend = rc.segs.first().unwrap().lend;
                let rend = rc.segs.last().unwrap().rend;
                out.push(EvidenceChain {
                    accession: rc.accession,
                    ev_type: rc.ev_type,
                    ev_class,
                    weight,
                    contig: rc.contig,
                    orient: rc.strand,
                    span: Coordset::new(lend, rend),
                    links: rc.segs,
                });
            }
            None => skipped += 1,
        }
    }
    Ok((out, skipped))
}

/// A grouped, unweighted chain straight from GFF3 rows.
struct RawChain {
    accession: String,
    ev_type: String,
    contig: String,
    strand: Strand,
    segs: Vec<Coordset>,
}

fn parse_chains(text: &str, file: &str, kind: EvidenceKind) -> Result<Vec<RawChain>> {
    let mut order: Vec<String> = Vec::new();
    let mut groups: HashMap<String, RawChain> = HashMap::new();

    for (lineno, raw) in text.lines().enumerate() {
        let Some(rec) = gff::record(raw, file, lineno + 1)? else {
            continue;
        };
        let attrs = gff::parse_attrs(rec.attrs, AttrSep::Eq);

        let key = match kind {
            EvidenceKind::Prediction => {
                if rec.ftype != "CDS" {
                    continue;
                }
                match attrs.get("Parent").or_else(|| attrs.get("ID")) {
                    Some(p) => p.clone(),
                    None => continue,
                }
            }
            EvidenceKind::Alignment => {
                let Some(target) = attrs.get("Target") else {
                    continue;
                };
                attrs
                    .get("ID")
                    .cloned()
                    .unwrap_or_else(|| first_token(target).to_string())
            }
        };

        let entry = groups.entry(key.clone()).or_insert_with(|| {
            order.push(key.clone());
            RawChain {
                accession: key.clone(),
                ev_type: rec.source.to_string(), // GFF column 2 = source
                contig: rec.contig.to_string(),
                strand: rec.strand,
                segs: Vec::new(),
            }
        });
        entry.segs.push(Coordset::new(rec.lend, rec.rend));
    }

    let mut out = Vec::with_capacity(order.len());
    for key in order {
        let mut rc = groups.remove(&key).expect("group present");
        rc.segs.sort_by_key(|c| c.lend);
        out.push(rc);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consensus::weights::Weights;

    fn weights() -> Weights {
        Weights::parse_str(
            "ABINITIO_PREDICTION fgenesh 1\nPROTEIN nap 5\nTRANSCRIPT gap2 1\n",
            "w",
        )
        .unwrap()
    }

    #[test]
    fn parses_prediction_cds_chain() {
        let gff = "\
chr1\tfgenesh\tgene\t100\t400\t.\t+\t.\tID=g1
chr1\tfgenesh\tmRNA\t100\t400\t.\t+\t.\tID=m1;Parent=g1
chr1\tfgenesh\tCDS\t300\t400\t.\t+\t2\tID=c1;Parent=m1
chr1\tfgenesh\tCDS\t100\t200\t.\t+\t0\tID=c1;Parent=m1";
        let chains = parse_chains(gff, "p.gff3", EvidenceKind::Prediction).unwrap();
        assert_eq!(chains.len(), 1);
        let c = &chains[0];
        assert_eq!(c.ev_type, "fgenesh");
        assert_eq!(c.accession, "m1");
        // sorted by lend even though rows were out of order
        assert_eq!(
            c.segs,
            vec![Coordset::new(100, 200), Coordset::new(300, 400)]
        );
    }

    #[test]
    fn parses_alignment_target_chain() {
        let gff = "\
chr1\tnap\tnucleotide_to_protein_match\t100\t200\t90\t+\t.\tID=a1;Target=P1 1 33
chr1\tnap\tnucleotide_to_protein_match\t300\t400\t88\t+\t.\tID=a1;Target=P1 34 66";
        let chains = parse_chains(gff, "pr.gff3", EvidenceKind::Alignment).unwrap();
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].ev_type, "nap");
        assert_eq!(chains[0].segs.len(), 2);
    }

    #[test]
    fn prediction_ignores_exon_rows_uses_cds() {
        // exon rows present but only CDS forms the chain
        let gff = "\
chr1\tfgenesh\texon\t90\t210\t.\t+\t.\tParent=m1
chr1\tfgenesh\tCDS\t100\t200\t.\t+\t0\tParent=m1";
        let chains = parse_chains(gff, "p", EvidenceKind::Prediction).unwrap();
        assert_eq!(chains[0].segs, vec![Coordset::new(100, 200)]);
    }

    #[test]
    fn ingest_attaches_class_and_skips_unweighted() {
        let gff = "\
chr1\tfgenesh\tCDS\t100\t200\t.\t+\t0\tParent=m1
chr1\tunknownsrc\tCDS\t500\t600\t.\t+\t0\tParent=m2";
        let (chains, skipped) =
            ingest_text(gff, "p.gff3", EvidenceKind::Prediction, &weights()).unwrap();
        assert_eq!(chains.len(), 1);
        assert_eq!(skipped, 1);
        assert_eq!(chains[0].ev_class, EvClass::AbinitioPrediction);
        assert_eq!(chains[0].weight, 1.0);
        assert_eq!(chains[0].span, Coordset::new(100, 200));
    }
}
