//! Render consensus genes to the format-neutral [`OutGene`] model (reused by the GFF3
//! and GTF writers). M3 emits gene/mRNA/exon with support attributes; the CDS/UTR
//! features land in M4.

use crate::consensus::engine::CalledGene;
use crate::io::out_model::{OutGene, OutTranscript};

/// Convert called genes to output genes, sorted deterministically by position.
pub fn to_out_genes(genes: &[CalledGene]) -> Vec<OutGene> {
    let mut order: Vec<&CalledGene> = genes.iter().collect();
    order.sort_by(|a, b| {
        let ak = (
            a.contig.as_str(),
            gene_lend(a),
            gene_rend(a),
            a.orient.to_char(),
        );
        let bk = (
            b.contig.as_str(),
            gene_lend(b),
            gene_rend(b),
            b.orient.to_char(),
        );
        ak.cmp(&bk)
    });

    order
        .into_iter()
        .enumerate()
        .map(|(i, g)| {
            let gene_id = format!("evm.{}.g{}", g.contig, i + 1);
            let lend = gene_lend(g);
            let rend = gene_rend(g);
            // Promoted transcript-ORF genes have no consensus noncoding baseline, so the
            // ratio is reported as NA rather than the placeholder infinity.
            let ratio = if g.promoted {
                "NA".to_string()
            } else if g.support.score_ratio.is_finite() {
                format!("{:.2}", g.support.score_ratio)
            } else {
                "inf".to_string()
            };
            let mut attrs = vec![
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
            // Tag only promoted genes, so default consensus output is unchanged.
            if g.promoted {
                attrs.insert(0, ("support".into(), vec!["transcript_orf".into()]));
            }
            OutGene {
                gene_id: gene_id.clone(),
                contig: g.contig.clone(),
                strand: g.orient,
                lend,
                rend,
                transcripts: vec![OutTranscript {
                    transcript_id: format!("{gene_id}.mRNA"),
                    contig: g.contig.clone(),
                    strand: g.orient,
                    exons: g.exons.clone(),
                    cds: g.cds.clone(),
                    five_utr: g.five_utr.clone(),
                    three_utr: g.three_utr.clone(),
                    attrs,
                }],
            }
        })
        .collect()
}

fn gene_lend(g: &CalledGene) -> i64 {
    g.exons.iter().map(|c| c.lend).min().unwrap_or(0)
}

fn gene_rend(g: &CalledGene) -> i64 {
    g.exons.iter().map(|c| c.rend).max().unwrap_or(0)
}
