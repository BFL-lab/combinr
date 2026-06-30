//! Render consensus genes to the format-neutral [`OutGene`] model (reused by the GFF3
//! and GTF writers). M3 emits gene/mRNA/exon with support attributes; the CDS/UTR
//! features land in M4.

use crate::consensus::engine::CalledGene;
use crate::io::out_model::{OutGene, OutTranscript};

/// The `score_ratio` attribute string: `NA` for a promoted transcript-ORF gene (no
/// consensus noncoding baseline to compare against), `{:.2}` for a finite ratio, else
/// `inf`. Shared by the default and `--alt-splice` consensus renderers.
pub(crate) fn format_ratio(g: &CalledGene) -> String {
    if g.promoted {
        "NA".to_string()
    } else if g.support.score_ratio.is_finite() {
        format!("{:.2}", g.support.score_ratio)
    } else {
        "inf".to_string()
    }
}

/// The six consensus attributes shared by both renderers, in their load-bearing GFF order:
/// `score, score_ratio, coding_length, low_support, partial5, partial3`. Callers prepend
/// the `support` tag (and the `--alt-splice` path appends `sources`/`contains`) themselves,
/// so the two modes stay intentionally — not accidentally — in sync.
pub(crate) fn consensus_core_attrs(g: &CalledGene) -> Vec<(String, Vec<String>)> {
    vec![
        ("score".into(), vec![format!("{:.1}", g.score)]),
        ("score_ratio".into(), vec![format_ratio(g)]),
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
    ]
}

/// Convert called genes to output genes, sorted deterministically by position.
pub fn to_out_genes(genes: &[CalledGene]) -> Vec<OutGene> {
    let mut order: Vec<&CalledGene> = genes.iter().collect();
    order.sort_by(|a, b| {
        let (al, ar) = a.span();
        let (bl, br) = b.span();
        (a.contig.as_str(), al, ar, a.orient.to_char()).cmp(&(
            b.contig.as_str(),
            bl,
            br,
            b.orient.to_char(),
        ))
    });

    order
        .into_iter()
        .enumerate()
        .map(|(i, g)| {
            let gene_id = format!("consensus.{}.g{}", g.contig, i + 1);
            let (lend, rend) = g.span();
            let mut attrs = consensus_core_attrs(g);
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
