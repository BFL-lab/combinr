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

/// The deterministic output order and gene IDs: sorted by `(contig, lend, rend, strand)`
/// (stable), with ID `consensus.{contig}.g{rank}` (1-based, numbered across all contigs).
/// The single source for the GFF3/GTF renderers (default and `--alt-splice`) and the
/// evidence report, so their gene IDs agree by construction.
pub fn ordered_with_ids(genes: &[CalledGene]) -> Vec<(String, &CalledGene)> {
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
        .map(|(i, g)| (format!("consensus.{}.g{}", g.contig, i + 1), g))
        .collect()
}

/// Convert called genes to output genes, sorted deterministically by position.
pub fn to_out_genes(genes: &[CalledGene]) -> Vec<OutGene> {
    ordered_with_ids(genes)
        .into_iter()
        .map(|(gene_id, g)| {
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
