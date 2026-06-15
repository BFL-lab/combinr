//! Recover transcript-only loci via de-novo ORF finding.
//!
//! `consensus` builds coding genes from weighted evidence, but a locus with only
//! transcript evidence yields no real consensus CDS — the trellis can only chain the
//! transcripts' *internal* exons into a both-ends-partial stub (no start/stop). When
//! `--promote-transcript-orfs` is set, this step assembles the transcripts (the PASA
//! path) and, at each locus with no genuine consensus gene, finds the longest ORF
//! (genetic-code aware) and emits the best isoform as a gene tagged
//! `support=transcript_orf`, **superseding** any both-partial stub there. Loci already
//! covered by a real (non-both-partial) consensus gene are left for the full alt-isoform
//! integration (Phase 2).

use crate::altsplice::{AltSpliceResult, Isoform, Locus};
use crate::consensus::engine::CalledGene;
use crate::consensus::filter::SupportFlags;
use crate::io::fasta::Fasta;
use crate::orf::GeneticCode;
use crate::orf::coords::SplicedTranscript;
use crate::orf::find_longest_orf;

/// Find transcript-ORF genes for uncovered loci and merge them with the trellis genes,
/// dropping the both-partial stubs they supersede.
pub fn promote_and_merge(
    mut genes: Vec<CalledGene>,
    asr: &AltSpliceResult,
    genome: &Fasta,
    code: &GeneticCode,
    min_coding_length: i64,
) -> Vec<CalledGene> {
    let promoted = recover_transcript_loci(&genes, asr, genome, code, min_coding_length);
    if promoted.is_empty() {
        return genes;
    }
    // Drop both-partial trellis stubs superseded by a promoted ORF gene.
    genes.retain(|g| {
        if !is_both_partial(g) {
            return true;
        }
        let (l, r) = (gene_lend(g), gene_rend(g));
        !promoted
            .iter()
            .any(|p| p.contig == g.contig && gene_lend(p) <= r && gene_rend(p) >= l)
    });
    genes.extend(promoted);
    genes
}

/// Build promoted ORF genes for loci not covered by a real consensus gene.
fn recover_transcript_loci(
    consensus_genes: &[CalledGene],
    asr: &AltSpliceResult,
    genome: &Fasta,
    code: &GeneticCode,
    min_coding_length: i64,
) -> Vec<CalledGene> {
    let mut out = Vec::new();
    for locus in &asr.loci {
        if locus.isoform_indices.is_empty() {
            continue;
        }
        let (lend, rend) = locus_span(&asr.isoforms, locus);
        // Only a *real* (not both-partial) consensus gene counts as coverage.
        let covered = consensus_genes.iter().any(|g| {
            g.contig == locus.contig
                && !is_both_partial(g)
                && gene_lend(g) <= rend
                && gene_rend(g) >= lend
        });
        if covered {
            continue;
        }
        let best = locus
            .isoform_indices
            .iter()
            .filter_map(|&i| orf_gene(&asr.isoforms[i], genome, code, min_coding_length))
            .max_by_key(|g| g.support.coding_length);
        if let Some(g) = best {
            out.push(g);
        }
    }
    out
}

/// Build a promoted gene from an isoform's longest ORF, or `None` if too short / no ORF.
fn orf_gene(
    iso: &Isoform,
    genome: &Fasta,
    code: &GeneticCode,
    min_coding_length: i64,
) -> Option<CalledGene> {
    let st = SplicedTranscript::new(&iso.exons, iso.strand);
    let seq = st.sequence(genome, &iso.contig)?;
    let orf = find_longest_orf(&seq, code)?;
    let coding_length = orf.len() as i64;
    if coding_length < min_coding_length {
        return None;
    }
    let cds = st.genomic_segments_for_tspan(orf.t_start, orf.t_end);
    if cds.is_empty() {
        return None;
    }
    Some(CalledGene {
        contig: iso.contig.clone(),
        orient: iso.strand,
        exons: iso.exons.clone(),
        cds,
        five_utr: st.genomic_segments_for_tspan(0, orf.t_start),
        three_utr: st.genomic_segments_for_tspan(orf.t_end, st.len()),
        partial5: !orf.has_start,
        partial3: !orf.has_stop,
        score: coding_length as f64,
        support: SupportFlags {
            score_ratio: f64::INFINITY, // no consensus noncoding to compare against
            coding_length,
            low_support: false,
        },
        promoted: true,
    })
}

fn is_both_partial(g: &CalledGene) -> bool {
    g.partial5 && g.partial3
}

fn gene_lend(g: &CalledGene) -> i64 {
    g.exons.first().map(|c| c.lend).unwrap_or(i64::MAX)
}

fn gene_rend(g: &CalledGene) -> i64 {
    g.exons.last().map(|c| c.rend).unwrap_or(i64::MIN)
}

fn locus_span(isoforms: &[Isoform], locus: &Locus) -> (i64, i64) {
    let lend = locus
        .isoform_indices
        .iter()
        .map(|&i| isoforms[i].exons.first().unwrap().lend)
        .min()
        .unwrap();
    let rend = locus
        .isoform_indices
        .iter()
        .map(|&i| isoforms[i].exons.last().unwrap().rend)
        .max()
        .unwrap();
    (lend, rend)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::altsplice::AltSpliceResult;
    use crate::model::{Coordset, Strand};
    use std::collections::BTreeSet;
    use std::io::Write;

    fn write_fasta(seq: &[u8]) -> Fasta {
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        writeln!(tmp, ">chr1").unwrap();
        tmp.write_all(seq).unwrap();
        writeln!(tmp).unwrap();
        Fasta::load(tmp.path()).unwrap()
    }

    fn iso(id: &str, segs: &[(i64, i64)]) -> Isoform {
        Isoform {
            id: id.into(),
            contig: "chr1".into(),
            strand: Strand::Plus,
            exons: segs.iter().map(|&(l, r)| Coordset::new(l, r)).collect(),
            contained_accs: vec![id.into()],
            source_set: BTreeSet::new(),
            is_fl: false,
        }
    }

    fn asr_one() -> (AltSpliceResult, Fasta) {
        // single-exon transcript 10..30 with a clean ATG..TAA ORF (21 nt)
        let mut g = vec![b'C'; 60];
        let put = |g: &mut Vec<u8>, p: i64, s: &[u8]| {
            for (k, &b) in s.iter().enumerate() {
                g[(p - 1) as usize + k] = b;
            }
        };
        put(&mut g, 10, b"ATG");
        put(&mut g, 28, b"TAA");
        let genome = write_fasta(&g);
        let isoforms = vec![iso("t1", &[(10, 30)])];
        let loci = vec![Locus {
            id: "L1".into(),
            contig: "chr1".into(),
            strand: Strand::Plus,
            isoform_indices: vec![0],
        }];
        (
            AltSpliceResult {
                isoforms,
                loci,
                events: vec![],
            },
            genome,
        )
    }

    fn gene(lend: i64, rend: i64, partial5: bool, partial3: bool) -> CalledGene {
        CalledGene {
            contig: "chr1".into(),
            orient: Strand::Plus,
            exons: vec![Coordset::new(lend, rend)],
            cds: vec![Coordset::new(lend, rend)],
            five_utr: vec![],
            three_utr: vec![],
            partial5,
            partial3,
            score: 100.0,
            support: SupportFlags {
                score_ratio: 2.0,
                coding_length: rend - lend + 1,
                low_support: false,
            },
            promoted: false,
        }
    }

    #[test]
    fn promotes_uncovered_transcript_locus() {
        let (asr, g) = asr_one();
        let out = promote_and_merge(vec![], &asr, &g, &GeneticCode::default(), 15);
        assert_eq!(out.len(), 1);
        assert!(out[0].promoted);
        assert_eq!(out[0].cds, vec![Coordset::new(10, 30)]);
    }

    #[test]
    fn real_consensus_gene_blocks_promotion() {
        let (asr, g) = asr_one();
        let real = gene(5, 35, false, false); // complete consensus gene over the locus
        let out = promote_and_merge(vec![real], &asr, &g, &GeneticCode::default(), 15);
        assert_eq!(out.len(), 1);
        assert!(!out[0].promoted, "real consensus gene blocks promotion");
    }

    #[test]
    fn both_partial_stub_is_superseded_by_promoted_orf() {
        let (asr, g) = asr_one();
        let stub = gene(12, 28, true, true); // internal-exon-only stub over the locus
        let out = promote_and_merge(vec![stub], &asr, &g, &GeneticCode::default(), 15);
        // the stub is dropped, replaced by the promoted ORF gene
        assert_eq!(out.len(), 1);
        assert!(out[0].promoted);
        assert_eq!(out[0].cds, vec![Coordset::new(10, 30)]);
    }

    #[test]
    fn short_orf_keeps_the_stub() {
        let (asr, g) = asr_one();
        let stub = gene(12, 28, true, true);
        // require 300 nt: the 21 nt ORF doesn't qualify, so nothing is promoted and the
        // stub stays.
        let out = promote_and_merge(vec![stub], &asr, &g, &GeneticCode::default(), 300);
        assert_eq!(out.len(), 1);
        assert!(!out[0].promoted);
    }
}
