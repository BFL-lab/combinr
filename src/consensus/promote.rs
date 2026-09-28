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
//!
//! Evidence report: a promoted gene's features are its CDS segments (exon rows, typed by
//! position, EVM frames) and the gaps between them (intron rows). The PASA assembly path
//! does not track per-feature attribution, so EVERY feature of a promoted gene lists all
//! of the isoform's contained accessions, each paired with its GFF column-2 source
//! (`"transcript"` when the accession is not among the loaded transcript chains).

use crate::altsplice::{AltSpliceResult, Isoform, Locus};
use crate::consensus::engine::{CalledGene, FeatureKind, FeatureSupport};
use crate::consensus::exon::{ExonType, end_frame_for};
use crate::consensus::filter::SupportFlags;
use crate::io::fasta::Fasta;
use crate::model::{Coordset, Strand};
use crate::orf::GeneticCode;
use crate::orf::coords::SplicedTranscript;
use crate::orf::find_longest_orf;
use std::collections::HashMap;

/// Find transcript-ORF genes for uncovered loci and merge them with the trellis genes,
/// dropping the both-partial stubs they supersede. `sources` maps a transcript accession
/// to its GFF column-2 source, for the evidence attribution.
pub fn promote_and_merge(
    mut genes: Vec<CalledGene>,
    asr: &AltSpliceResult,
    genome: &Fasta,
    code: &GeneticCode,
    min_coding_length: i64,
    sources: &HashMap<String, String>,
) -> Vec<CalledGene> {
    let promoted = recover_transcript_loci(&genes, asr, genome, code, min_coding_length, sources);
    if promoted.is_empty() {
        return genes;
    }
    // Drop both-partial trellis stubs superseded by a promoted ORF gene.
    genes.retain(|g| {
        if !is_both_partial(g) {
            return true;
        }
        let (l, r) = g.span();
        !promoted.iter().any(|p| {
            let (pl, pr) = p.span();
            p.contig == g.contig && pl <= r && pr >= l
        })
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
    sources: &HashMap<String, String>,
) -> Vec<CalledGene> {
    let mut out = Vec::new();
    for locus in &asr.loci {
        if locus.isoform_indices.is_empty() {
            continue;
        }
        let (lend, rend) = locus_span(&asr.isoforms, locus);
        // Only a *real* (not both-partial) consensus gene counts as coverage.
        let covered = consensus_genes.iter().any(|g| {
            let (gl, gr) = g.span();
            g.contig == locus.contig && !is_both_partial(g) && gl <= rend && gr >= lend
        });
        if covered {
            continue;
        }
        let best = locus
            .isoform_indices
            .iter()
            .filter_map(|&i| orf_gene(&asr.isoforms[i], genome, code, min_coding_length, sources))
            .max_by_key(|g| g.support.coding_length);
        if let Some(g) = best {
            out.push(g);
        }
    }
    out
}

/// Build a promoted gene from an isoform's longest ORF, or `None` if too short / no ORF.
///
/// `coding_length` (which gates `min_coding_length` and is the score) is the CDS length,
/// matching the trellis path, whose coding_length counts the leading partial codon
/// because the CDS starts at the first exon base. Genes with a start codon (phase 0)
/// are unaffected.
fn orf_gene(
    iso: &Isoform,
    genome: &Fasta,
    code: &GeneticCode,
    min_coding_length: i64,
    sources: &HashMap<String, String>,
) -> Option<CalledGene> {
    let st = SplicedTranscript::new(&iso.exons, iso.strand);
    let seq = st.sequence(genome, &iso.contig)?;
    let orf = find_longest_orf(&seq, code)?;
    // A start-less ORF begins at its frame offset (0-2) from the transcript's 5' end; as
    // on the trellis path its CDS starts at the first base, the leading partial codon
    // carried as the first CDS row's phase (GFF3), not as a 5'UTR.
    let (cds_t_start, cds_start_phase) = if orf.has_start {
        (orf.t_start, 0)
    } else {
        (0, orf.t_start as u8)
    };
    let (cds, five_utr, three_utr) = st.cds_and_utrs(cds_t_start, orf.t_end);
    if cds.is_empty() {
        return None;
    }
    let coding_length: i64 = cds.iter().map(|c| c.len()).sum();
    if coding_length < min_coding_length {
        return None;
    }
    let evidence: Vec<(String, String)> = iso
        .contained_accs
        .iter()
        .map(|acc| {
            let src = sources.get(acc).map_or("transcript", String::as_str);
            (acc.clone(), src.to_string())
        })
        .collect();
    let features = cds_features(
        &cds,
        iso.strand,
        cds_start_phase,
        orf.has_start,
        orf.has_stop,
        &evidence,
    );
    Some(CalledGene {
        contig: iso.contig.clone(),
        orient: iso.strand,
        exons: iso.exons.clone(),
        cds,
        five_utr,
        three_utr,
        partial5: !orf.has_start,
        partial3: !orf.has_stop,
        cds_start_phase,
        score: coding_length as f64,
        support: SupportFlags {
            // no consensus noncoding baseline: 0.0 placeholder (never printed — the
            // report writes NA for promoted genes; 0.0 keeps CalledGene: PartialEq sane)
            raw_noncoding: 0.0,
            noncoding_equivalent: 0.0,
            score_ratio: f64::INFINITY, // no consensus noncoding to compare against
            coding_length,
            low_support: false,
        },
        promoted: true,
        features,
    })
}

/// Evidence-report features of a promoted gene: the CDS segments as exon rows (typed by
/// [`structural_type`] in transcription order; EVM frames — the codon position of the segment's 5'
/// base, `cum % 3 + 1` over the coding bases 5' of it, exactly as the candidate builder
/// assigns them — plus 3 on the minus strand) and the gaps between them as intron rows,
/// every row carrying `evidence`. The walk starts at `-start_phase` (as `cds_phases`
/// does) so a 5'-partial CDS's leading partial codon keeps its codon position. Returned
/// in ascending `lend` order.
fn cds_features(
    cds: &[Coordset],
    strand: Strand,
    start_phase: u8,
    has_start: bool,
    has_stop: bool,
    evidence: &[(String, String)],
) -> Vec<FeatureSupport> {
    let mut segs = cds.to_vec();
    segs.sort_by_key(|c| c.lend);
    let minus = strand == Strand::Minus;
    let n = segs.len();
    let mut features = Vec::with_capacity(2 * n);
    let mut cum = -i64::from(start_phase);
    // walk 5'->3' so the type and frame follow transcription order
    let order: Vec<usize> = if minus {
        (0..n).rev().collect()
    } else {
        (0..n).collect()
    };
    for (k, &i) in order.iter().enumerate() {
        let seg = segs[i];
        let start_frame = (cum.rem_euclid(3) + 1) as u8;
        let end_frame = end_frame_for(start_frame, seg.len());
        let shift = if minus { 3 } else { 0 };
        features.push(FeatureSupport {
            coords: seg,
            kind: FeatureKind::Exon {
                exon_type: structural_type(k, n, has_start, has_stop),
                start_frame: start_frame + shift,
                end_frame: end_frame + shift,
            },
            evidence: evidence.to_vec(),
        });
        cum += seg.len();
    }
    for w in segs.windows(2) {
        features.push(FeatureSupport {
            coords: w[0].gap_to(&w[1]),
            kind: FeatureKind::Intron,
            evidence: evidence.to_vec(),
        });
    }
    features.sort_by_key(|f| f.coords.lend);
    features
}

/// Exon type of the `k`-th of `n` CDS segments (transcription order) of a gene with/without
/// a start and stop codon: the first segment is Initial only when the gene has a start,
/// the last Terminal only when it has a stop (a lone segment Single when both), else
/// Internal. Matches the trellis path, where a 5'-partial path starts at an Internal (or
/// Terminal) candidate and a 3'-partial one ends at an Internal (or Initial) one.
fn structural_type(k: usize, n: usize, has_start: bool, has_stop: bool) -> ExonType {
    let first = k == 0 && has_start;
    let last = k + 1 == n && has_stop;
    match (n == 1, first, last) {
        (true, true, true) => ExonType::Single,
        (_, true, _) => ExonType::Initial,
        (_, _, true) => ExonType::Terminal,
        _ => ExonType::Internal,
    }
}

fn is_both_partial(g: &CalledGene) -> bool {
    g.partial5 && g.partial3
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
            cds_start_phase: 0,
            score: 100.0,
            support: SupportFlags {
                raw_noncoding: 0.0,
                noncoding_equivalent: 0.0,
                score_ratio: 2.0,
                coding_length: rend - lend + 1,
                low_support: false,
            },
            promoted: false,
            features: vec![],
        }
    }

    #[test]
    fn promotes_uncovered_transcript_locus() {
        let (asr, g) = asr_one();
        let out = promote_and_merge(
            vec![],
            &asr,
            &g,
            &GeneticCode::default(),
            15,
            &HashMap::new(),
        );
        assert_eq!(out.len(), 1);
        assert!(out[0].promoted);
        assert_eq!(out[0].cds, vec![Coordset::new(10, 30)]);
    }

    #[test]
    fn promoted_gene_features_are_cds_segments_with_transcript_evidence() {
        let (asr, g) = asr_one();
        let sources: HashMap<String, String> = [("t1".to_string(), "PASA".to_string())].into();
        let out = promote_and_merge(vec![], &asr, &g, &GeneticCode::default(), 15, &sources);
        let ev = vec![("t1".to_string(), "PASA".to_string())];
        assert_eq!(
            out[0].features,
            vec![FeatureSupport {
                coords: Coordset::new(10, 30),
                kind: FeatureKind::Exon {
                    exon_type: crate::consensus::exon::ExonType::Single,
                    start_frame: 1,
                    end_frame: 3,
                },
                evidence: ev,
            }]
        );
    }

    #[test]
    fn startless_orf_cds_starts_at_the_transcript_base_with_its_phase() {
        // transcript 10..30 (offsets 0..20): TAA at offset 0 kills frame 0, TGA at
        // offsets 5..7 kills frame 2; frame 1 (offsets 1..19, no ATG, no stop) is the
        // longest ORF: start-less, first complete codon 1 base in.
        let mut g = vec![b'C'; 60];
        g[9..12].copy_from_slice(b"TAA"); // 10..12
        g[14..17].copy_from_slice(b"TGA"); // 15..17
        let genome = write_fasta(&g);
        let (mut asr, _) = asr_one();
        asr.isoforms = vec![iso("t1", &[(10, 30)])];
        let out = promote_and_merge(
            vec![],
            &asr,
            &genome,
            &GeneticCode::default(),
            15,
            &HashMap::new(),
        );
        assert_eq!(out.len(), 1);
        let gene = &out[0];
        assert!(gene.partial5 && gene.partial3);
        assert_eq!(
            gene.cds,
            vec![Coordset::new(10, 28)],
            "ORF 11..28 + partial codon"
        );
        assert!(gene.five_utr.is_empty());
        assert_eq!(gene.cds_start_phase, 1);
        assert_eq!(
            gene.support.coding_length, 19,
            "the CDS length, partial codon included"
        );
        // the reported frame stays the codon position of the CDS's first base; with no
        // start or stop codon the single segment is typed Internal
        assert!(matches!(
            gene.features[0].kind,
            FeatureKind::Exon {
                exon_type: ExonType::Internal,
                start_frame: 3,
                ..
            }
        ));
        // the min_coding_length gate uses the CDS length (19)
        let promote = |min| {
            promote_and_merge(
                vec![],
                &asr,
                &genome,
                &GeneticCode::default(),
                min,
                &HashMap::new(),
            )
        };
        assert_eq!(promote(19).len(), 1);
        assert!(promote(20).is_empty());
    }

    #[test]
    fn cds_features_minus_strand_types_frames_and_introns() {
        // minus CDS 5'->3': 200..209 (10 bp, initial), 100..108 (terminal, cum 10 -> 2)
        let segs = [Coordset::new(100, 108), Coordset::new(200, 209)];
        let ev = vec![("t9".to_string(), "transcript".to_string())];
        let f = cds_features(&segs, Strand::Minus, 0, true, true, &ev);
        let kinds: Vec<(i64, i64, FeatureKind)> = f
            .iter()
            .map(|x| (x.coords.lend, x.coords.rend, x.kind))
            .collect();
        use crate::consensus::exon::ExonType::{Initial, Terminal};
        assert_eq!(
            kinds,
            vec![
                (
                    100,
                    108,
                    FeatureKind::Exon {
                        exon_type: Terminal,
                        start_frame: 5,
                        end_frame: 4
                    }
                ),
                (109, 199, FeatureKind::Intron),
                (
                    200,
                    209,
                    FeatureKind::Exon {
                        exon_type: Initial,
                        start_frame: 4,
                        end_frame: 4
                    }
                ),
            ]
        );
        assert!(f.iter().all(|x| x.evidence == ev));
    }

    #[test]
    fn cds_features_types_ends_by_start_stop_presence() {
        let segs = [Coordset::new(10, 20), Coordset::new(40, 50)];
        let types = |has_start, has_stop| -> Vec<ExonType> {
            cds_features(&segs, Strand::Plus, 0, has_start, has_stop, &[])
                .iter()
                .filter_map(|f| match f.kind {
                    FeatureKind::Exon { exon_type, .. } => Some(exon_type),
                    FeatureKind::Intron => None,
                })
                .collect()
        };
        assert_eq!(
            types(true, false),
            vec![ExonType::Initial, ExonType::Internal]
        );
        assert_eq!(
            types(false, true),
            vec![ExonType::Internal, ExonType::Terminal]
        );
    }

    #[test]
    fn real_consensus_gene_blocks_promotion() {
        let (asr, g) = asr_one();
        let real = gene(5, 35, false, false); // complete consensus gene over the locus
        let out = promote_and_merge(
            vec![real],
            &asr,
            &g,
            &GeneticCode::default(),
            15,
            &HashMap::new(),
        );
        assert_eq!(out.len(), 1);
        assert!(!out[0].promoted, "real consensus gene blocks promotion");
    }

    #[test]
    fn both_partial_stub_is_superseded_by_promoted_orf() {
        let (asr, g) = asr_one();
        let stub = gene(12, 28, true, true); // internal-exon-only stub over the locus
        let out = promote_and_merge(
            vec![stub],
            &asr,
            &g,
            &GeneticCode::default(),
            15,
            &HashMap::new(),
        );
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
        let out = promote_and_merge(
            vec![stub],
            &asr,
            &g,
            &GeneticCode::default(),
            300,
            &HashMap::new(),
        );
        assert_eq!(out.len(), 1);
        assert!(!out[0].promoted);
    }
}
