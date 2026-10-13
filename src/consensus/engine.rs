//! Per-region consensus: run the trellis on both strands and recurse into large tails.
//!
//! The minus strand is handled by **reverse-complementing the region and running the
//! exact same forward machinery**, then transposing the resulting genes back to forward
//! coordinates — so the trellis and compatibility logic never need per-strand branches.
//! Forward and minus are explored independently across the whole region (leaning
//! anti-false-negative; see `prefer-flag-not-drop-genes`), rather than EVM's single
//! combined trellis that picks one strand per antisense locus.
//!
//! Recursion mirrors `generate_consensus_gene_predictions` (~483): after the best path,
//! re-search the left/right tails beyond the *retained* genes' span (≥ `research_size`),
//! which is also the fix for EVM's bug where an eliminated gene's span blocked re-search.

use crate::consensus::candidates::{
    CandidateParams, MAX_TERMINAL_EXTEND, RegionData, build_candidates,
};
use crate::consensus::evidence::EvidenceChain;
use crate::consensus::exon::{ExonType, frame_offset};
use crate::consensus::filter::{FilterParams, SupportFlags, assess};
use crate::consensus::region::ConsensusRegion;
use crate::consensus::trellis::{ConsensusGene, run_trellis, score_all_exons};
use crate::consensus::weights::Weights;
use crate::io::fasta::Fasta;
use crate::model::{Coordset, Strand};
use crate::orf::GeneticCode;
use crate::orf::coords::{OrfProjection, SplicedTranscript};
use crate::orf::translate::reverse_complement;

/// Trellis/recursion knobs.
#[derive(Clone, Copy, Debug)]
pub struct EngineParams {
    pub max_prev_exons: usize,
    /// Minimum tail size (bp) to re-search (EVM `terminal_intergenic_re_search`, 10000).
    pub research_size: i64,
    /// Minimum intergenic-gap size (bp) between called genes to re-search for additional
    /// genes (EVM `--re_search_intergenic`). 0 = off (the EVM default).
    pub research_intergenic: i64,
    /// Minimum intron length (bp) to re-search for nested genes (EVM
    /// `--search_long_introns`). 0 = off (the EVM default).
    pub search_long_introns: i64,
    pub intergenic_adjust: f64,
    /// The CDS of a called gene that begins at an initial/single exon must begin with ATG:
    /// trim in-frame to the first ATG (5' overhang becomes 5'UTR), else extend the first
    /// exon upstream to one, else flag the gene `partial5` (`--require-atg-start`).
    pub require_atg_start: bool,
    /// Extend an unterminated ORF's 3'-most exon downstream to the first in-frame stop
    /// when the path claims a terminal/single exon (`--extend-orf-stop`; EVM
    /// `extend_to_terminal`, applied after the trellis).
    pub extend_orf_stop: bool,
}

/// A consensus gene resolved to forward genomic coordinates, with its coding region
/// projected (CDS + flanking UTR).
#[derive(Debug, Clone, PartialEq)]
pub struct CalledGene {
    pub contig: String,
    pub orient: Strand,
    pub exons: Vec<Coordset>,
    pub cds: Vec<Coordset>,
    pub five_utr: Vec<Coordset>,
    pub three_utr: Vec<Coordset>,
    pub partial5: bool,
    pub partial3: bool,
    /// The ORF reached an in-frame stop followed by at least one more complete codon of
    /// the path, so the CDS stops short and the rest of the path is 3'UTR (a premature
    /// stop; 1-2 stray bases after a stop are not a lost codon). Independent of
    /// `partial3`, which keeps its meaning (no stop codon claimed or found).
    pub truncated: bool,
    /// GFF3 phase of the 5'-most CDS base: the number of bases before the first complete
    /// codon. Nonzero only for a 5'-partial gene, whose CDS starts at its first exon base
    /// (the partial codon is coding, not UTR); 0 for a gene with a start codon.
    pub cds_start_phase: u8,
    pub score: f64,
    pub support: SupportFlags,
    /// True for a transcript-ORF gene recovered at a locus with no consensus CDS
    /// (`--promote-transcript-orfs`), vs. a trellis-derived consensus gene.
    pub promoted: bool,
    /// Per-exon / per-intron evidence attribution for the evidence report, in ascending
    /// forward `lend` order.
    pub features: Vec<FeatureSupport>,
}

/// One exon or intron of a called gene with the evidence supporting it (forward genomic
/// coordinates). Rendered by [`crate::consensus::report`].
#[derive(Debug, Clone, PartialEq)]
pub struct FeatureSupport {
    pub coords: Coordset,
    pub kind: FeatureKind,
    /// `(accession, ev_type)` of each supporting evidence chain, in insertion order.
    pub evidence: Vec<(String, String)>,
}

/// An exon (typed, with EVM frames: 1-3 forward, 4-6 reverse) or an intron.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FeatureKind {
    Exon {
        exon_type: ExonType,
        start_frame: u8,
        end_frame: u8,
    },
    Intron,
}

impl CalledGene {
    /// Genomic span `(min lend, max rend)` over the gene's exons. Robust to ordering
    /// (the exons are lend-sorted in practice, but this does not rely on it).
    pub fn span(&self) -> (i64, i64) {
        let lend = self.exons.iter().map(|c| c.lend).min().unwrap_or(0);
        let rend = self.exons.iter().map(|c| c.rend).max().unwrap_or(0);
        (lend, rend)
    }
}

const MAX_RECURSION_DEPTH: usize = 32;

/// `--require-atg-start`: how far (in codons) past the projected CDS start the in-frame
/// ATG rescue may scan before the gene is flagged `partial5` instead.
const MAX_ATG_SCAN_CODONS: usize = 120;

/// Max bp to extend a path's 5'-most exon upstream looking for a stop-bounded ATG
/// (the evidence start is truncated); the extended 5'UTR is trimmed to just past the
/// bounding stop.
const MAX_START_EXTEND: i64 = 5_000;

/// First in-frame ATG at or after `t0`, within `MAX_ATG_SCAN_CODONS` codons, leaving at
/// least one codon before the ORF end `t_end`. `seq` is the spliced transcript 5'->3'.
fn first_in_frame_atg(seq: &[u8], t0: usize, t_end: usize) -> Option<usize> {
    let limit = (t0 + MAX_ATG_SCAN_CODONS * 3).min(t_end);
    let mut t = t0;
    while t + 3 < t_end && t < limit {
        if seq[t..t + 3].eq_ignore_ascii_case(b"ATG") {
            return Some(t);
        }
        t += 3;
    }
    None
}

/// Extend the 3'-most exon of a path downstream (in transcription orientation) by up to
/// MAX_TERMINAL_EXTEND bp and re-project the ORF; returns the extended exon set and the
/// new projection iff an in-frame stop is now reached (EVM `extend_to_terminal`).
fn extend_path_to_stop(
    exons: &[Coordset],
    orient: Strand,
    genome: &Fasta,
    contig: &str,
    cds_t_start: usize,
    code: &GeneticCode,
) -> Option<(Vec<Coordset>, OrfProjection)> {
    let mut ex = exons.to_vec(); // lend-sorted ascending
    let contig_len = genome.contig_len(contig)? as i64;
    if orient == Strand::Minus {
        // the 3' end of a minus-strand transcript is the smallest coordinate
        let first = ex.first_mut().unwrap();
        first.lend = (first.lend - MAX_TERMINAL_EXTEND).max(1);
    } else {
        let last = ex.last_mut().unwrap();
        last.rend = (last.rend + MAX_TERMINAL_EXTEND).min(contig_len);
    }
    let st = SplicedTranscript::new(&ex, orient);
    let seq = st.sequence(genome, contig)?;
    let proj = st.project_orf(&seq, cds_t_start.min(seq.len()), code);
    proj.hit_stop.then_some((ex, proj))
}

/// Extend the 5'-most exon of a path upstream by up to MAX_START_EXTEND bp and search it
/// for an in-frame ATG bounded by the nearest upstream in-frame stop (the maximal ORF
/// extension, for evidence whose 5' start is truncated). The retained 5'UTR is trimmed
/// to just past the bounding stop so the gene does not swallow intergenic sequence.
/// Returns the extended exons, transcript and the ORF re-projected
/// from the ATG, iff such an ATG exists.
fn extend_path_to_atg(
    exons: &[Coordset],
    orient: Strand,
    genome: &Fasta,
    contig: &str,
    cds_t_start: usize, // offset in the ORIGINAL spliced transcript
    code: &GeneticCode,
) -> Option<(Vec<Coordset>, SplicedTranscript, OrfProjection)> {
    let mut ex = exons.to_vec(); // lend-sorted ascending
    let contig_len = genome.contig_len(contig)? as i64;
    // extend the 5' end (smallest coordinate on '+', largest on '-')
    let shift;
    if orient == Strand::Minus {
        let last = ex.last_mut().unwrap();
        let old = last.rend;
        last.rend = (last.rend + MAX_START_EXTEND).min(contig_len);
        shift = (last.rend - old) as usize;
    } else {
        let first = ex.first_mut().unwrap();
        let old = first.lend;
        first.lend = (first.lend - MAX_START_EXTEND).max(1);
        shift = (old - first.lend) as usize;
    }
    if shift == 0 {
        return None;
    }
    let st = SplicedTranscript::new(&ex, orient);
    let seq = st.sequence(genome, contig)?;
    // original CDS start in the extended transcript; the in-frame codon grid is
    // { t : (start5 - t) % 3 == 0 }
    let start5 = shift + cds_t_start;
    // walk upstream in-frame to the bounding stop codon
    let mut stop_t: Option<usize> = None;
    let mut t = start5 as i64 - 3;
    while t >= 0 {
        if code.is_stop(&seq[t as usize..t as usize + 3]) {
            stop_t = Some(t as usize);
            break;
        }
        t -= 3;
    }
    // trim the 5'UTR to just past the bounding stop
    if let Some(ts) = stop_t {
        let trim = ts + 3;
        if trim >= shift {
            return None; // stop inside the original first exon: no room to rescue
        }
        if orient == Strand::Minus {
            ex.last_mut().unwrap().rend -= trim as i64;
        } else {
            ex.first_mut().unwrap().lend += trim as i64;
        }
        let st = SplicedTranscript::new(&ex, orient);
        let seq = st.sequence(genome, contig)?;
        return finish_atg_search(ex, st, seq, start5 - trim, code);
    }
    finish_atg_search(ex, st, seq, start5, code)
}
/// Scan the trimmed 5' extension of `seq` for the first in-frame ATG strictly upstream
/// of the original CDS start `start5`, and re-project the ORF from it.
fn finish_atg_search(
    ex: Vec<Coordset>,
    st: SplicedTranscript,
    seq: Vec<u8>,
    start5: usize,
    code: &GeneticCode,
) -> Option<(Vec<Coordset>, SplicedTranscript, OrfProjection)> {
    let mut t = start5 % 3;
    while t < start5 {
        if seq[t..t + 3].eq_ignore_ascii_case(b"ATG") {
            let proj = st.project_orf(&seq, t, code);
            // reject only an ATG with no complete codon after it in the transcript; a stop
            // counts, so an ATG immediately followed by a stop is accepted
            if proj.cds_t_end > t + 3 {
                return Some((ex, st, proj));
            }
            return None;
        }
        t += 3;
    }
    None
}

/// Find consensus genes in `region` on both strands.
#[allow(clippy::too_many_arguments)]
pub fn consensus_region(
    region: &ConsensusRegion,
    chains: &[EvidenceChain],
    genome: &Fasta,
    cand: &CandidateParams,
    filt: &FilterParams,
    weights: &Weights,
    eng: &EngineParams,
    mask: &[Coordset],
) -> Vec<CalledGene> {
    // clamp the region span to the contig
    let contig_len = genome
        .contig_len(&region.contig)
        .map(|n| n as i64)
        .unwrap_or(region.span.rend);
    let lo = region.span.lend.max(1);
    let hi = region.span.rend.min(contig_len);
    if hi < lo {
        return Vec::new();
    }
    let region = ConsensusRegion {
        contig: region.contig.clone(),
        span: Coordset::new(lo, hi),
        chain_indices: region.chain_indices.clone(),
    };

    let mut out = Vec::new();

    // ---- forward strand ----
    let fwd = build_candidates(&region, chains, genome, cand, mask);
    let fwd_base = score_all_exons(&fwd.exons, &fwd.vectors, weights);
    let mut fwd_genes = Vec::new();
    collect_genes(
        &fwd,
        &fwd_base,
        (lo, hi),
        &cand.code,
        filt,
        eng,
        &mut fwd_genes,
        0,
    );
    resolve_all(
        &fwd_genes,
        &fwd,
        &region.contig,
        Strand::Plus,
        None,
        genome,
        &cand.code,
        eng.require_atg_start,
        eng.extend_orf_stop,
        &mut out,
    );

    // ---- minus strand: reverse-complement the region, run forward, transpose back ----
    if let Some((rc_fasta, rc_region, rc_chains, rc_mask)) =
        build_rc_region(&region, chains, mask, genome, lo, hi)
    {
        let len = hi - lo + 1;
        let rc = build_candidates(&rc_region, &rc_chains, &rc_fasta, cand, &rc_mask);
        let rc_base = score_all_exons(&rc.exons, &rc.vectors, weights);
        let mut rc_genes = Vec::new();
        collect_genes(
            &rc,
            &rc_base,
            (1, len),
            &cand.code,
            filt,
            eng,
            &mut rc_genes,
            0,
        );
        resolve_all(
            &rc_genes,
            &rc,
            &region.contig,
            Strand::Minus,
            Some(hi),
            genome,
            &cand.code,
            eng.require_atg_start,
            eng.extend_orf_stop,
            &mut out,
        );
    }

    dedup_genes(out)
}

/// Reverse-complement transform for the minus strand: build the RC genome window over
/// `[lo, hi]`, transpose every minus-strand chain and the repeat mask into RC-local
/// coordinates (via [`Coordset::reflect`]), and return the RC region ready to run through
/// the same forward machinery. `None` when the genome window is unavailable.
fn build_rc_region(
    region: &ConsensusRegion,
    chains: &[EvidenceChain],
    mask: &[Coordset],
    genome: &Fasta,
    lo: i64,
    hi: i64,
) -> Option<(Fasta, ConsensusRegion, Vec<EvidenceChain>, Vec<Coordset>)> {
    let fwd_bytes = genome.subseq(&region.contig, lo, hi)?;
    let len = hi - lo + 1;
    let rc_fasta = Fasta::from_seq("rc", reverse_complement(fwd_bytes));

    let mut rc_chains = Vec::new();
    for &ci in &region.chain_indices {
        let c = &chains[ci];
        if c.orient != Strand::Minus {
            continue;
        }
        let mut links: Vec<Coordset> = c.links.iter().map(|s| s.reflect(hi)).collect();
        links.sort_by_key(|c| c.lend);
        let span = Coordset::new(links.first().unwrap().lend, links.last().unwrap().rend);
        rc_chains.push(EvidenceChain {
            accession: c.accession.clone(),
            ev_type: c.ev_type.clone(),
            ev_class: c.ev_class,
            weight: c.weight,
            contig: "rc".into(),
            orient: Strand::Plus,
            span,
            links,
            cds_start_phase: c.cds_start_phase,
        });
    }

    let rc_region = ConsensusRegion {
        contig: "rc".into(),
        span: Coordset::new(1, len),
        chain_indices: (0..rc_chains.len()).collect(),
    };
    let rc_mask: Vec<Coordset> = mask.iter().map(|iv| iv.reflect(hi)).collect();
    Some((rc_fasta, rc_region, rc_chains, rc_mask))
}

/// Resolve a strand's trellis genes to forward coordinates and append them to `out`.
#[allow(clippy::too_many_arguments)]
fn resolve_all(
    genes: &[(ConsensusGene, SupportFlags)],
    rd: &RegionData,
    contig: &str,
    orient: Strand,
    transpose_hi: Option<i64>,
    genome: &Fasta,
    code: &GeneticCode,
    require_atg: bool,
    extend_orf_stop: bool,
    out: &mut Vec<CalledGene>,
) {
    for (g, support) in genes {
        out.push(resolve(
            g,
            rd,
            contig,
            orient,
            transpose_hi,
            *support,
            genome,
            code,
            require_atg,
            extend_orf_stop,
        ));
    }
}

/// Drop genes with an identical (orientation, exon structure) — re-search can rediscover
/// the same gene via overlapping tail/gap/intron sub-ranges.
fn dedup_genes(mut genes: Vec<CalledGene>) -> Vec<CalledGene> {
    let mut seen = std::collections::HashSet::new();
    genes.retain(|g| seen.insert((g.orient, g.exons.clone())));
    genes
}

/// Resolve a trellis gene to forward genomic coordinates (transposing RC-local minus
/// coordinates when `transpose_hi` is set) and project its CDS/UTR.
///
/// The consensus exons ARE the coding structure: build a [`SplicedTranscript`] over them
/// and project the ORF from the gene's 5' start, in its reading frame, to the first
/// in-frame stop (reusing `orf::project_orf` — decision 6). For a clean model the CDS
/// equals the exon structure; a non-canonical terminus without a real stop is 3'-partial,
/// unless `extend_orf_stop` finds the stop downstream of a claimed terminal/single exon.
#[allow(clippy::too_many_arguments)]
fn resolve(
    g: &ConsensusGene,
    rd: &RegionData,
    contig: &str,
    orient: Strand,
    transpose_hi: Option<i64>,
    mut support: SupportFlags,
    genome: &Fasta,
    code: &GeneticCode,
    require_atg: bool,
    extend_orf_stop: bool,
) -> CalledGene {
    let mut exons: Vec<Coordset> = g
        .exon_indices
        .iter()
        .map(|&i| {
            let c = rd.exons[i].coords;
            match transpose_hi {
                Some(hi) => c.reflect(hi),
                None => c,
            }
        })
        .collect();
    exons.sort_by_key(|c| c.lend);

    // 5'/3'-most exons in the gene's reading order (path order is 5'→3' in each pass).
    let first = &rd.exons[g.exon_indices[0]];
    let last = &rd.exons[*g.exon_indices.last().unwrap()];
    // A start codon needs an initial/single exon that begins on a codon (frame 1). One
    // that begins mid-codon (a prediction whose 5'-most CDS row has a nonzero GFF3
    // phase) is a 5'-partial start, like a path that begins at an internal exon.
    let has_start =
        matches!(first.exon_type, ExonType::Initial | ExonType::Single) && first.start_frame == 1;
    let has_stop = matches!(last.exon_type, ExonType::Terminal | ExonType::Single);
    // first complete codon within the spliced transcript: a 5'-partial path begins with
    // a partial codon of `p` bases (0 when its first exon begins on a codon, frame 1)
    let p = frame_offset(first.start_frame);

    // set when --extend-orf-stop extends the 3'-most exon to a downstream stop
    let mut extended3 = false;
    // --require-atg-start may repair a gene's start or flag it 5'-partial
    let mut partial5 = !has_start;
    let st = SplicedTranscript::new(&exons, orient);
    let (cds, five_utr, three_utr, hit_stop, cds_start_phase) = match st.sequence(genome, contig) {
        Some(seq) => {
            // the ORF (and its stop) is read from the first complete codon ...
            let mut st = st;
            let mut proj = st.project_orf(&seq, p.min(seq.len()), code);
            // --require-atg-start: the CDS must begin with ATG. Trim in-frame to the
            // first ATG (the 5' overhang becomes 5'UTR); when the path holds no
            // in-frame ATG (truncated evidence start), extend the 5'-most exon
            // upstream to a stop-bounded ATG; otherwise flag the gene 5'-partial.
            // Only a start codon is trimmed: a first exon that begins mid-codon is only
            // extended, and a path that begins at an internal exon is left as it is.
            let mut atg_ok = true;
            let mut moved5 = false;
            if require_atg {
                let atg = if has_start {
                    first_in_frame_atg(&seq, proj.cds_t_start, proj.cds_t_end)
                } else {
                    None
                };
                match atg {
                    Some(t) => {
                        moved5 = t != proj.cds_t_start;
                        proj.cds_t_start = t;
                    }
                    None => {
                        if matches!(first.exon_type, ExonType::Initial | ExonType::Single)
                            && let Some((ex5, st5, proj5)) = extend_path_to_atg(
                                &exons,
                                orient,
                                genome,
                                contig,
                                proj.cds_t_start,
                                code,
                            )
                        {
                            exons = ex5;
                            st = st5;
                            proj = proj5;
                            moved5 = true;
                        } else {
                            atg_ok = false;
                        }
                    }
                }
                partial5 = !atg_ok;
            }
            // --extend-orf-stop: an ORF that never reaches a stop (3'-partial evidence)
            // extends its 3'-most exon downstream to the first in-frame stop, capped at
            // MAX_TERMINAL_EXTEND (EVM `extend_to_terminal`), only when the path claims a
            // terminal/single exon; a path ending in an internal exon stays 3'-partial
            if extend_orf_stop
                && has_stop
                && !proj.hit_stop
                && let Some((ex2, proj2)) =
                    extend_path_to_stop(&exons, orient, genome, contig, proj.cds_t_start, code)
            {
                exons = ex2;
                proj = proj2;
                extended3 = true;
            }
            // rebuild the transcript when the 3' exon was extended
            let st = if extended3 {
                SplicedTranscript::new(&exons, orient)
            } else {
                st
            };
            // ... but a 5'-partial CDS starts at the transcript's first base, the leading
            // partial codon carried as the first CDS row's phase (GFF3), not as a 5'UTR
            let (cds_t_start, phase) =
                if !has_start && proj.cds_t_end > proj.cds_t_start && !(require_atg && atg_ok) {
                    (0, p as u8)
                } else {
                    (proj.cds_t_start, 0)
                };
            let (cds, five_utr, mut three_utr) = st.cds_and_utrs(cds_t_start, proj.cds_t_end);
            if extended3 {
                // the stop lies in the cap-extended 3'-most exon (the unextended path held
                // none): cut that exon back to the end of the CDS, leaving no 3'UTR
                if orient == Strand::Minus {
                    exons[0].lend = cds[0].lend;
                } else {
                    exons.last_mut().unwrap().rend = cds.last().unwrap().rend;
                }
                three_utr = Vec::new();
            }
            if extended3 || moved5 {
                support.coding_length = cds.iter().map(|c| c.len()).sum();
            }
            (cds, five_utr, three_utr, proj.hit_stop, phase)
        }
        None => {
            // sequence unavailable (e.g. evidence past the contig end): emit the raw
            // path exons but never claim completeness — no stop scan was possible.
            let phase = if has_start { 0 } else { p as u8 };
            if require_atg {
                partial5 = true; // nor can the start codon be checked
            }
            (exons.clone(), Vec::new(), Vec::new(), false, phase)
        }
    };

    // a stop followed by at least one more complete codon of the path cut the CDS short
    // (premature); 1-2 trailing bases after the stop are not a lost codon
    let truncated = hit_stop && three_utr.iter().map(|c| c.len()).sum::<i64>() >= 3;
    CalledGene {
        contig: contig.into(),
        orient,
        exons,
        cds,
        five_utr,
        three_utr,
        partial5,
        // a path rescued by the 3' extension terminates at a real stop codon
        partial3: !(has_stop || extended3) || !hit_stop,
        truncated,
        cds_start_phase,
        score: g.score,
        support,
        promoted: false,
        features: path_features(g, rd, transpose_hi),
    }
}

/// The evidence attribution of a trellis gene's path: each path exon (its candidate's
/// type, frames and evidence) and each intron between consecutive path exons. Intron
/// evidence is looked up in the pass's WORKING coordinates (where `rd.introns` is keyed),
/// and only then are the coordinates transposed to forward space for the minus pass
/// (whose frames become EVM's reverse 4-6).
fn path_features(
    g: &ConsensusGene,
    rd: &RegionData,
    transpose_hi: Option<i64>,
) -> Vec<FeatureSupport> {
    let to_fwd = |c: Coordset| transpose_hi.map_or(c, |hi| c.reflect(hi));
    let frame_shift = if transpose_hi.is_some() { 3 } else { 0 };
    let mut features = Vec::with_capacity(2 * g.exon_indices.len());
    for (k, &i) in g.exon_indices.iter().enumerate() {
        let ex = &rd.exons[i];
        if k > 0 {
            let gap = rd.exons[g.exon_indices[k - 1]].coords.gap_to(&ex.coords);
            features.push(FeatureSupport {
                coords: to_fwd(gap),
                kind: FeatureKind::Intron,
                evidence: rd.introns.evidence((gap.lend, gap.rend)).to_vec(),
            });
        }
        features.push(FeatureSupport {
            coords: to_fwd(ex.coords),
            kind: FeatureKind::Exon {
                exon_type: ex.exon_type,
                start_frame: ex.start_frame + frame_shift,
                end_frame: ex.end_frame + frame_shift,
            },
            evidence: ex.evidence.clone(),
        });
    }
    features.sort_by_key(|f| f.coords.lend);
    features
}

/// Build the best path over `range`, flag low support, keep the retained genes, then
/// re-search the tails beyond their span. Flag-not-drop keeps everything by default;
/// `--strict` drops low-support genes, and only retained spans block re-search.
#[allow(clippy::too_many_arguments)]
fn collect_genes(
    rd: &RegionData,
    base: &[f64],
    range: (i64, i64),
    code: &GeneticCode,
    filt: &FilterParams,
    eng: &EngineParams,
    out: &mut Vec<(ConsensusGene, SupportFlags)>,
    depth: usize,
) {
    if depth >= MAX_RECURSION_DEPTH {
        return;
    }
    let genes = run_trellis(
        &rd.exons,
        base,
        range,
        &rd.introns,
        &rd.vectors,
        code,
        eng.intergenic_adjust,
        eng.max_prev_exons,
    );
    if genes.is_empty() {
        return;
    }

    let mut retained: Vec<(ConsensusGene, SupportFlags)> = Vec::new();
    for g in genes {
        let support = assess(&g, &rd.exons, &rd.vectors, filt);
        if filt.strict && support.low_support {
            continue; // EVM-style drop
        }
        retained.push((g, support));
    }
    if retained.is_empty() {
        return;
    }

    // Per-gene spans and (optionally) long introns, captured before `retained` moves.
    let mut gene_spans: Vec<(i64, i64)> = retained
        .iter()
        .map(|(g, _)| {
            let l = g
                .exon_indices
                .iter()
                .map(|&i| rd.exons[i].coords.lend)
                .min()
                .unwrap();
            let r = g
                .exon_indices
                .iter()
                .map(|&i| rd.exons[i].coords.rend)
                .max()
                .unwrap();
            (l, r)
        })
        .collect();
    let long_introns: Vec<(i64, i64)> = if eng.search_long_introns > 0 {
        retained
            .iter()
            .flat_map(|(g, _)| {
                g.exon_indices.windows(2).filter_map(|w| {
                    let intron = rd.exons[w[0]].coords.gap_to(&rd.exons[w[1]].coords);
                    (intron.len() >= eng.search_long_introns).then_some((intron.lend, intron.rend))
                })
            })
            .collect()
    } else {
        Vec::new()
    };

    let span_lend = gene_spans.iter().map(|s| s.0).min().unwrap();
    let span_rend = gene_spans.iter().map(|s| s.1).max().unwrap();
    out.extend(retained);

    // left / right tail re-search (EVM default behaviour, gated on `research_size`)
    if span_lend - range.0 >= eng.research_size {
        collect_genes(
            rd,
            base,
            (range.0, span_lend - 1),
            code,
            filt,
            eng,
            out,
            depth + 1,
        );
    }
    if range.1 - span_rend >= eng.research_size {
        collect_genes(
            rd,
            base,
            (span_rend + 1, range.1),
            code,
            filt,
            eng,
            out,
            depth + 1,
        );
    }

    // intergenic-gap re-search between called genes (EVM `--re_search_intergenic`, opt-in)
    if eng.research_intergenic > 0 {
        gene_spans.sort_unstable();
        for w in gene_spans.windows(2) {
            let (gap_lend, gap_rend) = (w[0].1 + 1, w[1].0 - 1);
            if gap_rend - gap_lend + 1 >= eng.research_intergenic {
                collect_genes(
                    rd,
                    base,
                    (gap_lend, gap_rend),
                    code,
                    filt,
                    eng,
                    out,
                    depth + 1,
                );
            }
        }
    }

    // long-intron re-search for nested genes (EVM `--search_long_introns`, opt-in)
    for &(il, ir) in &long_introns {
        collect_genes(rd, base, (il, ir), code, filt, eng, out, depth + 1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consensus::weights::EvClass;
    use std::io::Write;

    fn write_fasta(seq: &[u8]) -> Fasta {
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        writeln!(tmp, ">chr1").unwrap();
        tmp.write_all(seq).unwrap();
        writeln!(tmp).unwrap();
        Fasta::load(tmp.path()).unwrap()
    }

    fn chain(class: EvClass, weight: f64, orient: Strand, segs: &[(i64, i64)]) -> EvidenceChain {
        let links: Vec<Coordset> = segs.iter().map(|&(l, r)| Coordset::new(l, r)).collect();
        let span = Coordset::new(links.first().unwrap().lend, links.last().unwrap().rend);
        EvidenceChain {
            accession: "g".into(),
            ev_type: "fgenesh".into(),
            ev_class: class,
            weight,
            contig: "chr1".into(),
            orient,
            span,
            links,
            cds_start_phase: 0,
        }
    }

    fn config() -> (CandidateParams, FilterParams, Weights, EngineParams) {
        (
            CandidateParams {
                min_intron_length: 20,
                code: GeneticCode::default(),
                extend_terminal_stop: false,
                peak_augment: false,
                intron_score_scale: crate::consensus::IntronScoreScale::Linear,
                sum_genepred_weights: 0.0,
            },
            FilterParams {
                min_score_ratio: 0.75,
                min_coding_length: 150,
                intergenic_adjust: 1.0,
                strict: false,
            },
            Weights::parse_str("ABINITIO_PREDICTION fgenesh 1\n", "w").unwrap(),
            EngineParams {
                max_prev_exons: 500,
                research_size: 10_000,
                research_intergenic: 0,
                search_long_introns: 0,
                intergenic_adjust: 1.0,
                require_atg_start: false,
                extend_orf_stop: false,
            },
        )
    }

    #[test]
    fn forward_minus_strand_gene_is_transposed_back() {
        // A 200-bp genome of 'C' with a minus-strand single-exon prediction 50..120.
        // On the minus strand the coding reads revcomp(C-run)=G-run -> no stop, so the
        // single exon survives and must come back at its forward coordinates with '-'.
        let g = write_fasta(&[b'C'; 200]);
        let chains = vec![chain(
            EvClass::AbinitioPrediction,
            1.0,
            Strand::Minus,
            &[(50, 120)],
        )];
        let region = ConsensusRegion {
            contig: "chr1".into(),
            span: Coordset::new(1, 200),
            chain_indices: vec![0],
        };
        let (cand, filt, w, eng) = config();
        let genes = consensus_region(&region, &chains, &g, &cand, &filt, &w, &eng, &[]);
        let minus: Vec<_> = genes.iter().filter(|g| g.orient == Strand::Minus).collect();
        assert_eq!(minus.len(), 1, "the minus-strand gene is called");
        assert_eq!(
            minus[0].exons,
            vec![Coordset::new(50, 120)],
            "transposed back to forward coords"
        );
    }

    #[test]
    fn forward_strand_gene_called_in_place() {
        let g = write_fasta(&[b'C'; 200]);
        // forward single-exon prediction; revcomp not needed
        let chains = vec![chain(
            EvClass::AbinitioPrediction,
            1.0,
            Strand::Plus,
            &[(50, 120)],
        )];
        let region = ConsensusRegion {
            contig: "chr1".into(),
            span: Coordset::new(1, 200),
            chain_indices: vec![0],
        };
        let (cand, filt, w, eng) = config();
        let genes = consensus_region(&region, &chains, &g, &cand, &filt, &w, &eng, &[]);
        let plus: Vec<_> = genes.iter().filter(|g| g.orient == Strand::Plus).collect();
        assert_eq!(plus.len(), 1);
        assert_eq!(plus[0].exons, vec![Coordset::new(50, 120)]);
    }

    fn genome_with(len: usize, edits: &[(i64, &[u8])]) -> Fasta {
        let mut seq = vec![b'C'; len];
        for &(pos1, bytes) in edits {
            for (k, &b) in bytes.iter().enumerate() {
                seq[(pos1 - 1) as usize + k] = b;
            }
        }
        write_fasta(&seq)
    }

    fn call_single(g: &Fasta, segs: &[(i64, i64)], orient: Strand) -> Option<CalledGene> {
        let chains = vec![chain(EvClass::AbinitioPrediction, 1.0, orient, segs)];
        let region = ConsensusRegion {
            contig: "chr1".into(),
            span: Coordset::new(1, 200),
            chain_indices: vec![0],
        };
        let (cand, filt, w, eng) = config();
        consensus_region(&region, &chains, g, &cand, &filt, &w, &eng, &[])
            .into_iter()
            .find(|x| x.orient == orient)
    }

    fn ev_g() -> Vec<(String, String)> {
        vec![("g".to_string(), "fgenesh".to_string())]
    }

    fn exon_feat(l: i64, r: i64, t: ExonType, sf: u8, ef: u8) -> FeatureSupport {
        FeatureSupport {
            coords: Coordset::new(l, r),
            kind: FeatureKind::Exon {
                exon_type: t,
                start_frame: sf,
                end_frame: ef,
            },
            evidence: ev_g(),
        }
    }

    fn intron_feat(l: i64, r: i64) -> FeatureSupport {
        FeatureSupport {
            coords: Coordset::new(l, r),
            kind: FeatureKind::Intron,
            evidence: ev_g(),
        }
    }

    #[test]
    fn plus_gene_features_carry_exon_and_intron_evidence() {
        let g = genome_with(200, &[]);
        let gene = call_single(&g, &[(10, 18), (50, 60)], Strand::Plus).unwrap();
        assert_eq!(
            gene.features,
            vec![
                exon_feat(10, 18, ExonType::Initial, 1, 3),
                intron_feat(19, 49),
                exon_feat(50, 60, ExonType::Terminal, 1, 2),
            ]
        );
    }

    #[test]
    fn minus_gene_features_transposed_with_reverse_frames() {
        // Forward exons A=50..58 (3', terminal) and B=89..98 (5', initial), intron 59..88.
        // In the RC-local pass (hi = 200) the intron is keyed 113..142, so a lookup with
        // the forward coordinates would find no evidence. Frames: B initial 1->1 (len 10),
        // A terminal starts at cum 10 -> 2, ends 1; reverse encoding adds 3.
        let g = genome_with(200, &[]);
        let gene = call_single(&g, &[(50, 58), (89, 98)], Strand::Minus).unwrap();
        assert_eq!(
            gene.exons,
            vec![Coordset::new(50, 58), Coordset::new(89, 98)]
        );
        assert_eq!(
            gene.features,
            vec![
                exon_feat(50, 58, ExonType::Terminal, 5, 4),
                intron_feat(59, 88),
                exon_feat(89, 98, ExonType::Initial, 4, 4),
            ]
        );
    }

    #[test]
    fn clean_forward_orf_yields_cds_no_utr() {
        // 10..21 = ATG AAA CCC TAA (stop at 19..21)
        let g = genome_with(40, &[(10, b"ATGAAACCCTAA")]);
        let gene = call_single(&g, &[(10, 21)], Strand::Plus).unwrap();
        assert_eq!(gene.cds, vec![Coordset::new(10, 21)]);
        assert!(gene.five_utr.is_empty() && gene.three_utr.is_empty());
        assert!(!gene.partial5 && !gene.partial3);
        assert_eq!(
            gene.cds_start_phase, 0,
            "a gene with a start codon has phase 0"
        );
    }

    #[test]
    fn no_stop_codon_is_three_prime_partial() {
        // all-C exon: translation never hits a stop
        let g = genome_with(40, &[]);
        let gene = call_single(&g, &[(10, 21)], Strand::Plus).unwrap();
        assert!(gene.partial3, "no in-frame stop -> 3' partial");
    }

    #[test]
    fn minus_gene_cds_projected_in_forward_coords() {
        // minus single-exon, length a multiple of 3 (50..118 = 69 bp) so no remainder
        let g = genome_with(200, &[]);
        let gene = call_single(&g, &[(50, 118)], Strand::Minus).unwrap();
        assert_eq!(gene.orient, Strand::Minus);
        assert_eq!(gene.cds, vec![Coordset::new(50, 118)]);
    }

    #[test]
    fn long_intron_research_finds_nested_gene() {
        // Gene A: two CDS exons across a ~34 kb intron (its intron score dominates, so
        // the global best path is A alone). Gene B is a single exon nested inside that
        // intron, unreachable from A's exons by the grammar.
        let g = genome_with(40_000, &[]);
        let a = chain(
            EvClass::AbinitioPrediction,
            1.0,
            Strand::Plus,
            &[(1000, 1100), (35000, 35100)],
        );
        let b = chain(
            EvClass::AbinitioPrediction,
            1.0,
            Strand::Plus,
            &[(10000, 10300)],
        );
        let chains = vec![a, b];
        let region = ConsensusRegion {
            contig: "chr1".into(),
            span: Coordset::new(1, 40_000),
            chain_indices: vec![0, 1],
        };
        let (cand, filt, w, mut eng) = config();

        let has_b = |genes: &[CalledGene]| {
            genes
                .iter()
                .any(|x| x.exons.iter().any(|c| c.lend == 10000))
        };

        // off (default): B is hidden inside A's long intron
        let off = consensus_region(&region, &chains, &g, &cand, &filt, &w, &eng, &[]);
        assert!(
            !has_b(&off),
            "nested gene hidden when long-intron re-search is off"
        );

        // on: re-search introns >= 20 kb finds B
        eng.search_long_introns = 20_000;
        let on = consensus_region(&region, &chains, &g, &cand, &filt, &w, &eng, &[]);
        assert!(has_b(&on), "long-intron re-search recovers the nested gene");
    }

    /// A hand-built region (working coordinates) holding the given path exons, and the
    /// trellis gene that walks them in order.
    fn partial_path(len: usize, path: &[(i64, i64, ExonType, u8)]) -> (RegionData, ConsensusGene) {
        use crate::consensus::exon::{ExonCandidate, end_frame_for};
        use crate::consensus::vectors::{IntronScores, RegionVectors};
        let exons = path
            .iter()
            .map(|&(l, r, t, sf)| ExonCandidate {
                coords: Coordset::new(l, r),
                orient: Strand::Plus,
                exon_type: t,
                start_frame: sf,
                end_frame: end_frame_for(sf, r - l + 1),
                left_seq_boundary: *b"CC",
                right_seq_boundary: *b"CC",
                evidence: ev_g(),
            })
            .collect();
        let rd = RegionData {
            contig: "chr1".into(),
            span: Coordset::new(1, len as i64),
            exons,
            vectors: RegionVectors::new(1, len),
            introns: IntronScores::default(),
        };
        let g = ConsensusGene {
            exon_indices: (0..path.len()).collect(),
            orient: Strand::Plus,
            score: 1.0,
        };
        (rd, g)
    }

    fn resolve_partial(
        genome: &Fasta,
        len: usize,
        path: &[(i64, i64, ExonType, u8)],
        orient: Strand,
        transpose_hi: Option<i64>,
        require_atg: bool,
        extend_orf_stop: bool,
    ) -> CalledGene {
        let (rd, g) = partial_path(len, path);
        let support = SupportFlags {
            raw_noncoding: 0.0,
            noncoding_equivalent: 0.0,
            score_ratio: 0.0,
            coding_length: 0,
            low_support: false,
        };
        resolve(
            &g,
            &rd,
            "chr1",
            orient,
            transpose_hi,
            support,
            genome,
            &GeneticCode::default(),
            require_atg,
            extend_orf_stop,
        )
    }

    #[test]
    fn five_prime_partial_frame2_starts_cds_at_exon_with_phase_2() {
        // Internal 10..30 (start frame 2: its first base is codon position 2) -> Terminal
        // 50..60. The next codon starts 2 bases in, at 12: AGC CCC CCC TAA (stop 21..23).
        // Reading 1 base in (from 11) would hit TAG at 11..13 at once. The CDS starts at
        // the exon's first base (10) with phase 2; its 3' end (the stop) is unchanged.
        let g = genome_with(80, &[(11, b"TAG"), (21, b"TAA")]);
        let path = [
            (10, 30, ExonType::Internal, 2),
            (50, 60, ExonType::Terminal, 1),
        ];
        let gene = resolve_partial(&g, 80, &path, Strand::Plus, None, false, false);
        assert_eq!(
            gene.cds,
            vec![Coordset::new(10, 23)],
            "CDS from exon.lend, ORF read from exon.lend + 2 (stop at 21..23)"
        );
        assert!(
            gene.five_utr.is_empty(),
            "the partial codon is coding, not UTR"
        );
        assert_eq!(gene.cds_start_phase, 2);
        assert!(gene.partial5, "path starts at an internal exon");
        assert!(!gene.partial3, "terminal exon + in-frame stop");
    }

    #[test]
    fn five_prime_partial_frame3_starts_cds_at_exon_with_phase_1() {
        // start frame 3: the next codon starts 1 base in, at 11: CTA GCC CCC TAA (stop
        // 20..22). Reading 2 bases in (from 12) would hit TAG at 12..14 at once.
        let g = genome_with(80, &[(12, b"TAG"), (20, b"TAA")]);
        let path = [
            (10, 30, ExonType::Internal, 3),
            (50, 60, ExonType::Terminal, 1),
        ];
        let gene = resolve_partial(&g, 80, &path, Strand::Plus, None, false, false);
        assert_eq!(
            gene.cds,
            vec![Coordset::new(10, 22)],
            "CDS from exon.lend, ORF read from exon.lend + 1 (stop at 20..22)"
        );
        assert!(gene.five_utr.is_empty());
        assert_eq!(gene.cds_start_phase, 1);
        assert!(gene.partial5);
        assert!(!gene.partial3);
    }

    #[test]
    fn five_prime_partial_frame2_minus_strand() {
        // The frame-2 case in the RC-local pass (hi = 200): Internal 11..31 (sf 2) ->
        // Terminal 51..61, i.e. forward 170..190 (5') and 140..150. RC-local codons from
        // 13: AGC CCC CCC TAA (stop 22..24); from 12 it would hit TAG at 12..14 at once.
        // The forward genome is the reverse complement of that RC-local sequence.
        let mut rc = vec![b'C'; 200];
        rc[11..14].copy_from_slice(b"TAG"); // RC-local 12..14
        rc[21..24].copy_from_slice(b"TAA"); // RC-local 22..24
        let g = write_fasta(&reverse_complement(&rc));
        let path = [
            (11, 31, ExonType::Internal, 2),
            (51, 61, ExonType::Terminal, 1),
        ];
        let gene = resolve_partial(&g, 200, &path, Strand::Minus, Some(200), false, false);
        assert_eq!(
            gene.exons,
            vec![Coordset::new(140, 150), Coordset::new(170, 190)]
        );
        // ORF RC-local 13..24 -> forward 177..188 (first complete codon at forward
        // rend - 2 = 188); the CDS starts at the exon's 5' base (forward 190), phase 2
        assert_eq!(gene.cds, vec![Coordset::new(177, 190)]);
        assert!(gene.five_utr.is_empty());
        assert_eq!(gene.cds_start_phase, 2);
        assert!(gene.partial5);
        assert!(!gene.partial3);
    }

    #[test]
    fn mid_codon_initial_exon_is_a_five_prime_partial_start() {
        // A plus-strand prediction whose 5'-most CDS row has GFF3 phase 1: Initial 10..30
        // at start frame 3 -> Terminal 60..90. The first complete codon starts 1 base in,
        // at 11; the in-frame stop is TAA at 88..90 (decoys: reading from 10 would hit TAA
        // at 13..15, from 12 TAG at 18..20). Being typed initial does not make the partial
        // codon a start: the CDS starts at the first base with phase 1, no 5'UTR stub.
        let g = genome_with(200, &[(13, b"TAA"), (18, b"TAG"), (88, b"TAA")]);
        let path = [
            (10, 30, ExonType::Initial, 3),
            (60, 90, ExonType::Terminal, 3),
        ];
        let gene = resolve_partial(&g, 200, &path, Strand::Plus, None, false, false);
        assert_eq!(gene.cds, vec![Coordset::new(10, 30), Coordset::new(60, 90)]);
        assert!(
            gene.five_utr.is_empty(),
            "the partial codon is coding, not UTR"
        );
        assert_eq!(gene.cds_start_phase, 1);
        assert!(gene.partial5, "a mid-codon first base is not a start codon");
        assert!(!gene.partial3);
    }

    #[test]
    fn mid_codon_initial_exon_minus_strand() {
        // A minus-strand prediction whose 5'-most CDS row has GFF3 phase 2, in the RC-local
        // pass (hi = 200): Initial 21..41 (start frame 2) -> Terminal 71..90, i.e. forward
        // 160..180 (5') and 111..130. RC-local codons from 23 reach the stop TAA at 88..90
        // (forward 111..113); from 21 the ORF would hit TAA at 24..26, from 22 TAG at 28..30.
        let mut rc = vec![b'C'; 200];
        rc[23..26].copy_from_slice(b"TAA"); // RC-local 24..26
        rc[27..30].copy_from_slice(b"TAG"); // RC-local 28..30
        rc[87..90].copy_from_slice(b"TAA"); // RC-local 88..90
        let g = write_fasta(&reverse_complement(&rc));
        let path = [
            (21, 41, ExonType::Initial, 2),
            (71, 90, ExonType::Terminal, 2),
        ];
        let gene = resolve_partial(&g, 200, &path, Strand::Minus, Some(200), false, false);
        assert_eq!(
            gene.exons,
            vec![Coordset::new(111, 130), Coordset::new(160, 180)]
        );
        // the CDS runs from the exon's 5' base (forward 180, phase 2) to the stop
        assert_eq!(
            gene.cds,
            vec![Coordset::new(111, 130), Coordset::new(160, 180)]
        );
        assert!(gene.five_utr.is_empty());
        assert_eq!(gene.cds_start_phase, 2);
        assert!(gene.partial5);
        assert!(!gene.partial3);
    }

    #[test]
    fn missing_sequence_never_claims_a_stop() {
        // a single-exon prediction on a contig the FASTA does not hold: no stop scan is
        // possible, so the gene keeps its raw exons but is never marked 3'-complete
        let g = genome_with(200, &[]);
        let mut c = chain(EvClass::AbinitioPrediction, 1.0, Strand::Plus, &[(50, 120)]);
        c.contig = "chrX".into();
        let region = ConsensusRegion {
            contig: "chrX".into(),
            span: Coordset::new(1, 200),
            chain_indices: vec![0],
        };
        let (cand, filt, w, eng) = config();
        let genes = consensus_region(&region, &[c], &g, &cand, &filt, &w, &eng, &[]);
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].cds, vec![Coordset::new(50, 120)]);
        assert!(!genes[0].partial5);
        assert!(genes[0].partial3, "no sequence, no stop scan: 3'-partial");
    }

    #[test]
    fn one_or_two_bases_past_the_stop_are_not_a_truncation() {
        // Initial 10..30 -> Terminal 50..rend, frame 1, TAA at 68..70. With 1-2 bases
        // after the stop the path holds no further complete codon: the stop is the
        // gene's own, so the gene is neither truncated nor 3'-partial (the stray bases
        // stay 3'UTR). Three bases after the stop are one more codon: truncated.
        let g = genome_with(200, &[(10, b"ATG"), (68, b"TAA")]);
        for (rend, truncated) in [(71, false), (72, false), (73, true)] {
            let path = [
                (10, 30, ExonType::Initial, 1),
                (50, rend, ExonType::Terminal, 1),
            ];
            let gene = resolve_partial(&g, 200, &path, Strand::Plus, None, false, false);
            assert_eq!(gene.cds, vec![Coordset::new(10, 30), Coordset::new(50, 70)]);
            assert_eq!(gene.three_utr, vec![Coordset::new(71, rend)]);
            assert_eq!(gene.truncated, truncated, "terminal exon 50..{rend}");
            assert!(!gene.partial3);
        }
    }

    #[test]
    fn premature_stop_marks_the_gene_truncated() {
        // Initial 10..30 -> Internal 50..70 -> Terminal 100..120, all frame 1. An in-frame
        // TAA at 56..58 ends the ORF inside the internal exon: the CDS stops there and the
        // rest of the path becomes 3'UTR. The path claims a terminus and a stop was found,
        // so partial3 stays false; truncated reports the cut-short CDS.
        let path = [
            (10, 30, ExonType::Initial, 1),
            (50, 70, ExonType::Internal, 1),
            (100, 120, ExonType::Terminal, 1),
        ];
        let g = genome_with(200, &[(10, b"ATG"), (56, b"TAA")]);
        let gene = resolve_partial(&g, 200, &path, Strand::Plus, None, false, false);
        assert_eq!(gene.cds, vec![Coordset::new(10, 30), Coordset::new(50, 58)]);
        assert_eq!(
            gene.three_utr,
            vec![Coordset::new(59, 70), Coordset::new(100, 120)]
        );
        assert!(gene.truncated);
        assert!(!gene.partial3);

        // the same path with its stop at the end of the terminal exon is not truncated
        let g = genome_with(200, &[(10, b"ATG"), (118, b"TAA")]);
        let gene = resolve_partial(&g, 200, &path, Strand::Plus, None, false, false);
        assert!(gene.three_utr.is_empty());
        assert!(!gene.truncated);
        assert!(!gene.partial3);
    }

    #[test]
    fn extend_orf_stop_completes_a_runoff_orf() {
        // Initial 10..30 -> Terminal 50..69: the ORF runs off the end of the terminal exon.
        // The first in-frame stop downstream of the path is TAA at 71..73.
        let g = genome_with(300, &[(71, b"TAA")]);
        let path = [
            (10, 30, ExonType::Initial, 1),
            (50, 69, ExonType::Terminal, 1),
        ];
        let gene = resolve_partial(&g, 300, &path, Strand::Plus, None, false, true);
        assert!(!gene.partial3, "extended to the downstream stop");
        assert!(!gene.partial5);
        // CDS and exon structure end at the rescued stop codon
        assert_eq!(gene.cds, vec![Coordset::new(10, 30), Coordset::new(50, 73)]);
        assert!(gene.three_utr.is_empty());
        assert_eq!(
            gene.exons,
            vec![Coordset::new(10, 30), Coordset::new(50, 73)]
        );
        assert_eq!(gene.support.coding_length, 45, "the extended CDS length");
    }

    #[test]
    fn extend_orf_stop_off_leaves_runoff_partial() {
        let g = genome_with(300, &[(71, b"TAA")]);
        let path = [
            (10, 30, ExonType::Initial, 1),
            (50, 69, ExonType::Terminal, 1),
        ];
        let gene = resolve_partial(&g, 300, &path, Strand::Plus, None, false, false);
        assert!(gene.partial3, "no extension -> still 3'-partial");
        assert_eq!(
            gene.exons,
            vec![Coordset::new(10, 30), Coordset::new(50, 69)]
        );
        assert_eq!(gene.support.coding_length, 0, "untouched: the input value");
    }

    #[test]
    fn extend_orf_stop_minus_strand() {
        // the run-off case in the RC-local pass (hi = 300): Initial 10..30 -> Terminal
        // 50..69 with the stop TAA at RC-local 71..73, i.e. forward exons 271..291 (5')
        // and 232..251 and the stop at forward 228..230. On '-' the 3'-most exon is the
        // forward-first one, extended toward smaller coordinates.
        let mut rc = vec![b'C'; 300];
        rc[70..73].copy_from_slice(b"TAA"); // RC-local 71..73
        let g = write_fasta(&reverse_complement(&rc));
        let path = [
            (10, 30, ExonType::Initial, 1),
            (50, 69, ExonType::Terminal, 1),
        ];
        let gene = resolve_partial(&g, 300, &path, Strand::Minus, Some(300), false, true);
        assert!(!gene.partial3);
        assert_eq!(
            gene.exons,
            vec![Coordset::new(228, 251), Coordset::new(271, 291)]
        );
        assert_eq!(
            gene.cds,
            vec![Coordset::new(228, 251), Coordset::new(271, 291)]
        );
        assert!(gene.three_utr.is_empty());
        assert_eq!(gene.support.coding_length, 45);
    }

    #[test]
    fn extend_orf_stop_leaves_an_internal_ending_path_partial() {
        // Initial 10..30 -> Internal 50..69: the path claims no 3' terminus, so it is not
        // read past the internal exon's evidence donor into the intron, stop or no stop
        let g = genome_with(300, &[(71, b"TAA")]);
        let path = [
            (10, 30, ExonType::Initial, 1),
            (50, 69, ExonType::Internal, 1),
        ];
        let gene = resolve_partial(&g, 300, &path, Strand::Plus, None, false, true);
        assert!(gene.partial3);
        assert_eq!(
            gene.exons,
            vec![Coordset::new(10, 30), Coordset::new(50, 69)]
        );
        assert_eq!(gene.support.coding_length, 0);
    }

    #[test]
    fn require_atg_extends_five_prime_to_an_upstream_atg() {
        // Initial 100..120 -> Terminal 200..229 (frame 1), stop TAA at 227..229. The path
        // itself holds no ATG (all-C), but an in-frame ATG sits upstream at 70..72.
        let g = genome_with(300, &[(70, b"ATG"), (227, b"TAA")]);
        let path = [
            (100, 120, ExonType::Initial, 1),
            (200, 229, ExonType::Terminal, 1),
        ];
        let gene = resolve_partial(&g, 300, &path, Strand::Plus, None, true, false);
        assert!(!gene.partial5, "rescued by the upstream ATG");
        assert!(!gene.partial3);
        assert_eq!(
            gene.cds,
            vec![Coordset::new(70, 120), Coordset::new(200, 229)]
        );
        // no bounding stop upstream -> the 5'UTR runs to the contig start
        assert_eq!(gene.five_utr, vec![Coordset::new(1, 69)]);
        assert!(gene.three_utr.is_empty());
        assert_eq!(gene.support.coding_length, 81, "the repaired CDS length");
    }

    #[test]
    fn five_prime_extension_trims_utr_at_the_bounding_stop() {
        // as above, plus an in-frame stop TAA at 61..63 bounding the upstream ORF:
        // the 5'UTR must be trimmed to start just past it (genomic 64)
        let g = genome_with(300, &[(61, b"TAA"), (70, b"ATG"), (227, b"TAA")]);
        let path = [
            (100, 120, ExonType::Initial, 1),
            (200, 229, ExonType::Terminal, 1),
        ];
        let gene = resolve_partial(&g, 300, &path, Strand::Plus, None, true, false);
        assert!(!gene.partial5);
        assert_eq!(
            gene.cds,
            vec![Coordset::new(70, 120), Coordset::new(200, 229)]
        );
        assert_eq!(
            gene.five_utr,
            vec![Coordset::new(64, 69)],
            "UTR trimmed to just past the bounding stop"
        );
    }

    #[test]
    fn require_atg_trims_to_the_first_in_frame_atg() {
        // Initial 10..30 -> Terminal 50..61 on an all-C genome: no ATG at the exon start,
        // the first in-frame ATG is at 16..18; the terminal stop TAA is at 59..61.
        let g = genome_with(80, &[(16, b"ATG"), (59, b"TAA")]);
        let path = [
            (10, 30, ExonType::Initial, 1),
            (50, 61, ExonType::Terminal, 1),
        ];
        let gene = resolve_partial(&g, 80, &path, Strand::Plus, None, true, false);
        assert_eq!(
            gene.cds,
            vec![Coordset::new(16, 30), Coordset::new(50, 61)],
            "CDS starts at the ATG"
        );
        assert_eq!(gene.five_utr, vec![Coordset::new(10, 15)]);
        assert!(!gene.partial5, "rescued by the in-frame ATG");
        assert!(!gene.partial3, "stop at the very end of the path");
        assert!(gene.three_utr.is_empty());
        assert_eq!(gene.support.coding_length, 27);
    }

    #[test]
    fn require_atg_flags_partial5_when_no_in_frame_atg() {
        // same path but no ATG anywhere: stays 5'-partial, CDS unchanged
        let g = genome_with(80, &[(59, b"TAA")]);
        let path = [
            (10, 30, ExonType::Initial, 1),
            (50, 61, ExonType::Terminal, 1),
        ];
        let gene = resolve_partial(&g, 80, &path, Strand::Plus, None, true, false);
        assert!(gene.partial5, "no in-frame ATG -> 5'-partial");
        assert_eq!(gene.cds, vec![Coordset::new(10, 30), Coordset::new(50, 61)]);
        assert_eq!(gene.support.coding_length, 0, "unmodified: the input value");
    }

    #[test]
    fn require_atg_off_keeps_unchecked_start() {
        // flag off: no ATG check at all (combinr's non-canonical default)
        let g = genome_with(80, &[(16, b"ATG"), (59, b"TAA")]);
        let path = [
            (10, 30, ExonType::Initial, 1),
            (50, 61, ExonType::Terminal, 1),
        ];
        let gene = resolve_partial(&g, 80, &path, Strand::Plus, None, false, false);
        assert_eq!(gene.cds, vec![Coordset::new(10, 30), Coordset::new(50, 61)]);
        assert!(!gene.partial5);
    }

    #[test]
    fn require_atg_minus_strand_trims_to_the_first_in_frame_atg() {
        // the trim case in the RC-local pass (hi = 80): Initial 10..30 -> Terminal 50..61,
        // ATG at RC-local 16..18, stop at 59..61; forward exons 51..71 (5') and 20..31
        let mut rc = vec![b'C'; 80];
        rc[15..18].copy_from_slice(b"ATG"); // RC-local 16..18
        rc[58..61].copy_from_slice(b"TAA"); // RC-local 59..61
        let g = write_fasta(&reverse_complement(&rc));
        let path = [
            (10, 30, ExonType::Initial, 1),
            (50, 61, ExonType::Terminal, 1),
        ];
        let gene = resolve_partial(&g, 80, &path, Strand::Minus, Some(80), true, false);
        assert_eq!(
            gene.exons,
            vec![Coordset::new(20, 31), Coordset::new(51, 71)]
        );
        assert_eq!(gene.cds, vec![Coordset::new(20, 31), Coordset::new(51, 65)]);
        assert_eq!(gene.five_utr, vec![Coordset::new(66, 71)]);
        assert!(!gene.partial5);
        assert!(!gene.partial3);
    }

    #[test]
    fn require_atg_minus_strand_extends_to_a_stop_bounded_atg() {
        // the extension case in the RC-local pass (hi = 300): Initial 100..120 -> Terminal
        // 200..229 with no ATG in the path; upstream an in-frame ATG at RC-local 70..72,
        // bounded by the stop TAA at 61..63. Forward: exons 181..201 (5') and 72..101; the
        // 5'-most exon is the forward-last one, extended toward larger coordinates.
        let mut rc = vec![b'C'; 300];
        rc[60..63].copy_from_slice(b"TAA"); // RC-local 61..63
        rc[69..72].copy_from_slice(b"ATG"); // RC-local 70..72
        rc[226..229].copy_from_slice(b"TAA"); // RC-local 227..229
        let g = write_fasta(&reverse_complement(&rc));
        let path = [
            (100, 120, ExonType::Initial, 1),
            (200, 229, ExonType::Terminal, 1),
        ];
        let gene = resolve_partial(&g, 300, &path, Strand::Minus, Some(300), true, false);
        assert!(!gene.partial5);
        assert_eq!(
            gene.exons,
            vec![Coordset::new(72, 101), Coordset::new(181, 237)]
        );
        assert_eq!(
            gene.cds,
            vec![Coordset::new(72, 101), Coordset::new(181, 231)]
        );
        // the 5'UTR starts just past the bounding stop (RC-local 64 = forward 237)
        assert_eq!(gene.five_utr, vec![Coordset::new(232, 237)]);
        assert_eq!(gene.support.coding_length, 81);
    }

    #[test]
    fn require_atg_extends_but_never_trims_a_mid_codon_start() {
        // Initial 10..30 at start frame 3 (a 5'-partial prediction, first codon at 11) ->
        // Terminal 50..62, stop TAA at 60..62, and an in-frame ATG inside the path at
        // 17..19. That internal Met is not taken as the start: with no ATG upstream the
        // gene stays 5'-partial with its CDS from the first base, phase 1.
        let path = [
            (10, 30, ExonType::Initial, 3),
            (50, 62, ExonType::Terminal, 3),
        ];
        let g = genome_with(200, &[(17, b"ATG"), (60, b"TAA")]);
        let gene = resolve_partial(&g, 200, &path, Strand::Plus, None, true, false);
        assert!(gene.partial5);
        assert_eq!(gene.cds, vec![Coordset::new(10, 30), Coordset::new(50, 62)]);
        assert_eq!(gene.cds_start_phase, 1);

        // an in-frame ATG upstream (2..4) is reached by extending the first exon
        let g = genome_with(200, &[(2, b"ATG"), (17, b"ATG"), (60, b"TAA")]);
        let gene = resolve_partial(&g, 200, &path, Strand::Plus, None, true, false);
        assert!(!gene.partial5);
        assert_eq!(gene.cds, vec![Coordset::new(2, 30), Coordset::new(50, 62)]);
        assert_eq!(gene.five_utr, vec![Coordset::new(1, 1)]);
        assert_eq!(gene.cds_start_phase, 0);
    }

    #[test]
    fn require_atg_leaves_an_internal_first_exon_unextended() {
        // the extension case above, but the first exon is Internal 10..30 (5'-partial
        // from an evidence acceptor, first codon at 11) -> Terminal 50..62, stop TAA at
        // 60..62, with the in-frame ATG at 2..4 upstream. The path claims no start, so it
        // is neither trimmed nor extended: 5'-partial, CDS from its first base, phase 1.
        let path = [
            (10, 30, ExonType::Internal, 3),
            (50, 62, ExonType::Terminal, 3),
        ];
        let g = genome_with(200, &[(2, b"ATG"), (60, b"TAA")]);
        let gene = resolve_partial(&g, 200, &path, Strand::Plus, None, true, false);
        assert!(gene.partial5);
        assert!(!gene.partial3);
        assert_eq!(
            gene.exons,
            vec![Coordset::new(10, 30), Coordset::new(50, 62)]
        );
        assert_eq!(gene.cds, vec![Coordset::new(10, 30), Coordset::new(50, 62)]);
        assert!(gene.five_utr.is_empty());
        assert_eq!(gene.cds_start_phase, 1);
        assert_eq!(gene.support.coding_length, 0, "unmodified: the input value");
    }

    #[test]
    fn require_atg_keeps_a_single_exon_gene_without_atg() {
        // Single 50..121 on an all-C genome, stop TAA at 119..121, no ATG anywhere. The
        // trellis is unchanged by the flag, so the gene is still called; it is flagged
        // 5'-partial with its CDS untouched.
        let g = genome_with(200, &[(119, b"TAA")]);
        let chains = vec![chain(
            EvClass::AbinitioPrediction,
            1.0,
            Strand::Plus,
            &[(50, 121)],
        )];
        let region = ConsensusRegion {
            contig: "chr1".into(),
            span: Coordset::new(1, 200),
            chain_indices: vec![0],
        };
        let (cand, filt, w, mut eng) = config();
        eng.require_atg_start = true;
        let genes = consensus_region(&region, &chains, &g, &cand, &filt, &w, &eng, &[]);
        let plus: Vec<_> = genes.iter().filter(|x| x.orient == Strand::Plus).collect();
        assert_eq!(plus.len(), 1, "the no-ATG single-exon gene survives");
        assert_eq!(plus[0].cds, vec![Coordset::new(50, 121)]);
        assert!(plus[0].partial5);
        assert!(!plus[0].partial3);
    }

    #[test]
    fn require_atg_keeps_both_genes_of_a_region() {
        // Two plus-strand genes in one region, bridged by a protein match: A starts with
        // ATG, B does not. The flag only repairs starts after the trellis, so both genes
        // are called with it on, as with it off; B is flagged 5'-partial.
        let g = genome_with(1000, &[(101, b"ATG"), (228, b"TAA"), (628, b"TAA")]);
        let mut a = chain(
            EvClass::AbinitioPrediction,
            1.0,
            Strand::Plus,
            &[(101, 130), (201, 230)],
        );
        a.accession = "A".into();
        let mut b = chain(
            EvClass::AbinitioPrediction,
            1.0,
            Strand::Plus,
            &[(501, 530), (601, 630)],
        );
        b.accession = "B".into();
        let mut p = chain(
            EvClass::Protein,
            1.0,
            Strand::Plus,
            &[(211, 220), (511, 520)],
        );
        p.accession = "P".into();
        p.ev_type = "nap".into();
        let chains = vec![a, b, p];
        let region = ConsensusRegion {
            contig: "chr1".into(),
            span: Coordset::new(1, 1000),
            chain_indices: vec![0, 1, 2],
        };
        let (cand, filt, w, mut eng) = config();
        for require in [false, true] {
            eng.require_atg_start = require;
            let genes = consensus_region(&region, &chains, &g, &cand, &filt, &w, &eng, &[]);
            let gene_at = |lend: i64| {
                genes
                    .iter()
                    .find(|x| x.orient == Strand::Plus && x.span().0 == lend)
            };
            let (ga, gb) = (gene_at(101), gene_at(501));
            assert!(ga.is_some() && gb.is_some(), "require_atg_start={require}");
            assert!(!ga.unwrap().partial5);
            assert_eq!(gb.unwrap().partial5, require);
        }
    }

    #[test]
    fn require_atg_with_extend_orf_stop_keeps_the_start_exon_whole() {
        // Initial 10..30 -> Terminal 50..69: the ATG trim moves the start to 16 (5'UTR
        // 10..15) and the run-off ORF is extended to the stop TAA at 71..73. The start
        // exon stays one exon, 10..30, rather than splitting into 5'UTR and CDS pieces.
        let g = genome_with(300, &[(16, b"ATG"), (71, b"TAA")]);
        let path = [
            (10, 30, ExonType::Initial, 1),
            (50, 69, ExonType::Terminal, 1),
        ];
        let gene = resolve_partial(&g, 300, &path, Strand::Plus, None, true, true);
        assert_eq!(
            gene.exons,
            vec![Coordset::new(10, 30), Coordset::new(50, 73)]
        );
        assert_eq!(gene.cds, vec![Coordset::new(16, 30), Coordset::new(50, 73)]);
        assert_eq!(gene.five_utr, vec![Coordset::new(10, 15)]);
        assert!(gene.three_utr.is_empty());
        assert!(!gene.partial5 && !gene.partial3);
        assert_eq!(gene.support.coding_length, 39);
    }
}
