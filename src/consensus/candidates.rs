//! Turn weighted evidence chains into typed, framed candidate exons for a region,
//! populating the per-base coding vector, the intron-score map, and the site sets.
//!
//! Port of EVM's `load_prediction_data` / `load_evidence_data` /
//! `instantiate_evidence_based_exons` / `add_introns`, adapted to the non-canonical,
//! evidence-derived model:
//! - gene predictions (ABINITIO/OTHER): each CDS segment is a typed exon
//!   (initial/internal/terminal/single by position), its frame computed from the
//!   cumulative coding length (`cds_length % 3 + 1`) seeded by the GFF phase column
//!   of the chain's 5'-most CDS row (a nonzero phase marks a 5'-partial prediction);
//! - proteins: paint the coding vector for every segment, and seed internal exons for
//!   internal segments with all stop-free reading frames;
//! - transcripts: seed internal exons (no coding-vector painting — their contribution
//!   is exon-specific, applied in `score_exons`, M2);
//! - every chain's gaps become evidence-supported introns (no GT-AG gate).
//!
//! Frames are validated only against in-frame stops under the configurable
//! [`GeneticCode`] (see [`super::exon::determine_good_phases`]) — never canonical sites.
//!
//! Minus-strand handling (reverse-complement + transpose) lands in M3; this builds the
//! forward-strand candidates only.

use crate::consensus::evidence::EvidenceChain;
use crate::consensus::exon::{
    ExonCandidate, ExonType, determine_good_phases, end_frame_for, frame_offset,
};
use crate::consensus::peaks::PeakSignal;
use crate::consensus::region::ConsensusRegion;
use crate::consensus::vectors::{IntronScores, RegionVectors};
use crate::consensus::weights::EvClass;
use crate::io::fasta::Fasta;
use crate::model::{Coordset, Strand};
use crate::orf::GeneticCode;
use std::collections::HashMap;

/// Max bp to extend a terminal alignment segment downstream looking for a stop codon.
pub(crate) const MAX_TERMINAL_EXTEND: i64 = 10_000;

/// How the unmasked length of one intron maps to its score contribution
/// (`--intron-score-cap`), per intron per supporting chain.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IntronScoreScale {
    /// `weight x len` (EVM behaviour; the default).
    Linear,
    /// `weight x min(len, cap)`.
    Cap(i64),
    /// `weight x len` up to 1 kb, then logarithmic
    /// (`weight x 1 kb x (1 + ln(len / 1 kb))`): giant introns keep scoring more but
    /// can no longer outscore gene-splitting signals by orders of magnitude.
    Smooth,
}

/// Knee (bp) for [`IntronScoreScale::Smooth`]: introns up to this length score linearly.
pub(crate) const INTRON_SMOOTH_KNEE: i64 = 1_000;

impl IntronScoreScale {
    /// The length contributing to an intron's score under this scaling.
    pub fn scored_len(self, len: i64) -> i64 {
        match self {
            Self::Linear => len,
            Self::Cap(cap) => len.min(cap),
            Self::Smooth => {
                if len <= INTRON_SMOOTH_KNEE {
                    len
                } else {
                    let k = INTRON_SMOOTH_KNEE as f64;
                    (k * (1.0 + (len as f64 / k).ln())).round() as i64
                }
            }
        }
    }
}

impl std::str::FromStr for IntronScoreScale {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "0" | "off" => Ok(Self::Linear),
            "auto" => Ok(Self::Smooth),
            _ => match s.parse::<i64>() {
                Ok(n) if n > 0 => Ok(Self::Cap(n)),
                _ => Err(format!(
                    "invalid --intron-score-cap {s:?} (expected 0, auto, or a bp cap)"
                )),
            },
        }
    }
}

/// Tunables for candidate generation.
#[derive(Clone, Copy, Debug)]
pub struct CandidateParams {
    pub min_intron_length: i64,
    pub code: GeneticCode,
    /// Extend protein/transcript terminal segments downstream to a genetic-code stop,
    /// forming terminal exons (EVM `--extend_to_terminal`, stop half only — non-canonical).
    pub extend_terminal_stop: bool,
    /// Boost intergenic scores near start/stop evidence-density peaks
    /// (EVM `augment_intergenic_from_start_stop_peaks`).
    pub peak_augment: bool,
    /// Cap on the unmasked length contributing to an intron's score. Anti-fusion knob
    /// for repeat-rich genomes: without it, a giant intron from a misassembled
    /// transcript outweighs any intergenic signal and the trellis absorbs the
    /// neighbouring gene. `Cap(n)` hard-caps, `Smooth` scales logarithmically.
    pub intron_score_scale: IntronScoreScale,
    /// Sum of all gene-prediction weights (EVM `SUM_GENEPRED_WEIGHTS`), the peak threshold
    /// and the augmented intergenic level.
    pub sum_genepred_weights: f64,
}

/// Everything the trellis needs for one region.
pub struct RegionData {
    pub contig: String,
    pub span: Coordset,
    pub exons: Vec<ExonCandidate>,
    pub vectors: RegionVectors,
    pub introns: IntronScores,
}

/// Build candidate exons + scoring data for a region in its working coordinate space
/// (the caller reverse-complements the region for the minus strand). `mask` lists
/// repeat-masked intervals (in the same coordinate space) excluded from all scoring.
pub fn build_candidates(
    region: &ConsensusRegion,
    chains: &[EvidenceChain],
    genome: &Fasta,
    params: &CandidateParams,
    mask: &[Coordset],
) -> RegionData {
    let span = region.span;
    let len = (span.rend - span.lend + 1).max(0) as usize;
    let mut b = Builder {
        contig: region.contig.clone(),
        genome,
        code: &params.code,
        min_intron_length: params.min_intron_length,
        extend_terminal_stop: params.extend_terminal_stop,
        intron_score_scale: params.intron_score_scale,
        sum_genepred_weights: params.sum_genepred_weights,
        vectors: RegionVectors::new(span.lend, len),
        introns: IntronScores::default(),
        exons: Vec::new(),
        dedup: HashMap::new(),
        abinitio_spans: HashMap::new(),
        // Only allocate the per-base peak signal when the opt-in feature is on.
        peaks: params.peak_augment.then(|| PeakSignal::new(span.lend, len)),
    };

    // Mask repeats before painting any coverage, so masked bases contribute nothing.
    for iv in mask {
        b.vectors.set_masked(iv.lend, iv.rend);
    }

    for &ci in &region.chain_indices {
        let chain = &chains[ci];
        if chain.orient != Strand::Plus {
            continue; // this builder works in one (forward-oriented) coordinate space;
            // the minus strand is handled by reverse-complementing the region upstream
        }
        match chain.ev_class {
            EvClass::AbinitioPrediction | EvClass::OtherPrediction => b.add_prediction(chain),
            EvClass::Protein => b.add_protein(chain),
            EvClass::Transcript => b.add_transcript(chain),
        }
    }
    b.populate_intergenic();
    if let Some(peaks) = b.peaks.take() {
        peaks.augment(&b.exons, &mut b.vectors, b.sum_genepred_weights);
    }

    RegionData {
        contig: region.contig.clone(),
        span,
        exons: b.exons,
        vectors: b.vectors,
        introns: b.introns,
    }
}

struct Builder<'a> {
    contig: String,
    genome: &'a Fasta,
    code: &'a GeneticCode,
    min_intron_length: i64,
    extend_terminal_stop: bool,
    intron_score_scale: IntronScoreScale,
    sum_genepred_weights: f64,
    vectors: RegionVectors,
    introns: IntronScores,
    exons: Vec<ExonCandidate>,
    dedup: HashMap<(i64, i64, ExonType, u8), usize>,
    /// Abinitio prediction gene spans by ev_type, for intergenic-zone scoring.
    abinitio_spans: HashMap<String, (f64, Vec<Coordset>)>,
    /// Start/stop termini signal for the opt-in `--peak-augment` feature; `None` when off.
    peaks: Option<PeakSignal>,
}

impl Builder<'_> {
    fn add_prediction(&mut self, chain: &EvidenceChain) {
        let n = chain.links.len();
        if n == 0 {
            return;
        }
        let paints = chain.ev_class.paints_coding_vector(); // abinitio yes, other no
        let is_abinitio = chain.ev_class == EvClass::AbinitioPrediction;
        if is_abinitio {
            self.abinitio_spans
                .entry(chain.ev_type.clone())
                .or_insert_with(|| (chain.weight, Vec::new()))
                .1
                .push(chain.span);
        }
        // GFF3 phase of the 5'-most CDS row = bases skipped to the first complete codon.
        // Seed the cumulative length so the inherited frames match the predictor's own
        // reading frame: phase 0 -> frame 1 (cum_len 0), phase 1 -> frame 3 (cum_len 2),
        // phase 2 -> frame 2 (cum_len 1).
        let mut cum_len: i64 = match chain.cds_start_phase {
            1 => 2,
            2 => 1,
            _ => 0,
        };
        for (i, seg) in chain.links.iter().enumerate() {
            let exon_type = position_type(i, n);
            let start_frame = (cum_len.rem_euclid(3) + 1) as u8;
            if paints {
                self.vectors.add_coverage(seg.lend, seg.rend, chain.weight);
            }
            self.try_add_exon(*seg, exon_type, start_frame, chain);
            cum_len += seg.len();
        }
        self.add_introns(chain);
    }

    fn add_protein(&mut self, chain: &EvidenceChain) {
        let n = chain.links.len();
        // protein chain termini feed the start/stop peak signal (only built when peak
        // augmentation is enabled).
        if let Some(peaks) = &mut self.peaks {
            peaks.bump_begin(chain.span.lend, chain.weight);
            peaks.bump_end(chain.span.rend, chain.weight);
        }
        for (i, seg) in chain.links.iter().enumerate() {
            // proteins paint the coding vector for every segment, splice or not
            self.vectors.add_coverage(seg.lend, seg.rend, chain.weight);
            if n > 1 && i != 0 && i != n - 1 {
                self.add_internal_evidence_exon(*seg, chain);
            }
        }
        self.add_introns(chain);
        if self.extend_terminal_stop
            && let Some(last) = chain.links.last()
        {
            self.try_extend_terminal(*last, chain);
        }
    }

    fn add_transcript(&mut self, chain: &EvidenceChain) {
        let n = chain.links.len();
        for (i, seg) in chain.links.iter().enumerate() {
            if n > 1 && i != 0 && i != n - 1 {
                // transcripts do not paint the coding vector; exon-specific weight is
                // applied in score_exons (M2)
                self.add_internal_evidence_exon(*seg, chain);
            }
        }
        self.add_introns(chain);
        if self.extend_terminal_stop
            && let Some(last) = chain.links.last()
        {
            self.try_extend_terminal(*last, chain);
        }
    }

    /// Terminal-stop extension (EVM `extend_exon_downstream_to_stop`, non-canonical):
    /// extend a 3' alignment segment downstream, in each stop-free frame, to the first
    /// in-frame genetic-code stop, forming a terminal exon. Lets protein/transcript-only
    /// loci acquire a 3' terminus (anti-false-negative). EVM's 5'/start half is not ported
    /// at the candidate stage, where it would gate candidates on a genome ATG scan; the
    /// opt-in `--require-atg-start` repairs a called gene's start after the trellis
    /// instead (`engine::resolve`).
    fn try_extend_terminal(&mut self, seg: Coordset, chain: &EvidenceChain) {
        let phases = {
            let s = self
                .genome
                .subseq(&self.contig, seg.lend, seg.rend)
                .unwrap_or_default();
            determine_good_phases(s, self.code)
        };
        for phase in phases {
            if let Some(stop_end) = self.downstream_stop(seg, phase)
                && stop_end >= seg.rend
            {
                self.insert_exon(
                    Coordset::new(seg.lend, stop_end),
                    ExonType::Terminal,
                    phase,
                    chain,
                );
            }
        }
    }

    /// Genomic position of the last base of the first in-frame stop at or after `seg`
    /// (reading from `seg.lend` in the given phase), within the extension cap.
    fn downstream_stop(&self, seg: Coordset, phase: u8) -> Option<i64> {
        let offset = frame_offset(phase) as i64;
        let limit = seg.rend + MAX_TERMINAL_EXTEND;
        let mut cs = seg.lend + offset;
        while cs + 2 <= limit {
            let codon = self.genome.subseq(&self.contig, cs, cs + 2)?;
            if self.code.is_stop(codon) {
                return Some(cs + 2);
            }
            cs += 3;
        }
        None
    }

    /// Add internal exons in every reading frame with no in-frame stop.
    fn add_internal_evidence_exon(&mut self, seg: Coordset, chain: &EvidenceChain) {
        let phases = {
            let seq = self
                .genome
                .subseq(&self.contig, seg.lend, seg.rend)
                .unwrap_or_default();
            determine_good_phases(seq, self.code)
        };
        for phase in phases {
            self.insert_exon(seg, ExonType::Internal, phase, chain);
        }
    }

    /// Validate a single computed frame against in-frame stops, then add the exon.
    /// Terminal/single exons trim their trailing stop codon before the check.
    fn try_add_exon(
        &mut self,
        seg: Coordset,
        exon_type: ExonType,
        start_frame: u8,
        chain: &EvidenceChain,
    ) {
        let check_rend = match exon_type {
            ExonType::Terminal | ExonType::Single => seg.rend - 3,
            _ => seg.rend,
        };
        if check_rend >= seg.lend {
            let good = {
                let seq = self
                    .genome
                    .subseq(&self.contig, seg.lend, check_rend)
                    .unwrap_or_default();
                determine_good_phases(seq, self.code)
            };
            if !good.contains(&start_frame) {
                return; // computed frame has an in-frame stop -> skip this exon
            }
        }
        self.insert_exon(seg, exon_type, start_frame, chain);
    }

    fn insert_exon(
        &mut self,
        seg: Coordset,
        exon_type: ExonType,
        start_frame: u8,
        chain: &EvidenceChain,
    ) {
        let key = (seg.lend, seg.rend, exon_type, start_frame);
        let ev = (chain.accession.clone(), chain.ev_type.clone());
        if let Some(&idx) = self.dedup.get(&key) {
            self.exons[idx].evidence.push(ev);
            return;
        }
        let left = self.boundary(seg.lend, seg.lend + 1);
        let right = self.boundary(seg.rend - 1, seg.rend);
        let exon = ExonCandidate {
            coords: seg,
            orient: Strand::Plus,
            exon_type,
            start_frame,
            end_frame: end_frame_for(start_frame, seg.len()),
            left_seq_boundary: left,
            right_seq_boundary: right,
            evidence: vec![ev],
        };
        self.dedup.insert(key, self.exons.len());
        self.exons.push(exon);
    }

    /// Two genome bases `[lend, rend]`, padded with `N` if out of bounds.
    fn boundary(&self, lend: i64, rend: i64) -> [u8; 2] {
        let s = self
            .genome
            .subseq(&self.contig, lend, rend)
            .unwrap_or_default();
        let mut b = [b'N', b'N'];
        for (k, &c) in s.iter().take(2).enumerate() {
            b[k] = c;
        }
        b
    }

    /// `populate_intergenic_regions` (~3173): for each ab-initio prediction program,
    /// add its weight to every base of the gaps between its predicted genes.
    fn populate_intergenic(&mut self) {
        // `abinitio_spans` is not used again — take ownership instead of cloning it.
        let progs = std::mem::take(&mut self.abinitio_spans);
        for (weight, mut spans) in progs.into_values() {
            spans.sort_by_key(|c| c.lend);
            for w in spans.windows(2) {
                let gap = w[0].gap_to(&w[1]);
                if gap.rend >= gap.lend {
                    self.vectors.add_intergenic(gap.lend, gap.rend, weight);
                }
            }
        }
    }

    /// Every gap between consecutive links becomes an evidence-supported intron.
    /// EVM registers an intron even when repeat-masking zeroes its score
    /// (`$INTRONS_TO_SCORE{$key} += 0`), so the trellis can still bridge introns that
    /// are fully covered by repeats (TE-intron genes stay connected); masked bases
    /// simply contribute no score. Do the same.
    fn add_introns(&mut self, chain: &EvidenceChain) {
        for w in chain.links.windows(2) {
            let intron = w[0].gap_to(&w[1]);
            if intron.len() < self.min_intron_length {
                continue;
            }
            let unmasked = self.vectors.unmasked_len(intron.lend, intron.rend);
            let scored = self.intron_score_scale.scored_len(unmasked);
            let key = (intron.lend, intron.rend);
            self.introns.add(key, chain.weight, scored);
            self.introns
                .attribute(key, &chain.accession, &chain.ev_type);
        }
    }
}

pub(crate) fn position_type(i: usize, n: usize) -> ExonType {
    if n == 1 {
        ExonType::Single
    } else if i == 0 {
        ExonType::Initial
    } else if i == n - 1 {
        ExonType::Terminal
    } else {
        ExonType::Internal
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consensus::weights::EvClass;
    use std::io::Write;

    /// A FASTA of `len` 'C' bases with byte edits applied (1-based positions).
    fn genome(len: usize, edits: &[(i64, &[u8])]) -> Fasta {
        let mut seq = vec![b'C'; len];
        for &(pos1, bytes) in edits {
            for (k, &b) in bytes.iter().enumerate() {
                seq[(pos1 - 1) as usize + k] = b;
            }
        }
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        writeln!(tmp, ">chr1").unwrap();
        tmp.write_all(&seq).unwrap();
        writeln!(tmp).unwrap();
        Fasta::load(tmp.path()).unwrap()
    }

    fn chain(class: EvClass, weight: f64, segs: &[(i64, i64)]) -> EvidenceChain {
        let links: Vec<Coordset> = segs.iter().map(|&(l, r)| Coordset::new(l, r)).collect();
        let span = Coordset::new(links.first().unwrap().lend, links.last().unwrap().rend);
        EvidenceChain {
            accession: "acc".into(),
            ev_type: "src".into(),
            ev_class: class,
            weight,
            contig: "chr1".into(),
            orient: Strand::Plus,
            span,
            links,
            cds_start_phase: 0,
        }
    }

    fn region(chains: &[EvidenceChain]) -> ConsensusRegion {
        ConsensusRegion {
            contig: "chr1".into(),
            span: Coordset::new(1, 200),
            chain_indices: (0..chains.len()).collect(),
        }
    }

    fn params() -> CandidateParams {
        CandidateParams {
            min_intron_length: 20,
            code: GeneticCode::default(),
            extend_terminal_stop: false,
            peak_augment: false,
            intron_score_scale: IntronScoreScale::Linear,
            sum_genepred_weights: 0.0,
        }
    }

    #[test]
    fn prediction_yields_typed_framed_exons_and_intron_coverage() {
        // two CDS exons, intron 19..49 (31 bp >= min 20), terminal stop TAA at 56..58
        let g = genome(70, &[(56, b"TAA")]);
        let chains = vec![chain(
            EvClass::AbinitioPrediction,
            1.0,
            &[(10, 18), (50, 58)],
        )];
        let rd = build_candidates(&region(&chains), &chains, &g, &params(), &[]);

        assert_eq!(rd.exons.len(), 2);
        assert_eq!(rd.exons[0].exon_type, ExonType::Initial);
        assert_eq!(rd.exons[0].start_frame, 1);
        assert_eq!(rd.exons[1].exon_type, ExonType::Terminal);
        // cumulative coding length before exon 2 = 9 -> 9 % 3 + 1 = 1
        assert_eq!(rd.exons[1].start_frame, 1);

        // abinitio paints the coding vector
        assert_eq!(rd.vectors.coding_sum(10, 18), 9.0);
        assert_eq!(rd.vectors.coding_sum(50, 58), 9.0);

        // one evidence-supported intron, scored weight*unmasked_len = 1*31
        assert!(rd.introns.contains((19, 49)));
        assert_eq!(rd.introns.score((19, 49)), Some(31.0));
        assert_eq!(
            rd.introns.evidence((19, 49)),
            &[("acc".to_string(), "src".to_string())]
        );
    }

    #[test]
    fn intron_evidence_lists_each_chain_once() {
        // three chains share intron 19..49: "acc"/"src" twice (deduped) and a second
        // accession; the score still accumulates all three.
        let g = genome(70, &[(56, b"TAA")]);
        let mut other = chain(EvClass::AbinitioPrediction, 1.0, &[(10, 18), (50, 58)]);
        other.accession = "acc2".into();
        let chains = vec![
            chain(EvClass::AbinitioPrediction, 1.0, &[(10, 18), (50, 58)]),
            chain(EvClass::Transcript, 1.0, &[(1, 18), (50, 60)]),
            other,
        ];
        let rd = build_candidates(&region(&chains), &chains, &g, &params(), &[]);
        assert_eq!(rd.introns.score((19, 49)), Some(93.0));
        assert_eq!(
            rd.introns.evidence((19, 49)),
            &[
                ("acc".to_string(), "src".to_string()),
                ("acc2".to_string(), "src".to_string())
            ]
        );
        assert!(rd.introns.evidence((1, 5)).is_empty());
    }

    #[test]
    fn protein_paints_all_segments_and_seeds_internal_exon_per_frame() {
        // three segments, introns 19..39 and 49..69 (each 21 bp >= 20), no stops
        let g = genome(90, &[]);
        let chains = vec![chain(
            EvClass::Protein,
            5.0,
            &[(10, 18), (40, 48), (70, 78)],
        )];
        let rd = build_candidates(&region(&chains), &chains, &g, &params(), &[]);

        // internal exons only for the middle segment, one per stop-free frame (3)
        assert_eq!(rd.exons.len(), 3);
        assert!(rd.exons.iter().all(|e| e.exon_type == ExonType::Internal));
        let frames: std::collections::BTreeSet<u8> =
            rd.exons.iter().map(|e| e.start_frame).collect();
        assert_eq!(frames, [1, 2, 3].into_iter().collect());

        // proteins paint every segment (9 bases * weight 5 = 45)
        assert_eq!(rd.vectors.coding_sum(10, 18), 45.0);
        assert_eq!(rd.vectors.coding_sum(40, 48), 45.0);
        assert_eq!(rd.vectors.coding_sum(70, 78), 45.0);

        assert_eq!(rd.introns.len(), 2);
    }

    #[test]
    fn transcript_seeds_internal_exons_without_painting_coding() {
        let g = genome(90, &[]);
        let chains = vec![chain(
            EvClass::Transcript,
            7.0,
            &[(10, 18), (40, 48), (70, 78)],
        )];
        let rd = build_candidates(&region(&chains), &chains, &g, &params(), &[]);

        assert_eq!(rd.exons.len(), 3); // middle segment, 3 frames
        // transcripts do NOT paint the coding vector
        assert_eq!(rd.vectors.coding_sum(1, 90), 0.0);
        // ...but the evidence is recorded on the exons for exon-specific scoring (M2)
        assert!(
            rd.exons
                .iter()
                .all(|e| e.evidence == vec![("acc".into(), "src".into())])
        );
    }

    #[test]
    fn prediction_exon_with_in_frame_stop_is_skipped() {
        // single-exon prediction 10..21; premature stop TAA at 13..15 (offset 3 -> phase 1)
        let g = genome(40, &[(13, b"TAA")]);
        let chains = vec![chain(EvClass::AbinitioPrediction, 1.0, &[(10, 21)])];
        let rd = build_candidates(&region(&chains), &chains, &g, &params(), &[]);
        // computed frame for a single exon is 1, which the premature stop invalidates
        assert!(rd.exons.is_empty());
    }

    #[test]
    fn duplicate_exon_merges_evidence() {
        // two predictions with an identical first CDS exon -> one exon, two evidences
        let g = genome(70, &[(56, b"TAA")]);
        let mut c1 = chain(EvClass::AbinitioPrediction, 1.0, &[(10, 18), (50, 58)]);
        c1.accession = "p1".into();
        c1.ev_type = "fgenesh".into();
        let mut c2 = chain(EvClass::AbinitioPrediction, 1.0, &[(10, 18), (50, 58)]);
        c2.accession = "p2".into();
        c2.ev_type = "genemark".into();
        let chains = vec![c1, c2];
        let rd = build_candidates(&region(&chains), &chains, &g, &params(), &[]);
        // identical exon structure -> 2 unique exons (initial + terminal), each with 2 evidences
        assert_eq!(rd.exons.len(), 2);
        assert!(rd.exons.iter().all(|e| e.evidence.len() == 2));
    }

    #[test]
    fn prediction_frames_honor_the_gff_phase_column() {
        // same two-exon prediction as `prediction_yields_typed_framed_exons_and_intron_coverage`,
        // but the chain carries a 5' phase of 1 (a 5'-partial input prediction)
        let g = genome(70, &[(56, b"TAA")]);
        let mut c = chain(EvClass::AbinitioPrediction, 1.0, &[(10, 18), (50, 58)]);
        c.cds_start_phase = 1;
        let chains = vec![c];
        let rd = build_candidates(&region(&chains), &chains, &g, &params(), &[]);
        assert_eq!(rd.exons[0].start_frame, 3, "phase 1 -> skip 1 base -> frame 3");
        // cum_len = 2 (seed) + 9 = 11 -> 11 % 3 + 1 = 3
        assert_eq!(rd.exons[1].start_frame, 3);
    }

    #[test]
    fn prediction_exon_with_stops_in_all_frames_is_skipped() {
        // three-exon prediction whose internal exon 50..61 has a stop in every frame:
        // TAA at 50..52 (phase 1), TGA at 54..56 (phase 3), TAG at 58..60 (phase 2).
        // Internal exons are checked untrimmed, so no frame survives; an empty `good`
        // list used to slip past the guard and the exon was inserted anyway.
        let g = genome(200, &[(50, b"TAA"), (54, b"TGA"), (58, b"TAG")]);
        let chains = vec![chain(
            EvClass::AbinitioPrediction,
            1.0,
            &[(10, 18), (50, 61), (100, 108)],
        )];
        let rd = build_candidates(&region(&chains), &chains, &g, &params(), &[]);
        assert!(!rd.exons.iter().any(|e| e.coords == Coordset::new(50, 61)));
        assert_eq!(rd.exons.len(), 2, "initial and terminal exons stay");
    }

    #[test]
    fn intron_score_scale_parses_and_scales() {
        use std::str::FromStr;
        assert_eq!(
            IntronScoreScale::from_str("0"),
            Ok(IntronScoreScale::Linear)
        );
        assert_eq!(
            IntronScoreScale::from_str("auto"),
            Ok(IntronScoreScale::Smooth)
        );
        assert_eq!(
            IntronScoreScale::from_str("2000"),
            Ok(IntronScoreScale::Cap(2000))
        );
        assert!(IntronScoreScale::from_str("-5").is_err());
        assert!(IntronScoreScale::from_str("banana").is_err());

        // linear and cap
        assert_eq!(IntronScoreScale::Linear.scored_len(15_000), 15_000);
        assert_eq!(IntronScoreScale::Cap(2000).scored_len(15_000), 2000);
        // smooth: linear up to the 1 kb knee, logarithmic beyond
        assert_eq!(IntronScoreScale::Smooth.scored_len(500), 500);
        assert_eq!(IntronScoreScale::Smooth.scored_len(1_000), 1_000);
        let s15k = IntronScoreScale::Smooth.scored_len(15_000);
        assert!((3000..4500).contains(&s15k), "15 kb -> {s15k}"); // ~3708
        let s50k = IntronScoreScale::Smooth.scored_len(50_000);
        assert!((4500..5500).contains(&s50k), "50 kb -> {s50k}"); // ~4912
        // monotone: longer introns still score more, but barely
        assert!(s50k > s15k && s50k < 2 * s15k);
    }

    #[test]
    fn fully_masked_intron_is_still_registered() {
        // EVM parity: masking zeroes an intron's score but must not sever the trellis
        // connection (TE-filled introns are common in repeat-rich genomes)
        let g = genome(70, &[(56, b"TAA")]);
        let chains = vec![chain(EvClass::AbinitioPrediction, 1.0, &[(10, 18), (50, 58)])];
        let mask = [Coordset::new(19, 49)];
        let rd = build_candidates(&region(&chains), &chains, &g, &params(), &mask);
        assert!(
            rd.introns.contains((19, 49)),
            "fully masked intron must still be registered"
        );
        assert_eq!(rd.introns.score((19, 49)), Some(0.0), "but scores zero");
        assert_eq!(
            rd.introns.evidence((19, 49)),
            &[("acc".to_string(), "src".to_string())]
        );
    }

    #[test]
    fn terminal_stop_extension_creates_terminal_exon() {
        // all-C genome with a TAA at 52..54 (in-frame-1 reading from 40)
        let g = genome(60, &[(52, b"TAA")]);
        let chains = vec![chain(EvClass::Protein, 5.0, &[(10, 18), (40, 48)])];

        // off: a 2-segment protein makes no terminal exon
        let mut p = params();
        let off = build_candidates(&region(&chains), &chains, &g, &p, &[]);
        assert!(!off.exons.iter().any(|e| e.exon_type == ExonType::Terminal));

        // on: the 3' segment extends downstream to the stop, forming a terminal exon
        p.extend_terminal_stop = true;
        let on = build_candidates(&region(&chains), &chains, &g, &p, &[]);
        assert!(on.exons.iter().any(|e| {
            e.exon_type == ExonType::Terminal && e.coords.lend == 40 && e.coords.rend == 54
        }));
    }
}
