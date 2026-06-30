//! Turn weighted evidence chains into typed, framed candidate exons for a region,
//! populating the per-base coding vector, the intron-score map, and the site sets.
//!
//! Port of EVM's `load_prediction_data` / `load_evidence_data` /
//! `instantiate_evidence_based_exons` / `add_introns`, adapted to the non-canonical,
//! evidence-derived model:
//! - gene predictions (ABINITIO/OTHER): each CDS segment is a typed exon
//!   (initial/internal/terminal/single by position), its frame computed from the
//!   cumulative coding length (`cds_length % 3 + 1`), not the GFF phase column;
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
use crate::consensus::exon::{ExonCandidate, ExonType, determine_good_phases, end_frame_for};
use crate::consensus::peaks::PeakSignal;
use crate::consensus::region::ConsensusRegion;
use crate::consensus::vectors::{IntronScores, RegionVectors};
use crate::consensus::weights::EvClass;
use crate::io::fasta::Fasta;
use crate::model::{Coordset, Strand};
use crate::orf::GeneticCode;
use std::collections::HashMap;

/// Max bp to extend a terminal alignment segment downstream looking for a stop codon.
const MAX_TERMINAL_EXTEND: i64 = 10_000;

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
        let mut cum_len = 0i64;
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
    /// loci acquire a 3' terminus (anti-false-negative); the 5'/start half is omitted
    /// because it would require canonical ATG scanning (see `avoid-canonical-splice-bias`).
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
        let offset = match phase {
            1 => 0,
            2 => 2,
            _ => 1,
        };
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
            if !good.is_empty() && !good.contains(&start_frame) {
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
    fn add_introns(&mut self, chain: &EvidenceChain) {
        for w in chain.links.windows(2) {
            let intron = w[0].gap_to(&w[1]);
            if intron.len() < self.min_intron_length {
                continue;
            }
            let unmasked = self.vectors.unmasked_len(intron.lend, intron.rend);
            if unmasked <= 0 {
                continue;
            }
            self.introns
                .add((intron.lend, intron.rend), chain.weight, unmasked);
        }
    }
}

fn position_type(i: usize, n: usize) -> ExonType {
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
