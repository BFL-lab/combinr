//! The cDNA alignment assembler — faithful port of
//! `pasa_cpp/cdna_alignment_assembler.cpp`.
//!
//! Given a set of alignments (already grouped into one cluster and sharing a
//! contig), produce the minimal set of maximal non-redundant assemblies such
//! that every input alignment belongs to at least one assembly.

pub mod compat;
pub mod dp;
pub mod merge;
pub mod wrapper;

use crate::error::{CombinrError, Result};
use crate::model::{Alignment, Strand};
use compat::{can_merge, encapsulates};
use dp::Lobject;
use merge::merge_alignments;

pub use wrapper::{ClusterAssembly, assemble_cluster};

/// Default fuzz distance, matching `CDNA_alignment_assembler::fuzzlength`.
pub const DEFAULT_FUZZLENGTH: i64 = 20;

pub struct Assembler {
    /// Input alignments, sorted by `coords.lend`.
    pub alignments: Vec<Alignment>,
    pub fuzzlength: i64,
    n: usize,
    /// `compat[i*n + j]` — symmetric pairwise compatibility.
    compat: Vec<bool>,
    /// `encap[i*n + j]` — directional: alignment `i` encapsulates `j`.
    encap: Vec<bool>,
    lobjects: Vec<Lobject>,
    /// Resulting assemblies (merged structures).
    pub assemblies: Vec<Alignment>,
    /// For each assembly, the indices of its member alignments (ascending).
    pub containment: Vec<Vec<usize>>,
}

impl Assembler {
    /// Build an assembler over `alignments`, sorting by `lend` like the C++ ctor.
    pub fn new(mut alignments: Vec<Alignment>) -> Self {
        // Stable sort by lend. The C++ std::sort is unstable; stable preserves
        // input (file) order on tied lends, which is the documented tie risk.
        alignments.sort_by_key(|a| a.coords.lend);
        let n = alignments.len();
        Assembler {
            alignments,
            fuzzlength: DEFAULT_FUZZLENGTH,
            n,
            compat: vec![false; n * n],
            encap: vec![false; n * n],
            lobjects: Vec::new(),
            assemblies: Vec::new(),
            containment: Vec::new(),
        }
    }

    pub fn set_fuzzlength(&mut self, f: i64) {
        self.fuzzlength = f;
    }

    #[inline]
    fn compat(&self, i: usize, j: usize) -> bool {
        self.compat[i * self.n + j]
    }
    #[inline]
    fn encap(&self, i: usize, j: usize) -> bool {
        self.encap[i * self.n + j]
    }

    /// Port of `assembleAlignments` (lines 72–308). Populates `self.assemblies`
    /// and `self.containment`.
    pub fn assemble(&mut self) -> Result<()> {
        if self.n == 0 {
            return Ok(());
        }

        self.determine_compatibilities_and_encapsulations();
        self.populate_lobjects();
        self.do_full_fscan();

        // First pass: the single top-scoring forward assembly.
        let top_indices = self.get_top_scoring_alignment();
        self.push_assembly(top_indices.clone());

        if top_indices.len() == self.n {
            return Ok(()); // everything assembled into one.
        }

        // Otherwise discard it and re-derive everything from scratch (this
        // re-finds the top assembly and handles ties uniformly). Port literally.
        self.assemblies.clear();
        self.containment.clear();

        self.do_full_rscan();

        let mut accounted = vec![false; self.n];

        // Compute combined_score + nucleated trace for every index.
        for i in 0..self.n {
            let cs = self.lobjects[i].lscore_f + self.lobjects[i].lscore_r
                - self.lobjects[i].num_contained as i64;
            self.lobjects[i].combined_score = cs;
            let trace = self.nucleate(i);
            self.lobjects[i].trace = trace;
        }

        // Bin indices by combined_score. The C++ sorts ascending then walks from
        // the end, grouping consecutive-equal scores → bins ordered high→low,
        // with within-bin order = reverse of the ascending array.
        let bins = self.bin_by_combined_score();

        for mut bin in bins {
            loop {
                let max_missing = self.get_max_missing(&bin, &accounted);
                let lobj_idx = match max_missing {
                    Some(idx) => idx,
                    None => break, // nothing missing in this bin → next bin
                };

                let trace = self.lobjects[lobj_idx].trace.clone();
                let has_unconsumed = trace.iter().any(|&t| !accounted[t]);

                if has_unconsumed {
                    self.push_assembly(trace.clone());
                    for &t in &trace {
                        accounted[t] = true;
                    }
                    if accounted.iter().all(|&b| b) {
                        return Ok(());
                    }
                } else {
                    break; // (defensive; unreachable when max_missing > 0)
                }

                bin.retain(|&x| x != lobj_idx);
            }
        }

        // Reached only if some alignment was never covered.
        Err(CombinrError::UncoveredAlignments(self.n))
    }

    /// All-vs-all compatibility + encapsulation (`determine_..`, 660–689).
    fn determine_compatibilities_and_encapsulations(&mut self) {
        let n = self.n;
        for i in 0..n {
            for j in (i + 1)..n {
                if can_merge(&self.alignments[i], &self.alignments[j], self.fuzzlength) {
                    self.compat[i * n + j] = true;
                    self.compat[j * n + i] = true;
                    if encapsulates(&self.alignments[i], &self.alignments[j]) {
                        self.encap[i * n + j] = true;
                    }
                    if encapsulates(&self.alignments[j], &self.alignments[i]) {
                        self.encap[j * n + i] = true;
                    }
                }
            }
        }
    }

    /// `populateLobjects` (692–706): contained = self ∪ everything self encapsulates.
    fn populate_lobjects(&mut self) {
        let n = self.n;
        let mut lobjects = Vec::with_capacity(n);
        for i in 0..n {
            let mut l = Lobject::new(i, n);
            let mut contained = vec![i];
            for j in 0..n {
                if self.encap(i, j) {
                    contained.push(j);
                }
            }
            l.set_contained_indices(&contained);
            lobjects.push(l);
        }
        self.lobjects = lobjects;
    }

    /// `do_full_Fscan` (554–596).
    fn do_full_fscan(&mut self) {
        for i in 1..self.n {
            let mut top_score = 0i64;
            let mut top_idx: Option<usize> = None;
            for j in (0..i).rev() {
                let compatible = self.compat(i, j);
                let containment = self.encap(i, j) || self.encap(j, i);
                if compatible && !containment {
                    let curr = self.lobjects[j].lscore_f
                        + self.lobjects[i].num_unique_contained(&self.lobjects[j]);
                    if curr > top_score {
                        top_idx = Some(j);
                        top_score = curr;
                    }
                }
            }
            if let Some(j) = top_idx {
                self.lobjects[i].from = Some(j);
                self.lobjects[i].lscore_f = top_score;
            }
        }
    }

    /// `do_full_Rscan` (600–638).
    fn do_full_rscan(&mut self) {
        for i in (0..self.n.saturating_sub(1)).rev() {
            let mut top_score = 0i64;
            let mut top_idx: Option<usize> = None;
            for j in (i + 1)..self.n {
                let compatible = self.compat(i, j);
                let containment = self.encap(i, j) || self.encap(j, i);
                if compatible && !containment {
                    let curr = self.lobjects[j].lscore_r
                        + self.lobjects[i].num_unique_contained(&self.lobjects[j]);
                    if curr > top_score {
                        top_idx = Some(j);
                        top_score = curr;
                    }
                }
            }
            if let Some(j) = top_idx {
                self.lobjects[i].to = Some(j);
                self.lobjects[i].lscore_r = top_score;
            }
        }
    }

    /// Walk the `from` chain from `start`, unioning contained indices; ascending.
    fn back_trace(&self, start: usize) -> Vec<usize> {
        let mut tracker = vec![false; self.n];
        let mut cur = Some(start);
        while let Some(idx) = cur {
            for bit in self.lobjects[idx].contained.ones() {
                tracker[bit] = true;
            }
            cur = self.lobjects[idx].from;
        }
        (0..self.n).filter(|&i| tracker[i]).collect()
    }

    /// Walk the `to` chain from `start`, unioning contained indices; ascending.
    fn forward_trace(&self, start: usize) -> Vec<usize> {
        let mut tracker = vec![false; self.n];
        let mut cur = Some(start);
        while let Some(idx) = cur {
            for bit in self.lobjects[idx].contained.ones() {
                tracker[bit] = true;
            }
            cur = self.lobjects[idx].to;
        }
        (0..self.n).filter(|&i| tracker[i]).collect()
    }

    /// `get_top_scoring_alignment` (789–805): max `lscore_f` (first wins), backtrace.
    fn get_top_scoring_alignment(&self) -> Vec<usize> {
        let mut top_score = 0i64;
        let mut top_idx: Option<usize> = None;
        for i in 0..self.n {
            let score = self.lobjects[i].lscore_f;
            if score > top_score {
                top_idx = Some(i);
                top_score = score;
            }
        }
        // top_idx is None only when all lscore_f == 0, impossible (>=1 each).
        self.back_trace(top_idx.unwrap_or(0))
    }

    /// `get_alignment_assembly_nucleating_at_alignment_index`: union of back and
    /// forward trace, ascending unique.
    fn nucleate(&self, index: usize) -> Vec<usize> {
        let mut tracker = vec![false; self.n];
        for i in self.back_trace(index) {
            tracker[i] = true;
        }
        for i in self.forward_trace(index) {
            tracker[i] = true;
        }
        (0..self.n).filter(|&i| tracker[i]).collect()
    }

    /// Bin indices by `combined_score`, high→low, mirroring the C++ binning.
    fn bin_by_combined_score(&self) -> Vec<Vec<usize>> {
        // ascending stable sort by combined_score (ties keep index order).
        let mut order: Vec<usize> = (0..self.n).collect();
        order.sort_by_key(|&i| self.lobjects[i].combined_score);

        let mut bins: Vec<Vec<usize>> = Vec::new();
        let mut curr_bin: Vec<usize> = Vec::new();
        let mut curr_score: Option<i64> = None;
        // walk from the end (highest score) toward the start.
        for &idx in order.iter().rev() {
            let score = self.lobjects[idx].combined_score;
            match curr_score {
                Some(s) if s == score => curr_bin.push(idx),
                _ => {
                    if curr_score.is_some() {
                        bins.push(std::mem::take(&mut curr_bin));
                    }
                    curr_bin.push(idx);
                    curr_score = Some(score);
                }
            }
        }
        if !curr_bin.is_empty() {
            bins.push(curr_bin);
        }
        bins
    }

    /// `get_max_missing_Lobj`: index in `bin` whose trace has the most
    /// not-yet-accounted members (strict `>`, first wins); `None` if all 0.
    fn get_max_missing(&self, bin: &[usize], accounted: &[bool]) -> Option<usize> {
        let mut max_missing = 0usize;
        let mut result = None;
        for &lidx in bin {
            let num_missing = self.lobjects[lidx]
                .trace
                .iter()
                .filter(|&&t| !accounted[t])
                .count();
            if num_missing > max_missing {
                max_missing = num_missing;
                result = Some(lidx);
            }
        }
        result
    }

    /// `create_assembly` (766–786): fold member alignments left-to-right.
    fn create_assembly(&self, indices: &[usize]) -> Alignment {
        debug_assert!(!indices.is_empty());
        let mut idx = indices.to_vec();
        idx.sort_unstable();
        let mut assembly = self.alignments[idx[0]].clone();
        for &k in &idx[1..] {
            assembly = merge_alignments(&assembly, &self.alignments[k]);
        }
        assembly
    }

    /// Create + title an assembly and record its membership.
    fn push_assembly(&mut self, indices: Vec<usize>) {
        let mut assembly = self.create_assembly(&indices);
        assembly.acc = format!("assembly_{}", self.assemblies.len());
        self.assemblies.push(assembly);
        self.containment.push(indices);
    }

    /// Emit the `pasa`-compatible assembly lines (golden-diff target). One line
    /// per assembly:
    /// `assembly: (i) contains alignments: [accs] with structure [struct] score: (N)`.
    pub fn format_pasa_assemblies(&self) -> String {
        let mut out = String::new();
        for (i, asm) in self.assemblies.iter().enumerate() {
            let members = &self.containment[i];
            let accs: Vec<&str> = members
                .iter()
                .map(|&m| self.alignments[m].acc.as_str())
                .collect();
            out.push_str(&format!(
                "assembly: ({}) contains alignments: [{}] with structure [{}] score: ({})\n",
                i,
                accs.join(","),
                structure_string(asm),
                members.len()
            ));
        }
        out
    }
}

/// `title,orient,lend-rend,...` — the C++ `assembly->toString()` form.
fn structure_string(asm: &Alignment) -> String {
    let orient = match asm.aligned_orient {
        Strand::Plus => '+',
        Strand::Minus => '-',
        Strand::Unknown => '?',
    };
    let mut parts = vec![asm.acc.clone(), orient.to_string()];
    for s in &asm.segments {
        parts.push(format!("{}-{}", s.coords.lend, s.coords.rend));
    }
    parts.join(",")
}
