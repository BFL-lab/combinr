//! The `Lobject` dynamic-programming record.
//!
//! Port of `pasa_cpp/Lobject.{h,cpp}`. `contained` holds the alignment indices
//! this object subsumes (itself + everything it encapsulates); `lscore_f` /
//! `lscore_r` are the forward/reverse chain scores; `from` / `to` are the chain
//! links (indices, not raw pointers); `trace` caches the nucleated membership.

use fixedbitset::FixedBitSet;

#[derive(Debug, Clone)]
pub struct Lobject {
    pub index: usize,
    pub contained: FixedBitSet,
    pub num_contained: usize,
    pub lscore_f: i64,
    pub lscore_r: i64,
    pub combined_score: i64,
    pub from: Option<usize>,
    pub to: Option<usize>,
    pub trace: Vec<usize>,
}

impl Lobject {
    pub fn new(index: usize, num_alignments: usize) -> Self {
        Lobject {
            index,
            contained: FixedBitSet::with_capacity(num_alignments),
            num_contained: 0,
            lscore_f: 0,
            lscore_r: 0,
            combined_score: 0,
            from: None,
            to: None,
            trace: Vec::new(),
        }
    }

    /// Port of `setContainedIndices`: record the contained indices and set
    /// `lscore_f = lscore_r = num_contained = |indices|` (the triple increment).
    pub fn set_contained_indices(&mut self, indices: &[usize]) {
        self.num_contained = 0;
        for &i in indices {
            self.contained.insert(i);
            self.lscore_f += 1;
            self.lscore_r += 1;
            self.num_contained += 1;
        }
    }

    /// Port of `num_unique_contained`: number of indices contained in `self`
    /// but not in `other`.
    pub fn num_unique_contained(&self, other: &Lobject) -> i64 {
        self.contained
            .ones()
            .filter(|&i| !other.contained.contains(i))
            .count() as i64
    }
}
