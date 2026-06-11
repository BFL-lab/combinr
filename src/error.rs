//! Typed library errors. Variants map to the `die()`/`confess()` sites in the
//! PASA source so failures stay diagnosable.

use thiserror::Error;

#[derive(Debug, Error)]
pub enum CombinrError {
    #[error("parse error in {file}: line {line}: {msg}")]
    Parse {
        file: String,
        line: usize,
        msg: String,
    },

    #[error("conflicting spliced orientations in assembly (have {have}, {acc} has {found})")]
    ConflictingSplicedOrientation {
        have: char,
        acc: String,
        found: char,
    },

    /// The greedy cover failed to account for every input alignment — mirrors
    /// `assembleAlignments`'s terminal "Not all alignments were accounted for".
    #[error("assembly cover did not account for all {0} alignment(s) in a cluster")]
    UncoveredAlignments(usize),

    #[error("accession {0:?} contains a comma, which the token format disallows")]
    CommaInAccession(String),

    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
}

pub type Result<T> = std::result::Result<T, CombinrError>;
