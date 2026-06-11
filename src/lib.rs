//! combinr — PASA's core alignment-assembly and alternative-splicing algorithms,
//! reimplemented as a self-contained Rust library + binary.
//!
//! See the project plan and `combinr-project` notes. The two algorithms are:
//! 1. combine multiple transcript-alignment sources into one non-redundant
//!    assembly set (port of the C++ `pasa` assembler + its orientation wrapper);
//! 2. model alternative splicing (port of `Alternative_splice_comparer.pm`),
//!    with an optional CDS/UTR reconciliation step.

pub mod altsplice;
pub mod assemble;
pub mod cluster;
pub mod error;
pub mod filter;
pub mod io;
pub mod model;
pub mod orf;
pub mod pipeline;
pub mod token;

pub use assemble::Assembler;
pub use error::{CombinrError, Result};
pub use model::{Alignment, Coordset, Provenance, SegType, Segment, Strand};
