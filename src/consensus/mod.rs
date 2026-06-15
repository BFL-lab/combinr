//! EVidenceModeler-style weighted-evidence gene consensus (`combinr consensus`).
//!
//! A run mode distinct from the PASA-derived `assemble` path: per locus, integrate
//! weighted heterogeneous evidence (ab-initio gene predictions, protein alignments,
//! transcript alignments, other homology predictions) into a single best-scoring
//! consensus gene structure via a frame-aware gene-structure DP trellis. Ported from
//! EVidenceModeler's `EvmUtils/evidence_modeler.pl`, but with EVM's sliding-window
//! partitioning and canonical-site bias removed (see `consensus-evm-port` and
//! `avoid-canonical-splice-bias` notes, and the project plan).
//!
//! Milestone status: **M0** — scaffolding only (weights parsing, evidence ingestion,
//! per-contig region clustering). The trellis engine lands in later milestones.

pub mod candidates;
pub mod engine;
pub mod evidence;
pub mod exon;
pub mod filter;
pub mod grammar;
pub mod output;
pub mod promote;
pub mod region;
pub mod repeats;
pub mod sites;
pub mod trellis;
pub mod vectors;
pub mod weights;

pub use candidates::{CandidateParams, RegionData, build_candidates};
pub use engine::{CalledGene, EngineParams, consensus_region};
pub use evidence::EvidenceChain;
pub use exon::{ExonCandidate, ExonType};
pub use filter::{FilterParams, SupportFlags};
pub use output::to_out_genes;
pub use region::ConsensusRegion;
pub use sites::SiteSets;
pub use trellis::{ConsensusGene, run_trellis, score_all_exons};
pub use vectors::{IntronScores, RegionVectors};
pub use weights::{EvClass, Weights};
