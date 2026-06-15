# CLAUDE.md

Guidance for working in `combinr`. See `README.md` for user-facing docs.

## What this is

A self-contained Rust binary that reimplements evidence-combining genome-annotation
algorithms from PASA and EVidenceModeler — no Perl, no database. Two families, named by
the consolidation they perform:

- **`assemble` / `altsplice` / `orf` / `run`** — PASA: merge compatible transcript
  alignments into non-redundant isoforms (keeps alternatives), classify alt-splice
  events, and optionally graft an external CDS onto isoforms for CDS/UTR.
- **`consensus`** — EVidenceModeler: integrate *weighted* heterogeneous evidence
  (ab-initio predictions, protein/transcript alignments) into one best-scoring coding
  gene model per locus, via a frame-aware gene-structure DP trellis.

## Conventions (hold across the codebase)

- **Coordinates are 1-based, inclusive, `lend <= rend`** (`model::Coordset`). Strand is
  `Plus` / `Minus` / `Unknown`.
- **No canonical splice-site bias** (the project's defining principle). Splice junctions
  are positional, strand is verbatim from input, and `consensus` derives candidate
  splice/start/stop sites from the *evidence* (union of observed intron boundaries + CDS
  bounds), never a genome GT-AG/ATG scan. In-frame stops use a configurable `GeneticCode`
  (NCBI tables 1, 4, 6, 10, 12, 26).
- Idiomatic Rust plumbing (clap, thiserror, fixedbitset, rayon, hand-rolled parsers);
  faithful algorithm cores.

## Layout

- `src/model.rs` — `Coordset`, `Strand`, `Segment`, `Alignment` (the shared data model).
- `src/io/` — GFF3/GTF/BAM parsers (`load_sources`), FASTA, GFF3/GTF writers, `out_model`
  (`OutGene`/`OutTranscript`, consumed by both writers).
- `src/cluster.rs` — single-linkage overlap clustering (`cluster_alignments`;
  `cluster_spans_by_contig` for both-strand consensus regions).
- `src/assemble/`, `src/altsplice/`, `src/orf/` — the PASA path (`orf::SplicedTranscript`,
  `project_orf`, `GeneticCode`, `reconcile` are reused by `consensus`).
- `src/consensus/` — the EVM port (see below).
- `src/pipeline.rs` — top-level orchestration fns (`assemble_sources`, `analyze_sources`,
  `reconcile_sources`, `consensus_sources`); rayon-parallel per cluster/region.
- `src/cli.rs` + `src/main.rs` — clap CLI and dispatch.

### `src/consensus/` (the EVM port)

`weights` (EvClass + weights file) → `evidence` (self-contained GFF3 → `EvidenceChain`,
keyed on GFF **column-2 source**; predictions read from `CDS` rows, alignments from
`Target=` match chains) → `candidates` (per-base coding vector + evidence introns/sites +
frame-aware `ExonCandidate`s; intergenic population; opt-in terminal-stop extension +
peak augmentation) → `grammar` + `trellis` (the gene-structure DP: `classify`,
`score_exon`, `are_compatible_exons`, `run_trellis`) → `engine` (both strands via
reverse-complement-then-transpose; tail/intergenic/long-intron recursion; CDS/UTR
projection) → `filter` (flag-not-drop low-support) → `output` (→ `OutGene`). `repeats`
parses the optional mask.

## Deliberate divergences from EVM (don't "fix" these)

- Introns are keyed on the **actual gap span** `(A.rend+1, B.lend-1)`, not EVM's canonical
  `lend-2` dinucleotide convention — internally self-consistent for the non-canonical model.
- **Two independent per-strand trellises**, not EVM's single combined one, so an antisense
  locus can yield genes on both strands (anti-false-negative).
- **Flag-not-drop** low-support filtering by default (`--strict` for EVM's drop). Also
  fixes EVM's bug where an eliminated gene's span blocked re-search.
- **End-stitching and 5'/start terminal extension are intentionally NOT ported** — they
  rely on canonical splice/ATG scanning that conflicts with the project principle.
- The standalone `orf` subcommand is **kept** (it grafts a CDS onto PASA-assembled
  alt-isoforms — distinct from EVM consensus; removing it would regress).

## Testing

`cargo test`. Goldens are committed frozen snapshots under `tests/data/<category>/`; no
PASA/EVM code runs at test time (`tests/common/mod.rs` has the helpers). The consensus
self-golden (`tests/golden_consensus.rs`) runs on EVidenceModeler's `testing/` data
vendored under `tests/data/consensus/`; regenerate `smalltest.consensus.gff3.golden` from
the binary if a change to the consensus output is intentional. Run on real data:
`combinr consensus --weights tests/data/consensus/weights.txt --genome
tests/data/consensus/genome.fasta --gene-predictions tests/data/consensus/gene_predictions.gff3
--protein-alignments tests/data/consensus/protein_alignments.gff3
--transcript-alignments tests/data/consensus/transcript_alignments.gff3`.
