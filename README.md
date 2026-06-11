# combinr

`combinr` reimplements the two genuinely useful algorithms from
[PASA](https://github.com/PASApipeline/PASApipeline) as a single self-contained
Rust binary — no database, no Trinity requirement, and input from any source as
long as it is GTF or GFF3:

1. **Combine** multiple transcript-alignment sources into one **non-redundant**
   assembly set. A faithful port of PASA's C++ `pasa` assembler and its
   orientation wrapper.
2. **Model alternative splicing** between the resulting isoforms (retained
   intron, alternate donor/acceptor, transcription-start / polyA within an
   intron, exon skipping, alternate terminal exons), with an **optional** step
   that reconciles an external gene-prediction CDS onto the isoforms to derive
   CDS + 5'/3' UTRs.

## Design notes

- **No canonical splice-site bias.** combinr deliberately omits PASA's
  sequence-based splice validation, which matches intron dinucleotides against
  GT‑AG/GC‑AG/AT‑AC to infer strand and drop "non-consensus" alignments. Splice
  junctions are defined purely by exon position, and strand is taken verbatim
  from the input. Whatever introns and strand your upstream tools produce are
  honored as-is.
- **The genome FASTA is only needed for the ORF/UTR step.** Assembly and
  alt-splice classification are purely coordinate-based.
- **ORFs come from an external gene prediction, not from translating the
  transcripts.** Given a gene-prediction GFF3 with already-mapped CDS, an
  isoform that matches the prediction inherits the CDS verbatim; a divergent
  isoform is projected from the *same* start codon and translated along its own
  exon structure to the first in-frame stop (a premature stop becomes that
  isoform's ORF stop). Multiple coding models per locus are supported.

## Build

```sh
cargo build --release
# binary: target/release/combinr
```

## Usage

```sh
# 1) Combine sources into a non-redundant GFF3 (multiple -i, GTF and/or GFF3)
combinr assemble -i sampleA.gtf -i sampleB.gff3 > assemblies.gff3

# 2) Model alternative splicing (isoform GFF3 to stdout, events to a TSV)
combinr altsplice -i sampleA.gtf -i sampleB.gff3 --events events.tsv > isoforms.gff3

# 3) Reconcile an external CDS prediction into CDS/UTR annotations
combinr orf -i samples.gtf --gene-pred predictions.gff3 --genome genome.fa \
    --events events.tsv > annotated.gff3

# 4) Full pipeline (does the ORF step when --gene-pred and --genome are given)
combinr run -i samples.gtf [--gene-pred p.gff3 --genome genome.fa] --events events.tsv > out.gff3
```

Global options: `--format gff3|gtf`, `--fuzzlength <bp>` (default 20),
`--threads <n>`, and the off-by-default quality filters `--min-avg-per-id`,
`--min-intron`, `--max-intron`.

## Input

- **Transcript files** (`-i`, repeatable; GTF and/or GFF3, mixed allowed). GFF3
  is accepted as `cDNA_match`+`Target` alignment rows or as plain
  `gene`/`mRNA`/`exon` models; GTF as `exon` rows grouped by `transcript_id`.
  Each file's basename is recorded as provenance.
- **Gene-prediction GFF3** (`--gene-pred`, ORF step only): `CDS` rows grouped by
  mRNA `Parent`.
- **Genome FASTA** (`--genome`, ORF step only).

## Output

- **Transcript models** to stdout as GFF3 (`gene`/`mRNA`/`exon`, plus
  `CDS`/`five_prime_UTR`/`three_prime_UTR` when the ORF step runs) or GTF
  (`--format gtf`). Each transcript carries `sources=` (provenance) and
  `contains=` (contributing accessions).
- **Alt-splice event report** as a TSV (`--events`) with the event type, coords,
  the two isoforms, and — after the ORF step — the region (5'UTR / CDS / 3'UTR).

## Correctness

combinr is validated against the original PASA, using committed golden
references so the test suite is self-contained — **no PASA code runs at test
time**. Just `cargo test`.

- **Assembly (Algorithm 1)** — the original C++ `pasa` binary's output on every
  `pasa_cpp_sample_input*` is committed under `tests/data/assembler/`.
  `tests/golden_assembler.rs` runs three combinr code paths (raw assembler,
  orientation wrapper, full GFF3 pipeline) and asserts each reproduces the golden
  assembly set.
- **Alternative splicing (Algorithm 2)** — PASA's real
  `CDNA::Alternative_splice_comparer` was run on combinr's emitted isoforms and
  its events committed under `tests/data/altsplice/` (the `stringtie.gtf` fixture
  alone exercises all seven event types). `tests/golden_altsplice.rs` asserts
  combinr's classification reproduces them exactly.

The goldens are a frozen snapshot, generated once from a PASA checkout; the
one-off generation tooling is preserved in git history. combinr has also been run on the full
`sample_data` (a 6.5 MB cDNA_match GFF3 → 1334 assemblies; with the bundled gene
annotations + genome → CDS/UTR-annotated isoforms; alt-splice cross-check: 371/371
events identical to PASA).
