# combinr

`combinr` reimplements only the evidence-combining algorithms from
[PASA](https://github.com/PASApipeline/PASApipeline) and
[EVidenceModeler](https://github.com/EVidenceModeler/EVidenceModeler) as a
single self-contained Rust binary. This reimplementation was made to be agnostic
to input source, and also accepts a wide variety of input file types,
specifically GTF, GFF3, or BAM.

1. **Combine** multiple transcript-alignment sources into one **non-redundant**
   assembly set (`assemble`). A faithful port of PASA's C++ assembler and
   its orientation wrapper.
2. **Model alternative splicing** between the resulting isoforms
   (`assemble --alt-splice`: retained intron, alternate donor/acceptor,
   transcription-start / polyA within an intron, exon skipping, alternate terminal
   exons); `assemble --gene-pred --genome` additionally reconciles an external
   gene-prediction CDS onto the isoforms to derive CDS + 5'/3' UTRs.
3. **Build consensus gene models** (`consensus`) by integrating *weighted*
   heterogeneous evidence (ab-initio gene predictions, protein alignments,
   transcript alignments) into one best-scoring coding gene structure per locus,
   via a frame-aware gene-structure dynamic program. A reimplementation of
   EVidenceModeler's core algorithm.

The two families are complementary, named by what they consolidate: `assemble`
*merges compatible alignments and keeps alternatives* (transcript-level isoforms);
`consensus` *weighs competing evidence and selects one model* (a gene-level
coding structure).

## Design notes

- **No canonical splice-site bias.** combinr deliberately omits PASA's
  sequence-based splice validation, which matches intron dinucleotides against
  GT‑AG/GC‑AG/AT‑AC to infer strand and drop "non-consensus" alignments. Splice
  junctions are defined purely by exon position, and strand is taken verbatim
  from the input. Whatever introns and strand your upstream tools produce are
  honored as-is. **This extends to `consensus`:** EVM normally scans the genome
  for canonical GT‑AG splice sites and ATG/stop codons to gate candidate exons
  whereas combinr does not. The candidate splice/start/stop sites are *evidence-derived*
  (the union of intron boundaries and CDS bounds seen across the inputs), and
  in-frame stops are detected with the configurable genetic code, never a hard-coded
  table.
- **The genome FASTA is only needed for the ORF/UTR step.** Assembly and
  alt-splice classification are purely coordinate-based.
- **ORFs come from an external gene prediction, not from translating the
  transcripts.** Given a gene-prediction GFF3 with already-mapped CDS, an
  isoform that matches the prediction inherits the CDS verbatim; a divergent
  isoform is projected from the *same* start codon and translated along its own
  exon structure to the first in-frame stop (a premature stop becomes that
  isoform's ORF stop). Multiple coding models per locus are supported.
- **Genetic code.** Divergent-isoform stop detection uses the standard code
  (NCBI table 1) by default; pass `--genetic-code/-g <id>` to the `assemble`
  reconcile step or `consensus` for a non-standard table. Every assigned NCBI
  table is supported: `1-6`, `9-16`, `21-33` (ids 7, 8, and 17–20 were never
  assigned and are rejected). Codes are grouped by stop-codon set, since only stop
  assignments affect ORF bounds; matching isoforms inherit the predicted CDS and
  are unaffected by the code. The context-dependent tables (27/28/31, whose codons
  are dual sense/stop) use NCBI's declared stop set. Translation halts at the
  first such codon, the conservative choice for CDS bounds.
- **Consensus never blanks a locus (`consensus`).** EVM eliminates a gene whose
  coding/noncoding score ratio falls below a threshold, which can leave a locus empty.
  combinr computes the same metric but **flags** low-support genes (`low_support=true`)
  and keeps them by default. A false negative is treated as worse than a false
  positive. Pass `--strict` for EVM's drop-it behavior. (combinr also fixes EVM's bug
  where an *eliminated* gene's span still blocked the re-search for nested/adjacent
  genes.)
- **Both strands, independently.** `consensus` runs the trellis on each strand by
  reverse-complementing the locus, so the two strands are explored fully everywhere;
  an antisense locus can yield genes on both strands rather than EVM's single
  best-strand pick. The genome FASTA is required here (for the across-junction stop
  check and CDS projection) but, per the note above, never gates which sites are allowed.

## Install

`combinr` is a single self-contained binary with **no system dependencies**: no
Perl, no database, and no htslib (BAM reading comes from the pure-Rust
[`noodles`](https://github.com/zaeleus/noodles) crate). It runs on Linux and
macOS. Install a prebuilt binary (no toolchain required), or build from source.

### Prebuilt binary (recommended)

The quickest option, and it needs no Rust toolchain. This downloads a prebuilt,
statically-linked binary for your platform, verifies its checksum, installs it,
and updates your `PATH`:

```sh
curl --proto '=https' --tlsv1.2 -LsSf https://github.com/BFL-lab/combinr/releases/latest/download/combinr-installer.sh | sh
```

`releases/latest/download/…` always resolves to the newest release; pin a
specific version by using its tag instead:

```sh
curl --proto '=https' --tlsv1.2 -LsSf https://github.com/BFL-lab/combinr/releases/download/v0.1.1/combinr-installer.sh | sh
```

Prefer to do it by hand (or on Windows)? Grab the tarball for your target — or the
`.zip` on Windows — from the [Releases page](https://github.com/BFL-lab/combinr/releases),
check it against the bundled `.sha256`, and put `combinr` on your `PATH`.

### Build from source

Building needs a **Rust toolchain, 1.85 or newer** (the crate uses the 2024
edition). The easiest way to get one is [rustup](https://rustup.rs):

```sh
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh   # then: rustc --version
```

Then clone and build:

```sh
git clone https://github.com/BFL-lab/combinr
cd combinr
cargo build --release
# optimized binary lands at: target/release/combinr
./target/release/combinr --help
```

The first build compiles all dependencies and takes a few minutes; later builds
are incremental. Drop `--release` for a faster-compiling debug build under
`target/debug/`. To install it onto your `PATH` (into `~/.cargo/bin`):

```sh
cargo install --path .                                   # from a local checkout
cargo install --git https://github.com/BFL-lab/combinr   # or straight from git
combinr --version
```

### Run the tests (optional)

The suite is self-contained. It checks combinr against committed golden
references, so no PASA or EVidenceModeler installation is required:

```sh
cargo test
```

## Usage

```sh
# 1) Combine sources into a non-redundant GFF3 (multiple -i, GTF/GFF3/BAM, mixed)
combinr assemble -i sampleA.gtf -i sampleB.gff3 -i reads.bam > assemblies.gff3

# 2) Also classify alternative splicing (isoform GFF3 to stdout, events to a TSV)
combinr assemble -i sampleA.gtf -i sampleB.gff3 --alt-splice --events events.tsv > isoforms.gff3

# 3) Reconcile an external CDS into CDS/UTR (add --genetic-code/-g for a non-standard table)
combinr assemble -i samples.gtf --gene-pred p.gff3 --genome genome.fa --events events.tsv > out.gff3

# 4) Consensus gene models from weighted evidence (EvidenceModeler-style)
combinr consensus --weights weights.txt --genome genome.fa \
    --gene-predictions abinitio.gff3 \
    --protein-alignments proteins.gff3 \
    --transcript-alignments transcripts.gff3 > consensus.gff3
```

`consensus` options: `--strict` (drop low-support genes instead of flagging),
`--genetic-code/-g`, `--flank <bp>` (per-locus DP-window padding, default 10000.
Only widens the window for UTR/end placement, never merges loci), and the
off-by-default EVM heuristics `--research-intergenic <bp>` (re-search intergenic gaps
for missed genes), `--search-long-introns <bp>` (find nested genes in long introns),
`--extend-terminal-stop` (extend protein/transcript 3' ends to a stop to complete genes),
`--peak-augment` (boost intergenic scores near start/stop evidence peaks),
`--promote-transcript-orfs` (at loci with only transcript evidence and no consensus CDS,
find the longest ORF in the transcripts and emit it, tagged `support=transcript_orf`),
`--alt-splice` (also emit each consensus locus's alternative transcript isoforms as extra
mRNAs, with CDS derived from the consensus (inherited if matching, re-projected from the
consensus start if divergent) and write a region-tagged alt-splice events TSV to
`--events`), plus `--repeats <gff3>` (mask repeats from scoring).

`assemble` tuning: `--fuzzlength <bp>` (default 20) and the quality filters
`--min-avg-per-id` and `--min-intron` (off by default) plus `--max-intron`
(defaults to 100000 bp, matching PASA's `MAX_INTRON_LENGTH`, to drop spurious
long-range junctions; pass `--max-intron 0` to disable the cap).

Shared options, given *after* the subcommand name (e.g. `combinr consensus
--threads 4 ...`): `--output/-o <file>` (write models to a file instead of stdout,
the default), `--format gff3|gtf`, `--threads <n>` (default 4; `0` = all cores),
and `--verbose`. The `> out.gff3` redirections above are interchangeable with
`-o out.gff3`.

## Input

- **Transcript files** (`-i`, repeatable; GTF, GFF3, and/or BAM, mixed allowed).
  GFF3 is accepted as `cDNA_match`+`Target` alignment rows or as plain
  `gene`/`mRNA`/`exon` models; GTF as `exon` rows grouped by `transcript_id`.
  Each file's basename is recorded as provenance.
- **BAM** (auto-detected by the `.bam` extension): a cDNA/transcript-to-genome
  alignment from a splice-aware mapper (segemehl, `minimap2 -ax splice`, STAR,
  HISAT2, …). One primary mapped record becomes one transcript; exon blocks come
  from the CIGAR (a `N` skip is an intron; `D`/`I` stay within an exon); the
  transcribed strand is read from the `XS:A` tag (else left undetermined and
  resolved during assembly). CRAM and plain-text SAM are not supported.
- **Gene-prediction GFF3** (`--gene-pred`, the `assemble` reconcile step only): `CDS`
  rows grouped by mRNA `Parent`.
- **Genome FASTA** (`--genome`, the `assemble` reconcile step and `consensus`).
- **Consensus evidence** (`consensus`): a `--weights` file (three whitespace columns
  `CLASS TYPE WEIGHT`, where `CLASS` is `PROTEIN`/`TRANSCRIPT`/`ABINITIO_PREDICTION`/
  `OTHER_PREDICTION` and `TYPE` matches the GFF column-2 source) plus the evidence GFF3s:
  `--gene-predictions` (read as `CDS` rows grouped by `Parent`), `--protein-alignments`
  and `--transcript-alignments` (read as match chains carrying `Target=`). All repeatable.

## Output

- **Transcript models** to stdout as GFF3 (`gene`/`mRNA`/`exon`, plus
  `CDS`/`five_prime_UTR`/`three_prime_UTR` when the ORF step runs) or GTF
  (`--format gtf`). Each transcript carries `sources=` (provenance) and
  `contains=` (contributing accessions).
- **Alt-splice event report** as a TSV (`--events`) with the event type, coords,
  the two isoforms, and after the ORF step, the region (5'UTR / CDS / 3'UTR).
- **Consensus gene models** (`consensus`) to stdout as GFF3 `gene`/`mRNA`/`exon`/`CDS`
  with CDS phases (the exon structure is the coding structure). Each mRNA carries
  `score`, `score_ratio`, `coding_length`, `low_support`, and `partial5`/`partial3`.

## Correctness

combinr is validated against the original PASA, using committed golden
references so the test suite is self-contained, **no PASA code runs at test
time**. Just `cargo test`.

- **Assembly (Algorithm 1)**:  the original C++ PASA binary's output on every
  `pasa_cpp_sample_input*` is committed under `tests/data/assembler/`.
  `tests/golden_assembler.rs` runs three combinr code paths (raw assembler,
  orientation wrapper, full GFF3 pipeline) and asserts each reproduces the golden
  assembly set.
- **Alternative splicing (Algorithm 2)**:  PASA's real
  `CDNA::Alternative_splice_comparer` was run on combinr's emitted isoforms and
  its events committed under `tests/data/altsplice/` (the `stringtie.gtf` fixture
  alone exercises all seven event types). `tests/golden_altsplice.rs` asserts
  combinr's classification reproduces them exactly.
- **Consensus**: because `consensus` deliberately diverges from EVM (non-canonical
  sites + two independent per-strand trellises), there is no byte-for-byte EVM parity
  to assert. Instead `tests/golden_consensus.rs` is a **self-golden**: it runs
  `consensus` on EVidenceModeler's `testing/` data set (vendored under
  `tests/data/consensus/`) and asserts the GFF3 matches a committed frozen snapshot,
  and that the output is deterministic. This guards against regressions in combinr's
  own algorithm.

The PASA goldens are a frozen snapshot, generated once from a PASA checkout; the
one-off generation tooling is preserved in git history. combinr has also been run on the full
`sample_data` (a 6.5 MB cDNA_match GFF3 → 1334 assemblies; with the bundled gene
annotations + genome → CDS/UTR-annotated isoforms; alt-splice cross-check: 371/371
events identical to PASA).
