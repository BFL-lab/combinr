# Validation: regenerating the PASA golden references

The `tests/golden_*` suite compares combinr against **committed** PASA reference
outputs under `tests/data/` — so `cargo test` runs no PASA code. This directory
holds the tooling that *generated* those goldens from a PASA checkout:

- `regenerate_goldens.sh` — rebuilds every golden under `tests/data/` (the C++
  `pasa` assembler output for Algorithm 1, and the alt-splice events for
  Algorithm 2). Run `PASA=/path/to/PASApipeline ./validation/regenerate_goldens.sh`
  from the crate root.
- `pasa_altsplice_xcheck.pl` — drives PASA's real
  `CDNA::Alternative_splice_comparer` (the module that produces PASA's alt-splice
  report) on combinr's emitted isoforms; used both to generate the alt-splice
  golden and for ad-hoc cross-checks on larger datasets.

## How the alt-splice cross-check works

PASA's alt-splice analysis is normally database-bound (`.dbi` scripts), so there
is no standalone binary to diff against. But the classification *logic* lives in
a pure-Perl module that loads without a database, which the harness drives
directly:

1. Run `combinr altsplice` on some input, producing the isoform GFF3 and the
   event TSV.
2. The harness parses the isoform GFF3, rebuilds each isoform as a PASA
   `Gene_obj` (via `build_gene_obj_exons_n_cds_range`), and runs PASA's five
   finders (`find_unspliced_introns`, `find_conventional_alt_splice_isoforms`,
   `find_starts_and_ends_within_introns`, `find_exon_skipping_events`,
   `find_alternate_exons`) on every isoform pair **in the same order combinr
   uses**.
3. Both event sets are reduced to `(event_type, coords, isoform_a, isoform_b)`
   and compared.

Because both sides operate on identical isoform structures, this isolates the
alt-splice classification from any assembly differences.

## Running it

```sh
PL=/path/to/PASApipeline/PerlLib
S=/path/to/PASApipeline/sample_data

combinr altsplice -i "$S/stringtie.gtf" --events combinr_events.tsv > isoforms.gff3

perl validation/pasa_altsplice_xcheck.pl "$PL" isoforms.gff3 \
  | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4}' | sort -u > pasa.canon
tail -n +2 combinr_events.tsv \
  | awk -F'\t' '{print $3"\t"$5"\t"$6"\t"$7}' | sort -u > combinr.canon

comm -3 combinr.canon pasa.canon   # empty == perfect agreement
```

## Result

Run on PASA's bundled `sample_data` (`v2.5.3`):

| input | isoforms / loci | events | combinr-only | pasa-only |
|---|---|---|---|---|
| `stringtie.gtf` | 573 / 503 | 107 | 0 | 0 |
| `custom_alignments.gff3` | 1239 / 1034 | 371 | 0 | 0 |

**Exact agreement** across all seven event types (retained intron, alt
acceptor, alt donor, start/end-within-intron, exon skip, alternate exon) on
1,812 isoforms and 478 events — zero discrepancies.
