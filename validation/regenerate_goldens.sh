#!/usr/bin/env bash
# Regenerate the committed golden references under tests/data/ from a PASA
# checkout. Run from the crate root. Requires: the built C++ `pasa` binary, perl
# with the PASA PerlLib, and a release build of combinr.
#
#   PASA=/path/to/PASApipeline ./validation/regenerate_goldens.sh
#
# The committed goldens are produced by exactly these steps; CI/tests never run
# any of this — they compare combinr against the committed output.

set -euo pipefail

PASA="${PASA:-/home/matt/PASApipeline}"
PASA_CPP="$PASA/pasa_cpp"
PERLLIB="$PASA/PerlLib"

[ -x "$PASA_CPP/pasa" ] || { echo "build the pasa binary first: ( cd $PASA_CPP && make )"; exit 1; }

cargo build --release

# --- Algorithm 1: assembler golden (C++ pasa output per sample input) ---
mkdir -p tests/data/assembler
for f in "$PASA_CPP"/pasa_cpp_sample_input*; do
  case "$f" in *.o|*.output) continue;; esac
  name=$(basename "$f")
  cp "$f" "tests/data/assembler/$name"
  "$PASA_CPP/pasa" "$f" 2>/dev/null | grep '^assembly:' > "tests/data/assembler/$name.golden"
done

# --- Algorithm 2: alt-splice golden (PASA comparer on combinr's isoforms) ---
mkdir -p tests/data/altsplice
for base in tests/data/altsplice/*.gtf tests/data/altsplice/*.gff3; do
  [ -e "$base" ] || continue
  stem="${base%.*}"
  ./target/release/combinr altsplice -i "$base" --events /tmp/regen.events > /tmp/regen.iso.gff3 2>/dev/null
  perl validation/pasa_altsplice_xcheck.pl "$PERLLIB" /tmp/regen.iso.gff3 2>/dev/null \
    | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4}' | sort -u > "$stem.events.golden"
done

echo "Regenerated golden references under tests/data/."
