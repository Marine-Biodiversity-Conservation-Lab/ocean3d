#!/usr/bin/env bash
# Run bench.R for every reader x step on one dataset, one process each, and
# append "dataset reader step seconds peak_MB status result" rows to
# dev/data-model-trial/results/benchmarks.tsv.
#
# Usage (from the package root):
#   dev/data-model-trial/run_bench.sh <lib> <data_dir> <dataset> "<readers>" "<steps>"
set -u
lib="$1"; data_dir="$2"; dataset="$3"; readers="$4"; steps="$5"
here="dev/data-model-trial"
out="$here/results/benchmarks.tsv"
mkdir -p "$here/results"
[ -f "$out" ] || printf 'dataset\treader\tstep\tseconds\tpeak_MB\tstatus\tresult\n' > "$out"

for reader in $readers; do
  for step in $steps; do
    log="$(mktemp)"
    "$here/measure.sh" "$dataset/$reader/$step" \
      Rscript "$here/bench.R" "$lib" "$data_dir" "$dataset" "$reader" "$step" \
      > "$log" 2>&1
    m="$(grep '^MEASURE' "$log" | tail -n 1)"
    r="$(grep '^RESULT' "$log" | tail -n 1 | cut -f5-)"
    [ -n "$r" ] || r="FAILED: $(grep -v '^MEASURE' "$log" | tail -n 3 | tr '\n\t' '  ')"
    secs="$(cut -f3 <<< "$m")"; mb="$(cut -f4 <<< "$m")"; code="$(cut -f5 <<< "$m")"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$dataset" "$reader" "$step" \
      "$secs" "$mb" "$code" "$r" | tee -a "$out"
    rm -f "$log"
  done
done
