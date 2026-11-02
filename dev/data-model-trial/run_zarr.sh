#!/usr/bin/env bash
# Run each zarr.R step in its own process and append
# "zarr stars_mdim step seconds peak_MB status result" rows to
# dev/data-model-trial/results/benchmarks.tsv.
# Usage (from the package root): dev/data-model-trial/run_zarr.sh <lib> "<steps>"
set -u
lib="$1"; steps="$2"
here="dev/data-model-trial"
out="$here/results/benchmarks.tsv"
for step in $steps; do
  log="$(mktemp)"
  "$here/measure.sh" "zarr/$step" Rscript "$here/zarr.R" "$lib" "$step" > "$log" 2>&1
  m="$(grep '^MEASURE' "$log" | tail -n 1)"
  r="$(grep '^RESULT' "$log" | tail -n 1 | cut -f5-)"
  [ -n "$r" ] || r="FAILED: $(grep -v '^MEASURE' "$log" | tail -n 3 | tr '\n\t' '  ')"
  printf 'zarr\tstars_mdim\t%s\t%s\t%s\t%s\t%s\n' "$step" \
    "$(cut -f3 <<< "$m")" "$(cut -f4 <<< "$m")" "$(cut -f5 <<< "$m")" "$r" | tee -a "$out"
  rm -f "$log"
done
