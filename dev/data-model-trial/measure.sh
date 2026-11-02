#!/usr/bin/env bash
# Run a command in a fresh process and print its wall time (s) and peak
# resident memory (MB) as one tab-separated line on stderr:
#   MEASURE <label> <seconds> <peak_MB> <exit status>
# Usage: measure.sh <label> <command> [args...]
# The command runs with a 12 GB address-space limit, so a step that tries to
# read a whole large file fails instead of pushing the machine into swap.
label="$1"; shift
ulimit -v $((12 * 1024 * 1024))
tmp="$(mktemp)"
/usr/bin/time -f "%e %M %x" -o "$tmp" "$@"
status=$?
read -r secs kb code < <(tail -n 1 "$tmp")
rm -f "$tmp"
printf 'MEASURE\t%s\t%s\t%d\t%s\n' "$label" "$secs" "$((kb / 1024))" "$code" >&2
exit "$status"
