#!/usr/bin/env bash

# Stream a PVAR and emit variant IDs inside configured long-range/problem
# regions. Only the small region table is held in memory.

set -euo pipefail
export LC_ALL=C

usage() {
  echo "Usage: $0 --pvar PATH --regions PATH --out PATH" >&2
  exit 2
}

pvar=""
regions=""
out=""

while [[ $# -gt 0 ]]; do
  [[ $# -ge 2 ]] || usage
  case "$1" in
    --pvar) pvar=$2 ;;
    --regions) regions=$2 ;;
    --out) out=$2 ;;
    *) usage ;;
  esac
  shift 2
done

[[ -n "$pvar" && -n "$regions" && -n "$out" ]] || usage
for path in "$pvar" "$regions"; do
  [[ -r "$path" && -s "$path" ]] || {
    echo "ERROR: Phase 2 region-exclusion input is not readable or is empty: $path" >&2
    exit 2
  }
done

awk_bin=${GAWK_BIN:-gawk}
command -v "$awk_bin" >/dev/null 2>&1 || {
  echo "ERROR: GNU awk (gawk) is required for streaming Phase 2 region exclusions" >&2
  exit 127
}

mkdir -p "$(dirname "$out")"
unsorted_temp=$(mktemp "$(dirname "$out")/.$(basename "$out").unsorted.XXXXXX")
out_temp=$(mktemp "$(dirname "$out")/.$(basename "$out").XXXXXX")

cleanup() {
  [[ -z "${unsorted_temp:-}" ]] || rm -f -- "$unsorted_temp"
  [[ -z "${out_temp:-}" ]] || rm -f -- "$out_temp"
}
trap cleanup EXIT

"$awk_bin" \
    -v regions_file="$regions" \
    '
    function fail(message) {
      print "ERROR: " message > "/dev/stderr"
      failed = 1
      exit 2
    }

    function trim(value) {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
      return value
    }

    function normalized_chrom(value) {
      value = toupper(trim(value))
      if (substr(value, 1, 3) == "CHR") value = substr(value, 4)
      return value
    }

    function positive_integer(value, number) {
      value = trim(value)
      if (value !~ /^[0-9]+$/ || value + 0 <= 0) return ""
      number = value + 0
      return sprintf("%.0f", number)
    }

    function current_epoch(command, value) {
      command = "date +%s"
      command | getline value
      close(command)
      return value + 0
    }

    function report_progress(complete, elapsed, rate, suffix) {
      elapsed = current_epoch() - started
      if (elapsed < 1) elapsed = 1
      rate = rows_scanned / elapsed
      suffix = complete ? " (complete)" : ""
      printf "Phase 2 region PVAR scan: %d rows, %d excluded, %d s, %.0f rows/s%s\n", \
        rows_scanned, excluded, elapsed, rate, suffix > "/dev/stderr"
    }

    function load_regions(line, fields, field_count, i, name, required, chrom, start, end) {
      while ((getline line < regions_file) > 0) {
        sub(/\r$/, "", line)
        if (line ~ /^[[:space:]]*$/) continue
        field_count = split(line, fields, "\t")
        if (!region_header) {
          for (i = 1; i <= field_count; i++) {
            name = toupper(fields[i])
            sub(/^#/, "", name)
            if (name == "CHROM") region_chrom_column = i
            else if (name == "START") region_start_column = i
            else if (name == "END") region_end_column = i
          }
          if (!region_chrom_column || !region_start_column || !region_end_column) {
            fail("Phase 2 exclusion regions are missing CHROM, START, or END: " regions_file)
          }
          required = region_chrom_column
          if (region_start_column > required) required = region_start_column
          if (region_end_column > required) required = region_end_column
          region_required_columns = required
          region_header = 1
          continue
        }
        if (field_count < region_required_columns) {
          fail("Phase 2 exclusion-region row has fewer columns than its header: " regions_file)
        }
        chrom = normalized_chrom(fields[region_chrom_column])
        start = positive_integer(fields[region_start_column])
        end = positive_integer(fields[region_end_column])
        if (start == "" || end == "" || start + 0 > end + 0) {
          fail("Phase 2 exclusion region has invalid START/END coordinates: " line)
        }
        if (chrom ~ /^([1-9]|1[0-9]|2[0-2])$/) {
          region_count[chrom]++
          region_start[chrom, region_count[chrom]] = start + 0
          region_end[chrom, region_count[chrom]] = end + 0
        }
      }
      close(regions_file)
      if (!region_header) fail("Phase 2 exclusion-region header was not found: " regions_file)
    }

    BEGIN {
      FS = OFS = "\t"
      started = current_epoch()
      load_regions()
    }

    /^[[:space:]]*$/ || /^##/ { next }

    !pvar_header {
      for (i = 1; i <= NF; i++) {
        name = toupper($i)
        sub(/\r$/, "", name)
        sub(/^#/, "", name)
        if (name == "CHROM") chrom_column = i
        else if (name == "POS") pos_column = i
        else if (name == "ID") id_column = i
      }
      if (!chrom_column || !pos_column || !id_column) {
        fail("Phase 2 marker PVAR header is missing CHROM, POS, or ID: " FILENAME)
      }
      required_columns = chrom_column
      if (pos_column > required_columns) required_columns = pos_column
      if (id_column > required_columns) required_columns = id_column
      pvar_header = 1
      next
    }

    {
      if (NF < required_columns) {
        fail("Phase 2 marker PVAR row " FNR " has fewer columns than its header: " FILENAME)
      }
      rows_scanned++
      chrom = normalized_chrom($(chrom_column))
      position = positive_integer($(pos_column))
      if (position == "") {
        fail("Phase 2 marker PVAR row " FNR " has a non-integer or nonpositive position")
      }
      for (i = 1; i <= region_count[chrom]; i++) {
        if (position + 0 >= region_start[chrom, i] && position + 0 <= region_end[chrom, i]) {
          variant_id = $(id_column)
          sub(/\r$/, "", variant_id)
          print variant_id
          excluded++
          break
        }
      }
      if (rows_scanned % 5000000 == 0) report_progress(0)
    }

    END {
      if (failed) exit 2
      if (!pvar_header) {
        print "ERROR: Phase 2 marker PVAR header was not found: " FILENAME > "/dev/stderr"
        exit 2
      }
      report_progress(1)
    }
    ' \
    "$pvar" > "$unsorted_temp"

# Preserve the prior deterministic, unique ID artifact without using R memory.
sort -u "$unsorted_temp" > "$out_temp"

mv -f -- "$out_temp" "$out"
out_temp=""
