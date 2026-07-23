#!/usr/bin/env bash

# Stream PLINK metadata once and emit only genome-build marker candidates.

set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "Usage: $0 <pgen|bed> <metadata> <marker-ids> <marker-coordinates>" >&2
  exit 2
fi

kind=$1
metadata=$2
marker_ids=$3
marker_coordinates=$4

if [[ "$kind" != "pgen" && "$kind" != "bed" ]]; then
  echo "ERROR: unsupported genotype type '$kind'. Use 'pgen' or 'bed'." >&2
  exit 2
fi

for path in "$metadata" "$marker_ids" "$marker_coordinates"; do
  if [[ ! -r "$path" ]]; then
    echo "ERROR: genome-build scanner input is not readable: $path" >&2
    exit 2
  fi
done

if [[ ! -s "$metadata" ]]; then
  echo "ERROR: genome-build metadata file is empty: $metadata" >&2
  exit 2
fi

if ! command -v gawk >/dev/null 2>&1; then
  echo "ERROR: GNU awk (gawk) is required for fast genome-build metadata scanning" >&2
  exit 127
fi

exec gawk \
  -v kind="$kind" \
  -v marker_ids_file="$marker_ids" \
  -v marker_coordinates_file="$marker_coordinates" \
  '
  function fail(message) {
    print "ERROR: " message > "/dev/stderr"
    failed = 1
    exit 2
  }

  function clean_chrom(value) {
    if (tolower(substr(value, 1, 3)) == "chr") {
      value = substr(value, 4)
    }
    return value
  }

  function normalized_pos(value) {
    if (value !~ /^[0-9]+$/) {
      return ""
    }
    return sprintf("%.0f", value + 0)
  }

  function current_epoch(command, value) {
    command = "date +%s"
    command | getline value
    close(command)
    return value + 0
  }

  function report_progress(complete, elapsed, rate, suffix) {
    elapsed = current_epoch() - started
    if (elapsed < 1) {
      elapsed = 1
    }
    rate = rows_scanned / elapsed
    suffix = complete ? " (complete)" : ""
    printf "Genome-build metadata scan: %d rows, %d candidates, %d s, %.0f rows/s%s\n", \
      rows_scanned, candidates, elapsed, rate, suffix > "/dev/stderr"
  }

  BEGIN {
    FS = kind == "pgen" ? "\t" : " "
    OFS = "\t"
    started = current_epoch()

    while ((getline marker_id < marker_ids_file) > 0) {
      sub(/\r$/, "", marker_id)
      if (marker_id != "") {
        marker_ids[marker_id] = 1
      }
    }
    close(marker_ids_file)

    while ((getline coordinate < marker_coordinates_file) > 0) {
      sub(/\r$/, "", coordinate)
      count = split(coordinate, fields, "\t")
      position = count == 2 ? normalized_pos(fields[2]) : ""
      if (count != 2 || position == "") {
        fail("invalid marker coordinate: " coordinate)
      }
      marker_coordinates[clean_chrom(fields[1]) SUBSEP position] = 1
    }
    close(marker_coordinates_file)

    print "chrom", "pos", "variant_id"
  }

  /^[[:space:]]*$/ {
    next
  }

  kind == "pgen" && /^##/ {
    next
  }

  kind == "pgen" && !pvar_header {
    for (i = 1; i <= NF; i++) {
      name = $i
      sub(/\r$/, "", name)
      sub(/^#/, "", name)
      if (name == "CHROM") {
        chrom_column = i
      } else if (name == "POS") {
        pos_column = i
      } else if (name == "ID") {
        id_column = i
      }
    }
    if (!chrom_column || !pos_column || !id_column) {
      fail("PVAR primary header is missing required CHROM, POS, or ID columns")
    }
    pvar_header = 1
    next
  }

  {
    if (kind == "bed") {
      if (NF < 6) {
        fail("BIM row " FNR " must contain at least 6 columns")
      }
      chrom_column = 1
      id_column = 2
      pos_column = 4
    } else if (NF < chrom_column || NF < pos_column || NF < id_column) {
      fail("PVAR row " FNR " has fewer columns than its primary header")
    }

    rows_scanned++
    chrom = clean_chrom($(chrom_column))
    position = normalized_pos($(pos_column))
    variant_id = $(id_column)
    sub(/\r$/, "", variant_id)

    if (position == "") {
      fail(toupper(kind) " metadata row " FNR " has a nonnumeric position")
    }

    coordinate_key = chrom SUBSEP position
    if ((variant_id in marker_ids) || (coordinate_key in marker_coordinates)) {
      print chrom, position, variant_id
      candidates++
    }

    if (rows_scanned % 5000000 == 0) {
      report_progress(0)
    }
  }

  END {
    if (failed) {
      exit 2
    }
    if (kind == "pgen" && !pvar_header) {
      print "ERROR: PVAR primary header was not found" > "/dev/stderr"
      exit 2
    }
    report_progress(1)
  }
  ' \
  "$metadata"
