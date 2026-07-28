#!/usr/bin/env bash

# Inspect BIM/PVAR metadata with bounded memory. Clean inputs need one pass;
# only duplicate allele codes trigger a second pass for a PLINK-safe copy.

set -euo pipefail

usage() {
  echo "Usage: $0 --type bed|pgen --metadata PATH --summary PATH [--label TEXT] [--target-chromosomes CSV --stop-on-target] [--inspect-alleles] [--safe-metadata PATH --invalid-report PATH --invalid-exclude PATH]" >&2
  exit 2
}

kind=""
metadata=""
summary=""
label="genotype input"
targets=""
stop_on_target="false"
inspect_alleles="false"
safe_metadata=""
invalid_report=""
invalid_exclude=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --type|--metadata|--summary|--label|--target-chromosomes|--safe-metadata|--invalid-report|--invalid-exclude)
      [[ $# -ge 2 ]] || usage
      case "$1" in
        --type) kind=$2 ;;
        --metadata) metadata=$2 ;;
        --summary) summary=$2 ;;
        --label) label=$2 ;;
        --target-chromosomes) targets=$2 ;;
        --safe-metadata) safe_metadata=$2 ;;
        --invalid-report) invalid_report=$2 ;;
        --invalid-exclude) invalid_exclude=$2 ;;
      esac
      shift 2
      ;;
    --stop-on-target)
      stop_on_target="true"
      shift
      ;;
    --inspect-alleles)
      inspect_alleles="true"
      shift
      ;;
    *)
      usage
      ;;
  esac
done

[[ "$kind" == "bed" || "$kind" == "pgen" ]] || usage
[[ -n "$metadata" && -n "$summary" ]] || usage
format_label="BIM"
[[ "$kind" == "pgen" ]] && format_label="PVAR"
if [[ ! -r "$metadata" ]]; then
  echo "ERROR: $label metadata file is not readable: $metadata" >&2
  exit 2
fi
if [[ ! -s "$metadata" ]]; then
  echo "ERROR: $label $format_label file is empty: $metadata" >&2
  exit 2
fi

output_count=0
[[ -n "$safe_metadata" ]] && output_count=$((output_count + 1))
[[ -n "$invalid_report" ]] && output_count=$((output_count + 1))
[[ -n "$invalid_exclude" ]] && output_count=$((output_count + 1))
[[ "$output_count" -eq 0 || "$output_count" -eq 3 ]] || usage
if [[ "$output_count" -eq 3 ]]; then
  inspect_alleles="true"
fi

awk_bin=${GAWK_BIN:-gawk}
if ! command -v "$awk_bin" >/dev/null 2>&1; then
  echo "ERROR: GNU awk (gawk) is required for fast PLINK metadata inspection" >&2
  exit 127
fi

mkdir -p "$(dirname "$summary")"
if [[ "$output_count" -eq 3 ]]; then
  mkdir -p "$(dirname "$safe_metadata")" "$(dirname "$invalid_report")" "$(dirname "$invalid_exclude")"
fi

work_dir=$(mktemp -d "$(dirname "$summary")/.plink-metadata.XXXXXX")
summary_temp=$(mktemp "$(dirname "$summary")/.$(basename "$summary").XXXXXX")
invalid_index="$work_dir/invalid_rows.tsv"
if [[ "$output_count" -eq 3 ]]; then
  report_temp=$(mktemp "$(dirname "$invalid_report")/.$(basename "$invalid_report").XXXXXX")
  exclude_temp=$(mktemp "$(dirname "$invalid_exclude")/.$(basename "$invalid_exclude").XXXXXX")
else
  report_temp="$work_dir/invalid_report.tsv"
  exclude_temp="$work_dir/invalid_exclude.txt"
fi
safe_temp=""

cleanup() {
  if [[ -n "${safe_temp:-}" ]]; then
    rm -f -- "$safe_temp"
  fi
  rm -f -- "$summary_temp" "$report_temp" "$exclude_temp"
  if [[ -n "${work_dir:-}" && -d "$work_dir" ]]; then
    rm -rf -- "$work_dir"
  fi
}
trap cleanup EXIT

"$awk_bin" \
  -v kind="$kind" \
  -v label="$label" \
  -v targets_csv="$targets" \
  -v stop_on_target="$stop_on_target" \
  -v inspect_alleles="$inspect_alleles" \
  -v summary_file="$summary_temp" \
  -v invalid_index_file="$invalid_index" \
  -v report_file="$report_temp" \
  -v exclude_file="$exclude_temp" \
  '
  function fail(message) {
    print "ERROR: " label " " message > "/dev/stderr"
    failed = 1
    exit 2
  }

  function trim(value) {
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
    return value
  }

  function normalized_chrom(value) {
    value = toupper(trim(value))
    if (substr(value, 1, 3) == "CHR") {
      value = substr(value, 4)
    }
    if (value == "23") return "X"
    if (value == "24") return "Y"
    if (value == "25") return "XY"
    return value
  }

  function normalized_pos(value, number) {
    value = trim(value)
    if (value !~ /^[0-9]+$/ || value + 0 <= 0) {
      return ""
    }
    number = value + 0
    return sprintf("%.0f", number)
  }

  function split_row(raw, fields, line) {
    line = raw
    sub(/\r$/, "", line)
    if (kind == "pgen") {
      return split(line, fields, "\t")
    }
    line = trim(line)
    return split(line, fields, /[[:space:]]+/)
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
    suffix = stopped_early ? " (target found)" : (complete ? " (complete)" : "")
    printf "PLINK metadata scan: %d rows, %d invalid alleles, %d s, %.0f rows/s%s\n", \
      rows_scanned, invalid_count, elapsed, rate, suffix > "/dev/stderr"
  }

  function write_summary() {
    print "metric", "value" > summary_file
    print "metadata_type", kind >> summary_file
    print "rows_scanned", rows_scanned + 0 >> summary_file
    print "scan_complete", stopped_early ? "false" : "true" >> summary_file
    print "x_variants", x_count + 0 >> summary_file
    print "y_variants", y_count + 0 >> summary_file
    print "xy_variants", xy_count + 0 >> summary_file
    print "par1_variants", par1_count + 0 >> summary_file
    print "par2_variants", par2_count + 0 >> summary_file
    print "target_matches", target_matches + 0 >> summary_file
    print "invalid_duplicate_alleles", invalid_count + 0 >> summary_file
    print "pvar_info_column", info_column ? "true" : "false" >> summary_file
    close(summary_file)
  }

  BEGIN {
    OFS = "\t"
    started = current_epoch()
    report_header = kind == "bed" \
      ? "row_number\tvariant_id\tchrom\tpos\tallele1\tallele2\treplacement_variant_id\treason" \
      : "row_number\tvariant_id\tchrom\tpos\tref\talt\treplacement_variant_id\treason"
    print report_header > report_file
    close(report_file)
    printf "" > exclude_file
    close(exclude_file)
    printf "" > invalid_index_file
    close(invalid_index_file)

    count = split(targets_csv, requested, ",")
    for (i = 1; i <= count; i++) {
      target = normalized_chrom(requested[i])
      if (target != "") targets[target] = 1
    }
  }

  /^[[:space:]]*$/ { next }

  kind == "pgen" && /^##/ { next }

  {
    field_count = split_row($0, fields)

    if (kind == "pgen" && !pvar_header) {
      for (i = 1; i <= field_count; i++) {
        column_name = fields[i]
        if (column_name == "#CHROM" || column_name == "CHROM") chrom_column = i
        else if (column_name == "POS") pos_column = i
        else if (column_name == "ID") id_column = i
        else if (column_name == "REF") allele1_column = i
        else if (column_name == "ALT") allele2_column = i
        else if (column_name == "INFO") info_column = i
      }
      missing = ""
      if (!chrom_column) missing = missing " chrom"
      if (!id_column) missing = missing " variant_id"
      if (!pos_column) missing = missing " pos"
      if (inspect_alleles == "true" && !allele1_column) missing = missing " allele1"
      if (inspect_alleles == "true" && !allele2_column) missing = missing " allele2"
      if (missing != "") {
        gsub(/^ /, "", missing)
        gsub(/ /, ", ", missing)
        fail("PVAR header is missing required column(s): " missing ": " FILENAME)
      }
      pvar_header = 1
      next
    }

    if (kind == "bed") {
      if (field_count < 6) {
        fail("BIM row " FNR " must contain at least 6 columns: " FILENAME)
      }
      chrom_column = 1
      id_column = 2
      pos_column = 4
      allele1_column = 5
      allele2_column = 6
    } else {
      required_column = chrom_column
      if (pos_column > required_column) required_column = pos_column
      if (id_column > required_column) required_column = id_column
      if (inspect_alleles == "true" && allele1_column > required_column) required_column = allele1_column
      if (inspect_alleles == "true" && allele2_column > required_column) required_column = allele2_column
      if (field_count < required_column) {
        fail("PVAR metadata row " FNR " has fewer columns than its header: " FILENAME)
      }
    }

    rows_scanned++
    chrom = normalized_chrom(fields[chrom_column])
    position = normalized_pos(fields[pos_column])
    if (position == "") {
      fail(toupper(kind) " metadata row " FNR " has a non-integer or nonpositive position: " fields[pos_column])
    }

    if (chrom == "X") x_count++
    else if (chrom == "Y") y_count++
    else if (chrom == "XY") xy_count++
    else if (chrom == "PAR1") par1_count++
    else if (chrom == "PAR2") par2_count++

    if (chrom in targets) target_matches++

    if (inspect_alleles == "true") {
      allele1 = toupper(trim(fields[allele1_column]))
      allele2 = toupper(trim(fields[allele2_column]))
      invalid = allele1 != "" && allele2 != "" && allele1 == allele2
      if (kind == "pgen" && index(allele2, ",")) invalid = 0
      if (invalid) {
        invalid_count++
        replacement = kind == "bed" \
          ? "__stage1_excluded_invalid_bim_" rows_scanned \
          : "__stage1_excluded_invalid_pvar_" rows_scanned
        print rows_scanned, replacement >> invalid_index_file
        print replacement >> exclude_file
        print rows_scanned, fields[id_column], fields[chrom_column], position, \
          fields[allele1_column], fields[allele2_column], replacement, \
          "duplicate_allele_code" >> report_file
      }
    }

    if (rows_scanned % 5000000 == 0) report_progress(0)
    if (stop_on_target == "true" && target_matches > 0) {
      stopped_early = 1
      exit 0
    }
  }

  END {
    if (failed) exit 2
    if (kind == "pgen" && !pvar_header) {
      print "ERROR: " label " PVAR header not found: " FILENAME > "/dev/stderr"
      exit 2
    }
    write_summary()
    report_progress(1)
  }
  ' \
  "$metadata"

invalid_count=$("$awk_bin" -F $'\t' '$1 == "invalid_duplicate_alleles" { print $2 }' "$summary_temp")
invalid_count=${invalid_count:-0}

if [[ "$invalid_count" -gt 0 && "$output_count" -eq 3 ]]; then
  safe_temp=$(mktemp "$(dirname "$safe_metadata")/.$(basename "$safe_metadata").XXXXXX")

  # Invalid row numbers are naturally sorted, so the rewrite retains one row
  # from the index at a time instead of loading all exclusions into memory.
  "$awk_bin" \
    -v kind="$kind" \
    -v label="$label" \
    -v invalid_index_file="$invalid_index" \
    '
    function fail(message) {
      print "ERROR: " label " " message > "/dev/stderr"
      failed = 1
      exit 2
    }

    function trim(value) {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
      return value
    }

    function split_row(raw, fields, line) {
      line = raw
      sub(/\r$/, "", line)
      if (kind == "pgen") return split(line, fields, "\t")
      line = trim(line)
      return split(line, fields, /[[:space:]]+/)
    }

    function read_invalid(status, parts) {
      status = (getline invalid_line < invalid_index_file)
      if (status > 0) {
        split(invalid_line, parts, "\t")
        invalid_row = parts[1] + 0
        replacement = parts[2]
      } else {
        invalid_row = 0
        replacement = ""
        close(invalid_index_file)
      }
    }

    function print_fields(fields, count, i, line) {
      line = fields[1]
      for (i = 2; i <= count; i++) line = line OFS fields[i]
      print line
    }

    BEGIN {
      OFS = "\t"
      read_invalid()
    }

    /^[[:space:]]*$/ { print; next }

    kind == "pgen" && /^##/ { print; next }

    {
      field_count = split_row($0, fields)
      if (kind == "pgen" && !pvar_header) {
        for (i = 1; i <= field_count; i++) {
          if (fields[i] == "#CHROM" || fields[i] == "CHROM") chrom_column = i
          else if (fields[i] == "ID") id_column = i
          else if (fields[i] == "REF") allele1_column = i
          else if (fields[i] == "ALT") allele2_column = i
        }
        if (!chrom_column || !id_column || !allele1_column || !allele2_column) {
          fail("PVAR header changed during metadata sanitization: " FILENAME)
        }
        pvar_header = 1
        print
        next
      }
      if (kind == "bed") {
        id_column = 2
        allele1_column = 5
        allele2_column = 6
      }

      data_row++
      if (invalid_row && data_row > invalid_row) fail("duplicate-allele row index is inconsistent with " FILENAME)
      if (data_row == invalid_row) {
        fields[id_column] = replacement
        fields[allele1_column] = "A"
        fields[allele2_column] = "C"
        print_fields(fields, field_count)
        read_invalid()
      } else {
        print
      }
    }

    END {
      if (failed) exit 2
      if (invalid_row) {
        print "ERROR: " label " duplicate-allele row index extends beyond " FILENAME > "/dev/stderr"
        exit 2
      }
    }
    ' \
    "$metadata" > "$safe_temp"

  mv -f -- "$report_temp" "$invalid_report"
  mv -f -- "$exclude_temp" "$invalid_exclude"
  mv -f -- "$safe_temp" "$safe_metadata"
  safe_temp=""
elif [[ "$output_count" -eq 3 ]]; then
  # A successful clean scan invalidates stale sanitizer artifacts from a prior run.
  rm -f -- "$safe_metadata" "$invalid_report" "$invalid_exclude"
fi

mv -f -- "$summary_temp" "$summary"
summary_temp=""
