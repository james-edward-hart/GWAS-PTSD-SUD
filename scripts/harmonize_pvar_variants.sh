#!/usr/bin/env bash

# Harmonize ancestry-only PVAR IDs with bounded memory.
#
# Sorted PVARs take a streaming path. Unsorted PVARs are normalized first and
# then externally sorted; scientific matching and integrity checks are shared.

set -euo pipefail

usage() {
  echo "Usage: $0 --reference-pvar PATH --study-pvar PATH --output-dir DIR --exclude-palindromic true|false" >&2
  exit 2
}

reference_pvar=""
study_pvar=""
output_dir=""
exclude_palindromic=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --reference-pvar)
      [[ $# -ge 2 ]] || usage
      reference_pvar=$2
      shift 2
      ;;
    --study-pvar)
      [[ $# -ge 2 ]] || usage
      study_pvar=$2
      shift 2
      ;;
    --output-dir)
      [[ $# -ge 2 ]] || usage
      output_dir=$2
      shift 2
      ;;
    --exclude-palindromic)
      [[ $# -ge 2 ]] || usage
      exclude_palindromic=$2
      shift 2
      ;;
    *)
      usage
      ;;
  esac
done

[[ -n "$reference_pvar" && -n "$study_pvar" && -n "$output_dir" ]] || usage
[[ "$exclude_palindromic" == "true" || "$exclude_palindromic" == "false" ]] || usage

for path in "$reference_pvar" "$study_pvar"; do
  if [[ ! -r "$path" ]]; then
    echo "ERROR: ancestry PVAR is not readable: $path" >&2
    exit 2
  fi
  if [[ ! -s "$path" ]]; then
    echo "ERROR: ancestry PVAR is empty: $path" >&2
    exit 2
  fi
done

awk_bin=${GAWK_BIN:-gawk}
if ! command -v "$awk_bin" >/dev/null 2>&1; then
  echo "ERROR: GNU awk (gawk) is required for fast variant-ID harmonization" >&2
  exit 127
fi
if ! command -v sort >/dev/null 2>&1; then
  echo "ERROR: sort is required for variant-ID integrity checks" >&2
  exit 127
fi

mkdir -p "$output_dir"
work_dir=$(mktemp -d "$output_dir/.variant-harmonization.XXXXXX")
cleanup() {
  if [[ -n "${work_dir:-}" && -d "$work_dir" ]]; then
    rm -rf -- "$work_dir"
  fi
}
trap cleanup EXIT

reference_records="$work_dir/reference.records.tsv"
study_records="$work_dir/study.records.tsv"
reference_ids="$work_dir/reference.native_ids.tsv"
study_ids="$work_dir/study.native_ids.tsv"
reference_scan="$work_dir/reference.scan.tsv"
study_scan="$work_dir/study.scan.tsv"

mapping="$work_dir/variant_harmonization.tsv"
mismatches="$work_dir/shared_variant_mismatches.tsv"
shared="$work_dir/shared_variants.txt"
reference_extract="$work_dir/reference_native_variants.txt"
study_extract="$work_dir/study_native_variants.txt"
reference_update="$work_dir/reference_update_names.tsv"
study_update="$work_dir/study_update_names.tsv"
collision_index="$work_dir/harmonized_ids.tsv"
merge_summary="$work_dir/merge_summary.tsv"
summary="$work_dir/harmonization_summary.tsv"


# Normalize one PVAR and record whether eligible rows are coordinate-sorted.
normalize_pvar() {
  local side=$1
  local input=$2
  local records=$3
  local native_ids=$4
  local scan_summary=$5
  local started
  started=$(date +%s)

  "$awk_bin" \
    -v side="$side" \
    -v native_id_file="$native_ids" \
    -v summary_file="$scan_summary" \
    -v started="$started" \
    '
    function fail(message) {
      print "ERROR: " side " ancestry PVAR " message > "/dev/stderr"
      failed = 1
      exit 2
    }

    function trim(value) {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
      return value
    }

    function clean_chrom(value) {
      value = trim(value)
      if (tolower(substr(value, 1, 3)) == "chr") {
        value = substr(value, 4)
      }
      return value
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
      rate = scanned / elapsed
      suffix = complete ? " (complete)" : ""
      printf "%s ancestry PVAR scan: %d rows, %d eligible, %d s, %.0f rows/s%s\n", \
        side, scanned, eligible, elapsed, rate, suffix > "/dev/stderr"
    }

    BEGIN {
      FS = "\t"
      OFS = "\t"
      ordered = 1
    }

    /^[[:space:]]*$/ {
      next
    }

    /^##/ {
      next
    }

    !header_found {
      for (i = 1; i <= NF; i++) {
        name = trim($i)
        sub(/^#/, "", name)
        name = toupper(name)
        if (name == "CHROM") chrom_column = i
        else if (name == "POS") pos_column = i
        else if (name == "ID") id_column = i
        else if (name == "REF") ref_column = i
        else if (name == "ALT") alt_column = i
        else if (name == "INFO") info_column = i
      }
      if (!chrom_column || !pos_column || !id_column || !ref_column || !alt_column) {
        fail("header is missing required CHROM, POS, ID, REF, or ALT columns")
      }
      maximum_column = chrom_column
      if (pos_column > maximum_column) maximum_column = pos_column
      if (id_column > maximum_column) maximum_column = id_column
      if (ref_column > maximum_column) maximum_column = ref_column
      if (alt_column > maximum_column) maximum_column = alt_column
      if (info_column > maximum_column) maximum_column = info_column
      header_found = 1
      next
    }

    {
      scanned++
      if (NF < maximum_column) {
        fail("row " FNR " has fewer columns than its primary header")
      }

      position_text = trim($(pos_column))
      if (position_text !~ /^[0-9]+$/ || position_text + 0 <= 0) {
        fail("row " FNR " has a non-integer or nonpositive POS value")
      }
      position = position_text + 0

      chrom_text = clean_chrom($(chrom_column))
      if (chrom_text !~ /^[0-9]+$/) next
      chrom_number = chrom_text + 0
      if (chrom_number < 1 || chrom_number > 22) next

      native_id = trim($(id_column))
      ref = toupper(trim($(ref_column)))
      alt = toupper(trim($(alt_column)))
      if (native_id == "" || native_id == ".") next
      if (ref !~ /^[ACGT]+$/ || alt !~ /^[ACGT]+$/ || ref == alt) next

      chrom = sprintf("%.0f", chrom_number)
      canonical_position = sprintf("%.0f", position)
      chrom_order = sprintf("%02.0f", chrom_number)
      pos_order = sprintf("%012.0f", position)
      allele_low = ref < alt ? ref : alt
      allele_high = ref < alt ? alt : ref
      info = info_column ? trim($(info_column)) : ""
      provisional = (";" info ";") ~ /;PR;/ ? "True" : "False"

      coordinate_key = chrom_order OFS pos_order
      if (eligible > 0 && coordinate_key < previous_coordinate) ordered = 0
      previous_coordinate = coordinate_key
      eligible++

      print chrom_order, pos_order, allele_low, allele_high, native_id, \
        chrom, canonical_position, ref, alt, provisional
      print native_id, chrom_order, pos_order, allele_low, allele_high, \
        ref, alt, provisional > native_id_file

      if (scanned % 5000000 == 0) report_progress(0)
    }

    END {
      if (failed) exit 2
      if (!header_found) {
        print "ERROR: " side " ancestry PVAR primary header was not found" > "/dev/stderr"
        exit 2
      }
      if (!eligible) {
        print "ERROR: " side " ancestry PVAR contains no eligible autosomal biallelic ACGT variants" > "/dev/stderr"
        exit 2
      }
      print side, scanned, eligible, scanned - eligible, \
        (ordered ? "stream" : "fallback") > summary_file
      close(native_id_file)
      close(summary_file)
      report_progress(1)
    }
    ' \
    "$input" > "$records"
}


# Use a full normalized-record sort only when the input coordinate order fails.
prepare_merge_stream() {
  local side=$1
  local records=$2
  local scan_summary=$3
  local mode
  mode=$(cut -f 5 "$scan_summary")
  if [[ "$mode" == "fallback" ]]; then
    local sorted_records="${records}.sorted"
    echo "$side ancestry PVAR is not coordinate-sorted; using external-sort fallback" >&2
    LC_ALL=C TMPDIR="$work_dir" sort "$records" > "$sorted_records"
    mv "$sorted_records" "$records"
  elif [[ "$mode" != "stream" ]]; then
    echo "ERROR: internal $side PVAR scan mode is invalid: $mode" >&2
    exit 2
  fi
}


# Native IDs must be unique so each PLINK extraction and rename is unambiguous.
check_native_ids() {
  local side=$1
  local index=$2
  local sorted_index="${index}.sorted"
  LC_ALL=C TMPDIR="$work_dir" sort "$index" > "$sorted_index"
  "$awk_bin" -F $'\t' -v side="$side" '
    NR == 1 {
      previous_id = $1
      previous_line = $0
      next
    }
    $1 == previous_id {
      printf "ERROR: %s ancestry PVAR contains duplicate native variant ID '\''%s'\''. Records: %s | %s\n", \
        side, $1, previous_line, $0 > "/dev/stderr"
      exit 2
    }
    {
      previous_id = $1
      previous_line = $0
    }
  ' "$sorted_index"
  rm -f "$sorted_index"
}


normalize_pvar "reference" "$reference_pvar" "$reference_records" "$reference_ids" "$reference_scan"
normalize_pvar "study" "$study_pvar" "$study_records" "$study_ids" "$study_scan"
prepare_merge_stream "reference" "$reference_records" "$reference_scan"
prepare_merge_stream "study" "$study_records" "$study_scan"
check_native_ids "reference" "$reference_ids"
check_native_ids "study" "$study_ids"


# Merge normalized locus streams. Arrays are cleared after every locus, so
# memory is bounded by the explicit per-locus record ceiling.
"$awk_bin" \
  -v reference_path="$reference_records" \
  -v study_path="$study_records" \
  -v mapping_output="$mapping" \
  -v mismatch_output="$mismatches" \
  -v shared_output="$shared" \
  -v reference_extract_output="$reference_extract" \
  -v study_extract_output="$study_extract" \
  -v reference_update_output="$reference_update" \
  -v study_update_output="$study_update" \
  -v collision_output="$collision_index" \
  -v summary_output="$merge_summary" \
  -v exclude_palindromic="$exclude_palindromic" \
  '
  function fail(message) {
    print "ERROR: " message > "/dev/stderr"
    failed = 1
    exit 2
  }

  function clear_array(values, key) {
    for (key in values) delete values[key]
  }

  function load_group(path, state, rows, counts, first,    line, status, fields,
                      field_count, locus, next_locus, key, row_count) {
    clear_array(rows)
    clear_array(counts)
    clear_array(first)

    if (state["eof"]) return 0
    if (state["has_pending"]) {
      line = state["pending"]
      state["has_pending"] = 0
    } else {
      status = getline line < path
      if (status < 0) fail("could not read normalized PVAR stream: " path)
      if (status == 0) {
        state["eof"] = 1
        close(path)
        return 0
      }
    }

    field_count = split(line, fields, "\t")
    if (field_count != 10) fail("internal normalized PVAR record is malformed")
    locus = fields[1] "\t" fields[2]

    while (1) {
      row_count++
      if (row_count > 10000) {
        fail("ancestry PVAR has more than 10000 eligible records at locus " \
          fields[6] ":" fields[7])
      }

      key = fields[3] "/" fields[4]
      rows[row_count, "allele_low"] = fields[3]
      rows[row_count, "allele_high"] = fields[4]
      rows[row_count, "native_id"] = fields[5]
      rows[row_count, "chrom"] = fields[6]
      rows[row_count, "pos"] = fields[7]
      rows[row_count, "ref"] = fields[8]
      rows[row_count, "alt"] = fields[9]
      rows[row_count, "ref_provisional"] = fields[10]
      rows[row_count, "key"] = key
      rows[row_count, "line"] = line
      counts[key]++
      if (!(key in first)) first[key] = row_count

      status = getline line < path
      if (status < 0) fail("could not read normalized PVAR stream: " path)
      if (status == 0) {
        state["eof"] = 1
        close(path)
        break
      }
      field_count = split(line, fields, "\t")
      if (field_count != 10) fail("internal normalized PVAR record is malformed")
      next_locus = fields[1] "\t" fields[2]
      if (next_locus != locus) {
        state["pending"] = line
        state["has_pending"] = 1
        break
      }
    }

    state["key"] = locus
    return row_count
  }

  function collect_keys(counts, keys) {
    clear_array(keys)
    return asorti(counts, keys)
  }

  function collect_union_keys(reference_counts, study_counts, keys,
                              key) {
    clear_array(keys)
    clear_array(union_seen)
    for (key in reference_counts) union_seen[key] = 1
    for (key in study_counts) union_seen[key] = 1
    return asorti(union_seen, keys)
  }

  function row_selected(rows, row_index, mode, filter_key) {
    if (mode == "none") return 0
    if (mode == "all") return 1
    if (mode == "key") return rows[row_index, "key"] == filter_key
    if (mode == "valid") return rows[row_index, "valid"]
    if (mode == "unmatched") {
      return rows[row_index, "valid"] && !rows[row_index, "common"]
    }
    fail("internal mismatch row filter is invalid: " mode)
  }

  # Audit values follow normalized-record order so fast and fallback paths
  # produce identical mismatch reports.
  function collapse_field(rows, row_count, field, mode, filter_key,
                          i, selected_count, part_count, row_index, value,
                          result) {
    clear_array(collapse_order)
    clear_array(collapse_seen)
    for (i = 1; i <= row_count; i++) {
      if (row_selected(rows, i, mode, filter_key)) {
        collapse_order[++selected_count] = rows[i, "line"] SUBSEP i
      }
    }
    if (selected_count > 1) asort(collapse_order)

    for (i = 1; i <= selected_count; i++) {
      clear_array(collapse_parts)
      part_count = split(collapse_order[i], collapse_parts, SUBSEP)
      row_index = collapse_parts[part_count]
      value = rows[row_index, field]
      if (!(value in collapse_seen)) {
        collapse_seen[value] = 1
        result = result == "" ? value : result ";" value
      }
    }
    return result
  }

  function write_mismatch(reason, excluded_side,
                          reference_rows, reference_count, reference_mode, reference_key,
                          study_rows, study_count, study_mode, study_key) {
    print reason, excluded_side, \
      collapse_field(reference_rows, reference_count, "native_id", reference_mode, reference_key), \
      collapse_field(study_rows, study_count, "native_id", study_mode, study_key), \
      collapse_field(reference_rows, reference_count, "chrom", reference_mode, reference_key), \
      collapse_field(study_rows, study_count, "chrom", study_mode, study_key), \
      collapse_field(reference_rows, reference_count, "pos", reference_mode, reference_key), \
      collapse_field(study_rows, study_count, "pos", study_mode, study_key), \
      collapse_field(reference_rows, reference_count, "ref", reference_mode, reference_key), \
      collapse_field(reference_rows, reference_count, "alt", reference_mode, reference_key), \
      collapse_field(study_rows, study_count, "ref", study_mode, study_key), \
      collapse_field(study_rows, study_count, "alt", study_mode, study_key), \
      collapse_field(reference_rows, reference_count, "ref_provisional", reference_mode, reference_key), \
      collapse_field(study_rows, study_count, "ref_provisional", study_mode, study_key) \
      > mismatch_output
    mismatch_count++
  }

  function report_unshared_duplicates(side, row_count, counts,
                                      keys, key_count, i, key, reason) {
    key_count = collect_keys(counts, keys)
    for (i = 1; i <= key_count; i++) {
      key = keys[i]
      if (counts[key] <= 1) continue
      reason = "duplicate_locus_allele_key_" side
      if (side == "reference") {
        write_mismatch(reason, side,
          reference_rows, row_count, "key", key,
          study_rows, 0, "none", "")
      } else {
        write_mismatch(reason, side,
          reference_rows, 0, "none", "",
          study_rows, row_count, "key", key)
      }
    }
  }

  function retain_pair(reference_index, study_index, selected_id, selected_type,
                       selected_source, differing_rsid, orientation) {
    orientation = reference_rows[reference_index, "ref"] == study_rows[study_index, "ref"] && \
      reference_rows[reference_index, "alt"] == study_rows[study_index, "alt"] ? \
      "same" : "ref_alt_swap"

    print selected_id, selected_type, selected_source, \
      "locus_and_unordered_alleles", \
      reference_rows[reference_index, "native_id"], \
      study_rows[study_index, "native_id"], \
      reference_rows[reference_index, "chrom"], \
      reference_rows[reference_index, "pos"], \
      reference_rows[reference_index, "ref"], \
      reference_rows[reference_index, "alt"], \
      study_rows[study_index, "ref"], \
      study_rows[study_index, "alt"], \
      orientation, \
      reference_rows[reference_index, "ref_provisional"], \
      study_rows[study_index, "ref_provisional"], \
      differing_rsid > mapping_output

    print selected_id > shared_output
    print reference_rows[reference_index, "native_id"] > reference_extract_output
    print study_rows[study_index, "native_id"] > study_extract_output
    print reference_rows[reference_index, "native_id"], selected_id > reference_update_output
    print study_rows[study_index, "native_id"], selected_id > study_update_output
    print selected_id, reference_rows[reference_index, "chrom"], \
      reference_rows[reference_index, "pos"] > collision_output
    retained_count++
  }

  function process_shared_locus(reference_count, study_count,
                                keys, key_count, i, key, reference_index, study_index,
                                duplicate_reference, duplicate_study, reason,
                                excluded_side, reference_rsid, study_rsid,
                                selected_id, selected_type, selected_source,
                                differing_rsid, palindromic,
                                reference_valid, study_valid,
                                reference_unmatched, study_unmatched,
                                reference_mode, study_mode) {
    clear_array(duplicate_key)
    key_count = collect_union_keys(reference_key_counts, study_key_counts, keys)

    for (i = 1; i <= key_count; i++) {
      key = keys[i]
      duplicate_reference = reference_key_counts[key] > 1
      duplicate_study = study_key_counts[key] > 1
      if (!duplicate_reference && !duplicate_study) continue

      duplicate_key[key] = 1
      if (duplicate_reference && duplicate_study) {
        reason = "duplicate_locus_allele_key_both"
      } else if (duplicate_reference) {
        reason = "duplicate_locus_allele_key_reference"
      } else {
        reason = "duplicate_locus_allele_key_study"
      }
      excluded_side = reference_key_counts[key] > 0 && study_key_counts[key] > 0 ? \
        "both" : (reference_key_counts[key] > 0 ? "reference" : "study")
      write_mismatch(reason, excluded_side,
        reference_rows, reference_count,
        reference_key_counts[key] > 0 ? "key" : "none", key,
        study_rows, study_count,
        study_key_counts[key] > 0 ? "key" : "none", key)
    }

    for (i = 1; i <= reference_count; i++) {
      reference_rows[i, "valid"] = !(reference_rows[i, "key"] in duplicate_key)
      reference_rows[i, "common"] = 0
      if (reference_rows[i, "valid"]) reference_valid++
    }
    for (i = 1; i <= study_count; i++) {
      study_rows[i, "valid"] = !(study_rows[i, "key"] in duplicate_key)
      study_rows[i, "common"] = 0
      if (study_rows[i, "valid"]) study_valid++
    }

    for (i = 1; i <= key_count; i++) {
      key = keys[i]
      if (key in duplicate_key) continue
      if (reference_key_counts[key] != 1 || study_key_counts[key] != 1) continue

      reference_index = reference_first[key]
      study_index = study_first[key]
      reference_rows[reference_index, "common"] = 1
      study_rows[study_index, "common"] = 1

      palindromic = length(reference_rows[reference_index, "allele_low"]) == 1 && \
        (reference_rows[reference_index, "allele_low"] \
          reference_rows[reference_index, "allele_high"] == "AT" || \
        reference_rows[reference_index, "allele_low"] \
          reference_rows[reference_index, "allele_high"] == "CG")
      if (exclude_palindromic == "true" && palindromic) {
        write_mismatch("palindromic_snp_excluded", "both",
          reference_rows, reference_count, "key", key,
          study_rows, study_count, "key", key)
        continue
      }

      if (reference_rows[reference_index, "ref_provisional"] == "False" && \
          study_rows[study_index, "ref_provisional"] == "False" && \
          reference_rows[reference_index, "ref"] != study_rows[study_index, "ref"]) {
        write_mismatch("pvar_known_ref_disagreement", "both",
          reference_rows, reference_count, "key", key,
          study_rows, study_count, "key", key)
        continue
      }

      reference_rsid = reference_rows[reference_index, "native_id"] ~ /^rs[0-9]+$/
      study_rsid = study_rows[study_index, "native_id"] ~ /^rs[0-9]+$/
      selected_id = ""
      selected_type = ""
      selected_source = ""
      differing_rsid = "False"

      if (reference_rsid) {
        selected_id = reference_rows[reference_index, "native_id"]
        selected_type = "rsID"
        selected_source = "reference_rsid"
        if (study_rsid && reference_rows[reference_index, "native_id"] != \
            study_rows[study_index, "native_id"]) differing_rsid = "True"
      } else if (study_rsid) {
        selected_id = study_rows[study_index, "native_id"]
        selected_type = "rsID"
        selected_source = "study_rsid"
      } else if (reference_rows[reference_index, "ref_provisional"] == "False") {
        selected_id = reference_rows[reference_index, "chrom"] ":" \
          reference_rows[reference_index, "pos"] ":" \
          reference_rows[reference_index, "ref"] ":" \
          reference_rows[reference_index, "alt"]
        selected_type = "CPRA"
        selected_source = "reference_pvar_known_ref"
      } else if (study_rows[study_index, "ref_provisional"] == "False") {
        selected_id = study_rows[study_index, "chrom"] ":" \
          study_rows[study_index, "pos"] ":" \
          study_rows[study_index, "ref"] ":" \
          study_rows[study_index, "alt"]
        selected_type = "CPRA"
        selected_source = "study_pvar_known_ref"
      } else {
        write_mismatch("no_pvar_known_ref_for_cpra", "both",
          reference_rows, reference_count, "key", key,
          study_rows, study_count, "key", key)
        continue
      }

      retain_pair(reference_index, study_index, selected_id, selected_type, \
        selected_source, differing_rsid)
    }

    for (i = 1; i <= reference_count; i++) {
      if (reference_rows[i, "valid"] && !reference_rows[i, "common"]) {
        reference_unmatched++
      }
    }
    for (i = 1; i <= study_count; i++) {
      if (study_rows[i, "valid"] && !study_rows[i, "common"]) study_unmatched++
    }

    if (reference_unmatched || study_unmatched) {
      reference_mode = reference_unmatched ? "unmatched" : "valid"
      study_mode = study_unmatched ? "unmatched" : "valid"
      if (reference_valid && study_valid) {
        excluded_side = reference_unmatched && study_unmatched ? \
          "both" : (reference_unmatched ? "reference" : "study")
        write_mismatch("allele_mismatch", excluded_side,
          reference_rows, reference_count, reference_mode, "",
          study_rows, study_count, study_mode, "")
      }
    }
  }

  BEGIN {
    OFS = "\t"
    print "harmonized_id", "harmonized_id_type", "harmonized_id_source", \
      "match_method", "reference_native_id", "study_native_id", \
      "chrom", "pos", "reference_ref", "reference_alt", "study_ref", "study_alt", \
      "orientation", "reference_ref_provisional", "study_ref_provisional", \
      "differing_rsid" > mapping_output
    print "reason", "excluded_side", \
      "reference_native_id", "study_native_id", \
      "reference_chrom", "study_chrom", \
      "reference_pos", "study_pos", \
      "reference_ref", "reference_alt", "study_ref", "study_alt", \
      "reference_ref_provisional", "study_ref_provisional" > mismatch_output

    # Ensure headerless outputs exist even when no pair is retained.
    printf "%s", "" > shared_output
    printf "%s", "" > reference_extract_output
    printf "%s", "" > study_extract_output
    printf "%s", "" > reference_update_output
    printf "%s", "" > study_update_output
    printf "%s", "" > collision_output
    close(shared_output)
    close(reference_extract_output)
    close(study_extract_output)
    close(reference_update_output)
    close(study_update_output)
    close(collision_output)

    reference_count = load_group(reference_path, reference_state, reference_rows, \
      reference_key_counts, reference_first)
    study_count = load_group(study_path, study_state, study_rows, \
      study_key_counts, study_first)

    while (reference_count > 0 && study_count > 0) {
      if (reference_state["key"] < study_state["key"]) {
        report_unshared_duplicates("reference", reference_count, \
          reference_key_counts)
        reference_count = load_group(reference_path, reference_state, \
          reference_rows, reference_key_counts, reference_first)
      } else if (study_state["key"] < reference_state["key"]) {
        report_unshared_duplicates("study", study_count, \
          study_key_counts)
        study_count = load_group(study_path, study_state, study_rows, \
          study_key_counts, study_first)
      } else {
        process_shared_locus(reference_count, study_count)
        reference_count = load_group(reference_path, reference_state, \
          reference_rows, reference_key_counts, reference_first)
        study_count = load_group(study_path, study_state, study_rows, \
          study_key_counts, study_first)
      }
    }

    while (reference_count > 0) {
      report_unshared_duplicates("reference", reference_count, \
        reference_key_counts)
      reference_count = load_group(reference_path, reference_state, \
        reference_rows, reference_key_counts, reference_first)
    }
    while (study_count > 0) {
      report_unshared_duplicates("study", study_count, \
        study_key_counts)
      study_count = load_group(study_path, study_state, study_rows, \
        study_key_counts, study_first)
    }

    print retained_count + 0, mismatch_count + 0 > summary_output
    close(mapping_output)
    close(mismatch_output)
    close(shared_output)
    close(reference_extract_output)
    close(study_extract_output)
    close(reference_update_output)
    close(study_update_output)
    close(collision_output)
    close(summary_output)
  }
  ' </dev/null


# A selected operational ID must identify exactly one retained variant.
sorted_collisions="${collision_index}.sorted"
LC_ALL=C TMPDIR="$work_dir" sort "$collision_index" > "$sorted_collisions"
"$awk_bin" -F $'\t' '
  NR == 1 {
    previous_id = $1
    previous_locus = $2 ":" $3
    next
  }
  $1 == previous_id {
    locus = $2 ":" $3
    if (locus != previous_locus) {
      printf "ERROR: harmonized variant ID '\''%s'\'' maps to multiple loci: %s, %s\n", \
        $1, previous_locus, locus > "/dev/stderr"
    } else {
      printf "ERROR: harmonized variant ID '\''%s'\'' maps to multiple variants at locus %s\n", \
        $1, locus > "/dev/stderr"
    }
    exit 2
  }
  {
    previous_id = $1
    previous_locus = $2 ":" $3
  }
' "$sorted_collisions"
rm -f "$sorted_collisions"

IFS=$'\t' read -r _ reference_scanned reference_eligible reference_excluded reference_mode < "$reference_scan"
IFS=$'\t' read -r _ study_scanned study_eligible study_excluded study_mode < "$study_scan"
IFS=$'\t' read -r retained mismatched < "$merge_summary"

{
  printf "retained\tmismatches\treference_scanned\treference_eligible\tstudy_scanned\tstudy_eligible\treference_sort_mode\tstudy_sort_mode\n"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$retained" "$mismatched" \
    "$reference_scanned" "$reference_eligible" \
    "$study_scanned" "$study_eligible" \
    "$reference_mode" "$study_mode"
} > "$summary"

for artifact in \
  variant_harmonization.tsv \
  shared_variant_mismatches.tsv \
  shared_variants.txt \
  reference_native_variants.txt \
  study_native_variants.txt \
  reference_update_names.tsv \
  study_update_names.tsv \
  harmonization_summary.tsv
do
  mv -f "$work_dir/$artifact" "$output_dir/$artifact"
done

echo "Variant-ID harmonization: $retained retained, $mismatched exclusions; reference=$reference_mode, study=$study_mode" >&2
