#!/usr/bin/env bash

# Stream a PVAR to apply the Phase 2 Step 1 INFO/R2 threshold without loading
# millions of variant rows into R. Missing quality values intentionally pass.

set -euo pipefail
export LC_ALL=C

usage() {
  echo "Usage: $0 --pvar PATH --threshold NUMBER --pass-ids PATH --excluded-report PATH --summary PATH" >&2
  exit 2
}

pvar=""
threshold=""
pass_ids=""
excluded_report=""
summary=""

while [[ $# -gt 0 ]]; do
  [[ $# -ge 2 ]] || usage
  case "$1" in
    --pvar) pvar=$2 ;;
    --threshold) threshold=$2 ;;
    --pass-ids) pass_ids=$2 ;;
    --excluded-report) excluded_report=$2 ;;
    --summary) summary=$2 ;;
    *) usage ;;
  esac
  shift 2
done

[[ -n "$pvar" && -n "$threshold" && -n "$pass_ids" && -n "$excluded_report" && -n "$summary" ]] || usage
[[ -r "$pvar" && -s "$pvar" ]] || {
  echo "ERROR: Phase 2 Step 1 PVAR is not readable or is empty: $pvar" >&2
  exit 2
}

awk_bin=${GAWK_BIN:-gawk}
command -v "$awk_bin" >/dev/null 2>&1 || {
  echo "ERROR: GNU awk (gawk) is required for streaming Phase 2 INFO/R2 filtering" >&2
  exit 127
}

for path in "$pass_ids" "$excluded_report" "$summary"; do
  mkdir -p "$(dirname "$path")"
done

pass_temp=$(mktemp "$(dirname "$pass_ids")/.$(basename "$pass_ids").XXXXXX")
report_temp=$(mktemp "$(dirname "$excluded_report")/.$(basename "$excluded_report").XXXXXX")
summary_temp=$(mktemp "$(dirname "$summary")/.$(basename "$summary").XXXXXX")

cleanup() {
  [[ -z "${pass_temp:-}" ]] || rm -f -- "$pass_temp"
  [[ -z "${report_temp:-}" ]] || rm -f -- "$report_temp"
  [[ -z "${summary_temp:-}" ]] || rm -f -- "$summary_temp"
}
trap cleanup EXIT

# The first pass identifies the highest-priority usable quality field. A second
# pass is needed only when that field actually excludes variants; this preserves
# the existing pass-list semantics while keeping memory bounded.
"$awk_bin" \
  -v threshold="$threshold" \
  -v pass_file="$pass_temp" \
  -v report_file="$report_temp" \
  -v summary_file="$summary_temp" \
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

  function finite_number(raw, rendered) {
    raw = trim(raw)
    if (raw !~ /^[+-]?(([0-9]+([.][0-9]*)?)|([.][0-9]+))([eE][+-]?[0-9]+)?$/) {
      return 0
    }
    parsed_number = raw + 0
    rendered = tolower(sprintf("%.17g", parsed_number))
    return rendered !~ /inf|nan/
  }

  # INFO key precedence matches the former R parser.
  function info_key_number(raw, pieces, piece_count, key_index, piece_index, equals_at, key, value) {
    raw = trim(raw)
    if (raw == "" || raw == "." || raw == "NA") return 0
    piece_count = split(raw, pieces, ";")
    for (key_index = 1; key_index <= key_count; key_index++) {
      for (piece_index = 1; piece_index <= piece_count; piece_index++) {
        equals_at = index(pieces[piece_index], "=")
        if (!equals_at) continue
        key = toupper(substr(pieces[piece_index], 1, equals_at - 1))
        if (key != info_keys[key_index]) continue
        value = substr(pieces[piece_index], equals_at + 1)
        if (finite_number(value)) return 1
        break
      }
    }
    return 0
  }

  function current_epoch(command, value) {
    command = "date +%s"
    command | getline value
    close(command)
    return value + 0
  }

  function report_progress(pass_number, complete, elapsed, rate, suffix) {
    elapsed = current_epoch() - started[pass_number]
    if (elapsed < 1) elapsed = 1
    rate = pass_rows[pass_number] / elapsed
    suffix = complete ? " (complete)" : ""
    printf "Phase 2 INFO/R2 PVAR pass %d: %d rows, %d s, %.0f rows/s%s\n", \
      pass_number, pass_rows[pass_number], elapsed, rate, suffix > "/dev/stderr"
  }

  function select_metric(candidate) {
    for (candidate = 1; candidate <= candidate_count; candidate++) {
      if (finite_count[candidate] > 0) return candidate
    }
    return 0
  }

  function parse_header(field_count, fields, i, upper, candidate, required) {
    delete header_column
    delete header_name
    id_column = 0
    for (i = 1; i <= field_count; i++) {
      upper = toupper(fields[i])
      sub(/^#/, "", upper)
      if (upper == "ID" && !id_column) id_column = i
      if (!(upper in header_column)) {
        header_column[upper] = i
        header_name[upper] = fields[i]
      }
    }
    if (!id_column) fail("Phase 2 Step 1 PVAR header is missing the ID column: " FILENAME)

    delete metric_column
    delete metric_name
    for (candidate = 1; candidate <= candidate_count; candidate++) {
      upper = candidate_header[candidate]
      metric_column[candidate] = header_column[upper] + 0
      metric_name[candidate] = header_name[upper]
    }
    metric_column[3] = header_column["INFO"] + 0
    metric_name[3] = header_name["INFO"] ":INFO/R2"

    required = id_column
    for (candidate = 1; candidate <= candidate_count; candidate++) {
      if (metric_column[candidate] > required) required = metric_column[candidate]
    }
    required_columns = required
    header_seen[ARGIND] = 1
  }

  function metric_number(candidate, fields) {
    if (!metric_column[candidate]) return 0
    if (candidate == 3) return info_key_number(fields[metric_column[candidate]])
    return finite_number(fields[metric_column[candidate]])
  }

  BEGIN {
    FS = OFS = "\t"
    if (!finite_number(threshold) || parsed_number < 0 || parsed_number > 1) {
      fail("Phase 2 Step 1 INFO/R2 threshold must be between 0 and 1")
    }
    threshold += 0

    key_count = split("R2 INFO MACH_R2 MINIMAC3_R2 IMPUTE_INFO IMPUTE2_INFO INFO_SCORE RSQ", info_keys, " ")
    candidate_count = split("R2 INFO INFO MACH_R2 MINIMAC3_R2 IMPUTE_INFO IMPUTE2_INFO INFO_SCORE RSQ", candidate_header, " ")

    print "variant_id", "info_metric", "info_value", "info_min", "exclusion_reason" > report_file
    close(report_file)
    printf "" > pass_file
    close(pass_file)
  }

  BEGINFILE {
    started[ARGIND] = current_epoch()
    if (ARGIND == 2) {
      chosen = select_metric()
      skip_second = !chosen || low_count[chosen] == 0
    }
  }

  {
    if (ARGIND == 2 && skip_second) nextfile
    line = $0
    sub(/\r$/, "", line)
    if (line ~ /^[[:space:]]*$/ || line ~ /^##/) next

    field_count = split(line, fields, "\t")
    if (!header_seen[ARGIND]) {
      parse_header(field_count, fields)
      next
    }
    if (field_count < required_columns) {
      fail("Phase 2 Step 1 PVAR row " FNR " has fewer columns than its header: " FILENAME)
    }

    pass_rows[ARGIND]++
    if (ARGIND == 1) {
      for (candidate = 1; candidate <= candidate_count; candidate++) {
        if (!metric_number(candidate, fields)) continue
        finite_count[candidate]++
        if (parsed_number < threshold) low_count[candidate]++
      }
    } else if (metric_number(chosen, fields) && parsed_number < threshold) {
      print fields[id_column], metric_name[chosen], sprintf("%.15g", parsed_number), \
        sprintf("%.15g", threshold), "info_r2_below_min" >> report_file
    } else {
      id = fields[id_column]
      if (id != "" && id != "." && id != "NA") {
        print id >> pass_file
        pass_count++
      }
    }

    if (pass_rows[ARGIND] % 5000000 == 0) report_progress(ARGIND, 0)
  }

  ENDFILE {
    if (ARGIND == 1 && !header_seen[1]) {
      fail("Phase 2 Step 1 PVAR header was not found: " FILENAME)
    }
    if (!(ARGIND == 2 && skip_second)) report_progress(ARGIND, 1)
  }

  END {
    if (failed) exit 2
    chosen = select_metric()
    if (chosen && low_count[chosen] > 0 && pass_count == 0) {
      print "ERROR: Step 1 INFO/R2 marker filter removed all variants from " FILENAME > "/dev/stderr"
      exit 2
    }
    print "metric", "value" > summary_file
    print "rows_scanned", pass_rows[1] + 0 >> summary_file
    print "metric_name", chosen ? metric_name[chosen] : "" >> summary_file
    print "finite_values", chosen ? finite_count[chosen] : 0 >> summary_file
    print "below_min", chosen ? low_count[chosen] : 0 >> summary_file
    print "pass_variants", pass_count + 0 >> summary_file
    close(summary_file)
  }
  ' \
  "$pvar" "$pvar"

below_min=$("$awk_bin" -F $'\t' '$1 == "below_min" { print $2 }' "$summary_temp")
below_min=${below_min:-0}

if [[ "$below_min" -gt 0 ]]; then
  mv -f -- "$pass_temp" "$pass_ids"
  pass_temp=""
  mv -f -- "$report_temp" "$excluded_report"
  report_temp=""
else
  # Do not leave stale filter artifacts when no usable metric or low value exists.
  rm -f -- "$pass_ids" "$excluded_report"
fi

mv -f -- "$summary_temp" "$summary"
summary_temp=""
