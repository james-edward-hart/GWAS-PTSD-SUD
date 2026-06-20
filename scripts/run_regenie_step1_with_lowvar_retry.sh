#!/usr/bin/env bash
set -euo pipefail

die() {
  echo "ERROR: $*" >&2
  exit 1
}

usage() {
  cat >&2 <<'USAGE'
Usage: run_regenie_step1_with_lowvar_retry.sh \
  --out-prefix PATH \
  --base-extract PATH \
  --pred-list PATH \
  --done PATH \
  --max-low-variance-exclusions N \
  -- regenie [regenie step 1 arguments]
USAGE
  exit 2
}

out_prefix=""
base_extract=""
pred_list=""
done_file=""
max_low_variance_exclusions=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-prefix)
      [[ $# -ge 2 ]] || usage
      out_prefix="$2"
      shift 2
      ;;
    --base-extract)
      [[ $# -ge 2 ]] || usage
      base_extract="$2"
      shift 2
      ;;
    --pred-list)
      [[ $# -ge 2 ]] || usage
      pred_list="$2"
      shift 2
      ;;
    --done)
      [[ $# -ge 2 ]] || usage
      done_file="$2"
      shift 2
      ;;
    --max-low-variance-exclusions)
      [[ $# -ge 2 ]] || usage
      max_low_variance_exclusions="$2"
      shift 2
      ;;
    --)
      shift
      break
      ;;
    *)
      usage
      ;;
  esac
done

[[ -n "$out_prefix" ]] || usage
[[ -n "$base_extract" ]] || usage
[[ -n "$pred_list" ]] || usage
[[ -n "$done_file" ]] || usage
[[ "$max_low_variance_exclusions" =~ ^[0-9]+$ ]] || usage
[[ $# -gt 0 ]] || usage
[[ -s "$base_extract" ]] || die "base regenie Step 1 extract list is missing or empty: $base_extract"

regenie_cmd=("$@")
extract_idx=-1
out_idx=-1
lowmem_prefix="${out_prefix}.lowmem"
for i in "${!regenie_cmd[@]}"; do
  case "${regenie_cmd[$i]}" in
    --extract)
      [[ $((i + 1)) -lt ${#regenie_cmd[@]} ]] || die "regenie command has --extract without a value"
      [[ "$extract_idx" -eq -1 ]] || die "regenie command has multiple --extract options"
      extract_idx="$i"
      ;;
    --out)
      [[ $((i + 1)) -lt ${#regenie_cmd[@]} ]] || die "regenie command has --out without a value"
      [[ "$out_idx" -eq -1 ]] || die "regenie command has multiple --out options"
      out_idx="$i"
      ;;
    --lowmem-prefix)
      [[ $((i + 1)) -lt ${#regenie_cmd[@]} ]] || die "regenie command has --lowmem-prefix without a value"
      lowmem_prefix="${regenie_cmd[$((i + 1))]}"
      ;;
  esac
done

[[ "$extract_idx" -ge 0 ]] || die "regenie Step 1 command is missing --extract"
[[ "$out_idx" -ge 0 ]] || die "regenie Step 1 command is missing --out"
[[ "${regenie_cmd[$((extract_idx + 1))]}" == "$base_extract" ]] || {
  die "base extract does not match regenie command --extract value: $base_extract"
}
[[ "${regenie_cmd[$((out_idx + 1))]}" == "$out_prefix" ]] || {
  die "out prefix does not match regenie command --out value: $out_prefix"
}

runtime_extract="${out_prefix}.runtime_extract.snplist"
next_extract="${runtime_extract}.next"
observed_pred="${out_prefix}_pred.list"
low_variance_ids="${out_prefix}.low_variance_exclusions.txt"
low_variance_report="${out_prefix}.low_variance_exclusions.tsv"

mkdir -p "$(dirname "$out_prefix")" "$(dirname "$pred_list")" "$(dirname "$done_file")"
cp -f "$base_extract" "$runtime_extract"
: > "$low_variance_ids"
printf 'attempt\tvariant_id\tattempt_log\n' > "$low_variance_report"
rm -f "${out_prefix}.attempt_"*.stdout.log "${out_prefix}.attempt_"*.regenie.log

clean_partial_outputs() {
  rm -f \
    "$observed_pred" \
    "${out_prefix}.log" \
    "${lowmem_prefix}"_* \
    "${out_prefix}"_*.loco
}

extract_low_variance_snp() {
  local attempt_stdout="$1"
  {
    [[ -f "${out_prefix}.log" ]] &&
      sed -n -E 's/.*SNP ([^[:space:]]+) has low variance.*/\1/p' "${out_prefix}.log"
    [[ -f "$attempt_stdout" ]] &&
      sed -n -E 's/.*SNP ([^[:space:]]+) has low variance.*/\1/p' "$attempt_stdout"
  } | tail -n 1
}

attempt=1
while true; do
  clean_partial_outputs
  run_cmd=("${regenie_cmd[@]}")
  run_cmd[$((extract_idx + 1))]="$runtime_extract"
  attempt_stdout="${out_prefix}.attempt_${attempt}.stdout.log"
  attempt_regenie_log="${out_prefix}.attempt_${attempt}.regenie.log"
  variant_count="$(wc -l < "$runtime_extract" | tr -d '[:space:]')"
  echo "Phase 2 regenie Step 1 attempt ${attempt}: ${variant_count} variants"

  set +e
  "${run_cmd[@]}" 2>&1 | tee "$attempt_stdout"
  pipe_status=("${PIPESTATUS[@]}")
  set -e
  status="${pipe_status[0]}"
  tee_status="${pipe_status[1]}"
  [[ "$tee_status" -eq 0 ]] || die "could not write regenie attempt log: $attempt_stdout"

  if [[ "$status" -eq 0 ]]; then
    [[ -s "$observed_pred" ]] || die "regenie succeeded but prediction list is missing or empty: $observed_pred"
    if [[ "$observed_pred" != "$pred_list" ]]; then
      cp -f "$observed_pred" "$pred_list"
    fi
    printf '%s\n' 'ok' > "$done_file"
    exit 0
  fi

  if [[ -f "${out_prefix}.log" ]]; then
    cp -f "${out_prefix}.log" "$attempt_regenie_log"
  fi

  bad_snp="$(extract_low_variance_snp "$attempt_stdout")"
  if [[ -z "$bad_snp" ]]; then
    echo "ERROR: regenie Step 1 failed with exit status ${status}; no low-variance SNP was reported" >&2
    exit "$status"
  fi
  if grep -Fxq "$bad_snp" "$low_variance_ids"; then
    echo "ERROR: regenie reported a low-variance SNP that was already excluded: $bad_snp" >&2
    exit "$status"
  fi
  excluded_count="$(wc -l < "$low_variance_ids" | tr -d '[:space:]')"
  if [[ "$excluded_count" -ge "$max_low_variance_exclusions" ]]; then
    echo "ERROR: regenie Step 1 exceeded the low-variance SNP exclusion limit (${max_low_variance_exclusions})" >&2
    echo "Last low-variance SNP: $bad_snp" >&2
    exit "$status"
  fi
  if ! grep -Fxq "$bad_snp" "$runtime_extract"; then
    echo "ERROR: regenie reported a low-variance SNP absent from the runtime extract: $bad_snp" >&2
    exit "$status"
  fi

  echo "Excluding low-variance Step 1 SNP and retrying: $bad_snp" >&2
  printf '%s\n' "$bad_snp" >> "$low_variance_ids"
  printf '%s\t%s\t%s\n' "$attempt" "$bad_snp" "$attempt_stdout" >> "$low_variance_report"
  awk 'NR == FNR { exclude[$1] = 1; next } !($1 in exclude)' \
    "$low_variance_ids" "$base_extract" > "$next_extract"
  mv -f "$next_extract" "$runtime_extract"
  [[ -s "$runtime_extract" ]] || die "all Step 1 variants were removed by low-variance exclusions"
  attempt=$((attempt + 1))
done
