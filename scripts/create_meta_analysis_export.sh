#!/usr/bin/env bash

# Package only final review files and inputs required for downstream meta-analysis.

set -euo pipefail

die() {
  echo "ERROR: $*" >&2
  exit 1
}

[[ $# -eq 1 ]] || die "usage: bash scripts/create_meta_analysis_export.sh ARCHIVE.tar.gz"

requested_archive=$1
invocation_dir=$PWD
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
archive=$requested_archive
[[ $archive == /* ]] || archive="$invocation_dir/$archive"
[[ $archive == *.tar.gz ]] || die "archive name must end in .tar.gz"
[[ ! -e $archive ]] || die "archive already exists: $archive"

cd "$repo_root"
[[ -f config/config.yaml ]] || die "run this from a configured pipeline checkout"
[[ -d results ]] || die "results directory not found"

file_list=$(mktemp)
trap 'rm -f "$file_list"' EXIT

add_file() {
  [[ -f $1 ]] && printf '%s\n' "$1" >> "$file_list"
}

add_tree() {
  [[ -d $1 ]] && find "$1" -type f -print >> "$file_list"
}

add_file config/config.yaml
add_tree results/config
add_tree results/manifests
add_tree results/reports
add_tree results/plots

# This export tree is the complete HTP/LD handoff and must remain intact.
add_tree results/remeta/export

# Keep both common-variant meta-analysis streams and small reporting derivatives.
if [[ -d results/gwas ]]; then
  find results/gwas -type f \( \
    -name '*.plink2.glm.tsv' -o \
    -name '*.gwas_filter_summary.tsv' -o \
    -name '*.phase2_summary.tsv' -o \
    -name '*.association_metrics.tsv' -o \
    -name '*.top_hits.tsv' \
  \) -print >> "$file_list"
  find results/gwas -mindepth 3 -maxdepth 3 -type f \( \
    -name '*.PAN.*.regenie' -o \
    -name '*.PAN.*.regenie.gz' \
  \) -print >> "$file_list"
fi

# Aggregate QC is sufficient for review; sample-level and working data stay local.
for path in \
  results/qc/input_validation/validation.ok \
  results/qc/genome_build/genome_build.txt \
  results/qc/genome_build/genome_build_marker_matches.tsv \
  results/qc/strata/active_ancestries.tsv \
  results/qc/strata/excluded_ancestries.tsv \
  results/qc/strata/strata_counts.tsv \
  results/qc/sex/sex_check_summary.tsv \
  results/qc/relatedness/relatedness_summary.tsv \
  results/qc/ancestry/reference/reference_projection_validation.tsv \
  results/qc/ancestry/reference/reference_prep_report.md \
  results/qc/ancestry/production/popmad_population_counts.tsv \
  results/qc/ancestry/production/population_model_summary.tsv \
  results/qc/admixture/admixture_run_summary.tsv \
  results/qc/admixture/admixture_report.md \
  results/qc/phase2_regenie/pan_sample_summary.tsv \
  results/qc/phase2_regenie/trait_group_status.tsv
do
  add_file "$path"
done

LC_ALL=C sort -u "$file_list" -o "$file_list"
grep -q '\.plink2\.glm\.tsv$' "$file_list" || die "no ancestry summary statistics found"

tar -czf "$archive" -T "$file_list"
echo "Wrote meta-analysis export: $archive"
