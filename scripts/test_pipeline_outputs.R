#!/usr/bin/env Rscript

# Check that a local test run produced the expected Stage 1 outputs.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Parse expected results location and fixture labels.
args <- parse_args(defaults = list(results = "results", trait = "random_binary", build = "auto", ancestries = "HMAP_A,HMAP_B"))
results <- args$results


# Small assertions for required files and non-empty tables.
require_file <- function(path, message) {
  if (!file.exists(path)) die(message, ": ", path)
}
count_rows <- function(path) max(length(readLines(path, warn = FALSE)) - 1, 0)


# Check workflow-wide QC and provenance outputs.
build_file <- file.path(results, "qc", "genome_build", "genome_build.txt")
require_file(build_file, "missing inferred genome build")
require_file(file.path(results, "config", "effective_config.yaml"), "missing effective config snapshot")
require_file(file.path(results, "config", "resolved_config.yaml"), "missing resolved config snapshot")
require_file(file.path(results, "qc", "relatedness", "relatedness_summary.tsv"), "missing relatedness marker summary")
require_file(file.path(results, "qc", "sex", "sex_check_summary.tsv"), "missing sex-check summary")
build <- if (args$build == "auto") trimws(readLines(build_file, warn = FALSE)[[1]]) else args$build


# Check POP-MaD output tables.
assignments <- file.path(results, "qc", "ancestry", "popmad_assignments.tsv")
excluded <- file.path(results, "qc", "ancestry", "popmad_excluded.tsv")
counts <- file.path(results, "qc", "ancestry", "popmad_population_counts.tsv")
model_summary <- file.path(results, "qc", "ancestry", "population_model_summary.tsv")
require_file(assignments, "missing POP-MaD assignments")
if (count_rows(assignments) <= 0) die("POP-MaD assignments are empty")
require_file(excluded, "missing POP-MaD excluded file")
require_file(counts, "missing POP-MaD population counts")
require_file(model_summary, "missing POP-MaD model summary")


# Check each expected trait/ancestry report bundle.
for (ancestry in split_csv(args$ancestries)) {
  stats <- file.path(results, "gwas", args$trait, ancestry, paste0(args$trait, ".", ancestry, ".", build, ".plink2.glm.tsv"))
  report <- file.path(results, "reports", args$trait, paste0(args$trait, ".", ancestry, ".", build, ".report.md"))
  qq <- file.path(results, "plots", args$trait, ancestry, paste0(args$trait, ".", ancestry, ".", build, ".qq.png"))
  manhattan <- file.path(results, "plots", args$trait, ancestry, paste0(args$trait, ".", ancestry, ".", build, ".manhattan.png"))

  require_file(stats, "missing GWAS stats")
  if (count_rows(stats) <= 0) die("GWAS stats are empty: ", stats)
  header <- names(read_tsv(stats))
  for (column in c("a1_freq", "mac", "info", "test")) {
    if (!column %in% header) die("GWAS stats missing harmonized column ", column, ": ", stats)
  }
  require_file(report, "missing report")
  report_text <- readLines(report, warn = FALSE)
  if (!any(grepl("Covariates used", report_text))) die("report missing covariate section: ", report)
  if (!any(grepl("Relatedness LD-pruned variants", report_text))) die("report missing relatedness details: ", report)
  if (!any(grepl("Sex-check problems", report_text))) die("report missing sex-check details: ", report)
  if (!file.exists(qq) || file.info(qq)$size <= 0) die("missing QQ plot: ", qq)
  if (!file.exists(manhattan) || file.info(manhattan)$size <= 0) die("missing Manhattan plot: ", manhattan)
}


# Confirm the final run manifest exists.
require_file(file.path(results, "manifests", "run_manifest.tsv"), "missing run manifest")
cat("Pipeline output tests passed\n")
