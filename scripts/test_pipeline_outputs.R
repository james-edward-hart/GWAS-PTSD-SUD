#!/usr/bin/env Rscript

# Check that a production Stage 1 run produced expected output bundles.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Parse expected results location and optional trait/ancestry filters.
args <- parse_args(defaults = list(results = "results", config = "", trait = "", build = "auto", ancestries = ""))
results <- args$results
config_path <- args$config
if (blank(config_path)) config_path <- file.path(results, "config", "resolved_config.yaml")
require_existing_file(config_path, "resolved run config")
config <- load_config(config_path)

trait_ids <- if (blank(args$trait)) {
  traits <- read_tsv(config$inputs$trait_registry)
  require_columns(traits, "trait_id", "trait registry")
  traits$trait_id[nzchar(traits$trait_id)]
} else {
  split_csv(args$trait)
}
ancestry_labels <- if (blank(args$ancestries)) {
  as.character(unlist(config$analysis$ancestries, use.names = FALSE))
} else {
  split_csv(args$ancestries)
}
if (!length(trait_ids)) die("no traits configured for output validation")
if (!length(ancestry_labels)) die("no ancestries configured for output validation")


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
popmad_dir <- file.path(results, "qc", "ancestry", "production")
assignments <- file.path(popmad_dir, "popmad_assignments.tsv")
excluded <- file.path(popmad_dir, "popmad_excluded.tsv")
counts <- file.path(popmad_dir, "popmad_population_counts.tsv")
model_summary <- file.path(popmad_dir, "population_model_summary.tsv")
require_file(assignments, "missing POP-MaD assignments")
if (count_rows(assignments) <= 0) die("POP-MaD assignments are empty")
require_file(excluded, "missing POP-MaD excluded file")
require_file(counts, "missing POP-MaD population counts")
require_file(model_summary, "missing POP-MaD model summary")


# Check each expected trait/ancestry report bundle.
for (trait in trait_ids) {
  for (ancestry in ancestry_labels) {
    stats <- file.path(results, "gwas", trait, ancestry, paste0(trait, ".", ancestry, ".", build, ".plink2.glm.tsv"))
    summary <- file.path(results, "gwas", trait, ancestry, paste0(trait, ".", ancestry, ".", build, ".gwas_filter_summary.tsv"))
    report <- file.path(results, "reports", trait, paste0(trait, ".", ancestry, ".", build, ".report.md"))
    qq <- file.path(results, "plots", trait, ancestry, paste0(trait, ".", ancestry, ".", build, ".qq.png"))
    manhattan <- file.path(results, "plots", trait, ancestry, paste0(trait, ".", ancestry, ".", build, ".manhattan.png"))
    manhattan_pdf <- file.path(results, "plots", trait, ancestry, paste0(trait, ".", ancestry, ".", build, ".manhattan.pdf"))

    require_file(stats, "missing GWAS stats")
    require_file(summary, "missing GWAS filter summary")
    if (count_rows(stats) <= 0) die("GWAS stats are empty: ", stats)
    header <- names(read_tsv(stats))
    for (column in c("a1_freq", "mac", "info", "test")) {
      if (!column %in% header) die("GWAS stats missing harmonized column ", column, ": ", stats)
    }
    require_file(report, "missing report")
    report_text <- readLines(report, warn = FALSE)
    if (!any(grepl("Covariates used", report_text))) die("report missing covariate section: ", report)
    if (!any(grepl("Source genotype variants", report_text))) die("report missing variant-flow counts: ", report)
    if (!any(grepl("Genomic inflation factor", report_text))) die("report missing lambda GC: ", report)
    if (!any(grepl("Top Association Signals", report_text))) die("report missing top-signal section: ", report)
    if (!any(grepl("!\\[QQ plot\\]", report_text))) die("report missing embedded QQ plot: ", report)
    if (!any(grepl("!\\[Manhattan plot\\]", report_text))) die("report missing embedded Manhattan plot: ", report)
    if (!any(grepl("Relatedness LD-pruned variants", report_text))) die("report missing relatedness details: ", report)
    if (!any(grepl("Sex-check problems", report_text))) die("report missing sex-check details: ", report)
    if (!file.exists(qq) || file.info(qq)$size <= 0) die("missing QQ plot: ", qq)
    if (!file.exists(manhattan) || file.info(manhattan)$size <= 0) die("missing Manhattan plot: ", manhattan)
    if (!file.exists(manhattan_pdf) || file.info(manhattan_pdf)$size <= 0) die("missing Manhattan PDF plot: ", manhattan_pdf)
  }
}


# Confirm the final run manifest exists.
require_file(file.path(results, "manifests", "run_manifest.tsv"), "missing run manifest")
cat("Pipeline output tests passed\n")
