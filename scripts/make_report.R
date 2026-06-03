#!/usr/bin/env Rscript

# Write one compact markdown report per trait and ancestry stratum.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Parse all report inputs for one trait/ancestry pair.
args <- parse_args(defaults = list("reference-prep-report" = "NA"))
require_args(args, c(
  "config", "trait", "ancestry", "build", "stats", "qq", "manhattan",
  "strata-counts", "pheno", "covar", "keep", "relatedness-summary",
  "sex-check-summary", "genome-build-details", "ancestry-counts",
  "software", "reference", "plink-log", "out"
))


# Load report inputs.
config <- load_config(args$config)
stats <- read_tsv(args$stats)
keep <- read_tsv(args$keep)
pheno <- read_tsv(args$pheno)
covar <- read_tsv(args$covar)


# Count analyzed cases and controls from the final keep file.
keep_key <- paste(keep$FID, keep$IID, sep = "\t")
pheno <- pheno[paste(pheno$FID, pheno$IID, sep = "\t") %in% keep_key, , drop = FALSE]
cases <- sum(pheno$PHENO == "2")
controls <- sum(pheno$PHENO == "1")
n_samples <- cases + controls
underpowered <- n_samples < config$warnings$min_n ||
  cases < config$warnings$min_cases ||
  controls < config$warnings$min_controls


# Summarize harmonized GWAS result content.
p <- suppressWarnings(as.numeric(stats$p))
n_variants <- nrow(stats)
min_p <- if (any(is.finite(p))) min(p[is.finite(p)], na.rm = TRUE) else "NA"
tests <- if ("test" %in% names(stats)) paste(sort(unique(stats$test[nzchar(stats$test)])), collapse = ", ") else "NA"
if (!nzchar(tests)) tests <- "NA"
covariates <- paste(names(covar)[-(1:2)], collapse = ", ")


# Read metric/value summaries into named vectors.
kv <- function(path) {
  rows <- read_tsv(path)
  setNames(rows$value, rows$metric)
}
relatedness <- kv(args[["relatedness-summary"]])
sex_check <- kv(args[["sex-check-summary"]])


# Select the chosen genome-build detail row.
build_details <- read_tsv(args[["genome-build-details"]])
selected <- build_details[build_details$selected == "True", , drop = FALSE]
if (!nrow(selected)) selected <- build_details[build_details$build == args$build, , drop = FALSE]
selected <- selected[1, , drop = FALSE]


# Pull the total number of POP-MaD exclusions.
ancestry_counts <- read_tsv(args[["ancestry-counts"]])
excluded_total <- ancestry_counts$n[ancestry_counts$category == "excluded_total"]
if (!length(excluded_total)) excluded_total <- "NA"


# Extract useful one-line highlights from the PLINK log.
log_lines <- character()
if (file.exists(args[["plink-log"]])) {
  lines <- readLines(args[["plink-log"]], warn = FALSE)
  hits <- grepl("samples|variants|removed|remaining|covariate|phenotype", lines, ignore.case = TRUE)
  log_lines <- trimws(lines[hits & nzchar(trimws(lines))])
  log_lines <- head(log_lines, 12)
}
if (!length(log_lines)) log_lines <- "No highlights captured."


input_manifest_release_lines <- function(path) {
  if (blank(path) || !file.exists(path)) return(character())
  rows <- read_tsv(path)
  if (!all(c("file_role", "cohort_data_release") %in% names(rows))) return(character())
  rows <- rows[nzchar(rows$file_role), , drop = FALSE]
  if (!nrow(rows)) return(character())
  paste0("- Input release note for ", rows$file_role, ": ", rows$cohort_data_release)
}


reference_panel_lines <- function(config, section, label) {
  block <- config[[section]] %||% list()
  if (!truthy(block$enabled %||% FALSE)) return(character())
  c(
    paste0("- ", label, " panel: ", block$source_panel_id %||% "NA"),
    paste0("- ", label, " reference build: ", block$reference_genome_build %||% "NA"),
    paste0("- ", label, " reference metadata: `", block$metadata$path %||% "NA", "`")
  )
}


study_components <- genotype_component_paths(config$genotypes, "study")


# Keep reports compact while retaining enough QC evidence for review.
text <- c(
  paste0("# GWAS Report: ", args$trait, " / ", args$ancestry),
  "",
  "## Run",
  "",
  paste0("- Analysis: ", config$project$analysis_name),
  paste0("- Trait: ", args$trait),
  paste0("- Ancestry: ", args$ancestry),
  paste0("- Genome build label: ", args$build),
  "- Engine: PLINK2 `--glm`",
  paste0("- Effective config: `", args$config, "`"),
  paste0("- PLINK log: `", args[["plink-log"]], "`"),
  "",
  "## Input Files",
  "",
  paste0("- Cohort data release: ", config$project$cohort_data_release %||% "unspecified"),
  paste0("- Sample manifest: `", config$inputs$sample_manifest, "`"),
  paste0("- Trait registry: `", config$inputs$trait_registry, "`"),
  paste0("- Study genotype type: ", config$genotypes$type),
  paste0("- Study genotype prefix: `", config$genotypes$prefix, "`"),
  paste0("- Study genotype component (", names(study_components), "): `", unname(study_components), "`"),
  paste0("- Inferred genome build: ", args$build),
  paste0("- Genome-build marker details: `", args[["genome-build-details"]], "`"),
  paste0("- Reference package root: `", config$reference_package$root %||% "NA", "`"),
  paste0("- Reference package fingerprint: ", config$reference_package$observed_fingerprint %||% config$reference_package$fingerprint %||% "NA"),
  reference_panel_lines(config, "ancestry_reference", "POP-MaD"),
  reference_panel_lines(config, "admixture", "ADMIXTURE"),
  input_manifest_release_lines(config$resources$input_manifest %||% ""),
  "",
  "## Sample Counts",
  "",
  paste0("- N: ", n_samples),
  paste0("- Cases: ", cases),
  paste0("- Controls: ", controls),
  paste0("- Underpowered warning: ", ifelse(underpowered, "True", "False")),
  paste0("- POP-MaD excluded samples: ", excluded_total[[1]]),
  "",
  "## Covariates",
  "",
  paste0("- Covariate file: `", args$covar, "`"),
  paste0("- Covariates used: ", covariates),
  paste0("- PLINK covariate variance standardization: ", ifelse(truthy(config$gwas$covar_variance_standardize), "True", "False")),
  "",
  "## Variant Results",
  "",
  paste0("- Harmonized summary statistics: `", args$stats, "`"),
  paste0("- Variants in harmonized output: ", n_variants),
  paste0("- Minimum P value: ", min_p),
  paste0("- Retained PLINK2 test terms: ", tests),
  "",
  "## Plots",
  "",
  paste0("- QQ plot: `", args$qq, "`"),
  paste0("- Manhattan plot: `", args$manhattan, "`"),
  "",
  "## QC Settings",
  "",
  paste0("- INFO/R2 filter enabled: ", ifelse(truthy(config$qc$use_mach_r2_filter %||% FALSE), "True", "False")),
  paste0("- INFO/R2 minimum when enabled: ", config$qc$info_min),
  paste0("- MAF minimum: ", config$qc$maf_min),
  paste0("- HWE P minimum: ", config$qc$hwe_p_min),
  paste0("- Genotype missingness maximum: ", config$qc$geno_missing_max),
  paste0("- Sample missingness maximum: ", config$qc$sample_missing_max),
  paste0("- Relatedness mode: ", config$relatedness$mode),
  paste0("- KING cutoff: ", config$relatedness$king_cutoff),
  paste0("- Relatedness marker set: `", args[["relatedness-summary"]], "`"),
  paste0("- Relatedness LD-pruned variants: ", relatedness[["prune_in_variants"]] %||% "NA"),
  paste0("- Sex-check action: ", sex_check[["action"]] %||% "NA"),
  paste0("- Sex-check status: ", sex_check[["status"]] %||% "NA"),
  paste0("- Sex-check problems: ", sex_check[["n_problems"]] %||% "NA"),
  paste0("- Sex-check removed samples: ", sex_check[["n_removed"]] %||% "NA"),
  paste0("- Genome-build checked markers: ", selected$checked_markers %||% "NA"),
  paste0("- Genome-build matching markers: ", selected$matching_markers %||% "NA"),
  paste0("- Genome-build match fraction: ", selected$match_fraction %||% "NA"),
  paste0("- Reference prep report: `", args[["reference-prep-report"]], "`"),
  "",
  "## PLINK Log Highlights",
  "",
  paste0("- ", log_lines),
  "",
  "## Manifests",
  "",
  paste0("- Software manifest: `", args$software, "`"),
  paste0("- Reference data manifest: `", args$reference, "`")
)


# Write the markdown report.
ensure_parent(args$out)
writeLines(text, args$out)
cat("Wrote report:", args$out, "\n")
