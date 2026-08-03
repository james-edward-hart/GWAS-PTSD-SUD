#!/usr/bin/env Rscript

# Summarize one GWAS job for final reports without rerunning association tests.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


args <- parse_args()
require_args(args, c("config", "stats", "plink-log", "hwe-snplist", "hwe-log", "pheno", "covar", "keep", "out"))


count_records <- function(path, header = FALSE, drop_metadata = FALSE) {
  if (!file.exists(path)) return(NA_integer_)
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines)]
  if (drop_metadata) lines <- lines[!startsWith(lines, "##")]
  if (header && length(lines)) lines <- lines[-1]
  length(lines)
}

genotype_record_counts <- function(block) {
  kind <- tolower(block$type)
  prefix <- block$prefix
  if (kind == "bed") {
    return(list(
      source_samples = count_records(paste0(prefix, ".fam")),
      source_variants = count_records(paste0(prefix, ".bim"))
    ))
  }
  if (kind == "pgen") {
    return(list(
      source_samples = count_records(paste0(prefix, ".psam"), header = TRUE, drop_metadata = TRUE),
      source_variants = count_records(paste0(prefix, ".pvar"), header = TRUE, drop_metadata = TRUE)
    ))
  }
  die("unsupported genotype type '", kind, "'. Use 'pgen' or 'bed'.")
}

clean_log_number <- function(value) {
  if (is.null(value) || !length(value)) return("NA")
  value <- as.character(value)
  value[is.na(value) | !nzchar(value)] <- "NA"
  gsub(",", "", value)
}

extract_log_value <- function(lines, patterns) {
  for (pattern in patterns) {
    hit <- grep(pattern, lines, value = TRUE, perl = TRUE)
    if (length(hit)) {
      match <- regmatches(hit[[1]], regexec(pattern, hit[[1]], perl = TRUE))[[1]]
      if (length(match) >= 2) return(clean_log_number(match[[2]]))
    }
  }
  "NA"
}

extract_log_pair <- function(lines, patterns) {
  for (pattern in patterns) {
    hit <- grep(pattern, lines, value = TRUE, perl = TRUE)
    if (length(hit)) {
      match <- regmatches(hit[[1]], regexec(pattern, hit[[1]], perl = TRUE))[[1]]
      if (length(match) >= 3) return(clean_log_number(match[2:3]))
    }
  }
  c("NA", "NA")
}

config <- load_config(args$config)
genotype_counts <- genotype_record_counts(config$genotypes)
keep <- read_tsv(args$keep)
pheno <- read_tsv(args$pheno)
covar <- read_tsv(args$covar)
stats <- read_tsv(args$stats)
p <- suppressWarnings(as.numeric(stats$p))
valid_p <- is.finite(p) & p > 0 & p <= 1

log_lines <- if (file.exists(args[["plink-log"]])) readLines(args[["plink-log"]], warn = FALSE) else character()
hwe_log_lines <- if (file.exists(args[["hwe-log"]])) readLines(args[["hwe-log"]], warn = FALSE) else character()
initial_variant_filter <- extract_log_pair(log_lines, c(
  "([0-9,]+) excluded by .*?, ([0-9,]+) remaining\\.?$"
))
lambda_gc <- genomic_lambda(p)

summary <- data.frame(
  metric = c(
    "source_genotype_type",
    "source_genotype_samples",
    "source_genotype_variants",
    "sample_manifest_rows",
    "stratum_keep_samples",
    "phenotype_rows",
    "covariate_rows",
    "plink_loaded_samples",
    "plink_keep_remaining_samples",
    "plink_final_samples",
    "plink_loaded_variants",
    "plink_initial_filter_excluded_variants",
    "plink_initial_filter_remaining_variants",
    "control_hwe_controls",
    "control_hwe_passing_variants",
    "plink_geno_removed_variants",
    "plink_maf_removed_variants",
    "plink_hwe_removed_variants",
    "plink_info_removed_variants",
    "plink_final_variants",
    "harmonized_analysis_variants",
    "valid_p_value_variants",
    "genomewide_significant_variants",
    "suggestive_variants",
    "lambda_gc",
    "plink_log",
    "control_hwe_log"
  ),
  value = as.character(c(
    config$genotypes$type,
    genotype_counts$source_samples,
    genotype_counts$source_variants,
    count_records(config$inputs$sample_manifest, header = TRUE),
    nrow(keep),
    nrow(pheno),
    nrow(covar),
    extract_log_value(log_lines, c("^([0-9,]+) samples .* loaded from")),
    extract_log_value(log_lines, c("^--keep: ([0-9,]+) samples remaining")),
    extract_log_value(log_lines, c("^([0-9,]+) (?:samples|people).* remaining after main filters\\.?")),
    extract_log_value(log_lines, c("^([0-9,]+) variants (?:loaded from|in\\b)")),
    initial_variant_filter[[1]],
    initial_variant_filter[[2]],
    extract_log_value(hwe_log_lines, c("^--keep: ([0-9,]+) samples remaining")),
    count_records(args[["hwe-snplist"]]),
    extract_log_value(log_lines, c("^--geno: ([0-9,]+) variants? removed")),
    extract_log_value(log_lines, c("^--maf: ([0-9,]+) variants? removed", "^([0-9,]+) variants removed due to allele frequency threshold")),
    extract_log_value(hwe_log_lines, c("^--hwe: ([0-9,]+) variants? removed")),
    extract_log_value(log_lines, c("^--mach-r2-filter: ([0-9,]+) variants? removed")),
    extract_log_value(log_lines, c("^([0-9,]+) variants remaining after main filters\\.?")),
    nrow(stats),
    sum(valid_p),
    sum(valid_p & p <= 5e-8),
    sum(valid_p & p <= 1e-5),
    ifelse(is.finite(lambda_gc), sprintf("%.6f", lambda_gc), "NA"),
    args[["plink-log"]],
    args[["hwe-log"]]
  )),
  stringsAsFactors = FALSE
)

write_tsv(summary, args$out)
cat("Wrote GWAS filter summary:", args$out, "\n")
