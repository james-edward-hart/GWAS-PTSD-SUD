#!/usr/bin/env Rscript

# Validate config choices, input schemas, and required files before running GWAS.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Parse the config path and validation marker output.
args <- parse_args()
require_args(args, c("config", "out"))
config <- load_config(args$config)
genome_build <- config$project$inferred_genome_build %||% ""
if (!blank(args[["genome-build-file"]] %||% "")) {
  require_existing_file(args[["genome-build-file"]], "inferred genome-build file")
  genome_build <- trimws(readLines(args[["genome-build-file"]], warn = FALSE)[[1]])
}


# Fail with a labeled message when an input file is absent.
require_file <- function(path, label) {
  if (!file.exists(path)) die(label, " not found: ", path)
}


# Expand genotype prefixes into the required PLINK files.
genotype_files <- function(block, label) {
  kind <- tolower(block$type)
  prefix <- block$prefix
  files <- switch(kind,
    pgen = paste0(prefix, c(".pgen", ".pvar", ".psam")),
    bed = paste0(prefix, c(".bed", ".bim", ".fam")),
    die("unsupported ", label, " genotype type '", kind, "'. Use 'pgen' or 'bed'.")
  )
  for (path in files) require_file(path, label)
}


# When provided, the input manifest is only a compact cohort-release note. The
# config is the source of truth for paths; reports and the run manifest write
# those configured paths automatically.
validate_input_manifest <- function(path, config) {
  if (blank(path)) die("resources.input_manifest is blank")
  require_file(path, "input data manifest")
  rows <- read_tsv(path)
  if (!nrow(rows)) die("input data manifest is empty: ", path)
  require_columns(rows, c(
    "file_role", "cohort_data_release", "notes"
  ), "input data manifest")

  required_roles <- c("sample_manifest", "trait_registry", "study_genotype")
  for (role in required_roles) {
    hit <- rows$file_role == role
    if (!any(hit)) {
      die("input data manifest is missing file_role=", role)
    }
    if (sum(hit) > 1) {
      die("input data manifest has multiple rows for file_role=", role)
    }
    idx <- which(hit)
    for (field in c("cohort_data_release")) {
      if (!nzchar(rows[[field]][[idx]])) die("input data manifest row for ", role, " is missing ", field)
    }
  }
  invisible(TRUE)
}


# Confirm a configured executable is available when an optional branch needs it.
require_executable <- function(path, label, version_args = character()) {
  if (is.null(path) || !nzchar(path)) die(label, " path is empty")
  resolved <- if (grepl("/", path, fixed = TRUE)) path else Sys.which(path)
  if (!nzchar(resolved) || !file.exists(resolved)) die(label, " executable not found: ", path)
  if (file.access(resolved, mode = 1) != 0) die(label, " is not executable: ", resolved)
  if (length(version_args)) {
    status <- tryCatch(
      system2(resolved, version_args, stdout = TRUE, stderr = TRUE),
      warning = function(w) structure(character(), status = 1L),
      error = function(e) structure(character(), status = 1L)
    )
    code <- attr(status, "status") %||% 0L
    if (!identical(as.integer(code), 0L)) die(label, " did not run successfully with ", paste(version_args, collapse = " "))
  }
  resolved
}


validate_regenie_option_passthrough <- function(value, label) {
  value <- trimws(as.character(value %||% ""))
  if (!nzchar(value)) return(invisible(TRUE))
  blocked <- c(
    "--step", "--bed", "--pgen", "--bgen", "--sample", "--bgi", "--keep", "--remove",
    "--extract", "--exclude", "--extract-or", "--exclude-or", "--phenoFile", "--phenoCol", "--phenoColList",
    "--phenoExcludeList",
    "--covarFile", "--covarCol", "--covarColList", "--catCovarList", "--pred",
    "--covarExcludeList",
    "--bt", "--qt", "--out", "--bsize", "--lowmem", "--lowmem-prefix", "--threads",
    "--firth", "--approx", "--pThresh", "--minMAC", "--minINFO", "--apply-rint",
    "--gz", "--htp", "--no-split", "--chr", "--chrList", "--range", "--cc12", "--force-qt",
    "--strict", "--force-impute"
  )
  parts <- strsplit(value, "\\s+")[[1]]
  for (part in parts) {
    flag <- sub("=.*$", "", part)
    if (flag %in% blocked) {
      die(label, " cannot override pipeline-owned regenie option: ", flag)
    }
  }
  invisible(TRUE)
}


# Read sample IDs from BED/FAM or PGEN/PSAM input.
genotype_ids <- function(config) {
  kind <- tolower(config$genotypes$type)
  prefix <- config$genotypes$prefix
  if (kind == "bed") {
    fam <- read.table(paste0(prefix, ".fam"), stringsAsFactors = FALSE, quote = "", comment.char = "")
    return(paste(fam[[1]], fam[[2]], sep = "\t"))
  }
  psam <- read_tsv(paste0(prefix, ".psam"))
  ids <- table_sample_ids(psam, paste("genotype PSAM", paste0(prefix, ".psam")))
  paste(ids$FID, ids$IID, sep = "\t")
}


# Check whether the configured genotype dataset contains sex-chromosome markers.
genotype_has_sex_markers <- function(config) {
  kind <- tolower(config$genotypes$type)
  prefix <- config$genotypes$prefix
  if (kind == "bed") {
    chrom <- read.table(paste0(prefix, ".bim"), stringsAsFactors = FALSE, quote = "", comment.char = "")[[1]]
  } else {
    pvar <- read_tsv(paste0(prefix, ".pvar"))
    chrom_col <- if ("#CHROM" %in% names(pvar)) "#CHROM" else "CHROM"
    chrom <- pvar[[chrom_col]]
  }
  any(tolower(clean_chrom(chrom)) %in% c("23", "24", "x", "y"))
}


# Check required top-level config sections and simple scalar settings.
required_sections <- c("project", "analysis", "inputs", "tools", "genotypes", "genome_build",
  "popmad", "admixture", "ancestry_reference", "qc", "relatedness", "sex_check", "gwas",
  "warnings", "resources", "runtime")
missing_sections <- required_sections[vapply(required_sections, function(x) is.null(config[[x]]), logical(1))]
if (length(missing_sections)) die("config section missing: ", missing_sections[[1]])

if (!identical(config$project$genome_build, "auto")) die("project.genome_build must be 'auto'")
if (!is.null(config$project$run_mode) && !identical(config$project$run_mode, "production")) {
  die("project.run_mode is not part of the production config schema; remove project.run_mode or set it to 'production'")
}
if (!length(config$analysis$ancestries)) die("analysis.ancestries must list at least one ancestry")
min_stratum_n <- suppressWarnings(as.integer(config$analysis$min_stratum_n %||% 50))
if (!is.finite(min_stratum_n) || min_stratum_n < 1) die("analysis.min_stratum_n must be a positive integer")
max_unassigned_fraction <- suppressWarnings(as.numeric(config$popmad$max_unassigned_fraction %||% 0.07))
if (!is.finite(max_unassigned_fraction) || max_unassigned_fraction < 0 || max_unassigned_fraction > 1) {
  die("popmad.max_unassigned_fraction must be between 0 and 1")
}

deprecated_input_fields <- intersect(names(config$inputs), c(
  "ancestry_mode", "ancestry_file", "pcs_file", "projected_pcs_file", "reference_pcs_file"
))
if (length(deprecated_input_fields)) {
  die("unsupported input config field(s): ", paste(deprecated_input_fields, collapse = ", "),
    ". The workflow computes ancestry and PC covariates from the reference package; remove separate ancestry-label or PC-path fields.")
}

sex_action <- config$sex_check$action %||% "warn"
if (!sex_action %in% c("warn", "fail", "exclude")) die("sex_check.action must be warn, fail, or exclude")
invisible(sex_check_threshold_args(config$sex_check %||% list()))
if (truthy(config$relatedness$remove_sex_mismatches %||% FALSE) && sex_action != "exclude") {
  die("relatedness.remove_sex_mismatches: true requires sex_check.action: exclude")
}


# Confirm required files exist before loading their schemas.
genotype_files(config$genotypes, "genotype input")
invisible(require_executable(config$tools$plink2 %||% "plink2", "PLINK2", "--version"))
for (item in list(
  c(config$inputs$sample_manifest, "sample manifest"),
  c(config$inputs$trait_registry, "trait registry"),
  c(config$resources$software_manifest, "software manifest"),
  c(config$resources$reference_manifest, "reference manifest"),
  c(config$genome_build$marker_file, "genome-build marker file")
)) require_file(item[[1]], item[[2]])


# Load manifest tables and validate core schemas.
samples <- read_tsv(config$inputs$sample_manifest)
traits <- read_tsv(config$inputs$trait_registry)
pc_count <- as.integer(config$popmad$pcs %||% 10)
if (pc_count < 1 || pc_count > 20) die("popmad.pcs must be between 1 and 20")

require_columns(samples, c("FID", "IID", "age", "age2", "sex"), "sample manifest")
require_unique_ids(samples, "sample manifest")
require_columns(traits, c("trait_id", "phenotype_column", "case_value", "control_value", "missing_values"), "trait registry")


# Validate Phase 2 regenie settings and mixed binary/quantitative trait detection.
if (truthy(config$phase2_regenie$enabled %||% FALSE)) {
  pcs <- as.integer(config$phase2_regenie$global_pcs %||% NA)
  if (is.na(pcs) || pcs < 1 || pcs > 50) die("phase2_regenie.global_pcs must be an integer between 1 and 50")
  htp_cohort_name <- trimws(as.character(config$phase2_regenie$htp_cohort_name %||% ""))
  if (nzchar(htp_cohort_name) && !grepl("^[A-Za-z0-9._-]+$", htp_cohort_name)) {
    die("phase2_regenie.htp_cohort_name must be blank or contain only letters, numbers, dots, underscores, and hyphens")
  }
  min_mac <- as.numeric(config$phase2_regenie$min_mac %||% 1)
  if (!is.finite(min_mac) || min_mac < 1) die("phase2_regenie.min_mac must be a positive number")
  p_thresh <- as.numeric(config$phase2_regenie$p_thresh %||% 0.01)
  if (!is.finite(p_thresh) || p_thresh <= 0 || p_thresh > 1) die("phase2_regenie.p_thresh must be between 0 and 1")
  for (field in c("step1_bsize", "step2_bsize")) {
    value <- as.integer(config$phase2_regenie[[field]] %||% NA)
    if (is.na(value) || value < 1) die("phase2_regenie.", field, " must be a positive integer")
  }
  low_variance_limit <- suppressWarnings(as.integer(config$phase2_regenie$step1_low_variance_exclusion_limit %||% 25))
  if (is.na(low_variance_limit) || low_variance_limit < 0) {
    die("phase2_regenie.step1_low_variance_exclusion_limit must be a non-negative integer")
  }
  validate_regenie_option_passthrough(config$phase2_regenie$step1_options %||% "", "phase2_regenie.step1_options")
  validate_regenie_option_passthrough(config$phase2_regenie$step2_options %||% "", "phase2_regenie.step2_options")

  phase2_default_covars <- as.character(unlist(config$phase2_regenie$default_covariates %||% character(), use.names = FALSE))
  if (!length(phase2_default_covars)) phase2_default_covars <- c("age", "age2", "sex", paste0("PC", seq_len(pcs)))
  phase2_extra_covars <- as.character(unlist(config$phase2_regenie$extra_covariates %||% character(), use.names = FALSE))
  configured_phase2_covars <- unique(c(phase2_default_covars, phase2_extra_covars))
  configured_phase2_covars <- configured_phase2_covars[nzchar(configured_phase2_covars)]
  for (covar in configured_phase2_covars) {
    if (startsWith(covar, "PC")) {
      pc_idx <- suppressWarnings(as.integer(sub("^PC", "", covar)))
      if (is.na(pc_idx) || pc_idx < 1 || pc_idx > pcs) {
        die("Phase 2 PC covariate '", covar, "' is outside phase2_regenie.global_pcs=", pcs)
      }
    } else if (!covar %in% names(samples)) {
      die("Phase 2 covariate '", covar, "' is absent from sample manifest")
    }
  }

  for (branch in c("global_pca", "step1")) {
    settings <- config$phase2_regenie[[branch]]
    if (is.null(settings$filters)) die("phase2_regenie.", branch, ".filters section is required")
    if (is.null(settings$ld_prune)) die("phase2_regenie.", branch, ".ld_prune section is required")
    pruning <- settings$ld_prune
    if (is.na(suppressWarnings(as.numeric(pruning$step))) || as.numeric(pruning$step) < 1) {
      die("phase2_regenie.", branch, ".ld_prune.step must be a positive number")
    }
    if (is.na(suppressWarnings(as.numeric(pruning$r2))) || as.numeric(pruning$r2) <= 0 || as.numeric(pruning$r2) >= 1) {
      die("phase2_regenie.", branch, ".ld_prune.r2 must be between 0 and 1")
    }
    if (is.null(pruning$window) || !nzchar(as.character(pruning$window))) {
      die("phase2_regenie.", branch, ".ld_prune.window must be set")
    }
  }

  for (i in seq_len(nrow(traits))) {
    trait <- traits[i, , drop = FALSE]
    case_blank <- blank(trait$case_value[[1]])
    control_blank <- blank(trait$control_value[[1]])
    if (!case_blank && !control_blank) {
      next
    } else if (xor(case_blank, control_blank)) {
      die("Phase 2 trait autodetection requires both case_value and control_value, or both blank: ", trait$trait_id[[1]])
    }
    values <- samples[[trait$phenotype_column[[1]]]]
    values <- values[!values %in% c(split_csv(trait$missing_values[[1]]), "", "NA", "-9", ".")]
    if (!length(values)) die("Phase 2 quantitative trait has no nonmissing values: ", trait$trait_id[[1]])
    numeric_values <- suppressWarnings(as.numeric(values))
    if (any(is.na(numeric_values) | !is.finite(numeric_values))) {
      die("Phase 2 trait ", trait$trait_id[[1]], " has blank case/control values but nonnumeric phenotype values")
    }
  }
}


# Ensure sample manifest IDs exist in the genotype files.
genotype_key <- genotype_ids(config)
genotype_parts <- do.call(rbind, strsplit(genotype_key, "\t", fixed = TRUE))
genotype_table <- data.frame(FID = genotype_parts[, 1], IID = genotype_parts[, 2], stringsAsFactors = FALSE)
genotype_idx <- match_sample_rows(samples[c("FID", "IID")], sample_key_map(genotype_table, "genotype samples"))
missing_from_genotypes <- paste(samples$FID, samples$IID, sep = "\t")[is.na(genotype_idx)]
if (length(missing_from_genotypes)) {
  first <- paste(head(gsub("\t", " ", missing_from_genotypes), 5), collapse = ", ")
  die(length(missing_from_genotypes), " sample manifest IDs are absent from genotype files; first: ", first)
}


# Trait columns and non-PC covariates must be present in the manifest.
for (i in seq_len(nrow(traits))) {
  column <- traits$phenotype_column[[i]]
  if (!column %in% names(samples)) die("phenotype column '", column, "' is absent from sample manifest")
  trait_covars <- if ("covariates" %in% names(traits)) traits$covariates[[i]] else ""
  for (covar in split_csv(trait_covars)) {
    if (!covar %in% names(samples) && !startsWith(covar, "PC")) {
      die("trait covariate '", covar, "' is absent from sample manifest")
    }
  }
}


# Validate global covariate values and sex-code conventions.
manifest_covars <- c(unlist(config$gwas$default_covariates), unlist(config$gwas$extra_covariates %||% character()))
manifest_covars <- manifest_covars[!startsWith(manifest_covars, "PC")]
missing_global_covars <- setdiff(manifest_covars, names(samples))
if (length(missing_global_covars)) {
  die("global GWAS covariate(s) are absent from sample manifest: ", paste(missing_global_covars, collapse = ", "))
}
for (i in seq_len(nrow(samples))) {
  if (!samples$sex[[i]] %in% c("1", "2", "0", "NA", "-9", ".")) {
    die("sex code for ", samples$FID[[i]], " ", samples$IID[[i]], " must be 1, 2, 0, NA, -9, or .")
  }
  for (covar in intersect(manifest_covars, names(samples))) {
    value <- samples[[covar]][[i]]
    if (!value %in% c("", "NA", "-9", ".") && is.na(suppressWarnings(as.numeric(value)))) {
      die("covariate ", covar, " for ", samples$FID[[i]], " ", samples$IID[[i]], " must be numeric; observed '", value, "'")
    }
  }
}


# Validate an exclusion-region file against the inferred study build when possible.
validate_exclusion_regions <- function(path, label, build) {
  if (!nzchar(path %||% "")) return(invisible(TRUE))
  require_file(path, paste(label, "long-range LD/problem-region file"))
  rows <- read_tsv(path)
  require_columns(rows, c("chrom", "start", "end", "label"), paste(label, "exclusion regions"))
  if (nzchar(build)) {
    if ("build" %in% names(rows)) {
      observed <- unique(rows$build[nzchar(rows$build)])
      bad <- setdiff(observed, build)
      if (length(bad)) die(label, " exclusion-region build mismatch: inferred study build is ", build,
        ", file contains ", paste(bad, collapse = ", "))
    } else {
      known <- c("GRCh36", "GRCh37", "GRCh38")
      in_name <- known[vapply(known, function(x) grepl(x, basename(path), fixed = TRUE), logical(1))]
      if (length(in_name) && !build %in% in_name) {
        die(label, " exclusion-region build mismatch: inferred study build is ", build,
          ", file path suggests ", paste(in_name, collapse = ", "))
      }
      if (!length(in_name)) {
        die(label, " exclusion-region file must include a build column or build label in its filename: ", path)
      }
    }
  }
  invisible(TRUE)
}


# Validate production gates before optional branches are expanded.
ancestry_reference_enabled <- truthy(config$ancestry_reference$enabled %||% FALSE)
admixture_enabled <- truthy(config$admixture$enabled %||% FALSE)

if (!ancestry_reference_enabled) die("ancestry_reference.enabled: true is required")
if (!admixture_enabled) die("admixture.enabled: true is required")
if (blank(config$reference_package$root %||% "")) die("reference_package.root is required")
if (blank(config$reference_package$fingerprint %||% "")) die("reference_package.fingerprint is required")
if (blank(config$reference_package$observed_fingerprint %||% "")) die("a resolved reference package fingerprint is required")
if (!identical(tolower(config$reference_package$fingerprint), tolower(config$reference_package$observed_fingerprint))) {
  die("reference_package.fingerprint does not match observed_fingerprint")
}
if (truthy(config$gwas$allow_missing_pcs %||% FALSE)) die("gwas.allow_missing_pcs: true is not supported")
if (!truthy(config$sex_check$enabled %||% TRUE)) die("sex_check.enabled: true is required")
if ((config$sex_check$action %||% "warn") != "exclude") {
  die("sex_check.action: \"exclude\" is required in config/config.yaml")
}
if (!truthy(config$sex_check$allow_no_sex_markers %||% FALSE) && !genotype_has_sex_markers(config)) {
  die("sex-chromosome markers are required for sex_check.action: exclude, or set sex_check.allow_no_sex_markers: true with documented external sex QC")
}
if (!blank(config$resources$input_manifest %||% "")) {
  validate_input_manifest(config$resources$input_manifest, config)
}


# Validate the optional report-only supervised ADMIXTURE branch.
if (admixture_enabled) {
  # ADMIXTURE is report-only, but its labels still need to match the supervised
  # reference metadata so the Q-matrix columns can be interpreted safely.
  if (!identical(config$admixture$mode %||% "supervised", "supervised")) {
    die("admixture.mode must be 'supervised'")
  }
  admixture_k <- as.integer(config$admixture$k %||% NA)
  admixture_labels <- as.character(unlist(config$admixture$labels, use.names = FALSE))
  if (is.na(admixture_k) || admixture_k < 2) die("admixture.k must be an integer >= 2")
  if (length(admixture_labels) != admixture_k) {
    die("admixture.labels must contain exactly admixture.k labels")
  }
  if (any(!nzchar(admixture_labels))) die("admixture.labels cannot contain empty labels")
  if (any(duplicated(admixture_labels))) die("admixture.labels must be unique")

  invisible(require_executable(config$tools$admixture %||% "admixture", "ADMIXTURE"))
  genotype_files(config$admixture$reference_genotypes, "ADMIXTURE reference genotype input")
  if (nzchar(genome_build) && !identical(config$admixture$reference_genome_build %||% genome_build, genome_build)) {
    die("ADMIXTURE reference build mismatch: inferred study build is ", genome_build,
      ", reference build is ", config$admixture$reference_genome_build)
  }

  metadata_path <- config$admixture$metadata$path
  require_file(metadata_path, "ADMIXTURE reference metadata")
  metadata <- read_tsv(metadata_path)
  metadata_columns <- c(
    config$admixture$metadata$sample_id_column,
    config$admixture$metadata$population_column,
    config$admixture$metadata$super_population_column
  )
  fid_column <- config$admixture$metadata$fid_column %||% ""
  if (nzchar(fid_column)) metadata_columns <- c(metadata_columns, fid_column)
  require_columns(metadata, metadata_columns, "ADMIXTURE reference metadata")
  missing_labels <- setdiff(admixture_labels, unique(metadata[[config$admixture$metadata$super_population_column]]))
  if (length(missing_labels)) {
    die("ADMIXTURE reference metadata is missing configured super-population labels: ", paste(missing_labels, collapse = ", "))
  }

  invisible(require_executable(plink1_tool(config), "PLINK1", "--version"))
  min_pruned <- as.integer(config$admixture$min_pruned_variants %||% NA)
  if (is.na(min_pruned) || min_pruned < 1) die("admixture.min_pruned_variants must be a positive integer")
  pruning <- config$admixture$ld_prune
  if (is.null(pruning)) die("admixture.ld_prune section is required when ADMIXTURE is enabled")
  if (is.na(suppressWarnings(as.numeric(pruning$step))) || as.numeric(pruning$step) < 1) {
    die("admixture.ld_prune.step must be a positive number")
  }
  if (is.na(suppressWarnings(as.numeric(pruning$r2))) || as.numeric(pruning$r2) <= 0 || as.numeric(pruning$r2) >= 1) {
    die("admixture.ld_prune.r2 must be between 0 and 1")
  }
  if (is.null(pruning$window) || !nzchar(as.character(pruning$window))) {
    die("admixture.ld_prune.window must be set")
  }

  validate_exclusion_regions(config$admixture$exclusion_regions %||% "", "ADMIXTURE", genome_build)
}

# Validate the resolved package-backed ancestry reference panel.
variant_set <- config$ancestry_reference$variant_set %||% ""
if (nzchar(variant_set) && !variant_set %in% c("pre_ld_pruned", "workflow_pruned", "unpruned")) {
  die("ancestry_reference.variant_set must be pre_ld_pruned, workflow_pruned, unpruned, or empty")
}
genotype_files(config$ancestry_reference$reference_genotypes, "ancestry reference genotype input")
if (nzchar(genome_build) && !identical(config$ancestry_reference$reference_genome_build %||% genome_build, genome_build)) {
  die("ancestry reference build mismatch: inferred study build is ", genome_build,
    ", reference build is ", config$ancestry_reference$reference_genome_build)
}
metadata_path <- config$ancestry_reference$metadata$path
require_file(metadata_path, "ancestry reference metadata")
metadata <- read_tsv(metadata_path)
population_col <- config$ancestry_reference$metadata$population_column
super_col <- config$ancestry_reference$metadata$super_population_column
require_columns(metadata, c(
  config$ancestry_reference$metadata$sample_id_column,
  population_col,
  super_col
), "ancestry reference metadata")
missing_super <- setdiff(unlist(config$analysis$ancestries), unique(metadata[[super_col]]))
if (length(missing_super)) die("ancestry reference metadata is missing configured ancestry labels: ",
  paste(missing_super, collapse = ", "))
min_population_n <- as.integer(config$popmad$min_reference_population_n %||% 20)
coverage <- popmad_reference_model_coverage(
  metadata,
  population_col,
  super_col,
  config$analysis$ancestries,
  min_population_n
)
if (length(coverage$missing_super)) {
  die("ancestry reference metadata has no population with at least popmad.min_reference_population_n=", min_population_n,
    " samples for configured ancestry label(s): ", paste(coverage$missing_super, collapse = ", "))
}
if (nrow(coverage$low)) {
  labels <- paste0(coverage$low$population, " (n=", coverage$low$n, ")")
  shown <- head(labels, 30)
  if (length(labels) > length(shown)) {
    shown <- c(shown, paste0("... ", length(labels) - length(shown), " more"))
  }
  warning("fine-scale reference populations below popmad.min_reference_population_n=", min_population_n,
    " will be skipped during POP-MaD model fitting: ", paste(shown, collapse = ", "))
}
validate_exclusion_regions(config$ancestry_reference$exclusion_regions %||% "", "ancestry", genome_build)


# Write a small success marker for Snakemake.
ensure_parent(args$out)
writeLines("validation_passed", args$out)
cat("Validation passed\n")
