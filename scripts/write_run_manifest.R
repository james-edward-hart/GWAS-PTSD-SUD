#!/usr/bin/env Rscript

# Write a run-level provenance manifest with checksums.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Parse the effective config and genome-build artifact paths.
args <- parse_args()
require_args(args, c("config", "genome-build-file", "out"))


# Load run metadata needed in the manifest.
config <- load_config(args$config)
genome_build <- trimws(readLines(args[["genome-build-file"]], warn = FALSE)[[1]])


# Resolve a configured executable for version/hash recording.
resolve_executable <- function(path) {
  if (blank(path)) return("")
  resolved <- if (grepl("/", path, fixed = TRUE)) path else Sys.which(path)
  if (!nzchar(resolved) || !file.exists(resolved)) return("")
  normalizePath(resolved, mustWork = TRUE)
}


# Record executable path, hash, and best-effort version text.
tool_rows <- function(name, path, version_args = character()) {
  resolved <- resolve_executable(path)
  if (!nzchar(resolved)) {
    return(data.frame(key = paste0("tool_", name, "_status"), value = "not_found", stringsAsFactors = FALSE))
  }
  version <- "NA"
  if (length(version_args)) {
    output <- tryCatch(system2(resolved, version_args, stdout = TRUE, stderr = TRUE), error = function(e) character())
    if (length(output)) version <- paste(output, collapse = " | ")
  }
  data.frame(
    key = c(paste0("tool_", name, "_path"), paste0("tool_", name, "_sha256"), paste0("tool_", name, "_version")),
    value = c(resolved, sha256_file(resolved), version),
    stringsAsFactors = FALSE
  )
}


# Record final text/image output checksums while leaving raw genotype artifacts as path/count records.
final_output_rows <- function(results_dir = "results", analysis_name = "") {
  if (!dir.exists(results_dir)) return(data.frame(key = character(), value = character()))
  current_analysis <- function(paths) {
    if (!nzchar(analysis_name) || !length(paths)) return(paths)
    paths[startsWith(basename(paths), paste0(analysis_name, "."))]
  }
  candidates <- c(
    list.files(file.path(results_dir, "qc"), recursive = TRUE, full.names = TRUE),
    current_analysis(list.files(file.path(results_dir, "reports"), recursive = TRUE, full.names = TRUE)),
    list.files(file.path(results_dir, "gwas"), pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE),
    current_analysis(list.files(file.path(results_dir, "plots"), pattern = "\\.(png|pdf)$", recursive = TRUE, full.names = TRUE))
  )
  candidates <- candidates[file.exists(candidates) & !dir.exists(candidates)]
  raw_ext <- "\\.(bed|bim|fam|pgen|pvar|psam|eigenvec|eigenval|sscore|acount|log)$"
  raw <- candidates[grepl(raw_ext, candidates)]
  final <- sort(setdiff(candidates, raw))
  rows <- data.frame(key = character(), value = character(), stringsAsFactors = FALSE)
  if (length(final)) {
    rows <- rbind(rows, data.frame(
      key = paste0("output_sha256:", final),
      value = vapply(final, sha256_file, character(1)),
      stringsAsFactors = FALSE
    ))
  }
  rows <- rbind(rows, data.frame(
    key = c("raw_genotype_artifact_count", "raw_genotype_artifact_paths"),
    value = c(length(raw), paste(sort(raw), collapse = ";")),
    stringsAsFactors = FALSE
  ))
  rows
}


input_manifest_release_rows <- function(path) {
  if (blank(path) || !file.exists(path)) return(data.frame(key = character(), value = character()))
  rows <- read_tsv(path)
  if (!all(c("file_role", "cohort_data_release") %in% names(rows))) {
    return(data.frame(key = character(), value = character()))
  }
  rows <- rows[nzchar(rows$file_role), , drop = FALSE]
  if (!nrow(rows)) return(data.frame(key = character(), value = character()))
  data.frame(
    key = paste0("input_manifest_cohort_data_release:", rows$file_role),
    value = rows$cohort_data_release,
    stringsAsFactors = FALSE
  )
}


study_genotype_component_rows <- function(config) {
  components <- genotype_component_paths(config$genotypes, "study")
  data.frame(
    key = paste0("study_genotype_component:", names(components)),
    value = unname(components),
    stringsAsFactors = FALSE
  )
}


reference_panel_rows <- function(config, section, label) {
  block <- config[[section]] %||% list()
  if (!truthy(block$enabled %||% FALSE)) return(data.frame(key = character(), value = character()))
  values <- c(
    source_panel_id = block$source_panel_id %||% "",
    reference_genome_build = block$reference_genome_build %||% "",
    reference_genotype_type = block$reference_genotypes$type %||% "",
    reference_genotype_prefix = block$reference_genotypes$prefix %||% "",
    metadata_path = block$metadata$path %||% "",
    exclusion_regions = block$exclusion_regions %||% ""
  )
  data.frame(
    key = paste0(label, ":", names(values)),
    value = unname(values),
    stringsAsFactors = FALSE
  )
}


# Record key inputs and hashes for reproducibility.
rows <- data.frame(
  key = c(
    "created_utc",
    "analysis_name",
    "cohort_data_release",
    "genome_build",
    "config",
    "config_sha256",
    "genome_build_file",
    "genome_build_file_sha256",
    "software_manifest",
    "software_manifest_sha256",
    "reference_manifest",
    "reference_manifest_sha256",
    "input_manifest",
    "input_manifest_sha256",
    "sample_manifest",
    "sample_manifest_sha256",
    "trait_registry",
    "trait_registry_sha256",
    "study_genotype_type",
    "study_genotype_prefix",
    "reference_package_root",
    "reference_package_expected_fingerprint",
    "reference_package_observed_fingerprint"
  ),
  value = c(
    format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%S+00:00", tz = "UTC"),
    config$project$analysis_name,
    config$project$cohort_data_release %||% "",
    genome_build,
    args$config,
    sha256_file(args$config),
    args[["genome-build-file"]],
    sha256_file(args[["genome-build-file"]]),
    config$resources$software_manifest,
    sha256_file(config$resources$software_manifest),
    config$resources$reference_manifest,
    sha256_file(config$resources$reference_manifest),
    config$resources$input_manifest %||% "",
    if (!blank(config$resources$input_manifest %||% "") && file.exists(config$resources$input_manifest)) {
      sha256_file(config$resources$input_manifest)
    } else "",
    config$inputs$sample_manifest,
    sha256_file(config$inputs$sample_manifest),
    config$inputs$trait_registry,
    sha256_file(config$inputs$trait_registry),
    config$genotypes$type,
    config$genotypes$prefix,
    config$reference_package$root %||% "",
    config$reference_package$fingerprint %||% "",
    config$reference_package$observed_fingerprint %||% ""
  ),
  stringsAsFactors = FALSE
)

rows <- rbind(
  rows,
  study_genotype_component_rows(config),
  input_manifest_release_rows(config$resources$input_manifest %||% ""),
  reference_panel_rows(config, "ancestry_reference", "ancestry_reference"),
  reference_panel_rows(config, "admixture", "admixture"),
  tool_rows("plink2", config$tools$plink2 %||% "plink2", "--version"),
  tool_rows("plink1", plink1_tool(config), "--version"),
  tool_rows("admixture", config$tools$admixture %||% "admixture"),
  final_output_rows("results", analysis_output_name(config))
)


# Save the manifest.
write_tsv(rows, args$out)
cat("Wrote run manifest:", args$out, "\n")
