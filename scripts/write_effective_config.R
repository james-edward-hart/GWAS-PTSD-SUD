#!/usr/bin/env Rscript

# Save the merged Snakemake config used by the current run.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))

# Snakemake copies script: files, so fall back to the repository helper path.
stage1_path <- file.path(script_dir, "lib", "stage1.R")
if (!file.exists(stage1_path)) stage1_path <- file.path(getwd(), "scripts", "lib", "stage1.R")
source(stage1_path)


# Hash tracked pipeline files into one reproducibility fingerprint.
code_fingerprint <- function(paths) {
  paths <- sort(as.character(paths))
  hashes <- vapply(paths, function(path) paste(path, sha256_file(path), sep = "\t"), character(1))
  sha256_text(paste(hashes, collapse = "\n"))
}


# Write the resolved config with optional pipeline-code metadata.
write_config <- function(config, out, tracked_code = character()) {
  if (!requireNamespace("yaml", quietly = TRUE)) die("R package 'yaml' is required")
  config <- as.list(config)
  if (length(tracked_code)) {
    tracked_code <- sort(as.character(tracked_code))
    config[["_pipeline_code"]] <- list(
      fingerprint = code_fingerprint(tracked_code),
      tracked_files = tracked_code
    )
  }
  ensure_parent(out)
  yaml::write_yaml(config, out)
}


# Support both Snakemake script: execution and direct CLI use.
if (exists("snakemake", inherits = FALSE)) {
  # Snakemake provides the resolved config and the tracked code file list.
  write_config(snakemake@config, snakemake@output[["config"]], unlist(snakemake@input[["code"]], use.names = FALSE))
} else {
  args <- parse_args(repeated = "code")
  require_args(args, c("source", "out"))
  tracked <- args$code %||% character()
  write_config(load_config(args$source), args$out, tracked)
}
