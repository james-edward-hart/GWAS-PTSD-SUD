#!/usr/bin/env Rscript

# Resolve build-matched reference-package panels into the run config snapshot.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Parse the effective config and inferred build paths.
args <- parse_args()
require_args(args, c("config", "genome-build-file", "out"))

if (!requireNamespace("yaml", quietly = TRUE)) die("R package 'yaml' is required")


# Load the effective config and resolve any package-backed reference panels.
config <- load_config(args$config)
genome_build <- trimws(readLines(args[["genome-build-file"]], warn = FALSE)[[1]])
if (!nzchar(genome_build)) die("inferred genome-build file is empty: ", args[["genome-build-file"]])

resolved <- resolve_reference_package_config(config, genome_build, verify_hashes = TRUE)
resolved$project$inferred_genome_build <- genome_build


# Save the build-resolved config used by validation and downstream rules.
ensure_parent(args$out)
yaml::write_yaml(resolved, args$out)
cat("Resolved reference package for build", genome_build, "into", args$out, "\n")
