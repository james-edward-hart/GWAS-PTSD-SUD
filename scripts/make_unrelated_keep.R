#!/usr/bin/env Rscript

# Create the unrelated sample keep file with KING or all-samples mode.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Parse relatedness inputs and PLINK output paths.
args <- parse_args(defaults = list("pfile-prefix" = "", extract = "", threads = "1"))
require_args(args, c("config", "out", "plink-out-prefix"))


# Choose the configured relatedness strategy.
config <- load_config(args$config)
mode <- config$relatedness$mode

# Let PLINK2 choose unrelated samples with KING when requested.
if (mode == "plink2_king") {
  dataset <- if (nzchar(args[["pfile-prefix"]])) c("--pfile", args[["pfile-prefix"]]) else genotype_args(config$genotypes)
  command <- c(dataset)
  if (nzchar(args$extract)) command <- c(command, "--extract", args$extract)
  command <- c(command,
    "--king-cutoff", as.character(config$relatedness$king_cutoff),
    "--threads", args$threads,
    "--out", args[["plink-out-prefix"]]
  )
  ensure_parent(paste0(args[["plink-out-prefix"]], ".log"))
  run_command(plink_tool(config), command)

# Keep every sample when relatedness pruning is disabled.
} else if (mode == "all_samples") {
  samples <- read_tsv(config$inputs$sample_manifest)
  write_tsv(samples[c("FID", "IID")], args$out)
} else {
  die("unsupported relatedness mode: ", mode)
}


# Fail early if the keep file was not produced.
if (!file.exists(args$out) || file.info(args$out)$size == 0) {
  die("unrelated keep file was not created: ", args$out)
}

cat("Wrote unrelated keep file:", args$out, "\n")
