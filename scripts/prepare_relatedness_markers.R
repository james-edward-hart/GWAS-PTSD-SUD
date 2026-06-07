#!/usr/bin/env Rscript

# Build the QC'd LD-pruned variant set used for relatedness filtering.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Parse input keep files and output prefixes.
args <- parse_args(defaults = list(keep = "", threads = "1"))
require_args(args, c("config", "out-prefix", "prune-prefix", "summary"))


# Load relatedness-specific QC settings.
config <- load_config(args$config)
settings <- config$relatedness
qc <- config$qc


# Assemble the marker and sample filters for the relatedness subset.
filters <- character()
if (truthy(settings$autosome_only %||% qc$autosome_only)) filters <- c(filters, "--autosome")
if (truthy(settings$snps_only_acgt %||% qc$snps_only_acgt)) filters <- c(filters, "--snps-only", "just-acgt")
filters <- c(filters,
  "--max-alleles", "2",
  "--rm-dup", "exclude-all",
  "--maf", as.character(settings$maf_min %||% qc$maf_min),
  "--geno", as.character(settings$geno_missing_max %||% qc$geno_missing_max),
  "--mind", as.character(settings$sample_missing_max %||% qc$sample_missing_max)
)


# Create a QC'd temporary PGEN for relatedness pruning.
command <- c(plink_input_args(config$genotypes, args[["out-prefix"]], "relatedness genotype input"))
if (nzchar(args$keep)) command <- c(command, "--keep", args$keep)
command <- c(command, filters, "--make-pgen", "--threads", args$threads, "--out", args[["out-prefix"]])
ensure_parent(paste0(args[["out-prefix"]], ".pgen"))
run_command(plink_tool(config), command)


# LD-prune the relatedness marker set with PLINK2.
pruning <- settings$ld_prune
run_command(plink_tool(config), c(
  "--pfile", args[["out-prefix"]],
  "--indep-pairwise", as.character(pruning$window %||% "500kb"), as.character(pruning$step %||% 1), as.character(pruning$r2 %||% 0.2),
  "--threads", args$threads,
  "--out", args[["prune-prefix"]]
))

prune_in <- paste0(args[["prune-prefix"]], ".prune.in")
prune_out <- paste0(args[["prune-prefix"]], ".prune.out")
if (!file.exists(prune_in) || file.info(prune_in)$size == 0) die("relatedness LD pruning did not produce a non-empty prune.in file")


# Record marker-filter settings and pruning counts for reports.
summary <- data.frame(
  metric = c(
    "source_genotype_type", "maf_min", "geno_missing_max", "sample_missing_max",
    "ld_prune_window", "ld_prune_step", "ld_prune_r2",
    "prune_in_variants", "prune_out_variants", "qc_prefix", "prune_prefix"
  ),
  value = c(
    config$genotypes$type,
    as.character(settings$maf_min %||% qc$maf_min),
    as.character(settings$geno_missing_max %||% qc$geno_missing_max),
    as.character(settings$sample_missing_max %||% qc$sample_missing_max),
    as.character(pruning$window %||% "500kb"),
    as.character(pruning$step %||% 1),
    as.character(pruning$r2 %||% 0.2),
    count_lines(prune_in),
    if (file.exists(prune_out)) count_lines(prune_out) else 0,
    args[["out-prefix"]],
    args[["prune-prefix"]]
  ),
  stringsAsFactors = FALSE
)

# Save the summary consumed by reports.
write_tsv(summary, args$summary)
cat("Prepared relatedness marker set:", prune_in, "\n")
