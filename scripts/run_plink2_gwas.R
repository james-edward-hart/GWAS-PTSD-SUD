#!/usr/bin/env Rscript

# Run PLINK2 --glm for one trait/ancestry and harmonize the output.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Parse one trait/ancestry GWAS job from Snakemake.
args <- parse_args(defaults = list(threads = "1"))
require_args(args, c("config", "trait", "ancestry", "build", "pheno", "covar", "keep", "plink-prefix", "out"))


# Resolve config and covariate columns for the trait.
config <- load_config(args$config)
covars <- covariates_for_trait(config, args$trait)


# Build the PLINK2 --glm command.
command <- c(
  plink_input_args(config$genotypes, args[["plink-prefix"]], "GWAS genotype input"),
  "--keep", args$keep,
  "--pheno", args$pheno,
  "--pheno-name", "PHENO",
  "--covar", args$covar,
  "--covar-name", paste(covars, collapse = ",")
)
if (truthy(config$gwas$covar_variance_standardize %||% TRUE)) command <- c(command, "--covar-variance-standardize")

# Preserve free-form PLINK2 glm options from config.
glm_options <- trimws(config$gwas$glm_options %||% "")
glm_parts <- if (nzchar(glm_options)) strsplit(glm_options, "\\s+")[[1]] else character()
command <- c(command, gwas_filters(config), "--glm", glm_parts, "--threads", args$threads, "--out", args[["plink-prefix"]])


# Run PLINK2 and keep its raw output under the requested prefix.
ensure_parent(paste0(args[["plink-prefix"]], ".log"))
run_command(plink_tool(config), command)


# Harmonize PLINK2 output into the workflow summary-stat schema.
run_command("Rscript", c(
  file.path(script_dir, "harmonize_plink_glm.R"),
  "--plink-prefix", args[["plink-prefix"]],
  "--trait", args$trait,
  "--ancestry", args$ancestry,
  "--build", args$build,
  "--out", args$out,
  "--test", config$gwas$test %||% "ADD"
))
