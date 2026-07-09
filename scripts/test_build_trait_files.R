#!/usr/bin/env Rscript

# Focused checks for trait-file sample ID joins.

cmd <- commandArgs(FALSE)
script_path <- sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])
repo <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(repo, "scripts", "lib", "stage1.R"))

tmp <- tempfile("trait-files-")
dir.create(tmp, recursive = TRUE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

write_lines <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path)
}

samples <- file.path(tmp, "samples.tsv")
traits <- file.path(tmp, "traits.tsv")
pcs <- file.path(tmp, "pcs.tsv")
keep <- file.path(tmp, "keep.tsv")
config <- file.path(tmp, "config.yaml")
pheno <- file.path(tmp, "trait.pheno.tsv")
covar <- file.path(tmp, "trait.covar.tsv")

write_lines(c(
  "FID\tIID\tage\tsex\ttrait",
  "F1\tI1\t40\t1\t1",
  "F2\tI2\t42\t2\t0"
), samples)
write_lines(c(
  "trait_id\tphenotype_column\tcase_value\tcontrol_value\tmissing_values\tcovariates",
  "bt\ttrait\t1\t0\tNA\t"
), traits)
write_lines(c(
  "FID\tIID\tancestry\tPC1",
  "I1\tI1\tEUR\t0.11",
  "I2\tI2\tEUR\t0.12"
), pcs)
write_lines(c(
  "FID\tIID",
  "0\tI1",
  "0\tI2"
), keep)
write_lines(c(
  "inputs:",
  paste0("  sample_manifest: ", samples),
  paste0("  trait_registry: ", traits),
  "gwas:",
  "  default_covariates: [age, sex, PC1]",
  "  extra_covariates: []",
  "  allow_missing_pcs: false"
), config)

status <- system2("Rscript", c(
  file.path(repo, "scripts", "build_trait_files.R"),
  "--config", config,
  "--trait", "bt",
  "--pcs-file", pcs,
  "--keep", keep,
  "--pheno-out", pheno,
  "--covar-out", covar
))
stopifnot(identical(status, 0L))

pheno_rows <- read_tsv(pheno)
covar_rows <- read_tsv(covar)
stopifnot(identical(as.character(pheno_rows$PHENO), c("2", "1")))
stopifnot(identical(as.character(covar_rows$age), c("40", "42")))
stopifnot(identical(as.character(covar_rows$PC1), c("0.11", "0.12")))

cat("Trait file tests passed\n")
