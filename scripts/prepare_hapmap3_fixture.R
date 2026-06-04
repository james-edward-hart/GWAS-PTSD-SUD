#!/usr/bin/env Rscript

# Prepare the public HapMap3 fixture metadata and seeded binary phenotype.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Read FID/IID/sex from the public HapMap3 FAM file.
read_fam <- function(path) {
  fam <- read.table(path, stringsAsFactors = FALSE, comment.char = "", col.names = c("FID", "IID", "father", "mother", "sex", "pheno"))
  fam[c("FID", "IID", "sex")]
}


# Generate deterministic toy phenotype and covariates.
write_sample_manifest <- function(path, fam_rows) {
  set.seed(20260521)
  ages <- 25 + (((seq_len(nrow(fam_rows)) - 1) * 7) %% 45)
  centered <- ages - 47
  out <- data.frame(
    FID = fam_rows$FID,
    IID = fam_rows$IID,
    age = ages,
    age2 = centered * centered,
    sex = fam_rows$sex,
    ptsd_sud_case = ifelse(runif(nrow(fam_rows)) < 0.5, "1", "0"),
    stringsAsFactors = FALSE
  )
  write_tsv(out, path)
}


# Register the deterministic binary test trait.
write_traits <- function(path) {
  out <- data.frame(
    trait_id = "random_binary",
    phenotype_column = "ptsd_sud_case",
    case_value = "1",
    control_value = "0",
    missing_values = "NA,-9,.",
    description = "Random binary phenotype for real HapMap3 genotype fixture",
    stringsAsFactors = FALSE
  )
  write_tsv(out, path)
}


# Parse output locations for generated fixture files.
args <- parse_args(defaults = list(
  bfile = "data/example/hapmap3",
  "sample-manifest" = "config/example_sample_manifest.tsv",
  traits = "config/example_traits.tsv"
))


# Build the manifest and trait registry together.
fam_rows <- read_fam(paste0(args$bfile, ".fam"))
write_sample_manifest(args[["sample-manifest"]], fam_rows)
write_traits(args$traits)


# Print a short completion line for logs.
cat("Prepared real HapMap3 fixture with random binary phenotype\n")
