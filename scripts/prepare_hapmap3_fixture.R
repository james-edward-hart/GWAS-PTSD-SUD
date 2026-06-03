#!/usr/bin/env Rscript

# Prepare the public HapMap3 toy dataset and seeded binary phenotype.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# The fixture writes ten PCs for the local example workflow.
pcs <- paste0("PC", seq_len(10))


# Read FID/IID/sex from the public HapMap3 FAM file.
read_fam <- function(path) {
  fam <- read.table(path, stringsAsFactors = FALSE, comment.char = "", col.names = c("FID", "IID", "father", "mother", "sex", "pheno"))
  fam[c("FID", "IID", "sex")]
}


# Read PLINK2 PCA output and normalize the #FID header.
read_eigenvec <- function(path) {
  rows <- read.table(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
  if ("#FID" %in% names(rows)) names(rows)[names(rows) == "#FID"] <- "FID"
  rows
}


# Compute fixture PCs from the downloaded HapMap3 genotypes.
run_pca <- function(plink2, bfile, out_prefix) {
  ensure_parent(paste0(out_prefix, ".eigenvec"))
  run_command(plink2, c(
    "--bfile", bfile,
    "--maf", "0.05",
    "--geno", "0.05",
    "--mind", "0.05",
    "--pca", "10", "approx",
    "--out", out_prefix
  ))
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


# Save study PCs in the pipeline covariate format.
write_pcs <- function(path, eigen_rows) {
  write_tsv(eigen_rows[c("FID", "IID", pcs)], path)
}


# Split samples into two toy reference labels by PC1 median.
write_reference_pcs <- function(path, eigen_rows) {
  median_pc1 <- median(as.numeric(eigen_rows$PC1))
  label <- ifelse(as.numeric(eigen_rows$PC1) <= median_pc1, "HMAP_A", "HMAP_B")
  out <- data.frame(
    FID = eigen_rows$FID,
    IID = eigen_rows$IID,
    population = paste0(label, "_REF"),
    super_population = label,
    eigen_rows[pcs],
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
  plink2 = "software/local/plink2",
  bfile = "data/example/hapmap3",
  "pca-prefix" = "data/example/hapmap3_fixture_pca",
  "sample-manifest" = "config/hapmap3_sample_manifest.tsv",
  pcs = "config/hapmap3_pcs.tsv",
  "reference-pcs" = "config/hapmap3_reference_pcs.tsv",
  traits = "config/hapmap3_traits.tsv"
))


# Build PCs before deriving fixture metadata from PCA output.
run_pca(args$plink2, args$bfile, args[["pca-prefix"]])


# Build the manifest, study PCs, reference PCs, and trait registry together.
fam_rows <- read_fam(paste0(args$bfile, ".fam"))
eigen_rows <- read_eigenvec(paste0(args[["pca-prefix"]], ".eigenvec"))
write_sample_manifest(args[["sample-manifest"]], fam_rows)
write_pcs(args$pcs, eigen_rows)
write_reference_pcs(args[["reference-pcs"]], eigen_rows)
write_traits(args$traits)


# Print a short completion line for logs.
cat("Prepared real HapMap3 fixture with random binary phenotype\n")
