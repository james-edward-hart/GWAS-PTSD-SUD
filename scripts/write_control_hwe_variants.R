#!/usr/bin/env Rscript

# Write variants passing Hardy-Weinberg filtering in stratum controls only.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


args <- parse_args(defaults = list(threads = "1"))
require_args(args, c("config", "pheno", "keep", "plink-prefix", "out"))


config <- load_config(args$config)
pheno <- read_tsv(args$pheno)
keep <- read_tsv(args$keep)
require_columns(pheno, c("FID", "IID", "PHENO"), "phenotype file")
require_columns(keep, c("FID", "IID"), "GWAS keep file")

keep_key <- paste(keep$FID, keep$IID, sep = "\t")
control_key <- paste(pheno$FID, pheno$IID, sep = "\t")
controls <- pheno[pheno$PHENO == "1" & control_key %in% keep_key, c("FID", "IID"), drop = FALSE]
if (!nrow(controls)) die("no controls available for controls-only HWE filtering")

control_keep <- paste0(args[["plink-prefix"]], ".controls.keep.txt")
write_plink_id_file(controls, control_keep)

command <- c(
  plink_input_args(config$genotypes, args[["plink-prefix"]], "controls-only HWE genotype input"),
  "--keep", control_keep,
  control_hwe_filters(config),
  "--write-snplist",
  "--threads", args$threads,
  "--out", args[["plink-prefix"]]
)
run_command(plink_tool(config), command)

snplist <- paste0(args[["plink-prefix"]], ".snplist")
if (!file.exists(snplist)) die("PLINK2 did not write expected HWE variant list: ", snplist)
if (normalizePath(snplist, mustWork = TRUE) != normalizePath(args$out, mustWork = FALSE)) {
  ensure_parent(args$out)
  file.copy(snplist, args$out, overwrite = TRUE)
}
if (count_lines(args$out) == 0) die("controls-only HWE filter retained zero variants")

cat("Wrote", count_lines(args$out), "controls-only HWE-passing variants from", nrow(controls), "controls\n")
