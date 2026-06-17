#!/usr/bin/env Rscript

# Write variants passing Hardy-Weinberg filtering in the appropriate Stage 1
# sample subset: controls for binary traits, all nonmissing samples for
# quantitative traits.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


args <- parse_args(defaults = list(threads = "1"))
require_args(args, c("config", "trait", "pheno", "keep", "plink-prefix", "out"))


config <- load_config(args$config)
traits <- read_tsv(config$inputs$trait_registry)
pheno <- read_tsv(args$pheno)
keep <- read_tsv(args$keep)
require_columns(pheno, c("FID", "IID", "PHENO"), "phenotype file")
require_columns(keep, c("FID", "IID"), "GWAS keep file")

trait <- traits[traits$trait_id == args$trait, , drop = FALSE]
if (!nrow(trait)) die("unknown trait_id: ", args$trait)
case_value <- trait$case_value[[1]]
control_value <- trait$control_value[[1]]
is_binary_trait <- !blank(case_value) && !blank(control_value)
if (!is_binary_trait && (!blank(case_value) || !blank(control_value))) {
  die("trait ", args$trait, " must set both case_value and control_value for binary analysis, or leave both blank for quantitative analysis")
}

keep_key <- paste(keep$FID, keep$IID, sep = "\t")
control_key <- paste(pheno$FID, pheno$IID, sep = "\t")
if (is_binary_trait) {
  hwe_samples <- pheno[pheno$PHENO == "1" & control_key %in% keep_key, c("FID", "IID"), drop = FALSE]
  hwe_label <- "controls-only"
  if (!nrow(hwe_samples)) die("no controls available for controls-only HWE filtering")
} else {
  hwe_samples <- pheno[pheno$PHENO != "NA" & control_key %in% keep_key, c("FID", "IID"), drop = FALSE]
  hwe_label <- "quantitative nonmissing-sample"
  if (!nrow(hwe_samples)) die("no nonmissing quantitative samples available for HWE filtering")
}

control_keep <- paste0(args[["plink-prefix"]], ".controls.keep.txt")
write_plink_id_file(hwe_samples, control_keep)

command <- c(
  plink_input_args(config$genotypes, args[["plink-prefix"]], paste(hwe_label, "HWE genotype input")),
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
if (count_lines(args$out) == 0) die(hwe_label, " HWE filter retained zero variants")

cat("Wrote", count_lines(args$out), "HWE-passing variants from", nrow(hwe_samples), hwe_label, "samples\n")
