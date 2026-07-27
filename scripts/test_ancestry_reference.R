#!/usr/bin/env Rscript

# Focused regression checks for ancestry-reference harmonization thresholds and
# pre-LD-pruned package behavior.

tmp <- tempfile("ancestry_reference_test")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

write_pvar <- function(path, n) {
  rows <- data.frame(
    `#CHROM` = rep("1", n),
    POS = seq_len(n),
    ID = paste0("rs", seq_len(n)),
    REF = rep("A", n),
    ALT = rep("G", n),
    check.names = FALSE
  )
  write.table(rows, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

write_config <- function(path, min_shared = 10000, warn_below = 50000, variant_set = "pre_ld_pruned") {
  writeLines(c(
    "ancestry_reference:",
    paste0("  variant_set: \"", variant_set, "\""),
    paste0("  min_shared_variants: ", min_shared),
    paste0("  warn_shared_variants_below: ", warn_below),
    "  filters:",
    "    exclude_palindromic: true"
  ), path)
}

run_cmd <- function(args, log) {
  system2("Rscript", c("scripts/ancestry_reference.R", args), stdout = log, stderr = log)
}

harmonization_args <- function(stem) {
  c(
    "--mapping", file.path(tmp, paste0(stem, ".mapping.tsv")),
    "--reference-extract", file.path(tmp, paste0(stem, ".reference.extract.txt")),
    "--study-extract", file.path(tmp, paste0(stem, ".study.extract.txt")),
    "--reference-update", file.path(tmp, paste0(stem, ".reference.update.tsv")),
    "--study-update", file.path(tmp, paste0(stem, ".study.update.tsv"))
  )
}

config <- file.path(tmp, "config.yaml")
write_config(config)

ref_prefix_warn <- file.path(tmp, "ref_warn")
study_prefix_warn <- file.path(tmp, "study_warn")
write_pvar(paste0(ref_prefix_warn, ".pvar"), 49999)
write_pvar(paste0(study_prefix_warn, ".pvar"), 49999)

warn_log <- file.path(tmp, "warn.log")
status <- run_cmd(c(
  "shared-variants",
  "--config", config,
  "--reference-prefix", ref_prefix_warn,
  "--study-prefix", study_prefix_warn,
  "--out", file.path(tmp, "shared_warn.txt"),
  "--mismatch-report", file.path(tmp, "mismatch_warn.tsv"),
  harmonization_args("warn")
), warn_log)
stopifnot(identical(status, 0L))
warn_text <- readLines(warn_log)
stopifnot(any(grepl("POP-MaD study/reference overlapping variants: 49999", warn_text, fixed = TRUE)))
stopifnot(any(grepl("warning threshold is 50000", warn_text, fixed = TRUE)))
stopifnot(length(readLines(file.path(tmp, "shared_warn.txt"))) == 49999)

ref_prefix_fail <- file.path(tmp, "ref_fail")
study_prefix_fail <- file.path(tmp, "study_fail")
write_pvar(paste0(ref_prefix_fail, ".pvar"), 9999)
write_pvar(paste0(study_prefix_fail, ".pvar"), 9999)

fail_log <- file.path(tmp, "fail.log")
status <- run_cmd(c(
  "shared-variants",
  "--config", config,
  "--reference-prefix", ref_prefix_fail,
  "--study-prefix", study_prefix_fail,
  "--out", file.path(tmp, "shared_fail.txt"),
  "--mismatch-report", file.path(tmp, "mismatch_fail.tsv"),
  harmonization_args("fail")
), fail_log)
stopifnot(!identical(status, 0L))
stopifnot(any(grepl("minimum required is 10000", readLines(fail_log), fixed = TRUE)))

dup_config <- file.path(tmp, "dup_config.yaml")
write_config(dup_config, min_shared = 1, warn_below = 1)
allele_ref_prefix <- file.path(tmp, "ref_allele")
allele_study_prefix <- file.path(tmp, "study_allele")
allele_ref <- data.frame(
  `#CHROM` = c("1", "1"),
  POS = c(100, 200),
  ID = c("rs1", "rs2"),
  REF = "A",
  ALT = "G",
  check.names = FALSE
)
allele_study <- allele_ref
allele_study$ALT[[2]] <- "C"
write.table(allele_ref, paste0(allele_ref_prefix, ".pvar"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(allele_study, paste0(allele_study_prefix, ".pvar"), sep = "\t", quote = FALSE, row.names = FALSE)
allele_log <- file.path(tmp, "allele.log")
allele_mismatch <- file.path(tmp, "mismatch_allele.tsv")
status <- run_cmd(c(
  "shared-variants",
  "--config", dup_config,
  "--reference-prefix", allele_ref_prefix,
  "--study-prefix", allele_study_prefix,
  "--out", file.path(tmp, "shared_allele.txt"),
  "--mismatch-report", allele_mismatch,
  harmonization_args("allele")
), allele_log)
stopifnot(identical(status, 0L))
stopifnot(identical(readLines(file.path(tmp, "shared_allele.txt")), "rs1"))
allele_rows <- read.delim(allele_mismatch, sep = "\t", stringsAsFactors = FALSE)
stopifnot(identical(allele_rows$reason, "allele_mismatch"))

dup_ref_prefix <- file.path(tmp, "ref_dup")
dup_study_prefix <- file.path(tmp, "study_dup")
dup_rows <- data.frame(
  `#CHROM` = c("1", "1", "1"),
  POS = c(100, 200, 200),
  ID = c("rs1", "rs2", "rs3"),
  REF = "A",
  ALT = "G",
  check.names = FALSE
)
write.table(dup_rows, paste0(dup_ref_prefix, ".pvar"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(dup_rows, paste0(dup_study_prefix, ".pvar"), sep = "\t", quote = FALSE, row.names = FALSE)
dup_log <- file.path(tmp, "dup.log")
status <- run_cmd(c(
  "shared-variants",
  "--config", dup_config,
  "--reference-prefix", dup_ref_prefix,
  "--study-prefix", dup_study_prefix,
  "--out", file.path(tmp, "shared_dup.txt"),
  "--mismatch-report", file.path(tmp, "mismatch_dup.tsv"),
  harmonization_args("duplicate")
), dup_log)
stopifnot(identical(status, 0L))
stopifnot(identical(readLines(file.path(tmp, "shared_dup.txt")), "rs1"))
dup_mismatches <- read.delim(file.path(tmp, "mismatch_dup.tsv"), sep = "\t", stringsAsFactors = FALSE)
stopifnot(identical(dup_mismatches$reason, "duplicate_locus_allele_key_both"))

shared <- file.path(tmp, "shared_for_prune.txt")
writeLines(c("rs1", "rs2", "rs3"), shared)
prune_prefix <- file.path(tmp, "ld_prune", "ancestry_ld_prune")
prune_in <- paste0(prune_prefix, ".prune.in")
excluded <- paste0(prune_prefix, ".excluded_region_variants.txt")
prune_log <- file.path(tmp, "prune.log")
status <- run_cmd(c(
  "ld-prune",
  "--config", config,
  "--pfile-prefix", file.path(tmp, "unused_reference_shared"),
  "--shared-variants", shared,
  "--out-prefix", prune_prefix,
  "--prune-in", prune_in,
  "--excluded-regions", excluded,
  "--threads", "1"
), prune_log)
stopifnot(identical(status, 0L))
stopifnot(identical(readLines(prune_in), readLines(shared)))
stopifnot(file.exists(paste0(prune_prefix, ".prune.out")))
stopifnot(length(readLines(paste0(prune_prefix, ".prune.out"))) == 0)
stopifnot(file.exists(excluded))
stopifnot(length(readLines(excluded)) == 0)
stopifnot(any(grepl("pre-LD-pruned", readLines(prune_log), fixed = TRUE)))

pc_config <- file.path(tmp, "pc_config.yaml")
writeLines(c("popmad:", "  pcs: 2"), pc_config)
sscore_iid_only <- file.path(tmp, "study_projected.sscore")
writeLines(c(
  "#IID\tALLELE_CT\tPC1_AVG\tPC2_AVG",
  "I1\t100\t0.11\t0.21",
  "I2\t100\t0.12\t0.22"
), sscore_iid_only)
pcs_out <- file.path(tmp, "study_projected_pcs.tsv")
pcs_log <- file.path(tmp, "study_projected_pcs.log")
status <- run_cmd(c(
  "write-study-pcs",
  "--config", pc_config,
  "--sscore", sscore_iid_only,
  "--out", pcs_out
), pcs_log)
stopifnot(identical(status, 0L))
pcs_rows <- read.delim(pcs_out, sep = "\t", stringsAsFactors = FALSE)
stopifnot(identical(names(pcs_rows), c("FID", "IID", "PC1", "PC2")))
stopifnot(identical(pcs_rows$FID, pcs_rows$IID))
stopifnot(identical(pcs_rows$IID, c("I1", "I2")))

cat("Ancestry reference tests passed\n")
