#!/usr/bin/env Rscript

# Focused regression checks for ADMIXTURE QC parsing and POP-MaD comparison.


tmp <- tempfile("admixture_qc_test")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)


labels <- c("AFR", "AMR", "EAS", "EUR", "SAS")
config <- file.path(tmp, "config.yaml")
writeLines(c(
  "project:",
  "  analysis_name: test",
  "admixture:",
  "  enabled: true",
  "  mode: supervised",
  "  k: 5",
  "  labels:",
  paste0("    - ", labels),
  "  reference_genotypes:",
  "    type: pgen",
  "    prefix: reference-data/1000g/1000g",
  "  metadata:",
  "    path: metadata.tsv",
  "    sample_id_column: sample_id",
  "    fid_column: ''",
  "    population_column: population",
  "    super_population_column: super_population",
  "  exclusion_regions: resources/ancestry/long_range_ld_regions.GRCh38.tsv"
), config)

dup_ref_prefix <- file.path(tmp, "admixture_ref_dup")
dup_study_prefix <- file.path(tmp, "admixture_study_dup")
dup_rows <- data.frame(
  `#CHROM` = c("1", "1", "1"),
  POS = c(100, 200, 200),
  ID = c("rs_keep", "rs_dup_a", "rs_dup_b"),
  REF = "A",
  ALT = "G",
  check.names = FALSE
)
write.table(dup_rows, paste0(dup_ref_prefix, ".pvar"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(dup_rows, paste0(dup_study_prefix, ".pvar"), sep = "\t", quote = FALSE, row.names = FALSE)
dup_log <- file.path(tmp, "dup_shared.log")
status <- system2("Rscript", c(
  "scripts/admixture_qc.R",
  "shared-variants",
  "--config", config,
  "--reference-prefix", dup_ref_prefix,
  "--study-prefix", dup_study_prefix,
  "--out", file.path(tmp, "dup_shared.txt"),
  "--mismatch-report", file.path(tmp, "dup_mismatch.tsv")
), stdout = dup_log, stderr = dup_log)
stopifnot(identical(status, 0L))
stopifnot(identical(readLines(file.path(tmp, "dup_shared.txt")), "rs_keep"))
stopifnot(any(grepl("duplicated chromosome/position mappings", readLines(dup_log), fixed = TRUE)))

pos_ref_prefix <- file.path(tmp, "admixture_ref_pos")
pos_study_prefix <- file.path(tmp, "admixture_study_pos")
pos_ref <- data.frame(
  `#CHROM` = c("1", "1"),
  POS = c(100, 200),
  ID = c("rs_keep", "rs_pos_mismatch"),
  REF = "A",
  ALT = "G",
  check.names = FALSE
)
pos_study <- pos_ref
pos_study$POS[[2]] <- 250
write.table(pos_ref, paste0(pos_ref_prefix, ".pvar"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pos_study, paste0(pos_study_prefix, ".pvar"), sep = "\t", quote = FALSE, row.names = FALSE)
pos_log <- file.path(tmp, "pos_shared.log")
pos_mismatch <- file.path(tmp, "pos_mismatch.tsv")
status <- system2("Rscript", c(
  "scripts/admixture_qc.R",
  "shared-variants",
  "--config", config,
  "--reference-prefix", pos_ref_prefix,
  "--study-prefix", pos_study_prefix,
  "--out", file.path(tmp, "pos_shared.txt"),
  "--mismatch-report", pos_mismatch
), stdout = pos_log, stderr = pos_log)
stopifnot(identical(status, 0L))
stopifnot(identical(readLines(file.path(tmp, "pos_shared.txt")), "rs_keep"))
pos_rows <- read.delim(pos_mismatch, sep = "\t", stringsAsFactors = FALSE)
stopifnot(identical(pos_rows$reason, "position_mismatch"))
stopifnot(any(grepl("chromosome/position mismatches", readLines(pos_log), fixed = TRUE)))

fake_tool <- function(path, log, body) {
  writeLines(c("#!/bin/sh", paste0("echo \"$0 $@\" >> ", shQuote(log)), body), path)
  Sys.chmod(path, mode = "0755")
}

merge_log <- file.path(tmp, "fake_merge_tools.log")
fake_plink2 <- file.path(tmp, "fake_plink2")
fake_plink1 <- file.path(tmp, "fake_plink1")
fake_tool(fake_plink2, merge_log, c(
  "args=\"$*\"",
  "out=\"\"",
  "while [ \"$#\" -gt 0 ]; do",
  "  if [ \"$1\" = \"--out\" ]; then out=\"$2\"; shift 2; else shift; fi",
  "done",
  "[ -n \"$out\" ] || exit 2",
  "mkdir -p \"$(dirname \"$out\")\"",
  "case \" $args \" in *\" --make-bed \"*) : > \"$out.bed\"; : > \"$out.bim\"; : > \"$out.fam\";; esac",
  "case \" $args \" in *\" --make-pgen \"*) : > \"$out.pgen\"; : > \"$out.pvar\"; : > \"$out.psam\";; esac",
  "exit 0"
))
fake_tool(fake_plink1, merge_log, c(
  "args=\"$*\"",
  "out=\"\"",
  "while [ \"$#\" -gt 0 ]; do",
  "  if [ \"$1\" = \"--out\" ]; then out=\"$2\"; shift 2; else shift; fi",
  "done",
  "case \" $args \" in *\" --bmerge \"*) ;; *) exit 3;; esac",
  "[ -n \"$out\" ] || exit 2",
  "mkdir -p \"$(dirname \"$out\")\"",
  ": > \"$out.bed\"; : > \"$out.bim\"; : > \"$out.fam\"",
  "exit 0"
))
merge_config <- file.path(tmp, "merge_config.yaml")
writeLines(c(
  "tools:",
  paste0("  plink2: ", shQuote(fake_plink2)),
  paste0("  plink1: ", shQuote(fake_plink1))
), merge_config)
merge_prefix <- file.path(tmp, "merged")
status <- system2("Rscript", c(
  "scripts/admixture_qc.R",
  "merge",
  "--config", merge_config,
  "--reference-prefix", file.path(tmp, "reference_pruned"),
  "--study-prefix", file.path(tmp, "study_pruned"),
  "--out-prefix", merge_prefix,
  "--threads", "3"
), stdout = file.path(tmp, "merge.log"), stderr = file.path(tmp, "merge.log"))
stopifnot(identical(status, 0L))
stopifnot(all(file.exists(paste0(merge_prefix, c(".bed", ".bim", ".fam")))))
stopifnot(all(file.exists(paste0(merge_prefix, "_pmerge", c(".pgen", ".pvar", ".psam")))))
merge_calls <- readLines(merge_log)
stopifnot(any(grepl("--bmerge", merge_calls, fixed = TRUE)))
stopifnot(!any(grepl("--pmerge", merge_calls, fixed = TRUE)))


fam <- file.path(tmp, "merged.fam")
write.table(data.frame(
  FID = c("R_AFR", "R_AMR", "R_EAS", "R_EUR", "R_SAS", "S1", "S2"),
  IID = c("R_AFR", "R_AMR", "R_EAS", "R_EUR", "R_SAS", "S1", "S2"),
  PAT = 0,
  MAT = 0,
  SEX = 0,
  PHENO = -9
), fam, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

pop <- file.path(tmp, "merged.pop")
writeLines(c("AFR", "AMR", "EAS", "EUR", "SAS", "-", "-"), pop)

q <- file.path(tmp, "merged.5.Q")
write.table(rbind(
  c(0.01, 0.98, 0.00, 0.01, 0.00),
  c(0.01, 0.00, 0.00, 0.98, 0.01),
  c(0.00, 0.01, 0.00, 0.01, 0.98),
  c(0.98, 0.01, 0.00, 0.01, 0.00),
  c(0.01, 0.00, 0.98, 0.00, 0.01),
  c(0.04, 0.82, 0.04, 0.06, 0.04),
  c(0.10, 0.05, 0.75, 0.05, 0.05)
), q, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

p <- file.path(tmp, "merged.5.P")
write.table(rbind(
  c(0.10, 0.20, 0.30, 0.40, 0.50),
  c(0.11, 0.21, 0.31, 0.41, 0.51),
  c(0.12, 0.22, 0.32, 0.42, 0.52)
), p, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

bim <- file.path(tmp, "merged.bim")
write.table(data.frame(
  chrom = c(1, 1, 2),
  id = c("rs1", "rs2", "rs3"),
  cm = 0,
  pos = c(100, 200, 300),
  a1 = c("A", "C", "G"),
  a2 = c("G", "T", "A")
), bim, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

sample_populations <- file.path(tmp, "sample_populations.tsv")
write.table(data.frame(
  FID = c("R_AFR", "R_AMR", "R_EAS", "R_EUR", "R_SAS", "S1", "S2"),
  IID = c("R_AFR", "R_AMR", "R_EAS", "R_EUR", "R_SAS", "S1", "S2"),
  sample_set = c(rep("reference", 5), "study", "study"),
  population = c("YRI", "PEL", "CHB", "CEU", "GIH", "", ""),
  super_population = c(labels, "", ""),
  admixture_pop = c(labels, "-", "-")
), sample_populations, sep = "\t", quote = FALSE, row.names = FALSE)

metadata <- file.path(tmp, "metadata.tsv")
write.table(data.frame(
  sample_id = c("R_AFR", "R_AMR", "R_EAS", "R_EUR", "R_SAS"),
  population = c("YRI", "PEL", "CHB", "CEU", "GIH"),
  super_population = labels
), metadata, sep = "\t", quote = FALSE, row.names = FALSE)

popmad <- file.path(tmp, "popmad.tsv")
write.table(data.frame(
  FID = c("S1", "S2"),
  IID = c("S1", "S2"),
  ancestry = c("AFR", "EUR"),
  population = c("YRI", "CEU"),
  mahalanobis_distance = c(1.2, 1.3),
  method = "POP-MaD",
  confidence = c(0.9, 0.8),
  status = "assigned"
), popmad, sep = "\t", quote = FALSE, row.names = FALSE)


study_out <- file.path(tmp, "study.tsv")
reference_out <- file.path(tmp, "reference.tsv")
comparison_out <- file.path(tmp, "comparison.tsv")
summary_out <- file.path(tmp, "summary.tsv")
report_out <- file.path(tmp, "report.md")

status <- system2("Rscript", c(
  "scripts/admixture_qc.R",
  "parse-report",
  "--config", config,
  "--q", q,
  "--p", p,
  "--fam", fam,
  "--bim", bim,
  "--pop", pop,
  "--sample-populations", sample_populations,
  "--metadata", metadata,
  "--popmad", popmad,
  "--study-out", study_out,
  "--reference-out", reference_out,
  "--comparison-out", comparison_out,
  "--summary-out", summary_out,
  "--report-out", report_out
))
stopifnot(status == 0L)


study <- read.delim(study_out, stringsAsFactors = FALSE, check.names = FALSE)
reference <- read.delim(reference_out, stringsAsFactors = FALSE, check.names = FALSE)
comparison <- read.delim(comparison_out, stringsAsFactors = FALSE, check.names = FALSE)
summary <- read.delim(summary_out, stringsAsFactors = FALSE, check.names = FALSE)

stopifnot(nrow(study) == 2)
stopifnot(identical(names(study)[3:7], labels))
stopifnot(identical(study$top_component, c("AFR", "SAS")))
stopifnot(identical(study$popmad_ancestry, c("AFR", "EUR")))
stopifnot(identical(study$top_matches_popmad, c("True", "False")))

stopifnot(nrow(reference) == 5)
stopifnot(identical(reference$super_population, labels))
stopifnot(identical(comparison$comparison_status, c("match", "discordant")))
stopifnot(summary$value[summary$metric == "popmad_available"] == "True")
stopifnot(summary$value[summary$metric == "popmad_discordant"] == "1")

bad_sample_populations <- file.path(tmp, "bad_sample_populations.tsv")
bad_rows <- read.delim(sample_populations, stringsAsFactors = FALSE, check.names = FALSE)
bad_rows[c(1, 2), ] <- bad_rows[c(2, 1), ]
write.table(bad_rows, bad_sample_populations, sep = "\t", quote = FALSE, row.names = FALSE)
status <- system2("Rscript", c(
  "scripts/admixture_qc.R",
  "parse-report",
  "--config", config,
  "--q", q,
  "--p", p,
  "--fam", fam,
  "--bim", bim,
  "--pop", pop,
  "--sample-populations", bad_sample_populations,
  "--metadata", metadata,
  "--study-out", file.path(tmp, "bad_study.tsv"),
  "--reference-out", file.path(tmp, "bad_reference.tsv"),
  "--comparison-out", file.path(tmp, "bad_comparison.tsv"),
  "--summary-out", file.path(tmp, "bad_summary.tsv"),
  "--report-out", file.path(tmp, "bad_report.md")
), stdout = file.path(tmp, "bad.log"), stderr = file.path(tmp, "bad.log"))
stopifnot(!identical(status, 0L))
stopifnot(any(grepl("same order as the merged FAM", readLines(file.path(tmp, "bad.log")))))

cat("ADMIXTURE QC parser tests passed\n")
