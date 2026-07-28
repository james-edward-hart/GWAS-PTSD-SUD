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

harmonization_args <- function(stem) {
  c(
    "--mapping", file.path(tmp, paste0(stem, ".mapping.tsv")),
    "--reference-extract", file.path(tmp, paste0(stem, ".reference.extract.txt")),
    "--study-extract", file.path(tmp, paste0(stem, ".study.extract.txt")),
    "--reference-update", file.path(tmp, paste0(stem, ".reference.update.tsv")),
    "--study-update", file.path(tmp, paste0(stem, ".study.update.tsv"))
  )
}

dup_ref_prefix <- file.path(tmp, "admixture_ref_dup")
dup_study_prefix <- file.path(tmp, "admixture_study_dup")
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
dup_log <- file.path(tmp, "dup_shared.log")
status <- system2("Rscript", c(
  "scripts/admixture_qc.R",
  "shared-variants",
  "--config", config,
  "--reference-prefix", dup_ref_prefix,
  "--study-prefix", dup_study_prefix,
  "--out", file.path(tmp, "dup_shared.txt"),
  "--mismatch-report", file.path(tmp, "dup_mismatch.tsv"),
  harmonization_args("duplicate")
), stdout = dup_log, stderr = dup_log)
stopifnot(identical(status, 0L))
stopifnot(identical(readLines(file.path(tmp, "dup_shared.txt")), "rs1"))
dup_mismatches <- read.delim(file.path(tmp, "dup_mismatch.tsv"), sep = "\t", stringsAsFactors = FALSE)
stopifnot(identical(dup_mismatches$reason, "duplicate_locus_allele_key_both"))

allele_ref_prefix <- file.path(tmp, "admixture_ref_allele")
allele_study_prefix <- file.path(tmp, "admixture_study_allele")
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
allele_log <- file.path(tmp, "allele_shared.log")
allele_mismatch <- file.path(tmp, "allele_mismatch.tsv")
status <- system2("Rscript", c(
  "scripts/admixture_qc.R",
  "shared-variants",
  "--config", config,
  "--reference-prefix", allele_ref_prefix,
  "--study-prefix", allele_study_prefix,
  "--out", file.path(tmp, "allele_shared.txt"),
  "--mismatch-report", allele_mismatch,
  harmonization_args("allele")
), stdout = allele_log, stderr = allele_log)
stopifnot(identical(status, 0L))
stopifnot(identical(readLines(file.path(tmp, "allele_shared.txt")), "rs1"))
allele_rows <- read.delim(allele_mismatch, sep = "\t", stringsAsFactors = FALSE)
stopifnot(identical(allele_rows$reason, "allele_mismatch"))

fake_tool <- function(path, log, body) {
  writeLines(c("#!/bin/sh", paste0("echo \"$0 $@\" >> ", shQuote(log)), body), path)
  Sys.chmod(path, mode = "0755")
}

# TEMPORARY WORKAROUND TESTS: remove these with the GRCh38 reference-package
# workaround after the repaired package is installed.
temporary_reference_prefix <- file.path(tmp, "temporary_reference")
temporary_reference_ids <- sprintf("R%04d", seq_len(630))
stopifnot(file.create(paste0(temporary_reference_prefix, ".pgen")))
write.table(data.frame(
  `#CHROM` = "1", POS = 100, ID = "rs1", REF = "A", ALT = "G",
  check.names = FALSE
), paste0(temporary_reference_prefix, ".pvar"),
sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(`#IID` = temporary_reference_ids, check.names = FALSE),
  paste0(temporary_reference_prefix, ".psam"),
  sep = "\t", quote = FALSE, row.names = FALSE)

temporary_metadata <- file.path(tmp, "temporary_reference_metadata.tsv")
write.table(data.frame(
  sample_id = temporary_reference_ids[[1]],
  population = "YRI",
  super_population = "AFR"
), temporary_metadata, sep = "\t", quote = FALSE, row.names = FALSE)

temporary_convert_log <- file.path(tmp, "temporary_reference_convert.log")
temporary_plink2 <- file.path(tmp, "temporary_reference_plink2")
fake_tool(temporary_plink2, temporary_convert_log, c(
  "args=\"$*\"",
  "case \" $args \" in *\" --remove \"*) ;; *) exit 6;; esac",
  "out=\"\"",
  "while [ \"$#\" -gt 0 ]; do",
  "  if [ \"$1\" = \"--out\" ]; then out=\"$2\"; shift 2; else shift; fi",
  "done",
  "[ -n \"$out\" ] || exit 2",
  ": > \"$out.pgen\"; : > \"$out.pvar\"; : > \"$out.psam\"",
  "exit 0"
))
temporary_config <- file.path(tmp, "temporary_reference_config.yaml")
writeLines(c(
  "tools:",
  paste0("  plink2: ", shQuote(temporary_plink2)),
  "reference_package:",
  "  observed_fingerprint: 532a1e34598aa8c92ca72a8c3983dd06c82f19ea0562c45cd75d31054647976f",
  "admixture:",
  "  source_panel_id: admixture_grch38",
  "  reference_genome_build: GRCh38",
  "  reference_genotypes:",
  "    type: pgen",
  paste0("    prefix: ", shQuote(temporary_reference_prefix)),
  "  metadata:",
  paste0("    path: ", shQuote(temporary_metadata)),
  "    sample_id_column: sample_id",
  "    fid_column: ''",
  "  filters: {}"
), temporary_config)

temporary_remove <- file.path(tmp, "temporary_reference.remove.tsv")
status <- system2("Rscript", c(
  "scripts/admixture_qc.R",
  "temporary-reference-remove",
  "--config", temporary_config,
  "--out", temporary_remove
), stdout = file.path(tmp, "temporary_reference_remove.log"),
stderr = file.path(tmp, "temporary_reference_remove.log"))
stopifnot(identical(status, 0L))
temporary_remove_rows <- read.table(temporary_remove, sep = "\t", stringsAsFactors = FALSE)
stopifnot(nrow(temporary_remove_rows) == 629L)
stopifnot(identical(as.character(temporary_remove_rows[[2]]), temporary_reference_ids[-1]))

status <- system2("Rscript", c(
  "scripts/admixture_qc.R",
  "convert-reference",
  "--config", temporary_config,
  "--remove", temporary_remove,
  "--out-prefix", file.path(tmp, "temporary_reference_out")
), stdout = file.path(tmp, "temporary_reference_convert.stdout.log"),
stderr = file.path(tmp, "temporary_reference_convert.stdout.log"))
stopifnot(identical(status, 0L))
stopifnot(any(grepl(temporary_remove, readLines(temporary_convert_log), fixed = TRUE)))

write.table(data.frame(
  sample_id = temporary_reference_ids[1:2],
  population = "YRI",
  super_population = "AFR"
), temporary_metadata, sep = "\t", quote = FALSE, row.names = FALSE)
unexpected_remove_log <- file.path(tmp, "unexpected_temporary_reference_remove.log")
status <- system2("Rscript", c(
  "scripts/admixture_qc.R",
  "temporary-reference-remove",
  "--config", temporary_config,
  "--out", file.path(tmp, "unexpected_temporary_reference.remove.tsv")
), stdout = unexpected_remove_log, stderr = unexpected_remove_log)
stopifnot(!identical(status, 0L))
stopifnot(any(grepl("is not the known 629-sample", readLines(unexpected_remove_log), fixed = TRUE)))

write.table(data.frame(
  sample_id = temporary_reference_ids,
  population = "YRI",
  super_population = "AFR"
), temporary_metadata, sep = "\t", quote = FALSE, row.names = FALSE)
matched_remove <- file.path(tmp, "matched_temporary_reference.remove.tsv")
status <- system2("Rscript", c(
  "scripts/admixture_qc.R",
  "temporary-reference-remove",
  "--config", temporary_config,
  "--out", matched_remove
), stdout = file.path(tmp, "matched_temporary_reference_remove.log"),
stderr = file.path(tmp, "matched_temporary_reference_remove.log"))
stopifnot(identical(status, 0L))
stopifnot(file.info(matched_remove)$size == 0)

convert_log <- file.path(tmp, "fake_convert_tool.log")
fake_convert_plink2 <- file.path(tmp, "fake_convert_plink2")
fake_tool(fake_convert_plink2, convert_log, c(
  "args=\"$*\"",
  "case \" $args \" in *\" --keep \"*) ;; *) exit 5;; esac",
  "out=\"\"",
  "while [ \"$#\" -gt 0 ]; do",
  "  if [ \"$1\" = \"--out\" ]; then out=\"$2\"; shift 2; else shift; fi",
  "done",
  "[ -n \"$out\" ] || exit 2",
  "mkdir -p \"$(dirname \"$out\")\"",
  ": > \"$out.pgen\"; : > \"$out.pvar\"; : > \"$out.psam\"",
  "exit 0"
))
convert_prefix <- file.path(tmp, "study_source")
stopifnot(file.create(paste0(convert_prefix, ".pgen")))
write.table(data.frame(
  `#CHROM` = "1",
  POS = 100,
  ID = "rs1",
  REF = "A",
  ALT = "G",
  check.names = FALSE
), paste0(convert_prefix, ".pvar"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(FID = c("F1", "F2"), IID = c("I1", "I2")), paste0(convert_prefix, ".psam"),
  sep = "\t", quote = FALSE, row.names = FALSE)
convert_config <- file.path(tmp, "convert_config.yaml")
writeLines(c(
  "tools:",
  paste0("  plink2: ", shQuote(fake_convert_plink2)),
  "genotypes:",
  "  type: pgen",
  paste0("  prefix: ", shQuote(convert_prefix)),
  "admixture:",
  "  filters: {}"
), convert_config)

missing_keep_log <- file.path(tmp, "missing_convert_keep.log")
status <- system2("Rscript", c(
  "scripts/admixture_qc.R",
  "convert-study",
  "--config", convert_config,
  "--out-prefix", file.path(tmp, "missing_keep_out")
), stdout = missing_keep_log, stderr = missing_keep_log)
stopifnot(!identical(status, 0L))
stopifnot(any(grepl("missing required argument(s): keep", readLines(missing_keep_log), fixed = TRUE)))

convert_keep <- file.path(tmp, "convert.keep.tsv")
write.table(data.frame(FID = c("0", "I2"), IID = c("I1", "I2")), convert_keep,
  sep = "\t", quote = FALSE, row.names = FALSE)
convert_out <- file.path(tmp, "convert_out")
status <- system2("Rscript", c(
  "scripts/admixture_qc.R",
  "convert-study",
  "--config", convert_config,
  "--keep", convert_keep,
  "--out-prefix", convert_out,
  "--threads", "2"
), stdout = file.path(tmp, "convert.log"), stderr = file.path(tmp, "convert.log"))
stopifnot(identical(status, 0L))
canonical_keep <- paste0(convert_out, ".plink_keep.txt")
stopifnot(file.exists(canonical_keep))
canonical_rows <- read.table(canonical_keep, sep = "\t", stringsAsFactors = FALSE)
stopifnot(identical(canonical_rows[[1]], c("F1", "F2")))
stopifnot(identical(canonical_rows[[2]], c("I1", "I2")))
convert_calls <- readLines(convert_log)
stopifnot(any(grepl("--keep", convert_calls, fixed = TRUE)))
stopifnot(any(grepl(canonical_keep, convert_calls, fixed = TRUE)))

write_minimal_pgen <- function(prefix, psam_rows) {
  stopifnot(file.create(paste0(prefix, ".pgen")))
  write.table(data.frame(
    `#CHROM` = "1",
    POS = 100,
    ID = "rs1",
    REF = "A",
    ALT = "G",
    check.names = FALSE
  ), paste0(prefix, ".pvar"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(psam_rows, paste0(prefix, ".psam"), sep = "\t", quote = FALSE, row.names = FALSE)
}

fidless_cases <- list(
  hash_iid = data.frame(`#IID` = c("J1", "J2"), check.names = FALSE),
  iid = data.frame(IID = c("K1", "K2"), check.names = FALSE)
)
for (case_name in names(fidless_cases)) {
  fidless_prefix <- file.path(tmp, paste0("study_source_", case_name))
  write_minimal_pgen(fidless_prefix, fidless_cases[[case_name]])
  fidless_config <- file.path(tmp, paste0("convert_config_", case_name, ".yaml"))
  writeLines(c(
    "tools:",
    paste0("  plink2: ", shQuote(fake_convert_plink2)),
    "genotypes:",
    "  type: pgen",
    paste0("  prefix: ", shQuote(fidless_prefix)),
    "admixture:",
    "  filters: {}"
  ), fidless_config)
  fidless_ids <- fidless_cases[[case_name]][[1]]
  fidless_keep <- file.path(tmp, paste0("convert_", case_name, ".keep.tsv"))
  write.table(data.frame(IID = fidless_ids), fidless_keep, sep = "\t", quote = FALSE, row.names = FALSE)
  fidless_out <- file.path(tmp, paste0("convert_", case_name, "_out"))
  status <- system2("Rscript", c(
    "scripts/admixture_qc.R",
    "convert-study",
    "--config", fidless_config,
    "--keep", fidless_keep,
    "--out-prefix", fidless_out,
    "--threads", "2"
  ), stdout = file.path(tmp, paste0("convert_", case_name, ".log")),
    stderr = file.path(tmp, paste0("convert_", case_name, ".log")))
  stopifnot(identical(status, 0L))
  fidless_canonical <- read.table(paste0(fidless_out, ".plink_keep.txt"), sep = "\t", stringsAsFactors = FALSE)
  stopifnot(identical(as.character(fidless_canonical[[1]]), c("0", "0")))
  stopifnot(identical(as.character(fidless_canonical[[2]]), fidless_ids))
}

absent_keep <- file.path(tmp, "absent.keep.tsv")
write.table(data.frame(FID = "F9", IID = "I9"), absent_keep,
  sep = "\t", quote = FALSE, row.names = FALSE)
absent_keep_log <- file.path(tmp, "absent_convert_keep.log")
status <- system2("Rscript", c(
  "scripts/admixture_qc.R",
  "convert-study",
  "--config", convert_config,
  "--keep", absent_keep,
  "--out-prefix", file.path(tmp, "absent_keep_out"),
  "--threads", "2"
), stdout = absent_keep_log, stderr = absent_keep_log)
stopifnot(!identical(status, 0L))
stopifnot(any(grepl("contains samples absent from ADMIXTURE study genotype samples",
  readLines(absent_keep_log), fixed = TRUE)))

duplicate_stratum_keep <- file.path(tmp, "duplicate_stratum.keep.tsv")
write.table(data.frame(FID = c("0", "0"), IID = c("I1", "I1")), duplicate_stratum_keep,
  sep = "\t", quote = FALSE, row.names = FALSE)
duplicate_stratum_keep_log <- file.path(tmp, "duplicate_convert_keep.log")
status <- system2("Rscript", c(
  "scripts/admixture_qc.R",
  "convert-study",
  "--config", convert_config,
  "--keep", duplicate_stratum_keep,
  "--out-prefix", file.path(tmp, "duplicate_keep_out"),
  "--threads", "2"
), stdout = duplicate_stratum_keep_log, stderr = duplicate_stratum_keep_log)
stopifnot(!identical(status, 0L))
stopifnot(any(grepl("contains duplicate samples", readLines(duplicate_stratum_keep_log), fixed = TRUE)))

merge_log <- file.path(tmp, "fake_merge_tools.log")
fake_plink2 <- file.path(tmp, "fake_plink2")
fake_plink1 <- file.path(tmp, "fake_plink1")
fake_tool(fake_plink2, merge_log, c(
  "args=\"$*\"",
  "case \" $args \" in *\" --make-bed \"*\" --sort-vars \"*) exit 4;; esac",
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
plink2_calls <- merge_calls[grepl(basename(fake_plink2), merge_calls, fixed = TRUE)]
bed_calls <- plink2_calls[grepl("--make-bed", plink2_calls, fixed = TRUE)]
stopifnot(length(bed_calls) == 2L)
stopifnot(!any(grepl("--sort-vars", bed_calls, fixed = TRUE)))
stopifnot(all(grepl("--mind 0.999999", bed_calls, fixed = TRUE)))


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

write_pop_fam <- file.path(tmp, "write_pop.fam")
write.table(data.frame(
  FID = c("0", "S1"),
  IID = c("R_AFR", "S1"),
  PAT = 0,
  MAT = 0,
  SEX = 0,
  PHENO = -9
), write_pop_fam, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(data.frame(IID = c("R_AFR", "R_FILTERED")), paste0(file.path(tmp, "write_pop_ref"), ".psam"),
  sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(FID = c("S1", "S_FILTERED"), IID = c("S1", "S_FILTERED")), paste0(file.path(tmp, "write_pop_study"), ".psam"),
  sep = "\t", quote = FALSE, row.names = FALSE)
status <- system2("Rscript", c(
  "scripts/admixture_qc.R",
  "write-pop",
  "--config", config,
  "--fam", write_pop_fam,
  "--reference-prefix", file.path(tmp, "write_pop_ref"),
  "--study-prefix", file.path(tmp, "write_pop_study"),
  "--metadata", metadata,
  "--pop-out", file.path(tmp, "write_pop.pop"),
  "--sample-populations", file.path(tmp, "write_pop_samples.tsv")
), stdout = file.path(tmp, "write_pop.log"), stderr = file.path(tmp, "write_pop.log"))
stopifnot(identical(status, 0L))
stopifnot(identical(readLines(file.path(tmp, "write_pop.pop")), c("AFR", "-")))
write_pop_samples <- read.delim(file.path(tmp, "write_pop_samples.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(identical(write_pop_samples$sample_set, c("reference", "study")))

popmad <- file.path(tmp, "popmad.tsv")
write.table(data.frame(
  FID = c("0", "S2"),
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
  "--analysis-ancestry", "AFR",
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
stopifnot(identical(study$popmad_stratum, c("AFR", "AFR")))
stopifnot(identical(names(study)[4:8], labels))
stopifnot(identical(study$top_component, c("AFR", "SAS")))
stopifnot(identical(study$popmad_ancestry, c("AFR", "EUR")))
stopifnot(identical(study$top_matches_popmad, c("True", "False")))

stopifnot(nrow(reference) == 5)
stopifnot(identical(reference$popmad_stratum, rep("AFR", 5)))
stopifnot(identical(reference$super_population, labels))
stopifnot(identical(comparison$popmad_stratum, c("AFR", "AFR")))
stopifnot(identical(comparison$comparison_status, c("match", "discordant")))
stopifnot(identical(unique(summary$popmad_stratum), "AFR"))
summary_value <- function(metric) {
  value <- summary$value[summary$metric == metric]
  stopifnot(length(value) == 1)
  value
}
stopifnot(summary_value("popmad_available") == "True")
stopifnot(summary_value("popmad_discordant") == "1")
stopifnot(summary_value("popmad_comparable_samples") == "2")
stopifnot(summary_value("popmad_match_rate") == "0.500000")
stopifnot(summary_value("mean_study_proportion_AFR") == "0.435000")
stopifnot(summary_value("mean_study_proportion_SAS") == "0.395000")
stopifnot(summary_value("n_study_top_component_AFR") == "1")
stopifnot(summary_value("n_study_top_component_SAS") == "1")
report_text <- readLines(report_out, warn = FALSE)
stopifnot(any(grepl("Study Mean Proportions", report_text, fixed = TRUE)))

study_out_eur <- file.path(tmp, "study_eur.tsv")
reference_out_eur <- file.path(tmp, "reference_eur.tsv")
comparison_out_eur <- file.path(tmp, "comparison_eur.tsv")
summary_out_eur <- file.path(tmp, "summary_eur.tsv")
report_out_eur <- file.path(tmp, "report_eur.md")
status <- system2("Rscript", c(
  "scripts/admixture_qc.R",
  "parse-report",
  "--config", config,
  "--analysis-ancestry", "EUR",
  "--q", q,
  "--p", p,
  "--fam", fam,
  "--bim", bim,
  "--pop", pop,
  "--sample-populations", sample_populations,
  "--metadata", metadata,
  "--popmad", popmad,
  "--study-out", study_out_eur,
  "--reference-out", reference_out_eur,
  "--comparison-out", comparison_out_eur,
  "--summary-out", summary_out_eur,
  "--report-out", report_out_eur
))
stopifnot(status == 0L)

combined_study <- file.path(tmp, "combined_study.tsv")
combined_reference <- file.path(tmp, "combined_reference.tsv")
combined_comparison <- file.path(tmp, "combined_comparison.tsv")
combined_summary <- file.path(tmp, "combined_summary.tsv")
combined_report <- file.path(tmp, "combined_report.md")
status <- system2("Rscript", c(
  "scripts/admixture_qc.R",
  "combine-reports",
  "--config", config,
  "--study", study_out, study_out_eur,
  "--reference", reference_out, reference_out_eur,
  "--comparison", comparison_out, comparison_out_eur,
  "--summary", summary_out, summary_out_eur,
  "--report", report_out, report_out_eur,
  "--study-out", combined_study,
  "--reference-out", combined_reference,
  "--comparison-out", combined_comparison,
  "--summary-out", combined_summary,
  "--report-out", combined_report
))
stopifnot(status == 0L)
combined_study_rows <- read.delim(combined_study, stringsAsFactors = FALSE, check.names = FALSE)
combined_summary_rows <- read.delim(combined_summary, stringsAsFactors = FALSE, check.names = FALSE)
combined_counts <- table(combined_study_rows$popmad_stratum)
stopifnot(identical(names(combined_counts), c("AFR", "EUR")))
stopifnot(identical(as.integer(combined_counts), c(2L, 2L)))
stopifnot(identical(sort(unique(combined_summary_rows$popmad_stratum)), c("AFR", "EUR")))
combined_report_text <- readLines(combined_report, warn = FALSE)
stopifnot(any(grepl("run separately within each active POP-MaD stratum", combined_report_text, fixed = TRUE)))

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
