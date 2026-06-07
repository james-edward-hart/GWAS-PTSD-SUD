#!/usr/bin/env Rscript

# Focused regression checks for pre-PLINK genotype metadata sanitization.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))

tmp <- tempfile("genotype_input_sanitizer_test")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)


# BED/BIM/FAM inputs with same-allele BIM rows should be wrapped in a temporary
# parseable prefix, and those rows should be excluded from downstream PLINK2.
bed_prefix <- file.path(tmp, "study_bed")
writeLines("dummy-bed", paste0(bed_prefix, ".bed"))
write.table(data.frame(
  FID = c("F1", "F2"),
  IID = c("I1", "I2"),
  PAT = 0,
  MAT = 0,
  SEX = 0,
  PHENO = -9
), paste0(bed_prefix, ".fam"), sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(data.frame(
  chrom = c(1, 1, 2, 2),
  id = c("rs_valid", "AnBr12640909", "rs_lowercase", "rs_valid_2"),
  cm = 0,
  pos = c(100, 200, 300, 400),
  a1 = c("A", "A", "g", "C"),
  a2 = c("G", "A", "g", "T")
), paste0(bed_prefix, ".bim"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

bed_out <- file.path(tmp, "out", "study_bed_qc")
bed_args <- plink_input_args(list(type = "bed", prefix = bed_prefix), bed_out, "test BED input")
stopifnot(identical(bed_args[[1]], "--bfile"))
stopifnot(identical(bed_args[[3]], "--exclude"))
stopifnot(file.exists(paste0(bed_args[[2]], ".bim")))

bed_report <- read_tsv(paste0(bed_out, ".invalid_bim_alleles.tsv"))
stopifnot(identical(bed_report$variant_id, c("AnBr12640909", "rs_lowercase")))
stopifnot(identical(bed_report$reason, rep("duplicate_allele_code", 2)))
bed_exclude <- readLines(paste0(bed_out, ".invalid_bim_alleles.exclude.txt"), warn = FALSE)
stopifnot(identical(bed_exclude, bed_report$replacement_variant_id))

safe_bim <- read_bim_variants(paste0(bed_args[[2]], ".bim"))
stopifnot(!any(toupper(safe_bim$allele1) == toupper(safe_bim$allele2)))
stopifnot(identical(safe_bim$variant_id[[1]], "rs_valid"))
stopifnot(grepl("^__stage1_excluded_invalid_bim_", safe_bim$variant_id[[2]]))
stopifnot(identical(safe_bim$allele1[[2]], "A"))
stopifnot(identical(safe_bim$allele2[[2]], "C"))

clean_bed_prefix <- file.path(tmp, "clean_bed")
invisible(file.copy(paste0(bed_prefix, ".bed"), paste0(clean_bed_prefix, ".bed")))
invisible(file.copy(paste0(bed_prefix, ".fam"), paste0(clean_bed_prefix, ".fam")))
write.table(data.frame(
  chrom = 1,
  id = "rs_clean",
  cm = 0,
  pos = 100,
  a1 = "A",
  a2 = "G"
), paste0(clean_bed_prefix, ".bim"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
stopifnot(identical(
  plink_input_args(list(type = "bed", prefix = clean_bed_prefix), file.path(tmp, "clean_out"), "clean BED input"),
  c("--bfile", clean_bed_prefix)
))


# PGEN/PVAR/PSAM inputs should receive the same protection for REF == ALT rows.
pgen_prefix <- file.path(tmp, "study_pgen")
writeLines("dummy-pgen", paste0(pgen_prefix, ".pgen"))
writeLines(c("#FID\tIID", "F1\tI1", "F2\tI2"), paste0(pgen_prefix, ".psam"))
writeLines(c(
  "#CHROM\tPOS\tID\tREF\tALT",
  "1\t100\trs_valid\tA\tG",
  "1\t200\trs_same\tC\tC",
  "2\t300\trs_multiallelic\tA\tA,C"
), paste0(pgen_prefix, ".pvar"))

pgen_out <- file.path(tmp, "out", "study_pgen_qc")
pgen_args <- plink_input_args(list(type = "pgen", prefix = pgen_prefix), pgen_out, "test PGEN input")
stopifnot(identical(pgen_args[[1]], "--pfile"))
stopifnot(identical(pgen_args[[3]], "--exclude"))

pgen_report <- read_tsv(paste0(pgen_out, ".invalid_pvar_alleles.tsv"))
stopifnot(identical(pgen_report$variant_id, "rs_same"))
safe_pvar <- read_tsv(paste0(pgen_args[[2]], ".pvar"))
stopifnot(!any(safe_pvar$REF == safe_pvar$ALT & !grepl(",", safe_pvar$ALT, fixed = TRUE)))
stopifnot(identical(safe_pvar$REF[[2]], "A"))
stopifnot(identical(safe_pvar$ALT[[2]], "C"))
stopifnot(identical(safe_pvar$ID[[3]], "rs_multiallelic"))


# Headered pipeline keep files should be rewritten as headerless PLINK input.
keep_path <- file.path(tmp, "sex_checked.keep.tsv")
write_tsv(data.frame(FID = c("F1", "F2"), IID = c("I1", "I2")), keep_path)
keep_args <- plink_keep_args(keep_path, file.path(tmp, "plink_keep_out"), "test keep file")
stopifnot(identical(keep_args[[1]], "--keep"))
stopifnot(identical(readLines(keep_args[[2]], warn = FALSE), c("F1\tI1", "F2\tI2")))

empty_keep_path <- file.path(tmp, "empty_keep.tsv")
write_tsv(data.frame(FID = character(), IID = character()), empty_keep_path)
empty_keep_error <- tryCatch({
  plink_keep_args(empty_keep_path, file.path(tmp, "empty_plink_keep_out"), "empty keep file")
  ""
}, error = function(err) conditionMessage(err))
stopifnot(grepl("empty keep file is empty", empty_keep_error, fixed = TRUE))

empty_tsv_path <- file.path(tmp, "zero_byte.tsv")
invisible(file.create(empty_tsv_path))
empty_tsv_error <- tryCatch({
  read_tsv(empty_tsv_path)
  ""
}, error = function(err) conditionMessage(err))
stopifnot(grepl("tab-delimited file is empty", empty_tsv_error, fixed = TRUE))


# Blank sex-check thresholds should use conventional chrX defaults instead of
# PLINK2's strict no-threshold sanity-check defaults.
default_thresholds <- sex_check_threshold_args(list())
stopifnot(identical(as.character(default_thresholds), c("max-female-xf=0.2", "min-male-xf=0.8")))
stopifnot(isTRUE(attr(default_thresholds, "using_defaults")))

custom_thresholds <- sex_check_threshold_args(list(
  max_female_xf = "0.1",
  min_male_xf = "0.7",
  max_female_yrate = "0.01",
  min_male_yrate = "0.05"
))
stopifnot(identical(as.character(custom_thresholds), c(
  "max-female-xf=0.1",
  "min-male-xf=0.7",
  "max-female-yrate=0.01",
  "min-male-yrate=0.05"
)))
stopifnot(is.null(attr(custom_thresholds, "using_defaults")))

partial_threshold_error <- tryCatch({
  sex_check_threshold_args(list(max_female_xf = "0.1"))
  ""
}, error = function(err) conditionMessage(err))
stopifnot(grepl("custom sex_check thresholds must include both", partial_threshold_error, fixed = TRUE))

partial_yrate_error <- tryCatch({
  sex_check_threshold_args(list(max_female_xf = "0.1", min_male_xf = "0.7", max_female_yrate = "0.01"))
  ""
}, error = function(err) conditionMessage(err))
stopifnot(grepl("Y-rate thresholds must include both", partial_yrate_error, fixed = TRUE))

non_numeric_threshold_error <- tryCatch({
  sex_check_threshold_args(list(max_female_xf = "low", min_male_xf = "0.7"))
  ""
}, error = function(err) conditionMessage(err))
stopifnot(grepl("sex_check.max_female_xf must be numeric", non_numeric_threshold_error, fixed = TRUE))

cat("Genotype input sanitizer tests passed\n")
