#!/usr/bin/env Rscript

# Focused regression checks for pre-PLINK genotype metadata sanitization.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))

tmp <- tempfile("genotype_input_sanitizer_test")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)


# Invoke the scanner directly when checking its exact parser failures.
run_metadata_helper <- function(kind, metadata, stem, extra = character()) {
  summary <- file.path(tmp, paste0(stem, ".summary.tsv"))
  log <- file.path(tmp, paste0(stem, ".log"))
  helper <- file.path(script_dir, "inspect_plink_metadata.sh")
  helper_args <- c(
    helper,
    "--type", kind,
    "--metadata", metadata,
    "--summary", summary,
    "--inspect-alleles",
    extra
  )
  status <- system2(
    "bash",
    vapply(helper_args, shQuote, character(1), type = "sh"),
    stdout = log,
    stderr = log
  )
  list(status = status, summary = summary, log = readLines(log, warn = FALSE))
}


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
source_bim_sha256 <- sha256_file(paste0(bed_prefix, ".bim"))

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

safe_bim <- data.frame()
invisible(scan_plink_variant_metadata(list(type = "bed", prefix = bed_args[[2]]), function(rows) {
  safe_bim <<- rbind(safe_bim, rows)
  TRUE
}))
stopifnot(!any(toupper(safe_bim$allele1) == toupper(safe_bim$allele2)))
stopifnot(identical(safe_bim$variant_id[[1]], "rs_valid"))
stopifnot(grepl("^__stage1_excluded_invalid_bim_", safe_bim$variant_id[[2]]))
stopifnot(identical(safe_bim$allele1[[2]], "A"))
stopifnot(identical(safe_bim$allele2[[2]], "C"))
stopifnot(identical(sha256_file(paste0(bed_prefix, ".bim")), source_bim_sha256))

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
stopifnot(!file.exists(file.path(tmp, "clean_out.plink_safe_input.bim")))
stopifnot(!file.exists(file.path(tmp, "clean_out.invalid_bim_alleles.tsv")))


# PGEN/PVAR/PSAM inputs should receive the same protection for REF == ALT rows.
pgen_prefix <- file.path(tmp, "study_pgen")
writeLines("dummy-pgen", paste0(pgen_prefix, ".pgen"))
writeLines(c("#FID\tIID", "F1\tI1", "F2\tI2"), paste0(pgen_prefix, ".psam"))
noise <- paste(2, seq_len(10000), paste0("noise", seq_len(10000)), "A", "G", "noise", sep = "\t")
writeLines(c(
  "##fileformat=VCFv4.2",
  "#CHROM\tPOS\tID\tREF\tALT\tINFO",
  "1\t100\trs_valid\tA\tG\tPR;keep_valid",
  noise,
  "1\t200\trs_same\tC\tC\treplace",
  "chrX\t300\trs_multiallelic\tA\tA,C\tkeep_multiallelic"
), paste0(pgen_prefix, ".pvar"))
source_pvar_sha256 <- sha256_file(paste0(pgen_prefix, ".pvar"))

pgen_out <- file.path(tmp, "out", "study_pgen_qc")
pgen_args <- plink_input_args(list(type = "pgen", prefix = pgen_prefix), pgen_out, "test PGEN input")
stopifnot(identical(pgen_args[[1]], "--pfile"))
stopifnot(identical(pgen_args[[3]], "--exclude"))

pgen_report <- read_tsv(paste0(pgen_out, ".invalid_pvar_alleles.tsv"))
stopifnot(identical(pgen_report$variant_id, "rs_same"))
stopifnot(endsWith(pgen_report$replacement_variant_id, "_10002"))

safe_pvar <- data.frame()
invisible(scan_plink_variant_metadata(list(type = "pgen", prefix = pgen_args[[2]]), function(rows) {
  selected <- rows$row_number %in% c(1L, 10002L, 10003L)
  safe_pvar <<- rbind(safe_pvar, rows[selected, , drop = FALSE])
  TRUE
}, require_alleles = TRUE))
stopifnot(!any(safe_pvar$allele1 == safe_pvar$allele2 & !grepl(",", safe_pvar$allele2, fixed = TRUE)))
stopifnot(identical(safe_pvar$allele1[safe_pvar$row_number == 10002L], "A"))
stopifnot(identical(safe_pvar$allele2[safe_pvar$row_number == 10002L], "C"))
stopifnot(identical(safe_pvar$variant_id[safe_pvar$row_number == 10003L], "rs_multiallelic"))
stopifnot(isTRUE(safe_pvar$ref_provisional[safe_pvar$row_number == 1L]))
stopifnot(!any(safe_pvar$ref_provisional[safe_pvar$row_number != 1L]))
stopifnot(identical(sha256_file(paste0(pgen_prefix, ".pvar")), source_pvar_sha256))

scan_summary <- plink_metadata_summary(
  list(type = "pgen", prefix = pgen_prefix),
  "test PGEN input",
  require_alleles = TRUE
)
stopifnot(plink_metadata_count(scan_summary, "rows_scanned") == 10003)
stopifnot(plink_metadata_count(scan_summary, "x_variants") == 1)
stopifnot(plink_metadata_count(scan_summary, "invalid_duplicate_alleles") == 1)

safe_pvar_lines <- readLines(paste0(pgen_args[[2]], ".pvar"), warn = FALSE)
stopifnot(identical(safe_pvar_lines[[1]], "##fileformat=VCFv4.2"))
stopifnot(grepl("keep_multiallelic$", safe_pvar_lines[grepl("rs_multiallelic", safe_pvar_lines)]))
stopifnot(genotype_has_chromosomes(list(type = "pgen", prefix = pgen_args[[2]]), c("23", "X")))
stopifnot(!genotype_has_chromosomes(list(type = "bed", prefix = clean_bed_prefix), c("23", "X")))


# Headered pipeline keep files should be rewritten as headerless PLINK input.
keep_path <- file.path(tmp, "sex_checked.keep.tsv")
write_tsv(data.frame(FID = c("F1", "F2"), IID = c("I1", "I2")), keep_path)
keep_args <- plink_keep_args(keep_path, file.path(tmp, "plink_keep_out"), "test keep file")
stopifnot(identical(keep_args[[1]], "--keep"))
stopifnot(identical(readLines(keep_args[[2]], warn = FALSE), c("F1\tI1", "F2\tI2")))

king_keep_path <- file.path(tmp, "unrelated.king.cutoff.in.id")
writeLines(c("0\tI1", "I2\tI2"), king_keep_path)
canonical_keep_args <- plink_keep_args(
  king_keep_path,
  file.path(tmp, "plink_canonical_keep_out"),
  "KING unrelated keep file",
  reference_ids = data.frame(FID = c("F1", "F2"), IID = c("I1", "I2"), stringsAsFactors = FALSE),
  reference_label = "test PGEN samples"
)
stopifnot(identical(canonical_keep_args[[1]], "--keep"))
stopifnot(identical(readLines(canonical_keep_args[[2]], warn = FALSE), c("F1\tI1", "F2\tI2")))

iid_only_keep_path <- file.path(tmp, "unrelated.iid_only.in.id")
writeLines(c("I1", "I2"), iid_only_keep_path)
iid_only_keep_args <- plink_keep_args(
  iid_only_keep_path,
  file.path(tmp, "plink_iid_only_keep_out"),
  "IID-only unrelated keep file",
  reference_ids = data.frame(FID = c("F1", "F2"), IID = c("I1", "I2"), stringsAsFactors = FALSE),
  reference_label = "test PGEN samples"
)
stopifnot(identical(readLines(iid_only_keep_args[[2]], warn = FALSE), c("F1\tI1", "F2\tI2")))

psam_without_fid <- data.frame(`#IID` = c("I1", "I2"), check.names = FALSE)
plink_psam_ids <- table_sample_ids(psam_without_fid, "test PSAM", missing_fid = "zero")
stopifnot(identical(plink_psam_ids$FID, c("0", "0")))
stopifnot(identical(plink_psam_ids$IID, c("I1", "I2")))

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

empty_pvar_prefix <- file.path(tmp, "empty_pvar")
invisible(file.create(paste0(empty_pvar_prefix, ".pvar")))
empty_pvar_error <- tryCatch({
  scan_plink_variant_metadata(list(type = "pgen", prefix = empty_pvar_prefix), function(rows) TRUE)
  ""
}, error = function(err) conditionMessage(err))
stopifnot(grepl("PVAR file is empty", empty_pvar_error, fixed = TRUE))

missing_id_prefix <- file.path(tmp, "missing_id")
writeLines(c("#CHROM\tPOS\tREF\tALT", "1\t100\tA\tG"), paste0(missing_id_prefix, ".pvar"))
missing_id_error <- tryCatch({
  scan_plink_variant_metadata(list(type = "pgen", prefix = missing_id_prefix), function(rows) TRUE)
  ""
}, error = function(err) conditionMessage(err))
stopifnot(grepl("missing required column(s): variant_id", missing_id_error, fixed = TRUE))
missing_id_result <- run_metadata_helper("pgen", paste0(missing_id_prefix, ".pvar"), "missing_id")
stopifnot(missing_id_result$status != 0)
stopifnot(any(grepl("missing required column(s): variant_id", missing_id_result$log, fixed = TRUE)))

chrom_header_prefix <- file.path(tmp, "chrom_header")
writeLines(c("CHROM\tPOS\tID\tREF\tALT", "chrY\t100\trsY\tA\tG"), paste0(chrom_header_prefix, ".pvar"))
stopifnot(genotype_has_chromosomes(list(type = "pgen", prefix = chrom_header_prefix), c("24", "Y")))


# Reordered columns, repeated metadata lines, chromosome aliases, and optional
# INFO are all parsed from the primary header rather than fixed positions.
reordered_prefix <- file.path(tmp, "reordered")
writeLines(c(
  "##fileformat=VCFv4.2",
  "##source=test",
  "ID\tALT\tPOS\tINFO\t#CHROM\tREF",
  "rsX\tG\t100\t.\tchr23\tA",
  "rsY\tT\t200\t.\t24\tC",
  "rsXY\tC\t300\t.\t25\tA",
  "rsPAR1\tA\t400\t.\tchrPAR1\tG",
  "rsPAR2\tG\t500\t.\tPAR2\tT"
), paste0(reordered_prefix, ".pvar"))
reordered <- plink_metadata_summary(
  list(type = "pgen", prefix = reordered_prefix),
  "reordered PVAR",
  require_alleles = TRUE
)
stopifnot(plink_metadata_count(reordered, "rows_scanned") == 5)
stopifnot(plink_metadata_count(reordered, "x_variants") == 1)
stopifnot(plink_metadata_count(reordered, "y_variants") == 1)
stopifnot(plink_metadata_count(reordered, "xy_variants") == 1)
stopifnot(plink_metadata_count(reordered, "par1_variants") == 1)
stopifnot(plink_metadata_count(reordered, "par2_variants") == 1)
stopifnot(identical(reordered$pvar_info_column, "true"))

no_info_prefix <- file.path(tmp, "no_info")
writeLines(c("#CHROM\tPOS\tID\tREF\tALT", "X\t100\trsX\tA\tG"), paste0(no_info_prefix, ".pvar"))
no_info <- plink_metadata_summary(
  list(type = "pgen", prefix = no_info_prefix),
  "PVAR without INFO",
  require_alleles = TRUE
)
stopifnot(identical(no_info$pvar_info_column, "false"))

header_only_prefix <- file.path(tmp, "header_only")
writeLines("#CHROM\tPOS\tID\tREF\tALT", paste0(header_only_prefix, ".pvar"))
header_only <- plink_metadata_summary(
  list(type = "pgen", prefix = header_only_prefix),
  "empty eligible PVAR",
  require_alleles = TRUE
)
stopifnot(plink_metadata_count(header_only, "rows_scanned") == 0)


# Malformed rows and invalid positions fail before a summary is published.
short_pvar <- file.path(tmp, "short.pvar")
writeLines(c("#CHROM\tPOS\tID\tREF\tALT", "X\t100\trs_short\tA"), short_pvar)
short_result <- run_metadata_helper("pgen", short_pvar, "short")
stopifnot(short_result$status != 0)
stopifnot(any(grepl("has fewer columns than its header", short_result$log, fixed = TRUE)))
stopifnot(!file.exists(short_result$summary))

bad_position_pvar <- file.path(tmp, "bad_position.pvar")
writeLines(c("#CHROM\tPOS\tID\tREF\tALT", "X\tzero\trs_bad\tA\tG"), bad_position_pvar)
position_result <- run_metadata_helper("pgen", bad_position_pvar, "bad_position")
stopifnot(position_result$status != 0)
stopifnot(any(grepl("non-integer or nonpositive position", position_result$log, fixed = TRUE)))
stopifnot(!file.exists(position_result$summary))

malformed_bim_prefix <- file.path(tmp, "malformed_bim")
writeLines("1\trs1\t0\t100\tA", paste0(malformed_bim_prefix, ".bim"))
malformed_bim_error <- tryCatch({
  scan_plink_variant_metadata(list(type = "bed", prefix = malformed_bim_prefix), function(rows) TRUE)
  ""
}, error = function(err) conditionMessage(err))
stopifnot(grepl("must contain at least 6 columns", malformed_bim_error, fixed = TRUE))
malformed_bim_result <- run_metadata_helper(
  "bed",
  paste0(malformed_bim_prefix, ".bim"),
  "malformed_bim"
)
stopifnot(malformed_bim_result$status != 0)
stopifnot(any(grepl("must contain at least 6 columns", malformed_bim_result$log, fixed = TRUE)))

bad_position_bim <- file.path(tmp, "bad_position.bim")
writeLines("23\trs_bad\t0\t-1\tA\tG", bad_position_bim)
bad_bim_result <- run_metadata_helper("bed", bad_position_bim, "bad_bim_position")
stopifnot(bad_bim_result$status != 0)
stopifnot(any(grepl("non-integer or nonpositive position", bad_bim_result$log, fixed = TRUE)))


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
