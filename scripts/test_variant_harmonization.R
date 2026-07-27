#!/usr/bin/env Rscript

# Focused tests for ancestry-only rsID/CPRA harmonization.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))
source(file.path(script_dir, "lib", "variant_harmonization.R"))

tmp <- tempfile("variant_harmonization_test")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)


write_pvar <- function(prefix, rows) {
  lines <- c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tINFO",
    apply(rows, 1, paste, collapse = "\t")
  )
  writeLines(lines, paste0(prefix, ".pvar"))
  writeLines("dummy-pgen", paste0(prefix, ".pgen"))
  writeLines(c("#FID\tIID", "F1\tI1"), paste0(prefix, ".psam"))
}


artifact_paths <- function(directory) {
  list(
    shared = file.path(directory, "shared_variants.txt"),
    mismatch = file.path(directory, "shared_variant_mismatches.tsv"),
    mapping = file.path(directory, "variant_harmonization.tsv"),
    reference_extract = file.path(directory, "reference_native_variants.txt"),
    study_extract = file.path(directory, "study_native_variants.txt"),
    reference_update = file.path(directory, "reference_update_names.tsv"),
    study_update = file.path(directory, "study_update_names.tsv")
  )
}


run_harmonization <- function(reference_prefix, study_prefix, paths) {
  harmonize_pvar_variants(
    reference_prefix,
    study_prefix,
    paths$shared,
    paths$mismatch,
    paths$mapping,
    paths$reference_extract,
    paths$study_extract,
    paths$reference_update,
    paths$study_update,
    exclude_palindromic = TRUE
  )
}


run_helper <- function(reference_prefix, study_prefix, output_dir,
                       exclude_palindromic = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  log <- file.path(output_dir, "helper.log")
  helper <- file.path(script_dir, "harmonize_pvar_variants.sh")
  args <- c(
    helper,
    "--reference-pvar", paste0(reference_prefix, ".pvar"),
    "--study-pvar", paste0(study_prefix, ".pvar"),
    "--output-dir", output_dir,
    "--exclude-palindromic", if (exclude_palindromic) "true" else "false"
  )
  status <- system2(
    "bash", shQuote(args, type = "sh"),
    stdout = log, stderr = log
  )
  list(status = status, log = readLines(log, warn = FALSE))
}


harmonization_error <- function(reference_prefix, study_prefix, stem) {
  tryCatch({
    run_harmonization(
      reference_prefix,
      study_prefix,
      artifact_paths(file.path(tmp, stem))
    )
    ""
  }, error = function(err) conditionMessage(err))
}


# The combined fixture covers all supported ID combinations and exclusions.
reference_prefix <- file.path(tmp, "reference")
study_prefix <- file.path(tmp, "study")
reference_rows <- rbind(
  c("1", 100, "rs100", "A", "G", "."),
  c("1", 200, "rs200", "A", "C", "."),
  c("1", 300, "1:300:C:T", "C", "T", "."),
  c("1", 400, "1:400:G:T", "G", "T", "."),
  c("1", 500, "rs500", "A", "G", "."),
  c("1", 600, "rs600", "A", "C", "."),
  c("1", 700, "pal_ref", "A", "T", "."),
  c("1", 800, "mismatch_ref", "A", "G", "."),
  c("1", 900, "dup_ref_a", "A", "G", "."),
  c("1", 900, "dup_ref_b", "A", "G", "."),
  c("1", 1000, "provisional_ref", "A", "G", "PR"),
  c("1", 1100, "known_ref", "A", "G", "."),
  c("1", 1200, "study_cpra_ref", "A", "C", "PR")
)
study_rows <- rbind(
  c("chr01", 100, "rs100", "A", "G", "."),
  c("Chr1", 200, "1:200:A:C", "A", "C", "."),
  c("CHR1", 300, "rs300", "C", "T", "."),
  c("1", 400, "1:400:G:T", "G", "T", "."),
  c("1", 500, "rs501", "A", "G", "."),
  c("1", 600, "1:600:C:A", "C", "A", "PR"),
  c("1", 700, "pal_study", "T", "A", "PR"),
  c("1", 800, "mismatch_study", "A", "C", "."),
  c("1", 900, "dup_study", "A", "G", "."),
  c("1", 1000, "provisional_study", "A", "G", "PR"),
  c("1", 1100, "known_study", "G", "A", "."),
  c("1", 1200, "study_cpra_study", "C", "A", ".")
)
write_pvar(reference_prefix, reference_rows)
write_pvar(study_prefix, study_rows)
reference_sha <- sha256_file(paste0(reference_prefix, ".pvar"))
study_sha <- sha256_file(paste0(study_prefix, ".pvar"))

paths <- artifact_paths(file.path(tmp, "harmonized"))
result <- run_harmonization(reference_prefix, study_prefix, paths)
stopifnot(result$retained == 7L)
stopifnot(result$mismatches == 5L)

mapping <- read_tsv(paths$mapping)
stopifnot(identical(mapping$harmonized_id, c(
  "rs100", "rs200", "rs300", "1:400:G:T", "rs500", "rs600", "1:1200:C:A"
)))
stopifnot(identical(mapping$harmonized_id_type, c(
  "rsID", "rsID", "rsID", "CPRA", "rsID", "rsID", "CPRA"
)))
stopifnot(identical(mapping$harmonized_id_source, c(
  "reference_rsid",
  "reference_rsid",
  "study_rsid",
  "reference_pvar_known_ref",
  "reference_rsid",
  "reference_rsid",
  "study_pvar_known_ref"
)))
stopifnot(mapping$differing_rsid[mapping$harmonized_id == "rs500"] == "True")
stopifnot(mapping$orientation[mapping$harmonized_id == "rs600"] == "ref_alt_swap")
stopifnot(mapping$orientation[mapping$harmonized_id == "1:1200:C:A"] == "ref_alt_swap")

mismatches <- read_tsv(paths$mismatch)
stopifnot(setequal(mismatches$reason, c(
  "palindromic_snp_excluded",
  "allele_mismatch",
  "duplicate_locus_allele_key_reference",
  "no_pvar_known_ref_for_cpra",
  "pvar_known_ref_disagreement"
)))
stopifnot(all(c("reference_native_id", "study_native_id", "excluded_side") %in% names(mismatches)))
stopifnot(identical(readLines(paths$shared), mapping$harmonized_id))
stopifnot(identical(readLines(paths$reference_extract), mapping$reference_native_id))
stopifnot(identical(readLines(paths$study_extract), mapping$study_native_id))
stopifnot(identical(
  readLines(paths$reference_update),
  paste(mapping$reference_native_id, mapping$harmonized_id, sep = "\t")
))
stopifnot(identical(
  readLines(paths$study_update),
  paste(mapping$study_native_id, mapping$harmonized_id, sep = "\t")
))
stopifnot(sha256_file(paste0(reference_prefix, ".pvar")) == reference_sha)
stopifnot(sha256_file(paste0(study_prefix, ".pvar")) == study_sha)


# A tiny fake PLINK2 applies --extract and --update-name to text PVAR fixtures.
fake_plink2 <- file.path(tmp, "fake_plink2.R")
writeLines(c(
  "#!/usr/bin/env Rscript",
  "args <- commandArgs(trailingOnly = TRUE)",
  "value_after <- function(flag) {",
  "  index <- match(flag, args)",
  "  if (is.na(index) || index == length(args)) return('')",
  "  args[[index + 1L]]",
  "}",
  "input <- value_after('--pfile')",
  "output <- value_after('--out')",
  "extract <- value_after('--extract')",
  "update <- value_after('--update-name')",
  "dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)",
  "lines <- readLines(paste0(input, '.pvar'), warn = FALSE)",
  "header_index <- which(!startsWith(lines, '##'))[[1]]",
  "header <- strsplit(lines[[header_index]], '\\t', fixed = TRUE)[[1]]",
  "id_column <- match('ID', header)",
  "metadata <- if (header_index > 1L) lines[seq_len(header_index - 1L)] else character()",
  "body <- lines[(header_index + 1L):length(lines)]",
  "fields <- strsplit(body, '\\t', fixed = TRUE)",
  "if (nzchar(extract)) {",
  "  keep <- readLines(extract, warn = FALSE)",
  "  fields <- fields[vapply(fields, function(x) x[[id_column]] %in% keep, logical(1))]",
  "}",
  "if (nzchar(update)) {",
  "  names <- read.table(update, sep = '\\t', stringsAsFactors = FALSE, quote = '')",
  "  rename <- setNames(names[[2]], names[[1]])",
  "  fields <- lapply(fields, function(x) { x[[id_column]] <- rename[[x[[id_column]]]]; x })",
  "}",
  "body <- vapply(fields, paste, character(1), collapse = '\\t')",
  "writeLines(c(metadata, paste(header, collapse = '\\t'), body), paste0(output, '.pvar'))",
  "invisible(file.copy(paste0(input, '.pgen'), paste0(output, '.pgen'), overwrite = TRUE))",
  "invisible(file.copy(paste0(input, '.psam'), paste0(output, '.psam'), overwrite = TRUE))"
), fake_plink2)
Sys.chmod(fake_plink2, mode = "0755")
fake_config <- list(tools = list(plink2 = fake_plink2))

reference_renamed <- file.path(tmp, "reference_renamed")
study_renamed <- file.path(tmp, "study_renamed")
extract_and_rename_pgen(
  fake_config, reference_prefix, paths$reference_extract,
  paths$reference_update, reference_renamed, "1"
)
extract_and_rename_pgen(
  fake_config, study_prefix, paths$study_extract,
  paths$study_update, study_renamed, "1"
)
validation <- file.path(tmp, "harmonized_variants.ok")
validated_count <- validate_harmonized_pgen_variants(
  reference_renamed, study_renamed, validation
)
stopifnot(validated_count == result$retained)
stopifnot(grepl(paste0("\t", result$retained, "$"), readLines(validation)))
stopifnot(sha256_file(paste0(reference_prefix, ".pvar")) == reference_sha)
stopifnot(sha256_file(paste0(study_prefix, ".pvar")) == study_sha)

# Validation must reject even one ID/locus/allele-set difference.
bad_study <- file.path(tmp, "study_renamed_bad")
for (suffix in c(".pgen", ".pvar", ".psam")) {
  invisible(file.copy(paste0(study_renamed, suffix), paste0(bad_study, suffix)))
}
bad_lines <- readLines(paste0(bad_study, ".pvar"))
bad_header_index <- which(!startsWith(bad_lines, "##"))[[1]]
bad_header <- strsplit(bad_lines[[bad_header_index]], "\t", fixed = TRUE)[[1]]
bad_alt_column <- match("ALT", bad_header)
bad_first_data <- bad_header_index + 1L
bad_fields <- strsplit(bad_lines[[bad_first_data]], "\t", fixed = TRUE)[[1]]
bad_fields[[bad_alt_column]] <- "T"
bad_lines[[bad_first_data]] <- paste(bad_fields, collapse = "\t")
writeLines(bad_lines, paste0(bad_study, ".pvar"))
validation_error <- tryCatch({
  validate_harmonized_pgen_variants(
    reference_renamed, bad_study, file.path(tmp, "bad_harmonized_variants.ok")
  )
  ""
}, error = function(err) conditionMessage(err))
stopifnot(grepl("do not have identical ID/locus/allele sets", validation_error, fixed = TRUE))


# A selected operational ID at two loci is a hard integrity failure.
collision_reference <- file.path(tmp, "collision_reference")
collision_study <- file.path(tmp, "collision_study")
write_pvar(collision_reference, rbind(
  c("1", 2000, "rs999", "A", "G", "."),
  c("1", 2100, "reference_cpra", "A", "C", ".")
))
write_pvar(collision_study, rbind(
  c("1", 2000, "study_cpra", "A", "G", "."),
  c("1", 2100, "rs999", "A", "C", ".")
))
collision_error <- tryCatch({
  run_harmonization(
    collision_reference,
    collision_study,
    artifact_paths(file.path(tmp, "collision"))
  )
  ""
}, error = function(err) conditionMessage(err))
stopifnot(grepl("maps to multiple loci", collision_error, fixed = TRUE))


# Reordered inputs must take the fallback path without changing any artifact.
unsorted_reference <- file.path(tmp, "unsorted_reference")
unsorted_study <- file.path(tmp, "unsorted_study")
write_pvar(unsorted_reference, reference_rows[rev(seq_len(nrow(reference_rows))), , drop = FALSE])
write_pvar(unsorted_study, study_rows[rev(seq_len(nrow(study_rows))), , drop = FALSE])
unsorted_output <- file.path(tmp, "unsorted_helper")
unsorted_result <- run_helper(unsorted_reference, unsorted_study, unsorted_output)
stopifnot(identical(unsorted_result$status, 0L))
unsorted_summary <- read.delim(
  file.path(unsorted_output, "harmonization_summary.tsv"),
  sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
)
stopifnot(unsorted_summary$reference_sort_mode == "fallback")
stopifnot(unsorted_summary$study_sort_mode == "fallback")
stopifnot(any(grepl("external-sort fallback", unsorted_result$log, fixed = TRUE)))

expected_artifacts <- c(
  mapping = paths$mapping,
  mismatch = paths$mismatch,
  shared = paths$shared,
  reference_extract = paths$reference_extract,
  study_extract = paths$study_extract,
  reference_update = paths$reference_update,
  study_update = paths$study_update
)
fallback_artifacts <- c(
  mapping = file.path(unsorted_output, "variant_harmonization.tsv"),
  mismatch = file.path(unsorted_output, "shared_variant_mismatches.tsv"),
  shared = file.path(unsorted_output, "shared_variants.txt"),
  reference_extract = file.path(unsorted_output, "reference_native_variants.txt"),
  study_extract = file.path(unsorted_output, "study_native_variants.txt"),
  reference_update = file.path(unsorted_output, "reference_update_names.tsv"),
  study_update = file.path(unsorted_output, "study_update_names.tsv")
)
for (name in names(expected_artifacts)) {
  stopifnot(identical(
    readLines(expected_artifacts[[name]], warn = FALSE),
    readLines(fallback_artifacts[[name]], warn = FALSE)
  ))
}


# The parser accepts metadata, reordered columns, and an absent INFO column.
parser_reference <- file.path(tmp, "parser_reference")
parser_study <- file.path(tmp, "parser_study")
writeLines(c(
  "##fileformat=VCFv4.2",
  "##source=reordered-test",
  "#ID\tALT\tPOS\tCHROM\tREF\tINFO",
  "reference_native\tG\t3000\tchr01\tA\t."
), paste0(parser_reference, ".pvar"))
writeLines(c(
  "##fileformat=VCFv4.2",
  "#ALT\tID\tREF\tPOS\tCHROM",
  "G\tstudy_native\tA\t3000\t1"
), paste0(parser_study, ".pvar"))
parser_paths <- artifact_paths(file.path(tmp, "parser"))
parser_result <- run_harmonization(parser_reference, parser_study, parser_paths)
stopifnot(parser_result$retained == 1L)
parser_mapping <- read_tsv(parser_paths$mapping)
stopifnot(parser_mapping$harmonized_id == "1:3000:A:G")
stopifnot(parser_mapping$harmonized_id_source == "reference_pvar_known_ref")


# Malformed PVAR structures fail before any harmonized artifacts are published.
valid_parser_study <- parser_study

missing_header_prefix <- file.path(tmp, "missing_header")
writeLines(c(
  "##fileformat=VCFv4.2",
  "#CHROM\tPOS\tID\tREF",
  "1\t3100\trs3100\tA"
), paste0(missing_header_prefix, ".pvar"))
missing_header_error <- harmonization_error(
  missing_header_prefix, valid_parser_study, "missing_header_error"
)
stopifnot(grepl("header is missing required", missing_header_error, fixed = TRUE))

short_row_prefix <- file.path(tmp, "short_row")
writeLines(c(
  "#CHROM\tPOS\tID\tREF\tALT\tINFO",
  "1\t3100\trs3100\tA\tG"
), paste0(short_row_prefix, ".pvar"))
short_row_error <- harmonization_error(
  short_row_prefix, valid_parser_study, "short_row_error"
)
stopifnot(grepl("fewer columns", short_row_error, fixed = TRUE))

invalid_position_prefix <- file.path(tmp, "invalid_position")
writeLines(c(
  "#CHROM\tPOS\tID\tREF\tALT",
  "1\t3100.5\trs3100\tA\tG"
), paste0(invalid_position_prefix, ".pvar"))
invalid_position_error <- harmonization_error(
  invalid_position_prefix, valid_parser_study, "invalid_position_error"
)
stopifnot(grepl("non-integer or nonpositive POS", invalid_position_error, fixed = TRUE))

empty_eligible_prefix <- file.path(tmp, "empty_eligible")
writeLines(c(
  "#CHROM\tPOS\tID\tREF\tALT",
  "X\t3100\trs3100\tA\tG",
  "1\t3200\t.\tA\tG"
), paste0(empty_eligible_prefix, ".pvar"))
empty_eligible_error <- harmonization_error(
  empty_eligible_prefix, valid_parser_study, "empty_eligible_error"
)
stopifnot(grepl("contains no eligible", empty_eligible_error, fixed = TRUE))


# Native IDs must remain globally unique even when the loci differ.
duplicate_id_reference <- file.path(tmp, "duplicate_id_reference")
duplicate_id_study <- file.path(tmp, "duplicate_id_study")
write_pvar(duplicate_id_reference, rbind(
  c("1", 4000, "duplicate_native", "A", "G", "."),
  c("1", 4100, "duplicate_native", "A", "C", ".")
))
write_pvar(duplicate_id_study, rbind(
  c("1", 4000, "study_4000", "A", "G", "."),
  c("1", 4100, "study_4100", "A", "C", ".")
))
duplicate_id_error <- harmonization_error(
  duplicate_id_reference, duplicate_id_study, "duplicate_id_error"
)
stopifnot(grepl("duplicate native variant ID", duplicate_id_error, fixed = TRUE))


# Duplicate locus/allele keys are reported even when that locus is unshared.
unshared_duplicate_reference <- file.path(tmp, "unshared_duplicate_reference")
unshared_duplicate_study <- file.path(tmp, "unshared_duplicate_study")
write_pvar(unshared_duplicate_reference, rbind(
  c("1", 4150, "reference_dup_a", "A", "G", "."),
  c("1", 4150, "reference_dup_b", "A", "G", "."),
  c("1", 4175, "rs4175", "A", "C", ".")
))
write_pvar(unshared_duplicate_study, rbind(
  c("1", 4160, "study_dup_a", "A", "G", "."),
  c("1", 4160, "study_dup_b", "A", "G", "."),
  c("1", 4175, "study_4175", "A", "C", ".")
))
unshared_duplicate_paths <- artifact_paths(file.path(tmp, "unshared_duplicate"))
unshared_duplicate_result <- run_harmonization(
  unshared_duplicate_reference,
  unshared_duplicate_study,
  unshared_duplicate_paths
)
stopifnot(unshared_duplicate_result$retained == 1L)
unshared_duplicate_rows <- read_tsv(unshared_duplicate_paths$mismatch)
stopifnot(setequal(unshared_duplicate_rows$reason, c(
  "duplicate_locus_allele_key_reference",
  "duplicate_locus_allele_key_study"
)))


# A repeated selected ID at the same locus is also a hard integrity failure.
same_locus_reference <- file.path(tmp, "same_locus_reference")
same_locus_study <- file.path(tmp, "same_locus_study")
write_pvar(same_locus_reference, rbind(
  c("1", 4200, "rs4200", "A", "G", "."),
  c("1", 4200, "reference_ct", "C", "T", ".")
))
write_pvar(same_locus_study, rbind(
  c("1", 4200, "study_ag", "A", "G", "."),
  c("1", 4200, "rs4200", "C", "T", ".")
))
same_locus_error <- harmonization_error(
  same_locus_reference, same_locus_study, "same_locus_error"
)
stopifnot(grepl("multiple variants at locus", same_locus_error, fixed = TRUE))


# The explicit locus ceiling prevents pathological metadata from growing memory.
large_locus_reference <- file.path(tmp, "large_locus_reference")
large_locus_study <- file.path(tmp, "large_locus_study")
large_locus_rows <- c(
  "#CHROM\tPOS\tID\tREF\tALT",
  paste("1", "5000", paste0("large_", seq_len(10001)), "A", "G", sep = "\t")
)
writeLines(large_locus_rows, paste0(large_locus_reference, ".pvar"))
writeLines(c(
  "#CHROM\tPOS\tID\tREF\tALT",
  "1\t5000\tlarge_study\tA\tG"
), paste0(large_locus_study, ".pvar"))
large_locus_error <- harmonization_error(
  large_locus_reference, large_locus_study, "large_locus_error"
)
stopifnot(grepl("more than 10000 eligible records", large_locus_error, fixed = TRUE))


# Palindromic variants remain available only when the caller explicitly allows them.
pal_reference <- file.path(tmp, "pal_reference")
pal_study <- file.path(tmp, "pal_study")
write_pvar(pal_reference, rbind(c("1", 6000, "rs6000", "A", "T", ".")))
write_pvar(pal_study, rbind(c("1", 6000, "study_pal", "A", "T", ".")))
pal_paths <- artifact_paths(file.path(tmp, "pal_allowed"))
pal_result <- harmonize_pvar_variants(
  pal_reference, pal_study,
  pal_paths$shared, pal_paths$mismatch, pal_paths$mapping,
  pal_paths$reference_extract, pal_paths$study_extract,
  pal_paths$reference_update, pal_paths$study_update,
  exclude_palindromic = FALSE
)
stopifnot(pal_result$retained == 1L)
stopifnot(identical(readLines(pal_paths$shared), "rs6000"))

cat("Variant harmonization tests passed\n")
