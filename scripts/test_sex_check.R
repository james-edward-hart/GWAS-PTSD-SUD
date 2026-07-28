#!/usr/bin/env Rscript

# Focused real-PLINK integration checks for the compact sex-marker workflow.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))

repo <- normalizePath(file.path(script_dir, ".."))
local_plink <- file.path(repo, "software", "local", "plink2")
plink <- ""
for (candidate in unique(c(local_plink, Sys.which("plink2")))) {
  if (!nzchar(candidate) || !file.exists(candidate) || file.access(candidate, 1) != 0) next
  status <- suppressWarnings(system2(candidate, "--version", stdout = FALSE, stderr = FALSE))
  if (identical(status, 0L)) {
    plink <- candidate
    break
  }
}
if (!nzchar(plink)) {
  cat("Sex-check integration test skipped: PLINK2 is unavailable\n")
  quit(status = 0)
}
if (!nzchar(Sys.getenv("GAWK_BIN")) && !nzchar(Sys.which("gawk"))) {
  # The production environment supplies GNU awk. BSD awk is sufficient for
  # these small parser fixtures when running the macOS development test.
  Sys.setenv(GAWK_BIN = Sys.which("awk"))
}
if (!requireNamespace("yaml", quietly = TRUE)) {
  cat("Sex-check integration test skipped: R package 'yaml' is unavailable\n")
  quit(status = 0)
}

tmp <- tempfile("sex_check_test")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

source_prefix <- file.path(repo, "data", "example", "hapmap3")
manifest <- read_tsv(file.path(repo, "config", "example_sample_manifest.tsv"))
manifest <- manifest[seq_len(min(80L, nrow(manifest))), , drop = FALSE]
manifest_path <- file.path(tmp, "manifest.tsv")
write_tsv(manifest, manifest_path)


make_pgen_fixture <- function(prefix, filter_args) {
  run_command(plink, c(
    "--bfile", source_prefix, filter_args,
    "--make-pgen", "--sort-vars", "--out", prefix
  ))
}


metric_value <- function(summary, metric) {
  value <- summary$value[match(metric, summary$metric)]
  if (!length(value) || is.na(value)) die("test summary is missing metric: ", metric)
  as.character(value)
}


run_case <- function(name, block, build = "GRCh37", allow_no_markers = FALSE,
                     snps_only_acgt = FALSE, expect_skip = FALSE,
                     expected_error = "") {
  case_dir <- file.path(tmp, name)
  dir.create(case_dir)
  prefix <- file.path(case_dir, "plink_sex_check")
  config_path <- file.path(case_dir, "config.yaml")
  config <- list(
    project = list(inferred_genome_build = build),
    inputs = list(sample_manifest = manifest_path),
    tools = list(plink2 = plink),
    genotypes = block,
    qc = list(
      maf_min = 0.01,
      geno_missing_max = 0.05,
      snps_only_acgt = snps_only_acgt
    ),
    relatedness = list(ld_prune = list(window = "500kb", step = 1, r2 = 0.2)),
    sex_check = list(
      enabled = TRUE,
      action = "warn",
      allow_no_sex_markers = allow_no_markers,
      max_female_xf = "",
      min_male_xf = "",
      max_female_yrate = "",
      min_male_yrate = ""
    )
  )
  yaml::write_yaml(config, config_path)

  outputs <- list(
    detail = file.path(case_dir, "sexcheck.tsv"),
    remove = file.path(case_dir, "remove.tsv"),
    keep = file.path(case_dir, "keep.tsv"),
    summary = file.path(case_dir, "summary.tsv"),
    prune_in = paste0(prefix, ".sex_marker_prune.prune.in"),
    prune_out = paste0(prefix, ".sex_marker_prune.prune.out")
  )
  source_paths <- unname(genotype_component_paths(block))
  source_hashes <- vapply(source_paths, sha256_file, character(1))
  log <- file.path(case_dir, "run.log")
  error_log <- file.path(case_dir, "run.stderr.log")
  status <- system2(
    "Rscript",
    vapply(c(
      file.path(script_dir, "sex_check.R"),
      "--config", config_path,
      "--plink-out-prefix", prefix,
      "--sexcheck-out", outputs$detail,
      "--remove-out", outputs$remove,
      "--keep-out", outputs$keep,
      "--summary-out", outputs$summary,
      "--threads", "2"
    ), shQuote, character(1), type = "sh"),
    stdout = log,
    stderr = error_log
  )
  if (nzchar(expected_error)) {
    stopifnot(!identical(status, 0L))
    log_paths <- unique(c(
      log,
      error_log,
      list.files(case_dir, pattern = "\\.log$", full.names = TRUE)
    ))
    failure_log <- unlist(lapply(log_paths[file.exists(log_paths)], readLines, warn = FALSE))
    if (!any(grepl(expected_error, failure_log, fixed = TRUE))) {
      die("expected failure text '", expected_error, "' was absent in case ", name,
        "\n", paste(failure_log, collapse = "\n"))
    }
    stopifnot(identical(vapply(source_paths, sha256_file, character(1)), source_hashes))
    stopifnot(!file.exists(paste0(prefix, ".sex_markers.pgen")))
    return(list(prefix = prefix, log = log))
  }
  if (!identical(status, 0L)) {
    die("sex-check integration case failed: ", name, "\n",
      paste(c(readLines(log, warn = FALSE), readLines(error_log, warn = FALSE)), collapse = "\n"))
  }
  stopifnot(identical(vapply(source_paths, sha256_file, character(1)), source_hashes))

  summary <- read_tsv(outputs$summary)
  stopifnot(nrow(read_tsv(outputs$keep)) == nrow(manifest))
  stopifnot(file.exists(outputs$prune_in), file.exists(outputs$prune_out))
  if (expect_skip) {
    stopifnot(metric_value(summary, "status") == "skipped")
    stopifnot(metric_value(summary, "skipped_reason") == "no_sex_chromosome_markers")
    stopifnot(file.info(outputs$prune_in)$size == 0)
  } else {
    stopifnot(metric_value(summary, "status") == "completed")
    stopifnot(nrow(read_tsv(outputs$detail)) == nrow(manifest))
    stopifnot(as.numeric(metric_value(summary, "post_qc_x_variants")) > 0)
    stopifnot(as.numeric(metric_value(summary, "pruned_variants")) > 0)
    stopifnot(file.info(outputs$prune_in)$size > 0)
    stopifnot(!any(file.exists(paste0(prefix, ".sex_markers", c(".pgen", ".pvar", ".psam")))))
  }
  list(prefix = prefix, outputs = outputs, summary = summary, log = log, error_log = error_log)
}


# Main case: an unsplit BED takes the build-matched split path, filters only
# manifest samples, preserves numeric alleles, and retains non-PAR X markers.
main <- run_case("unsplit_b37", list(type = "bed", prefix = source_prefix))
stopifnot(metric_value(main$summary, "par_handling") == "split_b37")
stopifnot(as.numeric(metric_value(main$summary, "source_x_variants")) > 0)
stopifnot(as.numeric(metric_value(main$summary, "source_y_variants")) == 0)
main_run_log <- c(
  readLines(main$log, warn = FALSE),
  readLines(main$error_log, warn = FALSE)
)
stopifnot(sum(grepl("PLINK metadata scan: 13928 rows", main_run_log, fixed = TRUE)) == 1)

compact_log <- paste(readLines(paste0(main$prefix, ".sex_markers.log"), warn = FALSE), collapse = "\n")
prune_log <- paste(readLines(paste0(main$prefix, ".sex_marker_prune.log"), warn = FALSE), collapse = "\n")
check_log <- paste(readLines(paste0(main$prefix, ".log"), warn = FALSE), collapse = "\n")
for (flag in c("--chr X,Y,XY,PAR1,PAR2", "--split-par b37", "--maf 0.01",
               "--geno 0.05", "--max-alleles 2", "--rm-dup exclude-all",
               "--sort-vars")) {
  stopifnot(grepl(flag, compact_log, fixed = TRUE))
}
stopifnot(!grepl("--mind", compact_log, fixed = TRUE))
stopifnot(!grepl("--hwe", compact_log, fixed = TRUE))
stopifnot(!grepl("--snps-only", compact_log, fixed = TRUE))
stopifnot(grepl("--chr X,Y", prune_log, fixed = TRUE))
stopifnot(grepl("--indep-pairwise 500kb 1 0.2", prune_log, fixed = TRUE))
stopifnot(grepl("--extract", check_log, fixed = TRUE))
stopifnot(grepl("--check-sex", check_log, fixed = TRUE))

pruned_ids <- readLines(main$outputs$prune_in, warn = FALSE)
bim <- read.table(paste0(source_prefix, ".bim"), stringsAsFactors = FALSE)
names(bim)[c(1, 2, 4)] <- c("chrom", "id", "pos")
pruned_rows <- bim[match(pruned_ids, bim$id), , drop = FALSE]
stopifnot(!anyNA(pruned_rows$id))
stopifnot(all(pruned_rows$chrom %in% c(23, 24)))
stopifnot(all(
  pruned_rows$chrom != 23 |
    (pruned_rows$pos > 2699520 & pruned_rows$pos < 154931044)
))


# Build compact source fixtures once with real PLINK2.
x_prefix <- file.path(tmp, "source_x")
make_pgen_fixture(x_prefix, c("--chr", "X"))
x_only <- run_case("x_only", list(type = "pgen", prefix = x_prefix))
stopifnot(as.numeric(metric_value(x_only$summary, "source_y_variants")) == 0)

autosome_prefix <- file.path(tmp, "source_autosome")
make_pgen_fixture(autosome_prefix, "--autosome")
no_sex <- run_case(
  "no_sex_override",
  list(type = "pgen", prefix = autosome_prefix),
  allow_no_markers = TRUE,
  expect_skip = TRUE
)
stopifnot(metric_value(no_sex$summary, "par_handling") == "not_run_no_usable_x")
invisible(run_case(
  "no_sex_strict_failure",
  list(type = "pgen", prefix = autosome_prefix),
  expected_error = "requires usable non-PAR X markers"
))

# The source uses numeric allele codes, so enabling A/C/G/T-only filtering
# intentionally removes every marker and must fail instead of skipping QC.
invisible(run_case(
  "empty_after_qc",
  list(type = "pgen", prefix = x_prefix),
  snps_only_acgt = TRUE,
  expected_error = "--snps-only just-acgt"
))


# A PGEN with explicit PAR1 and Y labels must not be split again.
sex_prefix <- file.path(tmp, "source_explicit_par")
make_pgen_fixture(sex_prefix, c("--chr", "X,XY"))
pvar_lines <- readLines(paste0(sex_prefix, ".pvar"), warn = FALSE)
header_index <- which(grepl("^#CHROM\\t", pvar_lines))[[1]]
data_index <- seq.int(header_index + 1L, length(pvar_lines))
fields <- strsplit(pvar_lines[data_index], "\t", fixed = TRUE)
x_index <- which(vapply(fields, function(row) row[[1]] == "X", logical(1)))[[1]]
xy_index <- which(vapply(fields, function(row) row[[1]] == "XY", logical(1)))[[1]]
fields[[x_index]][1:2] <- c("PAR1", "100000")
fields[[xy_index]][[1]] <- "Y"
pvar_lines[data_index] <- vapply(fields, paste, collapse = "\t", character(1))
writeLines(pvar_lines, paste0(sex_prefix, ".pvar"))

explicit_par <- run_case("explicit_par", list(type = "pgen", prefix = sex_prefix))
stopifnot(metric_value(explicit_par$summary, "par_handling") == "already_split")
stopifnot(as.numeric(metric_value(explicit_par$summary, "source_y_variants")) > 0)
stopifnot(as.numeric(metric_value(explicit_par$summary, "source_par1_variants")) > 0)
explicit_log <- paste(
  readLines(paste0(explicit_par$prefix, ".sex_markers.log"), warn = FALSE),
  collapse = "\n"
)
stopifnot(!grepl("--split-par", explicit_log, fixed = TRUE))


# The same X-only fixture exercises the GRCh38 PAR boundary selection.
b38 <- run_case("unsplit_b38", list(type = "pgen", prefix = x_prefix), build = "GRCh38")
stopifnot(metric_value(b38$summary, "par_handling") == "split_b38")
b38_log <- paste(readLines(paste0(b38$prefix, ".sex_markers.log"), warn = FALSE), collapse = "\n")
stopifnot(grepl("--split-par b38", b38_log, fixed = TRUE))

cat("Sex-check integration tests passed\n")
