#!/usr/bin/env Rscript

# Focused tests for Phase 2 regenie helper behavior that does not require
# external PLINK2 or regenie execution.

cmd <- commandArgs(FALSE)
script_path <- sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])
repo <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(repo, "scripts", "lib", "stage1.R"))

if (!identical(qc_info_min(list(qc = list())), 0.9)) stop("shared INFO/R2 default is not 0.9")

tmp <- tempfile("phase2-reg-")
dir.create(tmp, recursive = TRUE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

write_lines <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path)
}

run_info_helper <- function(pvar, stem, expect_success = TRUE) {
  paths <- list(
    pass = paste0(stem, ".pass.snplist"),
    excluded = paste0(stem, ".excluded.tsv"),
    summary = paste0(stem, ".summary.tsv")
  )
  command <- c(
    file.path(repo, "scripts", "filter_pvar_info.sh"),
    "--pvar", pvar,
    "--threshold", "0.8",
    "--pass-ids", paths$pass,
    "--excluded-report", paths$excluded,
    "--summary", paths$summary
  )
  output <- suppressWarnings(system2("bash", command, stdout = TRUE, stderr = TRUE))
  code <- attr(output, "status") %||% 0L
  if (expect_success && !identical(as.integer(code), 0L)) {
    stop("filter_pvar_info.sh failed:\n", paste(output, collapse = "\n"), call. = FALSE)
  }
  if (!expect_success && identical(as.integer(code), 0L)) {
    stop("filter_pvar_info.sh unexpectedly succeeded", call. = FALSE)
  }
  c(paths, list(output = output))
}

run_phase2 <- function(args, expect_success = TRUE) {
  status <- system2("Rscript", c(file.path(repo, "scripts", "phase2_regenie.R"), args), stdout = TRUE, stderr = TRUE)
  code <- attr(status, "status") %||% 0L
  if (expect_success && !identical(as.integer(code), 0L)) {
    stop("phase2_regenie.R failed:\n", paste(status, collapse = "\n"), call. = FALSE)
  }
  if (!expect_success && identical(as.integer(code), 0L)) {
    stop("phase2_regenie.R unexpectedly succeeded", call. = FALSE)
  }
  status
}

# Reordered columns, metadata lines, direct-field precedence, and missing-value
# retention must match the former R implementation.
helper_direct_pvar <- file.path(tmp, "helper_direct.pvar")
write_lines(c(
  "##fileformat=VCFv4.3",
  "##INFO=<ID=R2,Number=1,Type=Float,Description=\"quality\">",
  "ALT\tID\tPOS\tINFO\tREF\t#CHROM\tR2",
  "G\trs_low\t100\tR2=0.99\tA\t1\t0.50",
  "T\trs_missing\t101\tR2=0.10\tC\t1\t.",
  "C\trs_high\t102\t.\tA\t1\t0.95"
), helper_direct_pvar)
helper_direct <- run_info_helper(helper_direct_pvar, file.path(tmp, "helper_direct"))
stopifnot(identical(readLines(helper_direct$pass), c("rs_missing", "rs_high")))
helper_direct_excluded <- read_tsv(helper_direct$excluded)
stopifnot(
  identical(helper_direct_excluded$variant_id, "rs_low"),
  identical(helper_direct_excluded$info_metric, "R2"),
  identical(helper_direct_excluded$exclusion_reason, "info_r2_below_min")
)

helper_keys_pvar <- file.path(tmp, "helper_keys.pvar")
write_lines(c(
  "#CHROM\tPOS\tID\tREF\tALT\tR2\tINFO",
  "1\t100\trs_key_low\tA\tG\t.\tPR;R2=0.50",
  "1\t101\trs_key_high\tC\tT\t.\tINFO=0.90",
  "1\t102\trs_key_missing\tA\tC\t.\tPR"
), helper_keys_pvar)
helper_keys <- run_info_helper(helper_keys_pvar, file.path(tmp, "helper_keys"))
stopifnot(identical(readLines(helper_keys$pass), c("rs_key_high", "rs_key_missing")))
helper_keys_excluded <- read_tsv(helper_keys$excluded)
stopifnot(
  identical(helper_keys_excluded$variant_id, "rs_key_low"),
  identical(helper_keys_excluded$info_metric, "INFO:INFO/R2")
)

# A PVAR with no quality annotation is scanned once and leaves no stale filter
# artifacts. Malformed short rows remain hard failures.
helper_none_pvar <- file.path(tmp, "helper_none.pvar")
write_lines(c("#CHROM\tPOS\tID\tREF\tALT", "1\t100\trs1\tA\tG"), helper_none_pvar)
helper_none_stem <- file.path(tmp, "helper_none")
write_lines("stale", paste0(helper_none_stem, ".pass.snplist"))
write_lines("stale", paste0(helper_none_stem, ".excluded.tsv"))
helper_none <- run_info_helper(helper_none_pvar, helper_none_stem)
stopifnot(!file.exists(helper_none$pass), !file.exists(helper_none$excluded))

helper_no_low_pvar <- file.path(tmp, "helper_no_low.pvar")
write_lines(c("#CHROM\tPOS\tID\tREF\tALT\tR2", "1\t100\trs1\tA\tG\t0.9"), helper_no_low_pvar)
helper_no_low <- run_info_helper(helper_no_low_pvar, file.path(tmp, "helper_no_low"))
helper_no_low_summary <- read_plink_metadata_summary(helper_no_low$summary)
stopifnot(
  identical(helper_no_low_summary$below_min, "0"),
  !file.exists(helper_no_low$pass),
  !file.exists(helper_no_low$excluded)
)

helper_short_pvar <- file.path(tmp, "helper_short.pvar")
write_lines(c("#CHROM\tPOS\tID\tREF\tALT\tR2", "1\t100\trs1\tA\tG"), helper_short_pvar)
invisible(run_info_helper(helper_short_pvar, file.path(tmp, "helper_short"), expect_success = FALSE))

helper_all_low_pvar <- file.path(tmp, "helper_all_low.pvar")
write_lines(c("#CHROM\tPOS\tID\tREF\tALT\tR2", "1\t100\trs1\tA\tG\t0.1"), helper_all_low_pvar)
invisible(run_info_helper(helper_all_low_pvar, file.path(tmp, "helper_all_low"), expect_success = FALSE))

# A larger fixture exercises both sequential passes and source preservation.
helper_large_pvar <- file.path(tmp, "helper_large.pvar")
helper_large_connection <- file(helper_large_pvar, open = "wt")
writeLines("#CHROM\tPOS\tID\tREF\tALT\tR2", helper_large_connection)
for (first in seq.int(1L, 100000L, by = 10000L)) {
  index <- seq.int(first, min(first + 9999L, 100000L))
  quality <- ifelse(index %% 10L == 0L, "0.5", "0.9")
  writeLines(paste("1", index, paste0("rs", index), "A", "G", quality, sep = "\t"),
    helper_large_connection)
}
close(helper_large_connection)
helper_large_sha <- sha256_file(helper_large_pvar)
helper_large <- run_info_helper(helper_large_pvar, file.path(tmp, "helper_large"))
helper_large_summary <- read_plink_metadata_summary(helper_large$summary)
if (as.integer(helper_large_summary$rows_scanned) != 100000L) {
  stop("large INFO/R2 fixture row count was ", helper_large_summary$rows_scanned)
}
stopifnot(
  identical(sha256_file(helper_large_pvar), helper_large_sha),
  as.integer(helper_large_summary$below_min) == 10000L,
  count_lines(helper_large$pass) == 90000L
)

samples <- file.path(tmp, "samples.tsv")
traits <- file.path(tmp, "traits.tsv")
config <- file.path(tmp, "config.yaml")

write_lines(c(
  "FID\tIID\tage\tage2\tsex\tbatch\tbt1\tbt2\tbt3\tqt1",
  "F1\tI1\t40\t1600\t1\t1\t1\t0\t1\t1.2",
  "F2\tI2\t42\t1764\t2\t1\t0\t1\t0\t2.4",
  "F3\tI3\t50\t2500\t1\t2\t1\t0\t.\tNA",
  "F4\tI4\t52\t2704\t2\t2\t0\t1\t0\t4.8"
), samples)

write_lines(c(
  "trait_id\tphenotype_column\tcase_value\tcontrol_value\tmissing_values\tcovariates",
  "bt1\tbt1\t1\t0\tNA\t",
  "bt2\tbt2\t1\t0\tNA\tbatch",
  "bt3\tbt3\t1\t0\tNA\t",
  "qt1\tqt1\t\t\tNA\t"
), traits)

write_lines(c(
  "project:",
  "  analysis_name: phase2_test",
  "  cohort_data_release: test",
  "inputs:",
  paste0("  sample_manifest: ", samples),
  paste0("  trait_registry: ", traits),
  "tools:",
  "  regenie: regenie",
  "  plink2: plink2",
  "genotypes:",
  "  type: pgen",
  "  prefix: dummy",
  "qc:",
  "  geno_missing_max: 0.05",
  "  sample_missing_max: 0.05",
  "  use_mach_r2_filter: false",
  "warnings:",
  "  min_n: 2",
  "  min_cases: 1",
  "  min_controls: 1",
  "ancestry_reference:",
  "  exclusion_regions: ''",
  "phase2_regenie:",
  "  enabled: true",
  "  global_pcs: 2",
  "  default_covariates: [age, age2, sex, PC1, PC2]",
  "  extra_covariates: []",
  "  apply_rint: false",
  "  htp_cohort_name: ''",
  "  min_mac: 1",
  "  p_thresh: 0.01",
  "  step1_bsize: 100",
  "  step2_bsize: 100",
  "  step1_options: ''",
  "  step2_options: ''",
  "  step2:",
  "    filters: {maf_min: 0.01}",
  "  global_pca:",
  "    filters: {maf_min: 0.01, geno_missing_max: 0.02, snps_only_acgt: true, autosome_only: true, max_alleles: 2, remove_duplicate_ids: true}",
  "    ld_prune: {window: 500kb, step: 1, r2: 0.2}",
  "  step1:",
  "    filters: {maf_min: 0.01, geno_missing_max: 0.02, snps_only_acgt: true, autosome_only: true, max_alleles: 2, remove_duplicate_ids: true}",
  "    ld_prune: {window: 1000kb, step: 1, r2: 0.2}",
  "remeta:",
  "  enabled: true"
), config)

groups <- file.path(tmp, "groups.tsv")
run_phase2(c("write-groups", "--config", config, "--out", groups))
group_rows <- read_tsv(groups)
if (nrow(group_rows) != 4) stop("expected one Phase 2 group per trait")
if (!setequal(group_rows$trait_type, c("bt", "qt"))) stop("expected binary and quantitative groups")
if (any(grepl(",", group_rows$traits, fixed = TRUE))) stop("Phase 2 groups are not singleton")
expected_groups <- c(bt1 = "bt__bt1", bt2 = "bt__bt2", bt3 = "bt__bt3", qt1 = "qt__qt1")
observed_groups <- setNames(group_rows$group, group_rows$traits)
if (!identical(observed_groups[names(expected_groups)], expected_groups)) stop("Phase 2 group IDs are not stable trait-derived IDs")
if (!any(group_rows$traits == "bt2" & grepl("batch", group_rows$covariates))) stop("trait-specific covariate was not retained")

whitespace_traits <- file.path(tmp, "whitespace_traits.tsv")
whitespace_config <- file.path(tmp, "whitespace_config.yaml")
write_lines(c(readLines(traits)[[1]], " bt1\tbt1\t1\t0\tNA\t"), whitespace_traits)
whitespace_config_lines <- readLines(config)
whitespace_config_lines[grepl("^  trait_registry:", whitespace_config_lines)] <-
  paste0("  trait_registry: ", whitespace_traits)
write_lines(whitespace_config_lines, whitespace_config)
invisible(run_phase2(c(
  "write-groups", "--config", whitespace_config, "--out", file.path(tmp, "whitespace_groups.tsv")
), expect_success = FALSE))

stats1 <- file.path(tmp, "stage1a.tsv")
stats2 <- file.path(tmp, "stage1b.tsv")
union <- file.path(tmp, "union.snplist")
union_summary <- file.path(tmp, "union.tsv")
write_lines(c("trait\tancestry\tvariant_id\tp", "bt1\tA\trs2\t0.2", "bt1\tA\trs1\t0.1"), stats1)
write_lines(c("ID\tP", "rs3\t0.3", "rs2\t0.2"), stats2)
run_phase2(c("stage1-pass-union", "--config", config, "--stage1-stats", stats1, stats2, "--out", union, "--summary-out", union_summary))
if (!identical(readLines(union), c("rs1", "rs2", "rs3"))) stop("Stage 1 pass-list union is wrong")

sscore <- file.path(tmp, "global_projected.sscore")
global_pcs <- file.path(tmp, "global_pcs.tsv")
write_lines(c(
  "#IID\tALLELE_CT\tPC1_AVG\tPC2_AVG",
  "I1\t100\t0.11\t0.21",
  "I2\t100\t0.12\t0.22"
), sscore)
run_phase2(c("write-global-pcs", "--config", config, "--sscore", sscore, "--out", global_pcs))
global_pc_rows <- read_tsv(global_pcs)
if (!identical(names(global_pc_rows), c("FID", "IID", "PC1", "PC2"))) stop("global PC parser wrote unexpected columns")
if (!identical(global_pc_rows$FID, global_pc_rows$IID)) stop("IID-only global PC file did not fall back to FID=IID")

default_pc_config <- file.path(tmp, "default_pc_config.yaml")
default_pc_sscore <- file.path(tmp, "default_pc_projected.sscore")
default_pc_out <- file.path(tmp, "default_pc_global_pcs.tsv")
write_lines(c(
  "phase2_regenie:",
  "  enabled: true"
), default_pc_config)
write_lines(c(
  paste(c("#IID", "ALLELE_CT", paste0("PC", 1:10, "_AVG")), collapse = "\t"),
  paste(c("I1", "100", sprintf("%.2f", seq(0.11, 0.20, by = 0.01))), collapse = "\t")
), default_pc_sscore)
run_phase2(c("write-global-pcs", "--config", default_pc_config, "--sscore", default_pc_sscore, "--out", default_pc_out))
default_pc_rows <- read_tsv(default_pc_out)
if (!identical(names(default_pc_rows), c("FID", "IID", paste0("PC", 1:10)))) stop("omitted phase2_regenie.global_pcs did not default to 10 PCs")

alias_keep <- file.path(tmp, "alias.keep.tsv")
alias_pcs <- file.path(tmp, "alias_pcs.tsv")
alias_pheno <- file.path(tmp, "alias.pheno.tsv")
alias_covar <- file.path(tmp, "alias.covar.tsv")
alias_summary <- file.path(tmp, "alias.summary.tsv")
alias_traits <- file.path(tmp, "alias.traits.txt")
alias_covars <- file.path(tmp, "alias.covars.txt")
alias_plink_keep <- file.path(tmp, "alias.plink.keep.txt")
bt1_group <- group_rows$group[group_rows$traits == "bt1"][[1]]
bt3_group <- group_rows$group[group_rows$traits == "bt3"][[1]]
write_lines(c("FID\tIID", "I1\tI1", "I2\tI2"), alias_keep)
write_lines(c("FID\tIID\tPC1\tPC2", "I1\tI1\t0.11\t0.21", "I2\tI2\t0.12\t0.22"), alias_pcs)
run_phase2(c(
  "build-group-inputs", "--config", config, "--group", bt3_group,
  "--keep", alias_keep, "--pcs", alias_pcs,
  "--pheno-out", alias_pheno, "--covar-out", alias_covar,
  "--summary-out", alias_summary, "--trait-list-out", alias_traits,
  "--covar-list-out", alias_covars, "--keep-plink-out", alias_plink_keep
))
alias_covar_rows <- read_tsv(alias_covar)
if (!identical(as.character(alias_covar_rows$age), c("40", "42"))) stop("Phase 2 group input builder did not match IID-only keep IDs to the manifest")
if (!identical(as.character(alias_covar_rows$PC1), c("0.11", "0.12"))) stop("Phase 2 group input builder did not match IID-only PC IDs")
alias_summary_rows <- read_tsv(alias_summary)
if (alias_summary_rows$model_sample_count != 2L || alias_summary_rows$usable_n != 2L) {
  stop("active singleton summary does not report the exact model keep")
}

status_skip_keep <- file.path(tmp, "status_skip.keep.tsv")
status_skip_pcs <- file.path(tmp, "status_skip.pcs.tsv")
status_skip_summary <- file.path(tmp, "status_skip.summary.tsv")
status_skip_traits <- file.path(tmp, "status_skip.traits.txt")
status_skip_model_keep <- file.path(tmp, "status_skip.model.keep.txt")
status_traits_registry <- file.path(tmp, "status_traits.tsv")
status_config <- file.path(tmp, "status_config.yaml")
write_lines(c(
  readLines(traits)[[1]],
  readLines(traits)[[4]],
  readLines(traits)[[2]]
), status_traits_registry)
status_config_lines <- readLines(config)
status_config_lines[grepl("^  trait_registry:", status_config_lines)] <-
  paste0("  trait_registry: ", status_traits_registry)
write_lines(status_config_lines, status_config)
write_lines(c("FID\tIID", "I1\tI1"), status_skip_keep)
write_lines(c("FID\tIID\tPC1\tPC2", "I1\tI1\t0.11\t0.21"), status_skip_pcs)
run_phase2(c(
  "build-group-inputs", "--config", config, "--group", bt1_group,
  "--keep", status_skip_keep, "--pcs", status_skip_pcs,
  "--pheno-out", file.path(tmp, "status_skip.pheno.tsv"),
  "--covar-out", file.path(tmp, "status_skip.covar.tsv"),
  "--summary-out", status_skip_summary, "--trait-list-out", status_skip_traits,
  "--covar-list-out", file.path(tmp, "status_skip.covars.txt"),
  "--keep-plink-out", status_skip_model_keep
))
status_out <- file.path(tmp, "trait_group_status.tsv")
run_phase2(c(
  "select-active-groups", "--config", status_config,
  "--status-summary", alias_summary, status_skip_summary,
  "--status-trait-list", alias_traits, status_skip_traits,
  "--status-keep", alias_plink_keep, status_skip_model_keep,
  "--out", status_out
))
status_rows <- read_tsv(status_out)
if (!identical(status_rows$remeta_eligible, c("True", "False"))) stop("active-group selection did not separate active and skipped models")
if (!identical(status_rows$keep_count, status_rows$model_sample_count)) stop("status keep and model counts differ")
if (any(!nzchar(status_rows$model_keep_sha256))) stop("status table omitted a model keep checksum")
if (file.info(status_skip_model_keep)$size != 0) stop("skipped singleton group wrote a nonempty model keep")

inconsistent_summary <- file.path(tmp, "inconsistent.summary.tsv")
inconsistent_rows <- read_tsv(alias_summary)
inconsistent_rows$usable_n <- inconsistent_rows$usable_n - 1L
write_tsv(inconsistent_rows, inconsistent_summary)
invisible(run_phase2(c(
  "select-active-groups", "--config", status_config,
  "--status-summary", inconsistent_summary, status_skip_summary,
  "--status-trait-list", alias_traits, status_skip_traits,
  "--status-keep", alias_plink_keep, status_skip_model_keep,
  "--out", file.path(tmp, "inconsistent_status.tsv")
), expect_success = FALSE))

partial_keep <- file.path(tmp, "partial.keep.tsv")
partial_pcs <- file.path(tmp, "partial_pcs.tsv")
partial_pheno <- file.path(tmp, "partial.pheno.tsv")
partial_covar <- file.path(tmp, "partial.covar.tsv")
partial_summary <- file.path(tmp, "partial.summary.tsv")
partial_traits <- file.path(tmp, "partial.traits.txt")
partial_covars <- file.path(tmp, "partial.covars.txt")
partial_plink_keep <- file.path(tmp, "partial.plink.keep.txt")
write_lines(c("FID\tIID", "I1\tI1", "I2\tI2", "I3\tI3"), partial_keep)
write_lines(c("FID\tIID\tPC1\tPC2", "I1\tI1\t0.11\t0.21", "I2\tI2\t0.12\t0.22"), partial_pcs)
run_phase2(c(
  "build-group-inputs", "--config", config, "--group", bt1_group,
  "--keep", partial_keep, "--pcs", partial_pcs,
  "--pheno-out", partial_pheno, "--covar-out", partial_covar,
  "--summary-out", partial_summary, "--trait-list-out", partial_traits,
  "--covar-list-out", partial_covars, "--keep-plink-out", partial_plink_keep
))
if (!identical(readLines(partial_plink_keep), c("I1\tI1", "I2\tI2"))) {
  stop("Phase 2 regenie keep file was not limited to covariate-complete samples")
}

qt_keep <- file.path(tmp, "qt.keep.tsv")
qt_pcs <- file.path(tmp, "qt_pcs.tsv")
qt_pheno <- file.path(tmp, "qt.pheno.tsv")
qt_covar <- file.path(tmp, "qt.covar.tsv")
qt_summary <- file.path(tmp, "qt.summary.tsv")
qt_traits <- file.path(tmp, "qt.traits.txt")
qt_covars <- file.path(tmp, "qt.covars.txt")
qt_plink_keep <- file.path(tmp, "qt.plink.keep.txt")
qt_group <- group_rows$group[group_rows$traits == "qt1"][[1]]
write_lines(c("FID\tIID", "I1\tI1", "I2\tI2", "I3\tI3", "I4\tI4"), qt_keep)
write_lines(c(
  "FID\tIID\tPC1\tPC2",
  "I1\tI1\t0.11\t0.21",
  "I2\tI2\t0.12\t0.22",
  "I3\tI3\t0.13\t0.23",
  "I4\tI4\t0.14\t0.24"
), qt_pcs)
run_phase2(c(
  "build-group-inputs", "--config", config, "--group", qt_group,
  "--keep", qt_keep, "--pcs", qt_pcs,
  "--pheno-out", qt_pheno, "--covar-out", qt_covar,
  "--summary-out", qt_summary, "--trait-list-out", qt_traits,
  "--covar-list-out", qt_covars, "--keep-plink-out", qt_plink_keep
))
if (!identical(readLines(qt_plink_keep), c("I1\tI1", "I2\tI2", "I4\tI4"))) {
  stop("Phase 2 regenie Step 1 keep file was not limited to phenotype-complete samples")
}

default_missing_keep <- file.path(tmp, "default_missing.keep.txt")
run_phase2(c(
  "build-group-inputs", "--config", config, "--group", bt3_group,
  "--keep", qt_keep, "--pcs", qt_pcs,
  "--pheno-out", file.path(tmp, "default_missing.pheno.tsv"),
  "--covar-out", file.path(tmp, "default_missing.covar.tsv"),
  "--summary-out", file.path(tmp, "default_missing.summary.tsv"),
  "--trait-list-out", file.path(tmp, "default_missing.traits.txt"),
  "--covar-list-out", file.path(tmp, "default_missing.covars.txt"),
  "--keep-plink-out", default_missing_keep
))
if (!identical(readLines(default_missing_keep), c("I1\tI1", "I2\tI2", "I4\tI4"))) {
  stop("binary default missing code was not excluded consistently")
}

skipped_group_summary <- file.path(tmp, "skipped_group_summary.tsv")
write_lines(c(
  "group\ttrait\ttrait_type\tcovariates\tphase2_pan_samples\tcomplete_covariate_samples\tusable_n\tcases\tcontrols\tmodel_sample_count\tmodel_cases\tmodel_controls\tmodel_keep_sha256\tskipped\tskip_reason",
  "bt__bt1\tbt1\tbt\tage,age2,sex,PC1,PC2\t4\t4\t1\t1\t0\t0\t\t\t\tTrue\tbelow_phase2_thresholds:n=1<min_n=2;controls=0<min_controls=1"
), skipped_group_summary)

skipped_stats <- file.path(tmp, "bt1.regenie")
skipped_summary <- file.path(tmp, "bt1.summary.tsv")
run_phase2(c(
  "stage-trait-output", "--config", config, "--trait", "bt1",
  "--group-summary", skipped_group_summary, "--raw-prefix", file.path(tmp, "raw"),
  "--out-stats", skipped_stats, "--out-summary", skipped_summary
))
if (!any(grepl("^## skipped:", readLines(skipped_stats)))) stop("skipped placeholder missing")

stale_raw_prefix <- file.path(tmp, "stale", "bt__bt1.GRCh38")
write_lines("stale", paste0(stale_raw_prefix, "_bt1.regenie"))
invisible(run_phase2(c(
  "stage-trait-output", "--config", config, "--trait", "bt1",
  "--group-summary", skipped_group_summary, "--raw-prefix", stale_raw_prefix,
  "--out-stats", file.path(tmp, "stale_skipped.regenie"),
  "--out-summary", file.path(tmp, "stale_skipped.summary.tsv")
), expect_success = FALSE))

raw <- paste0(file.path(tmp, "raw"), "_bt2.regenie")
write_lines(c(
  paste(c("Name", "Chr", "Pos", "Ref", "Alt", "Cohort", "Model", "Effect", "LCI_effect", "UCI_effect", "Pval", "AAF", "Num_Cases", "Cases_Ref", "Cases_Het", "Cases_Alt", "Num_Controls", "Controls_Ref", "Controls_Het", "Controls_Alt", "Info"), collapse = " "),
  "rs1 1 100 A G phase2_test ADD 1.1 1.0 1.2 0.01 0.2 2 1 1 0 2 1 1 0 NA"
), raw)
native_stats <- file.path(tmp, "bt2.regenie")
native_summary <- file.path(tmp, "bt2.summary.tsv")
native_keep <- file.path(tmp, "bt2.keep.txt")
native_ids <- file.path(tmp, "raw_bt2.regenie.ids")
write_lines(c("F1\tI1", "F2\tI2", "F3\tI3", "F4\tI4"), native_keep)
write_lines(c("F4\tI4", "F2\tI2", "F1\tI1", "F3\tI3"), native_ids)
group_summary <- file.path(tmp, "group_summary.tsv")
write_lines(c(
  "group\ttrait\ttrait_type\tcovariates\tphase2_pan_samples\tcomplete_covariate_samples\tusable_n\tcases\tcontrols\tmodel_sample_count\tmodel_cases\tmodel_controls\tmodel_keep_sha256\tskipped\tskip_reason",
  paste(c("bt__bt2", "bt2", "bt", "age,age2,sex,PC1,PC2,batch", "4", "4", "4", "2", "2", "4", "2", "2", sha256_file(native_keep), "False", ""), collapse = "\t")
), group_summary)
run_phase2(c(
  "stage-trait-output", "--config", config, "--trait", "bt2",
  "--group-summary", group_summary, "--raw-prefix", file.path(tmp, "raw"),
  "--sample-ids", native_ids, "--model-keep", native_keep,
  "--out-stats", native_stats, "--out-summary", native_summary
))
if (!identical(readLines(raw), readLines(native_stats))) stop("native regenie output was not preserved")

pan_summary <- file.path(tmp, "pan_summary.tsv")
ancestry <- file.path(tmp, "ancestry.tsv")
stage1_summary <- file.path(tmp, "stage1_summary.tsv")
remeta_validation <- file.path(tmp, "remeta_validation.tsv")
report <- file.path(tmp, "report.md")
stats_metrics <- file.path(tmp, "association_metrics.tsv")
top_hits <- file.path(tmp, "top_hits.tsv")
write_lines(c("metric\tvalue", "phase2_pan_samples\t4"), pan_summary)
write_lines(c("FID\tIID\tphase2_ancestry", "F1\tI1\tEUR", "F2\tI2\tUNKNOWN"), ancestry)
write_lines(c("metric\tvalue", "lambda_gc\t1.020000", "valid_p_value_variants\t100"), stage1_summary)
write_lines(c(
  "metric\tvalue",
  "total_variants\t12",
  "valid_p_value_variants\t12",
  "lambda_gc\t1.010000",
  "genomewide_significant_variants\t1",
  "suggestive_variants\t2",
  "qq_eligible_variants\t12",
  "qq_points_plotted\t7",
  "manhattan_eligible_variants\t12",
  "manhattan_points_plotted\t6",
  "large_plot_threshold\t10",
  "plot_thinning_applied\tTrue"
), stats_metrics)
write_lines(c(
  "chrom\tpos\tvariant_id\teffect\tse\tp",
  "1\t100\trs1\t1.1\tNA\t0.01"
), top_hits)
write_lines(c(
  "key\tvalue",
  "status\tvalidated",
  "sample_count\t4",
  "ordinary_regenie_sample_count\t4",
  "rare_regenie_sample_count\t4",
  "target_variant_count\t100",
  "unique_ld_target_variant_count\t98",
  "target_variants_not_indexed\t2",
  "target_variant_ld_coverage_pct\t98.000000",
  "reference_gene_count\t20",
  "indexed_gene_count\t18",
  "genes_without_indexed_variants\t2",
  "indexed_gene_coverage_pct\t90.000000",
  "ld_gene_variant_assignments\t110",
  "ld_assignments_within_gene_bounds\t109",
  "ld_assignments_outside_gene_bounds\t1",
  "ld_assignment_gene_bound_coverage_pct\t99.090909"
), remeta_validation)
run_phase2(c(
  "make-report", "--config", config, "--trait", "bt2", "--build", "GRCh38",
  "--stats", native_stats, "--stats-metrics", stats_metrics, "--top-hits", top_hits,
  "--summary", native_summary, "--group-summary", group_summary,
  "--union-summary", union_summary, "--pan-summary", pan_summary, "--ancestry-summary", ancestry,
  "--qq", "qq.png", "--manhattan", "mh.png", "--manhattan-pdf", "mh.pdf",
  "--stage1-summary", stage1_summary, "--remeta-validation", remeta_validation, "--out", report
))
text <- readLines(report)
if (!any(grepl("Stage 1 Lambda Comparison", text, fixed = TRUE))) stop("Phase 2 report missing Stage 1 comparison")
if (!any(grepl("Step 2 pooled MAF minimum: 0.01", text, fixed = TRUE))) stop("Phase 2 report missing pooled MAF threshold")
if (!any(grepl("Trait-specific group: bt__bt2", text, fixed = TRUE))) stop("Phase 2 report missing singleton group")
if (!any(grepl("Final REGENIE model samples: 4", text, fixed = TRUE))) stop("Phase 2 report missing final model count")
if (!any(grepl("REGENIE sample-ID validation: matched exact model keep", text, fixed = TRUE))) {
  stop("Phase 2 report missing exact sample-ID validation")
}
if (!any(grepl("![QQ plot](qq.png)", text, fixed = TRUE))) stop("Phase 2 report missing embedded QQ plot")
if (!any(grepl("![Manhattan plot](mh.png)", text, fixed = TRUE))) stop("Phase 2 report missing embedded Manhattan plot")
if (!any(grepl("QC-passing target-region variants represented in LD indexes | 98 | 100 | 98.00%", text, fixed = TRUE))) {
  stop("Phase 2 report missing ReMeta target-variant LD coverage")
}
if (!any(grepl("Conditional buffer variants: not included", text, fixed = TRUE))) {
  stop("Phase 2 report missing ReMeta buffer policy")
}
if (!any(grepl("Large PAN plot fallback", text, fixed = TRUE))) {
  stop("Phase 2 report missing large-plot fallback disclosure")
}
if (!any(grepl("Association metrics and top hits use all valid variants", text, fixed = TRUE))) {
  stop("Phase 2 report does not distinguish exact summaries from thinned plots")
}

skipped_report <- file.path(tmp, "skipped_report.md")
skipped_metrics <- file.path(tmp, "skipped_association_metrics.tsv")
skipped_top_hits <- file.path(tmp, "skipped_top_hits.tsv")
write_lines(c(
  "metric\tvalue",
  "total_variants\t0",
  "valid_p_value_variants\t0",
  "lambda_gc\tNA",
  "genomewide_significant_variants\t0",
  "suggestive_variants\t0",
  "qq_eligible_variants\t0",
  "qq_points_plotted\t0",
  "manhattan_eligible_variants\t0",
  "manhattan_points_plotted\t0",
  "large_plot_threshold\t10000000",
  "plot_thinning_applied\tFalse"
), skipped_metrics)
write_lines("chrom\tpos\tvariant_id\teffect\tse\tp", skipped_top_hits)
run_phase2(c(
  "make-report", "--config", config, "--trait", "bt1", "--build", "GRCh38",
  "--stats", skipped_stats, "--stats-metrics", skipped_metrics, "--top-hits", skipped_top_hits,
  "--summary", skipped_summary, "--group-summary", skipped_group_summary,
  "--union-summary", union_summary, "--pan-summary", pan_summary, "--ancestry-summary", ancestry,
  "--qq", "skipped_qq.png", "--manhattan", "skipped_mh.png", "--manhattan-pdf", "skipped_mh.pdf",
  "--stage1-summary", stage1_summary, "--out", skipped_report
))
skipped_text <- readLines(skipped_report)
if (!any(grepl("Skip reason: below_phase2_thresholds:n=1<min_n=2;controls=0<min_controls=1", skipped_text, fixed = TRUE))) {
  stop("skipped Phase 2 report missing its precise threshold reason")
}
if (!any(grepl("ReMeta export skipped", skipped_text, fixed = TRUE))) {
  stop("skipped Phase 2 report missing ReMeta skip disclosure")
}
if (any(grepl("rs1", skipped_text, fixed = TRUE))) stop("skipped report reused active-trait top hits")
pan_sections <- c(
  "## Model Overview",
  "## PAN Sample Set",
  "## Variant Sources and QC",
  "## ReMeta LD Target Coverage",
  "## REGENIE Run Settings",
  "## Association Results",
  "## Top Hits",
  "## Plots",
  "## Stage 1 Lambda Comparison"
)
pan_section_pos <- match(pan_sections, text)
if (any(is.na(pan_section_pos))) {
  stop("Phase 2 report missing section(s): ", paste(pan_sections[is.na(pan_section_pos)], collapse = ", "))
}
if (any(diff(pan_section_pos) <= 0)) stop("Phase 2 report sections are out of order")
if (!any(grepl("rs1", text, fixed = TRUE))) stop("HTP regenie parser did not expose top hit")

fake_regenie <- file.path(tmp, "fake regenie.sh")
write_lines(c(
  "#!/bin/sh",
  "echo fake regenie should not run here >&2",
  "exit 99"
), fake_regenie)
Sys.chmod(fake_regenie, "0755")
cmd_config <- file.path(tmp, "config_cmd.yaml")
cmd_lines <- readLines(config)
cmd_lines <- sub("analysis_name: phase2_test", "analysis_name: Phase 2 Test Cohort", cmd_lines, fixed = TRUE)
cmd_lines <- sub("regenie: regenie", paste0("regenie: '", fake_regenie, "'"), cmd_lines, fixed = TRUE)
write_lines(cmd_lines, cmd_config)

fake_plink2 <- file.path(tmp, "fake_plink2.sh")
write_lines(c(
  "#!/bin/sh",
  "args=\" $* \"",
  "case \"$args\" in *' fill-missing-from-dosage '* ) echo 'Step 2 genotype prep should not fill hardcalls from dosage' >&2; exit 20 ;; esac",
  "case \"$args\" in *' erase-dosage '* ) echo 'Step 2 genotype prep should not erase dosage' >&2; exit 21 ;; esac",
  "case \"$args\" in *' --maf 0.01 '* ) ;; * ) echo 'missing Step 2 --maf 0.01' >&2; exit 22 ;; esac",
  "out=''",
  "while [ \"$#\" -gt 0 ]; do",
  "  if [ \"$1\" = \"--out\" ]; then",
  "    shift",
  "    out=\"$1\"",
  "  fi",
  "  shift",
  "done",
  "if [ -z \"$out\" ]; then",
  "  echo 'missing --out' >&2",
  "  exit 2",
  "fi",
  "case \"$args\" in *' --write-snplist '* ) ;; * ) echo 'missing --write-snplist' >&2; exit 23 ;; esac",
  "printf 'rs1\\n' > \"$out.snplist\""
), fake_plink2)
Sys.chmod(fake_plink2, "0755")
plink_config <- file.path(tmp, "config_fake_plink.yaml")
plink_lines <- readLines(config)
plink_lines <- sub("plink2: plink2", paste0("plink2: '", fake_plink2, "'"), plink_lines, fixed = TRUE)
write_lines(plink_lines, plink_config)
assoc_extract <- file.path(tmp, "assoc_extract.txt")
assoc_prefix <- file.path(tmp, "assoc_norm")
write_lines("rs1", assoc_extract)
run_phase2(c(
  "prepare-assoc-variants", "--config", plink_config, "--pfile-prefix", file.path(tmp, "input"),
  "--extract", assoc_extract, "--out-prefix", assoc_prefix,
  "--summary-out", paste0(assoc_prefix, ".summary.tsv"), "--threads", "1"
))
if (!identical(readLines(paste0(assoc_prefix, ".snplist")), "rs1")) stop("pooled association variant filtering failed")
assoc_summary <- read_tsv(paste0(assoc_prefix, ".summary.tsv"))
if (assoc_summary$value[assoc_summary$metric == "pooled_filter_pass_variants"] != "1") stop("association filter summary has wrong pass count")

fake_pan_plink2 <- file.path(tmp, "fake_pan_plink2.sh")
write_lines(c(
  "#!/bin/sh",
  "out=''",
  "while [ \"$#\" -gt 0 ]; do",
  "  if [ \"$1\" = \"--out\" ]; then shift; out=\"$1\"; fi",
  "  shift",
  "done",
  "printf 'PGEN\\n' > \"$out.pgen\"",
  "printf '#CHROM\\tPOS\\tID\\tREF\\tALT\\n1\\t100\\trs1\\tA\\tG\\n' > \"$out.pvar\"",
  "printf '#IID\\tSEX\\nI1\\t1\\nI2\\t2\\n' > \"$out.psam\""
), fake_pan_plink2)
Sys.chmod(fake_pan_plink2, "0755")
pan_config <- file.path(tmp, "config_fake_pan.yaml")
pan_input <- file.path(tmp, "pan_input")
pan_config_lines <- readLines(config)
pan_config_lines <- sub("plink2: plink2", paste0("plink2: '", fake_pan_plink2, "'"), pan_config_lines, fixed = TRUE)
pan_config_lines <- sub("  prefix: dummy", paste0("  prefix: '", pan_input, "'"), pan_config_lines, fixed = TRUE)
write_lines(pan_config_lines, pan_config)
write_lines("PGEN", paste0(pan_input, ".pgen"))
write_lines(c("#CHROM\tPOS\tID\tREF\tALT", "1\t100\trs1\tA\tG"), paste0(pan_input, ".pvar"))
write_lines(c("#FID\tIID", "F1\tI1", "F2\tI2"), paste0(pan_input, ".psam"))
pan_prefix <- file.path(tmp, "pan")
pan_sex_keep <- file.path(tmp, "pan_sex.keep.txt")
pan_assignments <- file.path(tmp, "pan_assignments.tsv")
pan_excluded <- file.path(tmp, "pan_excluded.tsv")
write_lines(c("F1\tI1", "F2\tI2"), pan_sex_keep)
write_lines("FID\tIID", pan_assignments)
write_lines("FID\tIID", pan_excluded)
run_phase2(c(
  "prepare-pan-genotypes", "--config", pan_config,
  "--sex-keep", pan_sex_keep, "--assignments", pan_assignments, "--excluded", pan_excluded,
  "--out-prefix", pan_prefix, "--keep-out", file.path(tmp, "pan.keep.tsv"),
  "--ancestry-out", file.path(tmp, "pan.ancestry.tsv"),
  "--summary-out", file.path(tmp, "pan.summary.tsv"), "--threads", "1"
))
pan_psam <- read_tsv(paste0(pan_prefix, ".psam"))
if (!identical(names(pan_psam)[1:2], c("#FID", "IID")) || !identical(pan_psam[[1]], pan_psam$IID)) {
  stop("shared PAN PSAM was not normalized to two REGENIE ID columns")
}

fake_marker_plink2 <- file.path(tmp, "fake_marker_plink2.sh")
write_lines(c(
  "#!/bin/sh",
  "args=\" $* \"",
  "out=''",
  "while [ \"$#\" -gt 0 ]; do",
  "  if [ \"$1\" = \"--out\" ]; then",
  "    shift",
  "    out=\"$1\"",
  "  fi",
  "  shift",
  "done",
  "if [ -z \"$out\" ]; then echo 'missing --out' >&2; exit 2; fi",
  "case \"$args\" in",
  "  *' --make-pgen '* )",
  "    case \"$args\" in *' --extract '* ) echo 'zero-removal INFO/R2 filter unexpectedly passed --extract' >&2; exit 30 ;; esac",
  "    case \"$args\" in *' fill-missing-from-dosage '* ) ;; * ) echo 'Step 1 marker prep missing fill-missing-from-dosage' >&2; exit 22 ;; esac",
  "    case \"$args\" in *' erase-dosage '* ) ;; * ) echo 'Step 1 marker prep missing erase-dosage' >&2; exit 23 ;; esac",
  "    printf 'PGEN\\n' > \"$out.pgen\"",
  "    printf '#CHROM\\tPOS\\tID\\tREF\\tALT\\n1\\t100\\trs1\\tA\\tG\\n' > \"$out.pvar\"",
  "    printf '#IID\\nI1\\nI2\\n' > \"$out.psam\"",
  "    exit 0",
  "    ;;",
  "  *' --indep-pairwise '* )",
  "    case \"$args\" in *' --indep-pairwise 1000kb 1 0.2 '* ) ;; * ) echo 'Step 1 marker prep used wrong LD pruning settings' >&2; exit 28 ;; esac",
  "    printf 'rs1\\n' > \"$out.prune.in\"",
  "    exit 0",
  "    ;;",
  "esac",
  "echo 'unexpected fake marker PLINK2 command' >&2",
  "exit 24"
), fake_marker_plink2)
Sys.chmod(fake_marker_plink2, "0755")
marker_config <- file.path(tmp, "config_marker_plink.yaml")
marker_lines <- readLines(config)
marker_lines <- sub("plink2: plink2", paste0("plink2: '", fake_marker_plink2, "'"), marker_lines, fixed = TRUE)
write_lines(marker_lines, marker_config)
marker_input <- file.path(tmp, "marker_input")
write_lines(c(
  "#CHROM\tPOS\tID\tREF\tALT\tR2",
  "1\t100\trs1\tA\tG\t0.95"
), paste0(marker_input, ".pvar"))
marker_prefix <- file.path(tmp, "step1_marker_qc")
marker_prune_prefix <- file.path(tmp, "step1_marker_prune")
marker_output <- run_phase2(c(
  "prepare-marker-set", "--config", marker_config, "--branch", "step1",
  "--pfile-prefix", marker_input, "--out-prefix", marker_prefix,
  "--prune-prefix", marker_prune_prefix, "--prune-in", paste0(marker_prune_prefix, ".prune.in"),
  "--excluded-regions", file.path(tmp, "step1_marker.excluded_regions.txt"),
  "--threads", "1"
))
if (!file.exists(paste0(marker_prune_prefix, ".prune.in"))) stop("Step 1 marker prep did not produce prune.in")
if (!any(grepl("threshold=0.9; removed=0", marker_output, fixed = TRUE))) {
  stop("Step 1 INFO/R2 fallback was not 0.9 or zero-removal handling failed")
}

fake_marker_info_plink2 <- file.path(tmp, "fake_marker_info_plink2.sh")
write_lines(c(
  "#!/bin/sh",
  "args=\" $* \"",
  "out=''",
  "extract=''",
  "while [ \"$#\" -gt 0 ]; do",
  "  if [ \"$1\" = \"--out\" ]; then",
  "    shift",
  "    out=\"$1\"",
  "  elif [ \"$1\" = \"--extract\" ]; then",
  "    shift",
  "    extract=\"$1\"",
  "  fi",
  "  shift",
  "done",
  "if [ -z \"$out\" ]; then echo 'missing --out' >&2; exit 2; fi",
  "case \"$args\" in",
  "  *' --make-pgen '* )",
  "    if [ -z \"$extract\" ]; then echo 'Step 1 INFO/R2 filter did not pass --extract' >&2; exit 25; fi",
  "    expected=$(printf 'rs_high\\nrs_unimputed')",
  "    observed=$(cat \"$extract\")",
  "    if [ \"$observed\" != \"$expected\" ]; then echo 'Step 1 INFO/R2 pass list is wrong' >&2; cat \"$extract\" >&2; exit 26; fi",
  "    printf 'PGEN\\n' > \"$out.pgen\"",
  "    printf '##fileformat=VCFv4.3\\nID\\tALT\\tPOS\\tREF\\t#CHROM\\nrs_high\\tG\\t100\\tA\\tchr1\\nrs_unimputed\\tG\\t103\\tA\\tchr1\\n' > \"$out.pvar\"",
  "    printf '#IID\\nI1\\nI2\\n' > \"$out.psam\"",
  "    exit 0",
  "    ;;",
  "  *' --indep-pairwise '* )",
  "    case \"$args\" in *' --indep-pairwise 1000kb 1 0.2 '* ) ;; * ) echo 'Step 1 INFO marker prep used wrong LD pruning settings' >&2; exit 29 ;; esac",
  "    printf 'rs_high\\nrs_unimputed\\n' > \"$out.prune.in\"",
  "    exit 0",
  "    ;;",
  "esac",
  "echo 'unexpected fake marker INFO PLINK2 command' >&2",
  "exit 27"
), fake_marker_info_plink2)
Sys.chmod(fake_marker_info_plink2, "0755")
marker_info_config <- file.path(tmp, "config_marker_info_plink.yaml")
marker_info_regions <- file.path(tmp, "marker_info_regions.tsv")
write_lines(c(
  "chrom\tstart\tend\tlabel",
  "1\t103\t103\ttest_region"
), marker_info_regions)
marker_info_lines <- readLines(config)
marker_info_lines <- sub("plink2: plink2", paste0("plink2: '", fake_marker_info_plink2, "'"), marker_info_lines, fixed = TRUE)
marker_info_lines <- sub(
  "  exclusion_regions: ''",
  paste0("  exclusion_regions: '", marker_info_regions, "'"),
  marker_info_lines,
  fixed = TRUE
)
write_lines(marker_info_lines, marker_info_config)
marker_info_input <- file.path(tmp, "marker_info_input")
write_lines(c(
  "#CHROM\tPOS\tID\tREF\tALT\tR2",
  "1\t100\trs_high\tA\tG\t0.95",
  "1\t101\trs_between\tA\tG\t0.85",
  "1\t102\trs_low\tA\tG\t0.50",
  "1\t103\trs_unimputed\tA\tG\t."
), paste0(marker_info_input, ".pvar"))
marker_info_prefix <- file.path(tmp, "step1_marker_info_qc")
marker_info_prune_prefix <- file.path(tmp, "step1_marker_info_prune")
marker_info_region_excluded <- file.path(tmp, "step1_marker_info.excluded_regions.txt")
run_phase2(c(
  "prepare-marker-set", "--config", marker_info_config, "--branch", "step1",
  "--pfile-prefix", marker_info_input, "--out-prefix", marker_info_prefix,
  "--prune-prefix", marker_info_prune_prefix, "--prune-in", paste0(marker_info_prune_prefix, ".prune.in"),
  "--excluded-regions", marker_info_region_excluded,
  "--threads", "1"
))
info_excluded <- read_tsv(paste0(marker_info_prefix, ".info_r2.excluded.tsv"))
if (!identical(info_excluded$variant_id, c("rs_between", "rs_low"))) {
  stop("Step 1 INFO/R2 filter did not apply the default 0.9 threshold")
}
if (!identical(readLines(marker_info_region_excluded), "rs_unimputed")) {
  stop("Step 1 long-range-region filter excluded the wrong marker")
}

fake_filter_plink2 <- file.path(tmp, "fake_filter_plink2.sh")
write_lines(c(
  "#!/bin/sh",
  "args=\" $* \"",
  "mode=''",
  "case \"$args\" in *' --write-snplist '* ) mode='qc' ;; esac",
  "case \"$args\" in *' --glm '* ) echo 'old model-check --glm path should not run' >&2; exit 14 ;; esac",
  "case \"$args\" in *' --export Av '* ) echo 'old residual-variance export path should not run' >&2; exit 13 ;; esac",
  "if [ \"$mode\" = 'qc' ]; then",
  "  case \"$args\" in *' --mac 1 '* ) echo 'unexpected --mac 1 fallback' >&2; exit 3 ;; esac",
  "  case \"$args\" in *' --maf 0.01 '* ) ;; * ) echo 'missing configured --maf 0.01' >&2; exit 4 ;; esac",
  "  case \"$args\" in *' --mac 100 '* ) ;; * ) echo 'missing default Step 1 --mac 100' >&2; exit 9 ;; esac",
  "  case \"$args\" in *' --geno 0.02 '* ) ;; * ) echo 'missing configured --geno 0.02' >&2; exit 5 ;; esac",
  "  case \"$args\" in *' --snps-only just-acgt '* ) ;; * ) echo 'missing configured SNP allele filter' >&2; exit 6 ;; esac",
  "  case \"$args\" in *' --keep '* ) ;; * ) echo 'missing --keep' >&2; exit 7 ;; esac",
  "  case \"$args\" in *' --extract '* ) ;; * ) echo 'missing --extract' >&2; exit 8 ;; esac",
  "  case \"$args\" in *' --nonfounders '* ) ;; * ) echo 'missing --nonfounders' >&2; exit 15 ;; esac",
  "  case \"$args\" in *' --geno-counts '* ) ;; * ) echo 'missing --geno-counts' >&2; exit 16 ;; esac",
  "  case \"$args\" in *' cols=chrom,pos,ref,alt1,homref,refalt1,homalt1,missing,nobs '* ) ;; * ) echo 'missing expected --geno-counts columns' >&2; exit 17 ;; esac",
  "fi",
  "out=''",
  "while [ \"$#\" -gt 0 ]; do",
  "  if [ \"$1\" = \"--out\" ]; then",
  "    shift",
  "    out=\"$1\"",
  "  fi",
  "  shift",
  "done",
  "if [ -z \"$out\" ]; then",
  "  echo 'missing --out' >&2",
  "  exit 2",
  "fi",
  "if [ \"$mode\" = 'qc' ]; then",
  "  printf 'rs_pass\\nrs_low_mac\\nrs_all_het\\nrs_missing_count\\n' > \"$out.snplist\"",
  "  printf '#CHROM\\tPOS\\tID\\tREF\\tALT1\\tHOM_REF_CT\\tHET_REF_ALT1_CT\\tHOM_ALT1_CT\\tMISSING_CT\\tOBS_CT\\n' > \"$out.gcount\"",
  "  printf '1\\t100\\trs_pass\\tA\\tG\\t100\\t100\\t0\\t0\\t200\\n' >> \"$out.gcount\"",
  "  printf '1\\t101\\trs_low_mac\\tA\\tG\\t150\\t99\\t0\\t0\\t249\\n' >> \"$out.gcount\"",
  "  printf '1\\t102\\trs_all_het\\tA\\tG\\t0\\t200\\t0\\t0\\t200\\n' >> \"$out.gcount\"",
  "  exit 0",
  "fi",
  "echo 'unexpected fake PLINK2 mode' >&2",
  "exit 12"
), fake_filter_plink2)
Sys.chmod(fake_filter_plink2, "0755")
filter_config <- file.path(tmp, "config_filter_plink.yaml")
filter_lines <- readLines(config)
filter_lines <- sub("plink2: plink2", paste0("plink2: '", fake_filter_plink2, "'"), filter_lines, fixed = TRUE)
write_lines(filter_lines, filter_config)
filter_extract <- file.path(tmp, "step1_prune.in")
filter_keep <- file.path(tmp, "step1.keep.txt")
filter_traits <- file.path(tmp, "step1.traits.txt")
filter_out <- file.path(tmp, "step1.filtered.snplist")
filter_summary <- file.path(tmp, "step1.variant_qc.summary.tsv")
filter_excluded <- file.path(tmp, "step1.variant_qc.excluded.tsv")
write_lines(c("rs_pass", "rs_low_mac", "rs_all_het", "rs_missing_count"), filter_extract)
write_lines(c("I1\tI1", "I2\tI2", "I3\tI3", "I4\tI4"), filter_keep)
write_lines("bt1", filter_traits)
run_phase2(c(
  "filter-step1-variants", "--config", filter_config, "--pfile-prefix", file.path(tmp, "step1_qc"),
  "--extract", filter_extract, "--keep", filter_keep, "--trait-list", filter_traits,
  "--out", filter_out, "--summary-out", filter_summary, "--excluded-out", filter_excluded,
  "--threads", "1"
))
if (!identical(readLines(filter_out), "rs_pass")) {
  stop("Step 1 variant filter did not stage the hardcall-count QC snplist")
}
filter_summary_rows <- read_tsv(filter_summary)
expected_summary_cols <- c(
  "filter_method", "plink_nonfounders", "raw_plink_pass_snp_count",
  "hardcall_filter_pass_snp_count", "excluded_snp_count", "model_sample_count",
  "hardcall_mac_min", "hardcall_variance_min"
)
if (!identical(names(filter_summary_rows), expected_summary_cols)) stop("Step 1 hardcall-count summary has wrong columns")
if (!identical(filter_summary_rows$filter_method, "plink2_hardcall_count_qc")) stop("Step 1 hardcall-count summary has wrong method")
if (!identical(filter_summary_rows$plink_nonfounders, "True")) stop("Step 1 hardcall-count summary has wrong nonfounder flag")
if (!identical(as.integer(filter_summary_rows$raw_plink_pass_snp_count), 4L)) stop("Step 1 hardcall-count summary has wrong raw count")
if (!identical(as.integer(filter_summary_rows$hardcall_filter_pass_snp_count), 1L)) stop("Step 1 hardcall-count summary has wrong pass count")
if (!identical(as.integer(filter_summary_rows$excluded_snp_count), 3L)) stop("Step 1 hardcall-count summary has wrong excluded count")
if (!identical(as.integer(filter_summary_rows$model_sample_count), 4L)) stop("Step 1 hardcall-count summary has wrong sample count")
if (!identical(as.integer(filter_summary_rows$hardcall_mac_min), 100L)) stop("Step 1 hardcall-count summary has wrong MAC threshold")
if (!identical(as.numeric(filter_summary_rows$hardcall_variance_min), 0)) stop("Step 1 hardcall-count summary has wrong variance threshold")
filter_excluded_rows <- read_tsv(filter_excluded)
expected_excluded_cols <- c(
  "variant_id", "hardcall_ref_ct", "hardcall_alt_ct", "hardcall_mac",
  "hardcall_n", "hardcall_variance", "exclusion_reason"
)
if (!identical(names(filter_excluded_rows), expected_excluded_cols)) stop("Step 1 hardcall-count excluded table has wrong columns")
if (!identical(filter_excluded_rows$variant_id, c("rs_low_mac", "rs_all_het", "rs_missing_count"))) stop("Step 1 hardcall-count excluded table has wrong variants")
if (!identical(filter_excluded_rows$exclusion_reason, c("hardcall_mac_below_min", "zero_hardcall_variance", "missing_from_plink2_gcount"))) {
  stop("Step 1 hardcall-count excluded table has wrong exclusion reasons")
}

fake_dup_plink2 <- file.path(tmp, "fake_dup_snplist_plink2.sh")
write_lines(c(
  "#!/bin/sh",
  "out=''",
  "while [ \"$#\" -gt 0 ]; do",
  "  if [ \"$1\" = \"--out\" ]; then",
  "    shift",
  "    out=\"$1\"",
  "  fi",
  "  shift",
  "done",
  "if [ -z \"$out\" ]; then echo 'missing --out' >&2; exit 2; fi",
  "printf 'rs_dup\\nrs_dup\\n' > \"$out.snplist\"",
  "printf '#CHROM\\tPOS\\tID\\tREF\\tALT1\\tHOM_REF_CT\\tHET_REF_ALT1_CT\\tHOM_ALT1_CT\\tMISSING_CT\\tOBS_CT\\n' > \"$out.gcount\"",
  "printf '1\\t100\\trs_dup\\tA\\tG\\t100\\t100\\t0\\t0\\t200\\n' >> \"$out.gcount\""
), fake_dup_plink2)
Sys.chmod(fake_dup_plink2, "0755")
dup_filter_config <- file.path(tmp, "config_dup_filter_plink.yaml")
dup_filter_lines <- readLines(config)
dup_filter_lines <- sub("plink2: plink2", paste0("plink2: '", fake_dup_plink2, "'"), dup_filter_lines, fixed = TRUE)
write_lines(dup_filter_lines, dup_filter_config)
dup_failure <- run_phase2(c(
  "filter-step1-variants", "--config", dup_filter_config, "--pfile-prefix", file.path(tmp, "step1_qc"),
  "--extract", filter_extract, "--keep", filter_keep, "--trait-list", filter_traits,
  "--out", file.path(tmp, "step1.dup.filtered.snplist"), "--summary-out", file.path(tmp, "step1.dup.variant_qc.summary.tsv"),
  "--excluded-out", file.path(tmp, "step1.dup.variant_qc.excluded.tsv"),
  "--threads", "1"
), expect_success = FALSE)
if (!any(grepl("duplicate IDs", dup_failure))) stop("duplicate Step 1 snplist IDs did not fail as expected")

filter_empty_traits <- file.path(tmp, "step1.empty.traits.txt")
write_lines(character(), filter_empty_traits)
empty_filter_out <- file.path(tmp, "step1.empty.filtered.snplist")
empty_filter_summary <- file.path(tmp, "step1.empty.variant_qc.summary.tsv")
empty_filter_excluded <- file.path(tmp, "step1.empty.variant_qc.excluded.tsv")
fake_fail_plink2 <- file.path(tmp, "fake_fail_if_called_plink2.sh")
write_lines(c(
  "#!/bin/sh",
  "echo 'empty trait PLINK2 stub should not be called' >&2",
  "exit 99"
), fake_fail_plink2)
Sys.chmod(fake_fail_plink2, "0755")
empty_filter_config <- file.path(tmp, "config_empty_filter_plink.yaml")
empty_filter_lines <- readLines(config)
empty_filter_lines <- sub("plink2: plink2", paste0("plink2: '", fake_fail_plink2, "'"), empty_filter_lines, fixed = TRUE)
write_lines(empty_filter_lines, empty_filter_config)
invisible(run_phase2(c(
  "filter-step1-variants", "--config", empty_filter_config, "--pfile-prefix", file.path(tmp, "step1_qc"),
  "--extract", filter_extract, "--keep", filter_keep, "--trait-list", filter_empty_traits,
  "--out", empty_filter_out, "--summary-out", empty_filter_summary, "--excluded-out", empty_filter_excluded,
  "--threads", "1"
), expect_success = FALSE))
if (any(file.exists(c(empty_filter_out, empty_filter_summary, empty_filter_excluded)))) {
  stop("empty Step 1 trait list created fake variant-QC artifacts")
}

trait_list <- file.path(tmp, "bt1.traits.txt")
write_lines("bt1", trait_list)
bt_group <- bt1_group
step1_script <- file.path(tmp, "step1_command.sh")
run_phase2(c(
  "write-step1-command", "--config", cmd_config, "--group", bt_group, "--pfile-prefix", file.path(tmp, "step1_data"),
  "--extract", file.path(tmp, "extract_variants.txt"), "--pheno", file.path(tmp, "pheno.tsv"),
  "--covar", file.path(tmp, "covar.tsv"), "--keep", file.path(tmp, "keep_ids.txt"),
  "--trait-list", trait_list, "--pred-list", file.path(tmp, "step1_pred.list"),
  "--out-prefix", file.path(tmp, "step1"), "--script-out", step1_script, "--threads", "2"
))
step1_text <- paste(readLines(step1_script), collapse = "\n")
if (!grepl("'--step' '1'", step1_text, fixed = TRUE)) stop("Step 1 command script missing --step 1")
if (!grepl("'--lowmem'", step1_text, fixed = TRUE)) stop("Step 1 command script missing --lowmem")
if (!grepl("'--extract' '", step1_text, fixed = TRUE)) stop("Step 1 command script missing --extract")
if (!grepl("'--minCaseCount' '1'", step1_text, fixed = TRUE)) stop("Step 1 command script did not match the case threshold")
if (!grepl("step1_1.loco", step1_text, fixed = TRUE)) stop("Step 1 command script does not verify its prediction payload")
if (!grepl("'[^']*fake regenie.sh'", step1_text)) stop("Step 1 command script did not shell-quote the regenie tool path")

step2_script <- file.path(tmp, "step2_command.sh")
run_phase2(c(
  "write-step2-command", "--config", cmd_config, "--group", bt_group, "--trait", "bt1",
  "--pfile-prefix", file.path(tmp, "assoc_data"), "--extract", file.path(tmp, "assoc.snplist"),
  "--keep", file.path(tmp, "keep_ids.txt"),
  "--pheno", file.path(tmp, "pheno.tsv"), "--covar", file.path(tmp, "covar.tsv"),
  "--pred-list", file.path(tmp, "pred.list"), "--trait-list", trait_list,
  "--out-prefix", file.path(tmp, "step2"), "--done", file.path(tmp, "step2.done"),
  "--script-out", step2_script, "--threads", "2"
))
step2_text <- paste(readLines(step2_script), collapse = "\n")
if (!grepl("'--htp' 'Phase_2_Test_Cohort'", step2_text, fixed = TRUE)) stop("Step 2 command script missing the expected --htp cohort label")
if (!grepl("'--minMAC' '1'", step2_text, fixed = TRUE)) stop("Step 2 command script missing --minMAC 1")
if (!grepl("'--write-samples'", step2_text, fixed = TRUE)) stop("Step 2 command script missing --write-samples")
if (!grepl("'--minCaseCount' '1'", step2_text, fixed = TRUE)) stop("Step 2 command script did not match the case threshold")
if (!grepl("'--keep' '", step2_text, fixed = TRUE)) stop("Step 2 command script missing the exact model keep")
if (!grepl("'--extract' '", step2_text, fixed = TRUE)) stop("Step 2 command script missing the filtered variant list")
if (!grepl("'[^']*fake regenie.sh'", step2_text)) stop("Step 2 command script did not shell-quote the regenie tool path")

empty_traits <- file.path(tmp, "empty.traits.txt")
write_lines(character(), empty_traits)
run_phase2(c(
  "write-step1-command", "--config", cmd_config, "--group", bt_group, "--pfile-prefix", file.path(tmp, "step1"),
  "--extract", file.path(tmp, "extract.txt"), "--pheno", file.path(tmp, "pheno.tsv"),
  "--covar", file.path(tmp, "covar.tsv"), "--keep", file.path(tmp, "keep.txt"),
  "--trait-list", empty_traits, "--pred-list", file.path(tmp, "noop_pred.list"),
  "--out-prefix", file.path(tmp, "noop_step1"), "--script-out", file.path(tmp, "noop_step1.sh"), "--threads", "2"
), expect_success = FALSE)

run_phase2(c(
  "write-step2-command", "--config", cmd_config, "--group", bt_group, "--trait", "bt1",
  "--pfile-prefix", file.path(tmp, "assoc"), "--extract", file.path(tmp, "assoc.snplist"),
  "--keep", file.path(tmp, "keep.txt"),
  "--pheno", file.path(tmp, "pheno.tsv"), "--covar", file.path(tmp, "covar.tsv"),
  "--pred-list", file.path(tmp, "pred.list"), "--trait-list", empty_traits,
  "--out-prefix", file.path(tmp, "noop_step2"), "--done", file.path(tmp, "noop_step2.done"),
  "--script-out", file.path(tmp, "noop_step2.sh"), "--threads", "2"
), expect_success = FALSE)

python <- Sys.which("python3")
if (!nzchar(python)) python <- Sys.which("python")
if (!nzchar(python)) stop("python is required for record_regenie_tool.py test")
fake_versioned_regenie <- file.path(tmp, "fake_versioned_regenie.sh")
write_lines(c(
  "#!/bin/sh",
  "if [ \"$1\" = \"--version\" ]; then",
  "  echo 'regenie fake 4.1.2'",
  "  exit 0",
  "fi",
  "exit 0"
), fake_versioned_regenie)
Sys.chmod(fake_versioned_regenie, "0755")
tool_manifest <- file.path(tmp, "regenie_tool.tsv")
status <- system2(python, c(file.path(repo, "scripts", "record_regenie_tool.py"), "--tool", fake_versioned_regenie, "--out", tool_manifest), stdout = TRUE, stderr = TRUE)
if (!identical(as.integer(attr(status, "status") %||% 0L), 0L)) stop("record_regenie_tool.py failed:\n", paste(status, collapse = "\n"))
tool_rows <- read_tsv(tool_manifest)
if (!all(c("tool_regenie_path", "tool_regenie_sha256", "tool_regenie_version") %in% tool_rows$key)) stop("regenie tool provenance manifest is missing expected rows")
if (!any(tool_rows$key == "tool_regenie_version" & grepl("fake 4.1.2", tool_rows$value))) stop("regenie tool provenance did not record version output")

bad_options <- file.path(tmp, "bad_options.txt")
write_lines("--pgen x", bad_options)
blocked <- run_phase2(c("check-options", "--config", config, "--options-file", bad_options), expect_success = FALSE)
if (!any(grepl("cannot override", blocked))) stop("blocked pass-through option did not fail as expected")

cat("Phase 2 regenie helper tests passed\n")
