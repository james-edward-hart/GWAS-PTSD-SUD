#!/usr/bin/env Rscript

# Focused tests for Phase 2 regenie helper behavior that does not require
# external PLINK2 or regenie execution.

cmd <- commandArgs(FALSE)
script_path <- sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])
repo <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(repo, "scripts", "lib", "stage1.R"))

tmp <- tempfile("phase2-reg-")
dir.create(tmp, recursive = TRUE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

write_lines <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path)
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

samples <- file.path(tmp, "samples.tsv")
traits <- file.path(tmp, "traits.tsv")
config <- file.path(tmp, "config.yaml")

write_lines(c(
  "FID\tIID\tage\tage2\tsex\tbatch\tbt1\tbt2\tqt1",
  "F1\tI1\t40\t1600\t1\t1\t1\t0\t1.2",
  "F2\tI2\t42\t1764\t2\t1\t0\t1\t2.4",
  "F3\tI3\t50\t2500\t1\t2\t1\t0\tNA",
  "F4\tI4\t52\t2704\t2\t2\t0\t1\t4.8"
), samples)

write_lines(c(
  "trait_id\tphenotype_column\tcase_value\tcontrol_value\tmissing_values\tcovariates",
  "bt1\tbt1\t1\t0\tNA\t",
  "bt2\tbt2\t1\t0\tNA\tbatch",
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
  "  global_pca:",
  "    filters: {maf_min: 0.01, geno_missing_max: 0.02, snps_only_acgt: true, autosome_only: true, max_alleles: 2, remove_duplicate_ids: true}",
  "    ld_prune: {window: 500kb, step: 1, r2: 0.2}",
  "  step1:",
  "    filters: {maf_min: 0.01, geno_missing_max: 0.02, snps_only_acgt: true, autosome_only: true, max_alleles: 2, remove_duplicate_ids: true}",
  "    ld_prune: {window: 500kb, step: 1, r2: 0.5}"
), config)

groups <- file.path(tmp, "groups.tsv")
run_phase2(c("write-groups", "--config", config, "--out", groups))
group_rows <- read_tsv(groups)
if (nrow(group_rows) != 3) stop("expected three Phase 2 groups")
if (!setequal(group_rows$trait_type, c("bt", "qt"))) stop("expected binary and quantitative groups")
if (!any(group_rows$traits == "bt2" & grepl("batch", group_rows$covariates))) stop("trait-specific covariate did not split group")

stats1 <- file.path(tmp, "stage1a.tsv")
stats2 <- file.path(tmp, "stage1b.tsv")
union <- file.path(tmp, "union.snplist")
union_summary <- file.path(tmp, "union.tsv")
write_lines(c("trait\tancestry\tvariant_id\tp", "bt1\tA\trs2\t0.2", "bt1\tA\trs1\t0.1"), stats1)
write_lines(c("ID\tP", "rs3\t0.3", "rs2\t0.2"), stats2)
run_phase2(c("stage1-pass-union", "--config", config, "--stage1-stats", stats1, stats2, "--out", union, "--summary-out", union_summary))
if (!identical(readLines(union), c("rs1", "rs2", "rs3"))) stop("Stage 1 pass-list union is wrong")

group_summary <- file.path(tmp, "group_summary.tsv")
write_lines(c(
  "group\ttrait\ttrait_type\tcovariates\tphase2_pan_samples\tcomplete_covariate_samples\tusable_n\tcases\tcontrols\tskipped\tskip_reason",
  "bt_g1\tbt1\tbt\tage,age2,sex,PC1,PC2\t4\t4\t1\t1\t0\tTrue\tbelow_phase2_thresholds:n=1;cases=1;controls=0",
  "bt_g1\tbt2\tbt\tage,age2,sex,PC1,PC2,batch\t4\t4\t4\t2\t2\tFalse\t"
), group_summary)

skipped_stats <- file.path(tmp, "bt1.regenie")
skipped_summary <- file.path(tmp, "bt1.summary.tsv")
run_phase2(c(
  "stage-trait-output", "--config", config, "--trait", "bt1",
  "--group-summary", group_summary, "--raw-prefix", file.path(tmp, "raw"),
  "--out-stats", skipped_stats, "--out-summary", skipped_summary
))
if (!any(grepl("^## skipped:", readLines(skipped_stats)))) stop("skipped placeholder missing")

raw <- paste0(file.path(tmp, "raw"), "_bt2.regenie")
write_lines(c(
  paste(c("Name", "Chr", "Pos", "Ref", "Alt", "Cohort", "Model", "Effect", "LCI_effect", "UCI_effect", "Pval", "AAF", "Num_Cases", "Cases_Ref", "Cases_Het", "Cases_Alt", "Num_Controls", "Controls_Ref", "Controls_Het", "Controls_Alt", "Info"), collapse = " "),
  "rs1 1 100 A G phase2_test ADD 1.1 1.0 1.2 0.01 0.2 2 1 1 0 2 1 1 0 NA"
), raw)
native_stats <- file.path(tmp, "bt2.regenie")
native_summary <- file.path(tmp, "bt2.summary.tsv")
run_phase2(c(
  "stage-trait-output", "--config", config, "--trait", "bt2",
  "--group-summary", group_summary, "--raw-prefix", file.path(tmp, "raw"),
  "--out-stats", native_stats, "--out-summary", native_summary
))
if (!identical(readLines(raw), readLines(native_stats))) stop("native regenie output was not preserved")

pan_summary <- file.path(tmp, "pan_summary.tsv")
ancestry <- file.path(tmp, "ancestry.tsv")
stage1_summary <- file.path(tmp, "stage1_summary.tsv")
report <- file.path(tmp, "report.md")
write_lines(c("metric\tvalue", "phase2_pan_samples\t4"), pan_summary)
write_lines(c("FID\tIID\tphase2_ancestry", "F1\tI1\tEUR", "F2\tI2\tUNKNOWN"), ancestry)
write_lines(c("metric\tvalue", "lambda_gc\t1.020000", "valid_p_value_variants\t100"), stage1_summary)
run_phase2(c(
  "make-report", "--config", config, "--trait", "bt2", "--build", "GRCh38",
  "--stats", native_stats, "--summary", native_summary, "--group-summary", group_summary,
  "--union-summary", union_summary, "--pan-summary", pan_summary, "--ancestry-summary", ancestry,
  "--qq", "qq.png", "--manhattan", "mh.png", "--manhattan-pdf", "mh.pdf",
  "--stage1-summary", stage1_summary, "--out", report
))
text <- readLines(report)
if (!any(grepl("Stage 1 Lambda Comparison", text, fixed = TRUE))) stop("Phase 2 report missing Stage 1 comparison")
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
trait_list <- file.path(tmp, "bt1.traits.txt")
write_lines("bt1", trait_list)
bt_group <- group_rows$group[group_rows$traits == "bt1"][[1]]
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
if (!grepl("'[^']*fake regenie.sh'", step1_text)) stop("Step 1 command script did not shell-quote the regenie tool path")

step2_script <- file.path(tmp, "step2_command.sh")
run_phase2(c(
  "write-step2-command", "--config", cmd_config, "--group", bt_group, "--pfile-prefix", file.path(tmp, "assoc_data"),
  "--pheno", file.path(tmp, "pheno.tsv"), "--covar", file.path(tmp, "covar.tsv"),
  "--pred-list", file.path(tmp, "pred.list"), "--trait-list", trait_list,
  "--out-prefix", file.path(tmp, "step2"), "--done", file.path(tmp, "step2.done"),
  "--script-out", step2_script, "--threads", "2"
))
step2_text <- paste(readLines(step2_script), collapse = "\n")
if (!grepl("'--htp' 'Phase_2_Test_Cohort'", step2_text, fixed = TRUE)) stop("Step 2 command script missing the expected --htp cohort label")
if (!grepl("'--minMAC' '1'", step2_text, fixed = TRUE)) stop("Step 2 command script missing --minMAC 1")
if (!grepl("'[^']*fake regenie.sh'", step2_text)) stop("Step 2 command script did not shell-quote the regenie tool path")

empty_traits <- file.path(tmp, "empty.traits.txt")
write_lines(character(), empty_traits)
noop_step1 <- file.path(tmp, "noop_step1.sh")
noop_pred <- file.path(tmp, "noop_pred.list")
run_phase2(c(
  "write-step1-command", "--config", cmd_config, "--group", bt_group, "--pfile-prefix", file.path(tmp, "step1"),
  "--extract", file.path(tmp, "extract.txt"), "--pheno", file.path(tmp, "pheno.tsv"),
  "--covar", file.path(tmp, "covar.tsv"), "--keep", file.path(tmp, "keep.txt"),
  "--trait-list", empty_traits, "--pred-list", noop_pred,
  "--out-prefix", file.path(tmp, "noop_step1"), "--script-out", noop_step1, "--threads", "2"
))
status <- system2("bash", noop_step1, stdout = TRUE, stderr = TRUE)
if (!identical(as.integer(attr(status, "status") %||% 0L), 0L)) stop("Step 1 no-op script failed")
if (!file.exists(noop_pred) || file.info(noop_pred)$size != 0) stop("Step 1 no-op script did not create an empty prediction list")
if (!identical(readLines(file.path(tmp, "noop_step1.done")), "skipped_no_traits")) stop("Step 1 no-op script did not write skip sentinel")

noop_step2 <- file.path(tmp, "noop_step2.sh")
noop_done <- file.path(tmp, "noop_step2.done")
run_phase2(c(
  "write-step2-command", "--config", cmd_config, "--group", bt_group, "--pfile-prefix", file.path(tmp, "assoc"),
  "--pheno", file.path(tmp, "pheno.tsv"), "--covar", file.path(tmp, "covar.tsv"),
  "--pred-list", file.path(tmp, "pred.list"), "--trait-list", empty_traits,
  "--out-prefix", file.path(tmp, "noop_step2"), "--done", noop_done,
  "--script-out", noop_step2, "--threads", "2"
))
status <- system2("bash", noop_step2, stdout = TRUE, stderr = TRUE)
if (!identical(as.integer(attr(status, "status") %||% 0L), 0L)) stop("Step 2 no-op script failed")
if (!identical(readLines(noop_done), "skipped_no_traits")) stop("Step 2 no-op script did not write skip sentinel")

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
