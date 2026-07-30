#!/usr/bin/env Rscript

# Focused tests for cohort ReMeta invariants and command generation.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
remeta_script <- file.path(script_dir, "remeta_cohort.R")
rscript <- file.path(R.home("bin"), "Rscript")
tmp <- tempfile("test_remeta_cohort.")
dir.create(tmp, recursive = TRUE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)


write_lines <- function(lines, name) {
  path <- file.path(tmp, name)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path)
  path
}


write_gzip <- function(lines, name) {
  path <- file.path(tmp, name)
  con <- gzfile(path, "wt")
  writeLines(lines, con)
  close(con)
  path
}


run_task <- function(values) {
  system2(rscript, c(remeta_script, values), stdout = TRUE, stderr = TRUE)
}


htp_columns <- c(
  "Name", "Chr", "Pos", "Ref", "Alt", "Trait", "Cohort", "Model", "Effect", "LCI_Effect",
  "UCI_Effect", "Pval", "AAF", "Num_Cases", "Cases_Ref", "Cases_Het", "Cases_Alt",
  "Num_Controls", "Controls_Ref", "Controls_Het", "Controls_Alt", "Info"
)


htp_row <- function(name, chrom, pos, ref, alt, info) {
  paste(c(name, chrom, pos, ref, alt, "TRAIT1", "cohort_test", "ADD", "1.1", "1.0", "1.2",
    "0.5", "0.1", "10", "8", "2", "0", "10", "8", "2", "0", info), collapse = "\t")
}


target <- file.path(tmp, "target")
invisible(write_lines(c("#FID\tIID", "F1\tI1", "F2\tI2"), "target.psam"))
invisible(write_lines(c(
  "#CHROM\tPOS\tID\tREF\tALT",
  "1\t100\t1:100:A:G\tA\tG",
  "1\t120\t1:120:C:T\tC\tT"
), "target.pvar"))
invisible(write_lines(c("F1\tI1", "F2\tI2"), "keep.txt"))
gene_list <- write_lines(c(
  "ENSG000001\t1\t90\t130",
  "ENSG000002\t1\t200\t250"
), "gene_list.tsv")
htp <- write_gzip(c(
  paste(htp_columns, collapse = "\t"),
  htp_row("1:100:A:G", "1", "100", "A", "G", "SCORE=0.2;SKATV=0.1")
), "trait.regenie.gz")
index <- write_gzip("ENSG000001\t0\t0\t1:100:A:G,1:120:C:T\t", "chr1.remeta.ld.idx.gz")
ok <- file.path(tmp, "validation.ok")

output <- run_task(c(
  "validate-group", "--target-prefix", target, "--keep", file.path(tmp, "keep.txt"),
  "--gene-list", gene_list, "--htp", htp, "--index", index, "--out", ok
))
stopifnot(is.null(attr(output, "status")), file.exists(ok))
validation <- read.delim(ok, stringsAsFactors = FALSE)
stopifnot(validation$value[validation$key == "sample_count"] == "2")
stopifnot(validation$value[validation$key == "target_variant_count"] == "2")
stopifnot(validation$value[validation$key == "unique_ld_target_variant_count"] == "2")
stopifnot(validation$value[validation$key == "target_variant_ld_coverage_pct"] == "100.000000")
stopifnot(validation$value[validation$key == "reference_gene_count"] == "2")
stopifnot(validation$value[validation$key == "indexed_gene_count"] == "1")
stopifnot(validation$value[validation$key == "indexed_gene_coverage_pct"] == "50.000000")
stopifnot(validation$value[validation$key == "ld_gene_variant_assignments"] == "2")
stopifnot(validation$value[validation$key == "ld_assignments_within_gene_bounds"] == "2")
stopifnot(validation$value[validation$key == "ld_assignment_gene_bound_coverage_pct"] == "100.000000")


# A score-statistic variant without matching LD target data must fail.
bad_htp <- write_gzip(c(
  paste(htp_columns, collapse = "\t"),
  htp_row("1:999:A:G", "1", "999", "A", "G", "SCORE=0.2;SKATV=0.1")
), "bad.regenie.gz")
bad <- suppressWarnings(run_task(c(
  "validate-group", "--target-prefix", target, "--keep", file.path(tmp, "keep.txt"),
  "--gene-list", gene_list, "--htp", bad_htp, "--index", index,
  "--out", file.path(tmp, "bad.ok")
)))
stopifnot(!is.null(attr(bad, "status")), attr(bad, "status") != 0)

# A target-PVAR variant without an LD-index entry must also fail.
missing_ld_htp <- write_gzip(c(
  paste(htp_columns, collapse = "\t"),
  htp_row("1:120:C:T", "1", "120", "C", "T", "SCORE=0.2;SKATV=0.1")
), "missing_ld.regenie.gz")
index_one <- write_gzip("ENSG000001\t0\t0\t1:100:A:G\t", "chr1.one.remeta.ld.idx.gz")
missing_ld <- suppressWarnings(run_task(c(
  "validate-group", "--target-prefix", target, "--keep", file.path(tmp, "keep.txt"),
  "--gene-list", gene_list, "--htp", missing_ld_htp, "--index", index_one,
  "--out", file.path(tmp, "missing_ld.ok")
)))
stopifnot(!is.null(attr(missing_ld, "status")), attr(missing_ld, "status") != 0)

# A CPRA identifier is insufficient when PLINK still marks REF as provisional.
provisional <- file.path(tmp, "provisional")
invisible(file.copy(file.path(tmp, "target.psam"), paste0(provisional, ".psam")))
invisible(write_lines(c(
  "#CHROM\tPOS\tID\tREF\tALT\tINFO",
  "1\t100\t1:100:A:G\tA\tG\tPR"
), "provisional.pvar"))
bad_ref <- suppressWarnings(run_task(c(
  "validate-group", "--target-prefix", provisional, "--keep", file.path(tmp, "keep.txt"),
  "--gene-list", gene_list, "--htp", htp, "--index", index,
  "--out", file.path(tmp, "bad_ref.ok")
)))
stopifnot(!is.null(attr(bad_ref, "status")), attr(bad_ref, "status") != 0)


config <- write_lines(c(
  "project:",
  "  analysis_name: cohort_test",
  "  cohort_data_release: test_release",
  "phase2_regenie:",
  "  htp_cohort_name: cohort_test",
  "  step2_bsize: 400",
  "  p_thresh: 0.01",
  "  apply_rint: false",
  "  step2_options: ''",
  "remeta:",
  "  data_source: wes",
  "  genotype_mode: hardcall",
  "  input_variants_normalized: true",
  "  min_mac: 1",
  "  info_min: 0.8",
  "  target_r2: 0.0001",
  "tools:",
  "  regenie: regenie"
), "config.yaml")
summary <- write_lines(c(
  paste(c("group", "trait", "trait_type", "phase2_pan_samples", "usable_n", "skipped", "skip_reason"), collapse = "\t"),
  "bt_g1\tTRAIT1\tbt\t100\t95\tFalse\t"
), "summary.tsv")
trait_list <- write_lines("TRAIT1", "traits.txt")
covar_list <- write_lines("age,sex,PC1", "covars.txt")
command <- file.path(tmp, "step2.sh")
done <- file.path(tmp, "step2.done")
output <- run_task(c(
  "write-step2-command", "--config", config, "--group-summary", summary,
  "--pfile-prefix", target, "--pheno", "pheno.tsv", "--covar", "covar.tsv",
  "--pred-list", "step1_pred.list", "--trait-list", trait_list,
  "--covar-list", covar_list, "--out-prefix", file.path(tmp, "rare"),
  "--done", done, "--script-out", command, "--threads", "4"
))
stopifnot(is.null(attr(output, "status")), file.exists(command))
text <- paste(readLines(command), collapse = "\n")
stopifnot(grepl("--gz", text, fixed = TRUE))
stopifnot(grepl("--minMAC.*'1'", text))
stopifnot(grepl("--htp.*'cohort_test'", text))
stopifnot(!grepl("--minINFO", text, fixed = TRUE))

skipped_summary <- write_lines(c(
  paste(c("group", "trait", "trait_type", "phase2_pan_samples", "usable_n", "skipped", "skip_reason"), collapse = "\t"),
  "bt_g1\tSKIPPED\tbt\t100\t10\tTrue\tbelow_threshold"
), "skipped_summary.tsv")
skipped_htp <- file.path(tmp, "skipped.regenie.gz")
output <- run_task(c(
  "stage-trait", "--trait", "SKIPPED", "--group-summary", skipped_summary,
  "--raw-prefix", file.path(tmp, "absent"), "--out", skipped_htp
))
stopifnot(is.null(attr(output, "status")), file.exists(skipped_htp))
skipped_con <- gzfile(skipped_htp, "rt")
skipped_lines <- readLines(skipped_con, warn = FALSE)
close(skipped_con)
stopifnot(startsWith(skipped_lines[[1]], "## skipped:"))
stopifnot("Trait" %in% strsplit(skipped_lines[[2]], "\t", fixed = TRUE)[[1]])

dosage_step2_config <- write_lines(c(
  readLines(config),
  ""
), "dosage_step2_config.yaml")
dosage_lines <- readLines(dosage_step2_config)
dosage_lines[dosage_lines == "  data_source: wes"] <- "  data_source: imputed"
dosage_lines[dosage_lines == "  genotype_mode: hardcall"] <- "  genotype_mode: dosage"
writeLines(dosage_lines, dosage_step2_config)
dosage_command <- file.path(tmp, "dosage_step2.sh")
output <- run_task(c(
  "write-step2-command", "--config", dosage_step2_config, "--group-summary", summary,
  "--pfile-prefix", target, "--pheno", "pheno.tsv", "--covar", "covar.tsv",
  "--pred-list", "step1_pred.list", "--trait-list", trait_list,
  "--covar-list", covar_list, "--out-prefix", file.path(tmp, "rare_dosage"),
  "--done", file.path(tmp, "dosage.done"), "--script-out", dosage_command, "--threads", "4"
))
stopifnot(is.null(attr(output, "status")))
stopifnot(grepl("--minINFO.*'0.8'", paste(readLines(dosage_command), collapse = "\n")))


# Capture the target-preparation PLINK arguments and return a valid tiny PGEN
# trio. This checks the coordinate convention and literal CPRA template.
captured <- file.path(tmp, "plink_args.txt")
fake_plink <- write_lines(c(
  "#!/usr/bin/env bash",
  "set -euo pipefail",
  paste("printf '%s\\n' \"$@\" >", shQuote(captured)),
  "out=''",
  "while (($#)); do",
  "  if [[ \"$1\" == '--out' ]]; then out=\"$2\"; shift 2; else shift; fi",
  "done",
  "printf '#FID\\tIID\\nF1\\tI1\\nF2\\tI2\\n' > \"${out}.psam\"",
  "printf '#CHROM\\tPOS\\tID\\tREF\\tALT\\n1\\t100\\t1:100:A:G\\tA\\tG\\n' > \"${out}.pvar\"",
  ": > \"${out}.pgen\""
), "fake_plink.sh")
Sys.chmod(fake_plink, "0755")
prepare_config <- write_lines(c(
  "remeta:",
  "  data_source: wes",
  "  genotype_mode: hardcall",
  "  input_variants_normalized: true",
  "  min_mac: 1",
  "  geno_missing_max: 0.05",
  "tools:",
  paste0("  plink2: ", fake_plink)
), "prepare_config.yaml")
regions <- write_lines("1\t90\t130", "regions.bed")
prepared <- file.path(tmp, "prepared")
output <- run_task(c(
  "prepare-target", "--config", prepare_config, "--pfile-prefix", "source",
  "--keep", file.path(tmp, "keep.txt"), "--regions", regions,
  "--out-prefix", prepared, "--summary-out", file.path(tmp, "prepared.summary.tsv"),
  "--threads", "2"
))
stopifnot(is.null(attr(output, "status")))
plink_args <- readLines(captured)
stopifnot(any(plink_args == "bed0"))
stopifnot(any(plink_args == "26"))
stopifnot(any(plink_args == "@:#:$r:$a"))
stopifnot(any(plink_args == "erase-dosage"))

dosage_config <- write_lines(c(
  "remeta:",
  "  data_source: imputed",
  "  genotype_mode: dosage",
  "  input_variants_normalized: true",
  "  min_mac: 1",
  "  geno_missing_max: 0.05",
  "  info_min: 0.8",
  "tools:",
  paste0("  plink2: ", fake_plink)
), "dosage_config.yaml")
output <- run_task(c(
  "prepare-target", "--config", dosage_config, "--pfile-prefix", "source",
  "--keep", file.path(tmp, "keep.txt"), "--regions", regions,
  "--out-prefix", file.path(tmp, "prepared_dosage"),
  "--summary-out", file.path(tmp, "prepared_dosage.summary.tsv"), "--threads", "2"
))
stopifnot(is.null(attr(output, "status")))
plink_args <- readLines(captured)
stopifnot(any(plink_args == "dosage"))
stopifnot(any(plink_args == "--mach-r2-filter"))
stopifnot(!any(plink_args == "erase-dosage"))


provenance <- write_lines(c("key\tvalue", "genome_build\tGRCh38"), "provenance.tsv")
tool <- write_lines(c("key\tvalue", "tool_remeta_version\t0.11.2"), "tool.tsv")
manifest <- file.path(tmp, "manifest.tsv")
output <- run_task(c(
  "write-manifest", "--config", config, "--build", "GRCh38",
  "--gene-list", gene_list, "--provenance", provenance,
  "--target-summary", file.path(tmp, "prepared.summary.tsv"),
  "--trait-summary", summary,
  "--validation", ok, "--tool", tool, "--artifact", htp, index,
  "--out", manifest
))
stopifnot(is.null(attr(output, "status")), file.exists(manifest))
manifest_rows <- read.delim(manifest, stringsAsFactors = FALSE)
stopifnot(manifest_rows$value[manifest_rows$key == "analysis_scope"] == "marginal_gene_tests_only")
stopifnot(manifest_rows$value[manifest_rows$key == "conditional_buffer_included"] == "false")
stopifnot(manifest_rows$value[manifest_rows$key == "trait:TRAIT1:ld_group"] == "bt_g1")
stopifnot(any(grepl(":sha256$", manifest_rows$key)))

cat("ReMeta cohort helper tests passed\n")
