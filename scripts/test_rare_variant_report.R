#!/usr/bin/env Rscript

# Focused tests for active and skipped rare-variant reporting.

cmd <- commandArgs(FALSE)
script_path <- sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])
repo <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(repo, "scripts", "lib", "stage1.R"))

tmp <- tempfile("rare-report-")
dir.create(tmp, recursive = TRUE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

write_lines <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path)
  path
}

run_report <- function(values, expect_success = TRUE) {
  output <- suppressWarnings(system2(
    "Rscript", c(file.path(repo, "scripts", "make_rare_variant_report.R"), values),
    stdout = TRUE, stderr = TRUE
  ))
  success <- is.null(attr(output, "status"))
  if (success != expect_success) stop(paste(output, collapse = "\n"))
  output
}

config <- write_lines(c(
  "project:", "  analysis_name: rare_report_test"
), file.path(tmp, "config.yaml"))
summary_header <- paste(c(
  "group", "trait", "trait_type", "covariates", "usable_n", "cases", "controls",
  "model_sample_count", "model_cases", "model_controls", "skipped", "skip_reason"
), collapse = "\t")
active_row <- paste(c(
  "bt__BT", "BT", "bt", "age,sex,PC1", 30, 10, 20, 30, 10, 20, "False", ""
), collapse = "\t")
summary <- write_lines(c(summary_header, active_row), file.path(tmp, "summary.tsv"))
group_summary <- write_lines(c(summary_header, active_row), file.path(tmp, "group_summary.tsv"))
htp <- write_lines("validated HTP", file.path(tmp, "BT.PAN.regenie.gz"))
target <- write_lines(c(
  "key\tvalue", "data_source\twes", "genotype_mode\thardcall", "sample_count\t30",
  "variant_count\t5", "min_mac\t1", "geno_missing_max\t0.05", "info_min\tNA",
  "target_regions\tresources/remeta/GRCh38/target_regions.bed"
), file.path(tmp, "target.tsv"))
validation <- write_lines(c(
  "key\tvalue", "status\tvalidated", "group\tbt__BT", "trait\tBT", "sample_count\t30",
  "ordinary_regenie_sample_count\t30", "rare_regenie_sample_count\t30",
  paste0("htp_sha256\t", sha256_file(htp)), "htp_variant_rows\t5", "target_variant_count\t5",
  "unique_ld_target_variant_count\t4", "target_variants_not_indexed\t1",
  "target_variant_ld_coverage_pct\t80", "reference_gene_count\t3", "indexed_gene_count\t2",
  "genes_without_indexed_variants\t1", "indexed_gene_coverage_pct\t66.666667",
  "ld_gene_variant_assignments\t6", "ld_assignments_within_gene_bounds\t5",
  "ld_assignments_outside_gene_bounds\t1", "ld_assignment_gene_bound_coverage_pct\t83.333333"
), file.path(tmp, "validation.tsv"))
metrics <- write_lines(c(
  "metric\tvalue", "total_variants\t5", "valid_p_value_variants\t5", "lambda_gc\t0.9",
  "genomewide_significant_variants\t0", "suggestive_variants\t1", "plot_thinning_applied\tFalse",
  "genotype_mode\thardcall", "mac_definition\texact hardcall MAC", "mac_bin_1\t1",
  "mac_bin_2_5\t1", "mac_bin_6_10\t1", "mac_bin_11_20\t1", "mac_bin_gt20\t1",
  "qq_eligible_variants\t5", "qq_points_plotted\t5",
  "manhattan_eligible_variants\t5", "manhattan_points_plotted\t5",
  "mac_qq_eligible_variants\t5", "mac_qq_points_plotted\t5",
  "effect_frequency_eligible_variants\t5", "effect_frequency_points_plotted\t5"
), file.path(tmp, "metrics.tsv"))
hits <- write_lines(c(
  "chrom\tpos\tcpra\teffect\tlci\tuci\taaf\tmaf\tmac\tp",
  "1\t100\t1:100:A:G\t1.5\t1.1\t2.0\t0.01\t0.01\t1\t1e-6"
), file.path(tmp, "hits.tsv"))
plots <- setNames(file.path(tmp, paste0(c("qq", "manhattan", "manhattan_pdf", "mac_qq", "effect_frequency"), ".png")),
  c("qq", "manhattan", "manhattan-pdf", "mac-qq", "effect-frequency"))
invisible(vapply(plots, function(path) write_lines("plot", path), character(1)))

report <- file.path(tmp, "active.report.md")
active_values <- c(
  "--config", config, "--trait", "BT", "--build", "GRCh38", "--summary", summary,
  "--group-summary", group_summary, "--htp", htp, "--target-summary", target,
  "--validation", validation, "--metrics", metrics, "--top-hits", hits,
  unlist(Map(function(flag, path) c(paste0("--", flag), path), names(plots), plots), use.names = FALSE),
  "--out", report
)
invisible(run_report(active_values))
text <- readLines(report, warn = FALSE)
stopifnot(
  any(grepl("does not contain a local burden, SKAT, mask, gene, or", text, fixed = TRUE)),
  any(grepl("MAC-stratified QQ plot", text, fixed = TRUE)),
  any(grepl("Effect versus frequency plot", text, fixed = TRUE)),
  any(grepl("QC-passing target-region variants represented in LD indexes | 4 | 5 | 80.00%", text, fixed = TRUE)),
  any(grepl("all 22 autosomes passed", text, fixed = TRUE)),
  any(grepl("Conditional buffer variants: not included", text, fixed = TRUE)),
  any(grepl("ALT OR", text, fixed = TRUE))
)

# Metrics and validation must agree before a report can be emitted.
bad_metrics <- write_lines(sub("total_variants\\t5", "total_variants\\t4", readLines(metrics)),
  file.path(tmp, "bad_metrics.tsv"))
bad_values <- active_values
bad_values[match(metrics, bad_values)] <- bad_metrics
bad_values[match(report, bad_values)] <- file.path(tmp, "bad.report.md")
invisible(run_report(bad_values, expect_success = FALSE))

# Skipped traits receive one text-only audit artifact and no plot arguments.
skipped_row <- paste(c(
  "bt__SKIP", "SKIP", "bt", "age,sex,PC1", 200, 20, 180, 0, 0, 0,
  "True", "n=200<min_n=300; cases=20<min_cases=50"
), collapse = "\t")
skipped_summary <- write_lines(c(summary_header, skipped_row), file.path(tmp, "skipped.summary.tsv"))
skipped_group <- write_lines(c(summary_header, skipped_row), file.path(tmp, "skipped.group.tsv"))
skipped_report <- file.path(tmp, "skipped.report.md")
invisible(run_report(c(
  "--config", config, "--trait", "SKIP", "--build", "GRCh38",
  "--summary", skipped_summary, "--group-summary", skipped_group, "--out", skipped_report
)))
skipped_text <- readLines(skipped_report, warn = FALSE)
stopifnot(
  any(grepl("Export status: skipped", skipped_text, fixed = TRUE)),
  any(grepl("n=200<min_n=300; cases=20<min_cases=50", skipped_text, fixed = TRUE)),
  any(grepl("No rare-variant HTP file, target PGEN, LD matrices, metrics, or plots", skipped_text, fixed = TRUE)),
  !any(grepl("!\\[", skipped_text))
)

cat("Rare-variant report tests passed\n")
