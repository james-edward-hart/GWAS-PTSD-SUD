#!/usr/bin/env Rscript

# Focused tests for large-file-safe GWAS plotting and Phase 2 summary outputs.

cmd <- commandArgs(FALSE)
script_path <- sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])
repo <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(repo, "scripts", "lib", "stage1.R"))

tmp <- tempfile("plot-gwas-")
dir.create(tmp, recursive = TRUE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

write_lines <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path)
}

run_plot <- function(stats, name, threshold = "Inf") {
  stem <- file.path(tmp, name)
  outputs <- list(
    qq = paste0(stem, ".qq.png"),
    manhattan = paste0(stem, ".manhattan.png"),
    pdf = paste0(stem, ".manhattan.pdf"),
    metrics = paste0(stem, ".metrics.tsv"),
    hits = paste0(stem, ".hits.tsv")
  )
  command <- c(
    file.path(repo, "scripts", "plot_gwas.R"),
    "--stats", stats,
    "--qq", outputs$qq,
    "--manhattan", outputs$manhattan,
    "--manhattan-pdf", outputs$pdf,
    "--metrics-out", outputs$metrics,
    "--top-hits-out", outputs$hits,
    "--large-plot-threshold", threshold
  )
  log <- suppressWarnings(system2("Rscript", command, stdout = TRUE, stderr = TRUE))
  status <- attr(log, "status") %||% 0L
  if (!identical(as.integer(status), 0L)) {
    stop("plot_gwas.R failed:\n", paste(log, collapse = "\n"), call. = FALSE)
  }
  for (path in unlist(outputs, use.names = FALSE)) {
    if (!file.exists(path) || file.info(path)$size <= 0) stop("missing plotting output: ", path)
  }
  c(outputs, list(log = log))
}

metric <- function(path, key) {
  rows <- read_tsv(path)
  value <- rows$value[rows$metric == key]
  if (length(value) != 1L) stop("metric is missing or duplicated: ", key)
  as.character(value[[1]])
}

# Tab-delimited HTP input: metadata is skipped, invalid P values are excluded,
# and exact metrics/top hits use all valid variants.
tab_stats <- file.path(tmp, "tab.regenie")
write_lines(c(
  "##source=fixture",
  "Name\tChr\tPos\tTrait\tEffect\tPval\tInfo",
  "rs1\t1\t100\ttrait_a\t0.5\t1e-9\tSCORE=1",
  "rs2\t1\t200\ttrait_a\t0.2\t0.02\tSCORE=1",
  "rs3\t2\t300\ttrait_a\t-0.3\t1e-6\tSCORE=1",
  "rs4\t2\t400\ttrait_a\t0.1\t0\tSCORE=1",
  "rs5\t3\t500\ttrait_a\t0.1\t2\tSCORE=1",
  "rs6\t3\t600\ttrait_a\t0.1\tNA\tSCORE=1"
), tab_stats)
tab <- run_plot(tab_stats, "tab")
stopifnot(
  identical(metric(tab$metrics, "total_variants"), "6"),
  identical(metric(tab$metrics, "valid_p_value_variants"), "3"),
  identical(metric(tab$metrics, "genomewide_significant_variants"), "1"),
  identical(metric(tab$metrics, "suggestive_variants"), "2"),
  identical(metric(tab$metrics, "qq_points_plotted"), "3"),
  identical(metric(tab$metrics, "manhattan_points_plotted"), "3"),
  identical(metric(tab$metrics, "plot_thinning_applied"), "False")
)
tab_hits <- read_tsv(tab$hits)
stopifnot(identical(tab_hits$variant_id, c("rs1", "rs3", "rs2")))

# Whitespace-delimited LOG10P input exercises aliases and the built-in large
# plot fallback. All six rows still contribute to exact metrics and top hits.
log_stats <- file.path(tmp, "log10p.regenie")
write_lines(c(
  "ID CHROM GENPOS BETA SE LOG10P",
  "rs1 1 100 0.5 0.1 9",
  "rs2 1 200 0.4 0.1 6",
  "rs3 2 300 0.3 0.1 4",
  "rs4 2 400 0.2 0.1 3",
  "rs5 3 500 0.1 0.1 2",
  "rs6 3 600 0.0 0.1 1"
), log_stats)
fallback <- run_plot(log_stats, "fallback", threshold = "4")
stopifnot(
  identical(metric(fallback$metrics, "total_variants"), "6"),
  identical(metric(fallback$metrics, "valid_p_value_variants"), "6"),
  identical(metric(fallback$metrics, "qq_points_plotted"), "5"),
  identical(metric(fallback$metrics, "manhattan_points_plotted"), "3"),
  identical(metric(fallback$metrics, "plot_thinning_applied"), "True")
)
fallback_hits <- read_tsv(fallback$hits)
stopifnot(
  identical(fallback_hits$variant_id, paste0("rs", 1:6)),
  any(grepl("Applied large-plot fallback", fallback$log, fixed = TRUE))
)

# A skipped Phase 2 trait contains metadata plus a header but no association
# rows. It must still produce valid empty plots and machine-readable summaries.
skipped_stats <- file.path(tmp, "skipped.regenie")
write_lines(c(
  "## skipped: below_phase2_thresholds",
  "Name\tChr\tPos\tRef\tAlt\tCohort\tModel\tEffect\tLCI_effect\tUCI_effect\tPval\tAAF\tNum_Cases\tCases_Ref\tCases_Het\tCases_Alt\tNum_Controls\tControls_Ref\tControls_Het\tControls_Alt\tInfo"
), skipped_stats)
skipped <- run_plot(skipped_stats, "skipped", threshold = "4")
stopifnot(
  identical(metric(skipped$metrics, "total_variants"), "0"),
  identical(metric(skipped$metrics, "valid_p_value_variants"), "0"),
  identical(metric(skipped$metrics, "plot_thinning_applied"), "False"),
  nrow(read_tsv(skipped$hits)) == 0L
)

cat("GWAS plotting tests passed\n")
