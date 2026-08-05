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
  identical(metric(fallback$metrics, "qq_points_plotted"), "4"),
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


write_gzip <- function(lines, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con))
  writeLines(lines, con)
  path
}


run_rare_plot <- function(stats, name, config, summary, threshold = "Inf") {
  stem <- file.path(tmp, name)
  outputs <- list(
    qq = paste0(stem, ".qq.png"),
    manhattan = paste0(stem, ".manhattan.png"),
    pdf = paste0(stem, ".manhattan.pdf"),
    mac_qq = paste0(stem, ".mac_qq.png"),
    effect_frequency = paste0(stem, ".effect_frequency.png"),
    metrics = paste0(stem, ".metrics.tsv"),
    hits = paste0(stem, ".hits.tsv")
  )
  command <- c(
    file.path(repo, "scripts", "plot_gwas.R"),
    "--config", config, "--stats", stats, "--trait-summary", summary,
    "--qq", outputs$qq, "--manhattan", outputs$manhattan,
    "--manhattan-pdf", outputs$pdf, "--mac-qq", outputs$mac_qq,
    "--effect-frequency", outputs$effect_frequency,
    "--metrics-out", outputs$metrics, "--top-hits-out", outputs$hits,
    "--large-plot-threshold", threshold
  )
  log <- suppressWarnings(system2("Rscript", command, stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(log, "status"))) stop("rare plotting failed:\n", paste(log, collapse = "\n"))
  for (path in unlist(outputs, use.names = FALSE)) {
    if (!file.exists(path) || file.info(path)$size <= 0) stop("missing rare plotting output: ", path)
  }
  outputs
}


htp_header <- paste(c(
  "Name", "Chr", "Pos", "Ref", "Alt", "Trait", "Cohort", "Model", "Effect",
  "LCI_Effect", "UCI_Effect", "Pval", "AAF", "Num_Cases", "Cases_Ref",
  "Cases_Het", "Cases_Alt", "Num_Controls", "Controls_Ref", "Controls_Het",
  "Controls_Alt", "Info"
), collapse = "\t")
htp_row <- function(id, pos, trait, effect, p, aaf, cases, controls) paste(c(
  id, "1", pos, "A", "G", trait, "cohort", "ADD", effect, effect * 0.8,
  effect * 1.2, p, aaf, sum(cases), cases, sum(controls), controls, "SCORE=1"
), collapse = "\t")

# Hardcall HTP rows exercise exact MAC bins, gzip input, top-hit fields, and
# deterministic fallback without changing exact summaries.
hardcall_config <- file.path(tmp, "hardcall.yaml")
write_lines(c("project:", "  analysis_name: rare_test", "remeta:", "  genotype_mode: hardcall"), hardcall_config)
hardcall_summary <- file.path(tmp, "hardcall.summary.tsv")
write_lines(c(
  "group\ttrait\ttrait_type\tskipped",
  "bt__BT\tBT\tbt\tFalse"
), hardcall_summary)
hardcall_stats <- write_gzip(c(
  "##source=compressed_fixture", htp_header,
  htp_row("1:101:A:G", 101, "BT", 1.5, 1e-8, 1 / 60, c(14, 1, 0), c(15, 0, 0)),
  htp_row("1:102:A:G", 102, "BT", 1.4, 1e-6, 3 / 60, c(14, 0, 1), c(14, 1, 0)),
  htp_row("1:103:A:G", 103, "BT", 1.3, 1e-4, 8 / 60, c(13, 0, 2), c(13, 0, 2)),
  htp_row("1:104:A:G", 104, "BT", 1.2, 0.01, 15 / 60, c(10, 0, 5), c(12, 1, 2)),
  htp_row("1:105:A:G", 105, "BT", 1.1, 0.2, 0.5, c(7, 1, 7), c(7, 1, 7)),
  htp_row("1:106:A:G", 106, "BT", 1.0, 0, 2 / 60, c(13, 2, 0), c(15, 0, 0))
), file.path(tmp, "hardcall.regenie.gz"))
hardcall <- run_rare_plot(hardcall_stats, "hardcall", hardcall_config, hardcall_summary, "3")
stopifnot(
  metric(hardcall$metrics, "total_variants") == "6",
  metric(hardcall$metrics, "valid_p_value_variants") == "5",
  metric(hardcall$metrics, "mac_definition") == "exact hardcall MAC",
  all(vapply(c("mac_bin_1", "mac_bin_2_5", "mac_bin_6_10", "mac_bin_11_20", "mac_bin_gt20"),
    function(key) metric(hardcall$metrics, key) == "1", logical(1))),
  metric(hardcall$metrics, "qq_points_plotted") == "4",
  metric(hardcall$metrics, "manhattan_points_plotted") == "3",
  metric(hardcall$metrics, "mac_qq_points_plotted") == "5",
  metric(hardcall$metrics, "effect_frequency_points_plotted") == "4",
  metric(hardcall$metrics, "plot_thinning_applied") == "True"
)
hardcall_hits <- read_tsv(hardcall$hits)
stopifnot(
  identical(hardcall_hits$cpra, paste0("1:10", 1:5, ":A:G")),
  all(c("lci", "uci", "aaf", "maf", "mac") %in% names(hardcall_hits)),
  identical(as.numeric(hardcall_hits$mac), c(1, 3, 8, 15, 30))
)

# QT HTP uses Num_Cases as N; dosage MAC remains an expected, potentially
# non-integer count derived from the reported AAF.
dosage_config <- file.path(tmp, "dosage.yaml")
write_lines(c("project:", "  analysis_name: rare_test", "remeta:", "  genotype_mode: dosage"), dosage_config)
dosage_summary <- file.path(tmp, "dosage.summary.tsv")
write_lines(c("group\ttrait\ttrait_type\tskipped", "qt__QT\tQT\tqt\tFalse"), dosage_summary)
qt_row <- function(id, pos, effect, p, aaf) paste(c(
  id, "2", pos, "C", "T", "QT", "cohort", "ADD", effect, effect - 0.1,
  effect + 0.1, p, aaf, 10, 8, 2, 0, rep("NA", 4), "SCORE=1"
), collapse = "\t")
dosage_stats <- write_gzip(c(
  htp_header,
  qt_row("2:201:C:T", 201, -0.2, 0.001, 0.05),
  qt_row("2:202:C:T", 202, 0.3, 0.01, 0.15),
  qt_row("2:203:C:T", 203, 0.1, 0.1, 0.4)
), file.path(tmp, "dosage.regenie.gz"))
dosage <- run_rare_plot(dosage_stats, "dosage", dosage_config, dosage_summary)
dosage_hits <- read_tsv(dosage$hits)
stopifnot(
  metric(dosage$metrics, "mac_definition") == "dosage-based expected MAC",
  identical(as.numeric(dosage_hits$mac), c(1, 3, 8)),
  metric(dosage$metrics, "mac_bin_1") == "1",
  metric(dosage$metrics, "mac_bin_2_5") == "1",
  metric(dosage$metrics, "mac_bin_6_10") == "1"
)

# A validated active export may legitimately contain only the HTP header.
empty_rare <- run_rare_plot(
  write_gzip(htp_header, file.path(tmp, "empty.regenie.gz")),
  "empty_rare", hardcall_config, hardcall_summary
)
stopifnot(
  metric(empty_rare$metrics, "total_variants") == "0",
  metric(empty_rare$metrics, "valid_p_value_variants") == "0",
  nrow(read_tsv(empty_rare$hits)) == 0L
)

cat("GWAS plotting tests passed\n")
