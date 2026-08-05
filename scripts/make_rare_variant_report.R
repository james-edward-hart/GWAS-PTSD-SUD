#!/usr/bin/env Rscript

# Report existing target-region REGENIE statistics and validated ReMeta LD.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


active_args <- c(
  "htp", "target-summary", "validation", "metrics", "top-hits", "qq", "manhattan",
  "manhattan-pdf", "mac-qq", "effect-frequency"
)
args <- parse_args(defaults = as.list(setNames(rep("", length(active_args)), active_args)))
require_args(args, c("config", "trait", "build", "summary", "group-summary", "out"))


value_for <- function(rows, key, key_column = "key", default = "NA") {
  hit <- rows$value[rows[[key_column]] == key]
  if (length(hit) == 1L) as.character(hit[[1]]) else default
}


require_values <- function(rows, keys, key_column, label) {
  missing <- setdiff(keys, as.character(rows[[key_column]]))
  if (length(missing)) die(label, " is missing key(s): ", paste(missing, collapse = ", "))
}


format_number <- function(value, digits = 5L) {
  numeric_value <- suppressWarnings(as.numeric(value))
  if (!length(numeric_value) || !is.finite(numeric_value[[1]])) return("NA")
  format(signif(numeric_value[[1]], digits), scientific = TRUE, trim = TRUE)
}


format_percent <- function(value) {
  numeric_value <- suppressWarnings(as.numeric(value))
  if (!length(numeric_value) || !is.finite(numeric_value[[1]])) return("NA")
  sprintf("%.2f%%", numeric_value[[1]])
}


top_hit_lines <- function(rows, trait_type, mac_definition) {
  required <- c("chrom", "pos", "cpra", "effect", "lci", "uci", "aaf", "maf", "mac", "p")
  require_columns(rows, required, "rare-variant top hits")
  if (!nrow(rows)) return("No valid single-variant P values were available.")
  effect_label <- if (identical(trait_type, "bt")) "ALT OR" else "ALT beta"
  mac_label <- if (grepl("expected", mac_definition, fixed = TRUE)) "Expected MAC" else "MAC"
  c(
    paste0("| CHROM | POS | CPRA | ", effect_label, " | 95% CI | AAF | MAF | ", mac_label, " | P |"),
    "| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    vapply(seq_len(nrow(rows)), function(i) paste0(
      "| ", rows$chrom[[i]], " | ", rows$pos[[i]], " | ", rows$cpra[[i]], " | ",
      format_number(rows$effect[[i]]), " | ", format_number(rows$lci[[i]]), "–",
      format_number(rows$uci[[i]]), " | ", format_number(rows$aaf[[i]]), " | ",
      format_number(rows$maf[[i]]), " | ", format_number(rows$mac[[i]]), " | ",
      format_number(rows$p[[i]]), " |"
    ), character(1))
  )
}


config <- load_config(args$config)
summary <- read_tsv(args$summary)
group_summary <- read_tsv(args[["group-summary"]])
if (nrow(summary) != 1L || nrow(group_summary) != 1L ||
    summary$trait[[1]] != args$trait || group_summary$trait[[1]] != args$trait ||
    summary$group[[1]] != group_summary$group[[1]] ||
    summary$skipped[[1]] != group_summary$skipped[[1]]) {
  die("rare-variant report requires matching singleton public and group summaries")
}
skipped <- identical(summary$skipped[[1]], "True")
provided <- vapply(active_args, function(name) nzchar(args[[name]]), logical(1))


if (skipped) {
  if (any(provided)) die("skipped rare-variant report received active-only artifacts")
  lines <- c(
    paste0("# Rare-Variant REGENIE and ReMeta Export Report: ",
      config$project$analysis_name, " / ", args$trait), "",
    "## Scope and Export Status", "",
    "This trait did not enter the target-region REGENIE/ReMeta export branch.", "",
    "- Export status: skipped",
    paste0("- Trait-specific group: ", summary$group[[1]]),
    paste0("- Trait type: ", summary$trait_type[[1]]),
    paste0("- Genome build: ", args$build),
    paste0("- Candidate trait samples: ", summary$usable_n[[1]]),
    paste0("- Candidate cases: ", summary$cases[[1]]),
    paste0("- Candidate controls: ", summary$controls[[1]]),
    paste0("- Phase 2 skip reason: ", summary$skip_reason[[1]]), "",
    "No rare-variant HTP file, target PGEN, LD matrices, metrics, or plots were generated."
  )
  ensure_parent(args$out)
  writeLines(lines, args$out)
  quit(save = "no", status = 0L)
}


missing_args <- active_args[!provided]
if (length(missing_args)) die("active rare-variant report is missing argument(s): ", paste(missing_args, collapse = ", "))
for (path in unlist(args[active_args], use.names = FALSE)) {
  if (!file.exists(path) || file.info(path)$size == 0) die("rare-variant report input is empty: ", path)
}

target <- read_tsv(args[["target-summary"]])
validation <- read_tsv(args$validation)
metrics <- read_tsv(args$metrics)
hits <- read_tsv(args[["top-hits"]])
require_columns(target, c("key", "value"), "ReMeta target summary")
require_columns(validation, c("key", "value"), "ReMeta validation")
require_columns(metrics, c("metric", "value"), "rare-variant metrics")
require_values(target, c(
  "data_source", "genotype_mode", "sample_count", "variant_count", "min_mac",
  "geno_missing_max", "info_min", "target_regions"
), "key", "ReMeta target summary")
require_values(validation, c(
  "status", "group", "trait", "sample_count", "ordinary_regenie_sample_count",
  "rare_regenie_sample_count", "htp_sha256", "htp_variant_rows", "target_variant_count",
  "unique_ld_target_variant_count", "target_variants_not_indexed", "target_variant_ld_coverage_pct",
  "reference_gene_count", "indexed_gene_count", "genes_without_indexed_variants",
  "indexed_gene_coverage_pct", "ld_gene_variant_assignments", "ld_assignments_within_gene_bounds",
  "ld_assignments_outside_gene_bounds", "ld_assignment_gene_bound_coverage_pct"
), "key", "ReMeta validation")
required_metrics <- c(
  "total_variants", "valid_p_value_variants", "lambda_gc", "genomewide_significant_variants",
  "suggestive_variants", "plot_thinning_applied", "genotype_mode", "mac_definition",
  "mac_bin_1", "mac_bin_2_5", "mac_bin_6_10", "mac_bin_11_20", "mac_bin_gt20",
  "qq_eligible_variants", "qq_points_plotted", "manhattan_eligible_variants",
  "manhattan_points_plotted", "mac_qq_eligible_variants", "mac_qq_points_plotted",
  "effect_frequency_eligible_variants", "effect_frequency_points_plotted"
)
require_values(metrics, required_metrics, "metric", "rare-variant metrics")

metric <- function(key, default = "NA") value_for(metrics, key, "metric", default)
validated <- function(key, default = "NA") value_for(validation, key, "key", default)
target_value <- function(key, default = "NA") value_for(target, key, "key", default)
model_n <- as.integer(summary$model_sample_count[[1]])
sample_counts <- as.integer(c(
  target_value("sample_count"), validated("sample_count"),
  validated("ordinary_regenie_sample_count"), validated("rare_regenie_sample_count")
))
if (any(is.na(sample_counts)) || any(sample_counts != model_n)) {
  die("rare-variant report sample counts do not match the final singleton model")
}
if (validated("status") != "validated" || validated("group") != summary$group[[1]] ||
    validated("trait") != args$trait) {
  die("rare-variant validation identity does not match the requested trait")
}
if (validated("htp_sha256") != sha256_file(args$htp)) die("rare-variant HTP changed after validation")
if (as.integer(validated("htp_variant_rows")) != as.integer(metric("total_variants"))) {
  die("rare-variant metric row count does not match validated HTP rows")
}
if (as.integer(target_value("variant_count")) != as.integer(validated("target_variant_count"))) {
  die("target summary and ReMeta validation variant counts differ")
}
if (target_value("genotype_mode") != metric("genotype_mode")) {
  die("target summary and rare-variant plotting genotype modes differ")
}
mac_total <- sum(as.integer(vapply(
  c("mac_bin_1", "mac_bin_2_5", "mac_bin_6_10", "mac_bin_11_20", "mac_bin_gt20"),
  metric, character(1)
)))
if (mac_total != as.integer(metric("valid_p_value_variants")) ||
    as.integer(metric("mac_qq_eligible_variants")) != mac_total) {
  die("MAC strata do not reproduce all valid rare-variant P values")
}

trait_type <- summary$trait_type[[1]]
mac_definition <- metric("mac_definition")
mac_labels <- if (grepl("expected", mac_definition, fixed = TRUE)) {
  c("<=1", ">1–5", ">5–10", ">10–20", ">20")
} else {
  c("1", "2–5", "6–10", "11–20", ">20")
}
lambda_value <- suppressWarnings(as.numeric(metric("lambda_gc")))
lambda_label <- if (is.finite(lambda_value)) sprintf("%.6f", lambda_value) else "NA"
plot_links <- vapply(
  c("qq", "manhattan", "mac-qq", "effect-frequency"),
  function(name) report_relative_path(args[[name]], args$out),
  character(1)
)

lines <- c(
  paste0("# Rare-Variant REGENIE and ReMeta Export Report: ",
    config$project$analysis_name, " / ", args$trait), "",
  "## Scope and Export Status", "",
  paste0(
    "This report summarizes existing target-region single-variant REGENIE HTP statistics and their ",
    "sample-matched ReMeta LD export. It does not contain a local burden, SKAT, mask, gene, or ",
    "meta-analysis result."
  ), "",
  paste0(
    "For binary traits, HTP effects are ALT-allele odds ratios as defined by the ",
    "[ReMeta HTP specification](https://rgcgithub.github.io/remeta/file_formats/)."
  ), "",
  "- Export status: validated",
  paste0("- Trait-specific group: ", summary$group[[1]]),
  paste0("- Trait type: ", trait_type),
  paste0("- Genome build: ", args$build),
  paste0("- HTP file: `", args$htp, "`"), "",
  "## Sample and Model Validation", "",
  paste0("- Final model samples: ", model_n),
  paste0("- Final model cases: ", summary$model_cases[[1]]),
  paste0("- Final model controls: ", summary$model_controls[[1]]),
  "- Ordinary REGENIE IDs, rare-variant REGENIE IDs, model keep, and target PSAM: matched",
  paste0("- Covariates: ", summary$covariates[[1]]), "",
  "## Target Regions and Variant QC", "",
  paste0("- Data source: ", target_value("data_source")),
  paste0("- Genotype mode: ", target_value("genotype_mode")),
  paste0("- Target-region variants after QC: ", target_value("variant_count")),
  paste0("- Minimum MAC: ", target_value("min_mac")),
  paste0("- Maximum genotype missingness: ", target_value("geno_missing_max")),
  paste0("- Minimum INFO/MaCH R2: ", target_value("info_min")),
  paste0("- Target-region resource: `", target_value("target_regions"), "`"),
  "- Cohort AAF/MAC is descriptive; central ReMeta rarity and mask rules remain authoritative.", "",
  "## Single-Variant Association Summary", "",
  paste0("- HTP variant rows: ", metric("total_variants")),
  paste0("- Valid P-value variants: ", metric("valid_p_value_variants")),
  paste0("- Lambda GC: ", lambda_label),
  paste0("- Genome-wide significant variants (P <= 5e-8): ", metric("genomewide_significant_variants")),
  paste0("- Suggestive variants (P <= 1e-5): ", metric("suggestive_variants")),
  paste0("- MAC definition: ", mac_definition), "",
  "Overall lambda GC is descriptive for this target-region set; discrete low-MAC tests are better assessed in the MAC-stratified QQ plot.", "",
  "| MAC stratum | Valid variants |", "| --- | ---: |",
  paste0("| ", mac_labels[[1]], " | ", metric("mac_bin_1"), " |"),
  paste0("| ", mac_labels[[2]], " | ", metric("mac_bin_2_5"), " |"),
  paste0("| ", mac_labels[[3]], " | ", metric("mac_bin_6_10"), " |"),
  paste0("| ", mac_labels[[4]], " | ", metric("mac_bin_11_20"), " |"),
  paste0("| ", mac_labels[[5]], " | ", metric("mac_bin_gt20"), " |"), "",
  "## Top Single-Variant Signals", "",
  top_hit_lines(hits, trait_type, mac_definition), "",
  "## Diagnostic Plots", "",
  paste0("![Single-variant QQ plot](", plot_links[["qq"]], ")"), "",
  paste0("![Target-region Manhattan plot](", plot_links[["manhattan"]], ")"), "",
  paste0("![MAC-stratified QQ plot](", plot_links[["mac-qq"]], ")"), "",
  paste0("![Effect versus frequency plot](", plot_links[["effect-frequency"]], ")"), "",
  paste0("- Manhattan PDF: `", args[["manhattan-pdf"]], "`"),
  paste0("- Plot thinning applied: ", metric("plot_thinning_applied")),
  "- All metrics and top-hit rankings use every valid HTP row.", "",
  "| Plot | Eligible variants | Plotted points |", "| --- | ---: | ---: |",
  paste0("| QQ | ", metric("qq_eligible_variants"), " | ", metric("qq_points_plotted"), " |"),
  paste0("| Manhattan | ", metric("manhattan_eligible_variants"), " | ",
    metric("manhattan_points_plotted"), " |"),
  paste0("| MAC-stratified QQ | ", metric("mac_qq_eligible_variants"), " | ",
    metric("mac_qq_points_plotted"), " |"),
  paste0("| Effect versus frequency | ", metric("effect_frequency_eligible_variants"), " | ",
    metric("effect_frequency_points_plotted"), " |"), "",
  "## ReMeta LD Target Coverage", "",
  paste0("- Validated model samples: ", validated("sample_count")),
  paste0("- Ordinary REGENIE samples: ", validated("ordinary_regenie_sample_count")),
  paste0("- Rare-variant REGENIE samples: ", validated("rare_regenie_sample_count")),
  "- LD components: all 22 autosomes passed gene-LD, buffer-LD, and index validation",
  paste0("- LD prefix: `results/remeta/export/", args$build, "/ld/", summary$group[[1]], "/chr{1..22}.remeta.*`"), "",
  "| Coverage metric | Observed | Denominator | Coverage |",
  "| --- | ---: | ---: | ---: |",
  paste0("| Target genes with at least one indexed LD variant | ", validated("indexed_gene_count"),
    " | ", validated("reference_gene_count"), " | ", format_percent(validated("indexed_gene_coverage_pct")), " |"),
  paste0("| QC-passing target-region variants represented in LD indexes | ",
    validated("unique_ld_target_variant_count"), " | ", validated("target_variant_count"), " | ",
    format_percent(validated("target_variant_ld_coverage_pct")), " |"),
  paste0("| Indexed gene-variant assignments within declared gene spans | ",
    validated("ld_assignments_within_gene_bounds"), " | ", validated("ld_gene_variant_assignments"), " | ",
    format_percent(validated("ld_assignment_gene_bound_coverage_pct")), " |"), "",
  paste0("- Target genes without an indexed LD variant: ", validated("genes_without_indexed_variants")),
  paste0("- Target-region variants absent from every LD gene index: ", validated("target_variants_not_indexed")),
  paste0("- Indexed assignments outside declared gene spans: ", validated("ld_assignments_outside_gene_bounds")),
  "- Conditional buffer variants: not included (`--skip-buffer`; marginal LD export).", "",
  "## Central Handoff Boundary", "",
  paste0(
    "Gene masks, external rarity harmonization, ReMeta gene tests, and cross-cohort meta-analysis ",
    "must be performed centrally from the validated HTP and matching LD prefix; see the ",
    "[ReMeta gene workflow](https://rgcgithub.github.io/remeta/documentation/)."
  )
)

ensure_parent(args$out)
writeLines(lines, args$out)
cat("Wrote rare-variant report:", args$out, "\n")
