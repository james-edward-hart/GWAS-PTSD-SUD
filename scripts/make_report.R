#!/usr/bin/env Rscript

# Write one compact markdown report per trait and ancestry stratum.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Parse all report inputs for one trait/ancestry pair.
args <- parse_args(defaults = list("reference-prep-report" = "NA", "manhattan-pdf" = ""))
require_args(args, c(
  "config", "trait", "ancestry", "build", "stats", "gwas-summary", "qq", "manhattan",
  "popmad-plot",
  "strata-counts", "pheno", "covar", "keep", "relatedness-summary",
  "sex-check-summary", "genome-build-details", "ancestry-counts",
  "admixture-summary", "admixture-study", "admixture-comparison", "admixture-report",
  "software", "reference", "plink-log", "out"
))


# Load report inputs.
config <- load_config(args$config)
stats <- read_tsv(args$stats)
keep <- read_tsv(args$keep)
pheno <- read_tsv(args$pheno)
covar <- read_tsv(args$covar)


fmt <- function(value) {
  if (is.null(value) || !length(value)) return("NA")
  value <- as.character(value[[1]])
  if (!nzchar(value) || value %in% c("NA", "NaN")) "NA" else value
}

fmt_count <- function(value) {
  value <- fmt(value)
  number <- suppressWarnings(as.numeric(value))
  if (is.finite(number)) format(round(number), big.mark = ",", scientific = FALSE, trim = TRUE) else value
}

fmt_decimal <- function(value, digits = 3) {
  value <- fmt(value)
  number <- suppressWarnings(as.numeric(value))
  if (is.finite(number)) sprintf(paste0("%.", digits, "f"), number) else value
}

fmt_p <- function(value) {
  value <- fmt(value)
  number <- suppressWarnings(as.numeric(value))
  if (is.finite(number)) format(number, scientific = TRUE, digits = 3) else value
}

fmt_percent <- function(value) {
  value <- fmt(value)
  number <- suppressWarnings(as.numeric(value))
  if (is.finite(number)) paste0(sprintf("%.1f", 100 * number), "%") else value
}

markdown_escape <- function(value) {
  gsub("\\|", "\\\\|", fmt(value))
}

report_relative_path <- function(path, out_path) {
  path <- gsub("\\\\", "/", path)
  out_path <- gsub("\\\\", "/", out_path)
  if (startsWith(path, "results/") && startsWith(out_path, "results/reports/")) {
    return(file.path("..", "..", sub("^results/", "", path)))
  }
  path
}

genomic_lambda <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) return(NA_real_)
  chi <- suppressWarnings(qchisq(p, df = 1, lower.tail = FALSE))
  chi <- chi[is.finite(chi)]
  if (!length(chi)) return(NA_real_)
  median(chi, na.rm = TRUE) / qchisq(0.5, df = 1, lower.tail = FALSE)
}

top_signal_lines <- function(stats, n = 10) {
  if (!"p" %in% names(stats)) return("- No P-value column was available.")
  stats$p_num <- suppressWarnings(as.numeric(stats$p))
  rows <- stats[is.finite(stats$p_num) & stats$p_num > 0 & stats$p_num <= 1, , drop = FALSE]
  if (!nrow(rows)) return("- No valid P values were available.")
  rows <- head(rows[order(rows$p_num), , drop = FALSE], n)

  get_col <- function(name) if (name %in% names(rows)) rows[[name]] else rep("NA", nrow(rows))
  location <- paste0(get_col("chrom"), ":", get_col("pos"))
  table_rows <- c(
    "Variant | Location | Effect allele | A1 freq | Beta/log(OR) | SE | P",
    "--- | --- | --- | ---: | ---: | ---: | ---:"
  )
  for (i in seq_len(nrow(rows))) {
    table_rows <- c(table_rows, paste(
      markdown_escape(get_col("variant_id")[[i]]),
      markdown_escape(location[[i]]),
      markdown_escape(get_col("effect_allele")[[i]]),
      fmt_decimal(get_col("a1_freq")[[i]], 4),
      fmt_decimal(get_col("beta_or_log_or")[[i]], 4),
      fmt_decimal(get_col("se")[[i]], 4),
      fmt_p(rows$p_num[[i]]),
      sep = " | "
    ))
  }
  table_rows <- sub("^", "| ", table_rows)
  sub("$", " |", table_rows)
}


# Count analyzed cases and controls from the final keep file.
keep_key <- paste(keep$FID, keep$IID, sep = "\t")
pheno <- pheno[paste(pheno$FID, pheno$IID, sep = "\t") %in% keep_key, , drop = FALSE]
cases <- sum(pheno$PHENO == "2")
controls <- sum(pheno$PHENO == "1")
n_samples <- cases + controls
underpowered <- n_samples < config$warnings$min_n ||
  cases < config$warnings$min_cases ||
  controls < config$warnings$min_controls


# Summarize harmonized GWAS result content.
p <- suppressWarnings(as.numeric(stats$p))
valid_p <- is.finite(p) & p > 0 & p <= 1
n_variants <- nrow(stats)
min_p <- if (any(valid_p)) min(p[valid_p], na.rm = TRUE) else "NA"
lambda_gc <- genomic_lambda(p)
tests <- if ("test" %in% names(stats)) paste(sort(unique(stats$test[nzchar(stats$test)])), collapse = ", ") else "NA"
if (!nzchar(tests)) tests <- "NA"
covariates <- paste(names(covar)[-(1:2)], collapse = ", ")
top_signals <- top_signal_lines(stats)


# Read metric/value summaries into named vectors.
kv <- function(path) {
  rows <- read_tsv(path)
  setNames(rows$value, rows$metric)
}
metric_value <- function(values, name) {
  if (!is.null(names(values)) && name %in% names(values)) values[[name]] else "NA"
}
relatedness <- kv(args[["relatedness-summary"]])
sex_check <- kv(args[["sex-check-summary"]])
gwas_summary <- kv(args[["gwas-summary"]])
admixture_summary <- kv(args[["admixture-summary"]])


# Pull trait/ancestry counts before sex-check and relatedness intersections.
strata_counts <- read_tsv(args[["strata-counts"]])
stratum <- strata_counts[strata_counts$trait_id == args$trait & strata_counts$ancestry == args$ancestry, , drop = FALSE]
stratum <- stratum[1, , drop = FALSE]
trait_strata <- strata_counts[strata_counts$trait_id == args$trait, , drop = FALSE]
if (nrow(trait_strata)) {
  trait_strata <- trait_strata[order(trait_strata$ancestry), , drop = FALSE]
  strata_lines <- c(
    "| Ancestry | Phenotype-complete N | Cases | Controls | Active | Exclusion reason | Underpowered |",
    "| --- | ---: | ---: | ---: | --- | --- | --- |",
    vapply(seq_len(nrow(trait_strata)), function(i) {
      paste0(
        "| ", markdown_escape(trait_strata$ancestry[[i]]),
        " | ", fmt_count(trait_strata$n[[i]]),
        " | ", fmt_count(trait_strata$cases[[i]]),
        " | ", fmt_count(trait_strata$controls[[i]]),
        " | ", markdown_escape(trait_strata$active[[i]]),
        " | ", markdown_escape(trait_strata$excluded_reason[[i]]),
        " | ", markdown_escape(trait_strata$underpowered[[i]]),
        " |"
      )
    }, character(1))
  )
} else {
  strata_lines <- "- No trait-level strata rows were available."
}


# Select the chosen genome-build detail row.
build_details <- read_tsv(args[["genome-build-details"]])
selected <- build_details[build_details$selected == "True", , drop = FALSE]
if (!nrow(selected)) selected <- build_details[build_details$build == args$build, , drop = FALSE]
selected <- selected[1, , drop = FALSE]


# Pull the total number of POP-MaD exclusions.
ancestry_counts <- read_tsv(args[["ancestry-counts"]])
excluded_total <- ancestry_counts$n[ancestry_counts$category == "excluded_total"]
if (!length(excluded_total)) excluded_total <- "NA"

assigned_counts <- ancestry_counts[ancestry_counts$category == "assigned_ancestry", , drop = FALSE]
if (nrow(assigned_counts)) {
  assigned_counts <- assigned_counts[order(assigned_counts$ancestry), , drop = FALSE]
  assigned_total <- sum(suppressWarnings(as.numeric(assigned_counts$n)), na.rm = TRUE)
  popmad_assigned_lines <- c(
    "| Ancestry | POP-MaD assigned samples |",
    "| --- | ---: |",
    vapply(seq_len(nrow(assigned_counts)), function(i) {
      paste0("| ", markdown_escape(assigned_counts$ancestry[[i]]), " | ",
        fmt_count(assigned_counts$n[[i]]), " |")
    }, character(1))
  )
} else {
  assigned_total <- "NA"
  popmad_assigned_lines <- "- No POP-MaD assigned-count rows were available."
}

admixture_labels <- as.character(unlist(config$admixture$labels %||% character(), use.names = FALSE))
mean_metrics <- grep("^mean_study_proportion_", names(admixture_summary), value = TRUE)
mean_metric_labels <- sub("^mean_study_proportion_", "", mean_metrics)
if (!length(admixture_labels)) admixture_labels <- mean_metric_labels
admixture_labels <- unique(c(admixture_labels[admixture_labels %in% mean_metric_labels],
  mean_metric_labels[!mean_metric_labels %in% admixture_labels]))
if (length(admixture_labels)) {
  admixture_mean_lines <- c(
    "| ADMIXTURE ancestry | Mean study proportion | Study top-component samples |",
    "| --- | ---: | ---: |",
    vapply(admixture_labels, function(label) {
      paste0("| ", markdown_escape(label), " | ",
        fmt_decimal(metric_value(admixture_summary, paste0("mean_study_proportion_", label)), 3), " | ",
        fmt_count(metric_value(admixture_summary, paste0("n_study_top_component_", label))), " |")
    }, character(1))
  )
} else {
  admixture_mean_lines <- "- No ADMIXTURE mean-proportion rows were available."
}


# Extract useful one-line highlights from the PLINK log.
log_lines <- character()
if (file.exists(args[["plink-log"]])) {
  lines <- readLines(args[["plink-log"]], warn = FALSE)
  hits <- grepl("samples|variants|removed|remaining|covariate|phenotype", lines, ignore.case = TRUE)
  log_lines <- trimws(lines[hits & nzchar(trimws(lines))])
  log_lines <- head(log_lines, 12)
}
if (!length(log_lines)) log_lines <- "No highlights captured."


input_manifest_release_lines <- function(path) {
  if (blank(path) || !file.exists(path)) return(character())
  rows <- read_tsv(path)
  if (!all(c("file_role", "cohort_data_release") %in% names(rows))) return(character())
  rows <- rows[nzchar(rows$file_role), , drop = FALSE]
  if (!nrow(rows)) return(character())
  paste0("- Input release note for ", rows$file_role, ": ", rows$cohort_data_release)
}


reference_panel_lines <- function(config, section, label) {
  block <- config[[section]] %||% list()
  if (!truthy(block$enabled %||% FALSE)) return(character())
  c(
    paste0("- ", label, " panel: ", block$source_panel_id %||% "NA"),
    paste0("- ", label, " reference build: ", block$reference_genome_build %||% "NA"),
    paste0("- ", label, " reference metadata: `", block$metadata$path %||% "NA", "`")
  )
}


study_components <- genotype_component_paths(config$genotypes, "study")
qq_link <- report_relative_path(args$qq, args$out)
manhattan_link <- report_relative_path(args$manhattan, args$out)
popmad_plot_link <- report_relative_path(args[["popmad-plot"]], args$out)
plot_lines <- c(
  paste0("![POP-MaD projected PC space](", popmad_plot_link, ")"),
  "",
  paste0("![QQ plot](", qq_link, ")"),
  "",
  paste0("![Manhattan plot](", manhattan_link, ")"),
  "",
  paste0("- POP-MaD projection plot: `", args[["popmad-plot"]], "`"),
  paste0("- QQ plot file: `", args$qq, "`"),
  paste0("- Manhattan plot file (PNG): `", args$manhattan, "`"),
  if (nzchar(args[["manhattan-pdf"]])) paste0("- Manhattan plot (PDF): `", args[["manhattan-pdf"]], "`")
)


# Keep reports compact while retaining enough QC evidence for review.
text <- c(
  paste0("# GWAS Report: ", config$project$analysis_name, " / ", args$trait, " / ", args$ancestry),
  "",
  "## Run",
  "",
  paste0("- Analysis: ", config$project$analysis_name),
  paste0("- Trait: ", args$trait),
  paste0("- Ancestry: ", args$ancestry),
  paste0("- Genome build label: ", args$build),
  "- Engine: PLINK2 `--glm`",
  paste0("- Effective config: `", args$config, "`"),
  paste0("- PLINK log: `", args[["plink-log"]], "`"),
  "",
  "## Input Files",
  "",
  paste0("- Cohort data release: ", config$project$cohort_data_release %||% "unspecified"),
  paste0("- Sample manifest: `", config$inputs$sample_manifest, "`"),
  paste0("- Trait registry: `", config$inputs$trait_registry, "`"),
  paste0("- Study genotype type: ", config$genotypes$type),
  paste0("- Study genotype prefix: `", config$genotypes$prefix, "`"),
  paste0("- Study genotype component (", names(study_components), "): `", unname(study_components), "`"),
  paste0("- Inferred genome build: ", args$build),
  paste0("- Genome-build marker details: `", args[["genome-build-details"]], "`"),
  paste0("- Reference package root: `", config$reference_package$root %||% "NA", "`"),
  paste0("- Reference package fingerprint: ", config$reference_package$observed_fingerprint %||% config$reference_package$fingerprint %||% "NA"),
  reference_panel_lines(config, "ancestry_reference", "POP-MaD"),
  reference_panel_lines(config, "admixture", "ADMIXTURE"),
  input_manifest_release_lines(config$resources$input_manifest %||% ""),
  "",
  "## Ancestry and ADMIXTURE QC",
  "",
  paste0("- POP-MaD assigned samples: ", fmt_count(assigned_total)),
  paste0("- POP-MaD excluded samples: ", fmt_count(excluded_total[[1]])),
  paste0("- ADMIXTURE study samples: ", fmt_count(metric_value(admixture_summary, "n_study_samples"))),
  paste0("- ADMIXTURE reference samples: ", fmt_count(metric_value(admixture_summary, "n_reference_samples"))),
  paste0("- ADMIXTURE merged LD-pruned variants: ", fmt_count(metric_value(admixture_summary, "n_merged_variants"))),
  paste0("- Mean study top ADMIXTURE proportion: ", fmt_decimal(metric_value(admixture_summary, "mean_study_top_proportion"), 3)),
  paste0("- POP-MaD/ADMIXTURE comparable samples: ", fmt_count(metric_value(admixture_summary, "popmad_comparable_samples"))),
  paste0("- POP-MaD/ADMIXTURE matches: ", fmt_count(metric_value(admixture_summary, "popmad_matches")),
    "; discordant: ", fmt_count(metric_value(admixture_summary, "popmad_discordant")),
    "; match rate: ", fmt_percent(metric_value(admixture_summary, "popmad_match_rate"))),
  "",
  "### POP-MaD Assigned Counts",
  "",
  popmad_assigned_lines,
  "",
  "### ADMIXTURE Mean Study Proportions",
  "",
  admixture_mean_lines,
  "",
  paste0("- ADMIXTURE run summary: `", args[["admixture-summary"]], "`"),
  paste0("- ADMIXTURE study proportions: `", args[["admixture-study"]], "`"),
  paste0("- POP-MaD/ADMIXTURE comparison: `", args[["admixture-comparison"]], "`"),
  paste0("- ADMIXTURE QC report: `", args[["admixture-report"]], "`"),
  "",
  "### Trait GWAS Strata",
  "",
  strata_lines,
  "",
  "## Sample Filtering",
  "",
  paste0("- Source genotype samples: ", fmt_count(metric_value(gwas_summary, "source_genotype_samples"))),
  paste0("- Sample manifest rows: ", fmt_count(metric_value(gwas_summary, "sample_manifest_rows"))),
  paste0("- Ancestry stratum phenotype-complete samples before sex/relatedness filters: ", fmt_count(stratum$n),
    " (cases=", fmt_count(stratum$cases), ", controls=", fmt_count(stratum$controls), ")"),
  paste0("- Final GWAS keep samples after ancestry, sex-check, and relatedness filters: ", fmt_count(nrow(keep))),
  paste0("- Final analyzed phenotype count: ", fmt_count(n_samples), " (cases=", fmt_count(cases), ", controls=", fmt_count(controls), ")"),
  paste0("- PLINK loaded samples: ", fmt_count(metric_value(gwas_summary, "plink_loaded_samples"))),
  paste0("- PLINK samples after `--keep`: ", fmt_count(metric_value(gwas_summary, "plink_keep_remaining_samples"))),
  paste0("- PLINK samples after model/QC filters: ", fmt_count(metric_value(gwas_summary, "plink_final_samples"))),
  paste0("- Underpowered warning: ", ifelse(underpowered, "True", "False")),
  paste0("- POP-MaD excluded samples: ", excluded_total[[1]]),
  paste0("- Sex-check removed samples: ", metric_value(sex_check, "n_removed")),
  "",
  "## Covariates",
  "",
  paste0("- Covariate file: `", args$covar, "`"),
  paste0("- Covariates used: ", covariates),
  paste0("- PLINK covariate variance standardization: ", ifelse(truthy(config$gwas$covar_variance_standardize), "True", "False")),
  "",
  "## Variant Filtering and Results",
  "",
  paste0("- Harmonized summary statistics: `", args$stats, "`"),
  paste0("- GWAS filter summary: `", args[["gwas-summary"]], "`"),
  paste0("- Source genotype variants: ", fmt_count(metric_value(gwas_summary, "source_genotype_variants"))),
  paste0("- PLINK loaded variants: ", fmt_count(metric_value(gwas_summary, "plink_loaded_variants"))),
  paste0("- Initial PLINK variant exclusions: ", fmt_count(metric_value(gwas_summary, "plink_initial_filter_excluded_variants")),
    " excluded; ", fmt_count(metric_value(gwas_summary, "plink_initial_filter_remaining_variants")), " remaining"),
  paste0("- Controls used for HWE filtering: ", fmt_count(metric_value(gwas_summary, "control_hwe_controls"))),
  paste0("- Variants passing controls-only HWE prefilter: ", fmt_count(metric_value(gwas_summary, "control_hwe_passing_variants"))),
  paste0("- Variants removed by `--geno`: ", fmt_count(metric_value(gwas_summary, "plink_geno_removed_variants"))),
  paste0("- Variants removed by `--maf`: ", fmt_count(metric_value(gwas_summary, "plink_maf_removed_variants"))),
  paste0("- Variants removed by controls-only `--hwe`: ", fmt_count(metric_value(gwas_summary, "plink_hwe_removed_variants"))),
  paste0("- Variants removed by `--mach-r2-filter`: ", fmt_count(metric_value(gwas_summary, "plink_info_removed_variants"))),
  paste0("- PLINK variants after main filters: ", fmt_count(metric_value(gwas_summary, "plink_final_variants"))),
  paste0("- Variants in harmonized analysis output: ", fmt_count(n_variants)),
  paste0("- Variants with valid P values: ", fmt_count(sum(valid_p))),
  paste0("- Genome-wide significant variants (P <= 5e-8): ", fmt_count(metric_value(gwas_summary, "genomewide_significant_variants"))),
  paste0("- Suggestive variants (P <= 1e-5): ", fmt_count(metric_value(gwas_summary, "suggestive_variants"))),
  paste0("- Minimum P value: ", fmt_p(min_p)),
  paste0("- Genomic inflation factor (lambda GC): ", ifelse(is.finite(lambda_gc), sprintf("%.3f", lambda_gc), "NA")),
  paste0("- Retained PLINK2 test terms: ", tests),
  "",
  "## Top Association Signals",
  "",
  top_signals,
  "",
  "## Plots",
  "",
  plot_lines,
  "",
  "## QC Settings",
  "",
  paste0("- INFO/R2 filter enabled: ", ifelse(truthy(config$qc$use_mach_r2_filter %||% TRUE), "True", "False")),
  paste0("- INFO/R2 minimum when enabled: ", config$qc$info_min),
  paste0("- MAF minimum: ", config$qc$maf_min),
  paste0("- HWE P minimum: ", config$qc$hwe_p_min, " (calculated in controls only)"),
  paste0("- Genotype missingness maximum: ", config$qc$geno_missing_max),
  paste0("- Sample missingness maximum: ", config$qc$sample_missing_max),
  paste0("- Relatedness mode: ", config$relatedness$mode),
  paste0("- KING cutoff: ", config$relatedness$king_cutoff),
  paste0("- Relatedness marker set: `", args[["relatedness-summary"]], "`"),
  paste0("- Relatedness LD-pruned variants: ", metric_value(relatedness, "prune_in_variants")),
  paste0("- Sex-check action: ", metric_value(sex_check, "action")),
  paste0("- Sex-check status: ", metric_value(sex_check, "status")),
  paste0("- Sex-check problems: ", metric_value(sex_check, "n_problems")),
  paste0("- Sex-check removed samples: ", metric_value(sex_check, "n_removed")),
  paste0("- Genome-build checked markers: ", selected$checked_markers %||% "NA"),
  paste0("- Genome-build matching markers: ", selected$matching_markers %||% "NA"),
  paste0("- Genome-build match fraction: ", selected$match_fraction %||% "NA"),
  paste0("- Reference prep report: `", args[["reference-prep-report"]], "`"),
  "",
  "## PLINK Log Highlights",
  "",
  paste0("- ", log_lines),
  "",
  "## Manifests",
  "",
  paste0("- Software manifest: `", args$software, "`"),
  paste0("- Reference data manifest: `", args$reference, "`")
)


# Write the markdown report.
ensure_parent(args$out)
writeLines(text, args$out)
cat("Wrote report:", args$out, "\n")
