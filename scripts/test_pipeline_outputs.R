#!/usr/bin/env Rscript

# Check that a production Stage 1 run produced expected output bundles.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Parse expected results location and optional trait/ancestry filters.
args <- parse_args(defaults = list(results = "results", config = "", trait = "", build = "auto", ancestries = ""))
results <- args$results
config_path <- args$config
if (blank(config_path)) config_path <- file.path(results, "config", "resolved_config.yaml")
require_existing_file(config_path, "resolved run config")
config <- load_config(config_path)
analysis_name <- analysis_output_name(config)

trait_ids <- if (blank(args$trait)) {
  traits <- read_tsv(config$inputs$trait_registry)
  require_columns(traits, "trait_id", "trait registry")
  traits$trait_id[nzchar(traits$trait_id)]
} else {
  split_csv(args$trait)
}
ancestry_labels <- if (blank(args$ancestries)) {
  as.character(unlist(config$analysis$ancestries, use.names = FALSE))
} else {
  split_csv(args$ancestries)
}
if (!length(trait_ids)) die("no traits configured for output validation")
if (!length(ancestry_labels)) die("no ancestries configured for output validation")


# Small assertions for required files and non-empty tables.
require_file <- function(path, message) {
  if (!file.exists(path)) die(message, ": ", path)
}
require_nonempty_file <- function(path, message) {
  require_file(path, message)
  if (file.info(path)$size <= 0) die(message, " is empty: ", path)
}
count_rows <- function(path) max(length(readLines(path, warn = FALSE)) - 1, 0)
assert_section_order <- function(text, headings, report) {
  positions <- match(headings, text)
  if (any(is.na(positions))) {
    die("report missing section(s): ", paste(headings[is.na(positions)], collapse = ", "), ": ", report)
  }
  if (any(diff(positions) <= 0)) die("report sections are out of order: ", report)
}
id_keys <- function(path, label) {
  ids <- read_id_file(path, label)
  require_unique_ids(ids, label)
  sort(paste(ids$FID, ids$IID, sep = "\t"))
}


# Check workflow-wide QC and provenance outputs.
build_file <- file.path(results, "qc", "genome_build", "genome_build.txt")
require_file(build_file, "missing inferred genome build")
require_file(file.path(results, "config", "effective_config.yaml"), "missing effective config snapshot")
require_file(file.path(results, "config", "resolved_config.yaml"), "missing resolved config snapshot")
require_file(file.path(results, "qc", "relatedness", "relatedness_summary.tsv"), "missing relatedness marker summary")
sex_summary_path <- file.path(results, "qc", "sex", "sex_check_summary.tsv")
require_file(sex_summary_path, "missing sex-check summary")
sex_prune_prefix <- file.path(results, "qc", "sex", "plink_sex_check.sex_marker_prune.prune")
require_file(paste0(sex_prune_prefix, ".in"), "missing sex-check retained-marker audit")
require_file(paste0(sex_prune_prefix, ".out"), "missing sex-check pruned-marker audit")
sex_summary <- read_tsv(sex_summary_path)
required_sex_metrics <- c(
  "source_x_variants", "source_y_variants", "par_handling",
  "post_qc_variants", "pruned_variants"
)
missing_sex_metrics <- setdiff(required_sex_metrics, sex_summary$metric)
if (length(missing_sex_metrics)) {
  die("sex-check summary is missing marker audit metrics: ", paste(missing_sex_metrics, collapse = ", "))
}
build <- if (args$build == "auto") trimws(readLines(build_file, warn = FALSE)[[1]]) else args$build


# Check POP-MaD output tables.
popmad_dir <- file.path(results, "qc", "ancestry", "production")
assignments <- file.path(popmad_dir, "popmad_assignments.tsv")
excluded <- file.path(popmad_dir, "popmad_excluded.tsv")
counts <- file.path(popmad_dir, "popmad_population_counts.tsv")
model_summary <- file.path(popmad_dir, "population_model_summary.tsv")
require_file(assignments, "missing POP-MaD assignments")
if (count_rows(assignments) <= 0) die("POP-MaD assignments are empty")
require_file(excluded, "missing POP-MaD excluded file")
require_file(counts, "missing POP-MaD population counts")
require_file(model_summary, "missing POP-MaD model summary")


# Check ADMIXTURE QC outputs.
admixture_dir <- file.path(results, "qc", "admixture")
admixture_summary <- file.path(admixture_dir, "admixture_run_summary.tsv")
require_file(file.path(admixture_dir, "study_ancestry_proportions.tsv"), "missing ADMIXTURE study proportions")
require_file(file.path(admixture_dir, "reference_ancestry_proportions.tsv"), "missing ADMIXTURE reference proportions")
require_file(file.path(admixture_dir, "popmad_admixture_comparison.tsv"), "missing POP-MaD/ADMIXTURE comparison")
require_file(admixture_summary, "missing ADMIXTURE run summary")
require_file(file.path(admixture_dir, "admixture_report.md"), "missing ADMIXTURE QC report")
popmad_plot <- file.path(results, "plots", "ancestry", paste0(analysis_name, ".popmad_reference_study_pcs.png"))
require_file(popmad_plot, "missing POP-MaD projection plot")
if (file.info(popmad_plot)$size <= 0) die("POP-MaD projection plot is empty: ", popmad_plot)
admixture_summary_rows <- read_tsv(admixture_summary)
if (!any(grepl("^mean_study_proportion_", admixture_summary_rows$metric))) {
  die("ADMIXTURE summary missing per-ancestry mean study proportions")
}
if (!any(grepl("^n_study_top_component_", admixture_summary_rows$metric))) {
  die("ADMIXTURE summary missing per-ancestry top-component counts")
}


# Check each expected trait/ancestry report bundle.
for (trait in trait_ids) {
  for (ancestry in ancestry_labels) {
    stats <- file.path(results, "gwas", trait, ancestry, paste0(trait, ".", ancestry, ".", build, ".plink2.glm.tsv"))
    summary <- file.path(results, "gwas", trait, ancestry, paste0(trait, ".", ancestry, ".", build, ".gwas_filter_summary.tsv"))
    report <- file.path(results, "reports", trait, paste0(analysis_name, ".", trait, ".", ancestry, ".", build, ".report.md"))
    qq <- file.path(results, "plots", trait, ancestry, paste0(analysis_name, ".", trait, ".", ancestry, ".", build, ".qq.png"))
    manhattan <- file.path(results, "plots", trait, ancestry, paste0(analysis_name, ".", trait, ".", ancestry, ".", build, ".manhattan.png"))
    manhattan_pdf <- file.path(results, "plots", trait, ancestry, paste0(analysis_name, ".", trait, ".", ancestry, ".", build, ".manhattan.pdf"))

    require_file(stats, "missing GWAS stats")
    require_file(summary, "missing GWAS filter summary")
    if (count_rows(stats) <= 0) die("GWAS stats are empty: ", stats)
    header <- names(read_tsv(stats))
    for (column in c("a1_freq", "mac", "info", "test")) {
      if (!column %in% header) die("GWAS stats missing harmonized column ", column, ": ", stats)
    }
    require_file(report, "missing report")
    report_text <- readLines(report, warn = FALSE)
    if (!any(grepl("Covariates used", report_text))) die("report missing covariate section: ", report)
    if (!any(grepl("Source genotype variants", report_text))) die("report missing variant-flow counts: ", report)
    if (!any(grepl("Genomic inflation factor", report_text))) die("report missing lambda GC: ", report)
    if (!any(grepl("Ancestry Assignment and ADMIXTURE QC", report_text))) die("report missing ancestry/ADMIXTURE summary: ", report)
    if (!any(grepl("ADMIXTURE Mean Study Proportions", report_text))) die("report missing ADMIXTURE mean proportions: ", report)
    if (!any(grepl("Top Association Signals", report_text))) die("report missing top-signal section: ", report)
    if (!any(grepl("!\\[POP-MaD projected PC space\\]", report_text))) die("report missing embedded POP-MaD plot: ", report)
    if (!any(grepl("!\\[QQ plot\\]", report_text))) die("report missing embedded QQ plot: ", report)
    if (!any(grepl("!\\[Manhattan plot\\]", report_text))) die("report missing embedded Manhattan plot: ", report)
    if (!any(grepl("Relatedness LD-pruned variants", report_text))) die("report missing relatedness details: ", report)
    if (!any(grepl("Sex-check problems", report_text))) die("report missing sex-check details: ", report)
    assert_section_order(report_text, c(
      "## Run Summary",
      "## Inputs and Reference Provenance",
      "## QC Settings",
      "## Ancestry Assignment and ADMIXTURE QC",
      "## Trait Strata and Sample Filtering",
      "## Phenotype and Covariates",
      "## Variant Filtering",
      "## Association Results",
      "## Top Association Signals",
      "## Plots",
      "## PLINK Log Highlights"
    ), report)
    if (!file.exists(qq) || file.info(qq)$size <= 0) die("missing QQ plot: ", qq)
    if (!file.exists(manhattan) || file.info(manhattan)$size <= 0) die("missing Manhattan plot: ", manhattan)
    if (!file.exists(manhattan_pdf) || file.info(manhattan_pdf)$size <= 0) die("missing Manhattan PDF plot: ", manhattan_pdf)
  }
}


# Check singleton Phase 2 and exact-sample ReMeta outputs when enabled.
if (truthy(config$phase2_regenie$enabled %||% FALSE)) {
  phase2_dir <- file.path(results, "qc", "phase2_regenie")
  status_path <- file.path(phase2_dir, "trait_group_status.tsv")
  require_file(status_path, "missing Phase 2 trait-group status")
  status <- read_tsv(status_path)
  require_columns(status, c(
    "group", "trait", "trait_type", "usable_n", "cases", "controls",
    "model_sample_count", "keep_count", "model_keep_sha256", "skipped", "skip_reason"
  ), status_path)
  if (nrow(status) != length(trait_ids) || anyDuplicated(status$group) || anyDuplicated(status$trait) ||
      !setequal(status$trait, trait_ids)) {
    die("Phase 2 status is not a one-to-one mapping of configured traits")
  }

  expected_active <- data.frame(
    trait = c("co_ptsd_anysud", "co_ptsd_aud", "co_ptsd_cud", "co_ptsd_oud"),
    usable_n = c(625L, 885L, 768L, 1112L),
    cases = c(258L, 225L, 191L, 50L),
    controls = c(367L, 660L, 577L, 1062L),
    stringsAsFactors = FALSE
  )
  # These independently established counts are specific to the current ADAA production cohort.
  if (identical(analysis_name, "ADAA_Comorbid_GWAS") &&
      all(c(expected_active$trait, "co_ptsd_tud") %in% trait_ids)) {
    observed <- status[match(expected_active$trait, status$trait), , drop = FALSE]
    if (any(observed$skipped != "False") ||
        any(as.integer(observed$usable_n) != expected_active$usable_n) ||
        any(as.integer(observed$cases) != expected_active$cases) ||
        any(as.integer(observed$controls) != expected_active$controls)) {
      die("production Phase 2 singleton counts differ from the independently established acceptance counts")
    }
    tud <- status[status$trait == "co_ptsd_tud", , drop = FALSE]
    if (nrow(tud) != 1L || tud$skipped[[1]] != "True" || as.integer(tud$usable_n[[1]]) != 416L ||
        !grepl("below_phase2_thresholds:n=416<min_n=500", tud$skip_reason[[1]], fixed = TRUE)) {
      die("production TUD trait was not preserved as the expected 416-sample Phase 2 skip")
    }
  }

  for (trait in trait_ids) {
    status_row <- status[status$trait == trait, , drop = FALSE]
    expected_group <- paste0(status_row$trait_type[[1]], "__", trait)
    if (!identical(status_row$group[[1]], expected_group)) {
      die("Phase 2 trait has a non-stable group ID: ", trait)
    }
    group <- status_row$group[[1]]
    group_dir <- file.path(phase2_dir, "groups", group)
    group_summary_path <- file.path(group_dir, paste0(group, ".trait_summary.tsv"))
    keep_path <- file.path(group_dir, paste0(group, ".keep.txt"))
    require_file(group_summary_path, "missing singleton Phase 2 summary")
    require_file(keep_path, "missing singleton Phase 2 keep")
    group_summary <- read_tsv(group_summary_path)
    if (nrow(group_summary) != 1L || group_summary$trait[[1]] != trait || group_summary$group[[1]] != group) {
      die("Phase 2 group summary is not singleton: ", group_summary_path)
    }

    pan_dir <- file.path(results, "gwas", trait, "PAN")
    public_stats <- file.path(pan_dir, paste0(trait, ".PAN.", build, ".regenie"))
    public_summary_path <- file.path(pan_dir, paste0(trait, ".PAN.", build, ".phase2_summary.tsv"))
    metrics <- file.path(pan_dir, paste0(trait, ".PAN.", build, ".association_metrics.tsv"))
    top_hits <- file.path(pan_dir, paste0(trait, ".PAN.", build, ".top_hits.tsv"))
    report <- file.path(results, "reports", trait,
      paste0(analysis_name, ".", trait, ".PAN.", build, ".regenie.report.md"))
    require_nonempty_file(public_stats, "missing Phase 2 PAN statistics")
    require_file(public_summary_path, "missing Phase 2 public summary")
    require_file(metrics, "missing Phase 2 association metrics")
    require_file(top_hits, "missing Phase 2 top hits")
    require_nonempty_file(report, "missing Phase 2 report")
    report_text <- readLines(report, warn = FALSE)
    if (!any(grepl(paste0("Trait-specific group: ", group), report_text, fixed = TRUE))) {
      die("Phase 2 report lacks its singleton group: ", report)
    }

    skipped <- identical(status_row$skipped[[1]], "True")
    keep_count <- if (file.info(keep_path)$size > 0) length(id_keys(keep_path, paste(trait, "model keep"))) else 0L
    if (keep_count != as.integer(status_row$keep_count[[1]]) ||
        keep_count != as.integer(status_row$model_sample_count[[1]])) {
      die("Phase 2 status and keep counts differ for trait ", trait)
    }
    if (!identical(sha256_file(keep_path), status_row$model_keep_sha256[[1]])) {
      die("Phase 2 status and keep hashes differ for trait ", trait)
    }

    native_prefix <- file.path(results, "gwas", "PAN", "regenie", "groups", group,
      paste0(group, ".", build, "_", trait, ".regenie"))
    prediction <- file.path(results, "gwas", "PAN", "regenie", "groups", group,
      paste0(group, ".step1_1.loco"))
    if (skipped) {
      if (keep_count != 0L || !nzchar(status_row$skip_reason[[1]])) {
        die("skipped Phase 2 trait lacks an empty keep or precise reason: ", trait)
      }
      if (file.exists(prediction) || file.exists(native_prefix) || file.exists(paste0(native_prefix, ".ids"))) {
        die("skipped Phase 2 trait has fake native REGENIE artifacts: ", trait)
      }
    } else {
      native_ids <- paste0(native_prefix, ".ids")
      require_nonempty_file(prediction, "missing singleton REGENIE Step 1 prediction payload")
      require_nonempty_file(native_prefix, "missing native Phase 2 REGENIE statistics")
      require_nonempty_file(native_ids, "missing native Phase 2 REGENIE sample IDs")
      if (!identical(id_keys(native_ids, paste(trait, "ordinary REGENIE IDs")),
          id_keys(keep_path, paste(trait, "model keep")))) {
        die("ordinary REGENIE IDs differ from the exact model keep: ", trait)
      }
    }
  }

  if (truthy(config$remeta$enabled %||% FALSE)) {
    manifest_path <- file.path(results, "remeta", "export",
      paste0(analysis_name, ".", build, ".remeta_manifest.tsv"))
    require_file(manifest_path, "missing ReMeta cohort manifest")
    remeta_manifest <- read_tsv(manifest_path)
    if (!any(remeta_manifest$key == "manifest_schema" & remeta_manifest$value == "remeta_cohort_export_v2")) {
      die("ReMeta cohort manifest does not use schema v2")
    }
    if (any(grepl("(^|[[:space:]/_:])(bt|qt)_g[0-9]+($|[[:space:]/_:])",
        paste(remeta_manifest$key, remeta_manifest$value)))) {
      die("ReMeta cohort manifest references an obsolete shared Phase 2 group")
    }
    for (trait in trait_ids) {
      status_row <- status[status$trait == trait, , drop = FALSE]
      group <- status_row$group[[1]]
      skipped <- identical(status_row$skipped[[1]], "True")
      htp <- file.path(results, "remeta", "export", build, "htp", paste0(trait, ".PAN.regenie.gz"))
      target_prefix <- file.path(results, "remeta", "work", build, "groups", group, paste0(group, ".target"))
      validation <- file.path(results, "remeta", "work", build, "groups", group, paste0(group, ".validation.ok"))
      if (skipped) {
        if (file.exists(htp) || file.exists(paste0(target_prefix, ".pgen")) || file.exists(validation) ||
            dir.exists(file.path(results, "remeta", "export", build, "ld", group))) {
          die("skipped trait has fake ReMeta export artifacts: ", trait)
        }
        next
      }
      require_nonempty_file(htp, "missing active-trait ReMeta HTP")
      require_file(paste0(target_prefix, ".psam"), "missing active-trait ReMeta target PSAM")
      require_file(validation, "missing active-trait ReMeta validation")
      validation_rows <- read_tsv(validation)
      sample_metrics <- validation_rows$value[validation_rows$key %in% c(
        "sample_count", "ordinary_regenie_sample_count", "rare_regenie_sample_count"
      )]
      if (length(sample_metrics) != 3L || any(as.integer(sample_metrics) != as.integer(status_row$model_sample_count[[1]]))) {
        die("ReMeta validation sample counts differ from the Phase 2 model: ", trait)
      }
      for (chrom in 1:22) {
        prefix <- file.path(results, "remeta", "export", build, "ld", group, paste0("chr", chrom))
        for (suffix in c(".remeta.gene.ld", ".remeta.buffer.ld", ".remeta.ld.idx.gz")) {
          require_file(paste0(prefix, suffix), "missing ReMeta chromosome LD component")
        }
      }
    }
  }
}


# Confirm the final run manifest exists and contains no internal group work products.
run_manifest_path <- file.path(results, "manifests", "run_manifest.tsv")
require_file(run_manifest_path, "missing run manifest")
run_manifest <- read_tsv(run_manifest_path)
internal_pattern <- "results/(gwas/PAN/regenie/groups|qc/phase2_regenie/groups)/"
if (any(grepl(internal_pattern, paste(run_manifest$key, run_manifest$value)))) {
  die("run manifest contains internal or obsolete Phase 2 group artifacts")
}
production_traits <- c("co_ptsd_anysud", "co_ptsd_aud", "co_ptsd_cud", "co_ptsd_oud")
if (all(production_traits %in% trait_ids) && truthy(config$phase2_regenie$enabled %||% FALSE)) {
  output_paths <- c(
    file.path(results, "gwas", production_traits, "PAN", paste0(production_traits, ".PAN.", build, ".regenie")),
    file.path(results, "gwas", production_traits, "PAN", paste0(production_traits, ".PAN.", build, ".top_hits.tsv"))
  )
  hashes <- run_manifest$value[match(paste0("output_sha256:", output_paths), run_manifest$key)]
  if (any(is.na(hashes)) || anyDuplicated(hashes[seq_along(production_traits)]) ||
      anyDuplicated(hashes[length(production_traits) + seq_along(production_traits)])) {
    die("production PAN association or top-hit artifacts remain byte-identical across active traits")
  }
}
cat("Pipeline output tests passed\n")
