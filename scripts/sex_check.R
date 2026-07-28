#!/usr/bin/env Rscript

# Run PLINK2 genetic sex checks and write downstream keep/remove files.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Parse PLINK output prefixes and downstream QC paths.
args <- parse_args(defaults = list(threads = "1"))
require_args(args, c("config", "plink-out-prefix", "sexcheck-out", "remove-out", "keep-out", "summary-out"))


# Load config and manifest sample IDs.
config <- load_config(args$config)
settings <- config$sex_check %||% list()
action <- settings$action %||% "warn"
qc <- config$qc
pruning <- config$relatedness$ld_prune
snps_only_acgt <- truthy(qc$snps_only_acgt %||% TRUE)
source_chromosomes <- "X,Y,XY,PAR1,PAR2"
check_chromosomes <- "X,Y"
ld_prune_values <- c(
  window = as.character(pruning$window %||% "500kb"),
  step = as.character(pruning$step %||% 1),
  r2 = as.character(pruning$r2 %||% 0.2)
)
samples <- read_tsv(config$inputs$sample_manifest)
require_columns(samples, c("FID", "IID", "sex"), "sample manifest")
require_unique_ids(samples, "sample manifest")
sample_keys <- paste(samples$FID, samples$IID, sep = "\t")
plink_out_prefix <- args[["plink-out-prefix"]]
compact_prefix <- paste0(plink_out_prefix, ".sex_markers")
prune_prefix <- paste0(plink_out_prefix, ".sex_marker_prune")
prune_in <- paste0(prune_prefix, ".prune.in")
prune_out <- paste0(prune_prefix, ".prune.out")
manifest_keep <- paste0(plink_out_prefix, ".manifest.keep.txt")


# Extra summary fields are appended to the established metric/value schema.
sex_marker_metrics <- list(
  source_x_variants = "NA",
  source_y_variants = "NA",
  source_xy_variants = "NA",
  source_par1_variants = "NA",
  source_par2_variants = "NA",
  post_qc_variants = "NA",
  post_qc_x_variants = "NA",
  post_qc_y_variants = "NA",
  pruned_variants = "NA",
  par_handling = "not_run",
  marker_maf_min = as.character(qc$maf_min),
  marker_geno_missing_max = as.character(qc$geno_missing_max),
  marker_snps_only_acgt = ifelse(snps_only_acgt, "true", "false"),
  marker_max_alleles = "2",
  marker_rm_dup = "exclude-all",
  source_chromosomes = source_chromosomes,
  check_chromosomes = check_chromosomes,
  ld_prune_window = ld_prune_values[["window"]],
  ld_prune_step = ld_prune_values[["step"]],
  ld_prune_r2 = ld_prune_values[["r2"]]
)


# Read the small sample table so the manifest IDs can be rewritten to the exact
# FID/IID representation used by the source genotype.
read_genotype_ids <- function(block) {
  kind <- tolower(block$type)
  if (kind == "pgen") {
    rows <- read_tsv(paste0(block$prefix, ".psam"))
    return(table_sample_ids(rows, "sex-check genotype PSAM", missing_fid = "zero"))
  }
  if (kind == "bed") {
    rows <- read.table(paste0(block$prefix, ".fam"), stringsAsFactors = FALSE,
      quote = "", comment.char = "")
    if (ncol(rows) < 2 || !nrow(rows)) die("sex-check genotype FAM must contain sample IDs")
    return(data.frame(FID = rows[[1]], IID = rows[[2]], stringsAsFactors = FALSE))
  }
  die("unsupported genotype type '", kind, "'. Use 'pgen' or 'bed'.")
}


write_empty_prune_artifacts <- function() {
  ensure_parent(prune_in)
  for (path in c(prune_in, prune_out)) writeLines(character(), path)
}


# Read named scanner counts once, retaining numeric values for integrity checks.
metadata_counts <- function(summary, metric_map) {
  vapply(metric_map, function(metric) plink_metadata_count(summary, metric), numeric(1))
}


# Split sorted FID/IID keys into a two-column matrix, preserving empty outputs.
split_id_keys <- function(keys) {
  if (!length(keys)) return(matrix(character(), ncol = 2))
  do.call(rbind, strsplit(keys, "\t", fixed = TRUE))
}


# Write sex-check detail, remove, keep, and summary files together.
write_outputs <- function(rows, skipped_reason = "") {
  missing <- c("", "0", "NA", "-9", ".")

  # Build an empty detail file for skipped checks.
  if (!nrow(rows)) {
    out_rows <- data.frame(FID = character(), IID = character(), manifest_sex = character(),
      genotype_pedsex = character(), genetic_sex = character(), plink_status = character(),
      x_f = character(), y_count = character(), y_rate = character(),
      sex_check_status = character(), reason = character())
    problem_keys <- character()
  } else {

    # Compare manifest sex, PLINK PEDSEX, and inferred genetic sex.
    manifest <- setNames(samples$sex, sample_keys)
    keys <- paste(rows$FID, rows$IID, sep = "\t")
    manifest_sex <- unname(manifest[keys])
    manifest_sex[is.na(manifest_sex)] <- ""
    reasons <- character(nrow(rows))
    for (i in seq_len(nrow(rows))) {
      parts <- character()
      if (rows$plink_status[[i]] == "PROBLEM") parts <- c(parts, "plink_pedsex_genetic_mismatch")
      if (!manifest_sex[[i]] %in% missing && !rows$genetic_sex[[i]] %in% missing && manifest_sex[[i]] != rows$genetic_sex[[i]]) {
        parts <- c(parts, "manifest_genetic_mismatch")
      }
      reasons[[i]] <- if (length(parts)) paste(parts, collapse = ";") else "ok"
    }
    status <- ifelse(reasons == "ok", "ok", "problem")
    problem_keys <- keys[status == "problem"]
    out_rows <- data.frame(rows, manifest_sex = manifest_sex, sex_check_status = status, reason = reasons, check.names = FALSE)
    out_rows <- out_rows[c("FID", "IID", "manifest_sex", "genotype_pedsex", "genetic_sex", "plink_status", "x_f", "y_count", "y_rate", "sex_check_status", "reason")]
  }


  # Write detailed status rows and the remove list.
  write_tsv(out_rows, args[["sexcheck-out"]])
  problem_keys <- sort(unique(problem_keys))
  split_keys <- split_id_keys(problem_keys)
  remove <- data.frame(FID = split_keys[, 1], IID = split_keys[, 2], stringsAsFactors = FALSE)
  write_tsv(remove, args[["remove-out"]])


  # Keep all samples unless action=exclude.
  keep_keys <- if (action == "exclude") setdiff(sample_keys, problem_keys) else sample_keys
  keep_keys <- sort(keep_keys)
  keep_parts <- split_id_keys(keep_keys)
  keep <- data.frame(FID = keep_parts[, 1], IID = keep_parts[, 2], stringsAsFactors = FALSE)
  write_tsv(keep, args[["keep-out"]])


  # Summarize check status for reports.
  summary <- data.frame(
    metric = c(
      "enabled", "action", "status", "skipped_reason", "n_samples",
      "n_checked", "n_problems", "n_kept", "n_removed",
      names(sex_marker_metrics)
    ),
    value = c(
      ifelse(identical(skipped_reason, "disabled"), "false", "true"),
      action,
      ifelse(nzchar(skipped_reason), "skipped", "completed"),
      ifelse(nzchar(skipped_reason), skipped_reason, "NA"),
      length(sample_keys),
      nrow(rows),
      length(problem_keys),
      nrow(keep),
      ifelse(action == "exclude", length(problem_keys), 0),
      unlist(sex_marker_metrics, use.names = FALSE)
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(summary, args[["summary-out"]])
  if (action == "exclude" && nrow(keep) == 0 && length(sample_keys) > 0) {
    die("sex_check.action: exclude removed every sample. Review ", args[["sexcheck-out"]],
      " and configure sex_check thresholds or fix manifest/genotype sex coding.")
  }
  length(problem_keys)
}


# Skip cleanly when disabled. Enabled checks inspect and sanitize the source in
# one shared pass, then use those cached chromosome counts below.
if (!truthy(settings$enabled %||% TRUE)) {
  write_empty_prune_artifacts()
  problems <- write_outputs(data.frame(), "disabled")
} else {
  source_input_args <- plink_input_args(
    config$genotypes,
    plink_out_prefix,
    "sex-check genotype input"
  )
  source_summary <- plink_metadata_summary(
    config$genotypes,
    "sex-check genotype input",
    require_alleles = TRUE
  )
  source_counts <- metadata_counts(source_summary, c(
    source_x_variants = "x_variants",
    source_y_variants = "y_variants",
    source_xy_variants = "xy_variants",
    source_par1_variants = "par1_variants",
    source_par2_variants = "par2_variants"
  ))
  sex_marker_metrics[names(source_counts)] <- as.list(as.character(source_counts))

  # XF-based inference needs non-PAR chrX markers. A Y-only or PAR-only source
  # is therefore treated the same as an input with no usable sex markers.
  if (source_counts[["source_x_variants"]] == 0) {
    sex_marker_metrics$par_handling <- "not_run_no_usable_x"
    write_empty_prune_artifacts()
    if (!truthy(settings$allow_no_sex_markers %||% FALSE)) {
      die("sex_check requires usable non-PAR X markers, or explicit sex_check.allow_no_sex_markers: true")
    }
    problems <- write_outputs(data.frame(), "no_sex_chromosome_markers")
    cat("WARNING: sex_check.action=", action,
      " requested but no usable non-PAR X markers were found; keeping all samples\n", sep = "")
  } else {

    # Write only manifest samples, using the source genotype's exact FID/IID
    # representation so genotype-only samples cannot affect MAF or missingness.
    genotype_ids <- read_genotype_ids(config$genotypes)
    require_unique_ids(genotype_ids, "sex-check genotype samples")
    keep_ids <- canonicalize_sample_ids(
      samples[c("FID", "IID")],
      genotype_ids,
      "sample manifest",
      "sex-check genotype samples"
    )
    write_plink_id_file(keep_ids, manifest_keep)

    # Split pseudoautosomal regions when the input has not already encoded
    # PAR1/PAR2. The inferred study build, not a reference-panel label, decides
    # which PLINK boundary definition is used.
    source_par_count <- sum(source_counts[c("source_par1_variants", "source_par2_variants")])
    split_par_args <- character()
    if (source_par_count > 0) {
      sex_marker_metrics$par_handling <- "already_split"
    } else {
      inferred_build <- config$project$inferred_genome_build %||% ""
      split_build <- switch(
        as.character(inferred_build),
        GRCh37 = "b37",
        GRCh38 = "b38",
        ""
      )
      if (!nzchar(split_build)) {
        die("sex check needs project.inferred_genome_build GRCh37 or GRCh38 to split pseudoautosomal regions")
      }
      split_par_args <- c("--split-par", split_build)
      sex_marker_metrics$par_handling <- paste0("split_", split_build)
    }

    # Build a compact, sorted working PGEN. These marker filters are specific
    # to sex inference and do not alter association or ancestry inputs.
    compact_filters <- c(
      "--chr", source_chromosomes,
      split_par_args,
      "--maf", as.character(qc$maf_min),
      "--geno", as.character(qc$geno_missing_max),
      "--max-alleles", "2",
      "--rm-dup", "exclude-all"
    )
    if (snps_only_acgt) {
      compact_filters <- c(compact_filters, "--snps-only", "just-acgt")
    }
    ensure_parent(paste0(compact_prefix, ".pgen"))
    run_command(plink_tool(config), c(
      source_input_args,
      "--keep", manifest_keep,
      compact_filters,
      "--make-pgen", "--sort-vars",
      "--threads", args$threads,
      "--out", compact_prefix
    ))

    compact_summary <- plink_metadata_summary(
      list(type = "pgen", prefix = compact_prefix),
      "sex-check compact PGEN"
    )
    compact_counts <- metadata_counts(compact_summary, c(
      post_qc_variants = "rows_scanned",
      post_qc_x_variants = "x_variants",
      post_qc_y_variants = "y_variants"
    ))
    sex_marker_metrics[names(compact_counts)] <- as.list(as.character(compact_counts))
    if (compact_counts[["post_qc_x_variants"]] == 0) {
      die("sex-check MAF/missingness filters removed every usable non-PAR X marker")
    }

    # Prune only non-PAR X/Y markers, then use that exact audit list for sex
    # inference. No sample-missingness or generic HWE filter is applied here.
    run_command(plink_tool(config), c(
      "--pfile", compact_prefix,
      "--chr", check_chromosomes,
      "--indep-pairwise",
      unname(ld_prune_values),
      "--threads", args$threads,
      "--out", prune_prefix
    ))
    if (!file.exists(prune_in) || file.info(prune_in)$size == 0) {
      die("sex-check LD pruning produced an empty marker set: ", prune_in)
    }
    if (!file.exists(prune_out)) {
      die("sex-check LD pruning did not produce its prune.out audit file: ", prune_out)
    }
    sex_marker_metrics$pruned_variants <- as.character(count_lines(prune_in))

    # Use conventional chrX thresholds unless the config provides
    # cohort-specific values.
    thresholds <- sex_check_threshold_args(settings)
    if (isTRUE(attr(thresholds, "using_defaults"))) {
      cat("No sex_check thresholds configured; using max-female-xf=0.2 and min-male-xf=0.8\n")
    }
    run_command(plink_tool(config), c(
      "--pfile", compact_prefix,
      "--chr", check_chromosomes,
      "--extract", prune_in,
      "--check-sex", thresholds, "cols=fid,pedsex,status,xf,ycount,yrate",
      "--threads", args$threads,
      "--out", plink_out_prefix
    ))

    # Normalize PLINK2 output and require every manifest sample, with no
    # genotype-only extras, to have reached the check.
    native_path <- paste0(plink_out_prefix, ".sexcheck")
    native <- read_tsv(native_path)
    ids <- table_sample_ids(native, paste("PLINK sex-check output", native_path))
    sample_idx <- match_sample_rows(ids, sample_key_map(samples[c("FID", "IID")], "sample manifest"))
    if (anyNA(sample_idx)) {
      unexpected <- paste(ids$FID[is.na(sample_idx)], ids$IID[is.na(sample_idx)])
      die("PLINK sex-check output contains samples absent from the manifest: ",
        paste(head(unexpected, 5), collapse = ", "))
    }
    if (anyDuplicated(sample_idx)) {
      die("PLINK sex-check output contains duplicate manifest samples")
    }
    missing_sample_idx <- setdiff(seq_len(nrow(samples)), sample_idx)
    if (length(missing_sample_idx)) {
      missing_samples <- paste(samples$FID[missing_sample_idx], samples$IID[missing_sample_idx])
      die("PLINK sex-check output is missing manifest samples: ",
        paste(head(missing_samples, 5), collapse = ", "))
    }
    ids$FID <- samples$FID[sample_idx]
    ids$IID <- samples$IID[sample_idx]
    rows <- data.frame(
      FID = ids$FID,
      IID = ids$IID,
      genotype_pedsex = native$PEDSEX %||% "",
      genetic_sex = native$SNPSEX %||% "",
      plink_status = native$STATUS %||% "",
      x_f = if ("XF" %in% names(native)) native$XF else native$F %||% "",
      y_count = native$YCOUNT %||% "",
      y_rate = native$YRATE %||% "",
      stringsAsFactors = FALSE
    )
    problems <- write_outputs(rows)

    # Keep the prune lists and PLINK logs for audit, but remove the compact
    # genotype once the requested sex-check outputs have completed successfully.
    if (action != "fail" || problems == 0) {
      unlink(c(paste0(compact_prefix, c(".pgen", ".pvar", ".psam")), manifest_keep))
    }
  }
}


# Enforce fail mode after all QC files have been written.
if (action == "fail" && problems > 0) die("genetic sex check found ", problems, " problematic samples")
cat("Sex check action=", action, "; problematic samples=", problems, "\n", sep = "")
