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
samples <- read_tsv(config$inputs$sample_manifest)
require_columns(samples, c("FID", "IID", "sex"), "sample manifest")
sample_keys <- paste(samples$FID, samples$IID, sep = "\t")


# Split sorted FID/IID keys into a two-column matrix, preserving empty outputs.
split_id_keys <- function(keys) {
  if (!length(keys)) return(matrix(character(), ncol = 2))
  do.call(rbind, strsplit(keys, "\t", fixed = TRUE))
}


# Detect whether the genotype data includes X or Y markers.
has_sex_markers <- function(config) {
  kind <- tolower(config$genotypes$type)
  prefix <- config$genotypes$prefix
  if (kind == "bed") {
    chrom <- read.table(paste0(prefix, ".bim"), stringsAsFactors = FALSE, quote = "", comment.char = "")[[1]]
  } else {
    pvar <- read_tsv(paste0(prefix, ".pvar"))
    chrom_col <- if ("#CHROM" %in% names(pvar)) "#CHROM" else "CHROM"
    chrom <- pvar[[chrom_col]]
  }
  any(tolower(clean_chrom(chrom)) %in% c("23", "24", "x", "y"))
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
    metric = c("enabled", "action", "status", "skipped_reason", "n_samples", "n_checked", "n_problems", "n_kept", "n_removed"),
    value = c(
      ifelse(identical(skipped_reason, "disabled"), "false", "true"),
      action,
      ifelse(nzchar(skipped_reason), "skipped", "completed"),
      ifelse(nzchar(skipped_reason), skipped_reason, "NA"),
      length(sample_keys),
      nrow(rows),
      length(problem_keys),
      nrow(keep),
      ifelse(action == "exclude", length(problem_keys), 0)
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


# Skip cleanly when disabled or when no sex chromosomes are available.
if (!truthy(settings$enabled %||% TRUE)) {
  problems <- write_outputs(data.frame(), "disabled")
} else if (!has_sex_markers(config)) {
  if (!truthy(settings$allow_no_sex_markers %||% FALSE)) {
    die("sex_check requires sex-chromosome markers, or explicit sex_check.allow_no_sex_markers: true")
  }
  problems <- write_outputs(data.frame(), "no_sex_chromosome_markers")
  cat("WARNING: sex_check.action=", action, " requested but no sex-chromosome markers were found; keeping all samples\n", sep = "")
} else {

  # Use conventional chrX thresholds unless the config provides cohort-specific
  # values. This avoids PLINK2's very strict no-threshold sanity-check defaults.
  thresholds <- sex_check_threshold_args(settings)
  if (isTRUE(attr(thresholds, "using_defaults"))) {
    cat("No sex_check thresholds configured; using max-female-xf=0.2 and min-male-xf=0.8\n")
  }

  # Run PLINK2 native sex check.
  command <- c(
    plink_input_args(config$genotypes, args[["plink-out-prefix"]], "sex-check genotype input"),
    "--check-sex", thresholds, "cols=fid,pedsex,status,xf,ycount,yrate",
    "--threads", args$threads,
    "--out", args[["plink-out-prefix"]]
  )
  ensure_parent(paste0(args[["plink-out-prefix"]], ".sexcheck"))
  run_command(plink_tool(config), command)


  # Normalize PLINK2 sex-check output columns.
  native <- read_tsv(paste0(args[["plink-out-prefix"]], ".sexcheck"))
  ids <- table_sample_ids(native, paste("PLINK sex-check output", paste0(args[["plink-out-prefix"]], ".sexcheck")))
  sample_idx <- match_sample_rows(ids, sample_key_map(samples[c("FID", "IID")], "sample manifest"))
  matched <- !is.na(sample_idx)
  ids$FID[matched] <- samples$FID[sample_idx[matched]]
  ids$IID[matched] <- samples$IID[sample_idx[matched]]
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
}


# Enforce fail mode after all QC files have been written.
if (action == "fail" && problems > 0) die("genetic sex check found ", problems, " problematic samples")
cat("Sex check action=", action, "; problematic samples=", problems, "\n", sep = "")
