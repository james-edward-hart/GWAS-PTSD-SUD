#!/usr/bin/env Rscript

# Create ancestry-stratified keep files and sample-count summaries.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Parse the assignment file and output directory.
args <- parse_args()
require_args(args, c("config", "ancestry-file", "outdir"))


# Load sample, ancestry, and trait metadata.
config <- load_config(args$config)
samples <- read_tsv(config$inputs$sample_manifest)
ancestry <- read_tsv(args[["ancestry-file"]])
traits <- read_tsv(config$inputs$trait_registry)
labels <- unlist(config$analysis$ancestries, use.names = FALSE)
min_stratum_n <- as.integer(config$analysis$min_stratum_n %||% 50)
max_unassigned_fraction <- as.numeric(config$popmad$max_unassigned_fraction %||% 0.07)
if (!is.finite(min_stratum_n) || min_stratum_n < 1) die("analysis.min_stratum_n must be a positive integer")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)


# Restrict assignments to manifest samples and configured ancestries.
sample_key <- paste(samples$FID, samples$IID, sep = "\t")
ancestry_key <- paste(ancestry$FID, ancestry$IID, sep = "\t")
valid <- ancestry_key %in% sample_key & ancestry$ancestry %in% labels
assigned_n <- setNames(rep(0L, length(labels)), labels)
assigned_counts <- table(ancestry$ancestry[valid])
assigned_n[names(assigned_counts)] <- as.integer(assigned_counts)
active_labels <- labels[assigned_n[labels] >= min_stratum_n]
inactive_labels <- setdiff(labels, active_labels)
assigned_summary <- if (length(labels)) {
  paste(paste0(labels, "=", assigned_n[labels]), collapse = ", ")
} else {
  "none"
}
if (!length(active_labels)) {
  die("no ancestry strata meet analysis.min_stratum_n=", min_stratum_n,
    "; assigned configured ancestry counts: ", assigned_summary)
}

matched <- match(sample_key, ancestry_key)
assigned_label <- ifelse(is.na(matched), "", ancestry$ancestry[matched])
configured_assignment <- nzchar(assigned_label) & assigned_label %in% labels
active_assignment <- configured_assignment & assigned_label %in% active_labels
missing_keys <- sample_key[!active_assignment]
if (length(missing_keys)) {
  parts <- do.call(rbind, strsplit(missing_keys, "\t", fixed = TRUE))
  reasons <- ifelse(configured_assignment[!active_assignment], "below_min_stratum_n", "missing_or_unconfigured_ancestry")
  missing_rows <- data.frame(FID = parts[, 1], IID = parts[, 2], reason = reasons, stringsAsFactors = FALSE)
} else {
  missing_rows <- data.frame(FID = character(), IID = character(), reason = character())
}
write_tsv(missing_rows, file.path(args$outdir, "unassigned_ancestry.tsv"))

active_rows <- data.frame(ancestry = active_labels, n = assigned_n[active_labels], stringsAsFactors = FALSE)
write_tsv(active_rows, file.path(args$outdir, "active_ancestries.tsv"))
excluded_rows <- data.frame(
  ancestry = inactive_labels,
  n = assigned_n[inactive_labels],
  reason = ifelse(assigned_n[inactive_labels] == 0, "no_assigned_samples", "below_min_stratum_n"),
  stringsAsFactors = FALSE
)
write_tsv(excluded_rows, file.path(args$outdir, "excluded_ancestries.tsv"))

unassigned_fraction <- length(missing_keys) / max(length(sample_key), 1)
if (unassigned_fraction > max_unassigned_fraction) {
  first <- if (length(missing_keys)) paste(head(gsub("\t", " ", missing_keys), 20), collapse = ", ") else ""
  message <- paste0("inactive or unassigned ancestry fraction is ", sprintf("%.2f%%", 100 * unassigned_fraction),
    " (", length(missing_keys), "/", length(sample_key), "), above allowed ",
    sprintf("%.2f%%", 100 * max_unassigned_fraction),
    "; assigned configured ancestry counts: ", assigned_summary,
    if (length(inactive_labels)) paste0("; inactive labels: ", paste(inactive_labels, collapse = ", ")) else "",
    "; first samples: ", first,
    ". Raise popmad.max_unassigned_fraction only when these dropped samples are expected.")
  if (truthy(config$phase2_regenie$enabled %||% FALSE)) {
    warning(message, "; continuing because phase2_regenie.enabled is true and Phase 2 models excluded samples as UNKNOWN")
  } else {
    die(message)
  }
}
if (length(missing_keys)) {
  cat("WARNING: dropping ", length(missing_keys), " samples outside active ancestry strata (",
    sprintf("%.2f%%", 100 * unassigned_fraction), ")\n", sep = "")
}

ancestry <- ancestry[valid, , drop = FALSE]


# Write one PLINK keep file per configured ancestry label.
for (label in labels) {
  keep <- if (label %in% active_labels) {
    ancestry[ancestry$ancestry == label, c("FID", "IID"), drop = FALSE]
  } else {
    data.frame(FID = character(), IID = character())
  }
  write_tsv(keep, file.path(args$outdir, paste0(label, ".keep.tsv")))
}


# Summarize trait counts for every ancestry stratum.
counts <- list()
for (i in seq_len(nrow(traits))) {
  trait <- traits[i, ]
  is_binary_trait <- !blank(trait$case_value) && !blank(trait$control_value)
  if (blank(trait$case_value) != blank(trait$control_value)) {
    die("trait ", trait$trait_id, " must set both case_value and control_value for binary analysis, or leave both blank for quantitative analysis")
  }
  for (label in labels) {
    keep <- ancestry[ancestry$ancestry == label, , drop = FALSE]
    idx <- match(paste(keep$FID, keep$IID, sep = "\t"), sample_key)
    values <- samples[[trait$phenotype_column]][idx]
    if (is_binary_trait) {
      cases <- sum(values == trait$case_value)
      controls <- sum(values == trait$control_value)
      n <- cases + controls
      underpowered <- n < config$warnings$min_n ||
        cases < config$warnings$min_cases ||
        controls < config$warnings$min_controls
    } else {
      missing <- split_csv(trait$missing_values %||% "")
      observed <- values[!values %in% c(missing, "", "NA", "-9", ".")]
      numeric_observed <- suppressWarnings(as.numeric(observed))
      if (any(is.na(numeric_observed) | !is.finite(numeric_observed))) {
        die("nonnumeric quantitative phenotype values for trait ", trait$trait_id)
      }
      cases <- NA_integer_
      controls <- NA_integer_
      n <- length(observed)
      underpowered <- n < config$warnings$min_n
    }
    counts[[length(counts) + 1]] <- data.frame(
      trait_id = trait$trait_id,
      ancestry = label,
      n = n,
      cases = cases,
      controls = controls,
      active = ifelse(label %in% active_labels, "True", "False"),
      excluded_reason = ifelse(label %in% active_labels, "", excluded_rows$reason[match(label, excluded_rows$ancestry)]),
      underpowered = ifelse(underpowered, "True", "False"),
      stringsAsFactors = FALSE
    )
  }
}
counts[[length(counts) + 1]] <- data.frame(
  trait_id = "ALL",
  ancestry = "UNASSIGNED",
  n = length(missing_keys),
  cases = NA_integer_,
  controls = NA_integer_,
  active = "False",
  excluded_reason = "not_in_active_strata",
  underpowered = ifelse(length(missing_keys) > 0, "True", "False"),
  stringsAsFactors = FALSE
)


# Do not launch GWAS jobs for configured empty case/control cells.
count_rows <- do.call(rbind, counts)
checked <- count_rows[count_rows$trait_id != "ALL" & count_rows$active == "True", , drop = FALSE]
binary_checked <- checked[!is.na(checked$cases) & !is.na(checked$controls), , drop = FALSE]
quant_checked <- checked[is.na(checked$cases) | is.na(checked$controls), , drop = FALSE]
empty <- rbind(
  binary_checked[binary_checked$n == 0 | binary_checked$cases == 0 | binary_checked$controls == 0, , drop = FALSE],
  quant_checked[quant_checked$n == 0, , drop = FALSE]
)
if (nrow(empty)) {
  labels <- ifelse(
    is.na(empty$cases) | is.na(empty$controls),
    paste0(empty$trait_id, "/", empty$ancestry, " n=", empty$n),
    paste0(empty$trait_id, "/", empty$ancestry, " n=", empty$n, " cases=", empty$cases, " controls=", empty$controls)
  )
  die("empty GWAS strata or case/control cells: ", paste(labels, collapse = "; "))
}


# Save the combined sample-count summary.
write_tsv(count_rows, file.path(args$outdir, "strata_counts.tsv"))
cat("Wrote strata files to", args$outdir, "\n")
