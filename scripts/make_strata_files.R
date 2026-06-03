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
run_mode <- config$project$run_mode %||% "test"
max_unassigned_fraction <- as.numeric(config$popmad$max_unassigned_fraction %||% 0.02)
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)


# Restrict assignments to manifest samples and configured ancestries.
sample_key <- paste(samples$FID, samples$IID, sep = "\t")
ancestry_key <- paste(ancestry$FID, ancestry$IID, sep = "\t")
valid <- ancestry_key %in% sample_key & ancestry$ancestry %in% labels
valid_keys <- unique(ancestry_key[valid])
missing_keys <- setdiff(sample_key, valid_keys)
if (length(missing_keys)) {
  parts <- do.call(rbind, strsplit(missing_keys, "\t", fixed = TRUE))
  missing_rows <- data.frame(FID = parts[, 1], IID = parts[, 2], reason = "missing_or_unconfigured_ancestry", stringsAsFactors = FALSE)
} else {
  missing_rows <- data.frame(FID = character(), IID = character(), reason = character())
}
write_tsv(missing_rows, file.path(args$outdir, "unassigned_ancestry.tsv"))

unassigned_fraction <- length(missing_keys) / max(length(sample_key), 1)
if (unassigned_fraction > max_unassigned_fraction) {
  first <- if (length(missing_keys)) paste(head(gsub("\t", " ", missing_keys), 20), collapse = ", ") else ""
  die("unassigned ancestry fraction is ", sprintf("%.2f%%", 100 * unassigned_fraction),
    " (", length(missing_keys), "/", length(sample_key), "), above allowed ",
    sprintf("%.2f%%", 100 * max_unassigned_fraction), "; first samples: ", first)
}
if (length(missing_keys)) {
  cat("WARNING: dropping ", length(missing_keys), " samples with unassigned ancestry (",
    sprintf("%.2f%%", 100 * unassigned_fraction), ")\n", sep = "")
}

ancestry <- ancestry[valid, , drop = FALSE]


# Write one PLINK keep file per configured ancestry label.
for (label in labels) {
  keep <- ancestry[ancestry$ancestry == label, c("FID", "IID"), drop = FALSE]
  write_tsv(keep, file.path(args$outdir, paste0(label, ".keep.tsv")))
}


# Summarize trait counts for every ancestry stratum.
counts <- list()
for (i in seq_len(nrow(traits))) {
  trait <- traits[i, ]
  for (label in labels) {
    keep <- ancestry[ancestry$ancestry == label, , drop = FALSE]
    idx <- match(paste(keep$FID, keep$IID, sep = "\t"), sample_key)
    values <- samples[[trait$phenotype_column]][idx]
    cases <- sum(values == trait$case_value)
    controls <- sum(values == trait$control_value)
    n <- cases + controls
    underpowered <- n < config$warnings$min_n ||
      cases < config$warnings$min_cases ||
      controls < config$warnings$min_controls
    counts[[length(counts) + 1]] <- data.frame(
      trait_id = trait$trait_id,
      ancestry = label,
      n = n,
      cases = cases,
      controls = controls,
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
  underpowered = ifelse(length(missing_keys) > 0, "True", "False"),
  stringsAsFactors = FALSE
)


# Production should not launch GWAS jobs for configured empty case/control cells.
count_rows <- do.call(rbind, counts)
if (identical(run_mode, "production")) {
  checked <- count_rows[count_rows$trait_id != "ALL" & count_rows$ancestry != "UNASSIGNED", , drop = FALSE]
  empty <- checked[checked$n == 0 | checked$cases == 0 | checked$controls == 0, , drop = FALSE]
  if (nrow(empty)) {
    labels <- paste0(empty$trait_id, "/", empty$ancestry,
      " n=", empty$n, " cases=", empty$cases, " controls=", empty$controls)
    die("production has empty GWAS strata or case/control cells: ", paste(labels, collapse = "; "))
  }
}


# Save the combined sample-count summary.
write_tsv(count_rows, file.path(args$outdir, "strata_counts.tsv"))
if (identical(config$inputs$ancestry_mode %||% "", "precomputed")) {
  by_ancestry <- as.data.frame(table(ancestry$ancestry), stringsAsFactors = FALSE)
  ancestry_counts <- data.frame(
    category = "assigned_ancestry",
    population = "ALL",
    ancestry = by_ancestry$Var1,
    n = by_ancestry$Freq,
    stringsAsFactors = FALSE
  )
  ancestry_counts <- rbind(ancestry_counts,
    data.frame(category = "excluded_total", population = "ALL", ancestry = "ALL", n = length(missing_keys)))
  write_tsv(ancestry_counts, file.path(dirname(args$outdir), "ancestry", "precomputed_population_counts.tsv"))
} else {
  write_tsv(data.frame(category = character(), population = character(), ancestry = character(), n = integer()),
    file.path(dirname(args$outdir), "ancestry", "precomputed_population_counts.tsv"))
}
cat("Wrote strata files to", args$outdir, "\n")
