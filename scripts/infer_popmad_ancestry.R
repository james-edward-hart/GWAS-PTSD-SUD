#!/usr/bin/env Rscript

# Assign ancestry with POP-MaD-style Mahalanobis distances in projected PC space.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Parse input PC files, POP-MaD settings, and output paths.
args <- parse_args()
require_args(args, c(
  "study-pcs", "reference-pcs", "ancestries", "pcs", "assignments",
  "study-pcs-out", "distances", "reference-outliers", "model-summary",
  "within-pcs", "excluded", "counts", "min-confidence", "outlier-sd"
))


# Resolve requested PC columns and assignment thresholds.
pcs <- paste0("PC", seq_len(as.integer(args$pcs)))
allowed <- split_csv(args$ancestries)
min_confidence <- as.numeric(args[["min-confidence"]])
outlier_sd <- as.numeric(args[["outlier-sd"]])
min_reference_n <- as.integer(args[["min-reference-n"]] %||% 20)


# Load study/reference PCs. All reference super-populations are scored; GWAS
# strata are applied only after the nearest reference group is known.
study <- read_tsv(args[["study-pcs"]])
reference <- read_tsv(args[["reference-pcs"]])
require_columns(study, c("FID", "IID", pcs), "study PC file")
require_columns(reference, c("FID", "IID", "population", "super_population", pcs), "reference PC file")
require_unique_ids(study, "study PC file")
require_unique_ids(reference, "reference PC file")
if (any(!nzchar(reference$population)) || any(!nzchar(reference$super_population))) {
  die("reference PC file contains empty population or super_population labels")
}

mapping <- unique(reference[c("population", "super_population")])
conflicts <- unique(mapping$population[duplicated(mapping$population)])
if (length(conflicts)) die("reference PC file has conflicting population -> super_population mappings: ",
  paste(head(conflicts, 5), collapse = ", "))

empty_excluded <- function() {
  data.frame(FID = character(), IID = character(), best_ancestry = character(), best_population = character(),
    best_distance = numeric(), confidence = numeric(), reason = character())
}

validate_numeric_pcs <- function(rows, label, allow_invalid = FALSE) {
  invalid_rows <- rep(FALSE, nrow(rows))
  invalid_pc <- rep("", nrow(rows))
  for (pc in pcs) {
    value <- suppressWarnings(as.numeric(rows[[pc]]))
    invalid <- is.na(value) | !is.finite(value)
    if (any(invalid) && !allow_invalid) {
      bad <- which(invalid)[[1]]
      die(label, " has missing or non-finite ", pc, " for ", rows$FID[[bad]], " ", rows$IID[[bad]])
    }
    first_invalid <- invalid & !invalid_rows
    invalid_pc[first_invalid] <- pc
    invalid_rows <- invalid_rows | invalid
    rows[[pc]] <- value
  }
  if (!allow_invalid) return(rows)

  excluded <- empty_excluded()
  if (any(invalid_rows)) {
    excluded <- data.frame(
      FID = rows$FID[invalid_rows],
      IID = rows$IID[invalid_rows],
      best_ancestry = "",
      best_population = "",
      best_distance = NA_real_,
      confidence = NA_real_,
      reason = paste0("missing_or_nonfinite_", invalid_pc[invalid_rows]),
      stringsAsFactors = FALSE
    )
    warning("excluded ", nrow(excluded), " study sample(s) with missing or non-finite projected PCs")
  }
  list(valid = rows[!invalid_rows, , drop = FALSE], excluded = excluded)
}
study_check <- validate_numeric_pcs(study, "study PC file", allow_invalid = TRUE)
study <- study_check$valid
invalid_study <- study_check$excluded
reference <- validate_numeric_pcs(reference, "reference PC file")
if (!nrow(reference)) die("reference PC file contains no rows")


# Convert PC columns to a numeric matrix.
vec <- function(df) as.matrix(data.frame(df[pcs], check.names = FALSE))


# Estimate covariance with a small ridge term for stability.
covariance <- function(values, center) {
  if (nrow(values) <= 1) return(diag(length(center)))
  centered <- sweep(values, 2, center)
  cov(centered) + diag(1e-4, length(center))
}


# Compute Mahalanobis distance from each row to a population center.
mahalanobis_distance <- function(values, center, inverse_covariance) {
  diff <- sweep(values, 2, center)
  sqrt(pmax(rowSums((diff %*% inverse_covariance) * diff), 0))
}


# Use a median-plus-SD cutoff to flag distant samples.
distance_cutoff <- function(values, multiplier) {
  center <- median(values)
  variance <- sum((values - center)^2) / max(length(values) - 1, 1)
  center + multiplier * sqrt(variance)
}


# Fit one reference population model.
fit_model <- function(rows) {
  values <- vec(rows)
  center <- colMeans(values)
  inv <- solve(covariance(values, center))
  list(center = center, inverse = inv)
}


# Prepare model output containers.
populations <- unique(reference$population)
outliers <- data.frame(FID = character(), IID = character(), population = character(),
  super_population = character(), mahalanobis_distance = numeric(), cutoff = numeric())
summary_rows <- list()
retained <- list()


# First pass removes distant reference outliers when a population has enough samples.
for (population in populations) {
  rows <- reference[reference$population == population, , drop = FALSE]
  model <- fit_model(rows)
  distances <- mahalanobis_distance(vec(rows), model$center, model$inverse)
  cutoff <- distance_cutoff(distances, outlier_sd)
  drop <- distances > cutoff & nrow(rows) > 10
  if (any(drop)) {
    outliers <- rbind(outliers, data.frame(
      FID = rows$FID[drop], IID = rows$IID[drop],
      population = rows$population[drop], super_population = rows$super_population[drop],
      mahalanobis_distance = distances[drop], cutoff = cutoff,
      stringsAsFactors = FALSE
    ))
  }
  retained[[population]] <- rows[!drop, , drop = FALSE]
  summary_rows[[population]] <- data.frame(
    population = population,
    super_population = rows$super_population[[1]],
    n_input = nrow(rows),
    n_retained = sum(!drop),
    n_outliers = sum(drop),
    reference_outlier_cutoff = cutoff,
    assignment_cutoff = NA_real_,
    pc_count = length(pcs),
    model_status = "pending",
    stringsAsFactors = FALSE
  )
}


# Refit population models after outlier removal.
models <- list()
for (population in populations) {
  rows <- retained[[population]]
  if (nrow(rows) < min_reference_n) {
    summary_rows[[population]]$model_status <- "dropped_low_reference_n"
    summary_rows[[population]]$n_retained <- nrow(rows)
    next
  }
  model <- fit_model(rows)
  distances <- mahalanobis_distance(vec(rows), model$center, model$inverse)
  cutoff <- distance_cutoff(distances, outlier_sd)
  model$super_population <- rows$super_population[[1]]
  model$assignment_cutoff <- cutoff
  models[[population]] <- model
  summary_rows[[population]]$assignment_cutoff <- cutoff
  summary_rows[[population]]$model_status <- "retained"
}

if (!length(models)) die("no POP-MaD reference population models retained at min_reference_n=", min_reference_n)
valid_super <- unique(vapply(models, `[[`, character(1), "super_population"))
missing_super <- setdiff(allowed, valid_super)
if (length(missing_super)) {
  die("no valid POP-MaD reference population model retained for configured super-population(s): ",
    paste(missing_super, collapse = ", "))
}
populations <- names(models)


# Prepare assignment and diagnostic output tables.
assignments <- data.frame(FID = character(), IID = character(), ancestry = character(), population = character(),
  mahalanobis_distance = numeric(), method = character(), confidence = numeric(), status = character())
excluded <- invalid_study
distance_rows <- data.frame(FID = character(), IID = character(), population = character(),
  super_population = character(), mahalanobis_distance = numeric())
within <- data.frame(FID = character(), IID = character(), ancestry = character(), stringsAsFactors = FALSE)
for (pc in pcs) within[[pc]] <- character()


# Score every study sample against every reference population.
study_values <- vec(study)
for (i in seq_len(nrow(study))) {
  distances <- data.frame(
    population = populations,
    super_population = vapply(models, `[[`, character(1), "super_population"),
    mahalanobis_distance = vapply(models, function(model) {
      mahalanobis_distance(matrix(study_values[i, ], nrow = 1), model$center, model$inverse)
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
  distances <- distances[order(distances$mahalanobis_distance, distances$population), ]
  distance_rows <- rbind(distance_rows, data.frame(
    FID = study$FID[[i]], IID = study$IID[[i]], distances, stringsAsFactors = FALSE
  ))

  best <- distances[1, ]

  # Assign confidence between ancestry strata, not between fine populations.
  # Nearby populations from the same super-population should strengthen that
  # stratum rather than make the study sample look ambiguous.
  ancestry_distances <- aggregate(
    mahalanobis_distance ~ super_population,
    data = distances,
    FUN = min
  )
  ancestry_distances <- ancestry_distances[order(ancestry_distances$mahalanobis_distance, ancestry_distances$super_population), ]
  best_ancestry <- ancestry_distances[1, ]
  second_distance <- if (nrow(ancestry_distances) > 1) {
    ancestry_distances$mahalanobis_distance[[2]]
  } else {
    Inf
  }
  confidence <- if (!is.finite(second_distance) || second_distance == 0) 1 else 1 - best_ancestry$mahalanobis_distance / second_distance
  confidence <- max(0, min(confidence, 1))


  # Exclude samples that are nearest to an unconfigured stratum, distant, or ambiguous.
  reason <- ""
  if (!best$super_population %in% allowed) {
    reason <- "unconfigured_nearest_super_population"
  } else if (best$mahalanobis_distance > models[[best$population]]$assignment_cutoff) {
    reason <- "outside_reference_population_distance"
  } else if (confidence < min_confidence) {
    reason <- "ambiguous_nearest_ancestries"
  }

  if (nzchar(reason)) {
    excluded <- rbind(excluded, data.frame(
      FID = study$FID[[i]], IID = study$IID[[i]],
      best_ancestry = best$super_population,
      best_population = best$population,
      best_distance = best$mahalanobis_distance,
      confidence = confidence,
      reason = reason,
      stringsAsFactors = FALSE
    ))
  } else {

    # Assign samples that pass distance and confidence filters.
    assignments <- rbind(assignments, data.frame(
      FID = study$FID[[i]], IID = study$IID[[i]],
      ancestry = best$super_population,
      population = best$population,
      mahalanobis_distance = best$mahalanobis_distance,
      method = "POP-MaD",
      confidence = confidence,
      status = "assigned",
      stringsAsFactors = FALSE
    ))
    within <- rbind(within, data.frame(
      FID = study$FID[[i]], IID = study$IID[[i]], ancestry = best$super_population,
      study[i, pcs, drop = FALSE],
      check.names = FALSE
    ))
  }
}


# Write assignments, distances, and model diagnostics.
write_tsv(assignments, args$assignments)
write_tsv(study[c("FID", "IID", pcs)], args[["study-pcs-out"]])
write_tsv(distance_rows, args$distances)
write_tsv(outliers, args[["reference-outliers"]])
write_tsv(do.call(rbind, summary_rows), args[["model-summary"]])
write_tsv(within, args[["within-pcs"]])
write_tsv(excluded, args$excluded)


# Counts summarize assigned populations, assigned ancestries, and exclusion reasons.
counts <- data.frame(category = character(), population = character(), ancestry = character(), n = integer())
if (nrow(assignments)) {
  by_pop <- as.data.frame(table(assignments$population, assignments$ancestry), stringsAsFactors = FALSE)
  by_pop <- by_pop[by_pop$Freq > 0, ]
  counts <- rbind(counts, data.frame(category = "assigned_population", population = by_pop$Var1, ancestry = by_pop$Var2, n = by_pop$Freq))
  by_anc <- as.data.frame(table(assignments$ancestry), stringsAsFactors = FALSE)
  counts <- rbind(counts, data.frame(category = "assigned_ancestry", population = "ALL", ancestry = by_anc$Var1, n = by_anc$Freq))
}
if (nrow(excluded)) {
  by_reason <- as.data.frame(table(excluded$reason), stringsAsFactors = FALSE)
  counts <- rbind(counts, data.frame(category = "excluded_reason", population = "ALL", ancestry = by_reason$Var1, n = by_reason$Freq))
}
counts <- rbind(counts, data.frame(category = "excluded_total", population = "ALL", ancestry = "ALL", n = nrow(excluded)))
write_tsv(counts, args$counts)


# Print a compact summary for the Snakemake log.
cat("Assigned", nrow(assignments), "samples with POP-MaD; excluded", nrow(excluded), "samples\n")
