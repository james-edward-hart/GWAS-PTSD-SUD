#!/usr/bin/env Rscript

# Build a curated GRCh37/GRCh38 marker panel for genome-build inference.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Shared constants for autosomal marker selection and Ensembl lookups.
autosomes <- as.character(seq_len(22))
ensembl_endpoints <- c(
  GRCh37 = "https://grch37.rest.ensembl.org/variation/homo_sapiens",
  GRCh38 = "https://rest.ensembl.org/variation/homo_sapiens"
)


# Output schemas are fixed so regenerated resources remain comparable.
marker_cols <- c("variant_id", "chrom", "build", "pos", "source", "source_url", "generated_on", "notes")
qc_cols <- c("variant_id", "chrom", "grch37_pos", "grch38_pos", "delta_bp", "selected", "reason")


# Create an empty data frame with a named schema.
empty_df <- function(cols) {
  setNames(data.frame(matrix(nrow = 0, ncol = length(cols)), stringsAsFactors = FALSE), cols)
}


# Keep dbSNP rsIDs only.
is_rsid <- function(value) {
  grepl("^rs[0-9]+$", as.character(value))
}


# Read candidate rsIDs from a PLINK BIM file.
read_bim_candidates <- function(path) {
  rows <- read.table(path, stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (ncol(rows) < 6) die("BIM file has fewer than 6 columns: ", path)
  out <- data.frame(
    variant_id = rows[[2]],
    chrom = clean_chrom(rows[[1]]),
    pos = as.integer(rows[[4]]),
    allele1 = rows[[5]],
    allele2 = rows[[6]],
    stringsAsFactors = FALSE
  )
  out[out$chrom %in% autosomes & is_rsid(out$variant_id), , drop = FALSE]
}


# Read candidate rsIDs from a simple TSV table.
read_tsv_candidates <- function(path) {
  rows <- read_tsv(path)
  require_columns(rows, "variant_id", path)
  chrom <- clean_chrom(rows$chrom %||% "")
  out <- data.frame(
    variant_id = rows$variant_id,
    chrom = ifelse(chrom %in% autosomes, chrom, ""),
    pos = as.integer(rows$pos %||% 0),
    allele1 = rows$allele1 %||% "",
    allele2 = rows$allele2 %||% "",
    stringsAsFactors = FALSE
  )
  out[is_rsid(out$variant_id), , drop = FALSE]
}


# Sample candidates across each chromosome before remote lookup.
select_candidate_pool <- function(rows, markers_per_chrom, multiplier) {
  selected <- rows[FALSE, , drop = FALSE]
  per_chrom <- max(markers_per_chrom * multiplier, markers_per_chrom)
  for (chrom in sort(unique(rows$chrom[rows$chrom %in% autosomes]))) {
    chrom_rows <- rows[rows$chrom == chrom, , drop = FALSE]
    chrom_rows <- chrom_rows[order(chrom_rows$pos), , drop = FALSE]
    if (nrow(chrom_rows) <= per_chrom) {
      selected <- rbind(selected, chrom_rows)
      next
    }
    used <- character()
    min_pos <- min(chrom_rows$pos)
    max_pos <- max(chrom_rows$pos)
    for (index in seq_len(per_chrom)) {
      target <- min_pos + (index * (max_pos - min_pos) / (per_chrom + 1))
      eligible <- chrom_rows[!(chrom_rows$variant_id %in% used), , drop = FALSE]
      choice <- which.min(abs(eligible$pos - target))
      selected <- rbind(selected, eligible[choice, , drop = FALSE])
      used <- c(used, eligible$variant_id[[choice]])
    }
  }
  selected
}


# POST one JSON payload to Ensembl with basic retry behavior.
post_json <- function(url, payload, retries = 4) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) die("R package 'jsonlite' is required")
  curl <- Sys.which("curl")
  if (!nzchar(curl)) die("curl is required to query Ensembl REST")

  body_path <- tempfile("ensembl-payload-")
  out_path <- tempfile("ensembl-response-")
  err_path <- tempfile("ensembl-error-")
  on.exit(unlink(c(body_path, out_path, err_path)), add = TRUE)
  writeLines(jsonlite::toJSON(payload, auto_unbox = TRUE), body_path, useBytes = TRUE)

  for (attempt in seq_len(retries)) {
    code <- system2(curl, c(
      "-sS", "-o", out_path, "-w", "%{http_code}",
      "-X", "POST",
      "-H", "Content-Type: application/json",
      "-H", "Accept: application/json",
      "--data-binary", paste0("@", body_path),
      url
    ), stdout = TRUE, stderr = err_path)
    status <- attr(code, "status") %||% 0L
    code <- tail(code, 1)

    if (identical(status, 0L) && grepl("^2", code)) {
      text <- paste(readLines(out_path, warn = FALSE), collapse = "\n")
      return(jsonlite::fromJSON(text, simplifyVector = FALSE))
    }
    if (code == "429" && attempt < retries) {
      Sys.sleep(2)
      next
    }
    if (attempt < retries) {
      Sys.sleep(2 * attempt)
      next
    }
    details <- paste(readLines(err_path, warn = FALSE), collapse = "\n")
    die("failed POST request after ", retries, " attempts: ", url, " (HTTP ", code, ") ", details)
  }
}


# Fetch one batch of variation records for a genome build.
fetch_variation_batch <- function(ids, build) {
  post_json(ensembl_endpoints[[build]], list(ids = ids))
}


# Fetch all candidate positions from both target builds.
fetch_positions <- function(variant_ids, batch_size) {
  by_build <- list()
  for (build in c("GRCh37", "GRCh38")) {
    records <- list()
    total_batches <- ceiling(length(variant_ids) / batch_size)
    starts <- seq(1, length(variant_ids), by = batch_size)
    for (batch_number in seq_along(starts)) {
      start <- starts[[batch_number]]
      batch <- variant_ids[start:min(start + batch_size - 1, length(variant_ids))]
      cat("Fetching", build, "marker batch", paste0(batch_number, "/", total_batches), "\n")
      batch_records <- fetch_variation_batch(batch, build)
      records[names(batch_records)] <- batch_records
      Sys.sleep(0.2)
    }
    by_build[[build]] <- records
  }
  by_build
}


# Keep only unique autosomal SNP mappings for one build.
clean_mapping <- function(record, build) {
  if (is.null(record) || is.null(record$var_class) || record$var_class != "SNP") return(NULL)
  mappings <- list()
  for (mapping in record$mappings %||% list()) {
    chrom <- clean_chrom(mapping$seq_region_name %||% "")
    if (!(chrom %in% autosomes) || !identical(mapping$assembly_name %||% "", build)) next
    if (!identical(mapping$coord_system %||% "", "chromosome")) next
    if (as.integer(mapping$start %||% 0) != as.integer(mapping$end %||% -1)) next
    mappings[[length(mappings) + 1]] <- list(
      chrom = chrom,
      pos = as.integer(mapping$start),
      alleles = mapping$allele_string %||% ""
    )
  }
  if (length(mappings) != 1) return(NULL)
  mappings[[1]]
}


# Compare fetched mappings and retain markers that distinguish builds.
build_position_rows <- function(candidate_rows, fetched, min_build_delta_bp) {
  candidates <- data.frame()
  qc_rows <- empty_df(qc_cols[-6])

  for (i in seq_len(nrow(candidate_rows))) {
    row <- candidate_rows[i, , drop = FALSE]
    variant_id <- row$variant_id[[1]]
    grch37 <- clean_mapping(fetched$GRCh37[[variant_id]], "GRCh37")
    grch38 <- clean_mapping(fetched$GRCh38[[variant_id]], "GRCh38")

    reason <- "selected_candidate"
    if (is.null(grch37) || is.null(grch38)) {
      reason <- "missing_unique_autosomal_snp_mapping"
    } else if (grch37$chrom != grch38$chrom) {
      reason <- "chromosome_discordant_between_builds"
    } else if (abs(grch37$pos - grch38$pos) < min_build_delta_bp) {
      reason <- "insufficient_position_delta_between_builds"
    }

    qc_rows <- rbind(qc_rows, data.frame(
      variant_id = variant_id,
      chrom = if (!is.null(grch38)) grch38$chrom else row$chrom[[1]],
      grch37_pos = if (!is.null(grch37)) grch37$pos else "",
      grch38_pos = if (!is.null(grch38)) grch38$pos else "",
      delta_bp = if (!is.null(grch37) && !is.null(grch38)) abs(grch37$pos - grch38$pos) else "",
      reason = reason,
      stringsAsFactors = FALSE
    ))
    if (reason != "selected_candidate") next

    candidates <- rbind(candidates, data.frame(
      variant_id = variant_id,
      chrom = grch38$chrom,
      source_chrom = row$chrom[[1]],
      source_pos = row$pos[[1]],
      GRCh37 = grch37$pos,
      GRCh38 = grch38$pos,
      delta_bp = abs(grch37$pos - grch38$pos),
      stringsAsFactors = FALSE
    ))
  }
  list(candidates = candidates, qc_rows = qc_rows)
}


# Choose final markers balanced and spaced across chromosomes.
select_final_markers <- function(candidates, markers_per_chrom, min_spacing_bp) {
  if (!nrow(candidates)) return(candidates)
  selected <- candidates[FALSE, , drop = FALSE]
  for (chrom in as.character(sort(as.integer(unique(candidates$chrom))))) {
    rows <- candidates[candidates$chrom == chrom, , drop = FALSE]
    rows <- rows[order(rows$GRCh38), , drop = FALSE]
    if (nrow(rows) <= markers_per_chrom) {
      selected <- rbind(selected, rows)
      next
    }
    used_ids <- character()
    used_positions <- numeric()
    min_pos <- min(rows$GRCh38)
    max_pos <- max(rows$GRCh38)
    for (index in seq_len(markers_per_chrom)) {
      target <- min_pos + (index * (max_pos - min_pos) / (markers_per_chrom + 1))
      spaced <- vapply(rows$GRCh38, function(pos) !length(used_positions) || all(abs(pos - used_positions) >= min_spacing_bp), logical(1))
      eligible <- rows[!(rows$variant_id %in% used_ids) & spaced, , drop = FALSE]
      if (!nrow(eligible)) break
      choice <- which.min(abs(eligible$GRCh38 - target))
      selected <- rbind(selected, eligible[choice, , drop = FALSE])
      used_ids <- c(used_ids, eligible$variant_id[[choice]])
      used_positions <- c(used_positions, eligible$GRCh38[[choice]])
    }
  }
  selected
}


# Preserve existing marker rows when expanding a panel.
read_existing_markers <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(empty_df(marker_cols))
  rows <- read_tsv(path)
  out <- as.data.frame(setNames(lapply(marker_cols, function(name) rows[[name]] %||% rep("", nrow(rows))), marker_cols), stringsAsFactors = FALSE)
  out$chrom <- clean_chrom(out$chrom)
  out
}


# Write one marker table with build-specific position rows.
write_marker_table <- function(path, existing_rows, selected_rows, candidate_build) {
  rows <- existing_rows
  today <- Sys.Date()
  for (i in seq_len(nrow(selected_rows))) {
    row <- selected_rows[i, , drop = FALSE]
    if (nzchar(candidate_build) && !(candidate_build %in% c("GRCh37", "GRCh38")) && nzchar(as.character(row$source_pos[[1]]))) {
      rows <- rbind(rows, data.frame(
        variant_id = row$variant_id[[1]],
        chrom = ifelse(nzchar(row$source_chrom[[1]]), row$source_chrom[[1]], row$chrom[[1]]),
        build = candidate_build,
        pos = row$source_pos[[1]],
        source = "candidate BIM/PVAR position",
        source_url = "",
        generated_on = today,
        notes = "source-build position preserved from candidate genotype map",
        stringsAsFactors = FALSE
      ))
    }
    for (build in c("GRCh37", "GRCh38")) {
      rows <- rbind(rows, data.frame(
        variant_id = row$variant_id[[1]],
        chrom = row$chrom[[1]],
        build = build,
        pos = row[[build]][[1]],
        source = "Ensembl REST variation endpoint",
        source_url = ensembl_endpoints[[build]],
        generated_on = today,
        notes = "unique autosomal SNP mapping in both GRCh37 and GRCh38",
        stringsAsFactors = FALSE
      ))
    }
  }
  keep <- !duplicated(paste(rows$variant_id, rows$build, sep = "\t"))
  write_tsv(rows[keep, marker_cols, drop = FALSE], path)
}


# Write candidate-level QC and selected status.
write_qc <- function(path, qc_rows, selected_rows) {
  selected_ids <- selected_rows$variant_id
  qc_rows$selected <- qc_rows$variant_id %in% selected_ids
  write_tsv(qc_rows[qc_cols], path)
}


# Parse marker panel build inputs and thresholds.
args <- parse_args(defaults = list(
  "candidate-bim" = "",
  "candidate-tsv" = "",
  "include-existing" = "",
  "markers-per-chrom" = "20",
  "candidate-multiplier" = "6",
  "batch-size" = "100",
  "min-spacing-bp" = "1000000",
  "min-build-delta-bp" = "50",
  "candidate-build" = ""
))
require_args(args, c("out", "qc"))


# Require at least one candidate source.
if (!nzchar(args[["candidate-bim"]]) && !nzchar(args[["candidate-tsv"]])) {
  die("provide --candidate-bim or --candidate-tsv")
}


# Load candidate rsIDs from one or both supported sources.
candidate_rows <- data.frame()
if (nzchar(args[["candidate-bim"]])) candidate_rows <- rbind(candidate_rows, read_bim_candidates(args[["candidate-bim"]]))
if (nzchar(args[["candidate-tsv"]])) candidate_rows <- rbind(candidate_rows, read_tsv_candidates(args[["candidate-tsv"]]))
if (!nrow(candidate_rows)) die("no rsID candidates were found")


# Query enough candidates to survive mapping, spacing, and build-delta filters.
markers_per_chrom <- as.integer(args[["markers-per-chrom"]])
pool <- select_candidate_pool(candidate_rows, markers_per_chrom, as.integer(args[["candidate-multiplier"]]))
variant_ids <- sort(unique(pool$variant_id))
fetched <- fetch_positions(variant_ids, as.integer(args[["batch-size"]]))


# Build the final production marker panel and QC table.
position_rows <- build_position_rows(pool, fetched, as.integer(args[["min-build-delta-bp"]]))
selected <- select_final_markers(position_rows$candidates, markers_per_chrom, as.integer(args[["min-spacing-bp"]]))
if (!nrow(selected)) die("no production genome-build markers passed filters")


# Save marker and QC resources.
existing <- read_existing_markers(args[["include-existing"]])
write_marker_table(args$out, existing, selected, args[["candidate-build"]])
write_qc(args$qc, position_rows$qc_rows, selected)


# Print a compact build summary.
cat("Wrote", nrow(selected), "GRCh37/GRCh38 genome-build markers",
    paste0("(", nrow(selected) * 2, " build rows)"), "to", args$out, "\n")
