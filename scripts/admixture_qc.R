#!/usr/bin/env Rscript

# Prepare, run, and report supervised ADMIXTURE as a report-only QC branch.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Parse the first argument as the requested ADMIXTURE QC subtask.
raw <- commandArgs(trailingOnly = TRUE)
if (!length(raw)) die("missing ADMIXTURE QC subtask")
subtask <- raw[[1]]
args <- parse_args(defaults = list(threads = "1", popmad = ""), raw = raw[-1])


# ADMIXTURE harmonization is restricted to autosomal non-palindromic ACGT SNPs.
autosomes <- as.character(seq_len(22))
palindromic <- c("AT", "TA", "CG", "GC")


# Resolve configured ADMIXTURE labels in the exact output order.
admixture_labels <- function(config) {
  labels <- as.character(unlist(config$admixture$labels, use.names = FALSE))
  if (!length(labels)) die("admixture.labels is empty")
  labels
}


# Resolve the configured ADMIXTURE executable.
admixture_tool <- function(config) {
  config$tools$admixture %||% "admixture"
}


# Resolve relative or absolute executable paths before changing working dirs.
resolve_executable <- function(path) {
  if (grepl("/", path, fixed = TRUE)) return(normalizePath(path, mustWork = TRUE))
  path
}


# Translate ADMIXTURE QC settings into PLINK2 filters.
admixture_filters <- function(config) {
  filters <- config$admixture$filters
  out <- character()
  if (truthy(filters$autosome_only %||% TRUE)) out <- c(out, "--autosome")
  if (truthy(filters$snps_only_acgt %||% TRUE)) out <- c(out, "--snps-only", "just-acgt")
  out <- c(out, "--max-alleles", as.character(filters$max_alleles %||% 2))
  if (!is.null(filters$maf_min)) out <- c(out, "--maf", as.character(filters$maf_min))
  if (!is.null(filters$geno_missing_max)) out <- c(out, "--geno", as.character(filters$geno_missing_max))
  if (truthy(filters$remove_duplicate_ids %||% TRUE)) out <- c(out, "--rm-dup", "exclude-all")
  out
}


# Read a PVAR and keep simple autosomal biallelic ACGT variant rows.
read_pvar <- function(prefix_or_path) {
  path <- if (grepl("\\.pvar$", prefix_or_path)) prefix_or_path else paste0(prefix_or_path, ".pvar")
  rows <- read_tsv(path)
  chrom_col <- if ("#CHROM" %in% names(rows)) "#CHROM" else "CHROM"
  require_columns(rows, c(chrom_col, "POS", "ID", "REF", "ALT"), path)
  rows$chrom_clean <- clean_chrom(rows[[chrom_col]])
  rows$REF <- toupper(rows$REF)
  rows$ALT <- toupper(rows$ALT)
  keep <- rows$chrom_clean %in% autosomes &
    nzchar(rows$ID) & rows$ID != "." &
    grepl("^[ACGT]$", rows$REF) &
    grepl("^[ACGT]$", rows$ALT) &
    rows$REF != rows$ALT
  rows <- rows[keep, , drop = FALSE]
  duplicate_ids <- unique(rows$ID[duplicated(rows$ID)])
  if (length(duplicate_ids)) die("duplicate target variant IDs in ", path, ": ", paste(head(duplicate_ids, 5), collapse = ", "))
  rows <- drop_duplicate_variant_mappings(rows, "chrom_clean", "POS", path)
  data.frame(
    ID = rows$ID,
    chrom = rows$chrom_clean,
    pos = as.integer(rows$POS),
    ref = rows$REF,
    alt = rows$ALT,
    stringsAsFactors = FALSE
  )
}


# Read PLINK FAM sample IDs.
read_fam <- function(path) {
  rows <- read.table(path, stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (ncol(rows) < 2) die("FAM file must contain at least FID and IID columns: ", path)
  names(rows)[seq_len(min(6, ncol(rows)))] <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")[seq_len(min(6, ncol(rows)))]
  rows
}


# Read PLINK2 PSAM sample IDs.
read_psam_ids <- function(prefix_or_path) {
  path <- if (grepl("\\.psam$", prefix_or_path)) prefix_or_path else paste0(prefix_or_path, ".psam")
  rows <- read_tsv(path)
  iid_col <- if ("IID" %in% names(rows)) "IID" else if ("#IID" %in% names(rows)) "#IID" else ""
  if (!nzchar(iid_col)) die("PSAM file is missing IID column: ", path)
  fid_col <- if ("#FID" %in% names(rows)) "#FID" else if ("FID" %in% names(rows)) "FID" else ""
  fid <- if (nzchar(fid_col)) rows[[fid_col]] else rows[[iid_col]]
  out <- data.frame(FID = fid, IID = rows[[iid_col]], stringsAsFactors = FALSE)
  out$key <- paste(out$FID, out$IID, sep = "\t")
  out
}


sample_key_variants <- function(fid, iid) {
  unique(paste(c(fid, iid, "0"), iid, sep = "\t"))
}


sample_key_map <- function(ids, label) {
  variants <- mapply(sample_key_variants, ids$FID, ids$IID, SIMPLIFY = FALSE)
  out <- data.frame(
    key = unlist(variants, use.names = FALSE),
    row = rep(seq_len(nrow(ids)), lengths(variants)),
    stringsAsFactors = FALSE
  )
  conflict <- names(which(tapply(out$row, out$key, function(x) length(unique(x)) > 1)))
  if (length(conflict)) die(label, " has ambiguous sample IDs under FID/IID alias matching: ",
    paste(head(gsub("\t", " ", conflict), 5), collapse = ", "))
  out[!duplicated(out$key), , drop = FALSE]
}


match_sample_row <- function(fid, iid, key_map) {
  idx <- match(sample_key_variants(fid, iid), key_map$key)
  idx <- idx[!is.na(idx)]
  if (!length(idx)) return(NA_integer_)
  key_map$row[[idx[[1]]]]
}


# Read a BIM file for variant counts in reports.
read_bim <- function(path) {
  rows <- read.table(path, stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (ncol(rows) < 2) die("BIM file must contain at least chromosome and variant ID columns: ", path)
  rows
}


# Convert a configured genotype block to filtered, sorted PGEN.
convert_genotypes <- function(config, block, out_prefix, threads) {
  ensure_parent(paste0(out_prefix, ".pgen"))
  run_command(plink_tool(config), c(
    plink_input_args(block, out_prefix, "ADMIXTURE genotype input"), admixture_filters(config),
    "--make-pgen", "--sort-vars", "--threads", threads, "--out", out_prefix
  ))
}


# Intersect study/reference variants and report allele, position, and chromosome mismatches.
shared_variants <- function(config, reference_prefix, study_prefix, out, mismatch_report) {
  reference <- read_pvar(reference_prefix)
  study <- read_pvar(study_prefix)
  merged <- merge(reference, study, by = "ID", suffixes = c("_reference", "_study"))
  merged <- merged[order(merged$ID), , drop = FALSE]
  if (!nrow(merged)) die("no ADMIXTURE reference/study variants share IDs after basic SNP filtering")

  reason <- rep("", nrow(merged))
  reason[merged$chrom_reference != merged$chrom_study] <- "chromosome_mismatch"
  reason[reason == "" & merged$pos_reference != merged$pos_study] <- "position_mismatch"

  exclude_pal <- truthy(config$admixture$filters$exclude_palindromic %||% TRUE)
  if (exclude_pal) {
    pal <- paste0(merged$ref_reference, merged$alt_reference) %in% palindromic |
      paste0(merged$ref_study, merged$alt_study) %in% palindromic
    reason[reason == "" & pal] <- "palindromic_snp_excluded"
  }

  allele_match <- (merged$ref_reference == merged$ref_study & merged$alt_reference == merged$alt_study) |
    (merged$ref_reference == merged$alt_study & merged$alt_reference == merged$ref_study)
  reason[reason == "" & !allele_match] <- "allele_mismatch"

  keep <- merged$ID[reason == ""]
  if (!length(keep)) die("no ADMIXTURE variants survived reference/study harmonization")

  mismatches <- merged[reason != "", , drop = FALSE]
  mismatch_rows <- data.frame(
    variant_id = character(),
    reason = character(),
    reference_chrom = character(),
    study_chrom = character(),
    reference_pos = integer(),
    study_pos = integer(),
    reference_ref = character(),
    reference_alt = character(),
    study_ref = character(),
    study_alt = character(),
    stringsAsFactors = FALSE
  )
  if (nrow(mismatches)) {
    mismatch_rows <- data.frame(
      variant_id = mismatches$ID,
      reason = reason[reason != ""],
      reference_chrom = mismatches$chrom_reference,
      study_chrom = mismatches$chrom_study,
      reference_pos = mismatches$pos_reference,
      study_pos = mismatches$pos_study,
      reference_ref = mismatches$ref_reference,
      reference_alt = mismatches$alt_reference,
      study_ref = mismatches$ref_study,
      study_alt = mismatches$alt_study,
      stringsAsFactors = FALSE
    )
  }

  ensure_parent(out)
  writeLines(keep, out)
  write_tsv(mismatch_rows, mismatch_report)
  hard_mismatches <- mismatch_rows$reason %in% c("chromosome_mismatch", "position_mismatch")
  if (any(hard_mismatches)) {
    warning("excluded ", sum(hard_mismatches),
      " ADMIXTURE variants with chromosome/position mismatches; review ", mismatch_report)
  }
  cat("Wrote", length(keep), "shared ADMIXTURE markers;", nrow(mismatch_rows), "variants excluded\n")
}


# Extract a variant list from a PGEN dataset.
extract_variants <- function(config, input_prefix, variants, out_prefix, threads) {
  ensure_parent(paste0(out_prefix, ".pgen"))
  run_command(plink_tool(config), c(
    "--pfile", input_prefix,
    "--extract", variants,
    "--make-pgen", "--sort-vars",
    "--threads", threads,
    "--out", out_prefix
  ))
}


# Write variants that fall inside long-range LD/problem regions.
write_region_exclusions <- function(config, pfile_prefix, out) {
  regions_path <- config$admixture$exclusion_regions %||% ""
  variants <- read_pvar(paste0(pfile_prefix, ".pvar"))
  excluded <- character()
  if (nzchar(regions_path) && file.exists(regions_path)) {
    regions <- read_tsv(regions_path)
    missing <- setdiff(c("chrom", "start", "end"), names(regions))
    if (length(missing)) die("ADMIXTURE exclusion regions are missing columns: ", paste(missing, collapse = ", "))
    regions$chrom <- clean_chrom(regions$chrom)
    regions <- regions[regions$chrom %in% autosomes, , drop = FALSE]
    for (i in seq_len(nrow(regions))) {
      hit <- variants$chrom == regions$chrom[[i]] &
        variants$pos >= as.integer(regions$start[[i]]) &
        variants$pos <= as.integer(regions$end[[i]])
      excluded <- c(excluded, variants$ID[hit])
    }
  }
  ensure_parent(out)
  writeLines(sort(unique(excluded)), out)
  length(unique(excluded))
}


# LD-prune reference markers for supervised ADMIXTURE.
ld_prune <- function(config, pfile_prefix, out_prefix, prune_in, excluded_regions, threads) {
  n_excluded <- write_region_exclusions(config, pfile_prefix, excluded_regions)
  command <- c("--pfile", pfile_prefix)
  if (n_excluded > 0) command <- c(command, "--exclude", excluded_regions)
  pruning <- config$admixture$ld_prune
  command <- c(command,
    "--indep-pairwise",
    as.character(pruning$window %||% 50),
    as.character(pruning$step %||% 10),
    as.character(pruning$r2 %||% 0.1),
    "--threads", threads,
    "--out", out_prefix
  )
  run_command(plink_tool(config), command)
  if (!file.exists(prune_in) || file.info(prune_in)$size == 0) die("ADMIXTURE LD pruning did not produce a non-empty prune.in file")
  minimum <- as.integer(config$admixture$min_pruned_variants %||% 50000)
  observed <- count_lines(prune_in)
  if (observed < minimum) {
    die("only ", observed, " LD-pruned ADMIXTURE variants survived; minimum required is ", minimum)
  }
}


# Merge pruned reference and study samples, then export ADMIXTURE BED format.
merge_for_admixture <- function(config, reference_prefix, study_prefix, out_prefix, threads) {
  ensure_parent(paste0(out_prefix, ".bed"))
  reference_bed_prefix <- paste0(out_prefix, "_reference_pruned_bed")
  study_bed_prefix <- paste0(out_prefix, "_study_pruned_bed")
  pmerge_prefix <- paste0(out_prefix, "_pmerge")

  run_command(plink_tool(config), c(
    "--pfile", reference_prefix,
    "--mind", "0.999999",
    "--make-bed",
    "--indiv-sort", "none",
    "--threads", threads,
    "--out", reference_bed_prefix
  ))
  run_command(plink_tool(config), c(
    "--pfile", study_prefix,
    "--mind", "0.999999",
    "--make-bed",
    "--indiv-sort", "none",
    "--threads", threads,
    "--out", study_bed_prefix
  ))
  run_command(plink1_tool(config), c(
    "--bfile", reference_bed_prefix,
    "--bmerge", paste0(study_bed_prefix, ".bed"), paste0(study_bed_prefix, ".bim"), paste0(study_bed_prefix, ".fam"),
    "--make-bed",
    "--out", out_prefix
  ))
  run_command(plink_tool(config), c(
    "--bfile", out_prefix,
    "--make-pgen",
    "--indiv-sort", "none",
    "--sort-vars",
    "--threads", threads,
    "--out", pmerge_prefix
  ))
}


# Create a metadata lookup with multiple common FID/IID conventions.
metadata_lookup <- function(config, metadata_path) {
  if (blank(metadata_path)) metadata_path <- config$admixture$metadata$path %||% ""
  require_existing_file(metadata_path, "ADMIXTURE reference metadata")
  meta <- read_tsv(metadata_path)
  settings <- config$admixture$metadata
  sample_col <- settings$sample_id_column
  fid_col <- settings$fid_column %||% ""
  population_col <- settings$population_column
  super_col <- settings$super_population_column
  require_columns(meta, c(sample_col, population_col, super_col), "ADMIXTURE reference metadata")
  if (nzchar(fid_col) && !fid_col %in% names(meta)) die("ADMIXTURE reference metadata is missing configured FID column: ", fid_col)

  iid <- meta[[sample_col]]
  fid <- if (nzchar(fid_col)) meta[[fid_col]] else iid
  key_sets <- list(
    paste(fid, iid, sep = "\t"),
    paste(iid, iid, sep = "\t"),
    paste("0", iid, sep = "\t")
  )
  lookup <- data.frame(
    key = unlist(key_sets, use.names = FALSE),
    population = rep(meta[[population_col]], length(key_sets)),
    super_population = rep(meta[[super_col]], length(key_sets)),
    stringsAsFactors = FALSE
  )
  lookup[!duplicated(lookup$key), , drop = FALSE]
}


# Write the supervised ADMIXTURE .pop file and a readable sample annotation table.
write_supervised_pop <- function(config, fam_path, reference_prefix, study_prefix, metadata_path, pop_out, sample_populations) {
  if (blank(metadata_path)) metadata_path <- config$admixture$metadata$path %||% ""
  labels <- admixture_labels(config)
  fam <- read_fam(fam_path)
  fam$key <- paste(fam$FID, fam$IID, sep = "\t")
  reference <- read_psam_ids(reference_prefix)
  study <- read_psam_ids(study_prefix)
  if (any(duplicated(reference$key))) die("ADMIXTURE reference samples contain duplicate FID/IID rows")
  if (any(duplicated(study$key))) die("ADMIXTURE study samples contain duplicate FID/IID rows")
  reference_map <- sample_key_map(reference, "ADMIXTURE reference samples")
  study_map <- sample_key_map(study, "ADMIXTURE study samples")
  lookup <- metadata_lookup(config, metadata_path)

  out <- data.frame(
    FID = fam$FID,
    IID = fam$IID,
    sample_set = "",
    population = "",
    super_population = "",
    admixture_pop = "",
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(out))) {
    reference_row <- match_sample_row(out$FID[[i]], out$IID[[i]], reference_map)
    study_row <- match_sample_row(out$FID[[i]], out$IID[[i]], study_map)
    in_reference <- !is.na(reference_row)
    in_study <- !is.na(study_row)
    if (in_reference && in_study) die("sample appears in both ADMIXTURE reference and study sets: ", out$FID[[i]], " ", out$IID[[i]])
    if (in_reference) {
      idx <- match(sample_key_variants(out$FID[[i]], out$IID[[i]]), lookup$key)
      idx <- idx[!is.na(idx)]
      if (!length(idx)) die("ADMIXTURE reference metadata missing for sample: ", out$FID[[i]], " ", out$IID[[i]])
      idx <- idx[[1]]
      super_population <- lookup$super_population[[idx]]
      if (!super_population %in% labels) {
        die("ADMIXTURE reference sample ", out$FID[[i]], " ", out$IID[[i]],
          " has super-population label not in admixture.labels: ", super_population)
      }
      out$sample_set[[i]] <- "reference"
      out$population[[i]] <- lookup$population[[idx]]
      out$super_population[[i]] <- super_population
      out$admixture_pop[[i]] <- super_population
    } else if (in_study) {
      out$sample_set[[i]] <- "study"
      out$admixture_pop[[i]] <- "-"
    } else {
      die("merged ADMIXTURE FAM sample is absent from both reference and study PSAM files: ", out$FID[[i]], " ", out$IID[[i]])
    }
  }

  ensure_parent(pop_out)
  writeLines(out$admixture_pop, pop_out)
  write_tsv(out, sample_populations)
  cat("Wrote ADMIXTURE .pop for", sum(out$sample_set == "reference"), "reference and",
    sum(out$sample_set == "study"), "study samples\n")
}


# Run ADMIXTURE supervised mode from the BED directory so raw outputs land beside it.
run_admixture <- function(config, bed, k, threads) {
  if (!identical(config$admixture$mode %||% "supervised", "supervised")) die("only supervised ADMIXTURE mode is supported")
  command <- resolve_executable(admixture_tool(config))
  bed_dir <- dirname(bed)
  bed_base <- basename(bed)
  prefix <- sub("\\.bed$", "", bed_base)
  old <- getwd()
  setwd(bed_dir)
  on.exit(setwd(old), add = TRUE)
  run_command(command, c("--supervised", paste0("-j", threads), bed_base, as.character(k)))
  q <- paste0(prefix, ".", k, ".Q")
  p <- paste0(prefix, ".", k, ".P")
  if (!file.exists(q) || !file.exists(p)) die("ADMIXTURE did not create expected Q/P outputs for K=", k)
}


# Read ADMIXTURE Q/P matrices.
read_admixture_matrix <- function(path, expected_cols, label) {
  rows <- read.table(path, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(rows) != expected_cols) die(label, " has ", ncol(rows), " columns; expected ", expected_cols)
  for (i in seq_len(ncol(rows))) rows[[i]] <- as.numeric(rows[[i]])
  rows
}


# Attach POP-MaD comparison columns when assignments are available.
join_popmad <- function(study, popmad_path) {
  study$popmad_ancestry <- ""
  study$popmad_population <- ""
  study$popmad_confidence <- ""
  study$popmad_status <- ""
  study$top_matches_popmad <- ""
  study$comparison_status <- "popmad_not_available"

  available <- nzchar(popmad_path) && file.exists(popmad_path)
  if (!available) return(list(study = study, available = FALSE))

  popmad <- read_tsv(popmad_path)
  require_columns(popmad, c("FID", "IID", "ancestry", "population", "confidence", "status"), "POP-MaD assignments")
  idx <- match(paste(study$FID, study$IID, sep = "\t"), paste(popmad$FID, popmad$IID, sep = "\t"))
  matched <- !is.na(idx)
  if (any(matched)) {
    rows <- popmad[idx[matched], , drop = FALSE]
    study$popmad_ancestry[matched] <- rows$ancestry
    study$popmad_population[matched] <- rows$population
    study$popmad_confidence[matched] <- rows$confidence
    study$popmad_status[matched] <- rows$status
    same <- study$top_component[matched] == rows$ancestry
    study$top_matches_popmad[matched] <- ifelse(same, "True", "False")
    study$comparison_status[matched] <- ifelse(same, "match", "discordant")
  }
  study$comparison_status[!matched] <- "missing_popmad"
  list(study = study, available = TRUE)
}


# Infer ADMIXTURE Q columns from supervised reference samples and reorder to labels.
infer_and_order_q <- function(q, samples, labels) {
  # ADMIXTURE can emit components in arbitrary order. The supervised reference
  # samples identify which numeric component corresponds to each configured label.
  reference <- samples$sample_set == "reference"
  if (!any(reference)) die("ADMIXTURE sample populations contain no reference samples")
  observed_labels <- unique(samples$super_population[reference])
  missing <- setdiff(labels, observed_labels)
  if (length(missing)) die("ADMIXTURE reference samples are missing labels needed to infer Q columns: ",
    paste(missing, collapse = ", "))

  means <- matrix(NA_real_, nrow = length(labels), ncol = ncol(q), dimnames = list(labels, paste0("V", seq_len(ncol(q)))))
  for (label in labels) {
    rows <- reference & samples$super_population == label
    means[label, ] <- colMeans(as.matrix(q[rows, , drop = FALSE]))
  }

  best <- apply(means, 1, which.max)
  if (any(duplicated(best))) {
    ambiguous <- labels[duplicated(best) | duplicated(best, fromLast = TRUE)]
    die("ADMIXTURE Q-label mapping is ambiguous for labels: ", paste(ambiguous, collapse = ", "))
  }
  if (!setequal(as.integer(best), seq_len(ncol(q)))) {
    die("ADMIXTURE Q-label mapping did not use every component exactly once")
  }

  ordered <- q[, as.integer(best[labels]), drop = FALSE]
  names(ordered) <- labels
  ordered
}


# Parse ADMIXTURE output and write readable QC tables.
parse_report <- function(config, q_path, p_path, fam_path, bim_path, pop_path, sample_populations_path,
                         metadata_path, popmad_path, study_out, reference_out, comparison_out,
                         summary_out, report_out) {
  if (blank(metadata_path)) metadata_path <- config$admixture$metadata$path %||% ""
  require_existing_file(metadata_path, "ADMIXTURE reference metadata")
  labels <- admixture_labels(config)
  k <- as.integer(config$admixture$k %||% length(labels))
  if (length(labels) != k) die("admixture.labels length does not match admixture.k")

  q <- read_admixture_matrix(q_path, k, "ADMIXTURE Q file")
  names(q) <- labels
  p <- read_admixture_matrix(p_path, k, "ADMIXTURE P file")
  fam <- read_fam(fam_path)
  bim <- read_bim(bim_path)
  pop <- readLines(pop_path, warn = FALSE)
  samples <- read_tsv(sample_populations_path)
  metadata <- read_tsv(metadata_path)

  if (nrow(q) != nrow(fam)) die("ADMIXTURE Q rows do not match FAM sample count")
  if (length(pop) != nrow(fam)) die("ADMIXTURE .pop rows do not match FAM sample count")
  require_columns(samples, c("FID", "IID", "sample_set", "population", "super_population", "admixture_pop"), "ADMIXTURE sample populations")
  if (nrow(samples) != nrow(fam)) die("sample population rows do not match FAM sample count")
  if (!identical(paste(samples$FID, samples$IID, sep = "\t"), paste(fam$FID, fam$IID, sep = "\t"))) {
    die("ADMIXTURE sample population rows are not in the same order as the merged FAM")
  }
  if (!identical(pop, samples$admixture_pop)) die("ADMIXTURE .pop rows do not match sample population annotations")
  if (any(samples$sample_set == "study" & pop != "-")) die("study samples must be marked '-' in the ADMIXTURE .pop file")

  # Reorder Q after sanity checks so study/reference proportion tables use
  # semantic ancestry labels rather than ADMIXTURE's raw component order.
  q <- infer_and_order_q(q, samples, labels)
  props <- cbind(samples, q)
  top_index <- max.col(as.matrix(q), ties.method = "first")
  props$top_component <- labels[top_index]
  props$top_proportion <- apply(as.matrix(q), 1, max)

  reference <- props[props$sample_set == "reference", , drop = FALSE]
  reference <- reference[c("FID", "IID", "population", "super_population", labels, "top_component", "top_proportion")]

  study <- props[props$sample_set == "study", , drop = FALSE]
  study <- study[c("FID", "IID", labels, "top_component", "top_proportion")]
  joined <- join_popmad(study, popmad_path)
  study <- joined$study

  comparison <- data.frame(
    FID = study$FID,
    IID = study$IID,
    admixture_top_component = study$top_component,
    admixture_top_proportion = study$top_proportion,
    popmad_ancestry = study$popmad_ancestry,
    popmad_population = study$popmad_population,
    popmad_confidence = study$popmad_confidence,
    popmad_status = study$popmad_status,
    top_matches_popmad = study$top_matches_popmad,
    comparison_status = study$comparison_status,
    stringsAsFactors = FALSE
  )

  n_matches <- sum(comparison$comparison_status == "match")
  n_discordant <- sum(comparison$comparison_status == "discordant")
  n_missing_popmad <- sum(comparison$comparison_status == "missing_popmad")
  n_not_available <- sum(comparison$comparison_status == "popmad_not_available")
  n_comparable <- n_matches + n_discordant
  match_rate <- if (n_comparable) n_matches / n_comparable else NA_real_
  mean_study_top <- if (nrow(study)) mean(as.numeric(study$top_proportion)) else NA_real_
  mean_study_props <- if (nrow(study)) {
    vapply(labels, function(label) mean(as.numeric(study[[label]]), na.rm = TRUE), numeric(1))
  } else {
    stats::setNames(rep(NA_real_, length(labels)), labels)
  }
  study_top_counts <- table(factor(study$top_component, levels = labels))
  fmt_prop <- function(value) if (is.finite(value)) sprintf("%.6f", value) else "NA"

  summary <- data.frame(
    metric = c(
      "mode",
      "k",
      "labels",
      "n_total_samples",
      "n_reference_samples",
      "n_study_samples",
      "n_merged_variants",
      "n_admixture_p_rows",
      "n_metadata_rows",
      "mean_study_top_proportion",
      "popmad_available",
      "popmad_matches",
      "popmad_discordant",
      "popmad_missing",
      "popmad_not_available",
      "popmad_comparable_samples",
      "popmad_match_rate",
      paste0("mean_study_proportion_", labels),
      paste0("n_study_top_component_", labels)
    ),
    value = as.character(c(
      config$admixture$mode %||% "supervised",
      k,
      paste(labels, collapse = ","),
      nrow(props),
      nrow(reference),
      nrow(study),
      nrow(bim),
      nrow(p),
      nrow(metadata),
      fmt_prop(mean_study_top),
      ifelse(joined$available, "True", "False"),
      n_matches,
      n_discordant,
      n_missing_popmad,
      n_not_available,
      n_comparable,
      fmt_prop(match_rate),
      vapply(mean_study_props, fmt_prop, character(1)),
      as.integer(study_top_counts[labels])
    )),
    stringsAsFactors = FALSE
  )

  study_prop_lines <- c(
    "| Ancestry | Mean study proportion | Study top-component samples |",
    "| --- | ---: | ---: |",
    paste0("| ", labels, " | ",
      vapply(mean_study_props, fmt_prop, character(1)), " | ",
      as.integer(study_top_counts[labels]), " |")
  )

  report <- c(
    "# ADMIXTURE QC Report", "",
    "This branch is report-only QC. ADMIXTURE proportions do not replace POP-MaD ancestry labels, do not alter keep files, and are not used as GWAS covariates.", "",
    "## Inputs", "",
    paste0("- Mode: ", config$admixture$mode %||% "supervised"),
    paste0("- K: ", k),
    paste0("- Labels: ", paste(labels, collapse = ", ")),
    paste0("- Reference genotype prefix: `", config$admixture$reference_genotypes$prefix, "`"),
    paste0("- Reference metadata: `", metadata_path, "`"),
    paste0("- Exclusion regions: `", config$admixture$exclusion_regions %||% "", "`"), "",
    "## Counts", "",
    paste0("- Reference samples: ", nrow(reference)),
    paste0("- Study samples: ", nrow(study)),
    paste0("- Merged LD-pruned variants: ", nrow(bim)),
    paste0("- ADMIXTURE P rows: ", nrow(p)),
    paste0("- Mean study top ancestry proportion: ", summary$value[summary$metric == "mean_study_top_proportion"]), "",
    "## Study Mean Proportions", "",
    study_prop_lines, "",
    "## POP-MaD Comparison", "",
    paste0("- POP-MaD assignments available: ", ifelse(joined$available, "True", "False")),
    paste0("- Comparable study samples: ", n_comparable),
    paste0("- Matched study samples: ", n_matches),
    paste0("- Discordant study samples: ", n_discordant),
    paste0("- Match rate among comparable samples: ", fmt_prop(match_rate)),
    paste0("- Missing POP-MaD study samples: ", n_missing_popmad), "",
    "## Outputs", "",
    paste0("- Study proportions: `", study_out, "`"),
    paste0("- Reference proportions: `", reference_out, "`"),
    paste0("- POP-MaD comparison: `", comparison_out, "`"),
    paste0("- Run summary: `", summary_out, "`"),
    paste0("- Raw ADMIXTURE artifacts: `", dirname(q_path), "`")
  )

  write_tsv(study, study_out)
  write_tsv(reference, reference_out)
  write_tsv(comparison, comparison_out)
  write_tsv(summary, summary_out)
  ensure_parent(report_out)
  writeLines(report, report_out)
  cat("Wrote ADMIXTURE QC report for", nrow(study), "study samples and", nrow(reference), "reference samples\n")
}


# Load config shared by every subtask.
require_args(args, "config")
config <- load_config(args$config)
threads <- args$threads %||% "1"


# Dispatch to the requested ADMIXTURE QC subtask.
if (subtask == "convert-reference") {
  require_args(args, "out-prefix")
  convert_genotypes(config, config$admixture$reference_genotypes, args[["out-prefix"]], threads)
} else if (subtask == "convert-study") {
  require_args(args, "out-prefix")
  convert_genotypes(config, config$genotypes, args[["out-prefix"]], threads)
} else if (subtask == "shared-variants") {
  require_args(args, c("reference-prefix", "study-prefix", "out", "mismatch-report"))
  shared_variants(config, args[["reference-prefix"]], args[["study-prefix"]], args$out, args[["mismatch-report"]])
} else if (subtask == "extract-variants") {
  require_args(args, c("input-prefix", "variants", "out-prefix"))
  extract_variants(config, args[["input-prefix"]], args$variants, args[["out-prefix"]], threads)
} else if (subtask == "ld-prune") {
  require_args(args, c("pfile-prefix", "out-prefix", "prune-in", "excluded-regions"))
  ld_prune(config, args[["pfile-prefix"]], args[["out-prefix"]], args[["prune-in"]], args[["excluded-regions"]], threads)
} else if (subtask == "merge") {
  require_args(args, c("reference-prefix", "study-prefix", "out-prefix"))
  merge_for_admixture(config, args[["reference-prefix"]], args[["study-prefix"]], args[["out-prefix"]], threads)
} else if (subtask == "write-pop") {
  require_args(args, c("fam", "reference-prefix", "study-prefix", "pop-out", "sample-populations"))
  write_supervised_pop(config, args$fam, args[["reference-prefix"]], args[["study-prefix"]], args$metadata,
    args[["pop-out"]], args[["sample-populations"]])
} else if (subtask == "run-admixture") {
  require_args(args, c("bed", "k"))
  run_admixture(config, args$bed, as.integer(args$k), threads)
} else if (subtask == "parse-report") {
  require_args(args, c(
    "q", "p", "fam", "bim", "pop", "sample-populations",
    "study-out", "reference-out", "comparison-out", "summary-out", "report-out"
  ))
  parse_report(config, args$q, args$p, args$fam, args$bim, args$pop, args[["sample-populations"]],
    args$metadata, args$popmad %||% "", args[["study-out"]], args[["reference-out"]],
    args[["comparison-out"]], args[["summary-out"]], args[["report-out"]])
} else {
  die("unknown ADMIXTURE QC subtask: ", subtask)
}
