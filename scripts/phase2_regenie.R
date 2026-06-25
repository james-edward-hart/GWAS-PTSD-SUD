#!/usr/bin/env Rscript

# Phase 2 pooled pan-ancestry GWAS helpers for regenie.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


raw <- commandArgs(trailingOnly = TRUE)
if (!length(raw)) die("missing Phase 2 regenie subtask")
subtask <- raw[[1]]
args <- parse_args(
  defaults = list(threads = "1", keep = "", "stage1-summary" = character(), "stage1-stats" = character()),
  repeated = c("stage1-summary", "stage1-stats"),
  raw = raw[-1]
)


autosomes <- as.character(seq_len(22))
regenie_htp_columns <- c(
  "Name", "Chr", "Pos", "Ref", "Alt", "Cohort", "Model", "Effect",
  "LCI_effect", "UCI_effect", "Pval", "AAF", "Num_Cases", "Cases_Ref",
  "Cases_Het", "Cases_Alt", "Num_Controls", "Controls_Ref",
  "Controls_Het", "Controls_Alt", "Info"
)


read_tsv_no_metadata <- function(path) {
  if (!file.exists(path)) die("tab-delimited file not found: ", path)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^##", lines)]
  if (!length(lines) || !any(nzchar(trimws(lines)))) die("tab-delimited file is empty: ", path)
  header <- lines[nzchar(trimws(lines))][[1]]
  sep <- if (grepl("\t", header, fixed = TRUE)) "\t" else ""
  tryCatch(
    read.table(
      text = paste(lines, collapse = "\n"),
      sep = sep,
      header = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = "",
      na.strings = character()
    ),
    error = function(err) die("could not read tab-delimited file ", path, ": ", conditionMessage(err))
  )
}


read_psam_ids <- function(prefix_or_path, missing_fid = "iid") {
  path <- if (grepl("\\.psam$", prefix_or_path)) prefix_or_path else paste0(prefix_or_path, ".psam")
  rows <- read_tsv_no_metadata(path)
  table_sample_ids(rows, paste("PSAM file", path), missing_fid = missing_fid)
}


ensure_regenie_psam_has_fid <- function(prefix_or_path) {
  path <- if (grepl("\\.psam$", prefix_or_path)) prefix_or_path else paste0(prefix_or_path, ".psam")
  rows <- read_tsv_no_metadata(path)
  iid_col <- if ("IID" %in% names(rows)) "IID" else if ("#IID" %in% names(rows)) "#IID" else ""
  if (!nzchar(iid_col)) die("regenie PSAM file is missing IID/#IID sample ID column: ", path)
  fid_col <- if ("#FID" %in% names(rows)) "#FID" else if ("FID" %in% names(rows)) "FID" else ""
  fid <- if (nzchar(fid_col)) rows[[fid_col]] else rows[[iid_col]]

  out <- data.frame(`#FID` = as.character(fid), IID = as.character(rows[[iid_col]]),
    check.names = FALSE, stringsAsFactors = FALSE)
  extra_cols <- setdiff(names(rows), c("#FID", "FID", "IID", "#IID"))
  for (col in extra_cols) out[[col]] <- rows[[col]]
  write_tsv(out, path)
  invisible(path)
}


read_pvar_variants <- function(prefix_or_path) {
  path <- if (grepl("\\.pvar$", prefix_or_path)) prefix_or_path else paste0(prefix_or_path, ".pvar")
  rows <- read_tsv_no_metadata(path)
  chrom_col <- if ("#CHROM" %in% names(rows)) "#CHROM" else "CHROM"
  require_columns(rows, c(chrom_col, "POS", "ID"), path)
  data.frame(
    ID = rows$ID,
    chrom = clean_chrom(rows[[chrom_col]]),
    pos = suppressWarnings(as.integer(rows$POS)),
    stringsAsFactors = FALSE
  )
}


read_sscore <- function(path, pcs) {
  rows <- read_tsv_no_metadata(path)
  ids <- table_sample_ids(rows, paste("projected score file", path))
  pc_cols <- grep("_AVG$", names(rows), value = TRUE)
  if (length(pc_cols) < pcs) pc_cols <- grep("^PC[0-9]+$", names(rows), value = TRUE)
  if (length(pc_cols) < pcs) die("expected at least ", pcs, " projected PC columns in ", path)
  out <- ids
  for (i in seq_len(pcs)) out[[paste0("PC", i)]] <- rows[[pc_cols[[i]]]]
  out
}


regenie_tool <- function(config) {
  config$tools$regenie %||% "regenie"
}


phase2_htp_cohort_name <- function(config) {
  value <- trimws(as.character(config$phase2_regenie$htp_cohort_name %||% ""))
  if (nzchar(value)) {
    if (!grepl("^[A-Za-z0-9._-]+$", value)) {
      die("phase2_regenie.htp_cohort_name must be blank or contain only letters, numbers, dots, underscores, and hyphens")
    }
    return(value)
  }
  analysis_output_name(config)
}


phase2_enabled <- function(config) {
  truthy(config$phase2_regenie$enabled %||% FALSE)
}


phase2_pc_count <- function(config) {
  as.integer(config$phase2_regenie$global_pcs %||% 10)
}


phase2_default_covariates <- function(config) {
  covars <- as.character(unlist(config$phase2_regenie$default_covariates %||% character(), use.names = FALSE))
  if (length(covars)) return(covars[nzchar(covars)])
  c("age", "age2", "sex", paste0("PC", seq_len(phase2_pc_count(config))))
}


phase2_covariates_for_trait <- function(config, trait_row) {
  covars <- c(
    phase2_default_covariates(config),
    as.character(unlist(config$phase2_regenie$extra_covariates %||% character(), use.names = FALSE)),
    split_csv(if (blank(trait_row$covariates)) "" else trait_row$covariates)
  )
  unique(covars[nzchar(covars)])
}


phase2_step1_hardcall_mac_min <- function(config) {
  value <- config$phase2_regenie$step1$filters$mac_min
  if (blank(value)) return(100L)
  mac_min <- suppressWarnings(as.integer(value))
  if (is.na(mac_min) || mac_min < 1) die("phase2_regenie.step1.filters.mac_min must be a positive integer")
  mac_min
}


phase2_step1_info_min <- function(config) {
  value <- config$phase2_regenie$step1$filters$info_min %||% config$qc$info_min %||% 0.8
  if (blank(value)) return(NA_real_)
  threshold <- suppressWarnings(as.numeric(value))
  if (is.na(threshold) || !is.finite(threshold) || threshold < 0 || threshold > 1) {
    die("phase2_regenie.step1.filters.info_min/qc.info_min must be between 0 and 1")
  }
  threshold
}


trait_type <- function(samples, trait_row) {
  case_blank <- blank(trait_row$case_value)
  control_blank <- blank(trait_row$control_value)
  if (!case_blank && !control_blank) return("bt")
  if (xor(case_blank, control_blank)) {
    die("trait ", trait_row$trait_id, " has only one of case_value/control_value set")
  }
  phenotype_column <- trait_row$phenotype_column
  require_columns(samples, phenotype_column, "sample manifest")
  missing <- split_csv(trait_row$missing_values %||% "")
  values <- samples[[phenotype_column]]
  values <- values[!values %in% c(missing, "", "NA", "-9", ".")]
  if (!length(values)) die("quantitative trait ", trait_row$trait_id, " has no nonmissing values")
  numeric_values <- suppressWarnings(as.numeric(values))
  if (any(is.na(numeric_values) | !is.finite(numeric_values))) {
    bad <- values[which(is.na(numeric_values) | !is.finite(numeric_values))[[1]]]
    die("trait ", trait_row$trait_id, " has blank case/control values but nonnumeric phenotype value '", bad, "'")
  }
  "qt"
}


phase2_groups <- function(config) {
  if (!phase2_enabled(config)) return(list())
  samples <- read_tsv(config$inputs$sample_manifest)
  traits <- read_tsv(config$inputs$trait_registry)
  groups <- list()
  keys <- character()
  type_counts <- list(bt = 0L, qt = 0L)

  for (i in seq_len(nrow(traits))) {
    row <- traits[i, , drop = FALSE]
    id <- row$trait_id[[1]]
    if (!nzchar(id)) next
    type <- trait_type(samples, row)
    covars <- phase2_covariates_for_trait(config, row)
    key <- paste(type, paste(covars, collapse = "\r"), sep = "\t")
    idx <- match(key, keys)
    if (is.na(idx)) {
      type_counts[[type]] <- type_counts[[type]] + 1L
      group <- paste0(type, "_g", type_counts[[type]])
      keys <- c(keys, key)
      groups[[length(groups) + 1L]] <- list(group = group, trait_type = type, covariates = covars, traits = character())
      idx <- length(groups)
    }
    groups[[idx]]$traits <- c(groups[[idx]]$traits, id)
  }
  groups
}


group_info <- function(config, group) {
  for (item in phase2_groups(config)) {
    if (identical(item$group, group)) return(item)
  }
  die("unknown Phase 2 regenie group: ", group)
}


write_phase2_group_manifest <- function(config, out) {
  rows <- data.frame(group = character(), trait_type = character(), traits = character(), covariates = character())
  for (item in phase2_groups(config)) {
    rows <- rbind(rows, data.frame(
      group = item$group,
      trait_type = item$trait_type,
      traits = paste(item$traits, collapse = ","),
      covariates = paste(item$covariates, collapse = ","),
      stringsAsFactors = FALSE
    ))
  }
  write_tsv(rows, out)
}


phase2_filter_args <- function(config, branch) {
  filters <- config$phase2_regenie[[branch]]$filters
  out <- character()
  if (truthy(filters$autosome_only %||% TRUE)) out <- c(out, "--autosome")
  if (truthy(filters$snps_only_acgt %||% TRUE)) out <- c(out, "--snps-only", "just-acgt")
  out <- c(out, "--max-alleles", as.character(filters$max_alleles %||% 2))
  mac_min <- if (identical(branch, "step1")) phase2_step1_hardcall_mac_min(config) else filters$mac_min
  if (!blank(mac_min %||% "")) out <- c(out, "--mac", as.character(mac_min))
  if (!blank(filters$maf_min %||% "")) out <- c(out, "--maf", as.character(filters$maf_min))
  if (!blank(filters$geno_missing_max %||% "")) out <- c(out, "--geno", as.character(filters$geno_missing_max))
  if (truthy(filters$remove_duplicate_ids %||% TRUE)) out <- c(out, "--rm-dup", "exclude-all")
  out
}


parse_info_field_metric <- function(value, keys) {
  raw <- as.character(value)
  out <- rep(NA_real_, length(raw))
  missing <- raw %in% c("", ".", "NA")
  for (key in keys) {
    pattern <- paste0("(^|;)", key, "=([^;]+)")
    hit <- regexec(pattern, raw, ignore.case = TRUE)
    parts <- regmatches(raw, hit)
    matched <- lengths(parts) > 0
    if (any(matched)) {
      idx <- which(matched)
      parsed <- suppressWarnings(as.numeric(vapply(parts[matched], function(x) x[[3]], character(1))))
      fill <- is.na(out[idx])
      out[idx[fill]] <- parsed[fill]
    }
  }
  out[missing] <- NA_real_
  out
}


pvar_info_metric <- function(rows) {
  require_columns(rows, "ID", "PVAR file")
  names_upper <- toupper(names(rows))
  direct_candidates <- c("R2", "INFO", "MACH_R2", "MINIMAC3_R2", "IMPUTE_INFO", "IMPUTE2_INFO", "INFO_SCORE", "RSQ")
  info_keys <- c("R2", "INFO", "MACH_R2", "MINIMAC3_R2", "IMPUTE_INFO", "IMPUTE2_INFO", "INFO_SCORE", "RSQ")
  for (candidate in direct_candidates) {
    idx <- match(candidate, names_upper)
    if (is.na(idx)) next
    col <- names(rows)[[idx]]
    value <- suppressWarnings(as.numeric(as.character(rows[[col]])))
    if (any(is.finite(value))) return(list(name = col, value = value))
    if (identical(candidate, "INFO")) {
      value <- parse_info_field_metric(rows[[col]], info_keys)
      if (any(is.finite(value))) return(list(name = paste0(col, ":", "INFO/R2"), value = value))
    }
  }
  NULL
}


phase2_step1_info_filter_args <- function(config, pfile_prefix, out_prefix) {
  threshold <- phase2_step1_info_min(config)
  if (!is.finite(threshold)) return(character())

  pvar_path <- paste0(pfile_prefix, ".pvar")
  rows <- read_tsv_no_metadata(pvar_path)
  metric <- pvar_info_metric(rows)
  if (is.null(metric)) {
    cat("Step 1 INFO/R2 marker filter: skipped; no INFO/R2 column or INFO key found in ", pvar_path, "\n", sep = "")
    return(character())
  }

  value <- metric$value
  low <- is.finite(value) & value < threshold
  if (!any(low)) {
    cat("Step 1 INFO/R2 marker filter: active on ", metric$name,
      "; threshold=", threshold, "; removed=0\n", sep = "")
    return(character())
  }

  pass_ids <- as.character(rows$ID[!low])
  pass_ids <- pass_ids[nzchar(pass_ids) & pass_ids != "." & pass_ids != "NA"]
  if (!length(pass_ids)) die("Step 1 INFO/R2 marker filter removed all variants from ", pvar_path)

  pass_path <- paste0(out_prefix, ".info_r2.pass.snplist")
  excluded_path <- paste0(out_prefix, ".info_r2.excluded.tsv")
  ensure_parent(pass_path)
  writeLines(pass_ids, pass_path)
  write_tsv(data.frame(
    variant_id = as.character(rows$ID[low]),
    info_metric = metric$name,
    info_value = value[low],
    info_min = threshold,
    exclusion_reason = "info_r2_below_min",
    stringsAsFactors = FALSE
  ), excluded_path)
  cat("Step 1 INFO/R2 marker filter: active on ", metric$name,
    "; threshold=", threshold, "; removed=", sum(low), "\n", sep = "")
  c("--extract", pass_path)
}


write_region_exclusions <- function(config, pfile_prefix, out) {
  regions_path <- config$ancestry_reference$exclusion_regions %||% ""
  variants <- read_pvar_variants(pfile_prefix)
  excluded <- character()
  if (nzchar(regions_path) && file.exists(regions_path)) {
    regions <- read_tsv(regions_path)
    chrom_col <- if ("chrom" %in% names(regions)) "chrom" else "#chrom"
    require_columns(regions, c(chrom_col, "start", "end"), "Phase 2 exclusion regions")
    regions$chrom <- clean_chrom(regions[[chrom_col]])
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


blocked_regenie_flags <- c(
  "--step", "--bed", "--pgen", "--bgen", "--sample", "--bgi", "--keep", "--remove",
  "--extract", "--exclude", "--extract-or", "--exclude-or", "--phenoFile", "--phenoCol", "--phenoColList",
  "--phenoExcludeList",
  "--covarFile", "--covarCol", "--covarColList", "--catCovarList", "--pred",
  "--covarExcludeList",
  "--bt", "--qt", "--out", "--bsize", "--lowmem", "--lowmem-prefix", "--threads",
  "--firth", "--approx", "--pThresh", "--minMAC", "--minINFO", "--apply-rint",
  "--gz", "--htp", "--no-split", "--chr", "--chrList", "--range", "--cc12", "--force-qt",
  "--strict", "--force-impute"
)


extra_regenie_options <- function(value) {
  value <- trimws(as.character(value %||% ""))
  if (!nzchar(value)) return(character())
  parts <- strsplit(value, "\\s+")[[1]]
  for (part in parts) {
    flag <- sub("=.*$", "", part)
    if (flag %in% blocked_regenie_flags) {
      die("Phase 2 regenie option pass-through cannot override pipeline-owned option: ", flag)
    }
  }
  parts
}


regenie_type_args <- function(type) {
  if (identical(type, "bt")) return("--bt")
  "--qt"
}


regenie_common_inputs <- function(prefix, pheno, traits, covar, covars) {
  command <- c("--pgen", prefix, "--phenoFile", pheno)
  if (length(traits)) command <- c(command, "--phenoColList", paste(traits, collapse = ","))
  if (length(covars)) command <- c(command, "--covarFile", covar, "--covarColList", paste(covars, collapse = ","))
  command
}


regenie_p_values <- function(rows) {
  if ("P" %in% names(rows)) return(suppressWarnings(as.numeric(rows$P)))
  if ("p" %in% names(rows)) return(suppressWarnings(as.numeric(rows$p)))
  if ("Pval" %in% names(rows)) return(suppressWarnings(as.numeric(rows$Pval)))
  if ("LOG10P" %in% names(rows)) return(10 ^ -suppressWarnings(as.numeric(rows$LOG10P)))
  if ("log10p" %in% names(rows)) return(10 ^ -suppressWarnings(as.numeric(rows$log10p)))
  rep(NA_real_, nrow(rows))
}


genomic_lambda <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) return(NA_real_)
  chi <- suppressWarnings(qchisq(p, df = 1, lower.tail = FALSE))
  chi <- chi[is.finite(chi)]
  if (!length(chi)) return(NA_real_)
  median(chi, na.rm = TRUE) / qchisq(0.5, df = 1, lower.tail = FALSE)
}


metric_value <- function(df, key, default = "NA") {
  if (!nrow(df) || !"metric" %in% names(df) || !"value" %in% names(df)) return(default)
  hit <- df$value[df$metric == key]
  if (length(hit)) hit[[1]] else default
}


prepare_pan_genotypes <- function(config, sex_keep, assignments_path, excluded_path, out_prefix, keep_out, ancestry_out, summary_out, threads) {
  samples <- read_tsv(config$inputs$sample_manifest)
  require_columns(samples, c("FID", "IID"), "sample manifest")
  sex_ids <- read_id_file(sex_keep, "sex-check keep file")
  sex_idx <- match_sample_rows(samples[c("FID", "IID")], sample_key_map(sex_ids, "sex-check keep file"))
  initial <- samples[!is.na(sex_idx), c("FID", "IID"), drop = FALSE]
  if (!nrow(initial)) die("no samples passed sex-check for Phase 2")

  initial_keep <- paste0(out_prefix, ".sex_checked.keep.txt")
  write_plink_id_file(initial, initial_keep)
  run_command(plink_tool(config), c(
    plink_input_args(config$genotypes, out_prefix, "Phase 2 PAN genotype input"),
    "--keep", initial_keep,
    "--mind", as.character(config$qc$sample_missing_max %||% 0.05),
    "--make-pgen", "--sort-vars", "--threads", threads,
    "--out", out_prefix
  ))

  final <- read_psam_ids(out_prefix)
  write_tsv(final, keep_out)

  assignments <- if (file.exists(assignments_path)) read_tsv(assignments_path) else data.frame()
  excluded <- if (file.exists(excluded_path)) read_tsv(excluded_path) else data.frame()
  ancestry <- rep("UNKNOWN", nrow(final))
  population <- rep("", nrow(final))
  reason <- rep("missing_popmad_assignment", nrow(final))
  confidence <- rep("NA", nrow(final))
  if (nrow(assignments)) {
    require_columns(assignments, c("FID", "IID"), "POP-MaD assignments")
    idx <- match_sample_rows(final, sample_key_map(assignments[c("FID", "IID")], "POP-MaD assignments"))
    hit <- !is.na(idx)
    ancestry[hit] <- assignments$ancestry[idx[hit]]
    population[hit] <- assignments$population[idx[hit]]
    confidence[hit] <- as.character(assignments$confidence[idx[hit]])
    reason[hit] <- "assigned"
  }
  if (nrow(excluded)) {
    require_columns(excluded, c("FID", "IID"), "POP-MaD excluded samples")
    idx <- match_sample_rows(final, sample_key_map(excluded[c("FID", "IID")], "POP-MaD excluded samples"))
    hit <- !is.na(idx) & reason != "assigned"
    reason[hit] <- excluded$reason[idx[hit]]
  }
  ancestry_rows <- data.frame(
    FID = final$FID,
    IID = final$IID,
    phase2_ancestry = ancestry,
    popmad_population = population,
    popmad_confidence = confidence,
    phase2_ancestry_reason = reason,
    stringsAsFactors = FALSE
  )
  write_tsv(ancestry_rows, ancestry_out)

  summary <- data.frame(
    metric = c(
      "manifest_samples", "sex_checked_samples", "phase2_pan_samples",
      "sample_missing_max", "assigned_popmad_samples", "unknown_popmad_samples"
    ),
    value = as.character(c(
      nrow(samples), nrow(initial), nrow(final),
      config$qc$sample_missing_max %||% 0.05,
      sum(ancestry != "UNKNOWN"),
      sum(ancestry == "UNKNOWN")
    )),
    stringsAsFactors = FALSE
  )
  write_tsv(summary, summary_out)
}


prepare_marker_set <- function(config, branch, pfile_prefix, keep, out_prefix, prune_prefix, prune_in, excluded_regions, threads) {
  command <- c("--pfile", pfile_prefix)
  if (nzchar(keep)) {
    # PLINK2 matches FID-less PSAM samples as FID 0 in --keep files.
    command <- c(command, plink_keep_args(keep, out_prefix, paste("Phase 2", branch, "keep file"),
      reference_ids = read_psam_ids(pfile_prefix, missing_fid = "zero"),
      reference_label = paste("Phase 2", branch, "PGEN samples")))
  }
  make_pgen_args <- if (identical(branch, "step1")) {
    c("--make-pgen", "fill-missing-from-dosage", "erase-dosage")
  } else {
    "--make-pgen"
  }
  if (identical(branch, "step1")) {
    command <- c(command, phase2_step1_info_filter_args(config, pfile_prefix, out_prefix))
  }
  command <- c(command, phase2_filter_args(config, branch), make_pgen_args, "--sort-vars", "--threads", threads, "--out", out_prefix)
  ensure_parent(paste0(out_prefix, ".pgen"))
  run_command(plink_tool(config), command)

  n_excluded <- write_region_exclusions(config, out_prefix, excluded_regions)
  prune <- config$phase2_regenie[[branch]]$ld_prune
  default_r2 <- if (identical(branch, "step1")) 0.5 else 0.2
  command <- c("--pfile", out_prefix)
  if (n_excluded > 0) command <- c(command, "--exclude", excluded_regions)
  command <- c(command,
    "--indep-pairwise",
    as.character(prune$window %||% "500kb"),
    as.character(prune$step %||% 1),
    as.character(prune$r2 %||% default_r2),
    "--threads", threads,
    "--out", prune_prefix
  )
  run_command(plink_tool(config), command)
  if (!file.exists(prune_in) || file.info(prune_in)$size == 0) die("Phase 2 ", branch, " LD pruning did not produce a non-empty prune.in file")
  if (identical(branch, "step1")) ensure_regenie_psam_has_fid(out_prefix)
}


fit_global_pca <- function(config, pfile_prefix, variants, out_prefix, threads) {
  pcs <- phase2_pc_count(config)
  pca_args <- c(as.character(pcs), "allele-wts", "vcols=chrom,ref,alt", "approx")
  run_command(plink_tool(config), c(
    "--pfile", pfile_prefix,
    "--extract", variants,
    "--freq", "counts",
    "--pca", pca_args,
    "--threads", threads,
    "--out", out_prefix
  ))
}


score_global_pcs <- function(config, pfile_prefix, variants, weights, frequencies, out_prefix, threads) {
  pcs <- phase2_pc_count(config)
  run_command(plink_tool(config), c(
    "--pfile", pfile_prefix,
    "--extract", variants,
    "--read-freq", frequencies,
    "--score", weights, "2", "5", "header-read", "no-mean-imputation", "variance-standardize",
    "--score-col-nums", paste0("6-", 5 + pcs),
    "--threads", threads,
    "--out", out_prefix
  ))
}


build_group_inputs <- function(config, group, keep_path, pcs_path, pheno_out, covar_out, summary_out, trait_list_out, covar_list_out, keep_plink_out) {
  samples <- read_tsv(config$inputs$sample_manifest)
  traits <- read_tsv(config$inputs$trait_registry)
  keep <- read_id_file(keep_path, "Phase 2 PAN keep file")
  pcs <- read_tsv(pcs_path)
  info <- group_info(config, group)
  covars <- info$covariates

  sample_idx <- match_sample_rows(keep, sample_key_map(samples[c("FID", "IID")], "sample manifest"))
  if (any(is.na(sample_idx))) die("Phase 2 keep sample missing from sample manifest")
  samples <- samples[sample_idx, , drop = FALSE]

  pc_idx <- match_sample_rows(keep, sample_key_map(pcs[c("FID", "IID")], "Phase 2 global PC table"))

  covar <- data.frame(FID = keep$FID, IID = keep$IID, stringsAsFactors = FALSE)
  for (covar_name in covars) {
    if (covar_name %in% names(samples)) {
      value <- samples[[covar_name]]
    } else if (covar_name %in% names(pcs)) {
      value <- ifelse(is.na(pc_idx), "NA", pcs[[covar_name]][pc_idx])
    } else {
      die("Phase 2 covariate '", covar_name, "' is absent from sample manifest and global PC table")
    }
    value[!nzchar(value) | value %in% c(".", "-9")] <- "NA"
    covar[[covar_name]] <- value
  }

  covar_complete <- rep(TRUE, nrow(covar))
  for (covar_name in covars) {
    value <- suppressWarnings(as.numeric(covar[[covar_name]]))
    covar_complete <- covar_complete & is.finite(value)
  }

  pheno <- data.frame(FID = keep$FID, IID = keep$IID, stringsAsFactors = FALSE)
  summary_rows <- data.frame()
  analysis_traits <- character()
  for (trait_id in info$traits) {
    trait <- traits[traits$trait_id == trait_id, , drop = FALSE]
    if (!nrow(trait)) die("unknown trait in Phase 2 group ", group, ": ", trait_id)
    value <- samples[[trait$phenotype_column]]
    missing <- split_csv(trait$missing_values %||% "")
    if (identical(info$trait_type, "bt")) {
      encoded <- ifelse(value == trait$case_value[[1]], "1",
        ifelse(value == trait$control_value[[1]], "0",
          ifelse(value %in% missing, "NA", NA_character_)
        )
      )
      if (any(is.na(encoded))) die("unexpected binary phenotype value for Phase 2 trait ", trait_id)
      usable <- encoded != "NA" & covar_complete
      cases <- sum(encoded[usable] == "1")
      controls <- sum(encoded[usable] == "0")
      n <- cases + controls
      skip <- n < as.integer(config$warnings$min_n %||% 0) ||
        cases < as.integer(config$warnings$min_cases %||% 0) ||
        controls < as.integer(config$warnings$min_controls %||% 0)
      reason <- if (skip) paste0("below_phase2_thresholds:n=", n, ";cases=", cases, ";controls=", controls) else ""
    } else {
      encoded <- ifelse(value %in% c(missing, "", "NA", "-9", "."), "NA", value)
      numeric_value <- suppressWarnings(as.numeric(encoded))
      bad <- encoded != "NA" & (is.na(numeric_value) | !is.finite(numeric_value))
      if (any(bad)) die("nonnumeric quantitative phenotype value for Phase 2 trait ", trait_id)
      usable <- encoded != "NA" & covar_complete
      cases <- NA_integer_
      controls <- NA_integer_
      n <- sum(usable)
      skip <- n < as.integer(config$warnings$min_n %||% 0)
      reason <- if (skip) paste0("below_phase2_thresholds:n=", n) else ""
    }
    if (!skip) {
      pheno[[trait_id]] <- encoded
      analysis_traits <- c(analysis_traits, trait_id)
    }
    summary_rows <- rbind(summary_rows, data.frame(
      group = group,
      trait = trait_id,
      trait_type = info$trait_type,
      covariates = paste(covars, collapse = ","),
      phase2_pan_samples = nrow(keep),
      complete_covariate_samples = sum(covar_complete),
      usable_n = n,
      cases = cases,
      controls = controls,
      skipped = ifelse(skip, "True", "False"),
      skip_reason = reason,
      stringsAsFactors = FALSE
    ))
  }

  write_tsv(pheno, pheno_out)
  write_tsv(covar, covar_out)
  write_tsv(summary_rows, summary_out)
  writeLines(analysis_traits, trait_list_out)
  writeLines(paste(covars, collapse = ","), covar_list_out)
  step1_complete <- covar_complete
  if (length(analysis_traits)) {
    for (trait_id in analysis_traits) {
      step1_complete <- step1_complete & pheno[[trait_id]] != "NA"
    }
  } else {
    step1_complete <- rep(FALSE, nrow(keep))
  }
  write_plink_id_file(keep[step1_complete, , drop = FALSE], keep_plink_out)
}


stage1_pass_union <- function(paths, out, summary_out) {
  ids <- character()
  counts <- data.frame(file = character(), variants = integer())
  for (path in paths) {
    rows <- read_tsv_no_metadata(path)
    id_col <- if ("variant_id" %in% names(rows)) "variant_id" else if ("ID" %in% names(rows)) "ID" else ""
    if (!nzchar(id_col)) die("Stage 1 stats file is missing variant ID column: ", path)
    value <- rows[[id_col]]
    value <- value[nzchar(value) & value != "NA" & value != "."]
    ids <- c(ids, value)
    counts <- rbind(counts, data.frame(file = path, variants = length(unique(value))))
  }
  ids <- sort(unique(ids))
  if (!length(ids)) die("no Stage 1 variants available for Phase 2 union pass list")
  ensure_parent(out)
  writeLines(ids, out)
  totals <- data.frame(file = "UNION", variants = length(ids))
  write_tsv(rbind(counts, totals), summary_out)
}


prepare_assoc_genotypes <- function(config, pfile_prefix, extract, out_prefix, threads) {
  run_command(plink_tool(config), c(
    "--pfile", pfile_prefix,
    "--extract", extract,
    "--geno", as.character(config$qc$geno_missing_max %||% 0.05),
    "--make-pgen", "--sort-vars",
    "--threads", threads,
    "--out", out_prefix
  ))
  ensure_regenie_psam_has_fid(out_prefix)
}


step1_filter_report_paths <- function(out, summary_out = "", excluded_out = "") {
  stem <- sub("\\.snplist$", "", out)
  list(
    summary = if (blank(summary_out)) paste0(stem, ".variant_qc.summary.tsv") else summary_out,
    excluded = if (blank(excluded_out)) paste0(stem, ".variant_qc.excluded.tsv") else excluded_out
  )
}


write_step1_filter_summary <- function(path, raw_count, pass_count, model_sample_count,
                                       hardcall_mac_min, hardcall_variance_min = 0) {
  summary <- data.frame(
    filter_method = "plink2_hardcall_count_qc",
    plink_nonfounders = "True",
    raw_plink_pass_snp_count = as.integer(raw_count),
    hardcall_filter_pass_snp_count = as.integer(pass_count),
    excluded_snp_count = as.integer(raw_count - pass_count),
    model_sample_count = as.integer(model_sample_count),
    hardcall_mac_min = as.integer(hardcall_mac_min),
    hardcall_variance_min = as.numeric(hardcall_variance_min),
    stringsAsFactors = FALSE
  )
  write_tsv(summary, path)
}


empty_step1_excluded_table <- function(path) {
  write_tsv(data.frame(
    variant_id = character(),
    hardcall_ref_ct = numeric(),
    hardcall_alt_ct = numeric(),
    hardcall_mac = numeric(),
    hardcall_n = numeric(),
    hardcall_variance = numeric(),
    exclusion_reason = character(),
    stringsAsFactors = FALSE
  ), path)
}


read_step1_gcount_results <- function(path, raw_ids, hardcall_mac_min, hardcall_variance_min = 0) {
  if (!file.exists(path)) die("PLINK2 Step 1 hardcall-count QC did not produce expected .gcount output: ", path)
  duplicate_raw <- unique(raw_ids[duplicated(raw_ids)])
  if (length(duplicate_raw)) die("Step 1 raw PLINK-filtered variant list has duplicate IDs: ", paste(head(duplicate_raw, 5), collapse = ", "))

  rows <- read_tsv_no_metadata(path)
  required <- c("ID", "HOM_REF_CT", "HET_REF_ALT1_CT", "HOM_ALT1_CT", "MISSING_CT", "OBS_CT")
  missing_cols <- setdiff(required, names(rows))
  if (length(missing_cols)) die("PLINK2 Step 1 .gcount output is missing required columns: ", paste(missing_cols, collapse = ", "))
  rows$ID <- as.character(rows$ID)
  duplicate_gcount <- unique(rows$ID[duplicated(rows$ID)])
  if (length(duplicate_gcount)) die("PLINK2 Step 1 .gcount output has duplicate variant IDs: ", paste(head(duplicate_gcount, 5), collapse = ", "))

  numeric_col <- function(col) {
    value <- suppressWarnings(as.numeric(rows[[col]]))
    bad <- is.na(value) | !is.finite(value) | value < 0
    if (any(bad)) die("PLINK2 Step 1 .gcount column ", col, " contains missing, nonnumeric, or negative values")
    value
  }
  rows$HOM_REF_CT <- numeric_col("HOM_REF_CT")
  rows$HET_REF_ALT1_CT <- numeric_col("HET_REF_ALT1_CT")
  rows$HOM_ALT1_CT <- numeric_col("HOM_ALT1_CT")
  rows$MISSING_CT <- numeric_col("MISSING_CT")
  rows$OBS_CT <- numeric_col("OBS_CT")

  idx <- match(raw_ids, rows$ID)
  result <- data.frame(
    variant_id = raw_ids,
    hardcall_ref_ct = NA_real_,
    hardcall_alt_ct = NA_real_,
    hardcall_mac = NA_real_,
    hardcall_n = NA_real_,
    hardcall_variance = NA_real_,
    pass = FALSE,
    exclusion_reason = "missing_from_plink2_gcount",
    stringsAsFactors = FALSE
  )

  present <- !is.na(idx)
  if (any(present)) {
    hit <- idx[present]
    hom_ref <- rows$HOM_REF_CT[hit]
    het <- rows$HET_REF_ALT1_CT[hit]
    hom_alt <- rows$HOM_ALT1_CT[hit]
    hardcall_n <- hom_ref + het + hom_alt
    ref_ct <- 2 * hom_ref + het
    alt_ct <- het + 2 * hom_alt
    mac <- pmin(ref_ct, alt_ct)
    mean_alt <- rep(0, length(hardcall_n))
    mean_alt_sq <- rep(0, length(hardcall_n))
    observed <- hardcall_n > 0
    mean_alt[observed] <- alt_ct[observed] / hardcall_n[observed]
    mean_alt_sq[observed] <- (het[observed] + 4 * hom_alt[observed]) / hardcall_n[observed]
    variance <- mean_alt_sq - mean_alt^2
    variance[!is.finite(variance) | (variance < 0 & variance > -1e-12)] <- 0

    result$hardcall_ref_ct[present] <- ref_ct
    result$hardcall_alt_ct[present] <- alt_ct
    result$hardcall_mac[present] <- mac
    result$hardcall_n[present] <- hardcall_n
    result$hardcall_variance[present] <- variance

    pass <- is.finite(variance) & variance > hardcall_variance_min & mac >= hardcall_mac_min
    reason <- ifelse(!is.finite(variance) | variance <= hardcall_variance_min, "zero_hardcall_variance",
      ifelse(mac < hardcall_mac_min, "hardcall_mac_below_min", ""))
    result$pass[present] <- pass
    result$exclusion_reason[present] <- reason
  }

  result
}


filter_step1_variants <- function(config, pfile_prefix, extract, keep, trait_list, out, threads,
                                  summary_out = "", excluded_out = "") {
  traits <- readLines(trait_list, warn = FALSE)
  traits <- traits[nzchar(traits)]
  reports <- step1_filter_report_paths(out, summary_out, excluded_out)
  ensure_parent(out)
  hardcall_mac_min <- phase2_step1_hardcall_mac_min(config)
  hardcall_variance_min <- 0
  if (!length(traits)) {
    writeLines(character(), out)
    write_step1_filter_summary(reports$summary, 0L, 0L, 0L, hardcall_mac_min, hardcall_variance_min)
    empty_step1_excluded_table(reports$excluded)
    cat("No analysis traits; wrote empty Step 1 variant outputs and skipped PLINK hardcall-count QC\n")
    return(invisible(TRUE))
  }

  model_sample_count <- nrow(read_id_file(keep, "Phase 2 regenie Step 1 keep file"))
  tmp_prefix <- paste0(sub("\\.snplist$", "", out), ".variant_qc")
  # Reapply Step 1 marker filters after the final group keep file. The pooled
  # marker set can contain SNPs that become too rare for regenie's actual model
  # sample once phenotype/covariate-complete samples are selected.
  run_command(plink_tool(config), c(
    "--pfile", pfile_prefix,
    "--extract", extract,
    "--keep", keep,
    phase2_filter_args(config, "step1"),
    "--nonfounders",
    "--write-snplist",
    "--geno-counts", "cols=chrom,pos,ref,alt1,homref,refalt1,homalt1,missing,nobs",
    "--threads", threads,
    "--out", tmp_prefix
  ))
  snplist <- paste0(tmp_prefix, ".snplist")
  if (!file.exists(snplist) || file.info(snplist)$size == 0) {
    die(
      "Phase 2 regenie Step 1 group variant filter removed all markers after applying ",
      "the group keep file and configured Step 1 filters; check phenotype/covariate-complete ",
      "sample count or relax phase2_regenie.step1.filters"
    )
  }
  raw_ids <- readLines(snplist, warn = FALSE)
  raw_ids <- raw_ids[nzchar(raw_ids)]

  gcount <- paste0(tmp_prefix, ".gcount")
  all_results <- read_step1_gcount_results(gcount, raw_ids, hardcall_mac_min, hardcall_variance_min)
  final_ids <- all_results$variant_id[all_results$pass]
  excluded <- all_results[!all_results$pass, c(
    "variant_id", "hardcall_ref_ct", "hardcall_alt_ct", "hardcall_mac",
    "hardcall_n", "hardcall_variance", "exclusion_reason"
  ), drop = FALSE]

  writeLines(final_ids, out)
  if (nrow(excluded)) write_tsv(excluded, reports$excluded) else empty_step1_excluded_table(reports$excluded)
  write_step1_filter_summary(
    reports$summary, length(raw_ids), length(final_ids), model_sample_count,
    hardcall_mac_min, hardcall_variance_min
  )

  cat("Step 1 PLINK-pass variants:", length(raw_ids), "\n")
  cat("Step 1 PLINK2 hardcall-count QC: enabled\n")
  cat("Step 1 PLINK2 --nonfounders: enabled\n")
  cat("Step 1 model samples:", model_sample_count, "\n")
  cat("Step 1 hardcall MAC minimum:", hardcall_mac_min, "\n")
  cat("Step 1 hardcall variance minimum:", hardcall_variance_min, "\n")
  cat("Hardcall-count QC pass variants:", length(final_ids), "\n")
  cat("Hardcall-count QC excluded variants:", nrow(excluded), "\n")

  if (!file.exists(out) || file.info(out)$size == 0) die("staged Phase 2 regenie Step 1 variant list is empty: ", out)
}


shell_quote <- function(value) {
  shQuote(as.character(value), type = "sh")
}


shell_command_line <- function(command, args) {
  paste(c(shell_quote(command), vapply(args, shell_quote, character(1))), collapse = " ")
}


write_bash_script <- function(path, lines) {
  ensure_parent(path)
  writeLines(lines, path)
  Sys.chmod(path, mode = "0755")
}


mkdir_parent_line <- function(path) {
  paste("mkdir -p", shell_quote(dirname(path)))
}


printf_line <- function(value, path) {
  paste("printf '%s\\n'", shell_quote(value), ">", shell_quote(path))
}


regenie_step1_args <- function(config, group, pfile_prefix, extract, pheno, covar, keep, trait_list, out_prefix, threads) {
  info <- group_info(config, group)
  traits <- readLines(trait_list, warn = FALSE)
  traits <- traits[nzchar(traits)]
  covars <- info$covariates
  command <- c(
    "--step", "1",
    "--pgen", pfile_prefix,
    "--extract", extract,
    "--keep", keep,
    "--phenoFile", pheno,
    "--phenoColList", paste(traits, collapse = ","),
    "--covarFile", covar,
    "--covarColList", paste(covars, collapse = ","),
    regenie_type_args(info$trait_type),
    "--bsize", as.character(config$phase2_regenie$step1_bsize %||% 1000),
    "--lowmem",
    "--lowmem-prefix", paste0(out_prefix, ".lowmem"),
    "--threads", threads
  )
  if (identical(info$trait_type, "qt") && truthy(config$phase2_regenie$apply_rint %||% FALSE)) {
    command <- c(command, "--apply-rint")
  }
  list(
    traits = traits,
    args = c(command, extra_regenie_options(config$phase2_regenie$step1_options %||% ""), "--out", out_prefix)
  )
}


write_regenie_step1_command <- function(config, group, pfile_prefix, extract, pheno, covar, keep, trait_list,
                                        pred_list, out_prefix, script_out, threads) {
  command <- regenie_step1_args(config, group, pfile_prefix, extract, pheno, covar, keep, trait_list, out_prefix, threads)
  done <- paste0(out_prefix, ".done")
  if (!length(command$traits)) {
    write_bash_script(script_out, c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      mkdir_parent_line(pred_list),
      paste(": >", shell_quote(pred_list)),
      printf_line("skipped_no_traits", done)
    ))
    return(invisible(TRUE))
  }

  observed <- paste0(out_prefix, "_pred.list")
  lines <- c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    mkdir_parent_line(out_prefix),
    shell_command_line(regenie_tool(config), command$args),
    paste("test -s", shell_quote(observed))
  )
  if (!identical(observed, pred_list)) {
    lines <- c(lines, paste("cp -f", shell_quote(observed), shell_quote(pred_list)))
  }
  lines <- c(lines, printf_line("ok", done))
  write_bash_script(script_out, lines)
}


regenie_step2_args <- function(config, group, pfile_prefix, pheno, covar, pred_list, trait_list, out_prefix, threads) {
  info <- group_info(config, group)
  traits <- readLines(trait_list, warn = FALSE)
  traits <- traits[nzchar(traits)]
  covars <- info$covariates
  command <- c(
    "--step", "2",
    "--pgen", pfile_prefix,
    "--phenoFile", pheno,
    "--phenoColList", paste(traits, collapse = ","),
    "--covarFile", covar,
    "--covarColList", paste(covars, collapse = ","),
    "--pred", pred_list,
    "--htp", phase2_htp_cohort_name(config),
    regenie_type_args(info$trait_type),
    "--bsize", as.character(config$phase2_regenie$step2_bsize %||% 400),
    "--minMAC", as.character(config$phase2_regenie$min_mac %||% 1),
    "--threads", threads
  )
  if (truthy(config$qc$use_mach_r2_filter %||% FALSE)) {
    command <- c(command, "--minINFO", as.character(config$qc$info_min %||% 0.8))
  }
  if (identical(info$trait_type, "bt")) {
    command <- c(command, "--firth", "--approx", "--pThresh", as.character(config$phase2_regenie$p_thresh %||% 0.01))
  } else if (truthy(config$phase2_regenie$apply_rint %||% FALSE)) {
    command <- c(command, "--apply-rint")
  }
  list(
    traits = traits,
    args = c(command, extra_regenie_options(config$phase2_regenie$step2_options %||% ""), "--out", out_prefix)
  )
}


write_regenie_step2_command <- function(config, group, pfile_prefix, pheno, covar, pred_list, trait_list,
                                        out_prefix, done, script_out, threads) {
  command <- regenie_step2_args(config, group, pfile_prefix, pheno, covar, pred_list, trait_list, out_prefix, threads)
  if (!length(command$traits)) {
    write_bash_script(script_out, c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      mkdir_parent_line(done),
      printf_line("skipped_no_traits", done)
    ))
    return(invisible(TRUE))
  }

  expected_stats <- paste0(out_prefix, "_", command$traits, ".regenie")
  lines <- c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    mkdir_parent_line(out_prefix),
    shell_command_line(regenie_tool(config), command$args),
    paste("test -s", shell_quote(expected_stats))
  )
  lines <- c(lines, printf_line("ok", done))
  write_bash_script(script_out, lines)
}


stage_trait_output <- function(config, trait, group_summary, raw_prefix, out_stats, out_summary) {
  summary <- read_tsv(group_summary)
  row <- summary[summary$trait == trait, , drop = FALSE]
  if (!nrow(row)) die("trait ", trait, " is absent from Phase 2 group summary")
  skipped <- identical(row$skipped[[1]], "True")
  raw <- paste0(raw_prefix, "_", trait, ".regenie")
  if (skipped) {
    ensure_parent(out_stats)
    writeLines(c(paste0("## skipped: ", row$skip_reason[[1]]), paste(regenie_htp_columns, collapse = "\t")), out_stats)
  } else {
    if (!file.exists(raw)) die("expected native regenie output not found: ", raw)
    ensure_parent(out_stats)
    file.copy(raw, out_stats, overwrite = TRUE)
  }
  row$native_regenie <- ifelse(skipped, "", out_stats)
  write_tsv(row, out_summary)
}


summarize_native_stats <- function(path) {
  rows <- read_tsv_no_metadata(path)
  p <- regenie_p_values(rows)
  valid <- is.finite(p) & p > 0 & p <= 1
  list(
    n_variants = nrow(rows),
    valid_p = sum(valid),
    lambda_gc = genomic_lambda(p),
    genomewide = sum(valid & p <= 5e-8),
    suggestive = sum(valid & p <= 1e-5),
    rows = rows,
    p = p
  )
}


top_hit_lines <- function(rows, p, max_rows = 10) {
  if (!nrow(rows)) return(c("No variants available."))
  keep <- is.finite(p) & p > 0 & p <= 1
  if (!any(keep)) return(c("No valid P values available."))
  rows <- rows[keep, , drop = FALSE]
  p <- p[keep]
  ord <- order(p)
  ord <- ord[seq_len(min(length(ord), max_rows))]
  chrom <- if ("CHROM" %in% names(rows)) rows$CHROM else if ("Chr" %in% names(rows)) rows$Chr else if ("chrom" %in% names(rows)) rows$chrom else "NA"
  pos <- if ("GENPOS" %in% names(rows)) rows$GENPOS else if ("Pos" %in% names(rows)) rows$Pos else if ("pos" %in% names(rows)) rows$pos else "NA"
  id <- if ("ID" %in% names(rows)) rows$ID else if ("Name" %in% names(rows)) rows$Name else if ("variant_id" %in% names(rows)) rows$variant_id else "NA"
  effect_label <- if ("Effect" %in% names(rows)) "EFFECT" else "BETA"
  beta <- if ("BETA" %in% names(rows)) rows$BETA else if ("Effect" %in% names(rows)) rows$Effect else rep("NA", nrow(rows))
  se <- if ("SE" %in% names(rows)) rows$SE else rep("NA", nrow(rows))
  c(
    paste0("| CHROM | POS | ID | ", effect_label, " | SE | P |"),
    "| --- | ---: | --- | ---: | ---: | ---: |",
    paste0("| ", chrom[ord], " | ", pos[ord], " | ", id[ord], " | ", beta[ord], " | ", se[ord], " | ", signif(p[ord], 4), " |")
  )
}


make_phase2_report <- function(config, trait, build, stats, summary_path, group_summary, union_summary, pan_summary,
                               ancestry_summary, qq, manhattan, manhattan_pdf, stage1_summaries, out) {
  summary <- read_tsv(summary_path)
  skipped <- identical(summary$skipped[[1]], "True")
  stage1 <- data.frame()
  for (path in stage1_summaries) {
    rows <- read_tsv(path)
    ancestry <- basename(dirname(path))
    stage1 <- rbind(stage1, data.frame(
      ancestry = ancestry,
      lambda_gc = metric_value(rows, "lambda_gc"),
      valid_p_value_variants = metric_value(rows, "valid_p_value_variants"),
      stringsAsFactors = FALSE
    ))
  }
  stage1_lines <- if (nrow(stage1)) {
    c(
      "| Stage 1 ancestry | Lambda GC | Valid P variants |",
      "| --- | ---: | ---: |",
      paste0("| ", stage1$ancestry, " | ", stage1$lambda_gc, " | ", stage1$valid_p_value_variants, " |")
    )
  } else {
    "No Stage 1 comparison summaries were available."
  }

  stats_summary <- list(n_variants = 0, valid_p = 0, lambda_gc = NA_real_, genomewide = 0, suggestive = 0, rows = data.frame(), p = numeric())
  if (!skipped) stats_summary <- summarize_native_stats(stats)
  lambda <- if (is.finite(stats_summary$lambda_gc)) sprintf("%.6f", stats_summary$lambda_gc) else "NA"
  top_lines <- if (skipped) c("Trait skipped before regenie.", paste0("Reason: ", summary$skip_reason[[1]])) else {
    top_hit_lines(stats_summary$rows, stats_summary$p)
  }

  union <- read_tsv(union_summary)
  pan <- read_tsv(pan_summary)
  ancestry <- read_tsv(ancestry_summary)
  unknown <- sum(ancestry$phase2_ancestry == "UNKNOWN")
  assigned <- nrow(ancestry) - unknown
  union_n <- union$variants[union$file == "UNION"][[1]]

  lines <- c(
    paste0("# Phase 2 PAN Regenie Report: ", config$project$analysis_name, " / ", trait), "",
    "## Model", "",
    paste0("- Engine: regenie"),
    paste0("- Trait type: ", summary$trait_type[[1]]),
    paste0("- Genome build: ", build),
    paste0("- Covariates: ", summary$covariates[[1]]),
    paste0("- Global PCs: ", phase2_pc_count(config)),
    paste0("- Step 1 block size: ", config$phase2_regenie$step1_bsize %||% 1000),
    paste0("- Step 2 block size: ", config$phase2_regenie$step2_bsize %||% 400),
    paste0("- Regenie HTP cohort: ", phase2_htp_cohort_name(config)),
    paste0("- Binary approximate Firth pThresh: ", config$phase2_regenie$p_thresh %||% 0.01),
    paste0("- Quantitative RINT: ", ifelse(truthy(config$phase2_regenie$apply_rint %||% FALSE), "True", "False")),
    paste0("- Native regenie output: `", stats, "`"), "",
    "## Samples", "",
    paste0("- PAN samples after sex-check and sample missingness: ", metric_value(pan, "phase2_pan_samples")),
    paste0("- POP-MaD assigned samples in PAN set: ", assigned),
    paste0("- POP-MaD UNKNOWN samples in PAN set: ", unknown),
    paste0("- Complete covariate samples: ", summary$complete_covariate_samples[[1]]),
    paste0("- Usable trait samples: ", summary$usable_n[[1]]),
    paste0("- Cases: ", summary$cases[[1]]),
    paste0("- Controls: ", summary$controls[[1]]), "",
    "## Variant QC", "",
    paste0("- Stage 1 union-pass variants: ", union_n),
    paste0("- Pooled missingness threshold: ", config$qc$geno_missing_max %||% 0.05),
    paste0("- Regenie minMAC: ", config$phase2_regenie$min_mac %||% 1),
    paste0("- Regenie minINFO: ", ifelse(truthy(config$qc$use_mach_r2_filter %||% FALSE), as.character(config$qc$info_min %||% 0.8), "not_applied")), "",
    "## Results", "",
    paste0("- Skipped: ", summary$skipped[[1]]),
    paste0("- Valid P-value variants: ", stats_summary$valid_p),
    paste0("- Lambda GC: ", lambda),
    paste0("- Genome-wide significant variants (P <= 5e-8): ", stats_summary$genomewide),
    paste0("- Suggestive variants (P <= 1e-5): ", stats_summary$suggestive),
    paste0("- QQ plot: `", qq, "`"),
    paste0("- Manhattan PNG: `", manhattan, "`"),
    paste0("- Manhattan PDF: `", manhattan_pdf, "`"), "",
    "## Stage 1 Lambda Comparison", "",
    stage1_lines, "",
    "## Top Hits", "",
    top_lines
  )
  ensure_parent(out)
  writeLines(lines, out)
}


require_args(args, "config")
config <- load_config(args$config)
threads <- args$threads %||% "1"


if (subtask == "write-groups") {
  require_args(args, "out")
  write_phase2_group_manifest(config, args$out)
} else if (subtask == "prepare-pan-genotypes") {
  require_args(args, c("sex-keep", "assignments", "excluded", "out-prefix", "keep-out", "ancestry-out", "summary-out"))
  prepare_pan_genotypes(config, args[["sex-keep"]], args$assignments, args$excluded, args[["out-prefix"]],
    args[["keep-out"]], args[["ancestry-out"]], args[["summary-out"]], threads)
} else if (subtask == "prepare-marker-set") {
  require_args(args, c("branch", "pfile-prefix", "out-prefix", "prune-prefix", "prune-in", "excluded-regions"))
  prepare_marker_set(config, args$branch, args[["pfile-prefix"]], args$keep %||% "", args[["out-prefix"]],
    args[["prune-prefix"]], args[["prune-in"]], args[["excluded-regions"]], threads)
} else if (subtask == "fit-global-pca") {
  require_args(args, c("pfile-prefix", "variants", "out-prefix"))
  fit_global_pca(config, args[["pfile-prefix"]], args$variants, args[["out-prefix"]], threads)
} else if (subtask == "score-global-pcs") {
  require_args(args, c("pfile-prefix", "variants", "weights", "frequencies", "out-prefix"))
  score_global_pcs(config, args[["pfile-prefix"]], args$variants, args$weights, args$frequencies, args[["out-prefix"]], threads)
} else if (subtask == "write-global-pcs") {
  require_args(args, c("sscore", "out"))
  write_tsv(read_sscore(args$sscore, phase2_pc_count(config)), args$out)
} else if (subtask == "build-group-inputs") {
  require_args(args, c("group", "keep", "pcs", "pheno-out", "covar-out", "summary-out", "trait-list-out", "covar-list-out", "keep-plink-out"))
  build_group_inputs(config, args$group, args$keep, args$pcs, args[["pheno-out"]], args[["covar-out"]],
    args[["summary-out"]], args[["trait-list-out"]], args[["covar-list-out"]], args[["keep-plink-out"]])
} else if (subtask == "stage1-pass-union") {
  require_args(args, c("stage1-stats", "out", "summary-out"))
  stage1_pass_union(args[["stage1-stats"]], args$out, args[["summary-out"]])
} else if (subtask == "prepare-assoc-genotypes") {
  require_args(args, c("pfile-prefix", "extract", "out-prefix"))
  prepare_assoc_genotypes(config, args[["pfile-prefix"]], args$extract, args[["out-prefix"]], threads)
} else if (subtask == "filter-step1-variants") {
  require_args(args, c("pfile-prefix", "extract", "keep", "trait-list", "out"))
  filter_step1_variants(config, args[["pfile-prefix"]], args$extract, args$keep, args[["trait-list"]], args$out,
    threads, args[["summary-out"]] %||% "", args[["excluded-out"]] %||% "")
} else if (subtask == "write-step1-command") {
  require_args(args, c("group", "pfile-prefix", "extract", "pheno", "covar", "keep", "trait-list", "pred-list", "out-prefix", "script-out"))
  write_regenie_step1_command(config, args$group, args[["pfile-prefix"]], args$extract, args$pheno, args$covar,
    args$keep, args[["trait-list"]], args[["pred-list"]], args[["out-prefix"]], args[["script-out"]], threads)
} else if (subtask == "write-step2-command") {
  require_args(args, c("group", "pfile-prefix", "pheno", "covar", "pred-list", "trait-list", "out-prefix", "done", "script-out"))
  write_regenie_step2_command(config, args$group, args[["pfile-prefix"]], args$pheno, args$covar,
    args[["pred-list"]], args[["trait-list"]], args[["out-prefix"]], args$done, args[["script-out"]], threads)
} else if (subtask == "stage-trait-output") {
  require_args(args, c("trait", "group-summary", "raw-prefix", "out-stats", "out-summary"))
  stage_trait_output(config, args$trait, args[["group-summary"]], args[["raw-prefix"]], args[["out-stats"]], args[["out-summary"]])
} else if (subtask == "make-report") {
  require_args(args, c("trait", "build", "stats", "summary", "group-summary", "union-summary", "pan-summary", "ancestry-summary", "qq", "manhattan", "manhattan-pdf", "out"))
  make_phase2_report(config, args$trait, args$build, args$stats, args$summary, args[["group-summary"]],
    args[["union-summary"]], args[["pan-summary"]], args[["ancestry-summary"]], args$qq,
    args$manhattan, args[["manhattan-pdf"]], args[["stage1-summary"]], args$out)
} else if (subtask == "check-options") {
  if (!blank(args[["options-file"]] %||% "")) {
    value <- paste(readLines(args[["options-file"]], warn = FALSE), collapse = " ")
  } else {
    require_args(args, "options")
    value <- args$options
  }
  extra_regenie_options(value)
  cat("Phase 2 regenie option pass-through is valid\n")
} else {
  die("unknown Phase 2 regenie subtask: ", subtask)
}
