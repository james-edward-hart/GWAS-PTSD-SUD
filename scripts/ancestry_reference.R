#!/usr/bin/env Rscript

# Prepare reference and study PCs for production computed ancestry.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))
source(file.path(script_dir, "lib", "variant_harmonization.R"))


# Parse the first argument as the requested ancestry-reference subtask.
raw <- commandArgs(trailingOnly = TRUE)
if (!length(raw)) die("missing ancestry-reference subtask")
subtask <- raw[[1]]
args <- parse_args(defaults = list(threads = "1", eigenvec = character()), repeated = "eigenvec", raw = raw[-1])


# Translate ancestry-reference QC settings into PLINK2 filters.
ancestry_filters <- function(config) {
  filters <- config$ancestry_reference$filters
  out <- character()
  if (truthy(filters$autosome_only %||% TRUE)) out <- c(out, "--autosome")
  if (truthy(filters$snps_only_acgt %||% TRUE)) out <- c(out, "--snps-only", "just-acgt")
  out <- c(out, "--max-alleles", as.character(filters$max_alleles %||% 2))
  if (!is.null(filters$maf_min)) out <- c(out, "--maf", as.character(filters$maf_min))
  if (!is.null(filters$geno_missing_max)) out <- c(out, "--geno", as.character(filters$geno_missing_max))
  if (truthy(filters$remove_duplicate_ids %||% TRUE)) out <- c(out, "--rm-dup", "exclude-all")
  out
}


# Convert a configured genotype block to a filtered, sorted ancestry-only PGEN.
convert_genotypes <- function(config, block, out_prefix, threads) {
  ensure_parent(paste0(out_prefix, ".pgen"))
  run_command(plink_tool(config), c(
    plink_input_args(block, out_prefix, "ancestry genotype input"),
    ancestry_filters(config),
    "--make-pgen", "--sort-vars",
    "--threads", threads,
    "--out", out_prefix
  ))
}


# Harmonize source-native IDs only inside ancestry working copies.
shared_variants <- function(config, reference_prefix, study_prefix, out,
                            mismatch_report, mapping, reference_extract,
                            study_extract, reference_update, study_update) {
  result <- harmonize_pvar_variants(
    reference_prefix = reference_prefix,
    study_prefix = study_prefix,
    shared_variants = out,
    mismatch_report = mismatch_report,
    mapping = mapping,
    reference_extract = reference_extract,
    study_extract = study_extract,
    reference_update = reference_update,
    study_update = study_update,
    exclude_palindromic = truthy(
      config$ancestry_reference$filters$exclude_palindromic %||% TRUE
    )
  )
  minimum <- as.integer(config$ancestry_reference$min_shared_variants %||% 10000)
  warn_below <- as.integer(config$ancestry_reference$warn_shared_variants_below %||% 50000)
  if (result$retained < minimum) {
    die("only ", result$retained,
      " ancestry reference/study variants survived harmonization; minimum required is ",
      minimum)
  }
  if (result$retained < warn_below) {
    warning("only ", result$retained,
      " ancestry reference/study variants survived harmonization; warning threshold is ",
      warn_below)
  }
  cat("POP-MaD study/reference overlapping variants:", result$retained, "\n")
  cat("Wrote", result$retained, "shared ancestry markers;",
    result$mismatches, "exclusion records\n")
}


# Write variants that fall inside long-range LD/problem regions.
write_region_exclusions <- function(config, pfile_prefix, out) {
  regions_path <- config$ancestry_reference$exclusion_regions %||% ""
  regions <- data.frame(
    chrom = character(), start = integer(), end = integer(),
    stringsAsFactors = FALSE
  )
  if (nzchar(regions_path) && file.exists(regions_path)) {
    regions <- read_tsv(regions_path)
    chrom_col <- if ("chrom" %in% names(regions)) "chrom" else "#chrom"
    regions$chrom <- clean_chrom(regions[[chrom_col]])
    regions <- regions[regions$chrom %in% harmonization_autosomes, , drop = FALSE]
  }
  write_pvar_region_exclusions(
    pfile_prefix, regions, out, label = "shared ancestry reference PVAR"
  )
}


# Read PLINK2 score output and normalize projected PC columns.
read_sscore <- function(path, pcs) {
  rows <- read_tsv(path)
  fid_col <- intersect(c("#FID", "FID"), names(rows))
  iid_col <- intersect(c("IID", "#IID"), names(rows))
  if (!length(iid_col)) {
    die("projected score file is missing IID/#IID sample ID column: ", path)
  }
  iid_col <- iid_col[[1]]
  fid <- if (length(fid_col)) rows[[fid_col[[1]]]] else rows[[iid_col]]
  iid <- rows[[iid_col]]
  if (length(fid) != nrow(rows) || length(iid) != nrow(rows)) {
    die("projected score file has malformed sample ID columns: ", path)
  }
  pc_cols <- grep("_AVG$", names(rows), value = TRUE)
  if (length(pc_cols) < pcs) pc_cols <- grep("^PC[0-9]+$", names(rows), value = TRUE)
  if (length(pc_cols) < pcs) die("expected at least ", pcs, " projected PC columns in ", path)
  out <- data.frame(FID = fid, IID = iid, stringsAsFactors = FALSE)
  for (i in seq_len(pcs)) out[[paste0("PC", i)]] <- rows[[pc_cols[[i]]]]
  out
}


# Join projected reference PCs to population metadata.
write_reference_pcs <- function(config, sscore, metadata_path, out) {
  if (blank(metadata_path)) metadata_path <- config$ancestry_reference$metadata$path %||% ""
  require_existing_file(metadata_path, "ancestry reference metadata")
  pcs <- as.integer(config$popmad$pcs %||% 10)
  projected <- read_sscore(sscore, pcs)
  meta <- read_tsv(metadata_path)
  settings <- config$ancestry_reference$metadata
  sample_col <- settings$sample_id_column
  fid_col <- settings$fid_column %||% ""
  population_col <- settings$population_column
  super_col <- settings$super_population_column
  fid <- if (nzchar(fid_col) && fid_col %in% names(meta)) meta[[fid_col]] else meta[[sample_col]]
  iid <- meta[[sample_col]]
  meta_key <- c(paste(fid, iid, sep = "\t"), paste(iid, iid, sep = "\t"))
  population <- c(meta[[population_col]], meta[[population_col]])
  super_population <- c(meta[[super_col]], meta[[super_col]])
  idx <- match(paste(projected$FID, projected$IID, sep = "\t"), meta_key)
  if (any(is.na(idx))) die("reference metadata missing for ", sum(is.na(idx)), " PCA samples")
  projected$population <- population[idx]
  projected$super_population <- super_population[idx]
  write_tsv(projected[c("FID", "IID", "population", "super_population", paste0("PC", seq_len(pcs)))], out)
}


# Compare fitted reference PCs against projected reference PCs.
validate_projection <- function(config, eigenvec, projected_pcs, out) {
  pcs <- as.integer(config$popmad$pcs %||% 10)
  minimum <- as.numeric(config$ancestry_reference$pca$min_projection_pc_correlation %||% 0.95)
  original <- read_tsv(eigenvec)
  projected <- read_tsv(projected_pcs)
  original_ids <- table_sample_ids(original, paste("PCA eigenvec file", eigenvec))
  projected_ids <- table_sample_ids(projected, paste("projected PC file", projected_pcs))
  pi <- match_sample_rows(original_ids, sample_key_map(projected_ids, "projected PC samples"))
  shared <- !is.na(pi)
  if (sum(shared) < 3) die("fewer than 3 reference samples are shared between PCA eigenvec and projected PC files")
  oi <- which(shared)
  pi <- pi[shared]
  rows <- lapply(seq_len(pcs), function(i) {
    pc <- paste0("PC", i)
    corr <- cor(as.numeric(original[[pc]][oi]), as.numeric(projected[[pc]][pi]))
    data.frame(pc = pc, n_samples = length(shared), correlation = corr, abs_correlation = abs(corr),
      minimum_abs_correlation = minimum, status = ifelse(abs(corr) >= minimum, "pass", "fail"))
  })
  result <- do.call(rbind, rows)
  write_tsv(result, out)
  failed <- result[result$status == "fail", ]
  if (nrow(failed)) die("reference PCA projection validation failed; low absolute correlations: ",
    paste(paste0(failed$pc, "=", sprintf("%.4f", failed$correlation)), collapse = ", "))
}


# Combine per-ancestry PCA outputs into one covariate table.
combine_within_pcs <- function(config, eigenvecs, out) {
  pcs <- as.integer(config$popmad$pcs %||% 10)
  rows <- data.frame()
  for (item in eigenvecs) {
    parts <- strsplit(item, ":", fixed = TRUE)[[1]]
    ancestry <- parts[[1]]
    path <- paste(parts[-1], collapse = ":")
    eigenvec <- read_tsv(path)
    ids <- table_sample_ids(eigenvec, paste("within-ancestry PCA eigenvec file", path))
    block <- data.frame(FID = ids$FID, IID = ids$IID, ancestry = ancestry, stringsAsFactors = FALSE)
    for (i in seq_len(pcs)) block[[paste0("PC", i)]] <- eigenvec[[paste0("PC", i)]] %||% "NA"
    rows <- rbind(rows, block)
  }
  if (!nrow(rows)) die("no within-ancestry PCA rows were available")
  write_tsv(rows, out)
}


# Write a markdown report for reference preparation QC.
write_report <- function(config, shared_variants, mapping, prune_in,
                         mismatch_report, harmonized_validation,
                         projection_validation, reference_pcs, study_pcs, out) {
  projection <- read_tsv(projection_validation)
  projection_summary <- if (nrow(projection)) {
    paste(sprintf("%s abs(r)=%.3f", projection$pc, as.numeric(projection$abs_correlation)), collapse = ", ")
  } else "NA"
  shared_count <- count_lines(shared_variants)
  warn_below <- as.integer(config$ancestry_reference$warn_shared_variants_below %||% 50000)
  minimum <- as.integer(config$ancestry_reference$min_shared_variants %||% 10000)
  overlap_status <- if (shared_count < minimum) "fail" else if (shared_count < warn_below) "warn" else "pass"
  text <- c(
    "# Ancestry Reference Preparation Report", "",
    "This report is generated by the production computed-ancestry prep workflow.", "",
    "## Inputs", "",
    paste0("- Reference genotype type: `", config$ancestry_reference$reference_genotypes$type, "`"),
    paste0("- Reference variant set: `", config$ancestry_reference$variant_set %||% "", "`"),
    paste0("- Reference metadata: `", config$ancestry_reference$metadata$path, "`"),
    paste0("- Exclusion regions: `", config$ancestry_reference$exclusion_regions %||% "", "`"),
    paste0("- Reference manifest: `", config$resources$reference_manifest, "`"), "",
    "## Counts", "",
    paste0("- Shared harmonized variants: ", shared_count),
    paste0("- Shared variant status: ", overlap_status, " (warn below ", warn_below, "; fail below ", minimum, ")"),
    paste0("- LD-pruned PCA variants: ", count_lines(prune_in)),
    paste0("- Harmonization exclusion report rows: ", max(count_lines(mismatch_report) - 1, 0)),
    paste0("- Reference PC rows: ", max(count_lines(reference_pcs) - 1, 0)),
    paste0("- Study projected PC rows: ", max(count_lines(study_pcs) - 1, 0)), "",
    "## Variant-ID Harmonization", "",
    paste0("- Retained mapping table: `", mapping, "`"),
    paste0("- Exclusion report: `", mismatch_report, "`"),
    paste0("- Renamed-copy validation: `", harmonized_validation, "`"),
    "- Working-ID scope: reference and study ancestry copies only.",
    "- Association-ID scope: source study genotypes, Stage 1 GWAS, and ordinary Phase 2 outputs retain their native IDs.",
    "- ID precedence: reference rsID, then study rsID, then CPRA from a PVAR-known REF.",
    "- PVAR-known REF means the PVAR row lacks PLINK's INFO/PR provisional-REF flag; this workflow does not perform an independent FASTA check.", "",
    "## Projection Validation", "",
    paste0("- Validation file: `", projection_validation, "`"),
    paste0("- Reference eigenvec vs projected reference PC correlations: ", projection_summary), "",
    "## ADMIXTURE", "",
    paste0("- Enabled: ", truthy(config$admixture$enabled %||% FALSE)),
    "- Status: independent report-only QC branch when enabled. POP-MaD PCA/Mahalanobis assignment remains the active ancestry-label method.",
    "- Decision: ADMIXTURE outputs are reviewed for QC only and do not alter strata, keep files, or GWAS covariates.", "",
    "## Methods", "",
    "Study and reference genotypes were converted to sorted, filtered PGEN working copies and matched",
    "by normalized chromosome, position, and unordered allele pair. Each side was extracted with its",
    "source-native IDs and then renamed to the selected ancestry-only operational ID. Long-range LD/problem-region variants were",
    "excluded before PLINK2 LD pruning. PCA was fitted in the reference samples with allele weights,",
    "then both reference and study samples were scored with the same PLINK2 `--score` projection command",
    "so POP-MaD assignment uses PCs on the same projection scale.", "",
    "Review this report before using computed ancestry labels for production GWAS."
  )
  ensure_parent(out)
  writeLines(text, out)
}


# Load config shared by every subtask.
require_args(args, "config")
config <- load_config(args$config)
threads <- args$threads %||% "1"


# Dispatch to the requested reference-preparation subtask.
if (subtask == "convert-reference") {
  # Convert the external reference dataset to filtered PGEN.
  require_args(args, "out-prefix")
  convert_genotypes(config, config$ancestry_reference$reference_genotypes, args[["out-prefix"]], threads)
} else if (subtask == "convert-study") {
  # Convert the study dataset to filtered PGEN.
  require_args(args, "out-prefix")
  convert_genotypes(config, config$genotypes, args[["out-prefix"]], threads)
} else if (subtask == "shared-variants") {
  # Match by locus/alleles and write all extraction, rename, and audit files.
  require_args(args, c(
    "reference-prefix", "study-prefix", "out", "mismatch-report", "mapping",
    "reference-extract", "study-extract", "reference-update", "study-update"
  ))
  shared_variants(
    config,
    args[["reference-prefix"]],
    args[["study-prefix"]],
    args$out,
    args[["mismatch-report"]],
    args$mapping,
    args[["reference-extract"]],
    args[["study-extract"]],
    args[["reference-update"]],
    args[["study-update"]]
  )
} else if (subtask == "extract-shared") {
  # Extract with native IDs, then rename only the ancestry-local copy.
  require_args(args, c("input-prefix", "variants", "update-name", "out-prefix"))
  extract_and_rename_pgen(
    config,
    args[["input-prefix"]],
    args$variants,
    args[["update-name"]],
    args[["out-prefix"]],
    threads
  )
} else if (subtask == "validate-harmonized") {
  # Require identical operational ID/locus/allele sets before PCA.
  require_args(args, c("reference-prefix", "study-prefix", "out"))
  count <- validate_harmonized_pgen_variants(
    args[["reference-prefix"]], args[["study-prefix"]], args$out
  )
  cat("Validated", count, "identical harmonized ancestry markers\n")
} else if (subtask == "ld-prune") {
  # Exclude problem regions and LD-prune PCA markers.
  require_args(args, c("pfile-prefix", "out-prefix", "prune-in", "excluded-regions"))
  if (identical(tolower(config$ancestry_reference$variant_set %||% ""), "pre_ld_pruned")) {
    # POP-MaD package panels are already LD-pruned upstream. Re-pruning would
    # change the intended reference geometry, so reuse the harmonized shared set.
    require_args(args, "shared-variants")
    ensure_parent(args[["prune-in"]])
    file.copy(args[["shared-variants"]], args[["prune-in"]], overwrite = TRUE)
    writeLines(character(), args[["excluded-regions"]])
    writeLines(character(), paste0(args[["out-prefix"]], ".prune.out"))
    if (!file.exists(args[["prune-in"]]) || file.info(args[["prune-in"]])$size == 0) die("pre-LD-pruned shared marker copy did not produce a non-empty prune.in file")
    cat("POP-MaD reference panel is pre-LD-pruned; copied", count_lines(args[["prune-in"]]), "shared markers to", args[["prune-in"]], "\n")
    quit(save = "no", status = 0)
  }
  n_excluded <- write_region_exclusions(config, args[["pfile-prefix"]], args[["excluded-regions"]])
  command <- c("--pfile", args[["pfile-prefix"]])
  if (n_excluded > 0) command <- c(command, "--exclude", args[["excluded-regions"]])
  pruning <- config$ancestry_reference$ld_prune
  command <- c(command, "--indep-pairwise", as.character(pruning$window %||% 1500), as.character(pruning$step %||% 150), as.character(pruning$r2 %||% 0.1),
    "--threads", threads, "--out", args[["out-prefix"]])
  run_command(plink_tool(config), command)
  if (!file.exists(args[["prune-in"]]) || file.info(args[["prune-in"]])$size == 0) die("LD pruning did not produce a non-empty prune.in file")
} else if (subtask == "fit-reference-pca") {
  # Fit reference PCA and export allele weights.
  require_args(args, c("pfile-prefix", "variants", "out-prefix"))
  pcs <- as.integer(config$popmad$pcs %||% 10)
  pca_args <- c(as.character(pcs), "allele-wts", "vcols=chrom,ref,alt")
  if (truthy(config$ancestry_reference$pca$approx %||% TRUE)) pca_args <- c(pca_args, "approx")
  run_command(plink_tool(config), c("--pfile", args[["pfile-prefix"]], "--extract", args$variants, "--freq", "counts", "--pca", pca_args, "--threads", threads, "--out", args[["out-prefix"]]))
} else if (subtask == "score-pcs") {
  # Project samples onto reference PCA weights.
  require_args(args, c("pfile-prefix", "variants", "weights", "frequencies", "out-prefix"))
  pcs <- as.integer(config$popmad$pcs %||% 10)
  # The frequency file comes from the reference PCA fit; using it for both
  # reference and study scoring keeps projected PCs on the same scale.
  run_command(plink_tool(config), c("--pfile", args[["pfile-prefix"]], "--extract", args$variants, "--read-freq", args$frequencies,
    "--score", args$weights, "2", "5", "header-read", "no-mean-imputation", "variance-standardize",
    "--score-col-nums", paste0("6-", 5 + pcs), "--threads", threads, "--out", args[["out-prefix"]]))
} else if (subtask == "write-reference-pcs") {
  # Join projected reference PCs to population labels.
  require_args(args, c("sscore", "out"))
  write_reference_pcs(config, args$sscore, args$metadata %||% "", args$out)
} else if (subtask == "write-study-pcs") {
  # Write projected study PCs in POP-MaD input format.
  require_args(args, c("sscore", "out"))
  write_tsv(read_sscore(args$sscore, as.integer(config$popmad$pcs %||% 10)), args$out)
} else if (subtask == "validate-projection") {
  # Check that projected reference PCs match fitted PCs.
  require_args(args, c("eigenvec", "projected-pcs", "out"))
  validate_projection(config, args$eigenvec, args[["projected-pcs"]], args$out)
} else if (subtask == "within-ancestry-pca") {
  # Run within-ancestry PCA for final GWAS covariates.
  require_args(args, c("pfile-prefix", "keep", "variants", "out-prefix"))
  pcs <- as.integer(config$popmad$pcs %||% 10)
  pca_args <- as.character(pcs)
  if (truthy(config$ancestry_reference$pca$approx %||% TRUE)) pca_args <- c(pca_args, "approx")
  run_command(plink_tool(config), c("--pfile", args[["pfile-prefix"]], plink_keep_args(args$keep, args[["out-prefix"]], "within-ancestry PCA keep file"), "--extract", args$variants, "--pca", pca_args,
    "--threads", threads, "--out", args[["out-prefix"]]))
} else if (subtask == "combine-within-pcs") {
  # Combine per-ancestry PCA outputs into one table.
  require_args(args, c("eigenvec", "out"))
  combine_within_pcs(config, args$eigenvec, args$out)
} else if (subtask == "write-report") {
  # Write the reference-preparation QC report.
  require_args(args, c(
    "shared-variants", "mapping", "prune-in", "mismatch-report",
    "harmonized-validation", "projection-validation",
    "reference-pcs", "study-pcs", "out"
  ))
  write_report(
    config,
    args[["shared-variants"]],
    args$mapping,
    args[["prune-in"]],
    args[["mismatch-report"]],
    args[["harmonized-validation"]],
    args[["projection-validation"]],
    args[["reference-pcs"]],
    args[["study-pcs"]],
    args$out
  )
} else {
  die("unknown ancestry-reference subtask: ", subtask)
}
