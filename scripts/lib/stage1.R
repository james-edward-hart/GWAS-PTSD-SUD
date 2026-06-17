# Shared helpers for Stage 1 R task scripts.


# Stop with a consistent pipeline error prefix.
die <- function(...) {
  stop(paste0("ERROR: ", paste0(..., collapse = "")), call. = FALSE)
}


# Create an output directory just before writing a file.
ensure_parent <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}


# Load YAML config through the R yaml package.
load_config <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    die("R package 'yaml' is required")
  }
  yaml::read_yaml(path)
}


# Read tab-delimited files while preserving literal NA strings.
read_tsv <- function(path) {
  if (!file.exists(path)) die("tab-delimited file not found: ", path)
  if (file.info(path)$size == 0) die("tab-delimited file is empty: ", path)
  tryCatch(
    read.delim(
      path,
      sep = "\t",
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


# Write a standard tab-delimited file with no row names.
write_tsv <- function(x, path) {
  ensure_parent(path)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}


# Write a PLINK ID file with no header for direct --keep/--remove use.
write_plink_id_file <- function(ids, path) {
  ensure_parent(path)
  write.table(ids[c("FID", "IID")], path, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE, na = "")
}


# Parse simple --key value CLI arguments.
parse_args <- function(defaults = list(), repeated = character(), flags = character(), raw = commandArgs(trailingOnly = TRUE)) {
  out <- defaults
  i <- 1
  while (i <= length(raw)) {
    key <- raw[[i]]
    if (!startsWith(key, "--")) die("unknown argument: ", key)
    name <- sub("^--", "", key)
    if (name %in% flags) {
      out[[name]] <- TRUE
      i <- i + 1
      next
    }
    if (i == length(raw) || startsWith(raw[[i + 1]], "--")) {
      die("missing value for ", key)
    }
    if (name %in% repeated) {
      j <- i + 1
      values <- character()
      while (j <= length(raw) && !startsWith(raw[[j]], "--")) {
        values <- c(values, raw[[j]])
        j <- j + 1
      }
      out[[name]] <- c(out[[name]], values)
      i <- j
    } else {
      value <- raw[[i + 1]]
      out[[name]] <- value
      i <- i + 2
    }
  }
  out
}


# Require named arguments after parsing.
require_args <- function(args, names) {
  missing <- names[!vapply(names, function(x) isTRUE(length(args[[x]]) > 0) && !identical(args[[x]], ""), logical(1))]
  if (length(missing)) die("missing required argument(s): ", paste(missing, collapse = ", "))
}


# Interpret common config truth values.
truthy <- function(value) {
  isTRUE(value) || identical(value, "true") || identical(value, "True") || identical(value, "1")
}


# Normalize comma-separated strings to character vectors.
split_csv <- function(value) {
  if (is.null(value) || length(value) == 0) return(character())
  if (length(value) == 1 && (is.na(value) || identical(value, ""))) return(character())
  if (is.list(value) || length(value) > 1) return(as.character(value))
  items <- trimws(strsplit(as.character(value), ",", fixed = TRUE)[[1]])
  items[nzchar(items)]
}


# Merge global and trait-specific GWAS covariates.
covariates_for_trait <- function(config, trait_id) {
  covars <- c(
    unlist(config$gwas$default_covariates, use.names = FALSE),
    unlist(config$gwas$extra_covariates %||% character(), use.names = FALSE)
  )
  traits <- read_tsv(config$inputs$trait_registry)
  row <- traits[traits$trait_id == trait_id, , drop = FALSE]
  if (nrow(row)) covars <- c(covars, split_csv(row$covariates %||% ""))
  unique(covars[nzchar(covars)])
}


# Return a fallback when a value is NULL or empty.
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}


# Convert project.analysis_name to the filename-safe prefix used for final outputs.
analysis_output_name <- function(config_or_name) {
  name <- if (is.list(config_or_name)) {
    config_or_name$project$analysis_name %||% ""
  } else {
    config_or_name
  }
  name <- trimws(as.character(name[[1]] %||% ""))
  if (!nzchar(name)) die("project.analysis_name is required")
  safe <- gsub("[^A-Za-z0-9._-]+", "_", name)
  safe <- gsub("^[._-]+|[._-]+$", "", safe)
  if (!nzchar(safe)) die("project.analysis_name must contain at least one letter or number")
  safe
}


# Return the first non-missing value from a row.
first_value <- function(row, names) {
  for (name in names) {
    if (name %in% names(row)) {
      value <- as.character(row[[name]])
      if (!identical(value, "") && !identical(value, "NA") && !identical(value, ".")) return(value)
    }
  }
  "NA"
}


# Check that a table has required columns.
require_columns <- function(df, columns, label) {
  if (!nrow(df)) die(label, " is empty")
  missing <- setdiff(columns, names(df))
  if (length(missing)) die(label, " is missing required columns: ", paste(missing, collapse = ", "))
}


# Check that FID/IID rows are unique.
require_unique_ids <- function(df, label) {
  require_columns(df, c("FID", "IID"), label)
  keys <- paste(df$FID, df$IID, sep = "\t")
  duplicates <- unique(keys[duplicated(keys)])
  if (length(duplicates)) {
    die(label, " has duplicate FID/IID rows: ", paste(head(gsub("\t", " ", duplicates), 5), collapse = ", "))
  }
}


# Read a headered pipeline keep/remove file and normalize FID/IID names.
read_id_file <- function(path, label = path) {
  rows <- read_tsv(path)
  fid_col <- if ("FID" %in% names(rows)) "FID" else "#FID"
  require_columns(rows, c(fid_col, "IID"), label)
  data.frame(FID = rows[[fid_col]], IID = rows$IID, stringsAsFactors = FALSE)
}


# Convert a pipeline ID TSV to a headerless file for PLINK2 --keep.
plink_keep_args <- function(path, out_prefix, label = "keep file") {
  ids <- read_id_file(path, label)
  if (!nrow(ids)) die(label, " is empty: ", path)
  keep_path <- paste0(out_prefix, ".plink_keep.txt")
  write_plink_id_file(ids, keep_path)
  c("--keep", keep_path)
}


# Convert sex-check config thresholds to PLINK2 --check-sex modifiers.
sex_check_threshold_args <- function(settings) {
  threshold_map <- c(
    max_female_xf = "max-female-xf",
    min_male_xf = "min-male-xf",
    max_female_yrate = "max-female-yrate",
    min_male_yrate = "min-male-yrate"
  )
  configured <- vapply(names(threshold_map), function(name) {
    value <- settings[[name]] %||% ""
    nzchar(as.character(value[[1]]))
  }, logical(1))
  if (!any(configured)) {
    thresholds <- c("max-female-xf=0.2", "min-male-xf=0.8")
    attr(thresholds, "using_defaults") <- TRUE
    return(thresholds)
  }
  if (!all(configured[c("max_female_xf", "min_male_xf")])) {
    die("custom sex_check thresholds must include both max_female_xf and min_male_xf, or leave all threshold fields blank to use defaults")
  }
  if (xor(configured[["max_female_yrate"]], configured[["min_male_yrate"]])) {
    die("sex_check Y-rate thresholds must include both max_female_yrate and min_male_yrate, or leave both blank")
  }

  thresholds <- character()
  for (name in names(threshold_map)) {
    value <- settings[[name]] %||% ""
    if (configured[[name]] && is.na(suppressWarnings(as.numeric(value)))) {
      die("sex_check.", name, " must be numeric when nonblank")
    }
    if (nzchar(value)) thresholds <- c(thresholds, paste0(threshold_map[[name]], "=", value))
  }
  thresholds
}


# Summarize which fine-scale reference populations can contribute POP-MaD models.
popmad_reference_model_coverage <- function(metadata, population_col, super_col, ancestries, min_population_n) {
  require_columns(metadata, c(population_col, super_col), "ancestry reference metadata")
  if (any(!nzchar(metadata[[population_col]])) || any(!nzchar(metadata[[super_col]]))) {
    die("ancestry reference metadata contains empty population or super_population labels")
  }
  mapping <- unique(data.frame(
    population = metadata[[population_col]],
    super_population = metadata[[super_col]],
    stringsAsFactors = FALSE
  ))
  conflicts <- unique(mapping$population[duplicated(mapping$population)])
  if (length(conflicts)) {
    die("ancestry reference metadata has conflicting population -> super_population mappings: ",
      paste(head(conflicts, 5), collapse = ", "))
  }

  counts <- table(metadata[[population_col]])
  mapping$n <- as.integer(counts[mapping$population])
  retained <- mapping[mapping$n >= min_population_n, , drop = FALSE]
  low <- mapping[mapping$n < min_population_n, , drop = FALSE]
  missing_super <- setdiff(as.character(unlist(ancestries, use.names = FALSE)), unique(retained$super_population))
  list(summary = mapping, retained = retained, low = low, missing_super = missing_super)
}


# Run an external command and stop on failure.
run_command <- function(command, args) {
  status <- system2(command, args = args)
  if (!identical(status, 0L)) die("command failed: ", paste(c(command, args), collapse = " "))
}


# Normalize chromosome labels by dropping chr prefixes.
clean_chrom <- function(value) {
  sub("^(chr|CHR)", "", as.character(value))
}


# Drop variants whose chromosome/position mapping is not unique within a PVAR.
drop_duplicate_variant_mappings <- function(rows, chrom_col, pos_col, label) {
  coords <- paste(rows[[chrom_col]], rows[[pos_col]], sep = ":")
  duplicated_coord <- coords %in% coords[duplicated(coords)]
  if (!any(duplicated_coord)) return(rows)
  duplicate_coords <- sort(unique(coords[duplicated_coord]))
  warning("excluded ", sum(duplicated_coord), " variants at ", length(duplicate_coords),
    " duplicated chromosome/position mappings from ", label, "; first: ",
    paste(head(duplicate_coords, 5), collapse = ", "))
  rows[!duplicated_coord, , drop = FALSE]
}


# Resolve the configured PLINK2 executable.
plink_tool <- function(config) {
  config$tools$plink2 %||% "plink2"
}


# Resolve the configured PLINK1 executable for sample merges not supported by
# current PLINK2 pmerge builds.
plink1_tool <- function(config) {
  config$tools$plink1 %||% config$tools$plink %||% "plink"
}


# Convert genotype config to PLINK2 input arguments.
genotype_args <- function(block) {
  kind <- tolower(block$type)
  if (kind == "pgen") return(c("--pfile", block$prefix))
  if (kind == "bed") return(c("--bfile", block$prefix))
  die("unsupported genotype type '", kind, "'. Use 'pgen' or 'bed'.")
}


# Read PLINK1 BIM variant metadata with stable column names.
read_bim_variants <- function(path, label = "BIM file") {
  rows <- read.table(
    path,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    na.strings = character()
  )
  if (ncol(rows) < 6) die(label, " must contain at least 6 columns: ", path)
  names(rows)[1:6] <- c("chrom", "variant_id", "cm", "pos", "allele1", "allele2")
  rows$row_number <- seq_len(nrow(rows))
  rows
}


# Identify BIM rows PLINK2 rejects before normal variant filters can run.
duplicate_bim_allele_rows <- function(rows) {
  allele1 <- toupper(trimws(as.character(rows$allele1)))
  allele2 <- toupper(trimws(as.character(rows$allele2)))
  rows[nzchar(allele1) & nzchar(allele2) & allele1 == allele2, , drop = FALSE]
}


# Return TRUE for regular files and symlinks, including broken symlinks.
path_exists_or_symlink <- function(path) {
  file.exists(path) || nzchar(suppressWarnings(Sys.readlink(path)))
}


# Link a large genotype component into a sanitized temporary prefix.
link_plink_component <- function(src, dest, component) {
  src_abs <- normalizePath(src, mustWork = TRUE)
  ok <- suppressWarnings(file.symlink(src_abs, dest))
  if (!isTRUE(ok) || !file.exists(dest)) {
    if (path_exists_or_symlink(dest)) unlink(dest)
    ok <- suppressWarnings(file.link(src_abs, dest))
  }
  if (!isTRUE(ok) || !file.exists(dest)) {
    die("could not create ", component, " link for sanitized PLINK input at ", dest,
      ". Check filesystem permissions for ", dirname(dest))
  }
}


# Create a lightweight PLINK1 BED view whose BIM metadata can be parsed by
# PLINK2, while excluding malformed source rows from the downstream dataset.
plink_safe_bed_args <- function(block, out_prefix, label = "genotype input") {
  prefix <- block$prefix
  bim_path <- paste0(prefix, ".bim")
  rows <- read_bim_variants(bim_path, paste(label, "BIM file"))
  invalid <- duplicate_bim_allele_rows(rows)
  if (!nrow(invalid)) return(c("--bfile", prefix))

  if (blank(out_prefix)) die("out_prefix is required to sanitize malformed BIM alleles for ", label)
  safe_prefix <- paste0(out_prefix, ".plink_safe_input")
  safe_bim <- paste0(safe_prefix, ".bim")
  safe_bed <- paste0(safe_prefix, ".bed")
  safe_fam <- paste0(safe_prefix, ".fam")
  exclude_path <- paste0(out_prefix, ".invalid_bim_alleles.exclude.txt")
  report_path <- paste0(out_prefix, ".invalid_bim_alleles.tsv")

  ensure_parent(safe_bim)
  for (path in c(safe_bed, safe_fam, safe_bim, exclude_path, report_path)) {
    if (path_exists_or_symlink(path)) unlink(path)
  }

  safe_rows <- rows
  safe_ids <- paste0("__stage1_excluded_invalid_bim_", invalid$row_number)
  invalid_index <- match(invalid$row_number, safe_rows$row_number)
  safe_rows$variant_id[invalid_index] <- safe_ids
  safe_rows$allele1[invalid_index] <- "A"
  safe_rows$allele2[invalid_index] <- "C"

  report <- data.frame(
    row_number = invalid$row_number,
    variant_id = invalid$variant_id,
    chrom = invalid$chrom,
    pos = invalid$pos,
    allele1 = invalid$allele1,
    allele2 = invalid$allele2,
    replacement_variant_id = safe_ids,
    reason = "duplicate_allele_code",
    stringsAsFactors = FALSE
  )
  write_tsv(report, report_path)
  writeLines(safe_ids, exclude_path)
  write.table(
    safe_rows[, c("chrom", "variant_id", "cm", "pos", "allele1", "allele2")],
    safe_bim,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    na = ""
  )
  link_plink_component(paste0(prefix, ".bed"), safe_bed, "BED")
  link_plink_component(paste0(prefix, ".fam"), safe_fam, "FAM")

  cat("Excluded", nrow(invalid), "BIM row(s) with duplicate allele codes before PLINK2 conversion; report:", report_path, "\n")
  c("--bfile", safe_prefix, "--exclude", exclude_path)
}


# Create a lightweight PGEN view whose PVAR metadata can be parsed by PLINK2,
# while excluding malformed source rows from the downstream dataset.
plink_safe_pgen_args <- function(block, out_prefix, label = "genotype input") {
  prefix <- block$prefix
  pvar_path <- paste0(prefix, ".pvar")
  rows <- read_tsv(pvar_path)
  chrom_col <- if ("#CHROM" %in% names(rows)) "#CHROM" else "CHROM"
  require_columns(rows, c(chrom_col, "POS", "ID", "REF", "ALT"), paste(label, "PVAR file"))
  ref <- toupper(trimws(as.character(rows$REF)))
  alt <- toupper(trimws(as.character(rows$ALT)))
  invalid <- rows[nzchar(ref) & nzchar(alt) & !grepl(",", alt, fixed = TRUE) & ref == alt, , drop = FALSE]
  if (!nrow(invalid)) return(c("--pfile", prefix))

  if (blank(out_prefix)) die("out_prefix is required to sanitize malformed PVAR alleles for ", label)
  safe_prefix <- paste0(out_prefix, ".plink_safe_input")
  safe_pgen <- paste0(safe_prefix, ".pgen")
  safe_psam <- paste0(safe_prefix, ".psam")
  safe_pvar <- paste0(safe_prefix, ".pvar")
  exclude_path <- paste0(out_prefix, ".invalid_pvar_alleles.exclude.txt")
  report_path <- paste0(out_prefix, ".invalid_pvar_alleles.tsv")

  ensure_parent(safe_pvar)
  for (path in c(safe_pgen, safe_psam, safe_pvar, exclude_path, report_path)) {
    if (path_exists_or_symlink(path)) unlink(path)
  }

  safe_rows <- rows
  invalid_rows <- match(rownames(invalid), rownames(safe_rows))
  safe_ids <- paste0("__stage1_excluded_invalid_pvar_", invalid_rows)
  safe_rows$ID[invalid_rows] <- safe_ids
  safe_rows$REF[invalid_rows] <- "A"
  safe_rows$ALT[invalid_rows] <- "C"

  report <- data.frame(
    row_number = invalid_rows,
    variant_id = invalid$ID,
    chrom = invalid[[chrom_col]],
    pos = invalid$POS,
    ref = invalid$REF,
    alt = invalid$ALT,
    replacement_variant_id = safe_ids,
    reason = "duplicate_allele_code",
    stringsAsFactors = FALSE
  )
  write_tsv(report, report_path)
  writeLines(safe_ids, exclude_path)
  write_tsv(safe_rows, safe_pvar)
  link_plink_component(paste0(prefix, ".pgen"), safe_pgen, "PGEN")
  link_plink_component(paste0(prefix, ".psam"), safe_psam, "PSAM")

  cat("Excluded", nrow(invalid), "PVAR row(s) with duplicate allele codes before PLINK2 conversion; report:", report_path, "\n")
  c("--pfile", safe_prefix, "--exclude", exclude_path)
}


# Convert genotype config to PLINK2 input arguments, sanitizing BED metadata
# rows that PLINK2 cannot parse before normal --exclude/--snps-only filters.
plink_input_args <- function(block, out_prefix = "", label = "genotype input") {
  kind <- tolower(block$type)
  if (kind == "pgen") return(plink_safe_pgen_args(block, out_prefix, label))
  if (kind == "bed") return(plink_safe_bed_args(block, out_prefix, label))
  die("unsupported genotype type '", kind, "'. Use 'pgen' or 'bed'.")
}


# Expand a configured PLINK prefix into the component files used by reports.
genotype_component_paths <- function(block, label = "genotype") {
  kind <- tolower(block$type)
  prefix <- block$prefix
  if (kind == "pgen") {
    return(setNames(
      paste0(prefix, c(".pgen", ".pvar", ".psam")),
      paste0(label, c("_pgen", "_pvar", "_psam"))
    ))
  }
  if (kind == "bed") {
    return(setNames(
      paste0(prefix, c(".bed", ".bim", ".fam")),
      paste0(label, c("_bed", "_bim", "_fam"))
    ))
  }
  die("unsupported genotype type '", kind, "'. Use 'pgen' or 'bed'.")
}


# Build shared GWAS QC filters for PLINK2.
gwas_filters <- function(config, include_hwe = TRUE) {
  qc <- config$qc
  filters <- c(
    "--maf", as.character(qc$maf_min)
  )
  if (include_hwe) filters <- c(filters, "--hwe", as.character(qc$hwe_p_min))
  filters <- c(
    filters,
    "--geno", as.character(qc$geno_missing_max),
    "--mind", as.character(qc$sample_missing_max)
  )
  if (truthy(qc$snps_only_acgt %||% TRUE)) filters <- c(filters, "--snps-only", "just-acgt")
  if (truthy(qc$autosome_only %||% TRUE)) filters <- c(filters, "--autosome")
  if (truthy(qc$use_mach_r2_filter %||% TRUE)) filters <- c(filters, "--mach-r2-filter", as.character(qc$info_min))
  filters
}


# Filters used when deriving the controls-only HWE variant list.
control_hwe_filters <- function(config) {
  qc <- config$qc
  filters <- c("--hwe", as.character(qc$hwe_p_min))
  if (truthy(qc$snps_only_acgt %||% TRUE)) filters <- c(filters, "--snps-only", "just-acgt")
  if (truthy(qc$autosome_only %||% TRUE)) filters <- c(filters, "--autosome")
  filters
}


# Count non-empty lines in a small text file.
count_lines <- function(path) {
  if (!file.exists(path)) return(0L)
  sum(nzchar(readLines(path, warn = FALSE)))
}


# Compute a file SHA-256 using whichever common SHA tool is available.
sha256_file <- function(path) {
  tool <- Sys.which("sha256sum")
  if (nzchar(tool)) {
    output <- system2(tool, path, stdout = TRUE)
  } else {
    tool <- Sys.which("shasum")
    if (!nzchar(tool)) die("neither sha256sum nor shasum is available on PATH")
    output <- system2(tool, c("-a", "256", path), stdout = TRUE)
  }
  strsplit(output[[1]], "\\s+")[[1]][[1]]
}


# Compute a SHA-256 for generated text.
sha256_text <- function(text) {
  path <- tempfile("stage1-sha256-")
  on.exit(unlink(path), add = TRUE)
  writeLines(text, path, useBytes = TRUE)
  sha256_file(path)
}


# Return TRUE when a scalar config string is missing or empty.
blank <- function(value) {
  is.null(value) || length(value) == 0 || is.na(value[[1]]) || !nzchar(as.character(value[[1]]))
}


# Fail if a required path is not an existing file.
require_existing_file <- function(path, label) {
  if (blank(path) || !file.exists(path)) die(label, " not found: ", path)
}


# Resolve package-internal paths while keeping absolute site paths usable.
is_absolute_path <- function(path) {
  grepl("^/", path) | grepl("^[A-Za-z]:[\\\\/]", path)
}


package_path <- function(root, path) {
  if (blank(path)) return("")
  path <- as.character(path[[1]])
  if (is_absolute_path(path)) return(path)
  file.path(root, path)
}


# Fail if package panel paths escape the fingerprinted package root.
assert_package_relative_paths <- function(paths, label) {
  paths <- as.character(paths)
  paths <- paths[nzchar(paths)]
  bad <- paths[is_absolute_path(paths) | grepl("(^|/)\\.\\.(/|$)", paths)]
  if (length(bad)) die(label, " paths must be package-relative and cannot contain '..': ",
    paste(head(bad, 10), collapse = ", "))
  invisible(TRUE)
}


# Fail when raw Hail/VCF/BCF artifacts are present in a distributable package.
assert_no_raw_reference_artifacts <- function(root) {
  entries <- list.files(root, recursive = TRUE, all.files = TRUE, include.dirs = TRUE, full.names = FALSE)
  if (!length(entries)) return(invisible(TRUE))
  normalized <- gsub("\\\\", "/", entries)
  raw_hail <- grepl("(^|/)[^/]+\\.(mt|ht|vds)(/|$)", normalized, ignore.case = TRUE)
  raw_variant <- grepl("\\.(vcf|vcf\\.gz|vcf\\.bgz|bcf)(\\.tbi|\\.csi)?$", normalized, ignore.case = TRUE)
  bad <- sort(unique(normalized[raw_hail | raw_variant]))
  if (length(bad)) {
    die("reference package contains raw Hail/VCF/BCF artifacts; keep builder intermediates outside the package: ",
      paste(head(bad, 10), collapse = ", "))
  }
  invisible(TRUE)
}


# List all regular files in a package, including hidden sidecar files.
reference_package_files <- function(root) {
  entries <- list.files(root, recursive = TRUE, all.files = TRUE, include.dirs = FALSE, full.names = FALSE)
  if (!length(entries)) return(character())
  normalized <- gsub("\\\\", "/", entries)
  full <- file.path(root, normalized)
  sort(normalized[file.exists(full) & !dir.exists(full)])
}


# Identify filesystem/archive metadata that macOS can add during copy or
# extraction. These files are not reference-package content.
is_reference_package_sidecar <- function(path) {
  base <- basename(path)
  base %in% c(".DS_Store") |
    grepl("(^|/)__MACOSX(/|$)", path) |
    grepl("(^|/)\\._[^/]+$", path)
}


# Read one non-empty line from a small text file.
read_first_line <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines[nzchar(trimws(lines))])
  if (!length(lines)) return("")
  strsplit(lines[[1]], "\\s+")[[1]][[1]]
}


# Compute and verify the unpacked reference-package content fingerprint.
reference_package_fingerprint <- function(root, verify_hashes = TRUE) {
  if (blank(root)) die("reference_package.root is empty")
  assert_no_raw_reference_artifacts(root)
  file_manifest_path <- file.path(root, "file_manifest.tsv")
  fingerprint_path <- file.path(root, "content_fingerprint.sha256")
  require_existing_file(file_manifest_path, "reference package file manifest")
  require_existing_file(fingerprint_path, "reference package content fingerprint")

  manifest <- read_tsv(file_manifest_path)
  require_columns(manifest, c("path", "sha256", "size_bytes", "role"), "reference package file manifest")
  if (any(!nzchar(manifest$path))) die("reference package file manifest contains an empty path")
  if (any(is_absolute_path(manifest$path) | grepl("(^|/)\\.\\.(/|$)", manifest$path))) {
    die("reference package file manifest paths must be package-relative and cannot contain '..'")
  }
  duplicate_paths <- unique(manifest$path[duplicated(manifest$path)])
  if (length(duplicate_paths)) {
    die("reference package file manifest contains duplicate paths: ", paste(head(duplicate_paths, 10), collapse = ", "))
  }

  actual_files <- reference_package_files(root)
  actual_files <- actual_files[!is_reference_package_sidecar(actual_files)]
  allowed_unmanifested <- c("file_manifest.tsv", "content_fingerprint.sha256")
  extra_files <- setdiff(actual_files, c(manifest$path, allowed_unmanifested))
  if (length(extra_files)) {
    die("reference package contains files absent from file_manifest.tsv: ", paste(head(extra_files, 10), collapse = ", "))
  }

  normalized <- manifest[order(manifest$path), , drop = FALSE]
  lines <- character(nrow(normalized))
  for (i in seq_len(nrow(normalized))) {
    rel <- normalized$path[[i]]
    full <- file.path(root, rel)
    require_existing_file(full, paste0("reference package file '", rel, "'"))
    size <- as.character(file.info(full)$size)
    expected_size <- as.character(normalized$size_bytes[[i]])
    if (nzchar(expected_size) && expected_size != size) {
      die("reference package file size mismatch for ", rel, ": expected ", expected_size, ", observed ", size)
    }
    expected_sha <- tolower(as.character(normalized$sha256[[i]]))
    if (!nzchar(expected_sha)) die("reference package file manifest missing sha256 for ", rel)
    observed_sha <- if (verify_hashes) sha256_file(full) else expected_sha
    if (!identical(observed_sha, expected_sha)) {
      die("reference package file hash mismatch for ", rel, ": expected ", expected_sha, ", observed ", observed_sha)
    }
    lines[[i]] <- paste(rel, expected_sha, size, sep = "\t")
  }

  observed <- sha256_text(paste(lines, collapse = "\n"))
  recorded <- read_first_line(fingerprint_path)
  if (!blank(recorded) && !identical(tolower(recorded), observed)) {
    die("reference package content fingerprint file does not match unpacked contents: expected ",
      tolower(recorded), ", observed ", observed)
  }
  observed
}


# Return a single build-matched panel row from a package manifest.
select_reference_panel <- function(panel_manifest, role, genome_build) {
  require_columns(panel_manifest, c("panel_id", "role", "build", "genotype_type", "genotype_prefix",
    "metadata_path", "exclusion_regions"), "reference package panel manifest")
  rows <- panel_manifest[
    tolower(panel_manifest$role) == tolower(role) &
      toupper(panel_manifest$build) == toupper(genome_build),
    , drop = FALSE
  ]
  if (!nrow(rows)) die("reference package is missing a ", role, " panel for inferred build ", genome_build)
  if (nrow(rows) > 1) die("reference package has multiple ", role, " panels for inferred build ", genome_build)
  rows[1, , drop = FALSE]
}


# Return the manifest-covered files required by one panel row.
panel_required_paths <- function(row) {
  prefix <- panel_value(row, "genotype_prefix")
  kind <- tolower(panel_value(row, "genotype_type"))
  genotype_files <- switch(kind,
    pgen = paste0(prefix, c(".pgen", ".pvar", ".psam")),
    bed = paste0(prefix, c(".bed", ".bim", ".fam")),
    die("unsupported genotype_type in reference package panel manifest for ",
      panel_value(row, "panel_id", "<unknown>"), ": ", kind)
  )
  c(genotype_files, panel_value(row, "metadata_path"), panel_value(row, "exclusion_regions"))
}


# Ensure panel manifests only reference files covered by file_manifest.tsv.
validate_reference_panel_paths <- function(panel_manifest, file_manifest) {
  require_columns(panel_manifest, c("panel_id", "role", "build", "genotype_type", "genotype_prefix",
    "metadata_path", "exclusion_regions"), "reference package panel manifest")
  covered <- as.character(file_manifest$path)
  missing <- character()
  for (i in seq_len(nrow(panel_manifest))) {
    row <- panel_manifest[i, , drop = FALSE]
    required <- panel_required_paths(row)
    assert_package_relative_paths(required, paste0("reference package panel ", panel_value(row, "panel_id", "<unknown>")))
    missing <- c(missing, setdiff(required, covered))
  }
  missing <- unique(missing[nzchar(missing)])
  if (length(missing)) {
    die("reference package panel manifest references files absent from file_manifest.tsv: ",
      paste(head(missing, 10), collapse = ", "))
  }
  invisible(TRUE)
}


# Reject Hail array/stringified IDs that break PLINK ID-based harmonization.
validate_reference_panel_variant_ids <- function(root, panel_manifest) {
  for (i in seq_len(nrow(panel_manifest))) {
    row <- panel_manifest[i, , drop = FALSE]
    if (tolower(panel_value(row, "genotype_type")) != "pgen") next
    pvar <- package_path(root, paste0(panel_value(row, "genotype_prefix"), ".pvar"))
    variants <- read_tsv(pvar)
    require_columns(variants, "ID", paste0("PVAR for package panel ", panel_value(row, "panel_id", "<unknown>")))
    bad <- variants$ID[grepl("^\\[", variants$ID)]
    if (length(bad)) {
      die("PVAR for package panel ", panel_value(row, "panel_id", "<unknown>"),
        " contains Hail array-style variant IDs; rebuild with scalar rsIDs. First: ",
        paste(head(bad, 5), collapse = ", "))
    }
  }
  invisible(TRUE)
}


# Pull an optional scalar from a one-row manifest table.
panel_value <- function(row, name, default = "") {
  if (!name %in% names(row)) return(default)
  value <- as.character(row[[name]][[1]])
  if (is.na(value) || !nzchar(value)) default else value
}


# Install one resolved package panel into the active Stage 1 config tree.
apply_reference_panel <- function(config, section, row, root, genome_build) {
  if (is.null(config[[section]])) config[[section]] <- list()

  config[[section]]$reference_genotypes <- list(
    type = panel_value(row, "genotype_type"),
    prefix = package_path(root, panel_value(row, "genotype_prefix"))
  )

  existing_metadata <- config[[section]]$metadata %||% list()
  config[[section]]$metadata <- list(
    path = package_path(root, panel_value(row, "metadata_path")),
    sample_id_column = panel_value(row, "sample_id_column", existing_metadata$sample_id_column %||% "sample_id"),
    fid_column = panel_value(row, "fid_column", existing_metadata$fid_column %||% ""),
    population_column = panel_value(row, "population_column", existing_metadata$population_column %||% "population"),
    super_population_column = panel_value(row, "super_population_column", existing_metadata$super_population_column %||% "super_population")
  )

  config[[section]]$exclusion_regions <- package_path(root, panel_value(row, "exclusion_regions"))
  config[[section]]$reference_genome_build <- genome_build
  config[[section]]$source_panel_id <- panel_value(row, "panel_id")
  config[[section]]$variant_set <- panel_value(row, "variant_set", config[[section]]$variant_set %||% "")
  config[[section]]$source_uri <- panel_value(row, "source_uri", config[[section]]$source_uri %||% "")
  config[[section]]$source_sample_set <- panel_value(row, "source_sample_set", config[[section]]$source_sample_set %||% "")
  config[[section]]$source_variant_count <- panel_value(row, "source_variant_count", config[[section]]$source_variant_count %||% "")

  if (section == "admixture") {
    labels <- split_csv(panel_value(row, "labels", ""))
    if (length(labels)) {
      configured <- as.character(unlist(config$admixture$labels %||% character(), use.names = FALSE))
      if (length(configured) && !setequal(configured, labels)) {
        die("admixture.labels do not match package panel labels for ", panel_value(row, "panel_id"),
          ": config=", paste(configured, collapse = ","), "; package=", paste(labels, collapse = ","))
      }
      if (!length(configured)) config$admixture$labels <- labels
      if (is.null(config$admixture$k)) config$admixture$k <- length(labels)
    }
  }

  config
}


# Validate and resolve a custom reference package against the inferred study build.
resolve_reference_package_config <- function(config, genome_build, verify_hashes = TRUE) {
  package <- config$reference_package %||% list()
  root <- package$root %||% ""
  if (blank(root)) {
    die("reference_package.root is required")
  }

  if (!dir.exists(root)) die("reference_package.root does not exist: ", root)
  observed <- reference_package_fingerprint(root, verify_hashes = verify_hashes)
  expected <- package$fingerprint %||% ""
  if (blank(expected)) {
    die("reference_package.fingerprint is required")
  } else if (!identical(tolower(as.character(expected)), observed)) {
    die("reference package fingerprint mismatch: expected ", tolower(as.character(expected)), ", observed ", observed)
  }

  panel_manifest_path <- file.path(root, "panel_manifest.tsv")
  require_existing_file(panel_manifest_path, "reference package panel manifest")
  panels <- read_tsv(panel_manifest_path)
  file_manifest <- read_tsv(file.path(root, "file_manifest.tsv"))
  validate_reference_panel_paths(panels, file_manifest)
  validate_reference_panel_variant_ids(root, panels)

  if (truthy(config$ancestry_reference$enabled %||% FALSE)) {
    config <- apply_reference_panel(config, "ancestry_reference",
      select_reference_panel(panels, "popmad", genome_build), root, genome_build)
  }
  if (truthy(config$admixture$enabled %||% FALSE)) {
    config <- apply_reference_panel(config, "admixture",
      select_reference_panel(panels, "admixture", genome_build), root, genome_build)
  }

  config$reference_package <- as.list(package)
  config$reference_package$root <- root
  config$reference_package$observed_fingerprint <- observed
  config$reference_package$panel_manifest <- panel_manifest_path
  config$reference_package$file_manifest <- file.path(root, "file_manifest.tsv")
  config
}
