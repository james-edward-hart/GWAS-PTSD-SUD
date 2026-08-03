# Shared helpers for Stage 1 R task scripts.

# Resolve sibling shell helpers from this file, independent of the caller's
# current working directory.
stage1_library_dir <- local({
  source_path <- tryCatch(sys.frame(1)$ofile, error = function(err) "")
  if (is.null(source_path) || !nzchar(source_path)) {
    return(normalizePath(file.path("scripts", "lib"), mustWork = FALSE))
  }
  dirname(normalizePath(source_path, mustWork = FALSE))
})


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


# Return common PLINK ID aliases for matching outputs that may omit FID.
sample_key_variants <- function(fid, iid) {
  unique(paste(c(fid, iid, "0"), iid, sep = "\t"))
}


# Build an unambiguous lookup from FID/IID aliases to row numbers.
sample_key_map <- function(ids, label) {
  missing <- setdiff(c("FID", "IID"), names(ids))
  if (length(missing)) die(label, " is missing required columns: ", paste(missing, collapse = ", "))
  if (!nrow(ids)) return(data.frame(key = character(), row = integer(), stringsAsFactors = FALSE))
  variants <- mapply(sample_key_variants, ids$FID, ids$IID, SIMPLIFY = FALSE)
  out <- data.frame(
    key = unlist(variants, use.names = FALSE),
    row = rep(seq_len(nrow(ids)), lengths(variants)),
    stringsAsFactors = FALSE
  )
  conflict <- names(which(tapply(out$row, out$key, function(x) length(unique(x)) > 1)))
  if (length(conflict)) {
    die(label, " has ambiguous sample IDs under FID/IID alias matching: ",
      paste(head(gsub("\t", " ", conflict), 5), collapse = ", "))
  }
  out[!duplicated(out$key), , drop = FALSE]
}


# Match query FID/IID rows against a sample_key_map.
match_sample_row <- function(fid, iid, key_map) {
  idx <- match(sample_key_variants(fid, iid), key_map$key)
  idx <- idx[!is.na(idx)]
  if (!length(idx)) return(NA_integer_)
  key_map$row[[idx[[1]]]]
}


match_sample_rows <- function(ids, key_map) {
  missing <- setdiff(c("FID", "IID"), names(ids))
  if (length(missing)) die("sample ID table is missing required columns: ", paste(missing, collapse = ", "))
  vapply(seq_len(nrow(ids)), function(i) match_sample_row(ids$FID[[i]], ids$IID[[i]], key_map), integer(1))
}


# Normalize a table with PLINK-style ID columns to FID/IID.
table_sample_ids <- function(rows, label, missing_fid = c("iid", "zero")) {
  missing_fid <- match.arg(missing_fid)
  iid_col <- if ("IID" %in% names(rows)) "IID" else if ("#IID" %in% names(rows)) "#IID" else ""
  if (!nzchar(iid_col)) die(label, " is missing IID/#IID sample ID column")
  fid_col <- if ("#FID" %in% names(rows)) "#FID" else if ("FID" %in% names(rows)) "FID" else ""
  fid <- if (nzchar(fid_col)) {
    rows[[fid_col]]
  } else if (identical(missing_fid, "zero")) {
    rep("0", nrow(rows))
  } else {
    rows[[iid_col]]
  }
  data.frame(FID = fid, IID = rows[[iid_col]], stringsAsFactors = FALSE)
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


# Read a pipeline or PLINK keep/remove file and normalize FID/IID names.
read_id_file <- function(path, label = path) {
  if (!file.exists(path)) die("ID file not found: ", path)
  if (file.info(path)$size == 0) die("ID file is empty: ", path)
  first <- readLines(path, n = 1, warn = FALSE)
  tokens <- strsplit(trimws(first), "\\s+")[[1]]
  has_header <- any(tokens %in% c("FID", "#FID", "IID", "#IID"))
  rows <- tryCatch(
    read.table(
      path,
      header = has_header,
      sep = "",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = "",
      na.strings = character()
    ),
    error = function(err) die("could not read ID file ", path, ": ", conditionMessage(err))
  )
  if (has_header) {
    return(table_sample_ids(rows, label))
  }
  if (ncol(rows) == 1) {
    return(data.frame(FID = rows[[1]], IID = rows[[1]], stringsAsFactors = FALSE))
  }
  if (ncol(rows) < 2) die(label, " must contain at least one ID column")
  data.frame(FID = rows[[1]], IID = rows[[2]], stringsAsFactors = FALSE)
}


# Rewrite matched IDs to the exact FID/IID values used by a target genotype set.
canonicalize_sample_ids <- function(ids, reference_ids, label, reference_label = "target genotype samples") {
  if (is.null(reference_ids)) return(ids)
  require_columns(reference_ids, c("FID", "IID"), reference_label)
  if (!nrow(ids)) return(ids)
  idx <- match_sample_rows(ids, sample_key_map(reference_ids, reference_label))
  missing <- is.na(idx)
  if (any(missing)) {
    missing_labels <- paste(ids$FID[missing], ids$IID[missing])
    die(label, " contains samples absent from ", reference_label, ": ",
      paste(head(missing_labels, 5), collapse = ", "))
  }
  reference_ids[idx, c("FID", "IID"), drop = FALSE]
}


# Convert a pipeline ID TSV to a headerless file for PLINK2 --keep.
plink_keep_args <- function(path, out_prefix, label = "keep file", reference_ids = NULL,
                            reference_label = "target genotype samples") {
  ids <- read_id_file(path, label)
  if (!nrow(ids)) die(label, " is empty: ", path)
  ids <- canonicalize_sample_ids(ids, reference_ids, label, reference_label)
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


# Calculate genomic inflation from valid P values.
genomic_lambda <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) return(NA_real_)
  chi <- suppressWarnings(qchisq(p, df = 1, lower.tail = FALSE))
  chi <- chi[is.finite(chi)]
  if (!length(chi)) return(NA_real_)
  median(chi, na.rm = TRUE) / qchisq(0.5, df = 1, lower.tail = FALSE)
}


# Make result links relative to Markdown reports when both use standard paths.
report_relative_path <- function(path, out_path) {
  path <- gsub("\\\\", "/", path)
  out_path <- gsub("\\\\", "/", out_path)
  if (startsWith(path, "results/") && startsWith(out_path, "results/reports/")) {
    return(file.path("..", "..", sub("^results/", "", path)))
  }
  path
}


# Normalize chromosome labels by dropping case-insensitive chr prefixes.
clean_chrom <- function(value) {
  sub("^chr", "", as.character(value), ignore.case = TRUE)
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


# Cache completed scans only within the current R task. This lets a caller use
# chromosome counts after sanitization without rereading a large BIM/PVAR.
plink_metadata_cache <- new.env(parent = emptyenv())


plink_metadata_path <- function(block, label = "genotype input") {
  kind <- tolower(block$type)
  suffix <- if (kind == "bed") ".bim" else if (kind == "pgen") ".pvar" else {
    die("unsupported genotype type '", kind, "'. Use 'pgen' or 'bed'.")
  }
  path <- paste0(block$prefix, suffix)
  if (!file.exists(path)) die(label, " ", toupper(sub("^\\.", "", suffix)), " file not found: ", path)
  if (file.info(path)$size == 0) die(label, " ", toupper(sub("^\\.", "", suffix)), " file is empty: ", path)
  list(kind = kind, path = path)
}


plink_metadata_cache_key <- function(kind, path) {
  info <- file.info(path)
  paste(kind, normalizePath(path, mustWork = TRUE), info$size, as.numeric(info$mtime), sep = "\t")
}


# Parse the scanner's compact metric/value output into a named character list.
read_plink_metadata_summary <- function(path) {
  rows <- read_tsv(path)
  require_columns(rows, c("metric", "value"), "PLINK metadata scan summary")
  if (anyDuplicated(rows$metric)) die("PLINK metadata scan summary contains duplicate metrics: ", path)
  setNames(as.list(as.character(rows$value)), rows$metric)
}


# Run the shared GNU AWK scanner. Optional sanitizer paths cause a second pass
# only when duplicate allele codes are actually present.
inspect_plink_metadata <- function(block, label = "genotype input", targets = character(),
                                   stop_on_target = FALSE, inspect_alleles = FALSE,
                                   sanitizer_paths = NULL) {
  metadata <- plink_metadata_path(block, label)
  helper <- file.path(dirname(stage1_library_dir), "inspect_plink_metadata.sh")
  if (!file.exists(helper)) die("PLINK metadata helper not found: ", helper)

  summary_path <- tempfile("plink-metadata-summary-", fileext = ".tsv")
  on.exit(unlink(summary_path), add = TRUE)
  command_args <- c(
    helper,
    "--type", metadata$kind,
    "--metadata", metadata$path,
    "--summary", summary_path,
    "--label", label
  )
  if (length(targets)) {
    command_args <- c(command_args, "--target-chromosomes", paste(targets, collapse = ","))
  }
  if (isTRUE(stop_on_target)) command_args <- c(command_args, "--stop-on-target")
  if (isTRUE(inspect_alleles)) command_args <- c(command_args, "--inspect-alleles")
  if (!is.null(sanitizer_paths)) {
    command_args <- c(
      command_args,
      "--safe-metadata", sanitizer_paths$safe_metadata,
      "--invalid-report", sanitizer_paths$report,
      "--invalid-exclude", sanitizer_paths$exclude
    )
  }

  status <- suppressWarnings(system2(
    "bash",
    args = vapply(command_args, shQuote, character(1), type = "sh")
  ))
  if (!identical(status, 0L)) die("PLINK metadata inspection failed for ", metadata$path)

  summary <- read_plink_metadata_summary(summary_path)
  if (identical(summary$scan_complete, "true")) {
    key <- plink_metadata_cache_key(metadata$kind, metadata$path)
    assign(key, list(summary = summary, inspected_alleles = isTRUE(inspect_alleles)), envir = plink_metadata_cache)
  }
  summary
}


cached_plink_metadata_summary <- function(block, require_alleles = FALSE, label = "genotype input") {
  metadata <- plink_metadata_path(block, label)
  key <- plink_metadata_cache_key(metadata$kind, metadata$path)
  if (!exists(key, envir = plink_metadata_cache, inherits = FALSE)) return(NULL)
  cached <- get(key, envir = plink_metadata_cache, inherits = FALSE)
  if (isTRUE(require_alleles) && !isTRUE(cached$inspected_alleles)) return(NULL)
  cached$summary
}


# Return exact full-file counts, reusing a completed scan when one is available.
plink_metadata_summary <- function(block, label = "genotype input", require_alleles = FALSE) {
  cached <- cached_plink_metadata_summary(block, require_alleles, label)
  if (!is.null(cached)) return(cached)
  inspect_plink_metadata(block, label = label, inspect_alleles = require_alleles)
}


plink_metadata_count <- function(summary, metric) {
  value <- suppressWarnings(as.numeric(summary[[metric]]))
  if (length(value) != 1 || is.na(value)) {
    die("PLINK metadata scan summary is missing numeric metric: ", metric)
  }
  value
}


# Stream BIM/PVAR metadata in bounded chunks for ancestry validators that still
# need row-level R callbacks. Large input inspection/sanitization uses GNU AWK.
stream_plink_variant_chunks <- function(block, visit, label = "genotype input", output = "",
                                        chunk_size = 10000L, require_alleles = FALSE) {
  metadata <- plink_metadata_path(block, label)
  kind <- metadata$kind
  path <- metadata$path

  input_connection <- file(path, open = "rt")
  on.exit(close(input_connection), add = TRUE)

  output_connection <- NULL
  if (nzchar(output)) {
    ensure_parent(output)
    output_connection <- file(output, open = "wt")
    on.exit(close(output_connection), add = TRUE)
  }

  columns <- if (kind == "bed") {
    c(chrom = 1L, variant_id = 2L, pos = 4L, allele1 = 5L, allele2 = 6L, info = NA_integer_)
  } else {
    NULL
  }
  row_number <- 0L
  line_number <- 0L
  stopped_early <- FALSE

  repeat {
    lines <- readLines(input_connection, n = chunk_size, warn = FALSE)
    if (!length(lines)) break

    chunk_line_numbers <- line_number + seq_along(lines)
    line_number <- line_number + length(lines)
    capacity <- length(lines)
    data_line_index <- integer(capacity)
    data_line_number <- integer(capacity)
    data_row_number <- integer(capacity)
    chrom <- character(capacity)
    pos <- character(capacity)
    variant_id <- character(capacity)
    allele1 <- character(capacity)
    allele2 <- character(capacity)
    info <- character(capacity)
    fields <- vector("list", capacity)
    count <- 0L

    for (i in seq_along(lines)) {
      line <- lines[[i]]
      if (!nzchar(trimws(line))) next

      if (kind == "pgen" && startsWith(line, "##")) next
      parts <- if (kind == "pgen") {
        strsplit(line, "\t", fixed = TRUE)[[1]]
      } else {
        strsplit(trimws(line), "[[:space:]]+")[[1]]
      }

      if (kind == "pgen" && is.null(columns)) {
        chrom_index <- match("#CHROM", parts)
        if (is.na(chrom_index)) chrom_index <- match("CHROM", parts)
        columns <- c(
          chrom = chrom_index,
          variant_id = match("ID", parts),
          pos = match("POS", parts),
          allele1 = match("REF", parts),
          allele2 = match("ALT", parts),
          info = match("INFO", parts)
        )
        required <- c("chrom", "variant_id", "pos")
        if (require_alleles) required <- c(required, "allele1", "allele2")
        missing <- required[is.na(columns[required])]
        if (length(missing)) {
          die(label, " PVAR header is missing required column(s): ", paste(missing, collapse = ", "), ": ", path)
        }
        next
      }

      if (kind == "bed" && length(parts) < 6) {
        die(label, " BIM row ", chunk_line_numbers[[i]], " must contain at least 6 columns: ", path)
      }
      required_index <- max(columns[!is.na(columns)])
      if (length(parts) < required_index) {
        die(label, " ", toupper(kind), " metadata row ", chunk_line_numbers[[i]],
          " has fewer columns than its header: ", path)
      }

      count <- count + 1L
      row_number <- row_number + 1L
      data_line_index[[count]] <- i
      data_line_number[[count]] <- chunk_line_numbers[[i]]
      data_row_number[[count]] <- row_number
      chrom[[count]] <- parts[[columns[["chrom"]]]]
      pos[[count]] <- parts[[columns[["pos"]]]]
      variant_id[[count]] <- parts[[columns[["variant_id"]]]]
      allele1[[count]] <- if (is.na(columns[["allele1"]])) "" else parts[[columns[["allele1"]]]]
      allele2[[count]] <- if (is.na(columns[["allele2"]])) "" else parts[[columns[["allele2"]]]]
      info[[count]] <- if (is.na(columns[["info"]])) "" else parts[[columns[["info"]]]]
      fields[[count]] <- parts
    }

    result <- NULL
    if (count > 0) {
      keep <- seq_len(count)
      variants <- data.frame(
        row_number = data_row_number[keep],
        line_number = data_line_number[keep],
        chrom = chrom[keep],
        pos = suppressWarnings(as.integer(pos[keep])),
        variant_id = variant_id[keep],
        allele1 = allele1[keep],
        allele2 = allele2[keep],
        info = info[keep],
        ref_provisional = grepl("(^|;)PR(;|$)", info[keep]),
        stringsAsFactors = FALSE
      )
      result <- visit(list(
        variants = variants,
        lines = lines,
        data_line_index = data_line_index[keep],
        fields = fields[keep],
        columns = columns,
        kind = kind
      ))
      if (!is.null(result$lines)) {
        if (length(result$lines) != length(lines)) die("variant metadata callback returned the wrong number of lines")
        lines <- result$lines
      }
    }

    if (!is.null(output_connection)) writeLines(lines, output_connection, useBytes = TRUE)
    if (!is.null(result$continue) && identical(result$continue, FALSE)) {
      stopped_early <- TRUE
      break
    }
  }

  if (kind == "pgen" && is.null(columns)) die(label, " PVAR header not found: ", path)
  list(path = path, rows_scanned = row_number, stopped_early = stopped_early)
}


# Visit only normalized variant rows when raw metadata lines are not needed.
scan_plink_variant_metadata <- function(block, visit, label = "genotype input", chunk_size = 10000L,
                                        require_alleles = FALSE) {
  stream_plink_variant_chunks(
    block,
    function(chunk) {
      keep_going <- visit(chunk$variants)
      list(continue = !identical(keep_going, FALSE))
    },
    label = label,
    chunk_size = chunk_size,
    require_alleles = require_alleles
  )
}


# Normalize PLINK chromosome aliases used by BIM/PVAR and --chr.
normalize_plink_chromosome <- function(value) {
  value <- toupper(trimws(clean_chrom(value)))
  value[value == "23"] <- "X"
  value[value == "24"] <- "Y"
  value[value == "25"] <- "XY"
  value
}


# Stop the GNU AWK scan as soon as any requested chromosome is observed.
genotype_has_chromosomes <- function(block, chromosomes, label = "genotype input") {
  targets <- unique(normalize_plink_chromosome(chromosomes))
  cached <- cached_plink_metadata_summary(block, label = label)
  fixed_metrics <- c(
    X = "x_variants",
    Y = "y_variants",
    XY = "xy_variants",
    PAR1 = "par1_variants",
    PAR2 = "par2_variants"
  )
  if (!is.null(cached) && all(targets %in% names(fixed_metrics))) {
    return(any(vapply(fixed_metrics[targets], function(metric) {
      plink_metadata_count(cached, metric) > 0
    }, logical(1))))
  }

  summary <- inspect_plink_metadata(
    block,
    label = label,
    targets = targets,
    stop_on_target = TRUE
  )
  plink_metadata_count(summary, "target_matches") > 0
}


# Return TRUE for regular files and symlinks, including broken symlinks.
path_exists_or_symlink <- function(path) {
  file.exists(path) || nzchar(suppressWarnings(Sys.readlink(path)))
}


# Link a large genotype component into a sanitized temporary prefix.
link_plink_component <- function(src, dest, component) {
  src_abs <- normalizePath(src, mustWork = TRUE)
  ensure_parent(dest)
  temp <- tempfile(paste0(".", basename(dest), "-"), tmpdir = dirname(dest))
  on.exit(if (path_exists_or_symlink(temp)) unlink(temp), add = TRUE)
  ok <- suppressWarnings(file.symlink(src_abs, temp))
  if (!isTRUE(ok) || !file.exists(temp)) {
    if (path_exists_or_symlink(temp)) unlink(temp)
    ok <- suppressWarnings(file.link(src_abs, temp))
  }
  if (!isTRUE(ok) || !file.exists(temp)) {
    die("could not create ", component, " link for sanitized PLINK input at ", dest,
      ". Check filesystem permissions for ", dirname(dest))
  }
  if (path_exists_or_symlink(dest)) unlink(dest)
  if (!file.rename(temp, dest)) {
    die("could not finalize ", component, " link for sanitized PLINK input at ", dest)
  }
}


# Create a lightweight PLINK view around malformed metadata. The GNU AWK helper
# scans clean inputs once and performs a second streaming rewrite only when a
# duplicate allele code would otherwise prevent PLINK2 from loading the file.
sanitize_plink_input_args <- function(block, out_prefix, label = "genotype input") {
  kind <- tolower(block$type)
  if (!kind %in% c("bed", "pgen")) die("unsupported genotype type '", kind, "'. Use 'pgen' or 'bed'.")

  prefix <- block$prefix
  metadata_suffix <- if (kind == "bed") ".bim" else ".pvar"
  binary_suffix <- if (kind == "bed") ".bed" else ".pgen"
  sample_suffix <- if (kind == "bed") ".fam" else ".psam"
  format_label <- toupper(sub("^\\.", "", metadata_suffix))
  invalid_tag <- tolower(format_label)
  safe_prefix <- paste0(out_prefix, ".plink_safe_input")
  safe_metadata <- paste0(safe_prefix, metadata_suffix)
  safe_binary <- paste0(safe_prefix, binary_suffix)
  safe_samples <- paste0(safe_prefix, sample_suffix)
  exclude_path <- paste0(out_prefix, ".invalid_", invalid_tag, "_alleles.exclude.txt")
  report_path <- paste0(out_prefix, ".invalid_", invalid_tag, "_alleles.tsv")
  final_paths <- c(safe_metadata, safe_binary, safe_samples, exclude_path, report_path)

  sanitizer_paths <- if (blank(out_prefix)) NULL else list(
    safe_metadata = safe_metadata,
    report = report_path,
    exclude = exclude_path
  )
  summary <- inspect_plink_metadata(
    block,
    label = label,
    inspect_alleles = TRUE,
    sanitizer_paths = sanitizer_paths
  )
  invalid_count <- plink_metadata_count(summary, "invalid_duplicate_alleles")

  if (!invalid_count) {
    # The helper removes stale text artifacts only after a successful clean
    # scan; remove the corresponding stale component links at the same point.
    if (!blank(out_prefix)) {
      for (path in c(safe_binary, safe_samples)) {
        if (path_exists_or_symlink(path)) unlink(path)
      }
    }
    return(genotype_args(block))
  }
  if (blank(out_prefix)) {
    die("out_prefix is required to sanitize malformed ", format_label, " alleles for ", label)
  }

  successful <- FALSE
  on.exit({
    if (!successful) {
      for (path in final_paths) if (path_exists_or_symlink(path)) unlink(path)
    }
  }, add = TRUE)
  for (path in c(safe_binary, safe_samples)) {
    if (path_exists_or_symlink(path)) unlink(path)
  }
  link_plink_component(paste0(prefix, binary_suffix), safe_binary, toupper(sub("^\\.", "", binary_suffix)))
  link_plink_component(paste0(prefix, sample_suffix), safe_samples, toupper(sub("^\\.", "", sample_suffix)))

  successful <- TRUE
  cat(
    "Excluded", invalid_count, format_label,
    "row(s) with duplicate allele codes before PLINK2 conversion; report:", report_path, "\n"
  )
  c(if (kind == "bed") "--bfile" else "--pfile", safe_prefix, "--exclude", exclude_path)
}


# Convert genotype config to PLINK2 input arguments, sanitizing BIM/PVAR rows
# that PLINK2 cannot parse before normal --exclude/--snps-only filters.
plink_input_args <- function(block, out_prefix = "", label = "genotype input") {
  sanitize_plink_input_args(block, out_prefix, label)
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


qc_info_min <- function(config) config$qc$info_min %||% 0.9


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
  if (truthy(qc$use_mach_r2_filter %||% TRUE)) {
    filters <- c(filters, "--mach-r2-filter", as.character(qc_info_min(config)))
  }
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
