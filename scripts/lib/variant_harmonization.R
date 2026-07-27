# Shared variant-ID harmonization helpers for ancestry-only PGEN copies.
#
# The source genotype files are never rewritten. Large PVARs are streamed in
# chunks, sorted with the system sort utility, and compared one locus at a time.


harmonization_autosomes <- as.character(seq_len(22))


# Return a PGEN prefix whether the caller supplied a prefix or a .pvar path.
pvar_prefix <- function(prefix_or_path) {
  sub("\\.pvar$", "", prefix_or_path)
}


# Stream normalized, ancestry-eligible PVAR rows to a callback.
#
# REF is "PVAR-known" when INFO/PR is absent. This records PLINK metadata; it
# does not independently compare REF with a reference-genome FASTA.
scan_normalized_pvar <- function(prefix_or_path, visit, label = "ancestry PVAR",
                                 chunk_size = 10000L) {
  scanned <- 0L
  eligible <- 0L
  scan_plink_variant_metadata(
    list(type = "pgen", prefix = pvar_prefix(prefix_or_path)),
    function(rows) {
      scanned <<- scanned + nrow(rows)
      chrom_number <- suppressWarnings(as.integer(clean_chrom(trimws(rows$chrom))))
      chrom <- as.character(chrom_number)
      pos <- suppressWarnings(as.integer(rows$pos))
      native_id <- trimws(rows$variant_id)
      ref <- toupper(trimws(rows$allele1))
      alt <- toupper(trimws(rows$allele2))

      keep <- !is.na(chrom_number) & chrom_number >= 1L & chrom_number <= 22L &
        !is.na(pos) & pos > 0L &
        nzchar(native_id) & native_id != "." &
        grepl("^[ACGT]+$", ref) &
        grepl("^[ACGT]+$", alt) &
        ref != alt &
        !grepl(",", alt, fixed = TRUE)

      if (!any(keep)) return(TRUE)
      normalized <- data.frame(
        native_id = native_id[keep],
        chrom = chrom[keep],
        pos = pos[keep],
        ref = ref[keep],
        alt = alt[keep],
        ref_provisional = as.logical(rows$ref_provisional[keep]),
        stringsAsFactors = FALSE
      )
      eligible <<- eligible + nrow(normalized)
      visit(normalized)
    },
    label = label,
    chunk_size = chunk_size,
    require_alleles = TRUE
  )
  list(scanned = scanned, eligible = eligible, excluded_basic = scanned - eligible)
}


# Sort a text file without asking R to hold all records in memory.
sort_text_file <- function(input, output, label) {
  error_path <- paste0(output, ".sort.stderr")
  on.exit(unlink(error_path), add = TRUE)
  status <- suppressWarnings(system2(
    "sort",
    args = shQuote(input, type = "sh"),
    stdout = output,
    stderr = error_path,
    env = "LC_ALL=C"
  ))
  if (!identical(status, 0L)) {
    detail <- if (file.exists(error_path)) paste(readLines(error_path, warn = FALSE), collapse = " ") else ""
    die("could not sort ", label, if (nzchar(detail)) paste0(": ", detail) else "")
  }
}


# Read sorted text through a reusable buffer instead of one disk read per row.
new_buffered_line_reader <- function(path, chunk_size = 10000L) {
  connection <- file(path, open = "rt")
  buffer <- character()
  index <- 1L

  next_line <- function() {
    if (index > length(buffer)) {
      buffer <<- readLines(connection, n = chunk_size, warn = FALSE)
      index <<- 1L
      if (!length(buffer)) return(NULL)
    }
    line <- buffer[[index]]
    index <<- index + 1L
    line
  }

  list(next_line = next_line, close = function() close(connection))
}


# Copy a completed temporary artifact into its declared output path.
copy_harmonization_output <- function(source, target) {
  ensure_parent(target)
  if (!file.copy(source, target, overwrite = TRUE)) {
    die("could not write ancestry harmonization artifact: ", target)
  }
}


# Run the bounded-memory AWK engine and publish artifacts only after every
# parser, duplicate-ID, and harmonized-ID integrity check succeeds.
harmonize_pvar_variants <- function(reference_prefix, study_prefix,
                                    shared_variants, mismatch_report, mapping,
                                    reference_extract, study_extract,
                                    reference_update, study_update,
                                    exclude_palindromic = TRUE) {
  ensure_parent(shared_variants)
  temporary_dir <- tempfile(
    "variant-harmonization-",
    tmpdir = dirname(shared_variants)
  )
  dir.create(temporary_dir)
  on.exit(unlink(temporary_dir, recursive = TRUE), add = TRUE)

  helper <- file.path(script_dir, "harmonize_pvar_variants.sh")
  require_existing_file(helper, "variant harmonization helper")
  helper_args <- c(
    helper,
    "--reference-pvar", paste0(pvar_prefix(reference_prefix), ".pvar"),
    "--study-pvar", paste0(pvar_prefix(study_prefix), ".pvar"),
    "--output-dir", temporary_dir,
    "--exclude-palindromic", if (isTRUE(exclude_palindromic)) "true" else "false"
  )
  helper_log <- file.path(temporary_dir, "harmonization_helper.log")
  status <- suppressWarnings(system2(
    "bash",
    shQuote(helper_args, type = "sh"),
    stdout = "",
    stderr = helper_log
  ))
  helper_messages <- if (file.exists(helper_log)) {
    readLines(helper_log, warn = FALSE)
  } else {
    character()
  }
  if (length(helper_messages)) {
    cat(paste0(helper_messages, "\n"), file = stderr(), sep = "")
  }
  if (!identical(status, 0L)) {
    errors <- helper_messages[startsWith(helper_messages, "ERROR:")]
    detail <- if (length(errors)) paste(errors, collapse = " ") else {
      paste("helper exited with status", status)
    }
    die("variant harmonization failed: ", detail)
  }

  summary_path <- file.path(temporary_dir, "harmonization_summary.tsv")
  summary <- read_tsv(summary_path)
  required_summary <- c(
    "retained", "mismatches",
    "reference_scanned", "reference_eligible",
    "study_scanned", "study_eligible",
    "reference_sort_mode", "study_sort_mode"
  )
  missing_summary <- setdiff(required_summary, names(summary))
  if (length(missing_summary) || nrow(summary) != 1L) {
    die("variant harmonization helper wrote a malformed summary: ", summary_path)
  }

  read_count <- function(name) {
    value <- suppressWarnings(as.integer(summary[[name]][[1]]))
    if (is.na(value) || value < 0L) {
      die("variant harmonization helper wrote an invalid ", name, " count")
    }
    value
  }
  for (name in c("reference_sort_mode", "study_sort_mode")) {
    if (!summary[[name]][[1]] %in% c("stream", "fallback")) {
      die("variant harmonization helper wrote an invalid ", name)
    }
  }

  generated <- c(
    mapping = file.path(temporary_dir, "variant_harmonization.tsv"),
    mismatch_report = file.path(temporary_dir, "shared_variant_mismatches.tsv"),
    shared_variants = file.path(temporary_dir, "shared_variants.txt"),
    reference_extract = file.path(temporary_dir, "reference_native_variants.txt"),
    study_extract = file.path(temporary_dir, "study_native_variants.txt"),
    reference_update = file.path(temporary_dir, "reference_update_names.tsv"),
    study_update = file.path(temporary_dir, "study_update_names.tsv")
  )
  targets <- c(
    mapping = mapping, mismatch_report = mismatch_report,
    shared_variants = shared_variants,
    reference_extract = reference_extract, study_extract = study_extract,
    reference_update = reference_update, study_update = study_update
  )
  for (name in names(generated)) {
    require_existing_file(generated[[name]], paste(name, "harmonization artifact"))
  }
  for (name in names(generated)) {
    copy_harmonization_output(generated[[name]], targets[[name]])
  }

  list(
    retained = read_count("retained"),
    mismatches = read_count("mismatches"),
    reference_scanned = read_count("reference_scanned"),
    reference_eligible = read_count("reference_eligible"),
    study_scanned = read_count("study_scanned"),
    study_eligible = read_count("study_eligible")
  )
}


# Extract by source-native IDs, then rename only the ancestry-local PGEN copy.
extract_and_rename_pgen <- function(config, input_prefix, native_variants,
                                    update_names, out_prefix, threads) {
  ensure_parent(paste0(out_prefix, ".pgen"))
  native_prefix <- tempfile(".native-extract-", tmpdir = dirname(out_prefix))
  native_outputs <- paste0(native_prefix, c(
    ".pgen", ".pvar", ".psam", ".log", ".nosex"
  ))
  on.exit(unlink(native_outputs), add = TRUE)

  run_command(plink_tool(config), c(
    "--pfile", input_prefix,
    "--extract", native_variants,
    "--make-pgen", "--sort-vars",
    "--threads", threads,
    "--out", native_prefix
  ))
  run_command(plink_tool(config), c(
    "--pfile", native_prefix,
    "--update-name", update_names,
    "--make-pgen", "--sort-vars",
    "--threads", threads,
    "--out", out_prefix
  ))
}


# Write sortable ID/locus/allele signatures for validation.
write_harmonized_signatures <- function(prefix, unsorted_path, label) {
  connection <- file(unsorted_path, open = "wt")
  on.exit(close(connection), add = TRUE)
  counts <- scan_normalized_pvar(prefix, function(rows) {
    writeLines(paste(
      rows$native_id,
      sprintf("%02d", as.integer(rows$chrom)),
      sprintf("%012d", rows$pos),
      pmin(rows$ref, rows$alt),
      pmax(rows$ref, rows$alt),
      sep = "\t"
    ), connection)
    TRUE
  }, label = label)
  if (counts$eligible != counts$scanned) {
    die(label, " contains variants outside the harmonized autosomal biallelic ACGT set")
  }
  counts
}


# Confirm that both renamed PGEN copies expose the same operational marker set.
validate_harmonized_pgen_variants <- function(reference_prefix, study_prefix, out) {
  ensure_parent(out)
  temporary_dir <- tempfile("variant-set-validation-", tmpdir = dirname(out))
  dir.create(temporary_dir)
  on.exit(unlink(temporary_dir, recursive = TRUE), add = TRUE)

  reference_unsorted <- file.path(temporary_dir, "reference.unsorted.tsv")
  study_unsorted <- file.path(temporary_dir, "study.unsorted.tsv")
  reference_sorted <- file.path(temporary_dir, "reference.sorted.tsv")
  study_sorted <- file.path(temporary_dir, "study.sorted.tsv")
  reference_counts <- write_harmonized_signatures(
    reference_prefix, reference_unsorted, "renamed reference ancestry PVAR"
  )
  study_counts <- write_harmonized_signatures(
    study_prefix, study_unsorted, "renamed study ancestry PVAR"
  )
  sort_text_file(reference_unsorted, reference_sorted, "renamed reference variant set")
  sort_text_file(study_unsorted, study_sorted, "renamed study variant set")

  reference_reader <- new_buffered_line_reader(reference_sorted)
  study_reader <- new_buffered_line_reader(study_sorted)
  on.exit({
    reference_reader$close()
    study_reader$close()
  }, add = TRUE)

  count <- 0L
  previous_reference_id <- ""
  previous_study_id <- ""
  repeat {
    reference_line <- reference_reader$next_line()
    study_line <- study_reader$next_line()
    if (is.null(reference_line) && is.null(study_line)) break
    if (is.null(reference_line) || is.null(study_line) ||
        !identical(reference_line, study_line)) {
      die("renamed reference and study ancestry PVARs do not have identical ID/locus/allele sets; first difference: reference='",
        reference_line %||% "<EOF>",
        "', study='", study_line %||% "<EOF>", "'")
    }

    reference_id <- strsplit(reference_line, "\t", fixed = TRUE)[[1]][[1]]
    study_id <- strsplit(study_line, "\t", fixed = TRUE)[[1]][[1]]
    if (identical(reference_id, previous_reference_id) ||
        identical(study_id, previous_study_id)) {
      die("renamed ancestry PVAR contains duplicate harmonized variant ID: ", reference_id)
    }
    previous_reference_id <- reference_id
    previous_study_id <- study_id
    count <- count + 1L
  }
  if (!count) die("renamed ancestry PVARs contain no shared variants")
  if (reference_counts$eligible != study_counts$eligible || reference_counts$eligible != count) {
    die("renamed ancestry PVAR variant counts changed during validation")
  }

  writeLines(paste("harmonized_variant_sets_identical", count, sep = "\t"), out)
  count
}


# Stream a PVAR and write IDs inside any configured exclusion region.
write_pvar_region_exclusions <- function(prefix, regions, out,
                                         label = "ancestry PVAR") {
  ensure_parent(out)
  output_connection <- file(out, open = "wt")
  on.exit(close(output_connection), add = TRUE)
  count <- 0L

  scan_normalized_pvar(prefix, function(rows) {
    excluded <- rep(FALSE, nrow(rows))
    for (i in seq_len(nrow(regions))) {
      excluded <- excluded | (
        rows$chrom == regions$chrom[[i]] &
          rows$pos >= as.integer(regions$start[[i]]) &
          rows$pos <= as.integer(regions$end[[i]])
      )
    }
    if (any(excluded)) {
      writeLines(rows$native_id[excluded], output_connection)
      count <<- count + sum(excluded)
    }
    TRUE
  }, label = label)
  count
}
