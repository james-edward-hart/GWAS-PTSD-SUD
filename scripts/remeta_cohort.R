#!/usr/bin/env Rscript

# Cohort-side rare-variant statistics and marginal ReMeta LD helpers.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))

REMETA_HTP_COLUMNS <- c(
  "Name", "Chr", "Pos", "Ref", "Alt", "Trait", "Cohort", "Model", "Effect", "LCI_Effect",
  "UCI_Effect", "Pval", "AAF", "Num_Cases", "Cases_Ref", "Cases_Het", "Cases_Alt",
  "Num_Controls", "Controls_Ref", "Controls_Het", "Controls_Alt", "Info"
)
REMETA_LD_HEADER <- charToRaw("remetaLD.v1.1")


raw <- commandArgs(trailingOnly = TRUE)
if (!length(raw)) die("missing ReMeta cohort subtask")
subtask <- raw[[1]]
args <- parse_args(
  defaults = list(threads = "1", artifact = character(), "target-summary" = character(),
    "trait-summary" = character(), validation = character(), tool = character(),
    htp = character(), index = character(), ld = character(), "sample-ids" = "",
    "ordinary-sample-ids" = "", "group-status" = ""),
  repeated = c("artifact", "target-summary", "trait-summary", "validation", "tool", "htp", "index", "ld"),
  raw = raw[-1]
)


read_tsv_no_metadata <- function(path) {
  if (!file.exists(path)) die("tab-delimited file not found: ", path)
  con <- file(path, open = "rt")
  on.exit(close(con))
  skipped <- 0L
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) die("tab-delimited file is empty: ", path)
    if (!nzchar(trimws(line)) || startsWith(line, "##")) {
      skipped <- skipped + 1L
      next
    }
    break
  }
  seek(con, where = 0L, origin = "start")
  read.table(
    con, skip = skipped, sep = "\t", header = TRUE, check.names = FALSE,
    stringsAsFactors = FALSE, quote = "", comment.char = "", na.strings = character()
  )
}


read_psam_ids <- function(prefix_or_path) {
  path <- if (grepl("\\.psam$", prefix_or_path)) prefix_or_path else paste0(prefix_or_path, ".psam")
  table_sample_ids(read_tsv_no_metadata(path), paste("PSAM", path))
}


read_pvar <- function(prefix_or_path) {
  path <- if (grepl("\\.pvar$", prefix_or_path)) prefix_or_path else paste0(prefix_or_path, ".pvar")
  rows <- read_tsv_no_metadata(path)
  chrom_col <- if ("#CHROM" %in% names(rows)) "#CHROM" else "CHROM"
  require_columns(rows, c(chrom_col, "POS", "ID", "REF", "ALT"), path)
  rows$CHROM_CLEAN <- clean_chrom(rows[[chrom_col]])
  rows
}


assert_identical_samples <- function(keep_path, psam_path) {
  keep <- read_id_file(keep_path, "ReMeta/regenie group keep file")
  psam <- read_psam_ids(psam_path)
  require_unique_ids(keep, "ReMeta/regenie group keep file")
  require_unique_ids(psam, "ReMeta target PSAM")
  keep_key <- sort(paste(keep$FID, keep$IID, sep = "\t"))
  psam_key <- sort(paste(psam$FID, psam$IID, sep = "\t"))
  if (!identical(keep_key, psam_key)) {
    die("ReMeta target samples are not identical to the regenie Step 2 group keep samples")
  }
  length(keep_key)
}


assert_identical_id_files <- function(expected_path, observed_path, label) {
  expected <- read_id_file(expected_path, paste(label, "expected sample IDs"))
  observed <- read_id_file(observed_path, paste(label, "observed sample IDs"))
  require_unique_ids(expected, paste(label, "expected sample IDs"))
  require_unique_ids(observed, paste(label, "observed sample IDs"))
  expected_key <- sort(paste(expected$FID, expected$IID, sep = "\t"))
  observed_key <- sort(paste(observed$FID, observed$IID, sep = "\t"))
  if (!identical(expected_key, observed_key)) die(label, " samples do not match the exact Phase 2 keep")
  length(expected_key)
}


validate_cpra_pvar <- function(prefix_or_path) {
  rows <- read_pvar(prefix_or_path)
  if (!nrow(rows)) die("ReMeta target PVAR contains no variants")
  if ("INFO" %in% names(rows) && any(grepl("(^|;)PR(;|$)", rows$INFO))) {
    die("ReMeta target PVAR contains PLINK2 provisional REF alleles; CPRA requires reference-verified REF/ALT")
  }
  chrom <- rows$CHROM_CLEAN
  pos <- suppressWarnings(as.integer(rows$POS))
  ref <- toupper(as.character(rows$REF))
  alt <- toupper(as.character(rows$ALT))
  if (any(!chrom %in% as.character(seq_len(22)))) die("ReMeta target PVAR contains non-autosomal variants")
  if (any(is.na(pos) | pos < 1)) die("ReMeta target PVAR contains invalid positions")
  if (any(!grepl("^[ACGT]+$", ref)) || any(!grepl("^[ACGT]+$", alt))) {
    die("ReMeta target PVAR contains non-ACGT or symbolic alleles")
  }
  expected <- paste(chrom, pos, ref, alt, sep = ":")
  if (any(rows$ID != expected)) {
    idx <- which(rows$ID != expected)[[1]]
    die("ReMeta target ID is not normalized CPRA; observed ", rows$ID[[idx]], ", expected ", expected[[idx]])
  }
  duplicate <- unique(expected[duplicated(expected)])
  if (length(duplicate)) die("ReMeta target PVAR contains duplicate CPRA IDs: ", paste(head(duplicate, 5), collapse = ", "))
  rows
}


write_key_values <- function(rows, path) {
  write_tsv(data.frame(key = names(rows), value = unname(as.character(rows)), stringsAsFactors = FALSE), path)
}


prepare_target <- function(config, pfile_prefix, keep, regions, out_prefix, summary_out, threads) {
  settings <- config$remeta
  # This exact model keep also governs both REGENIE passes and the resulting LD matrices.
  ids <- read_id_file(keep, "ReMeta/regenie group keep file")
  if (!nrow(ids)) die("ReMeta group has no analyzable traits/samples after Phase 2 phenotype and covariate filtering")
  if (!file.exists(regions) || file.info(regions)$size == 0) die("ReMeta target region BED is empty: ", regions)

  mode <- tolower(settings$genotype_mode %||% "")
  filter_args <- c(
    "--mac", as.character(settings$min_mac %||% 1),
    "--geno", as.character(settings$geno_missing_max %||% 0.05)
  )
  make_args <- c("--make-pgen")
  if (identical(mode, "dosage")) {
    filter_args <- c(
      "--mac", as.character(settings$min_mac %||% 1),
      "--geno", as.character(settings$geno_missing_max %||% 0.05), "dosage",
      "--mach-r2-filter", as.character(settings$info_min %||% 0.8)
    )
  } else if (identical(mode, "hardcall")) {
    make_args <- c("--make-pgen", "erase-dosage")
  } else {
    die("remeta.genotype_mode must be hardcall or dosage")
  }

  run_command(plink_tool(config), c(
    "--pfile", pfile_prefix,
    "--keep", keep,
    "--extract", "bed0", regions,
    "--autosome",
    "--output-chr", "26",
    "--max-alleles", "2",
    filter_args,
    "--set-all-var-ids", shQuote("@:#:$r:$a", type = "sh"),
    "--new-id-max-allele-len", "1000", "missing",
    make_args,
    "--sort-vars",
    "--threads", threads,
    "--out", out_prefix
  ))

  sample_count <- assert_identical_samples(keep, paste0(out_prefix, ".psam"))
  variants <- validate_cpra_pvar(out_prefix)
  rows <- c(
    schema = "remeta_target_summary_v1",
    data_source = settings$data_source,
    genotype_mode = mode,
    input_variants_normalized = as.character(truthy(settings$input_variants_normalized %||% FALSE)),
    sample_count = sample_count,
    variant_count = nrow(variants),
    min_mac = settings$min_mac %||% 1,
    geno_missing_max = settings$geno_missing_max %||% 0.05,
    info_min = if (identical(mode, "dosage")) settings$info_min %||% 0.8 else "NA",
    target_regions = regions,
    target_regions_sha256 = sha256_file(regions),
    group_keep = keep,
    group_keep_sha256 = sha256_file(keep),
    pgen_sha256 = sha256_file(paste0(out_prefix, ".pgen")),
    pvar_sha256 = sha256_file(paste0(out_prefix, ".pvar")),
    psam_sha256 = sha256_file(paste0(out_prefix, ".psam"))
  )
  write_key_values(rows, summary_out)
}


shell_command_line <- function(command, values) {
  paste(c(shQuote(command, type = "sh"), vapply(values, shQuote, character(1), type = "sh")), collapse = " ")
}


write_bash_script <- function(path, lines) {
  ensure_parent(path)
  writeLines(lines, path)
  Sys.chmod(path, mode = "0755")
}


split_options <- function(value) {
  value <- trimws(as.character(value %||% ""))
  if (!nzchar(value)) character() else strsplit(value, "\\s+")[[1]]
}


write_step2_command <- function(config, expected_trait, group_summary, pfile_prefix, keep, pheno, covar,
                                pred_list, trait_list, covar_list, out_prefix, done, script_out, threads) {
  traits <- readLines(trait_list, warn = FALSE)
  traits <- traits[nzchar(traits)]
  if (length(traits) != 1L || !identical(traits, expected_trait)) {
    die("ReMeta Step 2 requires exactly the requested singleton trait: ", expected_trait)
  }
  covars <- trimws(readLines(covar_list, warn = FALSE))
  covars <- covars[nzchar(covars)]
  if (length(covars) != 1) die("Phase 2 group covariate list must contain exactly one non-empty line")
  summary <- read_tsv(group_summary)
  active <- summary[summary$trait %in% traits, , drop = FALSE]
  if (nrow(summary) != 1L || nrow(active) != 1L || identical(active$skipped[[1]], "True")) {
    die("ReMeta Step 2 received a skipped or non-singleton Phase 2 summary")
  }
  types <- unique(active$trait_type)
  if (length(types) != 1 || !types %in% c("bt", "qt")) die("could not determine one Phase 2 trait type for ReMeta group")

  settings <- config$remeta
  phase2 <- config$phase2_regenie
  command <- c(
    "--step", "2",
    "--pgen", pfile_prefix,
    "--keep", keep,
    "--phenoFile", pheno,
    "--phenoColList", paste(traits, collapse = ","),
    "--covarFile", covar,
    "--covarColList", covars,
    "--pred", pred_list,
    "--htp", trimws(as.character(phase2$htp_cohort_name %||% "")),
    if (identical(types, "bt")) "--bt" else "--qt",
    "--bsize", as.character(phase2$step2_bsize %||% 400),
    "--minMAC", as.character(settings$min_mac %||% 1),
    "--write-samples",
    "--threads", threads
  )
  cohort <- command[match("--htp", command) + 1]
  if (!nzchar(cohort)) command[match("--htp", command) + 1] <- analysis_output_name(config)
  if (identical(tolower(settings$genotype_mode), "dosage")) {
    command <- c(command, "--minINFO", as.character(settings$info_min %||% 0.8))
  }
  if (identical(types, "bt")) {
    command <- c(
      command,
      "--minCaseCount", as.character(config$warnings$min_cases %||% 10),
      "--firth", "--approx", "--pThresh", as.character(phase2$p_thresh %||% 0.01)
    )
  } else if (truthy(phase2$apply_rint %||% FALSE)) {
    command <- c(command, "--apply-rint")
  }
  command <- c(command, split_options(phase2$step2_options %||% ""), "--gz", "--out", out_prefix)
  expected <- paste0(out_prefix, "_", traits, ".regenie.gz")
  expected_ids <- paste0(out_prefix, "_", traits, ".regenie.ids")
  lines <- c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    paste("mkdir -p", shQuote(dirname(out_prefix), type = "sh")),
    shell_command_line(config$tools$regenie %||% "regenie", command),
    paste("test -s", vapply(expected, shQuote, character(1), type = "sh")),
    paste("test -s", vapply(expected_ids, shQuote, character(1), type = "sh")),
    paste("printf '%s\\n' ok >", shQuote(done, type = "sh"))
  )
  write_bash_script(script_out, lines)
}


stage_trait <- function(trait, group_summary, raw_stats, sample_ids, keep, out) {
  summary <- read_tsv(group_summary)
  row <- summary[summary$trait == trait, , drop = FALSE]
  if (nrow(row) != 1) die("trait ", trait, " is absent or duplicated in Phase 2 group summary")
  if (identical(row$skipped[[1]], "True")) die("cannot stage a fake ReMeta HTP file for skipped trait ", trait)
  sample_count <- assert_identical_id_files(keep, sample_ids, "rare-variant REGENIE Step 2")
  if (sample_count != as.integer(row$model_sample_count[[1]])) {
    die("rare-variant REGENIE sample count does not match the Phase 2 summary")
  }
  ensure_parent(out)
  if (!file.exists(raw_stats) || file.info(raw_stats)$size == 0) {
    die("expected compressed regenie HTP output not found: ", raw_stats)
  }
  if (!file.copy(raw_stats, out, overwrite = TRUE)) die("could not stage ReMeta regenie HTP output: ", out)
}


read_text_maybe_gzip <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  readLines(con, warn = FALSE)
}


validate_htp <- function(path, variant_ids, ld_variant_ids, expected_trait, trait_type, expected_cases,
                         expected_controls, chunk_size = 100000L) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))

  header <- character()
  while (!length(header)) {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) die("ReMeta HTP file has no header: ", path)
    if (!nzchar(line) || startsWith(line, "##")) next
    header <- strsplit(line, "\t", fixed = TRUE)[[1]]
  }
  if (!identical(header, REMETA_HTP_COLUMNS)) {
    die("ReMeta HTP header does not match the ordered 22-column HTP v4 schema: ", path)
  }
  if (identical(trait_type, "bt") &&
      (any(!is.finite(c(expected_cases, expected_controls))) ||
       any(c(expected_cases, expected_controls) < 0) ||
       any(c(expected_cases, expected_controls) != floor(c(expected_cases, expected_controls))))) {
    die("singleton Phase 2 model has invalid binary case/control counts: ", path)
  }
  column <- setNames(seq_along(header), header)
  variant_rows <- 0L

  repeat {
    lines <- readLines(con, n = chunk_size, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[nzchar(lines) & !startsWith(lines, "##")]
    if (!length(lines)) next
    data_rows <- strsplit(lines, "\t", fixed = TRUE)
    if (any(lengths(data_rows) < length(header))) die("ReMeta HTP file contains truncated data rows: ", path)
    names_seen <- vapply(data_rows, function(row) row[[column[["Name"]]]], character(1))
    traits_seen <- unique(vapply(data_rows, function(row) row[[column[["Trait"]]]], character(1)))
    if (!identical(traits_seen, expected_trait)) die("HTP trait does not match its singleton group: ", path)
    if (identical(trait_type, "bt")) {
      count_columns <- c(
        "Num_Cases", "Cases_Ref", "Cases_Het", "Cases_Alt",
        "Num_Controls", "Controls_Ref", "Controls_Het", "Controls_Alt"
      )
      counts <- do.call(rbind, lapply(data_rows, function(row) {
        suppressWarnings(as.numeric(row[column[count_columns]]))
      }))
      colnames(counts) <- count_columns
      if (any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts))) {
        die("HTP case/control genotype counts must be nonnegative integers: ", path)
      }
      # REGENIE reports nonmissing genotype counts per variant, which may be below the model N.
      if (any(counts[, "Num_Cases"] != rowSums(counts[, c("Cases_Ref", "Cases_Het", "Cases_Alt"), drop = FALSE])) ||
          any(counts[, "Num_Controls"] != rowSums(counts[, c("Controls_Ref", "Controls_Het", "Controls_Alt"), drop = FALSE]))) {
        die("HTP case/control totals do not equal their genotype-category counts: ", path)
      }
      if (any(counts[, "Num_Cases"] > expected_cases) ||
          any(counts[, "Num_Controls"] > expected_controls)) {
        die("HTP case/control counts exceed the singleton Phase 2 model counts: ", path)
      }
    }
    cpra_seen <- vapply(data_rows, function(row) {
      paste(clean_chrom(row[[column[["Chr"]]]]), row[[column[["Pos"]]]],
        toupper(row[[column[["Ref"]]]]), toupper(row[[column[["Alt"]]]]), sep = ":")
    }, character(1))
    if (any(names_seen != cpra_seen)) {
      idx <- which(names_seen != cpra_seen)[[1]]
      die("HTP Name is not consistent with its Chr/Pos/Ref/Alt columns: observed ",
        names_seen[[idx]], ", expected ", cpra_seen[[idx]])
    }
    info <- vapply(data_rows, function(row) row[[column[["Info"]]]], character(1))
    if (any(!grepl("(^|;)SCORE=", info))) {
      die("HTP data rows must contain the regenie SCORE field required by ReMeta gene tests: ", path)
    }
    missing <- setdiff(unique(names_seen), variant_ids)
    if (length(missing)) die("HTP variants are absent from the sample-matched LD target PVAR: ", paste(head(missing, 5), collapse = ", "))
    missing_ld <- setdiff(unique(names_seen), ld_variant_ids)
    if (length(missing_ld)) die("HTP variants are absent from every ReMeta LD index: ", paste(head(missing_ld, 5), collapse = ", "))
    variant_rows <- variant_rows + length(names_seen)
  }
  variant_rows
}


ld_component_table <- function(paths) {
  if (length(paths) != 66L || anyDuplicated(paths)) {
    die("ReMeta validation requires exactly 22 unique three-file LD component sets")
  }
  rows <- data.frame(path = paths, chrom = NA_integer_, component = "", stringsAsFactors = FALSE)
  pattern <- "^chr([0-9]+)\\.remeta\\.(gene\\.ld|buffer\\.ld|ld\\.idx\\.gz)$"
  for (i in seq_along(paths)) {
    match <- regmatches(basename(paths[[i]]), regexec(pattern, basename(paths[[i]])))[[1]]
    if (length(match) != 3L) die("unexpected ReMeta LD component filename: ", paths[[i]])
    rows$chrom[[i]] <- as.integer(match[[2]])
    rows$component[[i]] <- match[[3]]
    if (!file.exists(paths[[i]])) die("ReMeta LD component is missing: ", paths[[i]])
  }
  expected <- expand.grid(
    chrom = seq_len(22), component = c("gene.ld", "buffer.ld", "ld.idx.gz"),
    stringsAsFactors = FALSE
  )
  if (!setequal(paste(rows$chrom, rows$component), paste(expected$chrom, expected$component))) {
    die("ReMeta LD components do not contain exactly one three-file set for chromosomes 1-22")
  }
  rows
}


split_ld_index_fields <- function(line, path) {
  # Appending a sentinel preserves the required empty fifth field after the final tab.
  fields <- strsplit(paste0(line, "\t__END__"), "\t", fixed = TRUE)[[1]]
  fields <- fields[-length(fields)]
  if (length(fields) != 5L) die("malformed ReMeta LD index row in ", path)
  if (nzchar(fields[[5]])) {
    die("ReMeta LD index contains buffer variants despite the configured no-buffer policy: ", path)
  }
  fields
}


parse_ld_offset <- function(value, label, path) {
  if (!grepl("^(0|[1-9][0-9]*)$", value)) die("invalid ", label, " in ReMeta LD index: ", path)
  parsed <- suppressWarnings(as.numeric(value))
  if (!is.finite(parsed) || parsed > 2^53 || parsed != floor(parsed)) {
    die("invalid ", label, " in ReMeta LD index: ", path)
  }
  parsed
}


raw_uint_le <- function(value) {
  sum(as.numeric(value) * 256^(seq_along(value) - 1L))
}


read_raw_at <- function(con, offset, count, label, path) {
  seek(con, where = offset, origin = "start", rw = "read")
  value <- readBin(con, "raw", n = count)
  if (length(value) != count) die("truncated ", label, " in BGZF ReMeta LD file: ", path)
  value
}


scan_bgzf <- function(path) {
  physical_size <- file.info(path)$size
  con <- file(path, "rb")
  on.exit(close(con))
  block_rows <- list()
  physical_offset <- 0
  logical_offset <- 0
  while (physical_offset < physical_size) {
    base <- read_raw_at(con, physical_offset, 12L, "BGZF header", path)
    if (!identical(as.integer(base[1:4]), c(31L, 139L, 8L, 4L))) {
      die("invalid BGZF header in ReMeta LD file: ", path)
    }
    extra_length <- raw_uint_le(base[11:12])
    extra <- read_raw_at(con, physical_offset + 12, extra_length, "BGZF extra header", path)
    cursor <- 1L
    block_size <- NA_real_
    while (cursor + 3L <= length(extra)) {
      field_length <- raw_uint_le(extra[(cursor + 2L):(cursor + 3L)])
      field_end <- cursor + 3L + field_length
      if (field_end > length(extra)) die("malformed BGZF extra header in ReMeta LD file: ", path)
      if (rawToChar(extra[cursor:(cursor + 1L)]) == "BC" && field_length == 2L) {
        block_size <- raw_uint_le(extra[(cursor + 4L):(cursor + 5L)]) + 1
      }
      cursor <- field_end + 1L
    }
    if (!is.finite(block_size) || block_size < 12 + extra_length + 8 ||
        physical_offset + block_size > physical_size) {
      die("invalid BGZF block size in ReMeta LD file: ", path)
    }
    block <- read_raw_at(con, physical_offset, block_size, "BGZF block", path)
    logical_size <- raw_uint_le(tail(block, 4L))
    uncompressed <- tryCatch(memDecompress(block, type = "gzip"), error = function(err) raw())
    if (length(uncompressed) != logical_size || logical_size > 65536) {
      die("corrupt BGZF block in ReMeta LD file: ", path)
    }
    block_rows[[length(block_rows) + 1L]] <- data.frame(
      physical_offset = physical_offset, physical_size = block_size,
      logical_offset = logical_offset, logical_size = logical_size
    )
    physical_offset <- physical_offset + block_size
    logical_offset <- logical_offset + logical_size
  }
  blocks <- if (length(block_rows)) do.call(rbind, block_rows) else data.frame()
  if (!nrow(blocks) || physical_offset != physical_size || tail(blocks$logical_size, 1L) != 0) {
    die("ReMeta LD file is missing its BGZF terminator: ", path)
  }
  list(path = path, blocks = blocks, logical_size = logical_offset)
}


open_bgzf_reader <- function(bgzf) {
  reader <- new.env(parent = emptyenv())
  reader$path <- bgzf$path
  reader$blocks <- bgzf$blocks
  reader$logical_size <- bgzf$logical_size
  reader$con <- file(bgzf$path, "rb")
  reader$cached_row <- NA_integer_
  reader$cached_block <- raw()
  reader
}


read_bgzf_logical <- function(reader, offset, count, label) {
  if (offset < 0 || count < 0 || offset + count > reader$logical_size) {
    die("out-of-bounds ", label, " in BGZF ReMeta LD file: ", reader$path)
  }
  if (!count) return(raw())
  value <- raw()
  remaining <- count
  position <- offset
  while (remaining > 0) {
    matches <- which(
      reader$blocks$logical_size > 0 &
      position >= reader$blocks$logical_offset &
      position < reader$blocks$logical_offset + reader$blocks$logical_size
    )
    if (!length(matches)) die("could not map ", label, " to a BGZF block: ", reader$path)
    row <- matches[[1]]
    if (!identical(reader$cached_row, row)) {
      block <- read_raw_at(reader$con, reader$blocks$physical_offset[[row]],
        reader$blocks$physical_size[[row]], label, reader$path)
      reader$cached_block <- memDecompress(block, type = "gzip")
      reader$cached_row <- row
    }
    within <- position - reader$blocks$logical_offset[[row]]
    take <- min(remaining, length(reader$cached_block) - within)
    value <- c(value, reader$cached_block[within + seq_len(take)])
    position <- position + take
    remaining <- remaining - take
  }
  value
}


bgzf_virtual_to_logical <- function(bgzf, virtual_offset, label) {
  physical_offset <- floor(virtual_offset / 65536)
  within <- virtual_offset %% 65536
  row <- match(physical_offset, bgzf$blocks$physical_offset)
  if (is.na(row) || within >= bgzf$blocks$logical_size[[row]]) {
    die("invalid ", label, " in ReMeta LD index for ", bgzf$path)
  }
  bgzf$blocks$logical_offset[[row]] + within
}


read_bgzf_int32 <- function(reader, offset, label) {
  value <- raw_uint_le(read_bgzf_logical(reader, offset, 4L, label))
  if (value >= 2^31) value <- value - 2^32
  as.integer(value)
}


read_bgzf_float32 <- function(reader, offset, label) {
  con <- rawConnection(read_bgzf_logical(reader, offset, 4L, label), "rb")
  on.exit(close(con))
  value <- readBin(con, numeric(), n = 1L, size = 4L, endian = "little")
  if (length(value) != 1L || !is.finite(value)) die("invalid ", label, " in ReMeta LD file: ", reader$path)
  value
}


validate_ld_binary <- function(gene_path, buffer_path, index_records, chrom) {
  gene <- open_bgzf_reader(scan_bgzf(gene_path))
  buffer <- open_bgzf_reader(scan_bgzf(buffer_path))
  on.exit(close(gene$con), add = TRUE)
  on.exit(close(buffer$con), add = TRUE)
  if (!identical(read_bgzf_logical(gene, 0, length(REMETA_LD_HEADER), "file header"), REMETA_LD_HEADER)) {
    die("invalid remetaLD.v1.1 header in ", gene_path)
  }
  if (!identical(read_bgzf_logical(buffer, 0, length(REMETA_LD_HEADER), "file header"), REMETA_LD_HEADER)) {
    die("invalid remetaLD.v1.1 header in ", buffer_path)
  }
  float_size <- read_bgzf_int32(buffer, 13, "buffer float size")
  if (!float_size %in% c(1L, 2L, 4L)) die("unsupported correlation size in ReMeta buffer LD file: ", buffer_path)

  if (!length(index_records)) {
    if (gene$logical_size != 13 || buffer$logical_size != 17) {
      die("ReMeta LD files contain unindexed data for chromosome ", chrom)
    }
    return(invisible(NULL))
  }
  for (i in seq_along(index_records)) {
    index_records[[i]]$gene_logical_offset <- bgzf_virtual_to_logical(
      gene, index_records[[i]]$gene_offset, "gene offset"
    )
    index_records[[i]]$buffer_logical_offset <- bgzf_virtual_to_logical(
      buffer, index_records[[i]]$buffer_offset, "buffer offset"
    )
  }

  gene_order <- order(vapply(index_records, `[[`, numeric(1), "gene_logical_offset"))
  expected_gene_offset <- 13
  for (record in index_records[gene_order]) {
    offset <- record$gene_logical_offset
    if (offset != expected_gene_offset || offset >= gene$logical_size) {
      die("ReMeta gene LD offsets do not tile the indexed file for chromosome ", chrom)
    }
    variant_count <- read_bgzf_int32(gene, offset, "gene variant count")
    entry_count <- read_bgzf_int32(gene, offset + 4, "gene correlation count")
    read_bgzf_float32(gene, offset + 8, "gene sparsity threshold")
    if (variant_count != length(record$variant_ids) || entry_count < 0 ||
        entry_count > as.numeric(variant_count) * as.numeric(variant_count)) {
      die("ReMeta gene LD record disagrees with its index: ", record$gene)
    }

    cross_offset <- offset + 12 + 4 * variant_count + 12 * entry_count
    if (cross_offset + 20 > gene$logical_size) die("truncated ReMeta gene LD record: ", record$gene)
    cross_gene_count <- read_bgzf_int32(gene, cross_offset, "gene-buffer gene count")
    buffer_variant_count <- read_bgzf_int32(gene, cross_offset + 4, "buffer variant count")
    cross_entry_count <- read_bgzf_int32(gene, cross_offset + 8, "gene-buffer correlation count")
    read_bgzf_float32(gene, cross_offset + 12, "gene-buffer sparsity threshold")
    if (cross_gene_count != variant_count || buffer_variant_count != 0L || cross_entry_count != 0L) {
      die("ReMeta gene-buffer LD record violates the configured no-buffer policy: ", record$gene)
    }
    buffer_block_offset <- cross_offset + 16 + 8 * buffer_variant_count
    buffer_block_count <- read_bgzf_int32(gene, buffer_block_offset, "buffer block count")
    if (buffer_block_count < 0L) die("invalid buffer block count in ReMeta LD record: ", record$gene)
    expected_gene_offset <- buffer_block_offset + 4 + (8 + float_size) * cross_entry_count
    if (expected_gene_offset > gene$logical_size) die("truncated ReMeta gene LD record: ", record$gene)
  }
  if (expected_gene_offset != gene$logical_size) {
    die("ReMeta gene LD file contains trailing or unindexed bytes for chromosome ", chrom)
  }

  buffer_order <- order(vapply(index_records, `[[`, numeric(1), "buffer_logical_offset"))
  expected_buffer_offset <- 21
  for (record in index_records[buffer_order]) {
    offset <- record$buffer_logical_offset
    if (offset != expected_buffer_offset || offset + 4 > buffer$logical_size) {
      die("ReMeta buffer LD offsets do not tile the indexed file for chromosome ", chrom)
    }
    if (read_bgzf_int32(buffer, offset - 4, "buffer block end marker") != -3L) {
      die("invalid ReMeta buffer LD block-end marker for gene ", record$gene)
    }
    if (read_bgzf_int32(buffer, offset, "buffer block start marker") != -1L) {
      die("invalid ReMeta buffer LD block-start marker for gene ", record$gene)
    }
    expected_buffer_offset <- offset + 8
  }
  if (buffer$logical_size != expected_buffer_offset - 4) {
    die("ReMeta buffer LD file contains trailing or unindexed bytes for chromosome ", chrom)
  }
  invisible(NULL)
}


validate_group <- function(target_prefix, keep, ordinary_ids, rare_ids, group_summary, gene_list,
                           htp_paths, ld_paths, out) {
  summary <- read_tsv(group_summary)
  if (nrow(summary) != 1L || identical(summary$skipped[[1]], "True")) {
    die("ReMeta validation requires one active singleton Phase 2 summary")
  }
  require_columns(summary, c(
    "group", "trait", "trait_type", "model_sample_count", "model_cases", "model_controls"
  ), group_summary)
  if (length(htp_paths) != 1L) die("ReMeta singleton validation requires exactly one HTP file")
  # Compare sets, not row order: all four artifacts must represent the same analyzed people.
  ordinary_sample_count <- assert_identical_id_files(keep, ordinary_ids, "ordinary REGENIE Step 2")
  rare_sample_count <- assert_identical_id_files(keep, rare_ids, "rare-variant REGENIE Step 2")
  sample_count <- assert_identical_samples(keep, paste0(target_prefix, ".psam"))
  if (any(c(ordinary_sample_count, rare_sample_count, sample_count) != as.integer(summary$model_sample_count[[1]]))) {
    die("REGENIE/target sample counts do not match the singleton Phase 2 model summary")
  }
  variants <- validate_cpra_pvar(target_prefix)
  variant_ids <- variants$ID
  genes <- read.table(gene_list, sep = "", header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (ncol(genes) < 4 || !nrow(genes)) die("invalid ReMeta gene list: ", gene_list)
  gene_ids <- as.character(genes[[1]])
  if (anyDuplicated(gene_ids)) die("ReMeta gene list contains duplicate gene IDs")
  gene_chrom <- clean_chrom(genes[[2]])
  gene_start <- suppressWarnings(as.integer(genes[[3]]))
  gene_end <- suppressWarnings(as.integer(genes[[4]]))
  if (any(!gene_chrom %in% as.character(seq_len(22))) ||
      any(is.na(gene_start) | is.na(gene_end) | gene_start < 1L | gene_end < gene_start)) {
    die("ReMeta gene list contains invalid autosomal coordinates")
  }
  gene_row <- setNames(seq_along(gene_ids), gene_ids)
  variant_row <- setNames(seq_along(variant_ids), variant_ids)
  variant_pos <- suppressWarnings(as.integer(variants$POS))

  indexed_genes <- 0L
  indexed_gene_ids <- character()
  ld_variant_ids <- character()
  ld_gene_variant_assignments <- 0L
  ld_assignments_within_gene_bounds <- 0L
  ld_components <- ld_component_table(ld_paths)
  index_rows <- ld_components[ld_components$component == "ld.idx.gz", , drop = FALSE]
  index_rows <- index_rows[order(index_rows$chrom), , drop = FALSE]
  for (path_i in seq_len(nrow(index_rows))) {
    path <- index_rows$path[[path_i]]
    path_chrom <- as.character(index_rows$chrom[[path_i]])
    lines <- read_text_maybe_gzip(path)
    lines <- lines[nzchar(lines)]
    index_records <- list()
    for (line in lines) {
      fields <- split_ld_index_fields(line, path)
      if (!fields[[1]] %in% gene_ids) die("LD index contains gene absent from bundled gene list: ", fields[[1]])
      if (fields[[1]] %in% indexed_gene_ids) die("ReMeta LD indexes contain a duplicate gene row: ", fields[[1]])
      ids <- strsplit(fields[[4]], ",", fixed = TRUE)[[1]]
      ids <- ids[nzchar(ids)]
      if (!length(ids)) die("ReMeta LD index gene has no target variants: ", fields[[1]])
      if (anyDuplicated(ids)) die("ReMeta LD index gene contains duplicate target variants: ", fields[[1]])
      index_records[[length(index_records) + 1L]] <- list(
        gene = fields[[1]],
        gene_offset = parse_ld_offset(fields[[2]], "gene offset", path),
        buffer_offset = parse_ld_offset(fields[[3]], "buffer offset", path),
        variant_ids = ids
      )
      missing <- setdiff(ids, variant_ids)
      if (length(missing)) die("LD index variants are absent from target PVAR: ", paste(head(missing, 5), collapse = ", "))
      gene_idx <- unname(gene_row[[fields[[1]]]])
      if (!identical(gene_chrom[[gene_idx]], path_chrom)) {
        die("ReMeta LD index gene is stored under the wrong chromosome: ", fields[[1]])
      }
      variant_idx <- unname(variant_row[ids])
      if (any(variants$CHROM_CLEAN[variant_idx] != path_chrom)) {
        die("ReMeta LD index contains a variant under the wrong chromosome: ", fields[[1]])
      }
      within_gene <- variants$CHROM_CLEAN[variant_idx] == gene_chrom[[gene_idx]] &
        variant_pos[variant_idx] >= gene_start[[gene_idx]] &
        variant_pos[variant_idx] <= gene_end[[gene_idx]]
      ld_variant_ids <- c(ld_variant_ids, ids)
      indexed_gene_ids <- c(indexed_gene_ids, fields[[1]])
      indexed_genes <- indexed_genes + 1L
      ld_gene_variant_assignments <- ld_gene_variant_assignments + length(ids)
      ld_assignments_within_gene_bounds <- ld_assignments_within_gene_bounds + sum(within_gene)
    }
    gene_path <- ld_components$path[
      ld_components$chrom == as.integer(path_chrom) & ld_components$component == "gene.ld"
    ]
    buffer_path <- ld_components$path[
      ld_components$chrom == as.integer(path_chrom) & ld_components$component == "buffer.ld"
    ]
    validate_ld_binary(gene_path[[1]], buffer_path[[1]], index_records, path_chrom)
  }
  if (!indexed_genes) die("ReMeta LD indexes contain no genes")
  ld_variant_ids <- unique(ld_variant_ids)
  indexed_gene_count <- length(unique(indexed_gene_ids))
  target_variants_not_indexed <- length(setdiff(variant_ids, ld_variant_ids))
  ld_assignments_outside_gene_bounds <- ld_gene_variant_assignments - ld_assignments_within_gene_bounds
  percent <- function(numerator, denominator) {
    if (!denominator) return("NA")
    sprintf("%.6f", 100 * numerator / denominator)
  }

  htp_variants <- 0L
  for (path in htp_paths) {
    htp_variants <- htp_variants + validate_htp(
      path, variant_ids, ld_variant_ids, summary$trait[[1]], summary$trait_type[[1]],
      suppressWarnings(as.numeric(summary$model_cases[[1]])),
      suppressWarnings(as.numeric(summary$model_controls[[1]]))
    )
  }
  write_key_values(c(
    status = "validated",
    group = summary$group[[1]],
    trait = summary$trait[[1]],
    trait_type = summary$trait_type[[1]],
    sample_count = sample_count,
    ordinary_regenie_sample_count = ordinary_sample_count,
    rare_regenie_sample_count = rare_sample_count,
    group_keep_sha256 = sha256_file(keep),
    target_psam_sha256 = sha256_file(paste0(target_prefix, ".psam")),
    htp_sha256 = sha256_file(htp_paths[[1]]),
    target_variant_count = nrow(variants),
    unique_ld_target_variant_count = length(ld_variant_ids),
    target_variants_not_indexed = target_variants_not_indexed,
    target_variant_ld_coverage_pct = percent(length(ld_variant_ids), nrow(variants)),
    htp_variant_rows = htp_variants,
    reference_gene_count = length(gene_ids),
    indexed_gene_rows = indexed_genes,
    indexed_gene_count = indexed_gene_count,
    genes_without_indexed_variants = length(gene_ids) - indexed_gene_count,
    indexed_gene_coverage_pct = percent(indexed_gene_count, length(gene_ids)),
    ld_gene_variant_assignments = ld_gene_variant_assignments,
    ld_assignments_within_gene_bounds = ld_assignments_within_gene_bounds,
    ld_assignments_outside_gene_bounds = ld_assignments_outside_gene_bounds,
    ld_assignment_gene_bound_coverage_pct = percent(
      ld_assignments_within_gene_bounds, ld_gene_variant_assignments
    )
  ), out)
}


write_manifest <- function(config, config_path, build, gene_list, provenance, group_status, target_summaries,
                           trait_summaries, validations, tools, artifacts, out) {
  for (path in c(gene_list, provenance, group_status, target_summaries, trait_summaries, validations, tools, artifacts)) {
    if (!file.exists(path)) die("ReMeta manifest input not found: ", path)
  }
  status <- read_tsv(group_status)
  require_columns(status, c(
    "group", "trait", "trait_type", "usable_n", "cases", "controls", "model_sample_count",
    "model_cases", "model_controls", "keep_count", "model_keep_sha256",
    "skipped", "skip_reason", "remeta_eligible"
  ), group_status)
  if (!nrow(status) || anyDuplicated(status$group) || anyDuplicated(status$trait)) {
    die("ReMeta group status must contain unique group and trait rows")
  }
  if (any(!status$skipped %in% c("True", "False")) ||
      any(!status$remeta_eligible %in% c("True", "False"))) {
    die("ReMeta group status contains invalid active/skipped values")
  }
  status$skip_reason <- as.character(status$skip_reason)
  status$skip_reason[is.na(status$skip_reason)] <- ""
  skipped <- status$skipped == "True"
  eligible <- status$remeta_eligible == "True"
  if (any(status$keep_count != status$model_sample_count) || any(eligible == skipped) ||
      any(skipped != nzchar(status$skip_reason))) {
    die("ReMeta group status has inconsistent keep counts or active/skipped flags")
  }
  if (any(eligible & status$usable_n != status$model_sample_count)) {
    die("active ReMeta group status has different usable and model sample counts")
  }
  active <- status[status$remeta_eligible == "True" & status$skipped == "False", , drop = FALSE]
  target_groups <- vapply(target_summaries, function(path) basename(dirname(path)), character(1))
  validation_groups <- vapply(validations, function(path) basename(dirname(path)), character(1))
  if (anyDuplicated(target_groups) || anyDuplicated(validation_groups) ||
      length(target_groups) != nrow(active) || length(validation_groups) != nrow(active) ||
      !setequal(target_groups, active$group) || !setequal(validation_groups, active$group)) {
    die("ReMeta target summaries and validations must match the active singleton groups")
  }
  artifact_suffix <- function(path) {
    normalized <- gsub("\\\\", "/", path)
    marker <- paste0("/export/", build, "/")
    if (grepl(marker, normalized, fixed = TRUE)) return(sub(paste0("^.*", marker), "", normalized))
    sub(paste0("^results/remeta/export/", build, "/"), "", normalized)
  }
  expected_artifacts <- character()
  # A skipped trait is represented by manifest metadata only; fake empty exports are prohibited.
  for (i in seq_len(nrow(active))) {
    expected_artifacts <- c(expected_artifacts, paste0("htp/", active$trait[[i]], ".PAN.regenie.gz"))
    for (chrom in seq_len(22)) {
      prefix <- paste0("ld/", active$group[[i]], "/chr", chrom)
      expected_artifacts <- c(expected_artifacts, paste0(prefix, c(
        ".remeta.gene.ld", ".remeta.buffer.ld", ".remeta.ld.idx.gz"
      )))
    }
  }
  observed_artifacts <- vapply(artifacts, artifact_suffix, character(1))
  if (anyDuplicated(observed_artifacts) || length(observed_artifacts) != length(expected_artifacts) ||
      !setequal(observed_artifacts, expected_artifacts)) {
    die("ReMeta export artifacts do not exactly match the active trait-specific HTP and LD sets")
  }
  export_dir <- ""
  if (length(artifacts)) {
    normalized <- gsub("\\\\", "/", artifacts[[1]])
    marker <- paste0("/export/", build, "/")
    if (grepl(marker, normalized, fixed = TRUE)) {
      export_dir <- sub(paste0("^(.*?/export/", build, ").*$"), "\\1", normalized)
    } else if (startsWith(normalized, paste0("results/remeta/export/", build, "/"))) {
      export_dir <- file.path("results/remeta/export", build)
    }
  } else {
    normalized_out <- gsub("\\\\", "/", out)
    if (startsWith(normalized_out, "results/remeta/export/") ||
        grepl("/results/remeta/export/", normalized_out, fixed = TRUE)) {
      export_dir <- file.path(dirname(out), build)
    }
  }
  if (nzchar(export_dir) && dir.exists(export_dir)) {
    present <- list.files(export_dir, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    present <- present[file.exists(present) & !dir.exists(present)]
    # A prior copy of this manifest is safe to replace; no other undeclared
    # export may survive a trait-registry or eligibility change.
    present <- present[gsub("\\\\", "/", present) != gsub("\\\\", "/", out)]
    present <- substring(gsub("\\\\", "/", present), nchar(gsub("\\\\", "/", export_dir)) + 2L)
    if (!setequal(present, expected_artifacts)) {
      die("ReMeta export directory contains stale or unexpected artifacts; archive it before rerunning")
    }
  }
  if (nzchar(export_dir)) {
    work_groups_dir <- file.path(dirname(dirname(export_dir)), "work", build, "groups")
    if (dir.exists(work_groups_dir)) {
      work_files <- list.files(work_groups_dir, recursive = TRUE, full.names = TRUE,
        all.files = TRUE, no.. = TRUE)
      work_files <- work_files[file.exists(work_files) & !dir.exists(work_files)]
      relative <- substring(gsub("\\\\", "/", work_files),
        nchar(gsub("\\\\", "/", work_groups_dir)) + 2L)
      work_groups <- sub("/.*$", "", relative)
      if (any(!work_groups %in% active$group)) {
        die("ReMeta work directory contains skipped or unexpected group artifacts; archive it before rerunning")
      }
    }
  }
  metadata <- c(
    manifest_schema = "remeta_cohort_export_v2",
    grouping_schema = "trait_specific_v1",
    created_utc = format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%S+00:00", tz = "UTC"),
    analysis_name = config$project$analysis_name,
    cohort_data_release = config$project$cohort_data_release %||% "",
    genome_build = build,
    analysis_scope = "marginal_gene_tests_only",
    conditional_buffer_included = "false",
    data_source = config$remeta$data_source,
    genotype_mode = config$remeta$genotype_mode,
    input_variants_normalized = as.character(truthy(config$remeta$input_variants_normalized %||% FALSE)),
    target_r2 = config$remeta$target_r2,
    min_mac = config$remeta$min_mac,
    config = config_path,
    config_sha256 = sha256_file(config_path),
    gene_list = gene_list,
    gene_list_sha256 = sha256_file(gene_list),
    resource_provenance = provenance,
    resource_provenance_sha256 = sha256_file(provenance),
    group_status = group_status,
    group_status_sha256 = sha256_file(group_status),
    configured_trait_count = nrow(status),
    exported_trait_count = nrow(active)
  )
  rows <- data.frame(key = names(metadata), value = unname(as.character(metadata)), stringsAsFactors = FALSE)

  target_psam_hashes <- character()
  for (path in target_summaries) {
    values <- read_tsv(path)
    require_columns(values, c("key", "value"), path)
    group <- basename(dirname(path))
    required <- c("schema", "sample_count", "group_keep_sha256", "psam_sha256")
    if (anyDuplicated(values$key) || any(!required %in% values$key)) {
      die("ReMeta target summary is incomplete or duplicated: ", path)
    }
    value <- setNames(as.character(values$value), values$key)
    status_row <- active[active$group == group, , drop = FALSE]
    if (nrow(status_row) != 1L || value[["schema"]] != "remeta_target_summary_v1" ||
        suppressWarnings(as.integer(value[["sample_count"]])) != status_row$model_sample_count[[1]] ||
        value[["group_keep_sha256"]] != status_row$model_keep_sha256[[1]] ||
        !grepl("^[0-9a-f]{64}$", value[["psam_sha256"]])) {
      die("ReMeta target summary disagrees with its active singleton group: ", path)
    }
    target_psam_hashes[[group]] <- value[["psam_sha256"]]
    rows <- rbind(rows, data.frame(
      key = paste0("group:", group, ":", values$key), value = values$value, stringsAsFactors = FALSE
    ))
  }
  trait_tables <- lapply(trait_summaries, function(path) {
    values <- read_tsv(path)
    require_columns(values, c("group", "trait", "trait_type", "phase2_pan_samples",
      "usable_n", "cases", "controls", "model_sample_count", "model_cases", "model_controls",
      "model_keep_sha256", "skipped", "skip_reason"), path)
    if (nrow(values) != 1L) die("each ReMeta trait summary must contain exactly one singleton row: ", path)
    values$skip_reason <- as.character(values$skip_reason)
    values$skip_reason[is.na(values$skip_reason)] <- ""
    values
  })
  trait_rows <- if (length(trait_tables)) do.call(rbind, trait_tables) else data.frame()
  if (nrow(trait_rows) != nrow(status) || anyDuplicated(trait_rows$group) || anyDuplicated(trait_rows$trait) ||
      !setequal(paste(trait_rows$group, trait_rows$trait), paste(status$group, status$trait))) {
    die("ReMeta trait summaries must contain every configured singleton group exactly once")
  }
  for (values in trait_tables) {
    for (i in seq_len(nrow(values))) {
      status_row <- status[status$trait == values$trait[[i]] & status$group == values$group[[i]], , drop = FALSE]
      if (nrow(status_row) != 1L) die("trait summary does not match the Phase 2 group status: ", values$trait[[i]])
      shared_fields <- c(
        "trait_type", "usable_n", "cases", "controls", "model_sample_count", "model_cases",
        "model_controls", "model_keep_sha256", "skipped", "skip_reason"
      )
      if (!identical(as.character(values[i, shared_fields]), as.character(status_row[1, shared_fields]))) {
        die("trait summary disagrees with the Phase 2 group status: ", values$trait[[i]])
      }
      exported <- identical(status_row$remeta_eligible[[1]], "True") &&
        identical(status_row$skipped[[1]], "False")
      prefix <- paste0("trait:", values$trait[[i]], ":")
      trait_metadata <- c(
        phase2_group = values$group[[i]],
        ld_group = ifelse(exported, values$group[[i]], ""),
        trait_type = values$trait_type[[i]],
        phase2_pan_samples = values$phase2_pan_samples[[i]],
        candidate_n = values$usable_n[[i]],
        candidate_cases = values$cases[[i]],
        candidate_controls = values$controls[[i]],
        model_sample_count = values$model_sample_count[[i]],
        model_cases = values$model_cases[[i]],
        model_controls = values$model_controls[[i]],
        model_keep_sha256 = values$model_keep_sha256[[i]],
        skipped = values$skipped[[i]],
        skip_reason = values$skip_reason[[i]],
        remeta_exported = ifelse(exported, "true", "false"),
        htp = ifelse(exported, paste0("results/remeta/export/", build, "/htp/", values$trait[[i]], ".PAN.regenie.gz"), ""),
        ld_prefix = ifelse(exported, paste0("results/remeta/export/", build, "/ld/", values$group[[i]], "/chr{1-22}"), "")
      )
      rows <- rbind(rows, data.frame(
        key = paste0(prefix, names(trait_metadata)),
        value = unname(as.character(trait_metadata)), stringsAsFactors = FALSE
      ))
    }
  }
  for (path in validations) {
    values <- read_tsv(path)
    require_columns(values, c("key", "value"), path)
    group <- basename(dirname(path))
    required <- c(
      "status", "group", "trait", "sample_count", "ordinary_regenie_sample_count",
      "rare_regenie_sample_count", "group_keep_sha256", "target_psam_sha256"
    )
    if (anyDuplicated(values$key) || any(!required %in% values$key)) {
      die("ReMeta validation summary is incomplete or duplicated: ", path)
    }
    value <- setNames(as.character(values$value), values$key)
    status_row <- active[active$group == group, , drop = FALSE]
    expected_n <- if (nrow(status_row)) as.integer(status_row$model_sample_count[[1]]) else NA_integer_
    observed_n <- suppressWarnings(as.integer(value[c(
      "sample_count", "ordinary_regenie_sample_count", "rare_regenie_sample_count"
    )]))
    if (nrow(status_row) != 1L || value[["status"]] != "validated" ||
        value[["group"]] != group || value[["trait"]] != status_row$trait[[1]] ||
        any(is.na(observed_n)) || any(observed_n != expected_n) ||
        value[["group_keep_sha256"]] != status_row$model_keep_sha256[[1]] ||
        value[["target_psam_sha256"]] != target_psam_hashes[[group]]) {
      die("ReMeta validation summary disagrees with its active singleton group: ", path)
    }
    rows <- rbind(rows, data.frame(
      key = c(
        paste0("validation:", group, ":path"),
        paste0("validation:", group, ":sha256"),
        paste0("validation:", group, ":", values$key)
      ),
      value = c(path, sha256_file(path), values$value), stringsAsFactors = FALSE
    ))
  }
  for (path in tools) {
    values <- read_tsv(path)
    require_columns(values, c("key", "value"), path)
    rows <- rbind(rows, values[, c("key", "value"), drop = FALSE])
  }
  for (path in sort(unique(artifacts))) {
    label <- sub("^results/remeta/export/", "", path)
    rows <- rbind(rows, data.frame(
      key = c(paste0("artifact:", label, ":path"), paste0("artifact:", label, ":bytes"), paste0("artifact:", label, ":sha256")),
      value = c(path, as.character(file.info(path)$size), sha256_file(path)), stringsAsFactors = FALSE
    ))
  }
  write_tsv(rows, out)
}


if (subtask == "prepare-target") {
  require_args(args, c("config", "pfile-prefix", "keep", "regions", "out-prefix", "summary-out", "threads"))
  prepare_target(load_config(args$config), args[["pfile-prefix"]], args$keep, args$regions,
    args[["out-prefix"]], args[["summary-out"]], args$threads)
} else if (subtask == "write-step2-command") {
  require_args(args, c("config", "trait", "group-summary", "pfile-prefix", "keep", "pheno", "covar", "pred-list",
    "trait-list", "covar-list", "out-prefix", "done", "script-out", "threads"))
  write_step2_command(load_config(args$config), args$trait, args[["group-summary"]], args[["pfile-prefix"]],
    args$keep, args$pheno, args$covar, args[["pred-list"]], args[["trait-list"]], args[["covar-list"]],
    args[["out-prefix"]], args$done, args[["script-out"]], args$threads)
} else if (subtask == "stage-trait") {
  require_args(args, c("trait", "group-summary", "raw-stats", "sample-ids", "keep", "out"))
  stage_trait(args$trait, args[["group-summary"]], args[["raw-stats"]], args[["sample-ids"]], args$keep, args$out)
} else if (subtask == "validate-group") {
  require_args(args, c("target-prefix", "keep", "ordinary-sample-ids", "sample-ids", "group-summary", "gene-list", "out"))
  if (!length(args$htp) || !length(args$ld)) die("validate-group requires --htp and --ld inputs")
  validate_group(args[["target-prefix"]], args$keep, args[["ordinary-sample-ids"]], args[["sample-ids"]],
    args[["group-summary"]], args[["gene-list"]], args$htp, args$ld, args$out)
} else if (subtask == "write-manifest") {
  require_args(args, c("config", "build", "gene-list", "provenance", "group-status", "out"))
  if (!length(args[["trait-summary"]]) || !length(args$tool)) {
    die("write-manifest requires trait summaries and tool records")
  }
  write_manifest(load_config(args$config), args$config, args$build, args[["gene-list"]], args$provenance, args[["group-status"]],
    args[["target-summary"]], args[["trait-summary"]], args$validation, args$tool, args$artifact, args$out)
} else {
  die("unknown ReMeta cohort subtask: ", subtask)
}
