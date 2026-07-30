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


raw <- commandArgs(trailingOnly = TRUE)
if (!length(raw)) die("missing ReMeta cohort subtask")
subtask <- raw[[1]]
args <- parse_args(
  defaults = list(threads = "1", artifact = character(), "target-summary" = character(),
    "trait-summary" = character(), validation = character(), tool = character(),
    htp = character(), index = character()),
  repeated = c("artifact", "target-summary", "trait-summary", "validation", "tool", "htp", "index"),
  raw = raw[-1]
)


read_tsv_no_metadata <- function(path) {
  if (!file.exists(path)) die("tab-delimited file not found: ", path)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^##", lines)]
  if (!length(lines) || !any(nzchar(trimws(lines)))) die("tab-delimited file is empty: ", path)
  read.table(
    text = paste(lines, collapse = "\n"), sep = "\t", header = TRUE,
    check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "",
    na.strings = character()
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


write_step2_command <- function(config, group_summary, pfile_prefix, pheno, covar, pred_list,
                                trait_list, covar_list, out_prefix, done, script_out, threads) {
  traits <- readLines(trait_list, warn = FALSE)
  traits <- traits[nzchar(traits)]
  if (!length(traits)) die("ReMeta group has no analyzable traits after Phase 2 thresholds")
  covars <- trimws(readLines(covar_list, warn = FALSE))
  covars <- covars[nzchar(covars)]
  if (length(covars) != 1) die("Phase 2 group covariate list must contain exactly one non-empty line")
  summary <- read_tsv(group_summary)
  active <- summary[summary$trait %in% traits, , drop = FALSE]
  types <- unique(active$trait_type)
  if (length(types) != 1 || !types %in% c("bt", "qt")) die("could not determine one Phase 2 trait type for ReMeta group")

  settings <- config$remeta
  phase2 <- config$phase2_regenie
  command <- c(
    "--step", "2",
    "--pgen", pfile_prefix,
    "--phenoFile", pheno,
    "--phenoColList", paste(traits, collapse = ","),
    "--covarFile", covar,
    "--covarColList", covars,
    "--pred", pred_list,
    "--htp", trimws(as.character(phase2$htp_cohort_name %||% "")),
    if (identical(types, "bt")) "--bt" else "--qt",
    "--bsize", as.character(phase2$step2_bsize %||% 400),
    "--minMAC", as.character(settings$min_mac %||% 1),
    "--threads", threads
  )
  cohort <- command[match("--htp", command) + 1]
  if (!nzchar(cohort)) command[match("--htp", command) + 1] <- analysis_output_name(config)
  if (identical(tolower(settings$genotype_mode), "dosage")) {
    command <- c(command, "--minINFO", as.character(settings$info_min %||% 0.8))
  }
  if (identical(types, "bt")) {
    command <- c(command, "--firth", "--approx", "--pThresh", as.character(phase2$p_thresh %||% 0.01))
  } else if (truthy(phase2$apply_rint %||% FALSE)) {
    command <- c(command, "--apply-rint")
  }
  command <- c(command, split_options(phase2$step2_options %||% ""), "--gz", "--out", out_prefix)
  expected <- paste0(out_prefix, "_", traits, ".regenie.gz")
  lines <- c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    paste("mkdir -p", shQuote(dirname(out_prefix), type = "sh")),
    shell_command_line(config$tools$regenie %||% "regenie", command),
    paste("test -s", vapply(expected, shQuote, character(1), type = "sh")),
    paste("printf '%s\\n' ok >", shQuote(done, type = "sh"))
  )
  write_bash_script(script_out, lines)
}


stage_trait <- function(trait, group_summary, raw_prefix, out) {
  summary <- read_tsv(group_summary)
  row <- summary[summary$trait == trait, , drop = FALSE]
  if (nrow(row) != 1) die("trait ", trait, " is absent or duplicated in Phase 2 group summary")
  ensure_parent(out)
  if (identical(row$skipped[[1]], "True")) {
    con <- gzfile(out, "wt")
    on.exit(close(con), add = TRUE)
    writeLines(c(paste0("## skipped: ", row$skip_reason[[1]]), paste(REMETA_HTP_COLUMNS, collapse = "\t")), con)
    return(invisible(TRUE))
  }
  raw <- paste0(raw_prefix, "_", trait, ".regenie.gz")
  if (!file.exists(raw) || file.info(raw)$size == 0) die("expected compressed regenie HTP output not found: ", raw)
  if (!file.copy(raw, out, overwrite = TRUE)) die("could not stage ReMeta regenie HTP output: ", out)
}


read_text_maybe_gzip <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  readLines(con, warn = FALSE)
}


validate_htp <- function(path, variant_ids, ld_variant_ids, chunk_size = 100000L) {
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


validate_group <- function(target_prefix, keep, gene_list, htp_paths, index_paths, out) {
  sample_count <- assert_identical_samples(keep, paste0(target_prefix, ".psam"))
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
  for (path in index_paths) {
    lines <- read_text_maybe_gzip(path)
    lines <- lines[nzchar(lines)]
    for (line in lines) {
      fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (length(fields) < 4) die("malformed ReMeta LD index row in ", path)
      if (!fields[[1]] %in% gene_ids) die("LD index contains gene absent from bundled gene list: ", fields[[1]])
      ids <- strsplit(fields[[4]], ",", fixed = TRUE)[[1]]
      ids <- ids[nzchar(ids)]
      if (!length(ids)) die("ReMeta LD index gene has no target variants: ", fields[[1]])
      missing <- setdiff(ids, variant_ids)
      if (length(missing)) die("LD index variants are absent from target PVAR: ", paste(head(missing, 5), collapse = ", "))
      gene_idx <- unname(gene_row[[fields[[1]]]])
      variant_idx <- unname(variant_row[ids])
      within_gene <- variants$CHROM_CLEAN[variant_idx] == gene_chrom[[gene_idx]] &
        variant_pos[variant_idx] >= gene_start[[gene_idx]] &
        variant_pos[variant_idx] <= gene_end[[gene_idx]]
      ld_variant_ids <- c(ld_variant_ids, ids)
      indexed_gene_ids <- c(indexed_gene_ids, fields[[1]])
      indexed_genes <- indexed_genes + 1L
      ld_gene_variant_assignments <- ld_gene_variant_assignments + length(ids)
      ld_assignments_within_gene_bounds <- ld_assignments_within_gene_bounds + sum(within_gene)
    }
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
    htp_variants <- htp_variants + validate_htp(path, variant_ids, ld_variant_ids)
  }
  write_key_values(c(
    status = "validated",
    sample_count = sample_count,
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


write_manifest <- function(config, config_path, build, gene_list, provenance, target_summaries,
                           trait_summaries, validations, tools, artifacts, out) {
  for (path in c(gene_list, provenance, target_summaries, trait_summaries, validations, tools, artifacts)) {
    if (!file.exists(path)) die("ReMeta manifest input not found: ", path)
  }
  metadata <- c(
    manifest_schema = "remeta_cohort_export_v1",
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
    resource_provenance_sha256 = sha256_file(provenance)
  )
  rows <- data.frame(key = names(metadata), value = unname(as.character(metadata)), stringsAsFactors = FALSE)

  for (path in target_summaries) {
    values <- read_tsv(path)
    require_columns(values, c("key", "value"), path)
    group <- basename(dirname(path))
    rows <- rbind(rows, data.frame(
      key = paste0("group:", group, ":", values$key), value = values$value, stringsAsFactors = FALSE
    ))
  }
  for (path in trait_summaries) {
    values <- read_tsv(path)
    require_columns(values, c("group", "trait", "trait_type", "phase2_pan_samples",
      "usable_n", "skipped", "skip_reason"), path)
    for (i in seq_len(nrow(values))) {
      prefix <- paste0("trait:", values$trait[[i]], ":")
      trait_metadata <- c(
        ld_group = values$group[[i]],
        trait_type = values$trait_type[[i]],
        phase2_pan_samples = values$phase2_pan_samples[[i]],
        usable_n_before_group_intersection = values$usable_n[[i]],
        skipped = values$skipped[[i]],
        skip_reason = values$skip_reason[[i]]
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
    rows <- rbind(rows, data.frame(
      key = paste0("validation:", group, ":", values$key), value = values$value, stringsAsFactors = FALSE
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
  require_args(args, c("config", "group-summary", "pfile-prefix", "pheno", "covar", "pred-list",
    "trait-list", "covar-list", "out-prefix", "done", "script-out", "threads"))
  write_step2_command(load_config(args$config), args[["group-summary"]], args[["pfile-prefix"]],
    args$pheno, args$covar, args[["pred-list"]], args[["trait-list"]], args[["covar-list"]],
    args[["out-prefix"]], args$done, args[["script-out"]], args$threads)
} else if (subtask == "stage-trait") {
  require_args(args, c("trait", "group-summary", "raw-prefix", "out"))
  stage_trait(args$trait, args[["group-summary"]], args[["raw-prefix"]], args$out)
} else if (subtask == "validate-group") {
  require_args(args, c("target-prefix", "keep", "gene-list", "out"))
  if (!length(args$htp) || !length(args$index)) die("validate-group requires --htp and --index inputs")
  validate_group(args[["target-prefix"]], args$keep, args[["gene-list"]], args$htp, args$index, args$out)
} else if (subtask == "write-manifest") {
  require_args(args, c("config", "build", "gene-list", "provenance", "out"))
  if (!length(args[["target-summary"]]) || !length(args[["trait-summary"]]) ||
      !length(args$validation) || !length(args$tool) || !length(args$artifact)) {
    die("write-manifest requires target and trait summaries, validations, tool records, and artifacts")
  }
  write_manifest(load_config(args$config), args$config, args$build, args[["gene-list"]], args$provenance,
    args[["target-summary"]], args[["trait-summary"]], args$validation, args$tool, args$artifact, args$out)
} else {
  die("unknown ReMeta cohort subtask: ", subtask)
}
