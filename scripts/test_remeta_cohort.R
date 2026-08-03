#!/usr/bin/env Rscript

# Focused tests for cohort ReMeta invariants and command generation.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))
remeta_script <- file.path(script_dir, "remeta_cohort.R")
rscript <- file.path(R.home("bin"), "Rscript")
tmp <- tempfile("test_remeta_cohort.")
dir.create(tmp, recursive = TRUE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
bgzip <- Sys.which("bgzip")
if (!nzchar(bgzip)) stop("test_remeta_cohort.R requires bgzip from HTSlib")


write_lines <- function(lines, name) {
  path <- file.path(tmp, name)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path)
  path
}


write_gzip <- function(lines, name) {
  path <- file.path(tmp, name)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, "wt")
  writeLines(lines, con)
  close(con)
  path
}


rewrite_bgzf_payload <- function(path, mutate) {
  raw_path <- paste0(path, ".logical")
  on.exit(unlink(c(raw_path, paste0(raw_path, ".gz"))), add = TRUE)
  status <- system2(bgzip, c("-d", "-c", path), stdout = raw_path, stderr = FALSE)
  if (!identical(as.integer(status), 0L)) stop("could not decompress BGZF fixture")
  con <- file(raw_path, "rb")
  bytes <- readBin(con, "raw", n = file.info(raw_path)$size)
  close(con)
  con <- file(raw_path, "wb")
  writeBin(mutate(bytes), con)
  close(con)
  output <- suppressWarnings(system2(bgzip, c("--binary", "-f", raw_path), stdout = TRUE, stderr = TRUE))
  compressed <- paste0(raw_path, ".gz")
  if (!is.null(attr(output, "status")) || !file.copy(compressed, path, overwrite = TRUE)) {
    stop("could not recompress BGZF fixture")
  }
  invisible(path)
}


make_ld_paths <- function(name, index_lines, index_chrom = 1L, sparse_entries = FALSE) {
  raw_uint_le <- function(value) sum(as.numeric(value) * 256^(seq_along(value) - 1L))
  bgzip_binary <- function(raw_path, path, logical_offsets) {
    output <- suppressWarnings(system2(bgzip, c("--binary", "--index", "-f", raw_path),
      stdout = TRUE, stderr = TRUE))
    compressed <- paste0(raw_path, ".gz")
    index <- paste0(compressed, ".gzi")
    if (!is.null(attr(output, "status")) ||
        !file.exists(index)) {
      stop("could not create BGZF LD fixture: ", paste(output, collapse = "\n"))
    }
    bytes <- readBin(index, "raw", file.info(index)$size)
    entry_count <- raw_uint_le(bytes[1:8])
    if (length(bytes) != 8 + 16 * entry_count) stop("malformed bgzip fixture index")
    compressed_offsets <- logical_starts <- numeric(entry_count)
    for (i in seq_len(entry_count)) {
      start <- 9 + (i - 1L) * 16L
      compressed_offsets[[i]] <- raw_uint_le(bytes[start:(start + 7L)])
      logical_starts[[i]] <- raw_uint_le(bytes[(start + 8L):(start + 15L)])
    }
    virtual_offsets <- vapply(logical_offsets, function(offset) {
      block <- max(which(c(0, logical_starts) <= offset))
      c(0, compressed_offsets)[[block]] * 65536 + offset - c(0, logical_starts)[[block]]
    }, numeric(1))
    if (!file.rename(compressed, path)) stop("could not stage BGZF LD fixture")
    unlink(index)
    virtual_offsets
  }
  parse_specs <- function(lines) {
    if (!length(lines)) return(list())
    lapply(lines, function(line) {
      fields <- strsplit(paste0(line, "\t__END__"), "\t", fixed = TRUE)[[1]]
      fields <- fields[-length(fields)]
      list(gene = fields[[1]], variants = strsplit(fields[[4]], ",", fixed = TRUE)[[1]])
    })
  }
  write_gene_ld <- function(path, specs) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    raw_path <- paste0(path, ".raw")
    con <- file(raw_path, "wb")
    writeBin(charToRaw("remetaLD.v1.1"), con)
    offsets <- numeric(length(specs))
    for (i in seq_along(specs)) {
      offsets[[i]] <- seek(con)
      variant_count <- length(specs[[i]]$variants)
      entry_count <- if (sparse_entries && variant_count) 1L else 0L
      writeBin(as.integer(c(variant_count, entry_count)), con, size = 4L, endian = "little")
      writeBin(0.0001, con, size = 4L, endian = "little")
      writeBin(rep(1, variant_count), con, size = 4L, endian = "little")
      if (entry_count) {
        writeBin(as.integer(c(0L, 0L)), con, size = 4L, endian = "little")
        writeBin(1, con, size = 4L, endian = "little")
      }
      writeBin(as.integer(c(variant_count, 0L, 0L)), con, size = 4L, endian = "little")
      writeBin(0.0001, con, size = 4L, endian = "little")
      writeBin(0L, con, size = 4L, endian = "little")
    }
    close(con)
    bgzip_binary(raw_path, path, offsets)
  }
  write_buffer_ld <- function(path, specs) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    raw_path <- paste0(path, ".raw")
    con <- file(raw_path, "wb")
    writeBin(charToRaw("remetaLD.v1.1"), con)
    writeBin(4L, con, size = 4L, endian = "little")
    offsets <- numeric(length(specs))
    for (i in seq_along(specs)) {
      writeBin(-3L, con, size = 4L, endian = "little")
      offsets[[i]] <- seek(con)
      writeBin(-1L, con, size = 4L, endian = "little")
    }
    close(con)
    bgzip_binary(raw_path, path, offsets)
  }

  paths <- character()
  for (chrom in 1:22) {
    specs <- if (chrom == index_chrom) parse_specs(index_lines) else list()
    gene_path <- file.path(tmp, name, paste0("chr", chrom, ".remeta.gene.ld"))
    buffer_path <- file.path(tmp, name, paste0("chr", chrom, ".remeta.buffer.ld"))
    gene_offsets <- write_gene_ld(gene_path, specs)
    buffer_offsets <- write_buffer_ld(buffer_path, specs)
    rewritten_index <- vapply(seq_along(specs), function(i) paste(
      c(specs[[i]]$gene, gene_offsets[[i]], buffer_offsets[[i]],
        paste(specs[[i]]$variants, collapse = ","), ""),
      collapse = "\t"
    ), character(1))
    paths <- c(
      paths,
      gene_path,
      buffer_path,
      write_gzip(rewritten_index,
        file.path(name, paste0("chr", chrom, ".remeta.ld.idx.gz")))
    )
  }
  paths
}


run_task <- function(values) {
  system2(rscript, c(remeta_script, values), stdout = TRUE, stderr = TRUE)
}


htp_columns <- c(
  "Name", "Chr", "Pos", "Ref", "Alt", "Trait", "Cohort", "Model", "Effect", "LCI_Effect",
  "UCI_Effect", "Pval", "AAF", "Num_Cases", "Cases_Ref", "Cases_Het", "Cases_Alt",
  "Num_Controls", "Controls_Ref", "Controls_Het", "Controls_Alt", "Info"
)


htp_row <- function(name, chrom, pos, ref, alt, info, num_cases = 1L, case_genotypes = c(0L, 1L, 0L),
                    num_controls = 1L, control_genotypes = c(0L, 1L, 0L)) {
  paste(c(name, chrom, pos, ref, alt, "TRAIT1", "cohort_test", "ADD", "1.1", "1.0", "1.2",
    "0.5", "0.1", num_cases, case_genotypes, num_controls, control_genotypes, info), collapse = "\t")
}


target <- file.path(tmp, "target")
invisible(write_lines(c("#FID\tIID", "F1\tI1", "F2\tI2"), "target.psam"))
invisible(write_lines(c(
  "#CHROM\tPOS\tID\tREF\tALT",
  "1\t100\t1:100:A:G\tA\tG",
  "1\t120\t1:120:C:T\tC\tT"
), "target.pvar"))

# Terminal genes must never see the first variant on the next chromosome.
# The LD reader therefore receives an exact, nonempty target-ID partition.
partition_target <- file.path(tmp, "partition_target")
partition_rows <- c(
  "#CHROM\tPOS\tID\tREF\tALT",
  "1\t100\t1:100:A:G\tA\tG",
  "1\t200\t1:200:C:T\tC\tT",
  vapply(2:22, function(chrom) {
    paste(chrom, 100, paste(chrom, 100, "A", "G", sep = ":"), "A", "G", sep = "\t")
  }, character(1))
)
invisible(write_lines(partition_rows, "partition_target.pvar"))
extracts <- file.path(tmp, "extracts", paste0("chr", 22:1, ".target_variants.txt"))
output <- run_task(c(
  "write-ld-extracts", "--target-prefix", partition_target, "--extract-out", extracts
))
stopifnot(is.null(attr(output, "status")), all(file.exists(extracts)))
stopifnot(identical(readLines(file.path(tmp, "extracts/chr1.target_variants.txt")),
  c("1:100:A:G", "1:200:C:T")))
stopifnot(identical(readLines(file.path(tmp, "extracts/chr22.target_variants.txt")), "22:100:A:G"))

# An empty autosome is not silently represented by fabricated LD artifacts.
incomplete_target <- file.path(tmp, "incomplete_partition_target")
invisible(write_lines(partition_rows[!grepl("^22\\t", partition_rows)],
  "incomplete_partition_target.pvar"))
incomplete_extracts <- file.path(tmp, "incomplete_extracts",
  paste0("chr", 1:22, ".target_variants.txt"))
incomplete <- suppressWarnings(run_task(c(
  "write-ld-extracts", "--target-prefix", incomplete_target,
  "--extract-out", incomplete_extracts
)))
stopifnot(!is.null(attr(incomplete, "status")), attr(incomplete, "status") != 0,
  !any(file.exists(incomplete_extracts)))

invisible(write_lines(c("F1\tI1", "F2\tI2"), "keep.txt"))
keep_sha <- sha256_file(file.path(tmp, "keep.txt"))
ordinary_ids <- write_lines(c("F2\tI2", "F1\tI1"), "ordinary.regenie.ids")
rare_ids <- write_lines(c("F1\tI1", "F2\tI2"), "rare.regenie.ids")
group_summary <- write_lines(c(
  paste(c(
    "group", "trait", "trait_type", "phase2_pan_samples", "usable_n", "cases", "controls",
    "model_sample_count", "model_cases", "model_controls", "model_keep_sha256", "skipped", "skip_reason"
  ), collapse = "\t"),
  paste(c("bt__TRAIT1", "TRAIT1", "bt", "2", "2", "1", "1", "2", "1", "1",
    keep_sha, "False", ""), collapse = "\t")
), "bt__TRAIT1/summary.tsv")
gene_list <- write_lines(c(
  "ENSG000001\t1\t90\t130",
  "ENSG000002\t1\t200\t250"
), "gene_list.tsv")
htp <- write_gzip(c(
  paste(htp_columns, collapse = "\t"),
  htp_row("1:100:A:G", "1", "100", "A", "G", "SCORE=0.2;SKATV=0.1")
), "trait.regenie.gz")
ld_paths <- make_ld_paths(
  "ld", "ENSG000001\t0\t0\t1:100:A:G,1:120:C:T\t", sparse_entries = TRUE
)
ok <- file.path(tmp, "validation.ok")

output <- run_task(c(
  "validate-group", "--target-prefix", target, "--keep", file.path(tmp, "keep.txt"),
  "--ordinary-sample-ids", ordinary_ids, "--sample-ids", rare_ids,
  "--group-summary", group_summary, "--gene-list", gene_list, "--htp", htp,
  "--ld", ld_paths, "--out", ok
))
stopifnot(is.null(attr(output, "status")), file.exists(ok))
validation <- read.delim(ok, stringsAsFactors = FALSE)
stopifnot(validation$value[validation$key == "sample_count"] == "2")
stopifnot(validation$value[validation$key == "target_variant_count"] == "2")
stopifnot(validation$value[validation$key == "unique_ld_target_variant_count"] == "2")
stopifnot(validation$value[validation$key == "target_variant_ld_coverage_pct"] == "100.000000")
stopifnot(validation$value[validation$key == "reference_gene_count"] == "2")
stopifnot(validation$value[validation$key == "indexed_gene_count"] == "1")
stopifnot(validation$value[validation$key == "indexed_gene_coverage_pct"] == "50.000000")
stopifnot(validation$value[validation$key == "ld_gene_variant_assignments"] == "2")
stopifnot(validation$value[validation$key == "ld_assignments_within_gene_bounds"] == "2")
stopifnot(validation$value[validation$key == "ld_assignment_gene_bound_coverage_pct"] == "100.000000")

# The second record begins beyond one 64-KiB BGZF block, so its index address is
# a real HTSlib virtual offset rather than an ordinary byte position.
multiblock_count <- 16400L
multiblock_ids <- paste("1", seq_len(multiblock_count + 1L), "A", "G", sep = ":")
multiblock_target <- file.path(tmp, "multiblock_target")
invisible(write_lines(c("#FID\tIID", "F1\tI1", "F2\tI2"), "multiblock_target.psam"))
invisible(write_lines(c(
  "#CHROM\tPOS\tID\tREF\tALT",
  paste("1", seq_len(multiblock_count + 1L), multiblock_ids, "A", "G", sep = "\t")
), "multiblock_target.pvar"))
multiblock_genes <- write_lines(c(
  paste("ENSG_LARGE", 1, 1, multiblock_count, sep = "\t"),
  paste("ENSG_SECOND", 1, multiblock_count + 1L, multiblock_count + 1L, sep = "\t")
), "multiblock_genes.tsv")
multiblock_ld <- make_ld_paths("multiblock_ld", c(
  paste0("ENSG_LARGE\t0\t0\t", paste(multiblock_ids[seq_len(multiblock_count)], collapse = ","), "\t"),
  paste0("ENSG_SECOND\t0\t0\t", multiblock_ids[[multiblock_count + 1L]], "\t")
))
multiblock_index <- multiblock_ld[grepl("chr1\\.remeta\\.ld\\.idx\\.gz$", multiblock_ld)][[1]]
multiblock_index_lines <- readLines(gzfile(multiblock_index), warn = FALSE)
second_offset <- as.numeric(strsplit(multiblock_index_lines[[2]], "\t", fixed = TRUE)[[1]][[2]])
stopifnot(second_offset > 65536)
multiblock_htp <- write_gzip(c(
  paste(htp_columns, collapse = "\t"),
  htp_row(multiblock_ids[[1]], "1", "1", "A", "G", "SCORE=0.2;SKATV=0.1")
), "multiblock.regenie.gz")
multiblock <- run_task(c(
  "validate-group", "--target-prefix", multiblock_target, "--keep", file.path(tmp, "keep.txt"),
  "--ordinary-sample-ids", ordinary_ids, "--sample-ids", rare_ids,
  "--group-summary", group_summary, "--gene-list", multiblock_genes, "--htp", multiblock_htp,
  "--ld", multiblock_ld, "--out", file.path(tmp, "multiblock.ok")
))
stopifnot(is.null(attr(multiblock, "status")))

# HTP counts are variant-level nonmissing genotype counts, not repeated model Ns.
missing_genotypes_htp <- write_gzip(c(
  paste(htp_columns, collapse = "\t"),
  htp_row("1:100:A:G", "1", "100", "A", "G", "SCORE=0.2;SKATV=0.1",
    num_cases = 0L, case_genotypes = c(0L, 0L, 0L))
), "missing_genotypes.regenie.gz")
missing_genotypes <- run_task(c(
  "validate-group", "--target-prefix", target, "--keep", file.path(tmp, "keep.txt"),
  "--ordinary-sample-ids", ordinary_ids, "--sample-ids", rare_ids,
  "--group-summary", group_summary, "--gene-list", gene_list, "--htp", missing_genotypes_htp,
  "--ld", ld_paths, "--out", file.path(tmp, "missing_genotypes.ok")
))
stopifnot(is.null(attr(missing_genotypes, "status")))

bad_genotype_sum_htp <- write_gzip(c(
  paste(htp_columns, collapse = "\t"),
  htp_row("1:100:A:G", "1", "100", "A", "G", "SCORE=0.2;SKATV=0.1",
    num_cases = 1L, case_genotypes = c(0L, 0L, 0L))
), "bad_genotype_sum.regenie.gz")
bad_genotype_sum <- suppressWarnings(run_task(c(
  "validate-group", "--target-prefix", target, "--keep", file.path(tmp, "keep.txt"),
  "--ordinary-sample-ids", ordinary_ids, "--sample-ids", rare_ids,
  "--group-summary", group_summary, "--gene-list", gene_list, "--htp", bad_genotype_sum_htp,
  "--ld", ld_paths, "--out", file.path(tmp, "bad_genotype_sum.ok")
)))
stopifnot(!is.null(attr(bad_genotype_sum, "status")), attr(bad_genotype_sum, "status") != 0)

excess_genotypes_htp <- write_gzip(c(
  paste(htp_columns, collapse = "\t"),
  htp_row("1:100:A:G", "1", "100", "A", "G", "SCORE=0.2;SKATV=0.1",
    num_cases = 2L, case_genotypes = c(1L, 1L, 0L))
), "excess_genotypes.regenie.gz")
excess_genotypes <- suppressWarnings(run_task(c(
  "validate-group", "--target-prefix", target, "--keep", file.path(tmp, "keep.txt"),
  "--ordinary-sample-ids", ordinary_ids, "--sample-ids", rare_ids,
  "--group-summary", group_summary, "--gene-list", gene_list, "--htp", excess_genotypes_htp,
  "--ld", ld_paths, "--out", file.path(tmp, "excess_genotypes.ok")
)))
stopifnot(!is.null(attr(excess_genotypes, "status")), attr(excess_genotypes, "status") != 0)

bad_ids <- write_lines("F1\tI1", "bad.regenie.ids")
bad_samples <- suppressWarnings(run_task(c(
  "validate-group", "--target-prefix", target, "--keep", file.path(tmp, "keep.txt"),
  "--ordinary-sample-ids", ordinary_ids, "--sample-ids", bad_ids,
  "--group-summary", group_summary, "--gene-list", gene_list, "--htp", htp,
  "--ld", ld_paths, "--out", file.path(tmp, "bad_samples.ok")
)))
stopifnot(!is.null(attr(bad_samples, "status")), attr(bad_samples, "status") != 0)

wrong_trait_htp <- write_gzip(c(
  paste(htp_columns, collapse = "\t"),
  sub("\tTRAIT1\t", "\tOTHER\t", htp_row("1:100:A:G", "1", "100", "A", "G", "SCORE=0.2;SKATV=0.1"))
), "wrong_trait.regenie.gz")
wrong_trait <- suppressWarnings(run_task(c(
  "validate-group", "--target-prefix", target, "--keep", file.path(tmp, "keep.txt"),
  "--ordinary-sample-ids", ordinary_ids, "--sample-ids", rare_ids,
  "--group-summary", group_summary, "--gene-list", gene_list, "--htp", wrong_trait_htp,
  "--ld", ld_paths, "--out", file.path(tmp, "wrong_trait.ok")
)))
stopifnot(!is.null(attr(wrong_trait, "status")), attr(wrong_trait, "status") != 0)

incomplete_ld <- suppressWarnings(run_task(c(
  "validate-group", "--target-prefix", target, "--keep", file.path(tmp, "keep.txt"),
  "--ordinary-sample-ids", ordinary_ids, "--sample-ids", rare_ids,
  "--group-summary", group_summary, "--gene-list", gene_list, "--htp", htp,
  "--ld", ld_paths[-1], "--out", file.path(tmp, "incomplete_ld.ok")
)))
stopifnot(!is.null(attr(incomplete_ld, "status")), attr(incomplete_ld, "status") != 0)

duplicate_ld_paths <- ld_paths
duplicate_ld_paths[[2]] <- duplicate_ld_paths[[1]]
duplicate_ld <- suppressWarnings(run_task(c(
  "validate-group", "--target-prefix", target, "--keep", file.path(tmp, "keep.txt"),
  "--ordinary-sample-ids", ordinary_ids, "--sample-ids", rare_ids,
  "--group-summary", group_summary, "--gene-list", gene_list, "--htp", htp,
  "--ld", duplicate_ld_paths, "--out", file.path(tmp, "duplicate_ld.ok")
)))
stopifnot(!is.null(attr(duplicate_ld, "status")), attr(duplicate_ld, "status") != 0)

bad_header_ld <- make_ld_paths("bad_header_ld", "ENSG000001\t0\t0\t1:100:A:G,1:120:C:T\t")
bad_header_path <- bad_header_ld[grepl("chr1\\.remeta\\.gene\\.ld$", bad_header_ld)][[1]]
rewrite_bgzf_payload(bad_header_path, function(bytes) {
  bytes[[1]] <- charToRaw("X")
  bytes
})
bad_header <- suppressWarnings(run_task(c(
  "validate-group", "--target-prefix", target, "--keep", file.path(tmp, "keep.txt"),
  "--ordinary-sample-ids", ordinary_ids, "--sample-ids", rare_ids,
  "--group-summary", group_summary, "--gene-list", gene_list, "--htp", htp,
  "--ld", bad_header_ld, "--out", file.path(tmp, "bad_header.ok")
)))
stopifnot(!is.null(attr(bad_header, "status")), attr(bad_header, "status") != 0)

bad_offset_ld <- make_ld_paths("bad_offset_ld", "ENSG000001\t0\t0\t1:100:A:G,1:120:C:T\t")
bad_offset_index <- bad_offset_ld[grepl("chr1\\.remeta\\.ld\\.idx\\.gz$", bad_offset_ld)][[1]]
bad_offset_con <- gzfile(bad_offset_index, "wt")
writeLines("ENSG000001\t9999\t21\t1:100:A:G,1:120:C:T\t", bad_offset_con)
close(bad_offset_con)
bad_offset <- suppressWarnings(run_task(c(
  "validate-group", "--target-prefix", target, "--keep", file.path(tmp, "keep.txt"),
  "--ordinary-sample-ids", ordinary_ids, "--sample-ids", rare_ids,
  "--group-summary", group_summary, "--gene-list", gene_list, "--htp", htp,
  "--ld", bad_offset_ld, "--out", file.path(tmp, "bad_offset.ok")
)))
stopifnot(!is.null(attr(bad_offset, "status")), attr(bad_offset, "status") != 0)

truncated_buffer_ld <- make_ld_paths("truncated_buffer_ld", "ENSG000001\t0\t0\t1:100:A:G,1:120:C:T\t")
truncated_buffer_path <- truncated_buffer_ld[grepl("chr1\\.remeta\\.buffer\\.ld$", truncated_buffer_ld)][[1]]
rewrite_bgzf_payload(truncated_buffer_path, function(bytes) bytes[seq_len(17L)])
truncated_buffer <- suppressWarnings(run_task(c(
  "validate-group", "--target-prefix", target, "--keep", file.path(tmp, "keep.txt"),
  "--ordinary-sample-ids", ordinary_ids, "--sample-ids", rare_ids,
  "--group-summary", group_summary, "--gene-list", gene_list, "--htp", htp,
  "--ld", truncated_buffer_ld, "--out", file.path(tmp, "truncated_buffer.ok")
)))
stopifnot(!is.null(attr(truncated_buffer, "status")), attr(truncated_buffer, "status") != 0)

bad_marker_ld <- make_ld_paths("bad_marker_ld", "ENSG000001\t0\t0\t1:100:A:G,1:120:C:T\t")
bad_marker_path <- bad_marker_ld[grepl("chr1\\.remeta\\.buffer\\.ld$", bad_marker_ld)][[1]]
rewrite_bgzf_payload(bad_marker_path, function(bytes) {
  bytes[18:21] <- as.raw(c(254, 255, 255, 255))
  bytes
})
bad_marker <- suppressWarnings(run_task(c(
  "validate-group", "--target-prefix", target, "--keep", file.path(tmp, "keep.txt"),
  "--ordinary-sample-ids", ordinary_ids, "--sample-ids", rare_ids,
  "--group-summary", group_summary, "--gene-list", gene_list, "--htp", htp,
  "--ld", bad_marker_ld, "--out", file.path(tmp, "bad_marker.ok")
)))
stopifnot(!is.null(attr(bad_marker, "status")), attr(bad_marker, "status") != 0)

wrong_chrom_genes <- write_lines("ENSG_WRONG\t2\t90\t130", "wrong_chrom_gene_list.tsv")
wrong_chrom_ld <- make_ld_paths(
  "wrong_chrom_ld", "ENSG_WRONG\t0\t0\t1:100:A:G\t", index_chrom = 2L
)
wrong_chrom <- suppressWarnings(run_task(c(
  "validate-group", "--target-prefix", target, "--keep", file.path(tmp, "keep.txt"),
  "--ordinary-sample-ids", ordinary_ids, "--sample-ids", rare_ids,
  "--group-summary", group_summary, "--gene-list", wrong_chrom_genes, "--htp", htp,
  "--ld", wrong_chrom_ld, "--out", file.path(tmp, "wrong_chrom.ok")
)))
stopifnot(!is.null(attr(wrong_chrom, "status")), attr(wrong_chrom, "status") != 0)


# A score-statistic variant without matching LD target data must fail.
bad_htp <- write_gzip(c(
  paste(htp_columns, collapse = "\t"),
  htp_row("1:999:A:G", "1", "999", "A", "G", "SCORE=0.2;SKATV=0.1")
), "bad.regenie.gz")
bad <- suppressWarnings(run_task(c(
  "validate-group", "--target-prefix", target, "--keep", file.path(tmp, "keep.txt"),
  "--ordinary-sample-ids", ordinary_ids, "--sample-ids", rare_ids,
  "--group-summary", group_summary, "--gene-list", gene_list, "--htp", bad_htp, "--ld", ld_paths,
  "--out", file.path(tmp, "bad.ok")
)))
stopifnot(!is.null(attr(bad, "status")), attr(bad, "status") != 0)

# A target-PVAR variant without an LD-index entry must also fail.
missing_ld_htp <- write_gzip(c(
  paste(htp_columns, collapse = "\t"),
  htp_row("1:120:C:T", "1", "120", "C", "T", "SCORE=0.2;SKATV=0.1")
), "missing_ld.regenie.gz")
ld_paths_one <- make_ld_paths("ld_one", "ENSG000001\t0\t0\t1:100:A:G\t")
missing_ld <- suppressWarnings(run_task(c(
  "validate-group", "--target-prefix", target, "--keep", file.path(tmp, "keep.txt"),
  "--ordinary-sample-ids", ordinary_ids, "--sample-ids", rare_ids,
  "--group-summary", group_summary, "--gene-list", gene_list, "--htp", missing_ld_htp, "--ld", ld_paths_one,
  "--out", file.path(tmp, "missing_ld.ok")
)))
stopifnot(!is.null(attr(missing_ld, "status")), attr(missing_ld, "status") != 0)

# A CPRA identifier is insufficient when PLINK still marks REF as provisional.
provisional <- file.path(tmp, "provisional")
invisible(file.copy(file.path(tmp, "target.psam"), paste0(provisional, ".psam")))
invisible(write_lines(c(
  "#CHROM\tPOS\tID\tREF\tALT\tINFO",
  "1\t100\t1:100:A:G\tA\tG\tPR"
), "provisional.pvar"))
bad_ref <- suppressWarnings(run_task(c(
  "validate-group", "--target-prefix", provisional, "--keep", file.path(tmp, "keep.txt"),
  "--ordinary-sample-ids", ordinary_ids, "--sample-ids", rare_ids,
  "--group-summary", group_summary, "--gene-list", gene_list, "--htp", htp, "--ld", ld_paths,
  "--out", file.path(tmp, "bad_ref.ok")
)))
stopifnot(!is.null(attr(bad_ref, "status")), attr(bad_ref, "status") != 0)


config <- write_lines(c(
  "project:",
  "  analysis_name: cohort_test",
  "  cohort_data_release: test_release",
  "phase2_regenie:",
  "  htp_cohort_name: cohort_test",
  "  step2_bsize: 400",
  "  p_thresh: 0.01",
  "  apply_rint: false",
  "  step2_options: ''",
  "remeta:",
  "  data_source: wes",
  "  genotype_mode: hardcall",
  "  input_variants_normalized: true",
  "  min_mac: 1",
  "  info_min: 0.8",
  "  target_r2: 0.0001",
  "tools:",
  "  regenie: regenie"
), "config.yaml")
summary <- write_lines(c(
  paste(c("group", "trait", "trait_type", "phase2_pan_samples", "usable_n", "cases", "controls",
    "model_sample_count", "model_cases", "model_controls", "model_keep_sha256", "skipped", "skip_reason"), collapse = "\t"),
  "bt__TRAIT1\tTRAIT1\tbt\t100\t2\t1\t1\t2\t1\t1\ttesthash\tFalse\t"
), "summary.tsv")
trait_list <- write_lines("TRAIT1", "traits.txt")
covar_list <- write_lines("age,sex,PC1", "covars.txt")
command <- file.path(tmp, "step2.sh")
done <- file.path(tmp, "step2.done")
output <- run_task(c(
  "write-step2-command", "--config", config, "--trait", "TRAIT1", "--group-summary", summary,
  "--pfile-prefix", target, "--pheno", "pheno.tsv", "--covar", "covar.tsv",
  "--keep", file.path(tmp, "keep.txt"),
  "--pred-list", "step1_pred.list", "--trait-list", trait_list,
  "--covar-list", covar_list, "--out-prefix", file.path(tmp, "rare"),
  "--done", done, "--script-out", command, "--threads", "4"
))
stopifnot(is.null(attr(output, "status")), file.exists(command))
text <- paste(readLines(command), collapse = "\n")
stopifnot(grepl("--gz", text, fixed = TRUE))
stopifnot(grepl("--minMAC.*'1'", text))
stopifnot(grepl("--htp.*'cohort_test'", text))
stopifnot(grepl("--write-samples", text, fixed = TRUE))
stopifnot(grepl("--minCaseCount.*'10'", text))
stopifnot(grepl("--keep", text, fixed = TRUE))
stopifnot(!grepl("--minINFO", text, fixed = TRUE))

staged_htp <- file.path(tmp, "staged.regenie.gz")
output <- run_task(c(
  "stage-trait", "--trait", "TRAIT1", "--group-summary", summary,
  "--raw-stats", htp, "--sample-ids", rare_ids, "--keep", file.path(tmp, "keep.txt"),
  "--out", staged_htp
))
stopifnot(is.null(attr(output, "status")), identical(readBin(htp, "raw", file.info(htp)$size),
  readBin(staged_htp, "raw", file.info(staged_htp)$size)))

skipped_summary <- write_lines(c(
  paste(c("group", "trait", "trait_type", "phase2_pan_samples", "usable_n", "cases", "controls",
    "model_sample_count", "model_cases", "model_controls", "model_keep_sha256", "skipped", "skip_reason"), collapse = "\t"),
  "bt__SKIPPED\tSKIPPED\tbt\t100\t10\t5\t5\t0\t\t\ttesthash\tTrue\tbelow_threshold"
), "skipped_summary.tsv")
skipped_htp <- file.path(tmp, "skipped.regenie.gz")
output <- suppressWarnings(run_task(c(
  "stage-trait", "--trait", "SKIPPED", "--group-summary", skipped_summary,
  "--raw-stats", htp, "--sample-ids", rare_ids, "--keep", file.path(tmp, "keep.txt"),
  "--out", skipped_htp
)))
stopifnot(!is.null(attr(output, "status")), attr(output, "status") != 0, !file.exists(skipped_htp))

dosage_step2_config <- write_lines(c(
  readLines(config),
  ""
), "dosage_step2_config.yaml")
dosage_lines <- readLines(dosage_step2_config)
dosage_lines[dosage_lines == "  data_source: wes"] <- "  data_source: imputed"
dosage_lines[dosage_lines == "  genotype_mode: hardcall"] <- "  genotype_mode: dosage"
writeLines(dosage_lines, dosage_step2_config)
dosage_command <- file.path(tmp, "dosage_step2.sh")
output <- run_task(c(
  "write-step2-command", "--config", dosage_step2_config, "--trait", "TRAIT1", "--group-summary", summary,
  "--pfile-prefix", target, "--pheno", "pheno.tsv", "--covar", "covar.tsv",
  "--keep", file.path(tmp, "keep.txt"),
  "--pred-list", "step1_pred.list", "--trait-list", trait_list,
  "--covar-list", covar_list, "--out-prefix", file.path(tmp, "rare_dosage"),
  "--done", file.path(tmp, "dosage.done"), "--script-out", dosage_command, "--threads", "4"
))
stopifnot(is.null(attr(output, "status")))
stopifnot(grepl("--minINFO.*'0.8'", paste(readLines(dosage_command), collapse = "\n")))


# Capture the target-preparation PLINK arguments and return a valid tiny PGEN
# trio. This checks the coordinate convention and literal CPRA template.
captured <- file.path(tmp, "plink_args.txt")
fake_plink <- write_lines(c(
  "#!/usr/bin/env bash",
  "set -euo pipefail",
  paste("printf '%s\\n' \"$@\" >", shQuote(captured)),
  "out=''",
  "while (($#)); do",
  "  if [[ \"$1\" == '--out' ]]; then out=\"$2\"; shift 2; else shift; fi",
  "done",
  "printf '#FID\\tIID\\nF1\\tI1\\nF2\\tI2\\n' > \"${out}.psam\"",
  "printf '#CHROM\\tPOS\\tID\\tREF\\tALT\\n1\\t100\\t1:100:A:G\\tA\\tG\\n' > \"${out}.pvar\"",
  ": > \"${out}.pgen\""
), "fake_plink.sh")
Sys.chmod(fake_plink, "0755")
prepare_config <- write_lines(c(
  "remeta:",
  "  data_source: wes",
  "  genotype_mode: hardcall",
  "  input_variants_normalized: true",
  "  min_mac: 1",
  "  geno_missing_max: 0.05",
  "tools:",
  paste0("  plink2: ", fake_plink)
), "prepare_config.yaml")
regions <- write_lines("1\t90\t130", "regions.bed")
prepared <- file.path(tmp, "prepared")
output <- run_task(c(
  "prepare-target", "--config", prepare_config, "--pfile-prefix", "source",
  "--keep", file.path(tmp, "keep.txt"), "--regions", regions,
  "--out-prefix", prepared, "--summary-out", file.path(tmp, "prepared.summary.tsv"),
  "--threads", "2"
))
stopifnot(is.null(attr(output, "status")))
plink_args <- readLines(captured)
stopifnot(any(plink_args == "bed0"))
stopifnot(any(plink_args == "26"))
stopifnot(any(plink_args == "@:#:$r:$a"))
stopifnot(any(plink_args == "erase-dosage"))

dosage_config <- write_lines(c(
  "remeta:",
  "  data_source: imputed",
  "  genotype_mode: dosage",
  "  input_variants_normalized: true",
  "  min_mac: 1",
  "  geno_missing_max: 0.05",
  "  info_min: 0.8",
  "tools:",
  paste0("  plink2: ", fake_plink)
), "dosage_config.yaml")
output <- run_task(c(
  "prepare-target", "--config", dosage_config, "--pfile-prefix", "source",
  "--keep", file.path(tmp, "keep.txt"), "--regions", regions,
  "--out-prefix", file.path(tmp, "prepared_dosage"),
  "--summary-out", file.path(tmp, "prepared_dosage.summary.tsv"), "--threads", "2"
))
stopifnot(is.null(attr(output, "status")))
plink_args <- readLines(captured)
stopifnot(any(plink_args == "dosage"))
stopifnot(any(plink_args == "--mach-r2-filter"))
stopifnot(!any(plink_args == "erase-dosage"))


provenance <- write_lines(c("key\tvalue", "genome_build\tGRCh38"), "provenance.tsv")
tool <- write_lines(c("key\tvalue", "tool_remeta_version\t0.11.2"), "tool.tsv")
status <- write_lines(c(
  "group\ttrait\ttrait_type\tusable_n\tcases\tcontrols\tmodel_sample_count\tmodel_cases\tmodel_controls\tkeep_count\tmodel_keep_sha256\tskipped\tskip_reason\tremeta_eligible",
  paste(c("bt__TRAIT1", "TRAIT1", "bt", 2, 1, 1, 2, 1, 1, 2,
    keep_sha, "False", "", "True"), collapse = "\t")
), "group_status.tsv")
manifest_target <- write_lines(readLines(file.path(tmp, "prepared.summary.tsv")),
  "manifest_work/bt__TRAIT1/target.summary.tsv")
manifest_validation <- write_lines(readLines(ok), "manifest_work/bt__TRAIT1/validation.ok")
manifest_htp <- write_gzip(c(paste(htp_columns, collapse = "\t"),
  htp_row("1:100:A:G", "1", "100", "A", "G", "SCORE=0.2;SKATV=0.1")),
  "export/GRCh38/htp/TRAIT1.PAN.regenie.gz")
manifest_ld <- make_ld_paths(
  "export/GRCh38/ld/bt__TRAIT1", "ENSG000001\t0\t0\t1:100:A:G,1:120:C:T\t"
)
manifest <- file.path(tmp, "manifest.tsv")
output <- run_task(c(
  "write-manifest", "--config", config, "--build", "GRCh38",
  "--gene-list", gene_list, "--provenance", provenance, "--group-status", status,
  "--target-summary", manifest_target,
  "--trait-summary", group_summary,
  "--validation", manifest_validation, "--tool", tool, "--artifact", manifest_htp, manifest_ld,
  "--out", manifest
))
stopifnot(is.null(attr(output, "status")), file.exists(manifest))
manifest_rows <- read.delim(manifest, stringsAsFactors = FALSE)
stopifnot(manifest_rows$value[manifest_rows$key == "analysis_scope"] == "marginal_gene_tests_only")
stopifnot(manifest_rows$value[manifest_rows$key == "conditional_buffer_included"] == "false")
stopifnot(manifest_rows$value[manifest_rows$key == "manifest_schema"] == "remeta_cohort_export_v2")
stopifnot(manifest_rows$value[manifest_rows$key == "ld_target_partitioning"] ==
  "chromosome_specific_extract")
stopifnot(manifest_rows$value[manifest_rows$key == "trait:TRAIT1:ld_group"] == "bt__TRAIT1")
stopifnot(manifest_rows$value[manifest_rows$key == "trait:TRAIT1:ld_prefix"] ==
  "results/remeta/export/GRCh38/ld/bt__TRAIT1/chr{1-22}")
stopifnot(manifest_rows$value[manifest_rows$key == "trait:TRAIT1:remeta_exported"] == "true")
stopifnot(manifest_rows$value[manifest_rows$key == "validation:bt__TRAIT1:sha256"] ==
  sha256_file(manifest_validation))
stopifnot(any(grepl(":sha256$", manifest_rows$key)))

stale_export <- write_lines("stale", "export/GRCh38/ld/bt_g1/stale.txt")
stale_manifest <- suppressWarnings(run_task(c(
  "write-manifest", "--config", config, "--build", "GRCh38",
  "--gene-list", gene_list, "--provenance", provenance, "--group-status", status,
  "--target-summary", manifest_target, "--trait-summary", group_summary,
  "--validation", manifest_validation, "--tool", tool, "--artifact", manifest_htp, manifest_ld,
  "--out", file.path(tmp, "stale_manifest.tsv")
)))
stopifnot(!is.null(attr(stale_manifest, "status")), attr(stale_manifest, "status") != 0)
unlink(stale_export)

# ReMeta creates this native log beside --out; it is runtime provenance, not
# a central handoff artifact, and must never survive in the export directory.
native_log <- write_lines(
  "native log", "export/GRCh38/ld/bt__TRAIT1/chr1.compute_ref_ld.log"
)
native_log_manifest <- suppressWarnings(run_task(c(
  "write-manifest", "--config", config, "--build", "GRCh38",
  "--gene-list", gene_list, "--provenance", provenance, "--group-status", status,
  "--target-summary", manifest_target, "--trait-summary", group_summary,
  "--validation", manifest_validation, "--tool", tool, "--artifact", manifest_htp, manifest_ld,
  "--out", file.path(tmp, "native_log_manifest.tsv")
)))
stopifnot(!is.null(attr(native_log_manifest, "status")), attr(native_log_manifest, "status") != 0)
unlink(native_log)

skipped_status <- write_lines(c(
  "group\ttrait\ttrait_type\tusable_n\tcases\tcontrols\tmodel_sample_count\tmodel_cases\tmodel_controls\tkeep_count\tmodel_keep_sha256\tskipped\tskip_reason\tremeta_eligible",
  "bt__SKIPPED\tSKIPPED\tbt\t10\t5\t5\t0\t\t\t0\ttesthash\tTrue\tbelow_threshold\tFalse"
), "skipped_status.tsv")
skipped_manifest <- file.path(tmp, "skipped_manifest.tsv")
output <- run_task(c(
  "write-manifest", "--config", config, "--build", "GRCh38",
  "--gene-list", gene_list, "--provenance", provenance, "--group-status", skipped_status,
  "--trait-summary", skipped_summary, "--tool", tool, "--out", skipped_manifest
))
stopifnot(is.null(attr(output, "status")), file.exists(skipped_manifest))
skipped_manifest_rows <- read.delim(skipped_manifest, stringsAsFactors = FALSE)
stopifnot(skipped_manifest_rows$value[skipped_manifest_rows$key == "trait:SKIPPED:remeta_exported"] == "false")

# All-skipped runs still reject stale files in the canonical export tree.
canonical_stale <- write_lines(
  "stale", "canonical/results/remeta/export/GRCh38/htp/OLD.PAN.regenie.gz"
)
canonical_manifest <- suppressWarnings(run_task(c(
  "write-manifest", "--config", config, "--build", "GRCh38",
  "--gene-list", gene_list, "--provenance", provenance, "--group-status", skipped_status,
  "--trait-summary", skipped_summary, "--tool", tool,
  "--out", file.path(tmp, "canonical/results/remeta/export/cohort_test.GRCh38.remeta_manifest.tsv")
)))
stopifnot(!is.null(attr(canonical_manifest, "status")), attr(canonical_manifest, "status") != 0)
unlink(canonical_stale)

# A group that becomes skipped may not retain target or validation work files.
mixed_status <- write_lines(c(
  "group\ttrait\ttrait_type\tusable_n\tcases\tcontrols\tmodel_sample_count\tmodel_cases\tmodel_controls\tkeep_count\tmodel_keep_sha256\tskipped\tskip_reason\tremeta_eligible",
  paste(c("bt__TRAIT1", "TRAIT1", "bt", 2, 1, 1, 2, 1, 1, 2,
    keep_sha, "False", "", "True"), collapse = "\t"),
  "bt__SKIPPED\tSKIPPED\tbt\t10\t5\t5\t0\t\t\t0\ttesthash\tTrue\tbelow_threshold\tFalse"
), "mixed_status.tsv")
stale_work <- write_lines("stale", "work/GRCh38/groups/bt__SKIPPED/stale.txt")
stale_work_manifest <- suppressWarnings(run_task(c(
  "write-manifest", "--config", config, "--build", "GRCh38",
  "--gene-list", gene_list, "--provenance", provenance, "--group-status", mixed_status,
  "--target-summary", manifest_target,
  "--trait-summary", group_summary, "--trait-summary", skipped_summary,
  "--validation", manifest_validation, "--tool", tool, "--artifact", manifest_htp, manifest_ld,
  "--out", file.path(tmp, "stale_work_manifest.tsv")
)))
stopifnot(!is.null(attr(stale_work_manifest, "status")), attr(stale_work_manifest, "status") != 0)
unlink(stale_work)

cat("ReMeta cohort helper tests passed\n")
