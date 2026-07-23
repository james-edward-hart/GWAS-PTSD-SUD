#!/usr/bin/env Rscript

# Focused regression checks for genome-build inference.


# Use an isolated temporary workspace for generated fixtures.
tmp <- tempfile("genome_build_test")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)


# Write the minimal config needed by infer_genome_build.R.
write_case_config <- function(prefix, marker_file, config_file, kind = "bed") {
  writeLines(c(
    "genotypes:",
    paste0("  type: ", kind),
    paste0("  prefix: ", prefix),
    "genome_build:",
    paste0("  marker_file: ", marker_file),
    "  min_markers: 1",
    "  min_match_fraction: 0.5",
    "  min_marker_margin: 1",
    "  min_fraction_margin: 0.0"
  ), config_file)
}


# Run a fixture while safely preserving paths containing spaces.
run_inference <- function(case_dir, log = "") {
  command_args <- c(
    "scripts/infer_genome_build.R",
    "--config", file.path(case_dir, "config.yaml"),
    "--out", file.path(case_dir, "build.txt"),
    "--details", file.path(case_dir, "details.tsv")
  )
  command_args <- vapply(command_args, shQuote, character(1), type = "sh")
  if (nzchar(log)) {
    return(system2("Rscript", command_args, stdout = log, stderr = log))
  }
  system2("Rscript", command_args)
}


# Case 1: coordinate fallback should infer GRCh37.
case1 <- file.path(tmp, "case1")
dir.create(case1)
writeLines(c(
  "variant_id\tchrom\tbuild\tpos",
  "rs1\t1\tGRCh37\t100",
  "rs1\t1\tGRCh38\t110",
  "rs2\t2\tGRCh37\t200",
  "rs2\t2\tGRCh38\t210"
), file.path(case1, "markers.tsv"))
writeLines(c(
  "1\tchr1:100:A:G\t0\t100\tA\tG",
  "2\tchr2:200:C:T\t0\t200\tC\tT"
), file.path(case1, "study.bim"))
write_case_config(file.path(case1, "study"), file.path(case1, "markers.tsv"), file.path(case1, "config.yaml"))


# Run the positive case and check the selected build.
status <- run_inference(case1)
stopifnot(status == 0L)
stopifnot(readLines(file.path(case1, "build.txt")) == "GRCh37")
details <- read.delim(file.path(case1, "details.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
selected <- details[details$selected == "True", , drop = FALSE]
stopifnot(selected$matching_markers == 2L)
stopifnot(selected$checked_markers == 2L)


# Case 2: an rsID with mismatched coordinates should fail.
case2 <- file.path(tmp, "case2")
dir.create(case2)
writeLines(c(
  "variant_id\tchrom\tbuild\tpos",
  "rs1\t1\tGRCh37\t100",
  "rs1\t1\tGRCh38\t105",
  "rs2\t1\tGRCh37\t110",
  "rs2\t1\tGRCh38\t120"
), file.path(case2, "markers.tsv"))
writeLines("1\trs1\t0\t110\tA\tG", file.path(case2, "study.bim"))
write_case_config(file.path(case2, "study"), file.path(case2, "markers.tsv"), file.path(case2, "config.yaml"))


# Run the negative case and confirm the expected error.
case2_log <- file.path(case2, "run.log")
status <- run_inference(case2, case2_log)
stopifnot(!identical(status, 0L))
stopifnot(any(grepl("no genome-build marker positions matched", readLines(case2_log))))


# Case 3: scan reordered PVAR columns, mixed chr prefixes, and paths with spaces.
case3 <- file.path(tmp, "case 3 with spaces")
dir.create(case3)
writeLines(c(
  "variant_id\tchrom\tbuild\tpos",
  "rs1\t1\tGRCh37\t100",
  "rs1\t1\tGRCh38\t110"
), file.path(case3, "markers.tsv"))
noise_count <- 200005L
noise <- paste(
  paste0("noise", seq_len(noise_count)),
  "G",
  22,
  "PASS",
  seq_len(noise_count),
  "A",
  sep = "\t"
)
writeLines(c(
  "##fileformat=VCFv4.2",
  "ID\tALT\t#CHROM\tFILTER\tPOS\tREF",
  noise,
  "rs1\tG\tChr1\tPASS\t100\tA"
), file.path(case3, "study.pvar"))
write_case_config(
  file.path(case3, "study"),
  file.path(case3, "markers.tsv"),
  file.path(case3, "config.yaml"),
  kind = "pgen"
)

status <- run_inference(case3)
stopifnot(status == 0L)
stopifnot(readLines(file.path(case3, "build.txt")) == "GRCh37")
details <- read.delim(file.path(case3, "details.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
selected <- details[details$selected == "True", , drop = FALSE]
stopifnot(selected$matching_markers == 1L)
stopifnot(selected$checked_markers == 1L)


# Case 4: duplicate metadata candidates retain the original scoring denominator.
case4 <- file.path(tmp, "case4")
dir.create(case4)
writeLines(c(
  "variant_id\tchrom\tbuild\tpos",
  "rs1\t1\tGRCh37\t100",
  "rs1\t1\tGRCh38\t110"
), file.path(case4, "markers.tsv"))
writeLines(c(
  "chr1\trs1\t0\t100\tA\tG",
  "1\trs1\t0\t100\tA\tG"
), file.path(case4, "study.bim"))
write_case_config(file.path(case4, "study"), file.path(case4, "markers.tsv"), file.path(case4, "config.yaml"))

status <- run_inference(case4)
stopifnot(status == 0L)
details <- read.delim(file.path(case4, "details.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
selected <- details[details$selected == "True", , drop = FALSE]
stopifnot(selected$matching_markers == 2L)
stopifnot(selected$checked_markers == 2L)


# Case 5: malformed PVAR headers fail in the scanner and propagate to R.
case5 <- file.path(tmp, "case5")
dir.create(case5)
writeLines(c(
  "variant_id\tchrom\tbuild\tpos",
  "rs1\t1\tGRCh37\t100",
  "rs1\t1\tGRCh38\t110"
), file.path(case5, "markers.tsv"))
writeLines(c(
  "##fileformat=VCFv4.2",
  "#CHROM\tPOS\tREF\tALT",
  "1\t100\tA\tG"
), file.path(case5, "study.pvar"))
write_case_config(
  file.path(case5, "study"),
  file.path(case5, "markers.tsv"),
  file.path(case5, "config.yaml"),
  kind = "pgen"
)

case5_log <- file.path(case5, "run.log")
status <- run_inference(case5, case5_log)
stopifnot(!identical(status, 0L))
errors <- readLines(case5_log)
stopifnot(any(grepl("PVAR primary header is missing required CHROM, POS, or ID columns", errors)))
stopifnot(any(grepl("genome-build metadata scanner failed", errors)))


# Case 6: malformed BIM rows fail with their source row number.
case6 <- file.path(tmp, "case6")
dir.create(case6)
writeLines(c(
  "variant_id\tchrom\tbuild\tpos",
  "rs1\t1\tGRCh37\t100",
  "rs1\t1\tGRCh38\t110"
), file.path(case6, "markers.tsv"))
writeLines("1\trs1\t0\t100\tA", file.path(case6, "study.bim"))
write_case_config(file.path(case6, "study"), file.path(case6, "markers.tsv"), file.path(case6, "config.yaml"))

case6_log <- file.path(case6, "run.log")
status <- run_inference(case6, case6_log)
stopifnot(!identical(status, 0L))
stopifnot(any(grepl("BIM row 1 must contain at least 6 columns", readLines(case6_log))))


# Report test completion.
cat("Genome-build inference tests passed\n")
