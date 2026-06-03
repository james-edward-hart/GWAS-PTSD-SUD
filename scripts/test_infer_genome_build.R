#!/usr/bin/env Rscript

# Focused regression checks for genome-build inference.


# Use an isolated temporary workspace for generated fixtures.
tmp <- tempfile("genome_build_test")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)


# Write the minimal config needed by infer_genome_build.R.
write_case_config <- function(prefix, marker_file, config_file) {
  writeLines(c(
    "genotypes:",
    "  type: bed",
    paste0("  prefix: ", prefix),
    "genome_build:",
    paste0("  marker_file: ", marker_file),
    "  min_markers: 1",
    "  min_match_fraction: 0.5",
    "  min_marker_margin: 1",
    "  min_fraction_margin: 0.0"
  ), config_file)
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
status <- system2("Rscript", c("scripts/infer_genome_build.R", "--config", file.path(case1, "config.yaml"),
  "--out", file.path(case1, "build.txt"), "--details", file.path(case1, "details.tsv")))
stopifnot(status == 0L)
stopifnot(readLines(file.path(case1, "build.txt")) == "GRCh37")


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
status <- system2("Rscript", c("scripts/infer_genome_build.R", "--config", file.path(case2, "config.yaml"),
  "--out", file.path(case2, "build.txt"), "--details", file.path(case2, "details.tsv")),
  stdout = file.path(case2, "run.log"), stderr = file.path(case2, "run.log"))
stopifnot(!identical(status, 0L))
stopifnot(any(grepl("no genome-build marker positions matched", readLines(file.path(case2, "run.log")))))


# Report test completion.
cat("Genome-build inference tests passed\n")
