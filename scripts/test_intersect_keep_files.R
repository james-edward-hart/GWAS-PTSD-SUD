#!/usr/bin/env Rscript

# Regression checks for intersecting pipeline keep files with PLINK ID outputs.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))

repo <- normalizePath(file.path(script_dir, ".."))
tmp <- tempfile("intersect-keep-files-")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)


stratum <- file.path(tmp, "EUR.keep.tsv")
unrelated <- file.path(tmp, "unrelated.king.cutoff.in.id")
sex_keep <- file.path(tmp, "sex_checked.keep.tsv")
out <- file.path(tmp, "EUR.unrelated.keep.tsv")

write_tsv(data.frame(
  FID = c("FAM1", "FAM2", "FAM3"),
  IID = c("IID1", "IID2", "IID3"),
  stringsAsFactors = FALSE
), stratum)

writeLines(c(
  "0\tIID1",
  "IID2\tIID2",
  "0\tIID4"
), unrelated)

write_tsv(data.frame(
  FID = c("FAM1", "FAM2", "FAM3"),
  IID = c("IID1", "IID2", "IID3"),
  stringsAsFactors = FALSE
), sex_keep)

status <- system2("Rscript", c(
  file.path(repo, "scripts", "intersect_keep_files.R"),
  "--inputs", stratum, unrelated, sex_keep,
  "--out", out
), stdout = TRUE, stderr = TRUE)
if (!is.null(attr(status, "status"))) {
  stop("intersect_keep_files.R failed:\n", paste(status, collapse = "\n"), call. = FALSE)
}

rows <- read_tsv(out)
stopifnot(identical(rows$FID, c("FAM1", "FAM2")))
stopifnot(identical(rows$IID, c("IID1", "IID2")))

cat("intersect keep-file tests passed\n")
