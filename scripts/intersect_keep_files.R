#!/usr/bin/env Rscript

# Intersect PLINK keep files by FID/IID.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Accept either repeated --inputs or the older --left/--right flags.
args <- parse_args(defaults = list(inputs = character()), repeated = "inputs")
if (!is.null(args$left)) args$inputs <- c(args$inputs, args$left)
if (!is.null(args$right)) args$inputs <- c(args$inputs, args$right)
require_args(args, c("out"))
if (length(args$inputs) < 2) die("provide at least two keep files to intersect")


# Normalize PLINK keep files that may use FID or #FID.
read_ids <- function(path) {
  rows <- read_tsv(path)
  fid_col <- if ("FID" %in% names(rows)) "FID" else "#FID"
  require_columns(rows, c(fid_col, "IID"), path)
  paste(rows[[fid_col]], rows$IID, sep = "\t")
}


# Keep only samples present in every input file.
common <- Reduce(intersect, lapply(args$inputs, read_ids))
common <- sort(common)
parts <- do.call(rbind, strsplit(common, "\t", fixed = TRUE))

# Rebuild the shared sample IDs as a standard keep table.
out <- if (length(common)) {
  data.frame(FID = parts[, 1], IID = parts[, 2], stringsAsFactors = FALSE)
} else {
  data.frame(FID = character(), IID = character())
}


# Save the intersected keep file.
write_tsv(out, args$out)
cat("Wrote", nrow(out), "intersected samples\n")
