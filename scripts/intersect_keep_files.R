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


# Normalize PLINK keep files that may use headers, no headers, or IID aliases.
read_ids <- function(path) {
  read_id_file(path, path)
}


ids <- lapply(args$inputs, read_ids)
base <- ids[[1]]
keep <- rep(TRUE, nrow(base))
if (length(ids) > 1) {
  for (i in seq.int(2, length(ids))) {
    keep <- keep & !is.na(match_sample_rows(base, sample_key_map(ids[[i]], args$inputs[[i]])))
  }
}


# Save the intersected keep file using the first input's canonical sample IDs.
out <- base[keep, c("FID", "IID"), drop = FALSE]
write_tsv(out, args$out)
cat("Wrote", nrow(out), "intersected samples\n")
