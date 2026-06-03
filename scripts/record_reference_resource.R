#!/usr/bin/env Rscript

# Record an already-downloaded reference resource in the local TSV manifest.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Keep manifest columns stable for append/update operations.
fields <- c(
  "resource", "file_role", "version_or_build", "genome_build", "source_url",
  "local_path", "download_date", "checksum_algorithm", "sha256",
  "metadata_columns", "preprocessing_notes"
)


# Create an empty manifest with the expected schema.
empty_manifest <- function() {
  setNames(data.frame(matrix(nrow = 0, ncol = length(fields)), stringsAsFactors = FALSE), fields)
}


# Read any existing manifest and backfill missing columns.
read_manifest <- function(path) {
  if (!file.exists(path)) return(empty_manifest())
  rows <- read_tsv(path)
  missing <- setdiff(fields, names(rows))
  for (name in missing) rows[[name]] <- ""
  rows[fields]
}


# Parse the manifest row fields and optional existence check.
args <- parse_args(
  defaults = list(manifest = "resources/manifests/reference_data.tsv", sha256 = "", "metadata-columns" = "n/a", notes = "", "require-existing" = FALSE),
  flags = "require-existing"
)
require_args(args, c("resource", "file-role", "version-or-build", "genome-build", "source-url", "local-path", "download-date"))

if (isTRUE(args[["require-existing"]]) && !file.exists(args[["local-path"]])) {
  die("local path does not exist: ", args[["local-path"]])
}


# Compute a SHA-256 only when the local path is a readable file.
checksum <- args$sha256
checksum_algorithm <- if (nzchar(checksum)) "sha256" else ""
if (!nzchar(checksum) && file.exists(args[["local-path"]]) && !dir.exists(args[["local-path"]])) {
  checksum <- sha256_file(args[["local-path"]])
  checksum_algorithm <- "sha256"
}


# Build one replacement manifest row for this resource.
row <- data.frame(
  resource = args$resource,
  file_role = args[["file-role"]],
  version_or_build = args[["version-or-build"]],
  genome_build = args[["genome-build"]],
  source_url = args[["source-url"]],
  local_path = args[["local-path"]],
  download_date = args[["download-date"]],
  checksum_algorithm = checksum_algorithm,
  sha256 = checksum,
  metadata_columns = args[["metadata-columns"]],
  preprocessing_notes = args$notes,
  stringsAsFactors = FALSE
)


# Replace any prior row with the same resource key.
manifest <- read_manifest(args$manifest)
manifest <- manifest[manifest$resource != args$resource, , drop = FALSE]
write_tsv(rbind(manifest, row), args$manifest)
