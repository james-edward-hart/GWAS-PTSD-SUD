#!/usr/bin/env Rscript

# Focused regression checks for reference-package fingerprint and build resolution.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))

fixture <- "data/example/reference_package"
tmp <- tempfile("reference_package_test")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
pkg <- file.path(tmp, "pkg")
invisible(file.copy(fixture, tmp, recursive = TRUE))
invisible(file.rename(file.path(tmp, basename(fixture)), pkg))

observed <- reference_package_fingerprint(pkg, verify_hashes = TRUE)
stopifnot(identical(observed, read_first_line(file.path(pkg, "content_fingerprint.sha256"))))

config <- list(
  project = list(),
  reference_package = list(root = pkg, fingerprint = observed),
  ancestry_reference = list(enabled = TRUE),
  admixture = list(enabled = TRUE, labels = c("AFR", "AMR", "EAS", "EUR", "SAS"))
)
resolved <- resolve_reference_package_config(config, "GRCh37", verify_hashes = TRUE)
stopifnot(identical(resolved$ancestry_reference$reference_genome_build, "GRCh37"))
stopifnot(grepl("popmad.GRCh37", resolved$ancestry_reference$reference_genotypes$prefix, fixed = TRUE))
stopifnot(identical(resolved$ancestry_reference$variant_set, "pre_ld_pruned"))
stopifnot(identical(resolved$ancestry_reference$source_sample_set, "unrelateds_without_outliers"))
stopifnot(grepl("admixture.GRCh37", resolved$admixture$reference_genotypes$prefix, fixed = TRUE))

bad_fingerprint <- config
bad_fingerprint$reference_package$fingerprint <- paste0(substr(observed, 1, 63), ifelse(substr(observed, 64, 64) == "0", "1", "0"))
status <- tryCatch({
  resolve_reference_package_config(bad_fingerprint, "GRCh37", verify_hashes = TRUE)
  TRUE
}, error = function(e) FALSE)
stopifnot(!status)

status <- tryCatch({
  resolve_reference_package_config(config, "GRCh36", verify_hashes = TRUE)
  TRUE
}, error = function(e) FALSE)
stopifnot(!status)

pkg_sidecar <- file.path(tmp, "pkg_sidecar")
invisible(file.copy(fixture, tmp, recursive = TRUE))
invisible(file.rename(file.path(tmp, basename(fixture)), pkg_sidecar))
writeLines("macOS metadata", file.path(pkg_sidecar, "._file_manifest.tsv"))
writeLines("macOS metadata", file.path(pkg_sidecar, ".DS_Store"))
writeLines("macOS metadata", file.path(pkg_sidecar, "panels", "._admixture.GRCh37"))
dir.create(file.path(pkg_sidecar, "__MACOSX", "panels"), recursive = TRUE)
writeLines("macOS metadata", file.path(pkg_sidecar, "__MACOSX", "panels", "._popmad.GRCh37"))
sidecar_fingerprint <- reference_package_fingerprint(pkg_sidecar, verify_hashes = TRUE)
stopifnot(identical(sidecar_fingerprint, observed))

unlink(file.path(pkg, "panel_manifest.tsv"))
status <- tryCatch({
  reference_package_fingerprint(pkg, verify_hashes = TRUE)
  TRUE
}, error = function(e) FALSE)
stopifnot(!status)

pkg_raw <- file.path(tmp, "pkg_raw")
invisible(file.copy(fixture, tmp, recursive = TRUE))
invisible(file.rename(file.path(tmp, basename(fixture)), pkg_raw))
dir.create(file.path(pkg_raw, "raw_intermediate.mt"))
writeLines("raw hail marker", file.path(pkg_raw, "raw_intermediate.mt", "metadata.json"))
status <- tryCatch({
  reference_package_fingerprint(pkg_raw, verify_hashes = TRUE)
  TRUE
}, error = function(e) FALSE)
stopifnot(!status)

pkg_extra <- file.path(tmp, "pkg_extra")
invisible(file.copy(fixture, tmp, recursive = TRUE))
invisible(file.rename(file.path(tmp, basename(fixture)), pkg_extra))
writeLines("unmanifested", file.path(pkg_extra, "panels", "popmad.GRCh37", ".metadata.tsv.crc"))
status <- tryCatch({
  reference_package_fingerprint(pkg_extra, verify_hashes = TRUE)
  TRUE
}, error = function(e) FALSE)
stopifnot(!status)

pkg_escape <- file.path(tmp, "pkg_escape")
invisible(file.copy(fixture, tmp, recursive = TRUE))
invisible(file.rename(file.path(tmp, basename(fixture)), pkg_escape))
panels <- read_tsv(file.path(pkg_escape, "panel_manifest.tsv"))
files <- read_tsv(file.path(pkg_escape, "file_manifest.tsv"))
panels$metadata_path[[1]] <- "/tmp/metadata.tsv"
status <- tryCatch({
  validate_reference_panel_paths(panels, files)
  TRUE
}, error = function(e) FALSE)
stopifnot(!status)

pkg_bad_ids <- file.path(tmp, "pkg_bad_ids")
invisible(file.copy(fixture, tmp, recursive = TRUE))
invisible(file.rename(file.path(tmp, basename(fixture)), pkg_bad_ids))
pvar_path <- file.path(pkg_bad_ids, "panels", "popmad.GRCh37", "popmad.GRCh37.pvar")
pvar <- read_tsv(pvar_path)
pvar$ID[[1]] <- "[\"rs1\"]"
write_tsv(pvar, pvar_path)
status <- tryCatch({
  validate_reference_panel_variant_ids(pkg_bad_ids, read_tsv(file.path(pkg_bad_ids, "panel_manifest.tsv")))
  TRUE
}, error = function(e) FALSE)
stopifnot(!status)

cat("Reference package tests passed\n")
