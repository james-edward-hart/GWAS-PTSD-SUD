#!/usr/bin/env Rscript

# Focused checks for ancestry assignment coverage thresholds.

tmp <- tempfile("strata_test")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

sample_manifest <- file.path(tmp, "samples.tsv")
traits <- file.path(tmp, "traits.tsv")
config <- file.path(tmp, "config.yaml")
ancestry <- file.path(tmp, "ancestry.tsv")

samples <- data.frame(
  FID = paste0("F", seq_len(100)),
  IID = paste0("I", seq_len(100)),
  age = 50,
  age2 = 0,
  sex = "1",
  trait = rep(c("1", "0"), 50),
  stringsAsFactors = FALSE
)
write.table(samples, sample_manifest, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(
  trait_id = "t",
  phenotype_column = "trait",
  case_value = "1",
  control_value = "0",
  missing_values = "NA",
  stringsAsFactors = FALSE
), traits, sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(c(
  "analysis:",
  "  ancestries: [A]",
  "inputs:",
  paste0("  sample_manifest: ", sample_manifest),
  paste0("  trait_registry: ", traits),
  "popmad:",
  "  max_unassigned_fraction: 0.02",
  "warnings:",
  "  min_n: 1",
  "  min_cases: 1",
  "  min_controls: 1"
), config)

write_case <- function(n_assigned, path) {
  write.table(data.frame(
    FID = samples$FID[seq_len(n_assigned)],
    IID = samples$IID[seq_len(n_assigned)],
    ancestry = "A",
    stringsAsFactors = FALSE
  ), path, sep = "\t", quote = FALSE, row.names = FALSE)
}

write_case(98, ancestry)
status <- system2("Rscript", c(
  "scripts/make_strata_files.R",
  "--config", config,
  "--ancestry-file", ancestry,
  "--outdir", file.path(tmp, "pass")
))
stopifnot(status == 0L)
excluded_pass <- read.table(file.path(tmp, "pass", "excluded_ancestries.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
stopifnot(nrow(excluded_pass) == 0L)

write_case(97, ancestry)
status <- system2("Rscript", c(
  "scripts/make_strata_files.R",
  "--config", config,
  "--ancestry-file", ancestry,
  "--outdir", file.path(tmp, "fail")
), stdout = file.path(tmp, "fail.log"), stderr = file.path(tmp, "fail.log"))
stopifnot(!identical(status, 0L))
stopifnot(any(grepl("above allowed", readLines(file.path(tmp, "fail.log")))))

writeLines(c(
  "analysis:",
  "  min_stratum_n: 50",
  "  ancestries: [A, B]",
  "inputs:",
  paste0("  sample_manifest: ", sample_manifest),
  paste0("  trait_registry: ", traits),
  "popmad:",
  "  max_unassigned_fraction: 0.60",
  "warnings:",
  "  min_n: 1",
  "  min_cases: 1",
  "  min_controls: 1"
), config)
write.table(data.frame(
  FID = samples$FID[1:90],
  IID = samples$IID[1:90],
  ancestry = c(rep("A", 80), rep("B", 10)),
  stringsAsFactors = FALSE
), ancestry, sep = "\t", quote = FALSE, row.names = FALSE)
status <- system2("Rscript", c(
  "scripts/make_strata_files.R",
  "--config", config,
  "--ancestry-file", ancestry,
  "--outdir", file.path(tmp, "auto_filter")
))
stopifnot(status == 0L)
active <- read.table(file.path(tmp, "auto_filter", "active_ancestries.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
excluded <- read.table(file.path(tmp, "auto_filter", "excluded_ancestries.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
keep_b <- read.table(file.path(tmp, "auto_filter", "B.keep.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
stopifnot(identical(active$ancestry, "A"))
stopifnot(identical(excluded$ancestry, "B"))
stopifnot(identical(excluded$reason, "below_min_stratum_n"))
stopifnot(nrow(keep_b) == 0L)

cat("Strata assignment tests passed\n")
