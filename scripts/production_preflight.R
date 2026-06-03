#!/usr/bin/env Rscript

# Print production-readiness warnings before a SLURM pilot or full run.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


args <- parse_args(defaults = list(config = "config/config.yaml", profile = "profiles/slurm/config.yaml"),
  flags = c("allow-existing-results", "allow-test-run"))

config <- load_config(args$config)
run_mode <- config$project$run_mode %||% "test"
results_dir <- "results"
warnings <- character()
failures <- character()

add_warning <- function(...) warnings <<- c(warnings, paste0(...))
add_failure <- function(...) failures <<- c(failures, paste0(...))


# Production config gates duplicated here for operator-facing feedback.
if (run_mode != "production") {
  if (truthy(args[["allow-test-run"]] %||% FALSE)) {
    add_warning("project.run_mode is '", run_mode, "', not 'production'; continuing because --allow-test-run was set")
  } else {
    add_failure("project.run_mode must be 'production' for production preflight; use --allow-test-run only for rehearsal configs")
  }
}
if ((config$inputs$ancestry_mode %||% "") != "computed") add_failure("production requires inputs.ancestry_mode: computed")
if (!truthy(config$ancestry_reference$enabled %||% FALSE)) add_failure("production requires ancestry_reference.enabled: true")
if (!truthy(config$admixture$enabled %||% FALSE)) add_failure("production requires admixture.enabled: true")
if (!truthy(config$sex_check$enabled %||% TRUE)) add_failure("production requires sex_check.enabled: true")
if ((config$sex_check$action %||% "warn") != "exclude") add_failure("production requires sex_check.action: exclude")
if (truthy(config$sex_check$allow_no_sex_markers %||% FALSE)) {
  add_warning("sex_check.allow_no_sex_markers is true; confirm autosome-only study data have an external sex-QC record")
}
if (!blank(config$resources$input_manifest %||% "") && !file.exists(config$resources$input_manifest)) {
  add_failure("resources.input_manifest does not exist: ", config$resources$input_manifest)
}

root <- config$reference_package$root %||% ""
if (blank(root)) {
  add_failure("reference_package.root is empty")
} else if (!dir.exists(root)) {
  add_failure("reference_package.root does not exist: ", root)
} else {
  observed <- tryCatch(reference_package_fingerprint(root, verify_hashes = TRUE), error = function(e) e)
  if (inherits(observed, "error")) {
    add_failure(conditionMessage(observed))
  } else if (blank(config$reference_package$fingerprint %||% "")) {
    add_failure("reference_package.fingerprint is empty")
  } else if (!identical(tolower(config$reference_package$fingerprint), observed)) {
    add_failure("reference_package.fingerprint does not match observed content fingerprint")
  }
}


# SLURM profile checks are warnings because site policy varies.
if (!file.exists(args$profile)) {
  add_warning("SLURM profile not found: ", args$profile)
} else {
  profile_text <- readLines(args$profile, warn = FALSE)
  if (!any(grepl("^\\s*slurm_account\\s*:", profile_text))) {
    add_warning("SLURM profile has no active slurm_account")
  }
  if (!any(grepl("^\\s*conda-prefix\\s*:", profile_text))) {
    add_warning("SLURM profile has no active conda-prefix; --use-conda will use Snakemake defaults")
  }
  if (any(grepl("my_account|/path/to/shared", profile_text))) {
    add_warning("SLURM profile still contains placeholder account or conda-prefix text")
  }
}


# Disk and path hygiene checks.
if (dir.exists(results_dir) && length(list.files(results_dir, all.files = TRUE, no.. = TRUE)) &&
    !truthy(args[["allow-existing-results"]] %||% FALSE)) {
  add_failure("results/ is not clean; use a fresh results directory or rerun preflight with --allow-existing-results")
}

df <- tryCatch(system2("df", c("-Pk", "."), stdout = TRUE, stderr = TRUE), error = function(e) character())
if (length(df) >= 2) {
  fields <- strsplit(df[[2]], "\\s+")[[1]]
  available_kb <- suppressWarnings(as.numeric(fields[[4]]))
  if (!is.na(available_kb) && available_kb < 500 * 1024 * 1024) {
    add_warning("available workspace disk is below 500 GB by df -Pk; confirm production storage sizing")
  }
}

paths_to_review <- c(config$genotypes$prefix, root, results_dir, config$resources$input_manifest %||% "")
private_hits <- paths_to_review[grepl("protected|controlled|private|phi|pii", paths_to_review, ignore.case = TRUE)]
if (length(private_hits)) {
  add_warning("path names suggest protected data; confirm logs/manifests do not expose private identifiers: ",
    paste(private_hits, collapse = ", "))
}


# Print a compact operator summary and fail only for hard gates.
cat("Production preflight\n")
cat("- Config: ", args$config, "\n", sep = "")
cat("- Profile: ", args$profile, "\n", sep = "")
cat("- Run mode: ", run_mode, "\n", sep = "")
cat("- Reference package: ", ifelse(blank(root), "NA", root), "\n", sep = "")
if (length(warnings)) {
  cat("\nWarnings:\n")
  for (item in warnings) cat("- ", item, "\n", sep = "")
}
if (length(failures)) {
  cat("\nFailures:\n")
  for (item in failures) cat("- ", item, "\n", sep = "")
  quit(status = 1)
}
cat("\nPreflight passed\n")
