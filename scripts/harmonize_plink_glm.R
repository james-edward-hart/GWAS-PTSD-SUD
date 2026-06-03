#!/usr/bin/env Rscript

# Convert PLINK2 .glm output into the pipeline's stable summary-stat schema.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Parse the PLINK prefix and output metadata.
args <- parse_args(defaults = list(test = "ADD"))
require_args(args, c("plink-prefix", "trait", "ancestry", "build", "out"))


# Locate the first PLINK2 .glm result for this prefix.
matches <- sort(Sys.glob(paste0(args[["plink-prefix"]], "*.glm.*")))
matches <- matches[!grepl("\\.log$", matches)]
if (!length(matches)) die("no PLINK2 .glm output found for prefix: ", args[["plink-prefix"]])
glm_file <- matches[[1]]


# Read the result and normalize missing values from optional columns.
glm <- read_tsv(glm_file)
get_first <- function(df, names) {
  for (name in names) {
    if (name %in% colnames(df)) {
      value <- as.character(df[[name]])
      value[!nzchar(value) | value %in% c("NA", ".")] <- "NA"
      return(value)
    }
  }
  rep("NA", nrow(df))
}


# Keep the requested PLINK2 test term.
test <- get_first(glm, "TEST")
available_tests <- sort(unique(test[nzchar(test) & test != "NA"]))
keep <- if (nzchar(args$test)) test == args$test else rep(TRUE, nrow(glm))
glm <- glm[keep, , drop = FALSE]
test <- test[keep]
if (!nrow(glm)) {
  observed <- if (length(available_tests)) paste(available_tests, collapse = ", ") else "none"
  die("PLINK2 .glm output for ", args[["plink-prefix"]], " contains no rows for requested test '",
    args$test, "'; available TEST values: ", observed)
}


# Convert odds ratios to log-odds when PLINK did not report BETA.
or_value <- get_first(glm, "OR")
beta <- get_first(glm, "BETA")
needs_beta <- beta == "NA"
or_num <- suppressWarnings(as.numeric(or_value))
beta[needs_beta & is.finite(or_num) & or_num > 0] <- as.character(log(or_num[needs_beta & is.finite(or_num) & or_num > 0]))


# Emit a stable summary-statistics schema for downstream reports.
out <- data.frame(
  trait = args$trait,
  ancestry = args$ancestry,
  build = args$build,
  chrom = get_first(glm, c("#CHROM", "CHROM")),
  pos = get_first(glm, "POS"),
  variant_id = get_first(glm, "ID"),
  ref = get_first(glm, "REF"),
  alt = get_first(glm, "ALT"),
  effect_allele = get_first(glm, "A1"),
  a1_freq = get_first(glm, c("A1_FREQ", "ALT_FREQS", "A1_FREQS")),
  mac = get_first(glm, "MAC"),
  info = get_first(glm, c("MACH_R2", "INFO", "R2")),
  test = test,
  n = get_first(glm, "OBS_CT"),
  beta_or_log_or = beta,
  se = get_first(glm, c("SE", "LOG(OR)_SE")),
  or = or_value,
  p = get_first(glm, "P"),
  stringsAsFactors = FALSE
)


# Write the harmonized table.
write_tsv(out, args$out)
cat("Harmonized", glm_file, "->", args$out, "\n")
