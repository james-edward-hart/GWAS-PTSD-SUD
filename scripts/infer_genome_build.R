#!/usr/bin/env Rscript

# Infer the study dataset genome build from build-specific marker positions.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Parse effective config and output paths.
args <- parse_args()
require_args(args, c("config", "out", "details"))


# Load the offline build-marker table.
config <- load_config(args$config)
markers <- read_tsv(config$genome_build$marker_file)
require_columns(markers, c("variant_id", "chrom", "build", "pos"), "genome-build marker file")
if (!nrow(markers)) die("genome-build marker file is empty: ", config$genome_build$marker_file)
markers$chrom <- clean_chrom(markers$chrom)
markers$pos <- as.integer(markers$pos)


# Read study variants from the configured PLINK dataset.
kind <- tolower(config$genotypes$type)
prefix <- config$genotypes$prefix
if (kind == "bed") {
  bim <- read.table(paste0(prefix, ".bim"), stringsAsFactors = FALSE, quote = "", comment.char = "")
  variants <- data.frame(ID = bim[[2]], chrom = clean_chrom(bim[[1]]), pos = as.integer(bim[[4]]))
} else if (kind == "pgen") {
  pvar <- read_tsv(paste0(prefix, ".pvar"))
  chrom_col <- if ("#CHROM" %in% names(pvar)) "#CHROM" else "CHROM"
  variants <- data.frame(ID = pvar$ID, chrom = clean_chrom(pvar[[chrom_col]]), pos = as.integer(pvar$POS))
} else {
  die("unsupported genotype type '", kind, "'. Use 'pgen' or 'bed'.")
}


# Index markers by rsID and by coordinate for fallback matching.
builds <- sort(unique(markers$build))
scores <- setNames(rep(0L, length(builds)), builds)
checked <- 0L
by_id <- split(seq_len(nrow(markers)), markers$variant_id)
by_pos <- split(seq_len(nrow(markers)), paste(markers$chrom, markers$pos, sep = "\t"))


# rsID matches take precedence over coordinate fallback.
for (i in seq_len(nrow(variants))) {
  rows <- by_id[[variants$ID[[i]]]]
  if (is.null(rows)) rows <- by_pos[[paste(variants$chrom[[i]], variants$pos[[i]], sep = "\t")]]
  if (is.null(rows)) next
  checked <- checked + 1L
  matching <- rows[markers$chrom[rows] == variants$chrom[[i]] & markers$pos[rows] == variants$pos[[i]]]
  for (build in markers$build[matching]) scores[[build]] <- scores[[build]] + 1L
}


# Choose the build with the most matched markers.
if (!checked || max(scores) == 0) die("no genome-build marker positions matched the genotype data")

ordered <- order(-scores, names(scores))
build <- names(scores)[ordered[[1]]]
matches <- scores[[build]]
second <- if (length(scores) > 1) scores[[ordered[[2]]]] else 0L
fraction <- matches / checked
marker_margin <- matches - second
fraction_margin <- marker_margin / checked


# Enforce configured minimum support and winner margin.
minimum <- as.integer(config$genome_build$min_markers %||% 3)
min_fraction <- as.numeric(config$genome_build$min_match_fraction %||% 0.8)
min_marker_margin <- as.integer(config$genome_build$min_marker_margin %||% 1)
min_fraction_margin <- as.numeric(config$genome_build$min_fraction_margin %||% 0)

if (matches < minimum) die("only ", matches, " genome-build markers supported ", build, "; minimum required is ", minimum)
if (fraction < min_fraction) die("genome-build inference is weak for ", build, ": ", matches, "/", checked, " markers matched")
if (marker_margin < min_marker_margin) die("genome-build inference margin is weak for ", build)
if (fraction_margin < min_fraction_margin) die("genome-build inference margin fraction is weak for ", build)


# Build a per-build details table for reports.
fmt <- function(x) ifelse(abs(x - round(x)) < 1e-12, sprintf("%.1f", x), as.character(x))
details <- data.frame(
  build = names(scores),
  matching_markers = as.integer(scores),
  checked_markers = checked,
  match_fraction = fmt(as.numeric(scores) / checked),
  winner_marker_margin = ifelse(names(scores) == build, as.character(marker_margin), ""),
  winner_fraction_margin = ifelse(names(scores) == build, fmt(fraction_margin), ""),
  selected = ifelse(names(scores) == build, "True", "False"),
  stringsAsFactors = FALSE
)


# Save the selected build and scoring details.
ensure_parent(args$out)
writeLines(build, args$out)
write_tsv(details, args$details)
cat("Inferred genome build:", build, paste0("(", matches, "/", checked, " marker matches)"), "\n")
