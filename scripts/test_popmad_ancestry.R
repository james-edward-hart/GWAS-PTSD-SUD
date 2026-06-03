#!/usr/bin/env Rscript

# Focused regression checks for POP-MaD ancestry assignment.


# Use an isolated temporary workspace for generated PC files.
tmp <- tempfile("popmad_test")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)


# Define reference and study fixture paths.
reference <- file.path(tmp, "reference.tsv")
study <- file.path(tmp, "study.tsv")


# Build two compact reference populations in PC space.
ref <- data.frame(FID = character(), IID = character(), population = character(),
  super_population = character(), PC1 = character(), PC2 = character())
for (i in 0:11) {
  ref <- rbind(ref,
    data.frame(FID = paste0("A", i), IID = paste0("A", i), population = "POP_A", super_population = "ANC_A",
      PC1 = as.character(i * 0.01), PC2 = as.character(i * 0.01)),
    data.frame(FID = paste0("B", i), IID = paste0("B", i), population = "POP_B", super_population = "ANC_B",
      PC1 = as.character(10 + i * 0.01), PC2 = as.character(10 + i * 0.01))
  )
}

# Write one assignable study sample and one outlier.
write.table(ref, reference, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(FID = c("S1", "S2"), IID = c("S1", "S2"), PC1 = c("0.03", "100"), PC2 = c("0.03", "100")),
  study, sep = "\t", quote = FALSE, row.names = FALSE)


# Run the POP-MaD assignment script.
status <- system2("Rscript", c(
  "scripts/infer_popmad_ancestry.R",
  "--study-pcs", study,
  "--reference-pcs", reference,
  "--ancestries", "ANC_A,ANC_B",
  "--pcs", "2",
  "--assignments", file.path(tmp, "assignments.tsv"),
  "--study-pcs-out", file.path(tmp, "study_pcs.tsv"),
  "--distances", file.path(tmp, "distances.tsv"),
  "--reference-outliers", file.path(tmp, "outliers.tsv"),
  "--model-summary", file.path(tmp, "models.tsv"),
  "--within-pcs", file.path(tmp, "within.tsv"),
  "--excluded", file.path(tmp, "excluded.tsv"),
  "--counts", file.path(tmp, "counts.tsv"),
  "--min-confidence", "0.05",
  "--outlier-sd", "4",
  "--min-reference-n", "10"
))
stopifnot(status == 0L)


# Confirm the expected assignment and exclusion counts.
count_rows <- function(path) max(length(readLines(path)) - 1, 0)
stopifnot(count_rows(file.path(tmp, "assignments.tsv")) == 1)
stopifnot(count_rows(file.path(tmp, "excluded.tsv")) == 1)
stopifnot(count_rows(file.path(tmp, "models.tsv")) == 2)

status <- system2("Rscript", c(
  "scripts/infer_popmad_ancestry.R",
  "--study-pcs", study,
  "--reference-pcs", reference,
  "--ancestries", "ANC_A,ANC_B",
  "--pcs", "2",
  "--assignments", file.path(tmp, "low_assignments.tsv"),
  "--study-pcs-out", file.path(tmp, "low_study_pcs.tsv"),
  "--distances", file.path(tmp, "low_distances.tsv"),
  "--reference-outliers", file.path(tmp, "low_outliers.tsv"),
  "--model-summary", file.path(tmp, "low_models.tsv"),
  "--within-pcs", file.path(tmp, "low_within.tsv"),
  "--excluded", file.path(tmp, "low_excluded.tsv"),
  "--counts", file.path(tmp, "low_counts.tsv"),
  "--min-confidence", "0.05",
  "--outlier-sd", "4"
), stdout = file.path(tmp, "low.log"), stderr = file.path(tmp, "low.log"))
stopifnot(!identical(status, 0L))
stopifnot(any(grepl("no POP-MaD reference population models retained", readLines(file.path(tmp, "low.log")))))

dup_reference <- file.path(tmp, "duplicate_reference.tsv")
write.table(rbind(ref, ref[1, ]), dup_reference, sep = "\t", quote = FALSE, row.names = FALSE)
status <- system2("Rscript", c(
  "scripts/infer_popmad_ancestry.R",
  "--study-pcs", study,
  "--reference-pcs", dup_reference,
  "--ancestries", "ANC_A,ANC_B",
  "--pcs", "2",
  "--assignments", file.path(tmp, "dup_assignments.tsv"),
  "--study-pcs-out", file.path(tmp, "dup_study_pcs.tsv"),
  "--distances", file.path(tmp, "dup_distances.tsv"),
  "--reference-outliers", file.path(tmp, "dup_outliers.tsv"),
  "--model-summary", file.path(tmp, "dup_models.tsv"),
  "--within-pcs", file.path(tmp, "dup_within.tsv"),
  "--excluded", file.path(tmp, "dup_excluded.tsv"),
  "--counts", file.path(tmp, "dup_counts.tsv"),
  "--min-confidence", "0.05",
  "--outlier-sd", "4",
  "--min-reference-n", "10"
), stdout = file.path(tmp, "dup.log"), stderr = file.path(tmp, "dup.log"))
stopifnot(!identical(status, 0L))
stopifnot(any(grepl("duplicate FID/IID", readLines(file.path(tmp, "dup.log")))))


# A study sample closest to an unconfigured super-population should be excluded
# rather than forced into the nearest configured GWAS ancestry.
unconfigured_reference <- file.path(tmp, "unconfigured_reference.tsv")
ref_unconfigured <- ref
for (i in 0:11) {
  ref_unconfigured <- rbind(ref_unconfigured,
    data.frame(FID = paste0("X", i), IID = paste0("X", i), population = "POP_X", super_population = "ANC_X",
      PC1 = as.character(0.02 + i * 0.001), PC2 = as.character(0.02 + i * 0.001))
  )
}
write.table(ref_unconfigured, unconfigured_reference, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(FID = "SX", IID = "SX", PC1 = "0.025", PC2 = "0.025"),
  file.path(tmp, "unconfigured_study.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
status <- system2("Rscript", c(
  "scripts/infer_popmad_ancestry.R",
  "--study-pcs", file.path(tmp, "unconfigured_study.tsv"),
  "--reference-pcs", unconfigured_reference,
  "--ancestries", "ANC_A,ANC_B",
  "--pcs", "2",
  "--assignments", file.path(tmp, "unconfigured_assignments.tsv"),
  "--study-pcs-out", file.path(tmp, "unconfigured_study_pcs.tsv"),
  "--distances", file.path(tmp, "unconfigured_distances.tsv"),
  "--reference-outliers", file.path(tmp, "unconfigured_outliers.tsv"),
  "--model-summary", file.path(tmp, "unconfigured_models.tsv"),
  "--within-pcs", file.path(tmp, "unconfigured_within.tsv"),
  "--excluded", file.path(tmp, "unconfigured_excluded.tsv"),
  "--counts", file.path(tmp, "unconfigured_counts.tsv"),
  "--min-confidence", "0.05",
  "--outlier-sd", "4",
  "--min-reference-n", "10"
))
stopifnot(status == 0L)
excluded <- read.delim(file.path(tmp, "unconfigured_excluded.tsv"), sep = "\t", stringsAsFactors = FALSE)
stopifnot(nrow(excluded) == 1)
stopifnot(identical(excluded$reason[[1]], "unconfigured_nearest_super_population"))


# Ambiguity between two populations from the same super-population should not
# exclude a sample when that super-population is clearly separated from others.
same_super_reference <- file.path(tmp, "same_super_reference.tsv")
same_super <- data.frame(FID = character(), IID = character(), population = character(),
  super_population = character(), PC1 = character(), PC2 = character())
for (i in 0:11) {
  offset <- (i - 5.5) * 0.01
  same_super <- rbind(same_super,
    data.frame(FID = paste0("E1_", i), IID = paste0("E1_", i), population = "EUR_1", super_population = "EUR",
      PC1 = as.character(offset), PC2 = as.character(offset)),
    data.frame(FID = paste0("E2_", i), IID = paste0("E2_", i), population = "EUR_2", super_population = "EUR",
      PC1 = as.character(0.2 + offset), PC2 = as.character(0.2 + offset)),
    data.frame(FID = paste0("A_", i), IID = paste0("A_", i), population = "AFR_1", super_population = "AFR",
      PC1 = as.character(5 + offset), PC2 = as.character(5 + offset))
  )
}
write.table(same_super, same_super_reference, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(FID = "SE", IID = "SE", PC1 = "0.1", PC2 = "0.1"),
  file.path(tmp, "same_super_study.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
status <- system2("Rscript", c(
  "scripts/infer_popmad_ancestry.R",
  "--study-pcs", file.path(tmp, "same_super_study.tsv"),
  "--reference-pcs", same_super_reference,
  "--ancestries", "EUR,AFR",
  "--pcs", "2",
  "--assignments", file.path(tmp, "same_super_assignments.tsv"),
  "--study-pcs-out", file.path(tmp, "same_super_study_pcs.tsv"),
  "--distances", file.path(tmp, "same_super_distances.tsv"),
  "--reference-outliers", file.path(tmp, "same_super_outliers.tsv"),
  "--model-summary", file.path(tmp, "same_super_models.tsv"),
  "--within-pcs", file.path(tmp, "same_super_within.tsv"),
  "--excluded", file.path(tmp, "same_super_excluded.tsv"),
  "--counts", file.path(tmp, "same_super_counts.tsv"),
  "--min-confidence", "0.05",
  "--outlier-sd", "4",
  "--min-reference-n", "10"
))
stopifnot(status == 0L)
assigned <- read.delim(file.path(tmp, "same_super_assignments.tsv"), sep = "\t", stringsAsFactors = FALSE)
stopifnot(nrow(assigned) == 1)
stopifnot(identical(assigned$ancestry[[1]], "EUR"))


# Report test completion.
cat("POP-MaD ancestry tests passed\n")
