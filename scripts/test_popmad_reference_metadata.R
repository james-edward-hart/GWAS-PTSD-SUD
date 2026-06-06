#!/usr/bin/env Rscript

# Focused checks for POP-MaD reference-population coverage validation.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))

metadata <- rbind(
  data.frame(sample_id = paste0("A_big_", seq_len(25)), population = "A_big", super_population = "ANC_A"),
  data.frame(sample_id = paste0("A_small_", seq_len(5)), population = "A_small", super_population = "ANC_A"),
  data.frame(sample_id = paste0("B_big_", seq_len(22)), population = "B_big", super_population = "ANC_B")
)
coverage <- popmad_reference_model_coverage(metadata, "population", "super_population", c("ANC_A", "ANC_B"), 20)
stopifnot(!length(coverage$missing_super))
stopifnot(identical(coverage$low$population, "A_small"))
stopifnot(setequal(coverage$retained$population, c("A_big", "B_big")))

low_only <- rbind(
  metadata,
  data.frame(sample_id = paste0("C_small_", seq_len(8)), population = "C_small", super_population = "ANC_C")
)
coverage <- popmad_reference_model_coverage(low_only, "population", "super_population", c("ANC_A", "ANC_B", "ANC_C"), 20)
stopifnot(identical(coverage$missing_super, "ANC_C"))

conflicting <- rbind(
  metadata,
  data.frame(sample_id = "conflict_1", population = "A_big", super_population = "ANC_X")
)
status <- tryCatch({
  popmad_reference_model_coverage(conflicting, "population", "super_population", c("ANC_A"), 20)
  TRUE
}, error = function(e) FALSE)
stopifnot(!status)

cat("POP-MaD reference metadata tests passed\n")
