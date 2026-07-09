#!/usr/bin/env Rscript

# Build PLINK phenotype and covariate files for one configured trait.

cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))


# Parse the trait-specific file paths supplied by Snakemake.
args <- parse_args(repeated = "keep")
require_args(args, c("config", "trait", "pcs-file", "pheno-out", "covar-out"))


# Load config, sample metadata, PCs, and trait registry rows.
config <- load_config(args$config)
samples <- read_tsv(config$inputs$sample_manifest)
pcs <- read_tsv(args[["pcs-file"]])
traits <- read_tsv(config$inputs$trait_registry)


# Select the requested trait and resolve its covariate set.
trait <- traits[traits$trait_id == args$trait, , drop = FALSE]
if (!nrow(trait)) die("unknown trait_id: ", args$trait)

covars <- covariates_for_trait(config, args$trait)
allow_missing_pcs <- truthy(config$gwas$allow_missing_pcs %||% FALSE)

missing_values <- split_csv(trait$missing_values[[1]])
phenotype_column <- trait$phenotype_column[[1]]
case_value <- trait$case_value[[1]]
control_value <- trait$control_value[[1]]
is_binary_trait <- !blank(case_value) && !blank(control_value)


# PLINK case/control coding is 2=case, 1=control, NA=missing. Quantitative
# traits keep their numeric value and use NA for configured missing values.
value <- samples[[phenotype_column]]
if (is_binary_trait) {
  pheno <- ifelse(value == case_value, "2",
    ifelse(value == control_value, "1",
      ifelse(value %in% missing_values, "NA", NA_character_)
    )
  )
  bad <- which(is.na(pheno))
  if (length(bad)) {
    row <- samples[bad[[1]], ]
    die("unexpected phenotype value '", row[[phenotype_column]], "' for ", row$FID, " ", row$IID, " in trait ", args$trait)
  }
} else {
  if (!blank(case_value) || !blank(control_value)) {
    die("trait ", args$trait, " must set both case_value and control_value for binary analysis, or leave both blank for quantitative analysis")
  }
  pheno <- ifelse(value %in% c(missing_values, "", "NA", "-9", "."), "NA", value)
  numeric_pheno <- suppressWarnings(as.numeric(pheno))
  bad <- pheno != "NA" & (is.na(numeric_pheno) | !is.finite(numeric_pheno))
  if (any(bad)) {
    row <- samples[which(bad)[[1]], ]
    die("nonnumeric quantitative phenotype value '", row[[phenotype_column]], "' for ", row$FID, " ", row$IID, " in trait ", args$trait)
  }
}

# Write one phenotype row per sample.
write_tsv(data.frame(FID = samples$FID, IID = samples$IID, PHENO = pheno), args[["pheno-out"]])


# Match PC rows back to the manifest order.
pc_index <- match_sample_rows(samples[c("FID", "IID")], sample_key_map(pcs[c("FID", "IID")], "within-ancestry PC table"))


# PC completeness is enforced only for samples retained by final keep files.
keep_paths <- args$keep %||% character()
if (length(keep_paths)) {
  keep_rows <- data.frame(FID = character(), IID = character())
  for (path in keep_paths) {
    keep <- read_tsv(path)
    require_columns(keep, c("FID", "IID"), paste("keep file", path))
    keep_rows <- rbind(keep_rows, keep[c("FID", "IID")])
  }
  required_for_gwas <- !is.na(match_sample_rows(samples[c("FID", "IID")], sample_key_map(keep_rows, "GWAS keep files")))
} else {
  required_for_gwas <- rep(TRUE, nrow(samples))
}


# Build covariates from manifest columns or projected PC columns.
rows <- data.frame(FID = samples$FID, IID = samples$IID, stringsAsFactors = FALSE)
for (covar in covars) {
  if (covar %in% names(samples)) {
    rows[[covar]] <- samples[[covar]]
  } else if (covar %in% names(pcs)) {
    missing_required <- required_for_gwas & is.na(pc_index)
    if (any(missing_required) && !allow_missing_pcs) {
      first <- which(missing_required)[[1]]
      die("missing PCs for ", samples$FID[[first]], " ", samples$IID[[first]])
    }
    rows[[covar]] <- ifelse(is.na(pc_index), "NA", pcs[[covar]][pc_index])
    numeric_value <- suppressWarnings(as.numeric(rows[[covar]]))
    bad_required <- required_for_gwas & (is.na(numeric_value) | !is.finite(numeric_value))
    if (any(bad_required) && !allow_missing_pcs) {
      first <- which(bad_required)[[1]]
      die("missing or non-finite ", covar, " for kept sample ", samples$FID[[first]], " ", samples$IID[[first]])
    }
  } else if (allow_missing_pcs && startsWith(covar, "PC")) {
    rows[[covar]] <- "NA"
  } else if (startsWith(covar, "PC")) {
    first <- if (any(required_for_gwas)) which(required_for_gwas)[[1]] else 1
    die("required PC covariate '", covar, "' is missing for kept sample ",
      samples$FID[[first]], " ", samples$IID[[first]])
  } else {
    die("covariate '", covar, "' is missing for ", samples$FID[[1]], " ", samples$IID[[1]])
  }
}


# Save the final PLINK covariate table.
write_tsv(rows, args[["covar-out"]])
cat("Wrote trait files for", args$trait, "\n")
