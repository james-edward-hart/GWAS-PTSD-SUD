#!/usr/bin/env Rscript

# Real-PLINK checks for combined autosome/X analysis and downstream reporting.
cmd <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])))
source(file.path(script_dir, "lib", "stage1.R"))
repo <- normalizePath(file.path(script_dir, ".."))
plink <- c(file.path(repo, "software", "local", "plink2"), Sys.which("plink2"))
plink <- plink[nzchar(plink) & file.exists(plink) & file.access(plink, 1) == 0]
if (!length(plink)) stop("PLINK2 is required for the Stage 1 X integration test")
plink <- plink[[1]]
if (!nzchar(Sys.getenv("GAWK_BIN")) && !nzchar(Sys.which("gawk"))) Sys.setenv(GAWK_BIN = Sys.which("awk"))
args <- parse_args(defaults = list("out-dir" = ""))
tmp <- if (blank(args[["out-dir"]])) tempfile("plink2-gwas-x-") else normalizePath(args[["out-dir"]], mustWork = FALSE)
if (dir.exists(tmp)) stop("test output directory already exists: ", tmp)
dir.create(tmp, recursive = TRUE)
if (blank(args[["out-dir"]])) on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

invoke <- function(tool, arguments, error = "") {
  log <- suppressWarnings(system2(tool, vapply(arguments, shQuote, character(1)), stdout = TRUE, stderr = TRUE))
  status <- as.integer(attr(log, "status") %||% 0L)
  if (nzchar(error)) {
    stopifnot(status != 0L, any(grepl(error, log, fixed = TRUE)))
  } else if (status != 0L) {
    cat(paste(tail(log, 20), collapse = "\n"), "\n")
    stop("command failed: ", basename(tool), call. = FALSE)
  }
  invisible(log)
}
task <- function(script, arguments, ...) invoke("Rscript", c(file.path(script_dir, script), arguments), ...)
metric <- function(path, name) {
  rows <- read_tsv(path)
  value <- rows$value[rows$metric == name]
  stopifnot(length(value) == 1L)
  as.character(value)
}

# Independent genotypes avoid covariate collinearity; both sexes occur in both
# phenotype classes. PAR boundary markers exercise inclusive endpoints.
set.seed(420)
samples <- data.frame(FID = "0", IID = sprintf("I%03d", 1:160), sex = rep(1:2, 80),
  age = rnorm(160, 45, 8), binary = rep(c(1, 1, 0, 0), 40), quantitative = rnorm(160))
manifest <- file.path(tmp, "samples.tsv")
write_tsv(samples, manifest)
psam <- file.path(tmp, "samples.psam")
write_tsv(setNames(samples[c("FID", "IID", "sex")], c("#FID", "IID", "SEX")), psam)
keep <- file.path(tmp, "keep.tsv")
write_tsv(samples[c("FID", "IID")], keep)
covar <- file.path(tmp, "covar.tsv")
write_tsv(samples[c("FID", "IID", "age", "sex")], covar)
traits <- file.path(tmp, "traits.tsv")
write_tsv(data.frame(trait_id = c("binary", "quantitative"), phenotype_column = c("binary", "quantitative"),
  case_value = c("1", ""), control_value = c("0", ""), missing_values = "NA", covariates = ""), traits)
phenos <- setNames(file.path(tmp, paste0(c("binary", "quantitative"), ".pheno.tsv")), c("binary", "quantitative"))
for (trait in names(phenos)) {
  value <- if (trait == "binary") samples$binary + 1 else samples$quantitative
  write_tsv(data.frame(FID = samples$FID, IID = samples$IID, PHENO = value), phenos[[trait]])
}

make_source <- function(build, mono_x = FALSE) {
  stem <- file.path(tmp, paste0(build, if (mono_x) "_mono_x" else ""))
  bounds <- if (build == "GRCh37") c(2699520, 154931044) else c(2781479, 155701383)
  variants <- data.frame(chrom = c("1", "1", "22", rep("X", 7), "Y", "MT"),
    pos = c(100000000, 200000000, 50000000, 100000, bounds[1], bounds[1] + 1, 10000000,
      bounds[2] - 1, bounds[2], bounds[2] + 100000, 10000000, 1000),
    id = c("auto1", "auto2", "auto22", "par1", "par1_edge", "nonpar_head", "nonpar",
      "nonpar_tail", "par2_edge", "par2", "y", "mt"))
  header <- c("##fileformat=VCFv4.2", "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples$IID), collapse = "\t"))
  lines <- vapply(seq_len(nrow(variants)), function(i) {
    row <- variants[i, ]
    nonpar <- row$chrom == "X" && row$pos > bounds[1] && row$pos < bounds[2]
    haploid <- (nonpar & samples$sex == 1) | row$chrom %in% c("Y", "MT")
    dose <- rbinom(nrow(samples), ifelse(haploid, 1, 2), 0.3)
    if (mono_x && row$chrom == "X") dose[] <- 0
    gt <- ifelse(haploid, as.character(dose), c("0/0", "0/1", "1/1")[dose + 1])
    if (row$chrom == "Y") gt[samples$sex == 2] <- "."
    if (row$id == "par1") gt[1] <- "."
    paste(c(row$chrom, format(row$pos, scientific = FALSE, trim = TRUE), row$id,
      "A", "G", ".", "PASS", ".", "GT", gt), collapse = "\t")
  }, character(1))
  writeLines(c(header, lines), paste0(stem, ".vcf"))
  invoke(plink, c("--vcf", paste0(stem, ".vcf"), "--psam", psam, "--split-par",
    if (build == "GRCh37") "b37" else "b38", "--make-pgen", "--out", stem))
  list(type = "pgen", prefix = stem)
}
convert <- function(source, name, kind, flags = character()) {
  prefix <- file.path(tmp, name)
  invoke(plink, c(genotype_args(source), flags, if (kind == "bed") "--make-bed" else "--make-pgen", "--out", prefix))
  list(type = kind, prefix = prefix)
}

source37 <- make_source("GRCh37")
source38 <- make_source("GRCh38")
unsplit37 <- convert(source37, "unsplit37", "bed", "--merge-par")
unsplit38 <- convert(source38, "unsplit38", "pgen", "--merge-par")
split38 <- convert(source38, "split38", "bed")
legacy <- convert(source37, "legacy_unsorted", "bed", "--merge-par")
bim <- read.table(paste0(legacy$prefix, ".bim"), stringsAsFactors = FALSE)
bim[bim$V2 %in% c("par1", "par1_edge", "par2_edge", "par2"), 1] <- 25
write.table(bim, paste0(legacy$prefix, ".bim"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
legacy <- convert(legacy, "legacy_sorted", "pgen", "--sort-vars")
legacy <- convert(legacy, "legacy", "bed")
no_x <- convert(source37, "no_x", "pgen", "--autosome")
mono_x <- make_source("GRCh37", mono_x = TRUE)

base_config <- list(project = list(analysis_name = "X_integration", inferred_genome_build = "GRCh37"),
  inputs = list(sample_manifest = manifest, trait_registry = traits), tools = list(plink2 = plink),
  gwas = list(default_covariates = c("age", "sex"), covar_variance_standardize = TRUE,
    glm_options = "hide-covar firth-fallback", test = "ADD"),
  qc = list(autosome_only = FALSE, snps_only_acgt = FALSE, maf_min = 0.01, hwe_p_min = 1e-6,
    geno_missing_max = 0.05, sample_missing_max = 1, use_mach_r2_filter = FALSE),
  warnings = list(min_n = 1, min_cases = 1, min_controls = 1), relatedness = list(mode = "all_samples"))

run_case <- function(name, source, build = "GRCh37", autosome_only = FALSE, expect_x = TRUE, split = "") {
  case_dir <- file.path(tmp, name)
  dir.create(case_dir)
  config <- base_config
  config$genotypes <- source
  config$project$inferred_genome_build <- build
  config$qc$autosome_only <- autosome_only
  config_path <- file.path(case_dir, "config.yaml")
  yaml::write_yaml(config, config_path)
  for (trait in names(phenos)) {
    stem <- file.path(case_dir, trait)
    hwe <- paste0(stem, ".hwe")
    stats_path <- paste0(stem, ".tsv")
    task("write_control_hwe_variants.R", c("--config", config_path, "--trait", trait,
      "--pheno", phenos[[trait]], "--keep", keep, "--plink-prefix", hwe, "--out", paste0(hwe, ".snplist")))
    task("run_plink2_gwas.R", c("--config", config_path, "--trait", trait, "--ancestry", "TEST", "--build", build,
      "--pheno", phenos[[trait]], "--covar", covar, "--keep", keep, "--hwe-snplist", paste0(hwe, ".snplist"),
      "--plink-prefix", stem, "--out", stats_path))
    stats <- read_tsv(stats_path)
    p <- suppressWarnings(as.numeric(stats$p))
    stopifnot(all(stats$chrom %in% c("1", "22", "X")), sum(stats$chrom != "X") == 3L,
      sum(stats$chrom == "X") == if (expect_x) 7L else 0L, all(is.finite(p) & p > 0 & p <= 1))
    raw <- read_tsv(Sys.glob(paste0(stem, ".PHENO.glm.*"))[[1]])
    stopifnot(identical(as.character(stats$n), as.character(raw$OBS_CT)),
      identical(as.character(stats$pos), as.character(raw$POS)), identical(stats$variant_id, raw$ID))
    if (expect_x) stopifnot(all(c("par1", "par1_edge", "nonpar_head", "nonpar", "nonpar_tail", "par2_edge", "par2") %in% stats$variant_id))
    native_logs <- lapply(c(paste0(hwe, ".log"), paste0(stem, ".log")), readLines, warn = FALSE)
    for (log in native_logs) {
      if (nzchar(split)) stopifnot(any(grepl(paste("--split-par", split), log, fixed = TRUE)))
      else stopifnot(!any(grepl("--split-par", log, fixed = TRUE)))
    }
    stopifnot(any(grepl("--xchr-model 2", native_logs[[2]], fixed = TRUE)),
      any(grepl("no-x-sex", native_logs[[2]], fixed = TRUE)))
    summary <- paste0(stem, ".summary.tsv")
    task("summarize_gwas_run.R", c("--config", config_path, "--stats", stats_path, "--plink-log", paste0(stem, ".log"),
      "--hwe-snplist", paste0(hwe, ".snplist"), "--hwe-log", paste0(hwe, ".log"), "--pheno", phenos[[trait]],
      "--covar", covar, "--keep", keep, "--out", summary))
    stopifnot(metric(summary, "autosomal_analysis_variants") == "3", metric(summary, "autosomal_valid_p_value_variants") == "3",
      metric(summary, "x_analysis_variants") == as.character(if (expect_x) 7 else 0),
      metric(summary, "x_valid_p_value_variants") == as.character(if (expect_x) 7 else 0),
      metric(summary, "x_sample_count_min") == if (expect_x) "159" else "NA",
      metric(summary, "x_sample_count_max") == if (expect_x) "160" else "NA")
  }
  cat("Passed:", name, "(binary and quantitative)\n")
  invisible(list(dir = case_dir, config = config_path, stats = file.path(case_dir, "binary.tsv"),
    summary = file.path(case_dir, "binary.summary.tsv"), prefix = file.path(case_dir, "binary")))
}

combined <- run_case("bed_unsplit37", unsplit37, split = "b37")
run_case("pgen_split37", source37)
run_case("pgen_unsplit38", unsplit38, build = "GRCh38", split = "b38")
run_case("bed_split38", split38, build = "GRCh38")
run_case("bed_legacy_par", legacy, split = "b37")
autosomes <- run_case("bed_autosome_mode", unsplit37, autosome_only = TRUE, expect_x = FALSE)
absent <- run_case("pgen_no_x", no_x, expect_x = FALSE)
filtered <- run_case("pgen_all_x_filtered", mono_x, expect_x = FALSE)
stopifnot(identical(read_tsv(autosomes$stats), read_tsv(absent$stats)))

# Without a supplied sex covariate, retain PLINK's automatic X sex adjustment.
native_sex <- load_config(combined$config)
native_sex$gwas$default_covariates <- "age"
native_sex_path <- file.path(tmp, "native_sex.yaml")
yaml::write_yaml(native_sex, native_sex_path)
native_prefix <- file.path(tmp, "native_sex")
task("run_plink2_gwas.R", c("--config", native_sex_path, "--trait", "binary", "--ancestry", "TEST", "--build", "GRCh37",
  "--pheno", phenos[["binary"]], "--covar", covar, "--keep", keep, "--hwe-snplist", paste0(combined$prefix, ".hwe.snplist"),
  "--plink-prefix", native_prefix, "--out", paste0(native_prefix, ".tsv")))
stopifnot(!any(grepl("no-x-sex", readLines(paste0(native_prefix, ".log")), fixed = TRUE)),
  sum(read_tsv(paste0(native_prefix, ".tsv"))$chrom == "X") == 7L)

# Unknown build must not silently select incorrect PAR boundaries.
invalid <- base_config
invalid$genotypes <- unsplit37
invalid$project$inferred_genome_build <- ""
stopifnot(inherits(tryCatch(gwas_chromosome_args(invalid), error = identity), "error"))
# Conflicting explicit sex options must fail with an actionable message.
invalid <- load_config(combined$config)
invalid$gwas$glm_options <- "hide-covar sex"
invalid_path <- file.path(tmp, "duplicate_sex.yaml")
yaml::write_yaml(invalid, invalid_path)
task("run_plink2_gwas.R", c("--config", invalid_path, "--trait", "binary", "--ancestry", "TEST", "--build", "GRCh37",
  "--pheno", phenos[["binary"]], "--covar", covar, "--keep", keep, "--hwe-snplist", paste0(combined$prefix, ".hwe.snplist"),
  "--plink-prefix", file.path(tmp, "invalid_model"), "--out", file.path(tmp, "invalid.tsv")),
  error = "sex is already a configured covariate")

# Normalization preserves aliases, coordinates, row count, and missing P values.
aliases <- file.path(tmp, "aliases")
native <- data.frame(CHROM = c("1", "23", "chrX", "PAR1", "PAR2", "25", "chrXY"),
  POS = c(10, 3000000, 4000000, 100000, 155000000, 200000, 200001), ID = paste0("alias", 1:7),
  REF = "A", ALT = "G", A1 = "G", TEST = "ADD", OBS_CT = c(160, 159, 158, 157, 160, 159, 158),
  BETA = 0.1, SE = 0.2, P = c("0.1", "1e-9", "0.2", "0.3", "NA", ".", "0"))
write_tsv(native, paste0(aliases, ".PHENO.glm.linear"))
task("harmonize_plink_glm.R", c("--plink-prefix", aliases, "--trait", "binary", "--ancestry", "TEST",
  "--build", "GRCh37", "--out", paste0(aliases, ".tsv")))
normalized <- read_tsv(paste0(aliases, ".tsv"))
stopifnot(identical(normalized$chrom, c("1", rep("X", 6))), identical(as.character(normalized$pos), as.character(native$POS)),
  identical(as.character(normalized$n), as.character(native$OBS_CT)), identical(normalized$p[5:7], c("NA", "NA", "0")), nrow(normalized) == nrow(native))
task("summarize_gwas_run.R", c("--config", combined$config, "--stats", paste0(aliases, ".tsv"),
  "--plink-log", paste0(combined$prefix, ".log"), "--hwe-snplist", paste0(combined$prefix, ".hwe.snplist"),
  "--hwe-log", paste0(combined$prefix, ".hwe.log"), "--pheno", phenos[["binary"]], "--covar", covar,
  "--keep", keep, "--out", paste0(aliases, ".summary.tsv")))
stopifnot(metric(paste0(aliases, ".summary.tsv"), "x_analysis_variants") == "6",
  metric(paste0(aliases, ".summary.tsv"), "x_valid_p_value_variants") == "3")

# Reports and plots must agree with the same combined file, including zero-X
# cases. Small tables stand in for unrelated ancestry/provenance outputs.
stub <- file.path(tmp, "stub.tsv")
write_tsv(data.frame(metric = "n_removed", value = "0"), stub)
build_details <- file.path(tmp, "build.tsv")
write_tsv(data.frame(build = "GRCh37", selected = "True", checked_markers = 10, matching_markers = 10, match_fraction = 1), build_details)
ancestry_counts <- file.path(tmp, "ancestry.tsv")
write_tsv(data.frame(category = c("assigned_ancestry", "excluded_total"), ancestry = c("TEST", ""), n = c(160, 0)), ancestry_counts)
strata <- file.path(tmp, "strata.tsv")
write_tsv(data.frame(trait_id = "binary", ancestry = "TEST", n = 160, cases = 80, controls = 80,
  active = "True", excluded_reason = "", underpowered = "False"), strata)
alias_case <- list(dir = tmp, config = combined$config, stats = paste0(aliases, ".tsv"),
  summary = paste0(aliases, ".summary.tsv"), prefix = combined$prefix)
for (case in list(combined, autosomes, absent, filtered, alias_case)) {
  stem <- file.path(case$dir, "plot")
  task("plot_gwas.R", c("--stats", case$stats, "--qq", paste0(stem, ".qq.png"),
    "--manhattan", paste0(stem, ".manhattan.png"), "--manhattan-pdf", paste0(stem, ".manhattan.pdf"),
    "--metrics-out", paste0(stem, ".metrics.tsv"), "--top-hits-out", paste0(stem, ".hits.tsv")))
  stopifnot(metric(paste0(stem, ".metrics.tsv"), "valid_p_value_variants") == metric(case$summary, "valid_p_value_variants"),
    metric(paste0(stem, ".metrics.tsv"), "manhattan_points_plotted") == metric(case$summary, "valid_p_value_variants"))
  report <- paste0(stem, ".report.md")
  task("make_report.R", c("--config", case$config, "--trait", "binary", "--ancestry", "TEST", "--build", "GRCh37",
    "--stats", case$stats, "--gwas-summary", case$summary, "--pheno", phenos[["binary"]], "--covar", covar, "--keep", keep,
    "--qq", paste0(stem, ".qq.png"), "--manhattan", paste0(stem, ".manhattan.png"),
    "--manhattan-pdf", paste0(stem, ".manhattan.pdf"), "--popmad-plot", paste0(stem, ".qq.png"),
    "--strata-counts", strata, "--relatedness-summary", stub, "--sex-check-summary", stub, "--genome-build-details", build_details,
    "--ancestry-counts", ancestry_counts, "--admixture-summary", stub, "--admixture-study", stub, "--admixture-comparison", stub,
    "--admixture-report", stub, "--software", stub, "--reference", stub, "--plink-log", paste0(case$prefix, ".log"), "--out", report))
  text <- readLines(report, warn = FALSE)
  stopifnot(any(grepl(paste0("X results including PAR: ", metric(case$summary, "x_analysis_variants")), text, fixed = TRUE)),
    any(grepl("pool valid P values across all included chromosomes", text, fixed = TRUE)),
    any(grepl("MaCH R2 filtering does not evaluate non-PAR X", text, fixed = TRUE)))
  for (suffix in c(".qq.png", ".manhattan.png", ".manhattan.pdf")) stopifnot(file.info(paste0(stem, suffix))$size > 0)
}
cat("Stage 1 PLINK2 X integration tests passed\n")
