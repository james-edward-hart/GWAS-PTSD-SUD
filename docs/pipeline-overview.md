# Pipeline Overview

This Snakemake workflow runs Stage 1 ancestry-stratified GWAS and Phase 2 pooled pan-ancestry GWAS. It is designed to keep cohort inputs, reference package validation, QC decisions, association testing, and reports explicit and reproducible.

## What Stage 1 Does

- Validates config, input schemas, trait registry values, ancestry labels, genotype files, software manifests, and reference manifests.
- Writes `results/config/effective_config.yaml` from the merged Snakemake config and `results/config/resolved_config.yaml` after genome-build inference and reference-package resolution.
- Infers the genotype genome build from an offline GRCh36/GRCh37/GRCh38 marker panel before expanding build-labelled outputs.
- Resolves one fingerprinted, unpacked reference package.
- Computes POP-MaD ancestry and within-ancestry GWAS PCs from package-backed reference projection.
- Optionally runs supervised K=5 ADMIXTURE as report-only QC.
- Runs genetic sex checks on a manifest-only, PAR-aware, MAF/missingness-filtered,
  LD-pruned X/Y marker set and can warn, fail, or exclude mismatches based on config.
- Prepares one phenotype file and one covariate file per trait.
- Builds ancestry-stratum keep files, excludes ambiguous/unassigned POP-MaD samples, and keeps unrelated samples using the configured KING threshold.
- Runs PLINK2 `--glm` within each ancestry stratum.
- Harmonizes PLINK2 outputs into consistent summary statistics.
- Runs Phase 2 pooled pan-ancestry regenie GWAS under the `PAN` label when `phase2_regenie.enabled: true`.
- Optionally exports cohort rare-variant HTP statistics, matching marginal
  within-gene ReMeta LD, and a handoff manifest when `remeta.enabled: true`.
- Fits cohort-global PCs for Phase 2 by fitting PCA in sex-QC-passing unrelated samples and projecting all Phase 2 samples.
- Uses Stage 1 ancestry-stratified variant pass lists to build Phase 2 regenie association extract lists.
- Produces a POP-MaD reference/study projection plot, QQ plots with lambda GC, and Manhattan plots as PNG, with Manhattan plots also written as PDF.
- Writes one Markdown QC report per trait-by-ancestry Stage 1 GWAS and one slim Phase 2 `PAN` regenie report per trait.
- Writes a workflow run manifest.

## What Stage 1 Does Not Do

- It does not run imputation.
- It does not run METAL or trans-ancestry meta-analysis.
- It does not run central ReMeta gene tests or merge cohort results.
- It does not lift genotype coordinates or summary statistics.
- It does not download, create, rebuild, or modify the production reference package.

## Repository Map

```text
Snakefile                         Workflow entry point
workflow/rules/                   Snakemake rules
workflow/snake_helpers.py         Shared DAG helper functions
scripts/                          R and shell helper scripts
scripts/lib/                      Shared R helpers
config/                           Templates and public example config inputs
envs/                             Conda environment definitions
profiles/slurm/                   SLURM execution profile
resources/                        Marker panels, manifests, and resource docs
docs/                             User manuals and development notes
results/                          Generated outputs
```

`config/config.yaml` is intentionally ignored because production copies often contain private paths. Use `config/config.template.yaml` for production runs.

## Output Layout

```text
results/
  config/
  gwas/
  logs/
  manifests/
  plots/
  qc/
  reports/
```

Main per-analysis files:

```text
results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.plink2.glm.tsv
results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.gwas_filter_summary.tsv
results/gwas/{trait}/PAN/{trait}.PAN.{build}.regenie
results/gwas/{trait}/PAN/{trait}.PAN.{build}.phase2_summary.tsv
results/plots/ancestry/{analysis_name}.popmad_reference_study_pcs.png
results/plots/{trait}/{ancestry}/{analysis_name}.{trait}.{ancestry}.{build}.qq.png
results/plots/{trait}/{ancestry}/{analysis_name}.{trait}.{ancestry}.{build}.manhattan.png
results/plots/{trait}/{ancestry}/{analysis_name}.{trait}.{ancestry}.{build}.manhattan.pdf
results/plots/{trait}/PAN/{analysis_name}.{trait}.PAN.{build}.regenie.qq.png
results/plots/{trait}/PAN/{analysis_name}.{trait}.PAN.{build}.regenie.manhattan.png
results/reports/{trait}/{analysis_name}.{trait}.{ancestry}.{build}.report.md
results/reports/{trait}/{analysis_name}.{trait}.PAN.{build}.regenie.report.md
results/manifests/run_manifest.tsv
results/remeta/export/{build}/htp/{trait}.PAN.regenie.gz
results/remeta/export/{build}/ld/{group}/chr{chrom}.remeta.gene.ld
results/remeta/export/{analysis_name}.{build}.remeta_manifest.tsv
```

`{analysis_name}` is the filename-safe version of `project.analysis_name`.
The production review and archive sequence is in
[cohort-production-quickstart.md](cohort-production-quickstart.md).
