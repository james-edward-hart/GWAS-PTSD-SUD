# GWAS PTSD-SUD Stage 1 Pipeline

This repository contains a Snakemake workflow for Stage 1 ancestry-stratified GWAS.

The first implementation is intentionally minimal and readable. Helper scripts are small, checks fail early, and every required software/reference resource is tracked in manifests.

## What Stage 1 Does

- Validates config, phenotype/covariate manifests, trait registry, ancestry labels, and genotype inputs.
- Writes the active merged Snakemake config to `results/config/effective_config.yaml`, then writes the build-resolved runtime config to `results/config/resolved_config.yaml`.
- Infers genome build from an offline GRCh36/GRCh37/GRCh38 marker panel before naming outputs.
- Resolves and validates one fingerprinted reference package after build inference for production runs.
- Optionally prepares production computed ancestry from a reference panel: study/reference harmonization, long-range LD exclusion, LD pruning, reference PCA, study projection, POP-MaD assignment, and within-ancestry PCA.
- Optionally runs supervised K=5 ADMIXTURE as report-only QC with a 1000 Genomes reference; these proportions do not alter POP-MaD labels, strata, keep files, or GWAS covariates.
- Runs genetic sex checks and can warn, fail, or exclude mismatches based on config.
- Prepares one phenotype file and one covariate file per configured trait.
- Builds ancestry-stratum keep files and excludes ambiguous POP-MaD assignments.
- Keeps unrelated samples only using the configured KING threshold after QC and LD pruning of relatedness markers.
- Runs PLINK2 `--glm` within each ancestry stratum.
- Harmonizes PLINK2 output into consistent summary statistics.
- Produces QQ and Manhattan plots.
- Writes a Markdown QC report for each trait-by-ancestry GWAS.

## What Stage 1 Does Not Do

- It does not run imputation.
- It does not run pooled GWAS.
- It does not run METAL or trans-ancestry meta-analysis.
- It does not lift summary statistics to GRCh38.
- It does not hide reference-data preparation inside the GWAS rule.

## Quick Start

For real cohort runs on a cluster, use the production runbook first:

```text
docs/cohort-production-run-manual.md
```

Install Snakemake in a Linux or HPC-compatible environment.

```bash
# Create and enter the local test environment.
mamba create -n gwas-stage1 -c conda-forge -c bioconda snakemake plink2 r-base r-yaml r-jsonlite
mamba activate gwas-stage1
```

Download the public HapMap3 example data, then prepare the local fixture:

```bash
# Download public HapMap3 data and build the toy fixture.
bash scripts/download_test_data.sh
Rscript scripts/prepare_hapmap3_fixture.R
```

Run a dry run:

```bash
# Dry-run the workflow without creating outputs.
snakemake -n --use-conda
```

Run locally:

```bash
# Run the local workflow with four cores.
snakemake --cores 4 --shared-fs-usage input-output persistence software-deployment sources storage-local-copies
```

On this macOS test workstation, `config/config.yaml` points `tools.plink2` to `software/local/plink2`. On HPC, change this to a Linux PLINK2 binary, module path, or simply `plink2` if it is already on `PATH`.

Genome build is inferred from genotype marker positions before build-labelled output filenames are expanded. The marker-panel methods are documented in:

```text
# Marker-panel documentation.
resources/README.md
```

Check the example-run outputs:

```bash
# Validate the example outputs.
Rscript scripts/test_pipeline_outputs.R
```

Run on SLURM:

```bash
# Submit/run with the bundled SLURM profile.
snakemake --profile profiles/slurm
```

For production SLURM runs, create a driver environment with the executor plugin:

```bash
# Create the Snakemake driver environment for SLURM submission.
mamba env create -f envs/snakemake-driver.yaml
```

Tip for internet-insulated clusters: Conda needs channel access only while
creating or updating environments. If compute nodes cannot reach the internet,
create the driver environment and pre-build the workflow rule environments on a
login/build node or through your cluster's Conda mirror before production runs:

```bash
mamba env create -f envs/snakemake-driver.yaml
mamba activate gwas-stage1-driver
snakemake --profile profiles/slurm --conda-create-envs-only
```

After those environments exist under the SLURM profile's `conda-prefix`,
ordinary workflow jobs can reuse them without internet access unless the
environment YAML files change.

## Required Inputs

Edit `config/config.yaml` before running real data, or start from the production-oriented template:

```text
# Production config template.
config/config.template.yaml
```

Core inputs:

- `sample_manifest`: strict TSV with `FID`, `IID`, phenotype columns, `age`, `age2`, and `sex`.
- `trait_registry`: strict TSV defining trait IDs, phenotype columns, case/control values, and missing values.
- `ancestry_file`: strict TSV with `FID`, `IID`, and `ancestry` when using precomputed ancestry.
- `pcs_file`: strict TSV with `FID`, `IID`, and `PC1` through `PC10`.
- `genotypes`: one genome-wide PLINK2 `PGEN/PVAR/PSAM` dataset or one genome-wide PLINK1 `BED/BIM/FAM` dataset.

VCF/BCF is not accepted directly. Convert VCF/BCF to PGEN upstream.

Computed ancestry reference preparation is documented in:

```text
# Computed ancestry documentation.
docs/ancestry-reference-prep.md
```

The default example uses prebuilt toy PC tables. For production computed ancestry, set:

```yaml
# Production computed ancestry toggle; enabled options: true, false.
ancestry_reference:
  enabled: true
```

Then configure only `reference_package.root` and `reference_package.fingerprint`. Stage 1 consumes the prebuilt package, validates its fingerprint, and resolves the build-matched HGDP+1KG-compatible POP-MaD panel from the package manifest. It does not build or modify the reference package.

ADMIXTURE QC is disabled by default. To enable the report-only branch, set:

```yaml
# Report-only ADMIXTURE QC toggle; enabled options: true, false.
admixture:
  enabled: true
```

Then configure the same prebuilt reference package. Stage 1 resolves the build-matched 1000 Genomes ADMIXTURE panel from the package manifest.

Before a production SLURM pilot, run:

```bash
Rscript scripts/production_preflight.R --config config/config.yaml --profile profiles/slurm/config.yaml
```

## Software And Reference Manifests

Software is tracked in:

```text
# Software manifest path.
resources/manifests/software.tsv
```

Reference data is tracked in:

```text
# Reference manifest path.
resources/manifests/reference_data.tsv
```

Files downloaded by the helper scripts are tracked in:

```text
# Helper-download manifest path.
resources/manifests/downloaded_files.tsv
```

For every real analysis, record:

- tool name and version
- installation method
- source URL
- download date
- local path
- checksum, when available
- preprocessing steps

## HPC Notes

The workflow is developed for Linux/HPC first. The SLURM profile keeps scheduler details out of the workflow rules so the pipeline can be ported to another scheduler later.

Large production runs should use:

- shared read-only reference-data directories
- per-project writable `results/`
- conda or module environments pinned by version
- Snakemake dry runs before submission
- explicit resource requests per rule

## Output Layout

```text
# Main output directories.
results/
  gwas/
  logs/
  manifests/
  plots/
  qc/
  reports/
```

The main per-analysis files are:

- `results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.plink2.glm.tsv`
- `results/plots/{trait}/{ancestry}/`
- `results/reports/{trait}/{trait}.{ancestry}.{build}.report.md`
- `results/manifests/run_manifest.tsv`
