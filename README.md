# Stage 1 Ancestry-Stratified GWAS

This repository runs a Snakemake pipeline for Stage 1 ancestry-stratified GWAS.
It validates cohort inputs, infers the genotype genome build, resolves a
fingerprinted ancestry reference package, prepares ancestry/QC sample sets, runs
PLINK2 `--glm` within each ancestry stratum, and writes harmonized summary
statistics, plots, reports, and a run manifest.

For real cohort data, follow the production path. Do not run the HapMap3
example-data commands for production analyses.

## Choose Your Path

| Goal | Start here |
| --- | --- |
| Run a real cohort on SLURM | [Production Quickstart](#production-quickstart), then [docs/cohort-production-run-manual.md](docs/cohort-production-run-manual.md) |
| Run focused development checks | [Development Checks](#development-checks), then [docs/development-notes.md](docs/development-notes.md) |
| Understand or edit `config/config.yaml` | [docs/configuration.md](docs/configuration.md) |
| Review pipeline stages and output layout | [docs/pipeline-overview.md](docs/pipeline-overview.md) |
| Understand the ancestry reference package | [docs/ancestry-reference-prep.md](docs/ancestry-reference-prep.md) |
| Record software and reference provenance | [docs/resources-and-downloads.md](docs/resources-and-downloads.md) |

## What Stage 1 Does

- Accepts one genome-wide PLINK dataset: `PGEN/PVAR/PSAM` or `BED/BIM/FAM`.
- Validates config, sample manifest, trait registry, software manifests,
  reference manifests, genotype files, the ancestry reference package, and
  covariates.
- Infers the genotype build from marker positions and writes build-resolved
  config snapshots.
- Computes POP-MaD ancestry from a prebuilt reference package in production.
- Optionally runs supervised ADMIXTURE as report-only ancestry QC.
- Runs genetic sex checks, relatedness filtering, trait/covariate preparation,
  ancestry-stratified PLINK2 GWAS, plotting, reporting, and manifest generation.

Stage 1 does not run imputation, pooled GWAS, METAL, trans-ancestry
meta-analysis, or liftover. It also does not build, download, or modify the
production ancestry reference package.

## Production Quickstart

Use this checklist for real cohort data on a Linux/HPC system. The full runbook
with review points and troubleshooting is in
[docs/cohort-production-run-manual.md](docs/cohort-production-run-manual.md).

### 1. Set Up Environments

Load the site Conda/Mamba module, then create the driver environment for
Snakemake and the utility environment for standalone R checks.

```bash
<conda-or-mamba> env create -f envs/snakemake-driver.yaml
<conda-or-mamba> env create -f envs/gwas.yaml
<activate-command> gwas-stage1-driver
```

Verify the driver environment:

```bash
python -c "import snakemake_executor_plugin_slurm"
snakemake --version
```

### 2. Prepare Inputs

Production requires:

- one `pgen` or `bed` PLINK genotype prefix
- a sample manifest TSV
- a trait registry TSV
- the unpacked prebuilt reference package
- Linux-compatible PLINK2 and ADMIXTURE executables
- reviewed software and reference manifests

Input schemas are documented in [docs/configuration.md](docs/configuration.md).

### 3. Create `config/config.yaml`

```bash
cp config/config.template.yaml config/config.yaml
```

Set these fields first:

```yaml
project:
  analysis_name: "cohort_stage1_gwas"
  cohort_data_release: "cohort_freeze_or_release_label"
  genome_build: "auto"

reference_package:
  root: "/path/to/unpacked/stage1_reference_package"
  fingerprint: "value-from-content_fingerprint.sha256"

inputs:
  sample_manifest: "/path/to/sample_manifest.tsv"
  trait_registry: "/path/to/trait_registry.tsv"

genotypes:
  type: "pgen"
  prefix: "/path/to/study/genotypes_without_extension"
```

Production should normally keep:

```yaml
ancestry_reference:
  enabled: true

admixture:
  enabled: true

sex_check:
  enabled: true
  action: "exclude"
  allow_no_sex_markers: false

gwas:
  allow_missing_pcs: false
```

The reference fingerprint is the content fingerprint from the unpacked package,
not the checksum of a `.tar.gz` archive.

```bash
cat /path/to/unpacked/stage1_reference_package/content_fingerprint.sha256
```

### 4. Configure SLURM

Edit `profiles/slurm/config.yaml` for the cluster account, partition, QoS,
memory/runtime policy, job limits, and shared `conda-prefix`.

### 5. Run Preflight, Dry-Run, And Validation

```bash
<conda-or-mamba> run -n gwas-stage1 Rscript scripts/production_preflight.R \
  --config config/config.yaml \
  --profile profiles/slurm/config.yaml
```

```bash
snakemake -n --profile profiles/slurm
```

```bash
snakemake --profile profiles/slurm results/qc/input_validation/validation.ok
```

Production validation must run through Snakemake because the pipeline first
infers the genome build and writes `results/config/resolved_config.yaml`.

### 6. Run The Pipeline

If compute nodes cannot reach conda channels, pre-create rule environments on a
login/build node:

```bash
snakemake --profile profiles/slurm --conda-create-envs-only
```

Then run the full workflow:

```bash
snakemake --profile profiles/slurm
```

### 7. Review And Archive Results

Start review with:

```text
results/config/resolved_config.yaml
results/manifests/run_manifest.tsv
results/qc/strata/strata_counts.tsv
results/qc/sex/sex_check_summary.tsv
results/qc/relatedness/relatedness_summary.tsv
results/qc/ancestry/reference/reference_prep_report.md
results/qc/ancestry/production/popmad_population_counts.tsv
results/qc/admixture/admixture_report.md
results/reports/
```

Archive the final config, resolved config, manifest, QC reports, GWAS summary
statistics, and plots according to cohort policy.

## Development Checks

The production pipeline requires a real build-matched reference package and does
not support a separate `test` run mode. The old HapMap3 fixture remains useful
for focused script checks and marker-panel development, but it is no longer an
end-to-end pipeline example.

```bash
<conda-or-mamba> env create -f envs/snakemake-driver.yaml
<conda-or-mamba> env create -f envs/gwas.yaml
<activate-command> gwas-stage1-driver
snakemake -n --use-conda  # requires a production-style config/config.yaml
<conda-or-mamba> run -n gwas-stage1 Rscript scripts/test_infer_genome_build.R
<conda-or-mamba> run -n gwas-stage1 Rscript scripts/test_reference_package.R
```

Use [docs/development-notes.md](docs/development-notes.md) for development
checks. Use [docs/cohort-production-run-manual.md](docs/cohort-production-run-manual.md)
for real end-to-end runs.

## Required Inputs At A Glance

The sample manifest must include unique `FID`/`IID` rows, `age`, `age2`, `sex`,
all phenotype columns, and all non-PC covariates used by configured traits. Sex
codes must be `1`, `2`, `0`, `NA`, `-9`, or `.`.

Calculate `age2` as a centered quadratic term, not raw `age^2`:

```text
mean_age = mean(age) across non-missing analysis samples
age2 = (age - mean_age)^2
```

Use the same age units as `age`, usually years, and document the `mean_age`
value used for the cohort release.

The trait registry must include:

```text
trait_id	phenotype_column	case_value	control_value	missing_values
```

An optional `covariates` column can list trait-specific covariates as
comma-separated names.

Production ancestry is computed from the prebuilt reference package. Set only
`reference_package.root` and `reference_package.fingerprint`; package-derived
reference paths are resolved after genome-build inference.

## Main Outputs

```text
results/config/effective_config.yaml
results/config/resolved_config.yaml
results/qc/
results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.plink2.glm.tsv
results/plots/{trait}/{ancestry}/{trait}.{ancestry}.{build}.qq.png
results/plots/{trait}/{ancestry}/{trait}.{ancestry}.{build}.manhattan.png
results/plots/{trait}/{ancestry}/{trait}.{ancestry}.{build}.manhattan.pdf
results/reports/{trait}/{trait}.{ancestry}.{build}.report.md
results/manifests/run_manifest.tsv
```

## Repository Map

```text
Snakefile                         Workflow entry point
workflow/rules/                   Snakemake rules
workflow/snake_helpers.py         Shared DAG helpers
scripts/                          R and shell helper scripts
scripts/lib/                      Shared R helpers
config/                           Config templates and public fixture inputs
envs/                             Conda environment definitions
profiles/slurm/                   SLURM execution profile
resources/                        Marker panels and provenance manifests
docs/                             Stable user documentation
results/                          Generated outputs
```

`config/config.yaml` is intentionally untracked because production configs can
contain private paths. Start production configs from
`config/config.template.yaml`.

## Privacy And Security

Do not commit protected cohort data, credentials, or private external-drive
paths. Keep production inputs and `results/` on approved storage, and review
logs/manifests before sharing outputs because they intentionally record runtime
paths and provenance.
