# Stage 1 Ancestry-Stratified GWAS Instruction Manual

This repository runs Stage 1 ancestry-stratified GWAS with Snakemake. It validates cohort inputs, infers the genotype genome build, resolves a fingerprinted reference package, prepares ancestry/QC sample sets, runs PLINK2 `--glm` within each ancestry stratum, and writes harmonized summary statistics, plots, reports, and run manifests.

For real cohort data, start with the production path below. Do not run the HapMap3 example-data commands for production.

## Start Here

| Use case | Start with |
| --- | --- |
| Production cohort run on SLURM | [docs/cohort-production-run-manual.md](docs/cohort-production-run-manual.md) |
| Local toy/example run | [docs/local-example-run.md](docs/local-example-run.md) |
| Config fields and input schemas | [docs/configuration.md](docs/configuration.md) |
| Ancestry reference package behavior | [docs/ancestry-reference-prep.md](docs/ancestry-reference-prep.md) |
| Software, reference data, manifests | [docs/resources-and-downloads.md](docs/resources-and-downloads.md) |
| Pipeline scope and output overview | [docs/pipeline-overview.md](docs/pipeline-overview.md) |

## Production Run Sequence

Use a Linux/HPC-compatible environment. The SLURM profile requires Snakemake 8+ and the SLURM executor plugin.

```bash
mamba env create -f envs/snakemake-driver.yaml
mamba env create -f envs/gwas.yaml
mamba activate gwas-stage1-driver
python -c "import snakemake_executor_plugin_slurm"
```

Create the production config from the template:

```bash
cp config/config.template.yaml config/config.yaml
```

Edit `config/config.yaml` and set, at minimum:

- `project.analysis_name`
- `project.cohort_data_release`
- `project.run_mode: "production"`
- `reference_package.root`
- `reference_package.fingerprint`
- `inputs.sample_manifest`
- `inputs.trait_registry`
- `inputs.ancestry_mode: "computed"`
- `genotypes.type`
- `genotypes.prefix`
- `tools.plink2`
- `tools.admixture`

Production runs should usually keep:

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

Edit `profiles/slurm/config.yaml` for the cluster account, partition, QoS, and shared conda prefix. The bundled profile contains placeholders and is not cluster-ready until those values match your site.

Run preflight before submission:

```bash
mamba run -n gwas-stage1 Rscript scripts/production_preflight.R \
  --config config/config.yaml \
  --profile profiles/slurm/config.yaml
```

Dry-run the workflow:

```bash
snakemake -n --profile profiles/slurm
```

If compute nodes cannot reach conda channels, pre-create rule environments on a login/build node:

```bash
snakemake --profile profiles/slurm --conda-create-envs-only
```

Run the validation target:

```bash
snakemake --profile profiles/slurm results/qc/input_validation/validation.ok
```

Then run the full pipeline:

```bash
snakemake --profile profiles/slurm
```

Production validation is performed inside Snakemake after genome-build inference and reference-package resolution. Do not use `scripts/validate_config.R --config config/config.yaml` as the production validation step; the production workflow validates `results/config/resolved_config.yaml`.

## Required Inputs

Stage 1 accepts one genome-wide PLINK dataset:

- `PGEN/PVAR/PSAM`, or
- `BED/BIM/FAM`

VCF/BCF is not accepted directly. Convert VCF/BCF upstream.

The sample manifest TSV must include unique `FID`/`IID` rows, `age`, `age2`, `sex`, every phenotype column, and every non-PC covariate used by configured traits. Sex codes must be `1`, `2`, `0`, `NA`, `-9`, or `.`.

The trait registry TSV must include:

```text
trait_id	phenotype_column	case_value	control_value	missing_values
```

Optional trait-specific covariates can be added in a comma-separated `covariates` column.

Production ancestry is computed from the prebuilt, unpacked reference package. Set only `reference_package.root` and `reference_package.fingerprint`; the workflow validates the package fingerprint and resolves the build-matched POP-MaD and ADMIXTURE panels after genome-build inference.

## Main Outputs

The full workflow writes:

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

Review and archive the final config, resolved config, QC reports, GWAS summary statistics, plots, and run manifest according to cohort policy.

## Local Example Smoke Test

The local example uses public HapMap3 genotype data and a generated random binary phenotype. It is for workflow testing only.

```bash
mamba env create -f envs/snakemake-driver.yaml
mamba env create -f envs/gwas.yaml
mamba activate gwas-stage1-driver
cp config/config.hapmap3.example.yaml config/config.yaml
bash scripts/download_test_data.sh
mamba run -n gwas-stage1 Rscript scripts/prepare_hapmap3_fixture.R --plink2 plink2
snakemake -n --use-conda
snakemake --cores 4 --use-conda --shared-fs-usage input-output persistence software-deployment sources storage-local-copies
mamba run -n gwas-stage1 Rscript scripts/test_pipeline_outputs.R
```

If your local config points to `software/local/plink2`, make sure that path exists and matches your platform, or change `tools.plink2` to `plink2`.

## Scope

Stage 1 does not run imputation, pooled GWAS, METAL, trans-ancestry meta-analysis, or liftover. It also does not build or download the production reference package inside the GWAS workflow.
