# Cohort Production Quickstart

This is the short run manual for launching a production cohort GWAS. Use it
when you already have approved cohort data, a reference package, and cluster
access.

For full parameter details, see the
[Configuration Guide](configuration.md). If a run fails, use the
[Cohort Production Troubleshooting Guide](cohort-production-run-manual.md).

## 1. Clone The Repository

Work on cluster storage that is visible from login and compute nodes.

```bash
cd /path/to/project
git clone https://github.com/james-edward-hart/GWAS-PTSD-SUD.git
cd GWAS-PTSD-SUD
```

Keep protected cohort inputs outside the repository unless your data-use rules
explicitly allow them to be stored there.

## 2. Download The Reference Package

Download the approved Stage 1 reference package from Zenodo:

```bash
mkdir -p reference-data
curl -L \
  "https://zenodo.org/api/records/20615958/files/stage1_reference_package.tar.gz/content" \
  -o reference-data/stage1_reference_package.tar.gz
curl -L \
  "https://zenodo.org/api/records/20615958/files/stage1_reference_package.tar.gz.sha256/content" \
  -o reference-data/stage1_reference_package.tar.gz.sha256
sha256sum reference-data/stage1_reference_package.tar.gz
cat reference-data/stage1_reference_package.tar.gz.sha256
tar -xzf reference-data/stage1_reference_package.tar.gz -C reference-data
cat reference-data/stage1_reference_package/content_fingerprint.sha256
```

Confirm that the printed archive hash matches the hash in the `.sha256` file.
Use the unpacked package directory as `reference_package.root`. Use the printed
fingerprint as `reference_package.fingerprint`.

## 3. Prepare Cohort Inputs

Prepare one PLINK genotype dataset:

```text
PGEN/PVAR/PSAM
or
BED/BIM/FAM
```

Use the genotype prefix without the file extension in `config/config.yaml`.

Prepare a sample manifest TSV (phenotype + covariate file) with one row per
sample:

```text
FID	IID	age	age2	sex	<phenotype columns>	<non-PC covariates>
```

Required rules:

- `FID` and `IID` must match the genotype files exactly.
- `sex` must use PLINK-style codes: `1`, `2`, `0`, `NA`, `-9`, or `.`.
- Include every phenotype column and every non-PC covariate used in the GWAS.
- Export as tab-delimited text, not CSV or Excel.

Calculate `age2` as a centered quadratic age term:

```r
mean_age <- mean(manifest$age, na.rm = TRUE)
manifest$age2 <- (manifest$age - mean_age)^2
```

Prepare a trait registry TSV:

```text
trait_id	phenotype_column	case_value	control_value	missing_values
```

Optional trait-specific covariates can be listed in a `covariates` column as
comma-separated names.

## 4. Configure The Run

Create the editable config:

```bash
cp config/config.template.yaml config/config.yaml
```

Edit `config/config.yaml` first. These fields are the minimum production
starting point:

```yaml
project:
  analysis_name: "cohort_stage1_gwas"
  cohort_data_release: "cohort_freeze_or_release_label"
  genome_build: "auto"

reference_package:
  root: "/path/to/unpacked/stage1_reference_package"
  fingerprint: "value-from-content_fingerprint.sha256"

inputs:
  sample_manifest: "/path/to/cohort/sample_manifest.tsv"
  trait_registry: "/path/to/cohort/trait_registry.tsv"

genotypes:
  type: "pgen"
  prefix: "/path/to/cohort/genotypes_without_extension"

tools:
  plink2: "plink2"
  plink1: "plink"
  admixture: "admixture"

analysis:
  min_stratum_n: 50
  ancestries:
    - AFR
    - AMR
    - EAS
    - EUR
    - SAS
```

Keep these production safety settings unless the analysis plan says otherwise:

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

Edit the SLURM profile:

```text
profiles/slurm/config.yaml
```

Set your cluster account, partition, optional QOS, and Snakemake conda
environment directory:

```yaml
slurm-qos: "normal"
conda-prefix: "/path/to/shared/conda/envs"

default-resources:
  slurm_account: "your_account"
  slurm_partition: "standard"
  mem_mb: 4000
  runtime: 30
```

`conda-prefix` must be a writable directory accessible to compute nodes.

## 5. Set Up The Environment

Load the conda-providing module used on your cluster. On the project cluster:

```bash
module load miniforge3/23.3.1
```

Create the Snakemake driver environment:

```bash
mamba env create -f envs/snakemake-driver.yaml
conda activate gwas-stage1-driver
```

Create the utility R environment used for preflight:

```bash
mamba env create -f envs/gwas.yaml
```

If either environment already exists, update it instead:

```bash
mamba env update -n gwas-stage1-driver -f envs/snakemake-driver.yaml --prune
mamba env update -n gwas-stage1 -f envs/gwas.yaml --prune
```

## 6. Run Basic Checks

Run preflight:

```bash
mamba run -n gwas-stage1 Rscript scripts/production_preflight.R \
  --config config/config.yaml \
  --profile profiles/slurm/config.yaml
```

Run a Snakemake dry run:

```bash
snakemake -n --profile profiles/slurm
```

Run the validation target:

```bash
snakemake --profile profiles/slurm results/qc/input_validation/validation.ok
```

Continue only after these checks pass.

## 7. Run The Pipeline

Submit the full workflow:

```bash
snakemake --profile profiles/slurm
```

If a job fails, inspect the matching file under `results/logs/`, fix the cause,
and rerun the same command. Snakemake will continue from completed outputs.

## 8. Review Outputs

Start with the final GWAS reports. Each report is written per trait and active
ancestry stratum:

```text
results/reports/
```

The report is the primary review file. It includes the run inputs, reference
package fingerprint, ancestry and ADMIXTURE summaries, active/skipped strata,
sample filtering counts, covariates, variant filtering counts, lambda GC, top
association signals, QQ and Manhattan plots, and POP-MaD projection plot.

Use these supporting files when you need the underlying tables:

```text
results/manifests/run_manifest.tsv
results/config/resolved_config.yaml
results/qc/strata/strata_counts.tsv
results/qc/admixture/admixture_report.md
results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.plink2.glm.tsv
results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.gwas_filter_summary.tsv
results/plots/{trait}/{ancestry}/
```

## 9. Compress And Export Results

Create a run archive:

```bash
tar -czf cohort_stage1_gwas_results_YYYY-MM-DD.tar.gz \
  config/config.yaml \
  profiles/slurm/config.yaml \
  resources/manifests/software.tsv \
  resources/manifests/reference_data.tsv \
  results/
```

Export the archive to the approved destination for the cohort:

```bash
rsync -avP cohort_stage1_gwas_results_YYYY-MM-DD.tar.gz /path/to/export/location/
```

Review the archive contents before sharing. Reports, logs, manifests, and
config snapshots may contain private paths or cohort-specific metadata.
