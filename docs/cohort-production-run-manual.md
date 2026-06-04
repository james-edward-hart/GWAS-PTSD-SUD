# Cohort Production Run Manual

This is the step-by-step manual for running the Stage 1 ancestry-stratified GWAS pipeline on real cohort data.

Do not run the test-data setup commands. Do not build the reference package. The reference package is prebuilt; users only set its path and fingerprint.

## 1. Put The Pipeline On The Cluster

Copy or clone the repository into a project directory on the cluster.

```bash
cd /path/to/project
git clone https://github.com/james-edward-hart/GWAS-PTSD-SUD.git
cd GWAS-PTSD-SUD
```

Keep protected cohort data outside the repository when required by your data-use rules.

Use storage visible from login and compute nodes for the repository, input genotypes, prebuilt reference package, conda environments, and `results/`.

## 2. Set Up Snakemake

Use a login/build node for environment setup.

Load the cluster module or shell setup that provides `conda` or `mamba`. Use whichever package manager and activation command your site supports.

```bash
<conda-or-mamba> env create -f envs/snakemake-driver.yaml
<activate-command> gwas-stage1-driver
python -c "import snakemake_executor_plugin_slurm"
```

If `python -c "import snakemake_executor_plugin_slurm"` fails, the SLURM executor plugin is not installed in the active Snakemake environment.

For running the standalone R preflight script, use an existing R installation with `yaml` and `jsonlite`, or create the utility environment:

```bash
<conda-or-mamba> env create -f envs/gwas.yaml
```

Snakemake will create the rule-specific conda environments during the workflow run.

## 3. Prepare The Real Input Files

Prepare one genome-wide PLINK dataset:

- `PGEN/PVAR/PSAM`, or
- `BED/BIM/FAM`

VCF/BCF is not accepted directly. Convert it upstream before using this pipeline.

Prepare a sample manifest TSV with at least:

```text
FID	IID	age	age2	sex
```

Also include every phenotype column and every non-PC covariate used in the GWAS. Sex codes must be `1`, `2`, `0`, `NA`, `-9`, or `.`.

Do not use `M/F` sex codes. `FID/IID` rows must be unique.

Calculate `age2` before running the pipeline as a centered quadratic age term:

```text
mean_age = mean(age) across non-missing analysis samples
age2 = (age - mean_age)^2
```

Use the same age units as `age`, usually years. Do not use raw `age^2` unless
that is an explicitly approved analysis choice, and record the `mean_age` value
with the cohort release notes.

Prepare a trait registry TSV with:

```text
trait_id	phenotype_column	case_value	control_value	missing_values
```

Optional trait-specific covariates can be added in a `covariates` column as a comma-separated list.

Make sure `FID/IID` values match exactly between the sample manifest and genotype files.

## 4. Record Cohort Data Release

The pipeline records real input paths automatically from `config/config.yaml` into:

```text
results/config/resolved_config.yaml
results/manifests/run_manifest.tsv
results/reports/
```

Do not duplicate local input paths, genome build, dates, file sizes, or checksums in `resources/manifests/input_data.tsv`.

Set a clear cohort release label in the config:

```yaml
project:
  cohort_data_release: "cohort_freeze_or_release_label"
```

Optional: if the cohort wants a small free-text input note file, fill in:

```text
resources/manifests/input_data.tsv
```

with only these columns:

```text
file_role	cohort_data_release	notes
```

Leave `resources.input_manifest: ""` unless this file has the required rows below. A header-only `resources/manifests/input_data.tsv` is invalid if configured.

Use these `file_role` values if the optional file is used:

```text
sample_manifest
trait_registry
study_genotype
```

Example:

```text
sample_manifest	cohort_freeze_2026_06	phenotype/sample table from cohort freeze
trait_registry	cohort_freeze_2026_06	trait definitions approved for this run
study_genotype	cohort_freeze_2026_06	genotype prefix is set in config/config.yaml
```

Also review these manifests before sharing final results:

```text
resources/manifests/software.tsv
resources/manifests/reference_data.tsv
```

Record the Linux PLINK2 and ADMIXTURE versions/paths in `software.tsv`. Record cohort-approved reference provenance in `reference_data.tsv`; the prebuilt package fingerprint remains the runtime reference validation.

## 5. Create The Production Config

Start from the template.

```bash
cp config/config.template.yaml config/config.yaml
```

Edit:

```text
config/config.yaml
```

Set these fields first:

```yaml
project:
  analysis_name: "cohort_stage1_gwas"
  cohort_data_release: "cohort_freeze_or_release_label"
  genome_build: "auto"

reference_package:
  root: "/path/to/stage1_reference_package"
  fingerprint: "value-from-content_fingerprint.sha256"

inputs:
  sample_manifest: "/path/to/sample_manifest.tsv"
  trait_registry: "/path/to/trait_registry.tsv"

genotypes:
  type: "pgen"
  prefix: "/path/to/study/genotypes_without_extension"

resources:
  input_manifest: ""
```

Get the reference fingerprint with:

```bash
cat /path/to/stage1_reference_package/content_fingerprint.sha256
```

For this project reference package, the expected fingerprint may be:

```text
532a1e34598aa8c92ca72a8c3983dd06c82f19ea0562c45cd75d31054647976f
```

Use the value provided with your actual package.

Use the unpacked package directory as `reference_package.root`. Do not use the `.tar.gz` archive path, and do not use the archive checksum as the content fingerprint.

Do not manually add or edit `reference_package.observed_fingerprint`; Snakemake writes that into `results/config/resolved_config.yaml`.

The reference package handoff must include:

- `content_fingerprint.sha256`
- `file_manifest.tsv`
- `panel_manifest.tsv`
- package-relative paths only
- exactly one build-matched `popmad` panel for the inferred study build
- exactly one build-matched `admixture` panel for the inferred study build

The resolver rejects missing manifest files, unmanifested package files, absolute paths, `..` paths, file-size or SHA-256 mismatches, raw Hail/VCF/BCF artifacts, and package panels whose required genotype, metadata, or exclusion-region files are absent from `file_manifest.tsv`.

Keep these production settings unless there is a documented reason to change them:

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

If the genotype data are autosome-only, set `sex_check.allow_no_sex_markers: true` only after external sex QC has already been completed and documented.

Set `analysis.ancestries` to the ancestry strata you intend to run. Each configured trait/ancestry cell must have at least one case and one control.

For imputed dosage data with MACH_R2/INFO annotations, set:

```yaml
qc:
  use_mach_r2_filter: true
  info_min: 0.8
```

## 6. Configure The SLURM Profile

Edit:

```text
profiles/slurm/config.yaml
```

Set site-specific values such as:

```yaml
slurm_account: "your_account"
slurm_qos: "normal"
slurm_partition: "standard"
conda-prefix: "/path/to/shared/conda/envs"
```

Adjust memory, runtime, partition, and job limits if your cluster requires different settings.

`runtime` values in the profile are minutes. Put `conda-prefix` on shared writable storage, not node-local scratch.

## 7. Run Preflight Checks

From the repository root:

```bash
<activate-command> gwas-stage1-driver
```

Run the production preflight. Use the R environment available on your cluster, or run through the utility environment:

```bash
<conda-or-mamba> run -n gwas-stage1 Rscript scripts/production_preflight.R \
  --config config/config.yaml \
  --profile profiles/slurm/config.yaml
```

Preflight must pass before submission. It checks the reference package fingerprint, required ancestry settings, optional input-manifest path existence, SLURM profile placeholders, and whether `results/` is clean.

## 8. Dry-Run The Workflow

```bash
snakemake -n --profile profiles/slurm
```

Review the planned jobs. Do not start the real run until the dry-run completes without errors.

## 9. Create Conda Environments

If compute nodes cannot access conda channels, create environments before the main run:

```bash
snakemake --profile profiles/slurm --conda-create-envs-only
```

This step may take a while the first time.

## 10. Run A Validation Pilot

Run the input validation target before launching the full workflow:

```bash
snakemake --profile profiles/slurm results/qc/input_validation/validation.ok
```

Review:

```text
results/logs/validation/
results/qc/genome_build/
results/config/effective_config.yaml
results/config/resolved_config.yaml
```

Do not continue until validation passes.

Do not run `scripts/validate_config.R` directly on `config/config.yaml` for production. Production validation uses the build-resolved config written by Snakemake.

## 11. Run An Ancestry/QC Pilot

Before the full GWAS, it is useful to run the reference projection, POP-MaD, ADMIXTURE, strata, sex-check, and relatedness QC outputs:

```bash
snakemake --profile profiles/slurm \
  results/qc/ancestry/reference/reference_prep_report.md \
  results/qc/admixture/admixture_report.md \
  results/qc/strata/strata_counts.tsv \
  results/qc/ancestry/within_ancestry_pcs.tsv \
  results/qc/sex/sex_check_summary.tsv \
  results/qc/relatedness/relatedness_summary.tsv
```

Review these files before launching all GWAS jobs.

## 12. Run The Full Pipeline

```bash
snakemake --profile profiles/slurm
```

Monitor the scheduler queue and Snakemake logs. If a job fails, fix the cause and rerun the same command. Snakemake will continue from completed outputs.

## 13. Review The Results

Start with these files:

```text
results/manifests/run_manifest.tsv
results/qc/strata/strata_counts.tsv
results/qc/sex/sex_check_summary.tsv
results/qc/relatedness/relatedness_summary.tsv
results/qc/ancestry/reference/reference_prep_report.md
results/qc/ancestry/production/popmad_population_counts.tsv
results/qc/admixture/admixture_report.md
results/reports/
```

Main GWAS outputs are:

```text
results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.plink2.glm.tsv
results/plots/{trait}/{ancestry}/{trait}.{ancestry}.{build}.qq.png
results/plots/{trait}/{ancestry}/{trait}.{ancestry}.{build}.manhattan.png
results/plots/{trait}/{ancestry}/{trait}.{ancestry}.{build}.manhattan.pdf
results/reports/{trait}/{trait}.{ancestry}.{build}.report.md
```

Archive the final `config/config.yaml`, `results/config/resolved_config.yaml`, `results/manifests/run_manifest.tsv`, QC reports, GWAS summary statistics, and plots according to cohort policy.

## Helpful Tips

- If `snakemake` says `invalid choice: 'slurm'`, activate the driver environment and confirm `snakemake_executor_plugin_slurm` is installed.
- If PLINK2 or ADMIXTURE is not found, set `tools.plink2` or `tools.admixture` to a Linux executable path or cluster module shim.
- If the reference fingerprint fails, check `reference_package.root` and `reference_package.fingerprint`. Do not rebuild the package inside this pipeline.
- If the optional input manifest fails, either leave `resources.input_manifest: ""` or use only `file_role`, `cohort_data_release`, and `notes` columns with rows for `sample_manifest`, `trait_registry`, and `study_genotype`.
- If genome-build inference fails, check that the genotype prefix is correct and that BIM/PVAR marker positions match the intended genome build.
- If sex check fails because no X/Y markers exist, either provide genotype data with sex chromosomes or set `sex_check.allow_no_sex_markers: true` only with documented external sex QC.
- If production fails for an empty trait/ancestry cell, remove that ancestry or trait from the config/registry before rerunning.
- If conda fails on compute nodes, create environments on a login/build node with `snakemake --profile profiles/slurm --conda-create-envs-only`.
- If `results/` is not clean before a new production run, move or archive the old directory first.
- For any failed rule, inspect the matching file under `results/logs/` before rerunning.
- Snakemake engine logs are under `.snakemake/log/`; PLINK2 GWAS logs are under `results/gwas/{trait}/{ancestry}/plink2_raw/`.
- If a run is interrupted and Snakemake reports a lock, confirm no Snakemake process is still active, then run `snakemake --unlock --profile profiles/slurm`.
