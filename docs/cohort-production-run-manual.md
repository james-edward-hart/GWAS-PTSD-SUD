# Cohort Production Run Manual

This is the complete operating guide for running the Stage 1
ancestry-stratified GWAS pipeline on real cohort data. You should be able to run
a cohort analysis from this page without switching back to the main README.

Use the production reference package supplied for the analysis. Set its path and
fingerprint in `config/config.yaml`; do not run development-data setup commands
or rebuild the reference package as part of a cohort run.

## 1. Put The Pipeline On The Cluster

### 1.1 Copy Or Clone The Repository

Copy or clone the repository into a project directory on the cluster.

```bash
cd /path/to/project
git clone https://github.com/james-edward-hart/GWAS-PTSD-SUD.git
cd GWAS-PTSD-SUD
```

### 1.2 Keep Protected Data Outside The Repository

Keep protected cohort data outside the repository when required by your data-use
rules. Do not commit cohort data, credentials, or private storage paths.

### 1.3 Use Shared Storage

Use storage visible from login and compute nodes for the repository, input
genotypes, prebuilt reference package, conda environments, and `results/`.
If a path is only visible from the login node, Snakemake may dry-run
successfully and then fail when the SLURM job starts.

## 2. Set Up Snakemake

### 2.1 Use A Login Or Build Node

Use a login/build node for environment setup. Load the cluster module or shell
setup that provides `conda` or `mamba`. The exact module name is
site-specific; common choices are Miniforge, Mambaforge, Miniconda, or Anaconda.
Prefer a current conda/mamba module with conda-forge and bioconda access. Do
not use a Python-only module for environment creation.

On the cluster used for this project, load:

```bash
module load miniforge3/23.3.1
```

### 2.2 Create The Driver Environment

Create and activate the Snakemake driver environment. The repository includes
this environment in `envs/snakemake-driver.yaml`; it installs Snakemake and the
SLURM executor plugin needed by the bundled profile. It also installs a current
`conda`, which Snakemake uses to create rule-specific environments.

```bash
mamba env create -f envs/snakemake-driver.yaml
conda activate gwas-stage1-driver
```

If the environment already exists from an older checkout, update it from the
repository environment file:

```bash
mamba env update -n gwas-stage1-driver -f envs/snakemake-driver.yaml --prune
```

Use the driver environment for Snakemake commands. Do not keep the utility R
environment active at the same time; use `mamba run -n gwas-stage1 ...` for
preflight checks instead.

### 2.3 Verify The SLURM Executor

Confirm that the active environment includes the Snakemake SLURM executor
plugin.

```bash
echo "$CONDA_PREFIX"
"$CONDA_PREFIX/bin/conda" --version
python -c "import conda; print(conda.__version__)"
python -c "import snakemake_executor_plugin_slurm"
snakemake --version
```

`$CONDA_PREFIX` should point to the `gwas-stage1-driver` environment, and
`$CONDA_PREFIX/bin/conda --version` must report `24.7.1` or later. A plain
`conda --version` command may report the module or base conda version instead,
so use the explicit path above for this check. If Snakemake reports
`CreateCondaEnvironmentException: Conda must be version 24.7.1 or later`, update
`gwas-stage1-driver` from `envs/snakemake-driver.yaml` after loading the
intended conda/mamba module.

If the import command fails, the SLURM executor plugin is not installed in the
active Snakemake environment. The most common fix is to load the intended
conda/mamba module, activate `gwas-stage1-driver`, and rebuild or update that
environment from `envs/snakemake-driver.yaml`.

If `mamba env update` does not update the driver environment conda version, run:

```bash
mamba install -n gwas-stage1-driver -c conda-forge "conda>=24.7.1"
```

Snakemake may also warn if conda channel priority is not strict. This warning is
not the same as a failed dry run, but strict priority makes rule environment
creation more reproducible:

```bash
conda config --set channel_priority strict
```

### 2.4 Create The Utility R Environment

Create the utility R environment used by the production preflight script. This
environment is separate from the Snakemake driver environment; it provides R,
`yaml`, `jsonlite`, and the small command-line tools needed for standalone
checks.

```bash
mamba env create -f envs/gwas.yaml
```

If `gwas-stage1` already exists, update it from the repository environment file:

```bash
mamba env update -n gwas-stage1 -f envs/gwas.yaml --prune
```

If your cluster provides `conda` but not `mamba`, use:

```bash
conda env create -f envs/gwas.yaml
```

### 2.5 Let Snakemake Create Rule Environments

Snakemake creates rule-specific conda environments during the workflow run.
Create them ahead of time only if compute nodes cannot access conda channels;
that command is listed in Step 9. The rule environment definitions live in
`envs/`, so users should not manually install rule dependencies into the driver
environment.

## 3. Prepare The Real Input Files

### 3.1 Prepare One PLINK Genotype Dataset

Prepare one genome-wide PLINK dataset:

- `PGEN/PVAR/PSAM`, or
- `BED/BIM/FAM`

VCF/BCF is not accepted directly. Convert it upstream before using this
pipeline.

Set the config genotype prefix without the file extension. For example, use
`/path/to/cohort/genotypes` rather than `/path/to/cohort/genotypes.pgen`.

### 3.2 Prepare The Sample Manifest

Prepare a sample manifest TSV (phenotype + covariate file) with at least:

```text
FID	IID	age	age2	sex
```

**Also include every phenotype column** and every non-PC covariate used in the
GWAS. Sex codes must be `1`, `2`, `0`, `NA`, `-9`, or `.`. Do not use `M/F`
sex codes. `FID/IID` rows must be unique.

Export this file as tab-delimited text, not CSV or Excel. Header names are
matched exactly, and non-PC covariates used in the model should be numeric after
missing values are applied.

### 3.3 Calculate `age2`

Calculate `age2` before running the pipeline as a centered quadratic age term:

```r
mean_age <- mean(manifest$age, na.rm = TRUE)
manifest$age2 <- (manifest$age - mean_age)^2
```

Use the same age units as `age`, usually years. Do not use raw `age^2` unless
that is an explicitly approved analysis choice, and record the `mean_age` value
with the cohort release notes.

### 3.4 Prepare The Trait Registry

Prepare a trait registry TSV with:

```text
trait_id	phenotype_column	case_value	control_value	missing_values
```

Optional trait-specific covariates can be added in a `covariates` column as a
comma-separated list.

### 3.5 Confirm ID Matching

Make sure `FID/IID` values match exactly between the sample manifest and
genotype files before running preflight or Snakemake validation.
Order does not matter, but spelling, leading zeros, and whitespace do.

## 4. Record Cohort Data Release

### 4.1 Set A Cohort Release Label

Set a clear cohort release label in the config.

```yaml
project:
  cohort_data_release: "cohort_freeze_or_release_label"
```

### 4.2 Continue With The Config

After setting `project.cohort_data_release`, continue to Step 5. No additional
input log is required here. When the workflow runs, it automatically records the
resolved config and run-level metadata under `results/`.

### 4.3 Review Software And Reference Manifests

Review these repository provenance files during configuration and again before
sharing final results:

```text
resources/manifests/software.tsv
resources/manifests/reference_data.tsv
```

Record the Linux PLINK2 and ADMIXTURE versions/paths in `software.tsv`. Record
cohort-approved reference provenance in `reference_data.tsv`; the prebuilt
package fingerprint remains the runtime reference validation.

## 5. Create The Production Config

### 5.1 Know Which Files You Configure

For a standard production run, configure these repository files:

| File | Where it lives | What to do |
| --- | --- | --- |
| `config/config.yaml` | Copy from `config/config.template.yaml` in the repository root | Main run config. Set paths to cohort inputs, the reference package, tools, genotype files, ancestries, and QC options. |
| `profiles/slurm/config.yaml` | Bundled SLURM profile in `profiles/slurm/` | Set site scheduler values such as account, partition, QOS, job resources, and `conda-prefix`. |
| `resources/manifests/software.tsv` | Repository provenance file in `resources/manifests/` | Record the PLINK2, ADMIXTURE, Snakemake, and R versions used for the run. |
| `resources/manifests/reference_data.tsv` | Repository provenance file in `resources/manifests/` | Record the ancestry reference package provenance. Runtime validation still uses the package fingerprint in `config/config.yaml`. |

These cohort files do not have to live in the repository. For production, keep
protected cohort inputs on secure cluster-visible storage and point to them from
`config/config.yaml`:

| Cohort input | Config field |
| --- | --- |
| Sample manifest TSV (phenotype + covariate file) | `inputs.sample_manifest` |
| Trait registry TSV | `inputs.trait_registry` |
| Study genotype prefix without extension | `genotypes.prefix` |
| Unpacked ancestry reference package directory | `reference_package.root` |

Most runs should leave `resources.input_manifest: ""`. Only set it to
`resources/manifests/input_data.tsv` if your team wants an additional short
free-text note file about the cohort input release.

### 5.2 Copy The Template

Start from the template.

```bash
cp config/config.template.yaml config/config.yaml
```

Edit:

```text
config/config.yaml
```

Use full paths that are visible on compute nodes. Avoid `~` and shell variables
inside YAML values because they are not a reliable substitute for explicit
cluster-visible paths.

### 5.3 Set Required Fields First

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
  sample_manifest: "/path/to/cohort/sample_manifest.tsv"
  trait_registry: "/path/to/cohort/trait_registry.tsv"

tools:
  plink2: "plink2"
  plink1: "plink"
  admixture: "admixture"

genotypes:
  type: "pgen"
  prefix: "/path/to/cohort/genotypes_without_extension"

resources:
  input_manifest: ""
```

Use `plink2`, `plink1`, and `admixture` only if those commands resolve inside cluster
jobs. Otherwise, set each field to a full Linux executable path or a
site-approved module shim, and record the same tool versions in
`resources/manifests/software.tsv`.

### 5.4 Set The Reference Fingerprint

Get the reference fingerprint from the unpacked package.

```bash
cat /path/to/unpacked/stage1_reference_package/content_fingerprint.sha256
```

For this project reference package, the expected fingerprint may be:

```text
532a1e34598aa8c92ca72a8c3983dd06c82f19ea0562c45cd75d31054647976f
```

Use the value provided with your actual package. Use the unpacked package
directory as `reference_package.root`. Do not use the `.tar.gz` archive path,
and do not use the archive checksum as the content fingerprint.

Do not manually add or edit `reference_package.observed_fingerprint`;
Snakemake writes that into `results/config/resolved_config.yaml`.

### 5.5 Confirm The Reference Package Handoff

The reference package handoff must include:

- `content_fingerprint.sha256`
- `file_manifest.tsv`
- `panel_manifest.tsv`
- package-relative paths only
- exactly one build-matched `popmad` panel for the inferred study build
- exactly one build-matched `admixture` panel for the inferred study build

The resolver rejects missing manifest files, unmanifested package files,
absolute paths, `..` paths, file-size or SHA-256 mismatches, raw Hail/VCF/BCF
artifacts, and package panels whose required genotype, metadata, or
exclusion-region files are absent from `file_manifest.tsv`. macOS sidecar files
such as `._*`, `.DS_Store`, and `__MACOSX/` are ignored because they are copy
metadata, not reference-package content.

### 5.6 Keep Production Safety Settings

Keep these production settings unless there is a documented reason to change
them:

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

When the sex-check threshold fields are blank, the workflow uses chrX defaults
of `max-female-xf=0.2` and `min-male-xf=0.8`. If sex check removes every
sample, review `results/qc/sex/sexcheck.tsv` and adjust thresholds or fix
manifest/genotype sex coding before continuing. Custom thresholds must include
both chrX fields; Y-rate thresholds are optional but must be supplied as a pair.

If the genotype data are autosome-only, set
`sex_check.allow_no_sex_markers: true` only after external sex QC has already
been completed and documented.

### 5.7 Set Analysis Strata And QC Options

Set `analysis.ancestries` to the candidate GWAS ancestry strata you are willing
to analyze. The pipeline keeps only strata with at least
`analysis.min_stratum_n` assigned samples, writes the active list to
`results/qc/strata/active_ancestries.tsv`, and skips lower-count strata before
within-ancestry PCs or GWAS jobs are created. ADMIXTURE reference labels stay in
`admixture.labels` and can still include all report-only QC groups.

If POP-MaD excludes the expected samples but the run fails on
`popmad.max_unassigned_fraction`, raise that value deliberately in
`config/config.yaml`. Treat this as a documented QC decision, not a routine
workaround.

For imputed dosage data with MACH_R2/INFO annotations, set:

```yaml
qc:
  use_mach_r2_filter: true
  info_min: 0.8
```

## 6. Configure The SLURM Profile

### 6.1 Find Your Cluster Values

On the cluster login node, list available partitions and the account/QOS
combinations assigned to your user:

```bash
sinfo -o "%P %a %l %D %C"
sacctmgr -nP show assoc user=$USER format=Account,Partition,QOS,DefaultQOS
```

Use a partition that is available, has a sufficient time limit, and appears in
your user association. If `sinfo` prints a partition as `cpu*`, use `cpu` in
the profile; the `*` marks the cluster default and is not part of the partition
name. Set `slurm-qos` only to a QOS allowed for the selected
account and partition. If `sacctmgr` is not available on your cluster, use the
site HPC documentation or ask the cluster support team which account,
partition, and QOS should be used for batch jobs.

Choose `conda-prefix` rather than discovering it from SLURM. It should be a
writable directory for Snakemake-created rule environments, not the path to the
conda installation itself. A directory in your home folder is acceptable if it
is accessible to compute nodes and has enough quota.

### 6.2 Edit The Profile

Edit:

```text
profiles/slurm/config.yaml
```

### 6.3 Set Site Values

Set site-specific values such as:

```yaml
slurm-qos: "normal"
conda-prefix: "/path/to/shared/conda/envs"

default-resources:
  slurm_account: "your_account"
  slurm_partition: "standard"
  mem_mb: 4000
  runtime: 30
```

`default-resources` should stay active. It supplies the account, partition,
memory, and runtime defaults for each submitted job. `slurm-qos` is separate
because QOS is a SLURM executor option; leave it commented out or delete it if
your cluster does not use QOS.

For example, if your site values are account `hartj5`, QOS `normal`, partition
`cpu`, and environment directory `/lustre/home/hartj5/environments`, use:

```yaml
slurm-qos: "normal"
conda-prefix: "/lustre/home/hartj5/environments"

default-resources:
  slurm_account: "hartj5"
  slurm_partition: "cpu"
  mem_mb: 4000
  runtime: 30
```

### 6.4 Confirm Runtime And Environment Storage

Adjust memory, runtime, partition, and job limits if your cluster requires
different settings. `runtime` values in the profile are minutes. Put
`conda-prefix` on storage that compute nodes can access, not node-local scratch.

If `conda-prefix` is left unset, Snakemake may create environments under the
repository `.snakemake/` directory. That can work for small runs, but production
clusters often require an explicit shared environment directory to avoid home
quota and compute-node visibility problems.

## 7. Run Preflight Checks

### 7.1 Activate The Driver Environment

From the repository root:

```bash
<activate-command> gwas-stage1-driver
```

### 7.2 Run Production Preflight

Run the preflight script through the utility R environment from Step 2.4:

```bash
mamba run -n gwas-stage1 Rscript scripts/production_preflight.R \
  --config config/config.yaml \
  --profile profiles/slurm/config.yaml
```

If this command cannot find `gwas-stage1`, create the utility environment from
Step 2.4 after loading the conda/mamba module from Step 2.1. If your cluster
provides `conda` but not `mamba`, use `conda run -n gwas-stage1` instead.

If preflight reports `R package 'yaml' is required`, the `gwas-stage1`
environment is incomplete or older than `envs/gwas.yaml`. Update it with:

```bash
mamba env update -n gwas-stage1 -f envs/gwas.yaml --prune
```

Then confirm the required R packages are available:

```bash
mamba run -n gwas-stage1 Rscript -e 'library(yaml); library(jsonlite); cat("R utility environment OK\n")'
```

### 7.3 Confirm Preflight Passes

Preflight must pass before submission. It checks the reference package
fingerprint, required ancestry settings, optional input-manifest path existence,
SLURM profile placeholders, and whether `results/` is clean.
If it fails, fix the reported config or profile issue before continuing; do not
start the workflow and expect Snakemake to correct these settings.

## 8. Dry-Run The Workflow

### 8.1 Run The Dry-Run

```bash
snakemake -n --profile profiles/slurm
```

### 8.2 Review Planned Jobs

Review the planned jobs. Do not start the real run until the dry-run completes
without errors.

A successful dry-run confirms that Snakemake can build the workflow graph. It
does not prove that compute nodes can read every path or load every executable;
the validation and QC pilots below check those operational details.

## 9. Create Conda Environments

### 9.1 Pre-Create Environments If Needed

If compute nodes cannot access conda channels, create environments before the
main run:

```bash
snakemake --profile profiles/slurm --conda-create-envs-only
```

Run this from a login/build node that has conda channel access and can write to
the configured `conda-prefix`. After this step, compute jobs should be able to
use the prebuilt environments without reaching external channels.

### 9.2 Wait For First-Time Solves

This step may take a while the first time.

## 10. Run A Validation Pilot

### 10.1 Run Input Validation

Run the input validation target before launching the full workflow:

```bash
snakemake --profile profiles/slurm results/qc/input_validation/validation.ok
```

### 10.2 Review Validation Outputs

Review:

```text
results/logs/validation/
results/qc/genome_build/
results/config/effective_config.yaml
results/config/resolved_config.yaml
```

Do not continue until validation passes.
For failures, start with the matching files in `results/logs/validation/`; they
usually identify the missing column, inaccessible path, or build-resolution
problem directly.

### 10.3 Use Snakemake Validation For Production

Do not run `scripts/validate_config.R` directly on `config/config.yaml` for
production. Production validation uses the build-resolved config written by
Snakemake.

## 11. Run An Ancestry/QC Pilot

### 11.1 Run The QC Targets

Before the full GWAS, run the reference projection, POP-MaD, ADMIXTURE, strata,
sex-check, and relatedness QC outputs:

```bash
snakemake --profile profiles/slurm \
  results/qc/ancestry/reference/reference_prep_report.md \
  results/qc/admixture/admixture_report.md \
  results/qc/strata/strata_counts.tsv \
  results/qc/ancestry/within_ancestry_pcs.tsv \
  results/qc/sex/sex_check_summary.tsv \
  results/qc/relatedness/relatedness_summary.tsv
```

### 11.2 Review Before GWAS

Review these files before launching all GWAS jobs.
This pilot is the best place to catch reference-package, ADMIXTURE, PLINK2,
sample-strata, and sex-check problems before submitting many association jobs.

## 12. Run The Full Pipeline

### 12.1 Submit The Full Workflow

```bash
snakemake --profile profiles/slurm
```

### 12.2 Monitor And Resume

Monitor the scheduler queue and Snakemake logs. If a job fails, fix the cause
and rerun the same command. Snakemake will continue from completed outputs.

## 13. Review And Archive Results

### 13.1 Review Run-Level Files

Start with these files:

```text
results/manifests/run_manifest.tsv
results/config/resolved_config.yaml
results/qc/strata/strata_counts.tsv
results/qc/sex/sex_check_summary.tsv
results/qc/relatedness/relatedness_summary.tsv
results/qc/ancestry/reference/reference_prep_report.md
results/qc/ancestry/production/popmad_population_counts.tsv
results/qc/admixture/admixture_report.md
results/reports/
```

### 13.2 Review Main GWAS Outputs

Main GWAS outputs are:

```text
results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.plink2.glm.tsv
results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.gwas_filter_summary.tsv
results/plots/{trait}/{ancestry}/{trait}.{ancestry}.{build}.qq.png
results/plots/{trait}/{ancestry}/{trait}.{ancestry}.{build}.manhattan.png
results/plots/{trait}/{ancestry}/{trait}.{ancestry}.{build}.manhattan.pdf
results/reports/{trait}/{trait}.{ancestry}.{build}.report.md
```

The per-trait report embeds the QQ and Manhattan PNGs and summarizes sample
filtering, variant filtering, lambda GC, and the top association signals.

### 13.3 Archive The Run

Archive the final `config/config.yaml`, `results/config/resolved_config.yaml`,
`results/manifests/run_manifest.tsv`, QC reports, GWAS summary statistics,
GWAS filter summaries, and plots according to cohort policy.

## Helpful Tips

### Environment And Scheduler

- If `snakemake` says `invalid choice: 'slurm'`, activate the driver
  environment and confirm `snakemake_executor_plugin_slurm` is installed.
- If PLINK2, PLINK1, or ADMIXTURE is not found, set `tools.plink2`,
  `tools.plink1`, or `tools.admixture` to a Linux executable path or cluster
  module shim.
- If conda fails on compute nodes, create environments on a login/build node
  with `snakemake --profile profiles/slurm --conda-create-envs-only`.

### Inputs And Reference

- If the reference fingerprint fails, check `reference_package.root` and
  `reference_package.fingerprint`. Do not rebuild the package inside this
  pipeline.
- If a reference-package error lists files beginning with `._`, the package has
  macOS sidecar metadata. Current validation ignores those sidecars; real
  unmanifested analysis files still fail.
- If the optional input manifest fails, either leave
  `resources.input_manifest: ""` or use only `file_role`,
  `cohort_data_release`, and `notes` columns with rows for `sample_manifest`,
  `trait_registry`, and `study_genotype`.
- If genome-build inference fails, check that the genotype prefix is correct
  and that BIM/PVAR marker positions match the intended genome build.
- If the source BIM/PVAR contains rows with the same allele recorded twice
  (`A/A`, `C/C`, and similar), the workflow excludes those malformed variants
  before PLINK2 conversion and writes an `*.invalid_*_alleles.tsv` report next
  to the affected QC prefix.
- If validation warns that fine-scale reference populations are below
  `popmad.min_reference_population_n`, those populations will be skipped during
  POP-MaD model fitting. Fix the reference package or adjust the threshold only
  if a configured ancestry has no retained population model.
- If sex check fails because no X/Y markers exist, either provide genotype data
  with sex chromosomes or set `sex_check.allow_no_sex_markers: true` only with
  documented external sex QC.
- If production fails for an empty trait/ancestry cell, remove that ancestry or
  trait from the config/registry before rerunning.

### Reruns And Logs

- If `results/` is not clean before a new production run, move or archive the
  previous results directory first.
- For any failed rule, inspect the matching file under `results/logs/` before
  rerunning.
- Snakemake engine logs are under `.snakemake/log/`; PLINK2 GWAS logs are under
  `results/gwas/{trait}/{ancestry}/plink2_raw/`.
- If a run is interrupted and Snakemake reports a lock, confirm no Snakemake
  process is still active, then unlock the workflow:

```bash
snakemake --unlock --profile profiles/slurm
```
