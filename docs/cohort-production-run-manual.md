# Cohort Production Troubleshooting Guide

Use this guide when the
[Cohort Production Quickstart](cohort-production-quickstart.md) fails or when a
cohort run needs a targeted check before rerunning.

For any failure, start with the rule log named by Snakemake. Most logs are under
`results/logs/`. SLURM wrapper logs are under `.snakemake/slurm_logs/`, and
Snakemake engine logs are under `.snakemake/log/`.

After fixing the cause, rerun the same Snakemake command. Completed outputs are
reused automatically.

## Step 1: Clone The Repository

### SLURM Jobs Cannot See Files

**Symptom:** Dry-run works, but submitted jobs fail with missing input files.

**Likely cause:** The repository, input genotypes, reference package,
`conda-prefix`, or `results/` directory is only visible from the login node.

**Fix:** Use full paths on storage visible from login and compute nodes. Avoid
`~` and shell variables in YAML paths.

### Protected Data Appears In The Repository

**Symptom:** Cohort files, private paths, or credentials appear in `git status`
or are copied into repository folders.

**Likely cause:** Production inputs were placed under the cloned repository.

**Fix:** Keep protected cohort data outside the repository unless your
data-use rules explicitly allow it. Point `config/config.yaml` to the protected
data location.

## Step 2: Download The Reference Package

### Reference Fingerprint Fails

**Symptom:** Preflight or validation reports a reference fingerprint mismatch.

**Likely cause:** `reference_package.root` points to the `.tar.gz`, the
configured fingerprint is the archive checksum, the wrong package was unpacked,
or package contents changed after release.

**Fix:** Use the unpacked package directory as `reference_package.root`. Use the
value in:

```bash
cat /path/to/stage1_reference_package/content_fingerprint.sha256
```

Do not edit `reference_package.observed_fingerprint`; the workflow writes it in
`results/config/resolved_config.yaml`.

### Reference Package Content Fails Validation

**Symptom:** Validation reports missing manifests, extra package files,
absolute paths, `..` paths, size mismatches, SHA-256 mismatches, raw Hail/VCF/BCF
files, or missing panel files.

**Likely cause:** The package is incomplete, was edited, or was repackaged
incorrectly.

**Fix:** Re-extract the approved package and use that clean directory. macOS
sidecar files such as `._*`, `.DS_Store`, and `__MACOSX/` are ignored; real
unmanifested analysis files still fail.

### Reference Panel Build Does Not Match

**Symptom:** Validation reports that the package is missing a POP-MaD or
ADMIXTURE panel for the inferred study build.

**Likely cause:** The study build and reference package build do not align.

**Fix:** Inspect `results/qc/genome_build/genome_build_marker_matches.tsv`.
Use a reference package containing the inferred build, or fix genotype marker
IDs/coordinates if build inference is wrong.

## Step 3: Prepare Cohort Inputs

### Genotype Files Are Missing

**Symptom:** Validation cannot find `.pgen/.pvar/.psam` or `.bed/.bim/.fam`, or
the expected path looks like `.pgen.pgen`.

**Likely cause:** `genotypes.prefix` includes a file extension.

**Fix:** Set the prefix without an extension:

```yaml
genotypes:
  type: "pgen"
  prefix: "/path/to/cohort/genotypes_without_extension"
```

VCF/BCF is not accepted directly. Convert upstream to PLINK first.

### Sample Manifest IDs Do Not Match Genotypes

**Symptom:** Validation reports that sample manifest IDs are absent from the
genotype files.

**Likely cause:** `FID/IID` formatting differs between the sample manifest TSV
(phenotype + covariate file) and the FAM/PSAM files.

**Fix:** Compare IDs exactly. Preserve leading zeros, spelling, and whitespace.
Order does not matter.

### Manifest Columns Or Values Fail Validation

**Symptom:** Validation reports missing phenotype/covariate columns, nonnumeric
covariates, duplicate `FID/IID`, or invalid sex codes.

**Likely cause:** The manifest was exported as CSV/Excel, headers changed, or
categorical values were not recoded.

**Fix:** Export tab-delimited TSV with exact headers. Include `FID`, `IID`,
`age`, `age2`, `sex`, every phenotype column, and every non-PC covariate.
Sex must be `1`, `2`, `0`, `NA`, `-9`, or `.`.

### `age2` Looks Wrong

**Symptom:** `age2` values are very large or hard to interpret.

**Likely cause:** Raw `age^2` was used.

**Fix:** Use a centered quadratic term and record the mean age used:

```r
mean_age <- mean(manifest$age, na.rm = TRUE)
manifest$age2 <- (manifest$age - mean_age)^2
```

### Trait Registry Fails Before Jobs Start

**Symptom:** Snakemake fails while building the DAG with a missing or empty
`trait_id` message.

**Likely cause:** `inputs.trait_registry` is wrong, unreadable, or missing the
`trait_id` column.

**Fix:** Confirm the trait registry path and required columns:

```text
trait_id	phenotype_column	case_value	control_value	missing_values
```

Use simple, unique `trait_id` values because they are used in output paths.

### Unexpected Phenotype Value

**Symptom:** Trait-file generation reports an unexpected phenotype value.

**Likely cause:** The sample manifest contains a value not listed as case,
control, or missing for that trait.

**Fix:** Update the trait registry `case_value`, `control_value`, and
`missing_values`, or recode the phenotype column in the manifest.

## Step 4: Configure The Run

### Wrong File Is Being Edited

**Symptom:** Changes do not affect the run, or validation still sees old paths.

**Likely cause:** The configured production files were not edited.

**Fix:** For a standard run, configure:

```text
config/config.yaml
profiles/slurm/config.yaml
resources/manifests/software.tsv
resources/manifests/reference_data.tsv
```

Start `config/config.yaml` from `config/config.template.yaml`.

### Old Config Fields Are Still Present

**Symptom:** Validation reports unsupported fields such as `ancestry_file`,
`pcs_file`, `projected_pcs_file`, or `reference_pcs_file`.

**Likely cause:** The config still contains precomputed ancestry/PC settings.

**Fix:** Remove those fields. Production ancestry and PC covariates are computed
from the reference package.

### Production Safety Gates Fail

**Symptom:** Preflight or validation requires ancestry, ADMIXTURE, sex-check, or
PC settings.

**Likely cause:** Production safety settings were changed.

**Fix:** Keep these unless the analysis plan explicitly says otherwise:

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

### Tool Executables Are Not Found

**Symptom:** Validation or a cluster job cannot find `plink2`, `plink`, or
`admixture`.

**Likely cause:** The Snakemake rule environment was not created or is not
visible on compute nodes, or `tools.*` points to an invalid site-specific
binary.

**Fix:** Leave `tools.*` as `plink2`, `plink`, and `admixture` when using the
default conda rule environment. If using site-managed software instead, set full
Linux paths or approved module shims and record versions in
`resources/manifests/software.tsv`.

### Optional Input Manifest Fails

**Symptom:** Validation reports an error for `resources.input_manifest`.

**Likely cause:** The optional input manifest path is set but the file is
missing or has the wrong schema.

**Fix:** Most runs should leave `resources.input_manifest: ""`. If your team
uses the optional note file, keep it small and include only these columns:
`file_role`, `cohort_data_release`, and `notes`.

### Strata Are Skipped

**Symptom:** A configured ancestry does not launch GWAS jobs.

**Likely cause:** Fewer than `analysis.min_stratum_n` samples were assigned to
that ancestry, or the ancestry label is not present in POP-MaD output.

**Fix:** Review:

```text
results/qc/strata/active_ancestries.tsv
results/qc/strata/excluded_ancestries.tsv
results/qc/strata/strata_counts.tsv
```

Adjust `analysis.ancestries` or `analysis.min_stratum_n` only as a documented
analysis decision.

### Too Many Samples Are Unassigned

**Symptom:** Strata creation fails on `popmad.max_unassigned_fraction`.

**Likely cause:** POP-MaD excluded many samples, ancestry labels are not
configured, or low-count ancestries were removed from active strata.

**Fix:** Review:

```text
results/qc/ancestry/production/popmad_excluded.tsv
results/qc/strata/unassigned_ancestry.tsv
```

Raise `popmad.max_unassigned_fraction` only when the dropped samples are
expected and documented.

### Imputed Dosage Filters Are Not Applied

**Symptom:** Imputed variants are present but MACH_R2/INFO filtering is not used.

**Likely cause:** INFO filtering is disabled in the config.

**Fix:** For imputed dosage data with MACH_R2/INFO annotations, set:

```yaml
qc:
  use_mach_r2_filter: true
  info_min: 0.8
```

## Step 5: Set Up The Environment

### Conda Module Is Unclear

**Symptom:** `conda` or `mamba` is missing, or environment creation uses the
wrong Python-only module.

**Likely cause:** The conda-providing cluster module is not loaded.

**Fix:** Load the site conda/mamba module first. On the project cluster:

```bash
module load miniforge3/23.3.1
```

### SLURM Executor Is Missing

**Symptom:** Snakemake reports `invalid choice: 'slurm'`, or this import fails:

```bash
python -c "import snakemake_executor_plugin_slurm"
```

**Likely cause:** The driver environment is not active or is out of date.

**Fix:** Activate or update the Snakemake driver environment:

```bash
conda activate gwas-stage1-driver
mamba env update -n gwas-stage1-driver -f envs/snakemake-driver.yaml --prune
```

### Snakemake Rejects The Conda Version

**Symptom:** Snakemake reports `Conda must be version 24.7.1 or later`.

**Likely cause:** An old base/module conda is shadowing the driver environment.

**Fix:** Check the conda inside the active environment:

```bash
echo "$CONDA_PREFIX"
"$CONDA_PREFIX/bin/conda" --version
python -c "import conda; print(conda.__version__)"
```

If needed:

```bash
mamba install -n gwas-stage1-driver -c conda-forge "conda>=24.7.1"
```

### Utility R Environment Is Missing Packages

**Symptom:** Preflight fails with `R package 'yaml' is required` or
`gwas-stage1` cannot be found.

**Likely cause:** The workflow utility environment was not created or is stale.

**Fix:** Create or update it:

```bash
mamba env create -f envs/gwas.yaml
mamba env update -n gwas-stage1 -f envs/gwas.yaml --prune
```

Then check:

```bash
mamba run -n gwas-stage1 Rscript -e 'library(yaml); library(jsonlite); cat("R utility environment OK\n")'
```

### Rule Environments Fail On Compute Nodes

**Symptom:** Conda solve or download errors happen only after jobs start.

**Likely cause:** Compute nodes cannot access conda channels, or
`conda-prefix` is unset or inaccessible.

**Fix:** Set `conda-prefix` in `profiles/slurm/config.yaml` to a writable
directory accessible to compute nodes. Pre-create rule environments from a
login/build node:

```bash
snakemake --profile profiles/slurm --conda-create-envs-only
```

Strict channel priority warnings are not the same as failed runs, but strict
priority improves reproducibility:

```bash
conda config --set channel_priority strict
```

### SLURM Account, Partition, Or QOS Is Wrong

**Symptom:** Jobs stay pending or fail submission with invalid account,
partition, or QOS.

**Likely cause:** `profiles/slurm/config.yaml` still has placeholders or uses a
QOS/partition not assigned to your account.

**Fix:** Check available values:

```bash
sinfo -o "%P %a %l %D %C"
sacctmgr -nP show assoc user=$USER format=Account,Partition,QOS,DefaultQOS
```

If `sinfo` shows `cpu*`, use `cpu`; the `*` marks the default. Put account and
partition under `default-resources`. Set `slurm-qos` separately only if the
cluster uses QOS.

`conda-prefix` is not the conda install path. It is a writable directory for
Snakemake-created rule environments, and it may be in your home directory if
compute nodes can access it and quota is sufficient.

### Jobs Hit Time Or Memory Limits

**Symptom:** SLURM marks jobs as timed out, out of memory, or killed.

**Likely cause:** Profile resources are too small for the cohort or selected
partition.

**Fix:** Increase rule resources in `profiles/slurm/config.yaml`. `runtime`
values are minutes. Keep `default-resources` active; it supplies account,
partition, memory, and runtime defaults to submitted jobs.

## Step 6: Run Basic Checks

### Preflight Fails

**Symptom:** `production_preflight.R` reports failures.

**Likely cause:** A production gate failed before Snakemake submission.

**Fix:** Read the reported failures directly. Preflight checks the reference
fingerprint, required ancestry settings, optional input-manifest path, SLURM
profile placeholders, and whether `results/` is clean.

```bash
mamba run -n gwas-stage1 Rscript scripts/production_preflight.R \
  --config config/config.yaml \
  --profile profiles/slurm/config.yaml
```

If rerunning into an existing `results/` directory is intentional, use the
script's `--allow-existing-results` option.

### Dry-Run Passes But Real Jobs Fail

**Symptom:** `snakemake -n --profile profiles/slurm` succeeds, but submitted
jobs fail.

**Likely cause:** Dry-run only proves Snakemake can build the DAG. It does not
prove compute-node path visibility, executables, or conda access.

**Fix:** Run the validation target before the full workflow:

```bash
snakemake --profile profiles/slurm results/qc/input_validation/validation.ok
```

### Validation Fails

**Symptom:** The validation target fails.

**Likely cause:** The resolved config, genotype files, manifest, trait registry,
tools, build inference, or reference package failed a production check.

**Fix:** Start with:

```text
results/logs/validation/
results/qc/genome_build/
results/config/effective_config.yaml
results/config/resolved_config.yaml
```

Do not run `scripts/validate_config.R` directly on `config/config.yaml` for
production. Snakemake validation uses the build-resolved config.

### Genome Build Cannot Be Inferred

**Symptom:** Build inference reports no marker matches, weak support, or weak
margin.

**Likely cause:** The marker panel, variant IDs, or marker coordinates do not
match the study genotype build.

**Fix:** Inspect:

```text
results/qc/genome_build/genome_build_marker_matches.tsv
```

Then verify rsIDs and coordinates in the study BIM/PVAR.

### Smaller QC Pilot Is Needed

**Symptom:** You want to catch ancestry, ADMIXTURE, sex-check, or relatedness
problems before launching all GWAS jobs.

**Fix:** Run the QC targets:

```bash
snakemake --profile profiles/slurm \
  results/qc/ancestry/reference/reference_prep_report.md \
  results/qc/admixture/admixture_report.md \
  results/qc/strata/strata_counts.tsv \
  results/qc/ancestry/within_ancestry_pcs.tsv \
  results/qc/sex/sex_check_summary.tsv \
  results/qc/relatedness/relatedness_summary.tsv
```

## Step 7: Run The Pipeline

### Finding The Failed Rule

**Symptom:** A run failed after your SSH session closed.

**Fix:** Check the newest Snakemake engine logs and the rule logs:

```bash
ls -lt .snakemake/log/ | head
find results/logs -type f -name "*.log" -mtime -2
```

The Snakemake error usually names the failed rule and its log path.

### Snakemake Keeps Running After A Job Fails

**Symptom:** A rule fails, but other jobs continue.

**Likely cause:** The SLURM profile uses `keep-going: true`.

**Fix:** This is expected. Independent jobs may finish. Fix the failed rule and
rerun the same command.

### Workflow Is Locked

**Symptom:** Snakemake reports a lock after an interrupted run.

**Likely cause:** A previous Snakemake process was interrupted.

**Fix:** First confirm no Snakemake process is still active. Then unlock:

```bash
snakemake --unlock --profile profiles/slurm
```

### Malformed Duplicate Allele Codes

**Symptom:** PLINK2 reports a duplicate allele code such as `A/A` or `C/C`.

**Likely cause:** The source BIM/PVAR contains malformed variant rows.

**Fix:** The workflow excludes those rows before PLINK2 conversion and writes a
report next to the affected QC prefix, for example:

```text
*.invalid_bim_alleles.tsv
*.invalid_pvar_alleles.tsv
```

Review the report. If many rows are affected, investigate upstream genotype
preparation.

### Too Few Shared Ancestry Variants

**Symptom:** Ancestry reference preparation reports too few shared variants.

**Likely cause:** Build mismatch, rsID mismatch, allele mismatch, palindromic
SNP exclusion, or strict filters.

**Fix:** Review:

```text
results/qc/ancestry/reference/shared_variant_mismatches.tsv
results/qc/ancestry/reference/reference_prep_report.md
```

Fix build/ID problems before lowering thresholds.

### POP-MaD Model Or Projection Fails

**Symptom:** POP-MaD reports no valid model for an ancestry, or reference
projection validation has low PC correlations.

**Likely cause:** Reference metadata labels are wrong, too few fine-population
samples remain after thresholds, or allele harmonization/scoring does not
reproduce fitted reference PCs.

**Fix:** Check metadata `population` and `super_population` labels,
`popmad.min_reference_population_n`, PLINK2 version, and reference panel
integrity. Review:

```text
results/qc/ancestry/reference/reference_projection_validation.tsv
results/qc/ancestry/reference/reference_prep_report.md
```

### ADMIXTURE QC Fails

**Symptom:** ADMIXTURE validation fails on `k`, labels, harmonization, `.pop`
writing, or Q-label mapping.

**Likely cause:** `admixture.labels` does not match `admixture.k`, labels are
absent from reference metadata, study/reference builds differ, sample IDs are
ambiguous, or supervised components cannot be mapped cleanly.

**Fix:** Align labels exactly with reference super-populations and review:

```text
results/logs/admixture/
results/qc/admixture/raw/shared_variant_mismatches.tsv
results/qc/admixture/raw/ld_prune/admixture_ld_prune.prune.in
results/qc/admixture/admixture_report.md
```

ADMIXTURE is report-only QC; it does not replace POP-MaD strata.

### Sex Check Removes Samples

**Symptom:** Sex-check fails or removes every sample.

**Likely cause:** Manifest sex coding conflicts with genotype-inferred sex, or
thresholds are inappropriate for the dataset.

**Fix:** Review:

```text
results/qc/sex/sexcheck.tsv
results/qc/sex/sex_check_summary.tsv
results/logs/sample_prep/sex_check.log
```

Fix manifest sex coding first. Custom chrX thresholds must include both
`max_female_xf` and `min_male_xf`. Y-rate thresholds are optional but must be
supplied as a pair.

If genotype data are autosome-only, set `sex_check.allow_no_sex_markers: true`
only after external sex QC has been completed and documented.

### Relatedness Marker Set Is Empty

**Symptom:** Relatedness pruning does not produce a non-empty `prune.in`, or
the unrelated keep file is empty.

**Likely cause:** MAF, missingness, autosome, SNP, sex-check, or KING filters
removed all usable markers or samples.

**Fix:** Review:

```text
results/logs/sample_prep/prepare_relatedness_markers.log
results/logs/sample_prep/make_unrelated_keep.log
results/qc/relatedness/relatedness_summary.tsv
```

Relax relatedness-specific filters only if QC supports the change.

### Keep File Is Empty For One Ancestry

**Symptom:** Downstream GWAS has no samples for an ancestry.

**Likely cause:** The intersection of ancestry, sex-check, and unrelated keep
files removed the stratum.

**Fix:** Compare:

```text
results/qc/strata/{ancestry}.keep.tsv
results/qc/strata/{ancestry}.unrelated.keep.tsv
results/qc/sex/sex_checked.keep.tsv
results/qc/relatedness/unrelated.king.cutoff.in.id
```

### Empty Case/Control Cells

**Symptom:** The workflow stops before GWAS with empty trait/ancestry cells.

**Likely cause:** An active stratum has no cases, no controls, or no
phenotype-complete samples for a trait.

**Fix:** Inspect:

```text
results/qc/strata/strata_counts.tsv
```

Remove the trait or ancestry from the production run, or revise the phenotype
definition.

### PC Covariates Are Missing

**Symptom:** Trait-file generation reports missing or non-finite `PC` covariates
for kept samples.

**Likely cause:** Within-ancestry PCA did not produce rows for samples entering
GWAS.

**Fix:** Review within-ancestry PCA logs and final keep files. Production does
not support missing PCs.

### PLINK2 GWAS Fails

**Symptom:** A `run_plink2_gwas` job fails.

**Likely cause:** Empty keep file, bad phenotype/covariate file, invalid
genotype prefix, unavailable PLINK2, or model/QC filters removing data.

**Fix:** Inspect:

```text
results/logs/gwas/
results/gwas/{trait}/{ancestry}/plink2_raw/
```

If harmonization reports no `.glm` output or no requested `TEST` rows, check
the raw `.glm.*` files and confirm `gwas.test` matches a PLINK2 test value such
as `ADD`.

## Step 8: Review Outputs

### Final Report Is Missing

**Symptom:** No report appears under `results/reports/`.

**Likely cause:** One of the required upstream report inputs is missing.

**Fix:** Follow the missing input named by Snakemake back to its producing rule
and log under `results/logs/`.

### GWAS Plots Fail Or Look Empty

**Symptom:** QQ or Manhattan plotting fails, or plots contain very few points.

**Likely cause:** Harmonized GWAS stats are empty or missing required columns.

**Fix:** Inspect the harmonized stats and filter summary:

```text
results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.plink2.glm.tsv
results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.gwas_filter_summary.tsv
```

Do not diagnose this as a plotting issue until the summary statistics file is
confirmed non-empty.

### Report Shows High Lambda Or Heavy Variant Loss

**Symptom:** The final report shows high lambda GC, few valid P values, many
removed variants, or unexpected sample loss.

**Likely cause:** Model, covariate, phenotype, or QC-filter issues upstream of
reporting.

**Fix:** Start with the report, then review the filter summary and PLINK logs:

```text
results/reports/
results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.gwas_filter_summary.tsv
results/gwas/{trait}/{ancestry}/plink2_raw/
```

### POP-MaD Projection Plot Fails

**Symptom:** The POP-MaD 3D projection plot is missing or fails during
reporting.

**Likely cause:** PC1/PC2/PC3 columns are absent or projected PC rows are
non-finite.

**Fix:** Check:

```text
results/qc/ancestry/reference/reference_pcs.tsv
results/qc/ancestry/reference/study_projected_pcs.tsv
results/qc/ancestry/production/popmad_assignments.tsv
results/qc/ancestry/production/popmad_excluded.tsv
```

## Step 9: Compress And Export Results

### Archive Contains Old Or Extra Outputs

**Symptom:** The export archive contains outputs from earlier runs.

**Likely cause:** A previous `results/` directory was reused.

**Fix:** For a new production run, archive or move old `results/` first. Before
sharing, verify that the final report, run manifest, resolved config, QC tables,
plots, and GWAS files all belong to the intended cohort release.

### Archive Contains Private Runtime Details

**Symptom:** Governance review flags private paths or cohort metadata in the
archive.

**Likely cause:** Config snapshots, logs, reports, and manifests record runtime
paths and cohort release labels.

**Fix:** Review archive contents before export. Share only to approved
destinations under cohort policy.

### What To Review Before Export

Start with the final reports:

```text
results/reports/
```

Use these supporting files for audit:

```text
results/manifests/run_manifest.tsv
results/config/resolved_config.yaml
results/qc/strata/strata_counts.tsv
results/qc/sex/sex_check_summary.tsv
results/qc/relatedness/relatedness_summary.tsv
results/qc/ancestry/reference/reference_prep_report.md
results/qc/ancestry/production/popmad_population_counts.tsv
results/qc/admixture/admixture_report.md
results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.plink2.glm.tsv
results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.gwas_filter_summary.tsv
results/plots/{trait}/{ancestry}/
```
