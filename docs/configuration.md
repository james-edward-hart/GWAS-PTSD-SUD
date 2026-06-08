# Configuration Guide

This guide explains every field in `config/config.template.yaml`. Copy the
template to `config/config.yaml` for a real run, then edit the copy.

```bash
cp config/config.template.yaml config/config.yaml
```

Keep config values explicit. A failed validation step is preferable to running a
GWAS with the wrong samples, phenotype, ancestry labels, or build label.

## How Config Is Used

Snakemake reads `config/config.yaml` when building the DAG. The workflow then
writes two config snapshots:

```text
results/config/effective_config.yaml
results/config/resolved_config.yaml
```

`effective_config.yaml` records the merged config that Snakemake saw at parse
time. `resolved_config.yaml` is written after genome-build inference and
reference-package resolution. Downstream production validation and analysis
rules use `resolved_config.yaml`.

Do not edit `results/config/resolved_config.yaml` by hand. To change a run,
edit `config/config.yaml` and rerun Snakemake.

## Files To Configure

For a production run, configure these repository files:

| File | Where it lives | Purpose |
| --- | --- | --- |
| `config/config.yaml` | Copy from `config/config.template.yaml` in the repository root | Main run config. This is the file Snakemake reads. |
| `profiles/slurm/config.yaml` | Bundled scheduler profile in `profiles/slurm/` | Cluster account, partition, QOS, job resources, and shared conda location. |
| `resources/manifests/software.tsv` | Repository provenance file in `resources/manifests/` | Software names, versions, and paths used for the run. |
| `resources/manifests/reference_data.tsv` | Repository provenance file in `resources/manifests/` | Reference package provenance. Runtime validation uses the fingerprint in `config/config.yaml`. |

The cohort input files usually live outside the repository on secure,
cluster-visible storage. Do not copy protected cohort data into this repository
unless your data-governance plan explicitly allows it. Point to those files from
`config/config.yaml`:

| Cohort input | Config field |
| --- | --- |
| Sample manifest TSV (phenotype + covariate file) | `inputs.sample_manifest` |
| Trait registry TSV | `inputs.trait_registry` |
| Study genotype prefix without extension | `genotypes.prefix` |
| Unpacked ancestry reference package directory | `reference_package.root` |

`resources/manifests/input_data.tsv` is optional. Most production configs should
leave `resources.input_manifest: ""`. Use it only for short free-text notes about
the cohort input release, not as the source of truth for input paths.

## Production Minimum

Production runs should set these fields first:

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

Production mode also requires a clean, unpacked reference package whose
`content_fingerprint.sha256` matches `reference_package.fingerprint`.

## Input Table Schemas

### Sample Manifest

`inputs.sample_manifest` must be a TSV (phenotype + covariate file) with one
row per analysis sample.

Required columns:

```text
FID	IID	age	age2	sex
```

Also include every phenotype column named in the trait registry and every non-PC
covariate used by `gwas.default_covariates`, `gwas.extra_covariates`, or the
trait registry `covariates` column.

Rules:

- `FID`/`IID` pairs must be unique.
- Every sample manifest `FID`/`IID` pair must exist in the genotype files.
- `sex` must use PLINK-style codes: `1`, `2`, `0`, `NA`, `-9`, or `.`.
- Non-PC covariates must be numeric, except configured missing values.

### Trait Registry

`inputs.trait_registry` must be a TSV with these columns:

```text
trait_id	phenotype_column	case_value	control_value	missing_values
```

Optional column:

```text
covariates
```

Rules:

- `trait_id` values define the `{trait}` wildcard in GWAS outputs.
- `phenotype_column` must exist in the sample manifest.
- `case_value`, `control_value`, and `missing_values` define how the sample
  manifest phenotype is recoded for PLINK2.
- PLINK phenotype coding is written as control `1`, case `2`, missing `NA`.
- Trait-specific `covariates` are comma-separated and appended to the global
  GWAS covariates for that trait.

### Ancestry And PCs

Do not provide user-supplied ancestry labels, projected study PCs, reference
PCs, or GWAS PC covariates. The production workflow computes them from the
fingerprinted reference package:

```text
results/qc/ancestry/reference/reference_pcs.tsv
results/qc/ancestry/reference/study_projected_pcs.tsv
results/qc/ancestry/production/popmad_assignments.tsv
results/qc/ancestry/within_ancestry_pcs.tsv
results/plots/ancestry/popmad_reference_study_pcs.png
```

`results/qc/ancestry/within_ancestry_pcs.tsv` is the PC table used for GWAS
covariates. The POP-MaD projection plot is embedded in the final GWAS reports.

## Parameter Reference

### `project`

| Parameter | Required | Description |
| --- | --- | --- |
| `analysis_name` | Yes | Human-readable run name written to reports and manifests. |
| `cohort_data_release` | Recommended | Cohort freeze or release label written to provenance outputs. |
| `genome_build` | Yes | Must be `auto`. The workflow infers the build from genotype markers. |

`project.inferred_genome_build` is added to `resolved_config.yaml` after build
inference. Do not set it in the editable config.

### `reference_package`

| Parameter | Required | Description |
| --- | --- | --- |
| `root` | Yes | Path to the unpacked prebuilt reference package directory. |
| `fingerprint` | Production | Expected content fingerprint from `content_fingerprint.sha256` inside the unpacked package. |

The package must contain `content_fingerprint.sha256`, `file_manifest.tsv`, and
`panel_manifest.tsv`. The resolver validates package file hashes and injects
package-derived fields into `resolved_config.yaml`, including
`reference_package.observed_fingerprint` and build-matched reference paths for
`ancestry_reference` and `admixture`.

macOS sidecar files such as `._*`, `.DS_Store`, and `__MACOSX/` are ignored
during package validation because they are copy metadata, not reference content.
Other unmanifested files still fail validation.

Set only `root` and `fingerprint` for routine production runs. Do not manually
add `observed_fingerprint`.

### `analysis`

| Parameter | Required | Description |
| --- | --- | --- |
| `ancestries` | Yes | Candidate GWAS ancestry strata. Labels must exist in the resolved POP-MaD reference metadata `super_population` values. |
| `min_stratum_n` | Recommended | Minimum assigned sample count required for a candidate ancestry to become an active GWAS stratum. Template uses `50`; lower-count candidate strata are reported and excluded from downstream GWAS. |

Each listed ancestry is evaluated during strata creation. Only active strata
expand within-ancestry PCs, GWAS, plots, and reports. In production, each active
trait/ancestry cell must contain at least one case and one control.
This setting is separate from `admixture.labels`: ADMIXTURE can use all
reference ancestry labels for report-only QC even when GWAS is run for fewer
analysis strata.

### `inputs`

| Parameter | Required | Description |
| --- | --- | --- |
| `sample_manifest` | Yes | TSV (phenotype + covariate file) with `FID`, `IID`, demographics, phenotypes, and non-PC covariates. |
| `trait_registry` | Yes | TSV defining traits, phenotype coding, missing values, and optional trait covariates. |

The workflow generates ancestry assignments and PC covariates from the
fingerprinted reference package. Do not add ancestry-label or PC-path fields to
`inputs`, including `ancestry_mode`, `ancestry_file`, `pcs_file`,
`projected_pcs_file`, or `reference_pcs_file`.

### `tools`

| Parameter | Required | Description |
| --- | --- | --- |
| `plink2` | Yes | PLINK2 executable path or `PATH` command. Validation checks that it can run `--version`. |
| `plink1` | When ADMIXTURE enabled | PLINK 1.9 executable path or `PATH` command. Used for the supervised ADMIXTURE sample merge. |
| `admixture` | When ADMIXTURE enabled | ADMIXTURE executable path or `PATH` command. Required when `admixture.enabled: true`. |

On HPC, these must resolve to Linux executables, module shims, or active
environment commands. Do not reuse workstation-specific macOS binaries on the
cluster.

### `genotypes`

| Parameter | Required | Description |
| --- | --- | --- |
| `type` | Yes | `pgen` for `PGEN/PVAR/PSAM`, or `bed` for `BED/BIM/FAM`. |
| `prefix` | Yes | Genotype file prefix without extension. |

For `type: "pgen"`, the workflow requires:

```text
{prefix}.pgen
{prefix}.pvar
{prefix}.psam
```

For `type: "bed"`, the workflow requires:

```text
{prefix}.bed
{prefix}.bim
{prefix}.fam
```

VCF/BCF is not accepted directly. Convert upstream to PLINK format.

### `genome_build`

| Parameter | Required | Description |
| --- | --- | --- |
| `marker_file` | Yes | Marker panel used to infer genome build. The bundled file is `resources/build_markers.tsv`. |
| `min_markers` | Recommended | Minimum informative markers required before accepting a build call. Template uses `50`; focused development checks can use lower thresholds. |
| `min_match_fraction` | Recommended | Minimum fraction of informative markers that must match the winning build. Template uses `0.95`. |
| `min_marker_margin` | Recommended | Minimum marker-count lead of the winning build over the runner-up. Template uses `20`. |
| `min_fraction_margin` | Recommended | Minimum match-fraction lead of the winning build over the runner-up. Template uses `0.20`. |

Build inference writes:

```text
results/qc/genome_build/genome_build.txt
results/qc/genome_build/genome_build_marker_matches.tsv
```

The pipeline does not lift genotype coordinates or summary statistics.

### `popmad`

| Parameter | Required | Description |
| --- | --- | --- |
| `pcs` | Yes | Number of PCs used for POP-MaD assignment. Must be 1-20; template uses `10`. |
| `reference_outlier_sd` | Yes | Reference-population outlier cutoff in SD units before assignment. |
| `min_confidence` | Yes | Minimum assignment confidence. Lower-confidence samples are excluded as ambiguous. |
| `min_reference_population_n` | Recommended | Minimum reference samples required for a fine-scale population to contribute a POP-MaD model. Populations below this threshold are skipped; production fails only if a configured ancestry has no retained population model. |
| `max_unassigned_fraction` | Recommended | Maximum allowed fraction of samples dropped from active GWAS strata because they are unassigned, ambiguous/outlying, assigned to an ancestry not listed in `analysis.ancestries`, or assigned to a candidate stratum below `analysis.min_stratum_n`. Template uses `0.07`. |

POP-MaD ancestry outputs are written under
`results/qc/ancestry/production/`. Within-ancestry GWAS PCs are written to
`results/qc/ancestry/within_ancestry_pcs.tsv`.
Active and excluded GWAS strata are written to
`results/qc/strata/active_ancestries.tsv` and
`results/qc/strata/excluded_ancestries.tsv`.

### `admixture`

ADMIXTURE is a report-only QC branch. It does not define ancestry strata or GWAS
covariates. The run summary reports the mean study ADMIXTURE proportion for
each configured ancestry label.

| Parameter | Required | Description |
| --- | --- | --- |
| `enabled` | Production | `true` or `false`. Production requires `true`. |
| `mode` | When enabled | Must be `supervised`. |
| `k` | When enabled | Number of supervised ancestry labels. Must be an integer >= 2. |
| `labels` | When enabled | Ordered ancestry labels. Length must equal `k`, labels must be unique, and labels must exist in reference metadata. |
| `filters.maf_min` | Recommended | Minimum MAF for ADMIXTURE marker filtering. |
| `filters.geno_missing_max` | Recommended | Maximum marker missingness for ADMIXTURE marker filtering. |
| `filters.snps_only_acgt` | Recommended | Keep only A/C/G/T SNPs when `true`. |
| `filters.autosome_only` | Recommended | Keep autosomes only when `true`. |
| `filters.max_alleles` | Recommended | Maximum allele count; template uses `2`. |
| `filters.remove_duplicate_ids` | Recommended | Remove duplicate variant IDs before ADMIXTURE marker preparation. |
| `filters.exclude_palindromic` | Recommended | Remove strand-ambiguous A/T and C/G SNPs. |
| `min_pruned_variants` | Recommended | Minimum number of LD-pruned variants required for ADMIXTURE. |
| `exclusion_regions` | Optional | Long-range LD/problem-region TSV. Leave blank to skip. |
| `ld_prune.window` | Recommended | PLINK LD-pruning window. ADMIXTURE template uses variant count `50`. |
| `ld_prune.step` | Recommended | PLINK LD-pruning step. Must be positive. |
| `ld_prune.r2` | Recommended | PLINK LD-pruning `r2`. Must be between 0 and 1. |

Reference genotype, metadata, build, exclusion-region, and source fields are
resolved from the reference package. Do not add local reference paths to
editable production configs.

If `exclusion_regions` is set, the TSV must include:

```text
chrom	start	end	label
```

In production, include a `build` column or put a build label such as `GRCh37` or
`GRCh38` in the filename. The build must match the inferred study build.

### `ancestry_reference`

This branch prepares reference-projected PCs for POP-MaD. Production uses the
prebuilt reference package and treats it as read-only.

| Parameter | Required | Description |
| --- | --- | --- |
| `enabled` | Production | `true` or `false`. Production requires `true`. |
| `reference_genome_build` | Resolved in production | Reference build. Leave blank in the template; package resolution fills it. |
| `variant_set` | Optional/resolved | Empty, `pre_ld_pruned`, `workflow_pruned`, or `unpruned`. Package panels can set this. |
| `filters.maf_min` | Recommended | Minimum MAF for ancestry reference marker filtering. |
| `filters.geno_missing_max` | Recommended | Maximum marker missingness for ancestry reference marker filtering. |
| `filters.snps_only_acgt` | Recommended | Keep only A/C/G/T SNPs when `true`. |
| `filters.autosome_only` | Recommended | Keep autosomes only when `true`. |
| `filters.max_alleles` | Recommended | Maximum allele count; template uses `2`. |
| `filters.remove_duplicate_ids` | Recommended | Remove duplicate variant IDs before ancestry reference prep. |
| `filters.exclude_palindromic` | Recommended | Remove strand-ambiguous A/T and C/G SNPs. |
| `warn_shared_variants_below` | Recommended | Warn when the shared study/reference marker count is below this value. |
| `min_shared_variants` | Recommended | Hard fail when the shared study/reference marker count is below this value. |
| `exclusion_regions` | Optional/resolved | Long-range LD/problem-region TSV. Usually resolved from the package in production. |
| `ld_prune.window` | Recommended | PLINK LD-pruning window. Values can be variant counts or strings such as `500kb`. |
| `ld_prune.step` | Recommended | PLINK LD-pruning step. |
| `ld_prune.r2` | Recommended | PLINK LD-pruning `r2`. |
| `pca.approx` | Recommended | Whether to use approximate PCA behavior where supported. |
| `pca.min_projection_pc_correlation` | Recommended | Minimum correlation between original and reprojected reference PCs during projection validation. |

Production package resolution injects `ancestry_reference.reference_genotypes`,
`ancestry_reference.metadata`, `ancestry_reference.reference_genome_build`,
`ancestry_reference.exclusion_regions`, and source/provenance fields into
`resolved_config.yaml`.

Main outputs:

```text
results/qc/ancestry/reference/reference_pcs.tsv
results/qc/ancestry/reference/study_projected_pcs.tsv
results/qc/ancestry/reference/reference_prep_report.md
results/qc/ancestry/within_ancestry_pcs.tsv
```

### `qc`

These settings drive PLINK2 GWAS variant and sample filters.

| Parameter | Required | Description |
| --- | --- | --- |
| `info_min` | Optional | Minimum imputation INFO/MACH_R2 threshold when `use_mach_r2_filter: true`. |
| `maf_min` | Recommended | Minimum GWAS minor allele frequency. |
| `hwe_p_min` | Recommended | Minimum Hardy-Weinberg p-value filter. |
| `geno_missing_max` | Recommended | Maximum per-variant missingness. |
| `sample_missing_max` | Recommended | Maximum per-sample missingness. |
| `use_mach_r2_filter` | Recommended | Set `true` for imputed dosage data with MACH_R2/INFO annotations. |
| `snps_only_acgt` | Recommended | Keep only A/C/G/T SNPs when `true`. |
| `autosome_only` | Recommended | Restrict GWAS to autosomes when `true`. |

For imputed dosage data with MACH_R2/INFO annotations, use:

```yaml
qc:
  use_mach_r2_filter: true
  info_min: 0.8
```

### `relatedness`

| Parameter | Required | Description |
| --- | --- | --- |
| `mode` | Yes | `plink2_king` or `all_samples`. Production should usually use `plink2_king`. |
| `king_cutoff` | When `plink2_king` | KING relatedness cutoff. Template uses `0.0884`. |
| `remove_sex_mismatches` | Legacy | Prefer `sex_check.action: "exclude"`. If this is `true`, validation requires `sex_check.action: "exclude"`. |
| `maf_min` | Recommended | Minimum MAF for relatedness marker preparation. |
| `geno_missing_max` | Recommended | Maximum marker missingness for relatedness marker preparation. |
| `sample_missing_max` | Recommended | Maximum sample missingness for relatedness marker preparation. |
| `snps_only_acgt` | Recommended | Keep only A/C/G/T SNPs when `true`. |
| `autosome_only` | Recommended | Keep autosomes only when `true`. |
| `ld_prune.window` | Recommended | PLINK LD-pruning window, for example `500kb`. |
| `ld_prune.step` | Recommended | PLINK LD-pruning step. |
| `ld_prune.r2` | Recommended | PLINK LD-pruning `r2`. |

Key outputs:

```text
results/qc/relatedness/relatedness_qc.pgen
results/qc/relatedness/relatedness_ld_prune.prune.in
results/qc/relatedness/unrelated.king.cutoff.in.id
results/qc/relatedness/relatedness_summary.tsv
```

### `sex_check`

| Parameter | Required | Description |
| --- | --- | --- |
| `enabled` | Production | `true` or `false`. Production requires `true`. |
| `action` | Yes | `warn`, `fail`, or `exclude`. Production requires `exclude`. |
| `allow_no_sex_markers` | Production review | If `false`, production fails when genotype data lack X/Y markers. Set `true` only with documented external sex QC. |
| `max_female_xf` | Optional | Custom PLINK `--check-sex` female X inbreeding threshold when nonblank. |
| `min_male_xf` | Optional | Custom PLINK `--check-sex` male X inbreeding threshold when nonblank. |
| `max_female_yrate` | Optional | Custom PLINK `--check-sex` female Y-rate threshold when nonblank. |
| `min_male_yrate` | Optional | Custom PLINK `--check-sex` male Y-rate threshold when nonblank. |

If all four threshold fields are blank, the workflow uses `max-female-xf=0.2`
and `min-male-xf=0.8`. Set cohort-specific thresholds after reviewing the
`XF` and `YRATE` distributions in `results/qc/sex/sexcheck.tsv`. Custom
settings must include both `max_female_xf` and `min_male_xf`; Y-rate thresholds
are optional but must be supplied as a pair.

Outputs:

```text
results/qc/sex/sexcheck.tsv
results/qc/sex/sex_mismatches.remove.tsv
results/qc/sex/sex_checked.keep.tsv
results/qc/sex/sex_check_summary.tsv
```

`sex_checked.keep.tsv` contains all samples for `action: "warn"` and excludes
problematic samples for `action: "exclude"`.

### `gwas`

| Parameter | Required | Description |
| --- | --- | --- |
| `test` | Recommended | PLINK2 test term retained in harmonized output. Template uses `ADD`. |
| `default_covariates` | Yes | Covariates used for every trait unless removed from the config. Non-PC covariates must exist in the sample manifest. PC covariates must exist in the active PC table. |
| `extra_covariates` | Optional | Additional global covariates used for every trait. |
| `covar_variance_standardize` | Recommended | Adds PLINK2 `--covar-variance-standardize` when `true`. |
| `allow_missing_pcs` | Production false | Allows missing PC covariates when `true`. Production forbids `true`. |
| `glm_options` | Optional | Free-form options appended after PLINK2 `--glm`, for example `hide-covar firth-fallback`. |

Use a centered quadratic age term for `age2` when possible, for example
`(age - mean_age)^2`, to reduce collinearity.

Trait-specific covariates from the trait registry `covariates` column are
appended to `default_covariates` and `extra_covariates` for that trait only.

### `warnings`

These thresholds are report-only warnings; they do not stop the workflow.

| Parameter | Required | Description |
| --- | --- | --- |
| `min_n` | Optional | Warn when an analyzed trait/ancestry sample count is below this value. |
| `min_cases` | Optional | Warn when case count is below this value. |
| `min_controls` | Optional | Warn when control count is below this value. |

The workflow fails trait/ancestry cells with zero cases or zero controls.

### `resources`

| Parameter | Required | Description |
| --- | --- | --- |
| `software_manifest` | Yes | TSV documenting software provenance. |
| `reference_manifest` | Yes | TSV documenting reference data/package provenance. |
| `input_manifest` | Optional | Optional compact cohort-release note TSV. Leave blank unless needed. |

`resources.input_manifest` is not the source of truth for input paths. The
config is the source of truth, and the workflow records configured paths in
`resolved_config.yaml`, reports, and `run_manifest.tsv`.

If `resources.input_manifest` is nonblank in production, it must contain:

```text
file_role	cohort_data_release	notes
```

Required `file_role` values:

```text
sample_manifest
trait_registry
study_genotype
```

### `runtime`

These values become Snakemake resources. The bundled SLURM profile can override
or map them to cluster-specific resource requests.

| Parameter | Required | Description |
| --- | --- | --- |
| `threads_small` | Yes | Threads for small PLINK/R helper jobs. |
| `threads_gwas` | Yes | Threads for each PLINK2 GWAS job. |
| `threads_admixture` | Recommended | Threads for ADMIXTURE jobs. Falls back to small threads in some ADMIXTURE rules if omitted. |
| `mem_mb_small` | Yes | Memory in MB for small helper jobs. |
| `mem_mb_gwas` | Yes | Memory in MB for each GWAS job. |
| `mem_mb_admixture` | Recommended | Memory in MB for ADMIXTURE jobs. |
| `time_min_small` | Yes | Runtime in minutes for small helper jobs. |
| `time_min_gwas` | Yes | Runtime in minutes for each GWAS job. |
| `time_min_admixture` | Recommended | Runtime in minutes for ADMIXTURE jobs. |

For HPC runs, also review `profiles/slurm/config.yaml`. To identify candidate
SLURM values on the cluster login node, run:

```bash
sinfo -o "%P %a %l %D %C"
sacctmgr -nP show assoc user=$USER format=Account,Partition,QOS,DefaultQOS
```

Use a partition that is available, has enough runtime for the job, and appears
in your user association. If `sinfo` marks a partition as `cpu*`, write `cpu`
in the profile; the `*` only marks the cluster default. Use only a QOS allowed
for the selected account and partition.

Set `conda-prefix` to a writable directory for Snakemake-created rule
environments. It is not the conda installation path. A directory in your home
folder is acceptable if it is accessible to compute nodes and has enough quota.

In `profiles/slurm/config.yaml`, set account and partition under
`default-resources` so they are attached to every submitted job:

```yaml
slurm-qos: "normal"
conda-prefix: "/path/to/shared/conda/envs"

default-resources:
  slurm_account: "your_account"
  slurm_partition: "standard"
  mem_mb: 4000
  runtime: 30
```

Leave `slurm-qos` commented out or remove it if your cluster does not use QOS.

## Production Validation Gates

The production-only workflow requires:

- `project.genome_build: "auto"`
- `ancestry_reference.enabled: true`
- `admixture.enabled: true`
- `reference_package.root` and `reference_package.fingerprint`
- a resolved package fingerprint that matches the expected fingerprint
- `sex_check.enabled: true`
- `sex_check.action: "exclude"`
- sex-chromosome markers unless `sex_check.allow_no_sex_markers: true`
- `gwas.allow_missing_pcs: false`

Run production validation through Snakemake:

```bash
snakemake --profile profiles/slurm results/qc/input_validation/validation.ok
```

Do not run `scripts/validate_config.R` directly against editable
`config/config.yaml`. Production validation must use the Snakemake target above
so the workflow can pass a build-resolved config. For script development, use
focused `scripts/test_*.R` tests; if manually debugging validation, pass a
build-resolved config.

## Main Config-Derived Outputs

```text
results/qc/genome_build/genome_build.txt
results/qc/genome_build/genome_build_marker_matches.tsv
results/config/effective_config.yaml
results/config/resolved_config.yaml
results/qc/input_validation/validation.ok
results/qc/traits/{trait}.pheno.tsv
results/qc/traits/{trait}.covar.tsv
results/qc/strata/{ancestry}.unrelated.keep.tsv
results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.plink2.glm.tsv
results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.gwas_filter_summary.tsv
results/reports/{trait}/{trait}.{ancestry}.{build}.report.md
results/manifests/run_manifest.tsv
```
