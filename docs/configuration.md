# Configuration Guide

Edit `config/config.yaml` for each project.

Keep config values explicit. A failed validation step is preferable to running a GWAS with the wrong samples, phenotype, or build label.

## Genome Build

`project.genome_build` must stay set to `auto`.

The workflow infers the build before expanding build-labelled outputs. It compares genotype variant IDs and positions to `genome_build.marker_file`, with a coordinate-position fallback for non-rsID BIM/PVAR variant IDs, and writes:

```text
# Genome-build inference outputs.
results/qc/genome_build/genome_build.txt
results/qc/genome_build/genome_build_marker_matches.tsv
```

The active merged Snakemake config is also written to:

```text
# Active config snapshot used by scripts.
results/config/effective_config.yaml
results/config/resolved_config.yaml
```

Every script receives this snapshot, so runs launched with `snakemake --configfile project.yaml` use the same config during DAG construction and script execution.
`effective_config.yaml` is written before build-matched package resolution. `resolved_config.yaml` is written after genome-build inference and is used by validation and downstream rules.

Input paths are also reported from the resolved config. `resources.input_manifest` is optional and is only for lightweight cohort-release notes; it is not the source of truth for sample, trait, or genotype paths.

The pipeline does not lift genotype coordinates or summary statistics.

The bundled production marker table contains 1,057 rsID markers with GRCh36, GRCh37, and GRCh38 positions. GRCh37 and GRCh38 positions come from Ensembl REST; GRCh36 positions are preserved from the local HapMap fixture for regression testing. Rebuild methods and QC output are documented in:

```text
# Marker-panel methods and QC files.
resources/README.md
resources/build_marker_panel_qc.tsv
```

Production defaults require:

```yaml
# Genome-build thresholds; project.genome_build options: auto only.
genome_build:
  min_markers: 50
  min_match_fraction: 0.95
  min_marker_margin: 20
  min_fraction_margin: 0.20
```

When rsIDs are unavailable, the inference falls back to chromosome-position matching against the same marker table. Sparse datasets or datasets with too few informative positions should still fail loudly rather than guessing a build.

## Genotype Inputs

Supported values:

- `pgen`
- `bed`

For `pgen`, set:

```yaml
# PLINK2 genotype input; type options: pgen, bed.
genotypes:
  type: "pgen"
  prefix: "data/my_project/my_data"
```

For `bed`, set:

```yaml
# PLINK1 genotype input; type options: pgen, bed.
genotypes:
  type: "bed"
  prefix: "data/my_project/my_data"
```

## Phenotypes

The sample manifest must contain:

- `FID`
- `IID`
- `age`
- `age2`
- `sex`
- all configured phenotype columns
- all configured non-PC covariates

Case/control traits are converted to PLINK2 coding:

- control: `1`
- case: `2`
- missing: `NA`

`FID/IID` rows must be unique and must exist in the genotype files. Non-PC covariates must be numeric except for missing values. Sex codes must be `1`, `2`, `0`, `NA`, `-9`, or `.`. Phenotype values must match each trait registry's `case_value`, `control_value`, or `missing_values`.

## Covariates

The default covariates are:

- `age`
- `age2`
- `sex`
- `PC1` through `PC10`

Use a centered quadratic age term for `age2` when possible, for example `(age - mean_age)^2`, to avoid unnecessary collinearity.

Add project-specific covariates in `gwas.extra_covariates`.

Trait-specific covariates can be added with an optional comma-separated `covariates` column in the trait registry. These are appended to the default and extra covariates for that trait only.

By default, PLINK2 covariates are variance-standardized before regression:

```yaml
# Covariate standardization toggle; options: true, false.
gwas:
  covar_variance_standardize: true
```

## Ancestry

The pipeline supports `computed` and `precomputed` ancestry modes.

Production runs must use the prebuilt unpacked reference package:

```yaml
project:
  run_mode: "production"

inputs:
  ancestry_mode: "computed"

reference_package:
  root: "/path/to/stage1_reference_package"
  fingerprint: "sha256-from-content_fingerprint.sha256"
```

Production mode forbids precomputed ancestry, requires computed POP-MaD ancestry, requires report-only ADMIXTURE, checks the package content fingerprint, and requires reference/exclusion-region builds to match the inferred study build. The pipeline consumes this package as read-only data; it does not create or rebuild it.

For computed ancestry, provide projected study PCs and HGDP+1KG-style reference PCs:

```yaml
# Preprojected computed ancestry inputs; ancestry_mode options: computed, precomputed.
inputs:
  ancestry_mode: "computed"
  projected_pcs_file: "config/study_projected_pcs.tsv"
  reference_pcs_file: "config/hgdp_1kg_reference_pcs.tsv"
```

These two files are consumed only for non-package/test computed ancestry. In production package-backed computed ancestry, set `ancestry_reference.enabled: true`; the workflow creates run-specific projected study PCs and reference PCs from the resolved reference package instead of consuming the placeholder `projected_pcs_file` and `reference_pcs_file` values.

The reference PC file must include:

- `FID`
- `IID`
- `population`
- `super_population`
- `PC1` through `PC10`

POP-MaD assigns each study sample to the nearest reference population, collapses that to the configured ancestry label, and excludes ambiguous or outlying samples from ancestry-stratified GWAS. Assignment counts are written to:

```text
# POP-MaD assignment count outputs.
results/qc/ancestry/popmad_population_counts.tsv
results/qc/ancestry/popmad_excluded.tsv
results/qc/ancestry/population_model_summary.tsv
```

For production computed ancestry, enable the reference projection branch:

```yaml
# Production ancestry reference projection; reference paths are resolved from reference_package.
ancestry_reference:
  enabled: true
  warn_shared_variants_below: 50000
  min_shared_variants: 10000
```

When enabled, the pipeline creates run-specific projection and QC outputs instead of consuming precomputed projected PCs:

```text
# Production computed-ancestry outputs.
results/qc/ancestry/reference/reference_pcs.tsv
results/qc/ancestry/reference/study_projected_pcs.tsv
results/qc/ancestry/production/popmad_assignments.tsv
results/qc/ancestry/production/popmad_population_counts.tsv
results/qc/ancestry/reference/reference_prep_report.md
results/qc/ancestry/within_ancestry_pcs.tsv
```

The reference build must match the inferred study build. Reference metadata column names are configurable because HGDP+1KG releases and local manifests may use different labels. Package POP-MaD panels marked `variant_set=pre_ld_pruned` skip workflow-level LD pruning; the harmonized shared marker list is copied to the downstream PCA marker path.

In production computed mode, within-ancestry PCs are fitted on the final sex-checked unrelated samples for each ancestry stratum. These PCs are the GWAS covariates used in that stratum.

ADMIXTURE can be enabled as a separate report-only QC branch:

```yaml
# Report-only supervised ADMIXTURE QC.
admixture:
  enabled: true
  mode: "supervised"
  k: 5
  labels: [AFR, AMR, EAS, EUR, SAS]
```

When enabled, the workflow writes:

```text
results/qc/admixture/study_ancestry_proportions.tsv
results/qc/admixture/reference_ancestry_proportions.tsv
results/qc/admixture/popmad_admixture_comparison.tsv
results/qc/admixture/admixture_run_summary.tsv
results/qc/admixture/admixture_report.md
```

ADMIXTURE outputs are for QC review only. POP-MaD remains the ancestry-label source for strata and GWAS covariates.

For precomputed ancestry labels:

```yaml
# Precomputed ancestry input; ancestry_mode options: computed, precomputed.
inputs:
  ancestry_mode: "precomputed"
  ancestry_file: "config/my_ancestry.tsv"
```

The ancestry file must include:

- `FID`
- `IID`
- `ancestry`

Ancestry labels must match `analysis.ancestries`.

Samples with missing, ambiguous, or unconfigured ancestry labels are excluded from stratum keep files. The workflow fails when the excluded/unassigned fraction exceeds `popmad.max_unassigned_fraction`. In production mode, each configured trait-by-ancestry cell must contain at least one case and one control.

## Exclusion Region Files

`ancestry_reference.exclusion_regions` and `admixture.exclusion_regions` are TSV files used to remove long-range LD or other problem regions before ancestry-reference projection or ADMIXTURE QC.

Required columns:

```text
chrom	start	end	label
```

An optional `build` column can be included:

```text
chrom	start	end	label	build
6	25000000	34000000	MHC	GRCh38
```

In production, the file must either include a `build` column or include a build label such as `GRCh37` or `GRCh38` in the filename. The build must match the inferred study genotype build.

## Relatedness

The default is production pruning with PLINK2 KING:

```yaml
# Relatedness mode options: plink2_king, all_samples.
relatedness:
  mode: "plink2_king"
  king_cutoff: 0.0884
```

PLINK2 writes unrelated samples to:

```text
# KING unrelated keep output.
results/qc/relatedness/unrelated.king.cutoff.in.id
```

Before KING, the workflow creates an autosomal, biallelic, common, LD-pruned marker set:

```text
# Relatedness QC marker outputs.
results/qc/relatedness/relatedness_qc.pgen
results/qc/relatedness/relatedness_ld_prune.prune.in
results/qc/relatedness/relatedness_summary.tsv
```

Tune these filters under `relatedness.maf_min`, `relatedness.geno_missing_max`, and `relatedness.ld_prune`.

## Genetic Sex Check

The pipeline runs PLINK2 `--check-sex` when sex-chromosome markers are present. By default it reports problems but does not remove samples:

```yaml
# Genetic sex check; action options: warn, fail, exclude.
sex_check:
  enabled: true
  action: "warn"
  allow_no_sex_markers: false
  max_female_xf: ""
  min_male_xf: ""
  max_female_yrate: ""
  min_male_yrate: ""
```

Outputs are:

```text
# Genetic sex-check outputs.
results/qc/sex/sexcheck.tsv
results/qc/sex/sex_mismatches.remove.tsv
results/qc/sex/sex_checked.keep.tsv
results/qc/sex/sex_check_summary.tsv
```

Use `action: fail` to stop when mismatches are detected. Use `action: exclude` to remove problematic samples from relatedness pruning and final GWAS keep files. Production mode requires `action: exclude` and sex-chromosome markers unless `sex_check.allow_no_sex_markers: true` is set explicitly after external sex QC has already been completed and documented.
