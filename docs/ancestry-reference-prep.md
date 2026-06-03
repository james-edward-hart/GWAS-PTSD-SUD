# Ancestry Reference Preparation

Stage 1 supports two ancestry modes:

- `precomputed`: validated ancestry labels are supplied by the analyst.
- `computed`: POP-MaD assignment from projected study PCs and an HGDP + 1000 Genomes style reference PC table.

The runnable example uses computed mode with a small HapMap-derived fixture. Production use requires the reviewed, fingerprinted reference package containing POP-MaD and ADMIXTURE panels for the inferred study build.

Reference-package creation is intentionally outside the GWAS rules. The workflow consumes one local prebuilt package, validates its fingerprint, resolves build-matched panel paths, and does not download or build HGDP+1KG or 1000 Genomes reference data.

## Maintainer Notes

The `reference_package_builder/` directory is maintainer tooling for rebuilding the package when needed. Snakemake does not call it, and collaborators running Stage 1 should not need it. Routine users only need:

```text
reference_package.root: "/path/to/stage1_reference_package"
reference_package.fingerprint: "sha256-from-content_fingerprint.sha256"
```

The current production package was built upstream in a cluster environment and carries its own `panel_manifest.tsv`, `file_manifest.tsv`, `content_fingerprint.sha256`, QC reports, methods, and provenance files. The runtime pipeline treats those as read-only inputs.

Keep raw reference downloads and builder caches outside collaborator run directories. The unpacked package and its fingerprinted archive are the handoff artifacts.

## Preparation Steps

1. Validate study genotype build and reference build metadata.
2. Convert reference genotypes to one filtered PLINK2 dataset.
3. Intersect study and reference variants by ID, build, chromosome, position, and allele.
4. Keep autosomal biallelic SNPs and remove strand-ambiguous palindromic variants unless a reviewed allele-frequency check supports rescue.
5. Exclude long-range LD and other cohort-approved problematic regions.
6. For pre-LD-pruned POP-MaD package panels, copy the harmonized shared marker list to the stable PCA marker path; otherwise LD-prune variants in the reference panel.
7. Fit reference PCA with PLINK2 allele weights.
8. Project reference and study samples with the same PLINK2 `--score` command.
9. Assign ancestry using POP-MaD across 10 PCs: remove reference population outliers, calculate Mahalanobis distance to each population, assign the nearest reviewed label, and exclude ambiguous or outlying samples.
10. Run within-ancestry PCA in final sex-checked unrelated study strata for final GWAS covariates.
11. Export labels, excluded sample counts, population counts, assignment confidence metrics, PCs, and the reference-prep report.

## Workflow Inputs

Enable production reference use with:

```yaml
project:
  run_mode: "production"

reference_package:
  root: "/path/to/stage1_reference_package"
  fingerprint: "sha256-from-content_fingerprint.sha256"

ancestry_reference:
  enabled: true

admixture:
  enabled: true
```

Stage 1 resolves `ancestry_reference` and `admixture` reference paths from the package `panel_manifest.tsv` after genome-build inference. Package creation remains outside Snakemake.

Main generated files:

```text
results/qc/ancestry/reference/reference_qc.pgen
results/qc/ancestry/reference/study_qc.pgen
results/qc/ancestry/reference/shared_variants.txt
results/qc/ancestry/reference/shared_variant_mismatches.tsv
results/qc/ancestry/reference/ld_prune/ancestry_ld_prune.prune.in
results/qc/ancestry/reference/reference_pca/reference.eigenvec.allele
results/qc/ancestry/reference/reference_pcs.tsv
results/qc/ancestry/reference/study_projected_pcs.tsv
results/qc/ancestry/reference/reference_projection_validation.tsv
results/qc/ancestry/production/popmad_assignments.tsv
results/qc/ancestry/production/popmad_population_counts.tsv
results/qc/ancestry/production/population_model_summary.tsv
results/qc/ancestry/reference/reference_prep_report.md
results/qc/ancestry/within_ancestry_pcs.tsv
```

The pipeline is conservative: package fingerprint mismatch fails validation, build mismatch fails validation, raw Hail/VCF/BCF artifacts fail package validation, chromosome/position mismatches fail harmonization, strand-ambiguous palindromic SNPs are excluded by default, allele mismatches are reported rather than flipped silently, fewer than 10,000 shared POP-MaD variants fails harmonization, 10,000-49,999 shared variants warns and continues, reference projection is checked by comparing original reference PCs with reprojected reference samples, and POP-MaD model cutoffs are written for review.

ADMIXTURE is available as an independent report-only QC branch when `admixture.enabled: true`. It uses a separate 1000 Genomes reference configuration, writes proportions and POP-MaD comparison tables under `results/qc/admixture/`, and does not alter POP-MaD labels, strata, keep files, within-ancestry PCs, or GWAS covariates.

Do not run production GWAS from computed ancestry labels until the reference-preparation report and manifest have been reviewed.
