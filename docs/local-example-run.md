# HapMap3 Development Data

The HapMap3 files provide a small public dataset for focused development
checks, such as genome-build inference, script-level tests, ADMIXTURE QC checks,
and marker-panel maintenance.

Do not use these files as a template for cohort analysis. Full Stage 1 runs
require one cohort genotype dataset, a sample manifest (phenotype + covariate
file), a trait registry, a build-matched fingerprinted reference package, and
the production workflow described in
[cohort-production-quickstart.md](cohort-production-quickstart.md).

## Fixture Setup

Download the public HapMap3 PLINK files:

```bash
bash scripts/download_test_data.sh
```

Build the example sample manifest and trait registry:

```bash
mamba run -n gwas-stage1 Rscript scripts/prepare_hapmap3_fixture.R
```

The generated phenotype is random and not analytically meaningful.

## Focused Checks

Run focused tests through the utility environment:

```bash
mamba run -n gwas-stage1 Rscript scripts/test_infer_genome_build.R
mamba run -n gwas-stage1 Rscript scripts/test_reference_package.R
mamba run -n gwas-stage1 Rscript scripts/test_ancestry_reference.R
mamba run -n gwas-stage1 Rscript scripts/test_admixture_qc.R
```

Use the production quickstart for real cohort analyses.
