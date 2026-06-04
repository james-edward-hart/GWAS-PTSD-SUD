# HapMap3 Development Fixture

The old HapMap3 end-to-end example has been retired. The production pipeline no
longer supports `project.run_mode: "test"`, precomputed ancestry labels, or
user-supplied projected PC files.

Use this fixture only for focused development checks, such as genome-build
inference, script-level tests, and marker-panel maintenance. Real end-to-end
pipeline runs require a build-matched, fingerprinted reference package and the
production workflow described in
[cohort-production-run-manual.md](cohort-production-run-manual.md).

## Fixture Setup

Download the public HapMap3 PLINK files:

```bash
bash scripts/download_test_data.sh
```

Build the sample manifest, toy trait registry, and local PCA fixture files:

```bash
mamba run -n gwas-stage1 Rscript scripts/prepare_hapmap3_fixture.R --plink2 plink2
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

Do not use this fixture as production evidence and do not expect it to exercise
the full package-backed POP-MaD/GWAS DAG.
