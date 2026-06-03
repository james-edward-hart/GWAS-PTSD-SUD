# Local Example Run

This smoke test uses public HapMap3 genotype data and a generated random binary phenotype. It is useful for checking the DAG, rule environments, QC reports, GWAS harmonization, and plots. It is not a production analysis.

## 1. Create Environments

The driver environment runs Snakemake. The utility environment runs standalone R helper scripts used before or after the workflow.

```bash
mamba env create -f envs/snakemake-driver.yaml
mamba env create -f envs/gwas.yaml
mamba activate gwas-stage1-driver
```

## 2. Create The Local Config

`config/config.yaml` is ignored because production configs can contain private paths. For the public fixture, copy the tracked example config:

```bash
cp config/config.hapmap3.example.yaml config/config.yaml
```

The example config uses:

- `project.run_mode: "test"`
- `inputs.ancestry_mode: "computed"`
- `ancestry_reference.enabled: false`
- `admixture.enabled: false`
- trait `random_binary`
- ancestries `HMAP_A` and `HMAP_B`
- genotype prefix `data/example/hapmap3`

## 3. Download And Prepare The Fixture

Download the public HapMap3 PLINK files:

```bash
bash scripts/download_test_data.sh
```

Build the sample manifest, toy trait registry, study PCs, and toy reference PC table:

```bash
mamba run -n gwas-stage1 Rscript scripts/prepare_hapmap3_fixture.R --plink2 plink2
```

The script defaults to `software/local/plink2` for this workstation checkout. Passing `--plink2 plink2` uses the PLINK2 installed in the `gwas-stage1` environment. On HPC or Linux, use `--plink2 plink2`, `--plink2 software/bin/plink2`, or another executable path that matches the platform.

## 4. Dry-Run And Run Locally

```bash
snakemake -n --use-conda
```

```bash
snakemake --cores 4 --use-conda --shared-fs-usage input-output persistence software-deployment sources storage-local-copies
```

`--use-conda` lets Snakemake build and reuse rule-specific environments from `envs/`. If you omit it, the active shell environment must provide every dependency used by the rules.

## 5. Check The Outputs

```bash
mamba run -n gwas-stage1 Rscript scripts/test_pipeline_outputs.R
```

Expected main files include:

```text
results/qc/genome_build/genome_build.txt
results/qc/ancestry/popmad_assignments.tsv
results/qc/strata/strata_counts.tsv
results/qc/relatedness/relatedness_summary.tsv
results/qc/sex/sex_check_summary.tsv
results/gwas/random_binary/HMAP_A/
results/gwas/random_binary/HMAP_B/
results/plots/random_binary/HMAP_A/
results/plots/random_binary/HMAP_B/
results/reports/random_binary/
results/manifests/run_manifest.tsv
```

To restart a clean local example, move or remove `results/` and rerun Snakemake. Do not reuse local example outputs as production evidence.
