# Repository Guidelines

## Project Structure & Module Organization

This repository is a Snakemake pipeline for Stage 1 ancestry-stratified GWAS. The workflow entrypoint is `Snakefile`; rule files live in `workflow/rules/`, and shared Snakemake helpers live in `workflow/snake_helpers.py`. R and shell helper scripts are in `scripts/`, with reusable R code in `scripts/lib/`. Configuration and fixture inputs are in `config/` and `data/example/`. Conda environment definitions are in `envs/`, SLURM settings are in `profiles/slurm/`, reference manifests are in `resources/manifests/`, and generated outputs go under `results/`.

## Build, Test, and Development Commands

- `mamba env create -f envs/gwas.yaml`: create the local analysis environment.
- `bash scripts/download_test_data.sh`: download public HapMap3 example data.
- `Rscript scripts/prepare_hapmap3_fixture.R`: build the toy test fixture.
- `snakemake -n --use-conda`: dry-run the workflow and inspect planned jobs.
- `snakemake --cores 4 --shared-fs-usage input-output persistence software-deployment sources storage-local-copies`: run locally.
- `snakemake --profile profiles/slurm`: submit with the bundled SLURM profile.
- `Rscript scripts/test_pipeline_outputs.R`: validate expected example-run outputs.

## Coding Style & Naming Conventions

Keep helper scripts small and single-purpose. Prefer explicit command-line inputs and outputs, clear error messages, and no hidden writes outside declared paths. Use Snakemake rules and YAML config instead of long shell scripts when workflow behavior changes. R scripts use readable snake_case names such as `infer_genome_build.R`; Snakemake outputs follow `{trait}/{ancestry}/{trait}.{ancestry}.{build}.*` patterns. Use comments only for non-obvious GWAS, QC, or scheduler behavior.

## Testing Guidelines

Before changing workflow behavior, run:

```bash
snakemake -n --use-conda
Rscript scripts/validate_config.R --config config/config.yaml --out /tmp/stage1-validation.ok
```

Use existing `scripts/test_*.R` files for focused checks and keep new tests close to the script they validate. Test names should describe the behavior under test, for example `test_popmad_ancestry.R`.

## Commit & Pull Request Guidelines

This checkout does not include git history, so no existing commit convention can be inferred. Use short, imperative commit messages such as `Add sex-check exclusion report` or `Fix SLURM runtime resource mapping`. Pull requests should summarize pipeline behavior changes, list commands run, note config or resource-manifest updates, and include representative output paths when results change.

## Security & Configuration Tips

Do not commit protected cohort data, credentials, or external-drive paths containing private identifiers. For real analyses, start from `config/config.template.yaml`, record software and reference provenance in `resources/manifests/`, and keep production `results/` separate from public example outputs.
