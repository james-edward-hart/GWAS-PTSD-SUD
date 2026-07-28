# Development Notes

Keep scripts small.

Each helper script should:

- do one job
- accept explicit inputs and outputs
- print clear errors
- avoid hidden writes
- leave enough white space to scan quickly

Prefer Snakemake rules and config over long shell scripts.

Use comments only where they explain non-obvious GWAS or workflow behavior.

Before changing workflow behavior:

```bash
snakemake -n --use-conda  # requires a production-style config/config.yaml
Rscript scripts/test_reference_package.R
Rscript scripts/test_infer_genome_build.R
Rscript scripts/test_genotype_input_sanitizer.R
Rscript scripts/test_sex_check.R
Rscript scripts/test_remeta_cohort.R
python scripts/test_remeta_resources.py
```

Full validation must run through Snakemake so the workflow can infer the genome
build, resolve the reference package, and validate
`results/config/resolved_config.yaml`:

```bash
snakemake --profile profiles/slurm results/qc/input_validation/validation.ok
```
