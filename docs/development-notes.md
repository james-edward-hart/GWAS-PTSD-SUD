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
# Smoke-test workflow behavior after edits.
snakemake -n --use-conda
Rscript scripts/validate_config.R --config config/config.yaml --out /tmp/stage1-validation.ok
```
