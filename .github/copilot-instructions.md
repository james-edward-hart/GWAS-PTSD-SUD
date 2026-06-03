# Copilot Instructions — GWAS-PTSD-SUD

## Project Overview

Cross-ancestry GWAS of PTSD × substance use disorder (SUD) comorbidity across five PGC (Psychiatric Genomics Consortium) cohorts of predominantly African-American participants. The pipeline goes from raw PGC phenotype files → harmonized phenotype files → genotype QC → logistic GWAS → post-GWAS diagnostics.

## Cohorts

| Code | Study | Genotyping array | Notes |
|------|-------|-----------------|-------|
| AAND | African-American Nicotine Dependence | Illumina OmniExpressExome | No opioid data |
| ADAA | Alcohol Dependence in African Americans | Custom Illumina OmniExpressExome | Has duplicate samples to remove |
| CGNDNICO | COGEND — Omni2.5M chip | Illumina Omni2.5M | 28 subjects overlap with CGNDSAGE (kept here) |
| CGNDSAGE | COGEND SAGE subset — Human1M chip | Illumina Human1M | Must remove 28 overlapping CGNDNICO subjects |
| FSCDSAGE | Family Study of Cocaine Dependence | Illumina Human1M | — |

## Pipeline Architecture

### Phase 1: Phenotype Preparation (`Data/{COHORT}/{COHORT}-pgc-prep.r`)

Each R script reads raw PGC phenotype files from an external drive, harmonizes SUD case/control definitions, creates comorbidity variables, and writes a standardized `{COHORT}_gwas_pheno.txt`. **Each cohort has idiosyncratic source-data quirks** — do not assume the recoding logic is identical across cohorts. Key differences:

- AAND: No opioid data at all; `EvrAlc` uses 1=never/2=ever coding
- ADAA: 40 duplicate FIDs to remove; ALC PHENO=0 is truly missing (no abuse-only split); opioid data uses `EverUseOpi` 0/1 coding
- CGNDSAGE: Must cross-reference CGNDNICO to remove 28 overlapping subjects
- FSCDSAGE: Sex is coded 2=female/3=male in source (not the usual 1/2)

### Phase 2: Ancestry Estimation (`Scripts/run_ancestry.sh`)

Supervised ADMIXTURE (K=5) against 1000 Genomes super-populations (AFR/AMR/EAS/EUR/SAS). Merges each cohort with all five reference panels, handles strand flips iteratively, LD-prunes, then runs ADMIXTURE in supervised mode. Output: `{COHORT}_ancestry_proportions.csv`.

### Phase 3: GWAS (`Scripts/gwas_template.sh`)

Single script for all cohort × phenotype combinations. Usage:
```bash
bash Scripts/gwas_template.sh COHORT PHENOTYPE
# Example: bash Scripts/gwas_template.sh AAND co_ptsd_aud
```

Steps: genotype QC → sample QC (heterozygosity ±3SD, sex check) → PCA (10 components) → logistic regression (covariates: sex, age, PC1–PC10) → genomic inflation λ + QQ plot.

## Phenotype Coding Conventions

- **Case/control**: `2` = case, `1` = control, `-9` = missing (PLINK convention)
- **PTSD**: `ptsd_dx` — `1` = yes, `0` = no (note: different scale than SUD variables)
- **Comorbidity** (`co_ptsd_*`): case = both PTSD + SUD; control = neither; all other combinations = `-9`
- **anySUD**: case if any individual SUD is case; control only if all assessed SUDs are control
- **`make_comorbid()`** helper function is defined identically in each prep script

## Phenotype Column Names

All `_gwas_pheno.txt` files share this schema (tab-delimited):
`FID, IID, ptsd_dx, trauma_exposed, aud_casecon, cud_casecon, tud_casecon, oud_casecon, anysud_casecon, co_ptsd_aud, co_ptsd_cud, co_ptsd_tud, co_ptsd_oud, co_ptsd_anysud, sex, age`

## Key Tools & Paths

- **PLINK2** for GWAS and QC; **PLINK 1.9** for ancestry (merge operations)
- **ADMIXTURE 1.3.0** for supervised ancestry estimation
- **R/Rscript** for phenotype prep, covariate building, sex checks, and QQ plots
- Genotype data lives on external drive (`/Volumes/Amstadter-GCPS-1/`), with optional local cache at `~/tmp_gwas/`
- Results are written locally then copied back to the external drive under `GWAS-Results/{COHORT}/`

## Conventions

- All bash scripts use `set -euo pipefail`
- Sex from `.fam` files (always 1=male, 2=female) is preferred over phenotype-file sex columns due to inconsistent coding across cohorts
- Ambiguous A/T and C/G SNPs are excluded before ancestry merges
- HWE filtering is applied to controls only
- GWAS covariates are always: sex, age, PC1–PC10 (variance-standardized)
