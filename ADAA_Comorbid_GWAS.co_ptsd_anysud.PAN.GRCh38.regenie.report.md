# Phase 2 PAN Regenie Report: ADAA_Comorbid_GWAS / co_ptsd_anysud

## Model Overview

- Engine: regenie
- Trait type: bt
- Genome build: GRCh38
- Covariates: age,age2,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10

## PAN Sample Set

- PAN samples after sex-check and sample missingness: 1947
- POP-MaD assigned samples in PAN set: 1905
- POP-MaD UNKNOWN samples in PAN set: 42
- Complete covariate samples: 1947
- Usable trait samples: 625
- Cases: 258
- Controls: 367

## Variant Sources and QC

- Stage 1 union-pass variants: 14266172
- Pooled missingness threshold: 0.05
- Step 2 pooled MAF minimum: 0.01
- Regenie minMAC: 1
- Regenie minINFO: 0.9

## ReMeta LD Target Coverage

Coverage is calculated from the group-specific LD indexes and their matching QC-filtered target PGEN. A target variant can be assigned to more than one overlapping gene, so unique variants and gene-variant assignments are reported separately.

| Coverage metric | Observed | Denominator | Coverage |
| --- | ---: | ---: | ---: |
| Target genes with at least one indexed LD variant | 18728 | 19166 | 97.71% |
| QC-passing target-region variants represented in LD indexes | 1015489 | 1015489 | 100.00% |
| Indexed gene-variant assignments within declared gene spans | 1183144 | 1183179 | 100.00% |

- Target genes without an indexed LD variant: 438
- Target-region variants absent from every LD gene index: 0
- Indexed gene-variant assignments outside the declared gene span: 35
- Conditional buffer variants: not included (`--skip-buffer`; marginal LD export).

## REGENIE Run Settings

- Global PCs: 10
- Step 1 block size: 1000
- Step 2 block size: 400
- Regenie HTP cohort: ADAA_Comorbid_GWAS
- Binary approximate Firth pThresh: 0.01
- Quantitative RINT: False

## Association Results

- Native regenie output: `results/gwas/co_ptsd_anysud/PAN/co_ptsd_anysud.PAN.GRCh38.regenie`
- Skipped: False
- Native regenie variant rows: 13819687
- Valid P-value variants: 13819687
- Lambda GC: 0.855438
- Genome-wide significant variants (P <= 5e-8): 0
- Suggestive variants (P <= 1e-5): 0

## Top Hits

| CHROM | POS | ID | EFFECT | SE | P |
| --- | ---: | --- | ---: | ---: | ---: |
| 13 | 29736781 | 13:29736781:G:A | 0.0677836 | NA | 2.321e-05 |
| 16 | 63963509 | 16:63963509:C:G | 15.7835 | NA | 2.416e-05 |
| 4 | 56884691 | 4:56884691:T:G | 3.992 | NA | 2.488e-05 |
| 9 | 93935234 | 9:93935234:G:A | 20.224 | NA | 2.612e-05 |
| 9 | 118335011 | 9:118335011:T:C | 40.2181 | NA | 3.057e-05 |
| 9 | 118310119 | 9:118310119:C:T | 40.2111 | NA | 3.06e-05 |
| 21 | 22436995 | 21:22436995:A:C | 46.6612 | NA | 3.112e-05 |
| 17 | 6478316 | 17:6478316:T:C | 13.9217 | NA | 3.229e-05 |
| 12 | 23936388 | 12:23936388:T:C | 3.05754 | NA | 3.773e-05 |
| 12 | 97395899 | 12:97395899:C:T | 6.95493 | NA | 3.82e-05 |

## Plots

![QQ plot](../../plots/co_ptsd_anysud/PAN/ADAA_Comorbid_GWAS.co_ptsd_anysud.PAN.GRCh38.regenie.qq.png)

![Manhattan plot](../../plots/co_ptsd_anysud/PAN/ADAA_Comorbid_GWAS.co_ptsd_anysud.PAN.GRCh38.regenie.manhattan.png)

- QQ plot: `results/plots/co_ptsd_anysud/PAN/ADAA_Comorbid_GWAS.co_ptsd_anysud.PAN.GRCh38.regenie.qq.png`
- Manhattan PNG: `results/plots/co_ptsd_anysud/PAN/ADAA_Comorbid_GWAS.co_ptsd_anysud.PAN.GRCh38.regenie.manhattan.png`
- Manhattan PDF: `results/plots/co_ptsd_anysud/PAN/ADAA_Comorbid_GWAS.co_ptsd_anysud.PAN.GRCh38.regenie.manhattan.pdf`
- Large PAN plot fallback: QQ displayed 6909844 of 13819687 ordered P-value points (every second rank plus all P <= 1e-5); Manhattan displayed 6909844 of 13819687 eligible variants (the most significant 50%) because the 10000000-variant threshold was exceeded. Association metrics and top hits use all valid variants.

## Stage 1 Lambda Comparison

| Stage 1 ancestry | Lambda GC | Valid P variants |
| --- | ---: | ---: |
| AFR | 1.028676 | 14266169 |
