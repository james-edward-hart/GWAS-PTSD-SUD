# Cohort ReMeta Export

## Purpose

This optional branch prepares the minimum cohort evidence needed for central,
marginal ReMeta gene tests:

1. Rare-variant single-variant score statistics in regenie HTP format.
2. Cohort-specific marginal LD among variants that can enter a protein-coding
   or splice mask.
3. A checksummed manifest that binds the data source, genome build, samples,
   parameters, resources, software, and artifacts.

The cohort pipeline does not compute gene p-values. Harmonized annotation,
external population-frequency filtering, masks, ReMeta `gene`, and ReMeta
`merge` remain central so every cohort uses one definition. The final central
output is one p-value per gene/mask/test combination.

## Scientific Contract

- Eligible data are WES hardcalls or imputed dosages. Basic array and WGS
  inputs are rejected by configuration validation.
- Every active trait has its own target PGEN and LD prefix, subset to that
  trait's exact Phase 2 model keep. Ordinary regenie IDs, rare-variant regenie
  IDs, the keep, and the target PSAM must be identical.
- Variant IDs are rewritten as build-specific `CHR:POS:REF:ALT`. Enabling the
  branch requires an explicit attestation that the source variants were split,
  left-normalized, and aligned to the inferred reference build before entering
  the workflow. The pipeline cannot scientifically infer that history from a
  PVAR alone.
- Explicit PLINK2 provisional-REF flags are rejected; the REF allele in CPRA
  must be reference-verified.
- The target resource is the union of autosomal GENCODE v50 basic
  protein-coding exons plus two splice bases. It is a superset, not a mask.
- Cohort filtering uses `MAC >= 1`, missingness, and—for imputed data—MaCH R2.
  It deliberately does not apply cohort MAF, because a variant common in one
  ancestry may still satisfy the centrally harmonized external rarity rule.
- ReMeta is always called with `--skip-buffer`. The LD files therefore contain
  marginal within-gene target LD only: no flanking variants and no conditional
  analysis payload.
- Every chromosome job receives an exact target-ID extract generated from the
  trait-specific PVAR. This prevents a terminal gene from acquiring the first
  variant on the next chromosome at a multichromosome reader boundary.
- WES runs erase dosages and calculate LD from hardcalls. Imputed runs retain
  dosages and use ReMeta `--use-dosages`.
- WES call-level quality control (for example PASS, depth, genotype quality,
  allele balance, contamination, and batch checks) must be completed before
  conversion to PGEN; those fields cannot be reconstructed here. This branch
  adds target missingness and MAC filters but is not a variant-calling QC
  workflow.

ReMeta models each cohort's observed LD; it does not itself correct ancestry
differences across cohorts. The central model and cohort inclusion strategy
remain responsible for ancestry heterogeneity.

ReMeta can reuse one cohort-wide LD matrix across traits. This pipeline instead
uses exact-trait samples as a deliberate, stricter provenance policy. That
policy increases computation and can make LD estimates noisier for small trait
samples; the manifest makes the choice and sample set explicit.

## Configuration

Enable the branch only after Phase 2 regenie is configured:

```yaml
tools:
  regenie: regenie
  remeta: remeta

phase2_regenie:
  enabled: true

remeta:
  enabled: true
  data_source: wes             # wes | imputed
  genotype_mode: hardcall      # hardcall for WES; dosage for imputed
  input_variants_normalized: true
  resource_root: resources/remeta
  min_mac: 1
  geno_missing_max: 0.05
  info_min: 0.8                # imputed only
  target_r2: 0.0001
```

For imputed cohorts use:

```yaml
  data_source: imputed
  genotype_mode: dosage
```

The input must be PGEN/PVAR/PSAM. Do not set the normalization attestation to
true until multiallelic variants have been split and alleles have been
left-normalized and reference-aligned against the inferred GRCh37 or GRCh38
FASTA. A standard upstream implementation is `bcftools norm -f BUILD.fa -m
-any`, followed by conversion that preserves REF/ALT and dosages.

## What Runs

For each active Phase 2 singleton trait group, the branch:

1. Intersects the PAN genotypes with the trait's phenotype/covariate-complete
   keep file and the bundled target intervals.
2. Applies source-specific variant QC and writes stable CPRA IDs.
3. Runs rare-variant regenie Step 2 with the existing group Step 1 predictions.
4. Partitions target IDs by chromosome, then runs ReMeta independently on
   chromosomes 1-22 with `--target-extract` and `--skip-buffer`.
5. Verifies ordinary/rare regenie sample IDs, target PSAM identity, CPRA
   consistency, all 22 three-file LD sets, HTP identity and counts,
   HTP-to-PVAR/LD membership, and gene/index consistency.
6. Writes the export manifest with SHA-256 hashes and byte counts.

When ReMeta is enabled, each final Phase 2 PAN report includes a **ReMeta LD
Target Coverage** section derived from that trait group's validated artifacts.
It reports target genes with indexed variants, unique QC-passing target-region
variants represented in the LD indexes, gene-variant assignments within their
declared gene spans, uncovered counts, and the no-buffer policy.

The group keep contains complete covariates and a valid phenotype for exactly
one trait. A status checkpoint schedules REGENIE and ReMeta only for active
traits. Skipped traits remain documented in reports and the manifest but never
receive fake target, HTP, LD, or validation artifacts.

All 22 target-ID partitions must be nonempty and together reproduce the target
PVAR exactly. A cohort with no retained target variant on an autosome fails
explicitly instead of emitting a synthetic LD set. ReMeta's native
`compute_ref_ld` logs are moved to `results/logs/remeta/`; the export tree is
reserved for declared, checksummed handoff artifacts.

The 22 LD jobs per active trait are chromosome-parallel, so sufficient cluster
capacity makes their wall time approximate the slowest chromosome rather than
their sum.
Runtime scales mainly with group sample count and retained coding/splice
variants. The template requests 4 threads, 32 GB, and 4 hours per chromosome;
benchmark chromosome 1 in the first cohort and tune the SLURM profile from its
observed peak memory and elapsed time. The bundled SLURM profile caps the custom
`remeta_ld_jobs` resource at 8, while every LD rule consumes one slot. Override
that capacity only after benchmarking the shared filesystem; this prevents an
88-job I/O burst when four traits are active.

## Export Layout

```text
results/remeta/export/{build}/htp/{trait}.PAN.regenie.gz
results/remeta/export/{build}/ld/{group}/chr1.remeta.gene.ld
results/remeta/export/{build}/ld/{group}/chr1.remeta.buffer.ld
results/remeta/export/{build}/ld/{group}/chr1.remeta.ld.idx.gz
...
results/remeta/export/{analysis_name}.{build}.remeta_manifest.tsv
```

The chromosome partitions are reusable work products at
`results/remeta/work/{build}/groups/{group}/ld_extracts/`. They are not part of
the central handoff.

Group names are stable type-and-trait IDs such as `bt__co_ptsd_aud`. The
`.remeta.buffer.ld` file is part of ReMeta's required three-file format but
contains no conditional buffer payload because `--skip-buffer` was used.

Share the HTP files, all three LD files for every chromosome/group, and the
manifest. The central analyst should reject an export if hashes, build,
normalization attestation, genotype mode, resource identity, or sample-group
metadata do not agree with the analysis specification. The manifest records
each trait's matching LD group, trait type, active/skip state, exact sample
metadata, and hashes, so
no separate cohort config is needed to pair an HTP file with its LD prefix.

## Central Handoff Boundary

Central analysis is deliberately not implemented in this repository. Its
external analysis specification must:

1. Normalize and harmonize CPRA/alleles across cohort HTP files.
2. Annotate once with the centrally pinned annotation source.
3. Define rarity from the agreed external population rule (for example,
   ancestry-aware maximum frequency), not independently within each cohort.
4. Build deterministic gene masks and keep a mask-version identifier.
5. Run ReMeta `gene` with each cohort's matching HTP and LD prefix.
6. Run ReMeta `merge` and report gene/mask/test p-values plus cohort and variant
   counts.

Because this cohort branch intentionally does not include conditional buffer
variants, central conditional gene tests must not be requested from these LD
exports.
