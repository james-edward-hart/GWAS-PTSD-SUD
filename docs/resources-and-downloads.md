# Resources And Downloads

This pipeline keeps software and reference downloads separate from analysis rules.

Production HPC deployments should prefer site-managed modules, conda
environments, and approved data-transfer procedures. The workflow does not
download HGDP+1KG reference data. For a cohort run, use the approved unpacked
reference package and its `content_fingerprint.sha256`; record package-level
provenance in `resources/manifests/reference_data.tsv`.

## Software

### PLINK2

Official source:

```text
https://www.cog-genomics.org/plink/2.0/
```

Use a Linux x86_64 build on HPC systems. The workflow uses PLINK2 for genotype input, filtering, relatedness pruning, PCA, projection, and `--glm` GWAS.

### ADMIXTURE

Official source:

```text
https://dalexander.github.io/admixture/download.html
```

The workflow can run supervised ADMIXTURE as a report-only QC branch when `admixture.enabled: true`. Stage 1 ancestry labels, strata, keep files, and GWAS covariates still come from POP-MaD PCA/Mahalanobis assignment.

### Snakemake

Official source:

```text
https://snakemake.readthedocs.io/
```

Install through conda/mamba on HPC unless the cluster already provides a tested module.

The bundled SLURM profile expects the Snakemake 8 SLURM executor plugin in the driver environment:

```bash
mamba env create -f envs/snakemake-driver.yaml
mamba activate gwas-stage1-driver
```

## Reference Data

### gnomAD HGDP + 1000 Genomes

Primary POP-MaD MatrixTable source:

```text
gs://gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/pca_results/unrelateds_without_outliers.mt
```

Optional metadata fallback:

```text
https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz
```

For a standard production run, configure only the unpacked package path and the
value in `content_fingerprint.sha256`. Package creation and source-object
provenance are managed separately from each cohort run.

The reference manifest uses these columns:

```text
resource	file_role	version_or_build	genome_build	source_url	local_path	download_date	checksum_algorithm	sha256	metadata_columns	preprocessing_notes
```

The manifest is still useful for documenting package-level provenance, but runtime validation relies on the package's own `file_manifest.tsv` and `content_fingerprint.sha256`.

To record a local file without any network access when maintaining provenance:

```bash
Rscript scripts/record_reference_resource.R \
  --resource gnomAD_HGDP_1KG_prepared_pgen \
  --file-role prepared_pgen \
  --version-or-build v3.1.2_prepared \
  --genome-build GRCh38 \
  --source-url https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg \
  --local-path reference-data/hgdp_1kg/hgdp_1kg.pgen \
  --download-date YYYY-MM-DD \
  --notes 'Converted to PLINK2 after allele/build QC and LD-region exclusions.'
```

Genome-build inference uses `resources/build_markers.tsv`. The bundled table
includes 1,057 rsID markers with GRCh36, GRCh37, and GRCh38 positions. GRCh37
and GRCh38 positions were fetched from Ensembl REST; GRCh36 positions were
preserved from the local HapMap example data. Rebuild methods and QC are in:

```text
resources/README.md
resources/build_marker_panel_qc.tsv
```

Maintainer rebuilds should archive the exact marker-panel build command, source URLs, and QC table with the reference package.

### Long-Range LD Regions

Starter long-range LD/problem-region files are stored in:

```text
resources/ancestry/long_range_ld_regions.GRCh37.tsv
resources/ancestry/long_range_ld_regions.GRCh38.tsv
```

These are used when `ancestry_reference.enabled: true` or `admixture.enabled: true`. Replace them with a cohort-approved region list when the production reference bundle has its own documented exclusions.

## Development Data

The helper script `scripts/download_test_data.sh` downloads the official ADMIXTURE HapMap3 sample archive:

```text
https://dalexander.github.io/admixture/hapmap3-files.tar.gz
```

The development data use public HapMap3 genome-wide array genotypes with a
generated random binary phenotype:

```bash
bash scripts/download_test_data.sh
Rscript scripts/prepare_hapmap3_fixture.R
```

The HapMap3 dataset is real public genotype data, but it is not real imputed
dosage data and is not a cohort analysis template. Public, individual-level,
genome-wide imputed human array datasets are usually controlled-access. For a
cohort run, point `config/config.yaml` at one genome-wide PGEN/PVAR/PSAM or
BED/BIM/FAM dataset and use the production workflow with a build-matched
reference package. The workflow records configured cohort input paths in the
resolved config, final reports, and run manifest.
