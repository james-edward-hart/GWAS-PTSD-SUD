# Resource Methods

## Genome Build Marker Panel

`resources/build_markers.tsv` is the offline marker-position table used by `scripts/infer_genome_build.R`.

The shipped panel has 1,057 unique rsID markers. Each marker has one position row for
GRCh36, GRCh37, and GRCh38, for 3,171 build-specific rows total.

The current panel was generated on 2026-05-22 with:

```bash
# Rebuild the genome-build marker panel from a reviewed candidate BIM.
Rscript scripts/build_genome_marker_panel.R \
  --candidate-bim data/example/hapmap3.bim \
  --candidate-build GRCh36 \
  --out resources/build_markers.tsv \
  --qc resources/build_marker_panel_qc.tsv \
  --markers-per-chrom 50 \
  --candidate-multiplier 3 \
  --batch-size 80
```

Method summary:

- Candidate rsIDs were sampled across autosomes from a real genome-wide HapMap3 array BIM.
- GRCh37 positions were pulled from `https://grch37.rest.ensembl.org/variation/homo_sapiens`.
- GRCh38 positions were pulled from `https://rest.ensembl.org/variation/homo_sapiens`.
- Markers were retained only when Ensembl returned one unique autosomal SNP mapping in both builds.
- Markers with near-identical GRCh37 and GRCh38 coordinates were removed because they do not help distinguish builds.
- Final markers were spread across all autosomes.
- GRCh36 positions were preserved from the candidate HapMap3 BIM for regression testing.

The generated QC table is `resources/build_marker_panel_qc.tsv`.

The production inference defaults require at least 50 matching markers, at least a 95%
match fraction, at least a 20-marker lead over the next-best build, and at least a 20%
lead fraction over the next-best build.

For a new production reference, rebuild this table from that project reference or a reviewed common-SNP candidate list, then archive the command, source URLs, date, and QC table with the analysis.
