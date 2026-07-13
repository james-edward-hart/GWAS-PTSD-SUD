#!/usr/bin/env python3

import csv
import hashlib
import os
import re


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_provenance(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return {row["key"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}


for build in ("GRCh37", "GRCh38"):
    directory = os.path.join(ROOT, "resources", "remeta", build)
    provenance = read_provenance(os.path.join(directory, "provenance.tsv"))
    assert provenance["genome_build"] == build
    assert provenance["gencode_release"] == "50"
    assert provenance["splice_padding_bp"] == "2"
    assert provenance["source"].startswith("https://ftp.ebi.ac.uk/")
    assert 18000 <= int(provenance["gene_count"]) <= 21000
    assert provenance["gene_list_sha256"] == sha256_file(os.path.join(directory, "gene_list.tsv"))
    assert provenance["genes_sha256"] == sha256_file(os.path.join(directory, "genes.tsv"))
    assert provenance["target_regions_sha256"] == sha256_file(os.path.join(directory, "target_regions.bed"))

    gene_rows = []
    with open(os.path.join(directory, "gene_list.tsv"), encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            gene, chrom, start, end = line.rstrip("\n").split("\t")
            assert re.fullmatch(r"ENSG[0-9]+", gene), (build, line_number, gene)
            assert 1 <= int(chrom) <= 22
            assert 1 <= int(start) <= int(end)
            gene_rows.append((gene, int(chrom), int(start), int(end)))
    assert len(gene_rows) == int(provenance["gene_count"])
    assert len({row[0] for row in gene_rows}) == len(gene_rows)
    assert gene_rows == sorted(gene_rows, key=lambda row: (row[1], row[2], row[0]))

    with open(os.path.join(directory, "genes.tsv"), newline="", encoding="utf-8") as handle:
        gene_metadata = {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert set(gene_metadata) == {row[0] for row in gene_rows}
    for gene, chrom, start, end in gene_rows:
        row = gene_metadata[gene]
        assert (int(row["chrom"]), int(row["start"]), int(row["end"])) == (chrom, start, end)
        assert row["build"] == build

    previous = None
    region_count = 0
    with open(os.path.join(directory, "target_regions.bed"), encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            chrom, start, end = line.rstrip("\n").split("\t")
            region = (int(chrom), int(start), int(end))
            assert 1 <= region[0] <= 22 and 0 <= region[1] < region[2]
            if previous is not None:
                assert (region[0], region[1]) >= (previous[0], previous[1])
                assert region[0] != previous[0] or region[1] > previous[2]
            previous = region
            region_count += 1
    assert region_count == int(provenance["merged_target_region_count"])

manifest_path = os.path.join(ROOT, "resources", "manifests", "reference_data.tsv")
with open(manifest_path, newline="", encoding="utf-8") as handle:
    manifest = {row["genome_build"]: row for row in csv.DictReader(handle, delimiter="\t")
                if row["file_role"] == "remeta_gene_target_resource"}
for build in ("GRCh37", "GRCh38"):
    assert build in manifest
    row = manifest[build]
    local_path = os.path.join(ROOT, row["local_path"])
    assert row["checksum_algorithm"] == "sha256"
    assert row["sha256"] == sha256_file(local_path)

print("ReMeta resource tests passed")
