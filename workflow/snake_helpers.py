# Helper functions used by Snakemake while expanding the Stage 1 DAG.

import csv
import os


def workflow_error(message):
    try:
        from snakemake.exceptions import WorkflowError
    except Exception:
        raise RuntimeError(message)
    raise WorkflowError(message)


# Read small TSV config/manifests used while constructing the Snakemake DAG.
def read_tsv(path):
    if not os.path.exists(path):
        workflow_error(f"TSV input not found while building the DAG: {path}")
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def trait_ids(config):
    path = config.get("inputs", {}).get("trait_registry", "")
    if not path:
        workflow_error("config inputs.trait_registry is empty")
    rows = read_tsv(path)
    if not rows or "trait_id" not in rows[0]:
        workflow_error(f"trait registry is missing required column 'trait_id': {path}")
    ids = [row["trait_id"] for row in rows if row.get("trait_id", "")]
    if not ids:
        workflow_error(f"trait registry contains no non-empty trait_id values: {path}")
    return ids


# Track every manifest-listed package file so edits invalidate the resolved config.
def reference_package_inputs(config):
    root = config.get("reference_package", {}).get("root", "")
    if not root:
        return []

    required = [
        f"{root}/panel_manifest.tsv",
        f"{root}/file_manifest.tsv",
        f"{root}/content_fingerprint.sha256",
    ]
    manifest_path = f"{root}/file_manifest.tsv"
    if not os.path.exists(manifest_path):
        return required

    paths = list(required)
    for row in read_tsv(manifest_path):
        rel = row.get("path", "")
        normalized = rel.replace("\\", "/")
        if not rel or os.path.isabs(rel) or normalized.startswith("../") or "/../" in normalized:
            continue
        paths.append(os.path.join(root, rel))
    return sorted(set(paths))


# Add production ancestry-reference report only when that branch is enabled.
def ancestry_reference_targets(config):
    if config.get("ancestry_reference", {}).get("enabled", False) and config["inputs"]["ancestry_mode"] == "computed":
        return ["results/qc/ancestry/reference/reference_prep_report.md"]
    return []


# Add the report-only ADMIXTURE branch only when explicitly enabled.
def admixture_targets(config):
    if not config.get("admixture", {}).get("enabled", False):
        return []
    return [
        "results/qc/admixture/study_ancestry_proportions.tsv",
        "results/qc/admixture/reference_ancestry_proportions.tsv",
        "results/qc/admixture/popmad_admixture_comparison.tsv",
        "results/qc/admixture/admixture_run_summary.tsv",
        "results/qc/admixture/admixture_report.md",
    ]


# Provide validation inputs for the enabled ancestry-reference branch.
def ancestry_reference_validation_inputs(config, _wildcards):
    if not config.get("ancestry_reference", {}).get("enabled", False):
        return []
    if config.get("reference_package", {}).get("root", ""):
        return []

    settings = config.get("ancestry_reference", {})
    paths = []

    metadata = settings.get("metadata", {}).get("path")
    if metadata:
        paths.append(metadata)

    exclusion_regions = settings.get("exclusion_regions")
    if exclusion_regions:
        paths.append(exclusion_regions)

    reference = settings.get("reference_genotypes", {})
    prefix = reference.get("prefix", "")
    kind = reference.get("type", "").lower()
    if kind == "pgen":
        paths.extend([f"{prefix}.pgen", f"{prefix}.pvar", f"{prefix}.psam"])
    elif kind == "bed":
        paths.extend([f"{prefix}.bed", f"{prefix}.bim", f"{prefix}.fam"])

    return [path for path in paths if path]


# Provide validation inputs for the enabled ADMIXTURE QC branch.
def admixture_validation_inputs(config, _wildcards):
    if not config.get("admixture", {}).get("enabled", False):
        return []
    if config.get("reference_package", {}).get("root", ""):
        return []

    settings = config.get("admixture", {})
    paths = []

    metadata = settings.get("metadata", {}).get("path")
    if metadata:
        paths.append(metadata)

    exclusion_regions = settings.get("exclusion_regions")
    if exclusion_regions:
        paths.append(exclusion_regions)

    reference = settings.get("reference_genotypes", {})
    prefix = reference.get("prefix", "")
    kind = reference.get("type", "").lower()
    if kind == "pgen":
        paths.extend([f"{prefix}.pgen", f"{prefix}.pvar", f"{prefix}.psam"])
    elif kind == "bed":
        paths.extend([f"{prefix}.bed", f"{prefix}.bim", f"{prefix}.fam"])

    return [path for path in paths if path]


# Read the build chosen by the genome-build checkpoint.
def inferred_build(checkpoints):
    build_file = checkpoints.infer_genome_build.get().output.build
    with open(build_file) as handle:
        return handle.read().strip()


# Expand final report targets after the genome-build checkpoint completes.
def report_targets(checkpoints, traits, ancestries, wildcards):
    build = inferred_build(checkpoints)
    return [
        f"results/reports/{trait}/{trait}.{ancestry}.{build}.report.md"
        for trait in traits
        for ancestry in ancestries
    ]


# Keep production computed-ancestry outputs separate from toy/local outputs.
def popmad_output_dir(config):
    if config["inputs"]["ancestry_mode"] == "computed" and config.get("ancestry_reference", {}).get("enabled", False):
        return "results/qc/ancestry/production"
    return "results/qc/ancestry"


# Build paths inside the active POP-MaD output directory.
def popmad_path(config, name):
    return f"{popmad_output_dir(config)}/{name}"


# Resolve the ancestry assignment file for precomputed or computed mode.
def ancestry_file(config):
    if config["inputs"]["ancestry_mode"] == "computed":
        return popmad_path(config, "popmad_assignments.tsv")
    return config["inputs"]["ancestry_file"]


# Resolve the projected study PC file used by POP-MaD.
def ancestry_study_pcs_file(config, _wildcards=None):
    if config["inputs"]["ancestry_mode"] == "computed" and config.get("ancestry_reference", {}).get("enabled", False):
        return "results/qc/ancestry/reference/study_projected_pcs.tsv"
    return config["inputs"]["projected_pcs_file"]


# Resolve the reference PC file used by POP-MaD.
def ancestry_reference_pcs_file(config, _wildcards=None):
    if config["inputs"]["ancestry_mode"] == "computed" and config.get("ancestry_reference", {}).get("enabled", False):
        return "results/qc/ancestry/reference/reference_pcs.tsv"
    return config["inputs"]["reference_pcs_file"]


# Resolve the PC table that includes assigned ancestry labels.
def popmad_within_file(config):
    if config["inputs"]["ancestry_mode"] == "computed" and config.get("ancestry_reference", {}).get("enabled", False):
        return f"{popmad_output_dir(config)}/popmad_projected_pcs_with_ancestry.tsv"
    return "results/qc/ancestry/within_ancestry_pcs.tsv"


# Resolve the GWAS covariate PC file for precomputed, toy, or production mode.
def pcs_file(config):
    if config["inputs"]["ancestry_mode"] == "computed" and config.get("ancestry_reference", {}).get("enabled", False):
        return "results/qc/ancestry/within_ancestry_pcs.tsv"
    if config["inputs"]["ancestry_mode"] == "computed":
        return "results/qc/ancestry/study_pcs.tsv"
    return config["inputs"]["pcs_file"]
