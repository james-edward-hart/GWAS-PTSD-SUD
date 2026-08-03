# Helper functions used by Snakemake while expanding the Stage 1 DAG.

import csv
import math
import os
import re


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
    ids = [str(row.get("trait_id", "")) for row in rows]
    if any(trait_id != trait_id.strip() for trait_id in ids):
        workflow_error(f"trait registry contains trait_id values with surrounding whitespace: {path}")
    if not ids or any(not trait_id for trait_id in ids):
        workflow_error(f"trait registry contains no non-empty trait_id values: {path}")
    duplicates = sorted({trait_id for trait_id in ids if ids.count(trait_id) > 1})
    if duplicates:
        workflow_error(f"trait registry contains duplicate trait_id values: {', '.join(duplicates)}")
    unsafe = [trait_id for trait_id in ids if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", trait_id)]
    if unsafe:
        workflow_error(
            "trait_id values must start with a letter or number and contain only letters, "
            f"numbers, dots, underscores, and hyphens: {', '.join(unsafe)}"
        )
    return ids


def blank(value):
    return value is None or str(value).strip() == ""


def split_csv(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [str(item).strip() for item in value if str(item).strip()]
    items = [item.strip() for item in str(value).split(",")]
    return [item for item in items if item]


def config_list(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [str(item) for item in value if str(item)]
    value = str(value)
    return [value] if value else []


def truthy(value):
    return value is True or str(value) in {"true", "True", "1"}


def phase2_enabled(config):
    return truthy(config.get("phase2_regenie", {}).get("enabled", False))


def remeta_enabled(config):
    return truthy(config.get("remeta", {}).get("enabled", False))


def phase2_default_covariates(config):
    settings = config.get("phase2_regenie", {})
    covars = config_list(settings.get("default_covariates"))
    if covars:
        return covars
    pc_count = int(settings.get("global_pcs", 10))
    return ["age", "age2", "sex"] + [f"PC{i}" for i in range(1, pc_count + 1)]


def phase2_covariates_for_trait(config, trait_row):
    settings = config.get("phase2_regenie", {})
    covars = phase2_default_covariates(config)
    covars.extend(config_list(settings.get("extra_covariates")))
    covars.extend(split_csv(trait_row.get("covariates", "")))
    out = []
    for covar in covars:
        if covar and covar not in out:
            out.append(covar)
    return out


def phase2_trait_type(samples, trait_row):
    case_value = str(trait_row.get("case_value", "")).strip()
    control_value = str(trait_row.get("control_value", "")).strip()
    if case_value and control_value:
        return "bt"
    if case_value or control_value:
        workflow_error(
            f"trait {trait_row.get('trait_id', '<unknown>')} has only one of case_value/control_value set"
        )

    phenotype_column = trait_row.get("phenotype_column", "")
    if not samples or phenotype_column not in samples[0]:
        workflow_error(f"phenotype column '{phenotype_column}' is absent from sample manifest")
    missing = set(split_csv(trait_row.get("missing_values", "")))
    values = []
    for sample in samples:
        value = str(sample.get(phenotype_column, "")).strip()
        if value in missing or value in {"", "NA", "-9", "."}:
            continue
        values.append(value)
    if not values:
        workflow_error(f"quantitative trait {trait_row.get('trait_id', '<unknown>')} has no nonmissing values")
    for value in values:
        try:
            numeric = float(value)
        except ValueError:
            workflow_error(
                f"trait {trait_row.get('trait_id', '<unknown>')} has blank case/control values "
                f"but nonnumeric phenotype value '{value}'"
            )
        if not math.isfinite(numeric):
            workflow_error(
                f"trait {trait_row.get('trait_id', '<unknown>')} has non-finite phenotype value '{value}'"
            )
    return "qt"


def phase2_trait_groups(config):
    if not phase2_enabled(config):
        return []
    trait_path = config.get("inputs", {}).get("trait_registry", "")
    sample_path = config.get("inputs", {}).get("sample_manifest", "")
    traits = read_tsv(trait_path)
    samples = read_tsv(sample_path)
    expected_ids = trait_ids(config)
    rows = []
    for trait in traits:
        trait_id = str(trait.get("trait_id", ""))
        trait_type = phase2_trait_type(samples, trait)
        covars = phase2_covariates_for_trait(config, trait)
        # Group identity is trait-derived so registry ordering cannot change a model's sample set.
        group_id = f"{trait_type}__{trait_id}"
        rows.append(
            {
                "group": group_id,
                "trait_type": trait_type,
                "covariates": list(covars),
                "traits": [trait_id],
            }
        )
    if [row["traits"][0] for row in rows] != expected_ids:
        workflow_error("Phase 2 singleton groups do not match the trait registry")
    group_ids = [row["group"] for row in rows]
    if len(group_ids) != len(set(group_ids)):
        workflow_error("Phase 2 singleton group IDs are not unique")
    return rows


def phase2_group_ids(config):
    return [row["group"] for row in phase2_trait_groups(config)]


def phase2_group_row(config, group):
    for row in phase2_trait_groups(config):
        if row["group"] == group:
            return row
    workflow_error(f"unknown Phase 2 regenie trait group: {group}")


def phase2_trait_group(config, trait):
    for row in phase2_trait_groups(config):
        if trait in row["traits"]:
            return row["group"]
    workflow_error(f"trait is not assigned to a Phase 2 regenie group: {trait}")


def phase2_group_traits(config, group):
    return phase2_group_row(config, group)["traits"]


def analysis_output_name(config):
    name = str(config.get("project", {}).get("analysis_name", "")).strip()
    if not name:
        workflow_error("config project.analysis_name is empty")
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._-")
    if not safe:
        workflow_error("config project.analysis_name must contain at least one letter or number")
    return safe


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


# Add ancestry-reference report only when that branch is enabled.
def ancestry_reference_targets(config):
    if config.get("ancestry_reference", {}).get("enabled", False):
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


# Read active GWAS strata after strata creation has reviewed assigned ancestry counts.
def active_ancestries(checkpoints, wildcards):
    active_path = checkpoints.make_strata_files.get().output.active
    rows = read_tsv(str(active_path))
    if not rows or "ancestry" not in rows[0]:
        workflow_error(f"active ancestry table is empty or missing column 'ancestry': {active_path}")
    ancestries = [row["ancestry"] for row in rows if row.get("ancestry", "")]
    if not ancestries:
        workflow_error(f"no active ancestry strata were written: {active_path}")
    return ancestries


# Build active final-GWAS keep paths after sex-check and relatedness filters.
def active_unrelated_keep_files(checkpoints, wildcards):
    return [
        f"results/qc/strata/{ancestry}.unrelated.keep.tsv"
        for ancestry in active_ancestries(checkpoints, wildcards)
    ]


def active_within_ancestry_eigenvecs(checkpoints, wildcards):
    return [
        f"results/qc/ancestry/within/{ancestry}.eigenvec"
        for ancestry in active_ancestries(checkpoints, wildcards)
    ]


def active_admixture_stratum_outputs(checkpoints, wildcards, filename):
    return [
        f"results/qc/admixture/by_ancestry/{ancestry}/{filename}"
        for ancestry in active_ancestries(checkpoints, wildcards)
    ]


# Expand final report targets after genome-build and strata checkpoints complete.
def report_targets(checkpoints, traits, wildcards, config):
    build = inferred_build(checkpoints)
    ancestries = active_ancestries(checkpoints, wildcards)
    analysis_name = analysis_output_name(config)
    return [
        f"results/reports/{trait}/{analysis_name}.{trait}.{ancestry}.{build}.report.md"
        for trait in traits
        for ancestry in ancestries
    ]


def phase2_report_targets(checkpoints, traits, wildcards, config):
    if not phase2_enabled(config):
        return []
    build = inferred_build(checkpoints)
    # Force strata checkpoint completion because Phase 2 consumes Stage 1 QC.
    active_ancestries(checkpoints, wildcards)
    analysis_name = analysis_output_name(config)
    return [
        f"results/reports/{trait}/{analysis_name}.{trait}.PAN.{build}.regenie.report.md"
        for trait in traits
    ]


def phase2_group_stage1_stats(checkpoints, wildcards, config):
    build = inferred_build(checkpoints)
    ancestries = active_ancestries(checkpoints, wildcards)
    traits = phase2_group_traits(config, wildcards.group)
    return [
        f"results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.plink2.glm.tsv"
        for trait in traits
        for ancestry in ancestries
    ]


def phase2_trait_stage1_summaries(checkpoints, wildcards, config):
    build = inferred_build(checkpoints)
    ancestries = active_ancestries(checkpoints, wildcards)
    return [
        f"results/gwas/{wildcards.trait}/{ancestry}/{wildcards.trait}.{ancestry}.{build}.gwas_filter_summary.tsv"
        for ancestry in ancestries
    ]


def phase2_status_rows(checkpoints):
    path = str(checkpoints.select_phase2_active_groups.get().output.status)
    rows = read_tsv(path)
    required = {
        "group", "trait", "trait_type", "usable_n", "model_sample_count",
        "keep_count", "model_keep_sha256", "skipped", "remeta_eligible",
    }
    if not rows or not required.issubset(rows[0]):
        workflow_error(f"Phase 2 trait-group status table is empty or incomplete: {path}")
    groups = [row["group"] for row in rows]
    traits = [row["trait"] for row in rows]
    if len(groups) != len(set(groups)) or len(traits) != len(set(traits)):
        workflow_error(f"Phase 2 trait-group status table contains duplicate groups or traits: {path}")
    if any(
        row["skipped"] not in {"True", "False"}
        or row["remeta_eligible"] not in {"True", "False"}
        or row["skipped"] == row["remeta_eligible"]
        for row in rows
    ):
        workflow_error(f"Phase 2 trait-group status table has inconsistent active/skipped flags: {path}")
    return rows


def phase2_trait_status(checkpoints, trait):
    matches = [row for row in phase2_status_rows(checkpoints) if row["trait"] == trait]
    if len(matches) != 1:
        workflow_error(f"trait is absent or duplicated in the Phase 2 status table: {trait}")
    return matches[0]


def phase2_active_group_rows(checkpoints):
    return [
        row for row in phase2_status_rows(checkpoints)
        if row.get("remeta_eligible") == "True" and row.get("skipped") == "False"
    ]


def phase2_status_file(checkpoints):
    return str(checkpoints.select_phase2_active_groups.get().output.status)


def phase2_trait_regenie_outputs(checkpoints, wildcards, config):
    status = phase2_trait_status(checkpoints, wildcards.trait)
    # Empty dynamic input is intentional: skipped traits stage an honest public placeholder only.
    if status["skipped"] == "True":
        return []
    group = phase2_trait_group(config, wildcards.trait)
    prefix = f"results/gwas/PAN/regenie/groups/{group}/{group}.{wildcards.build}_{wildcards.trait}"
    return [f"{prefix}.regenie", f"{prefix}.regenie.ids", f"{prefix}.step2.done"]


def remeta_resource_file(config, build, filename):
    root = str(config.get("remeta", {}).get("resource_root", "resources/remeta")).rstrip("/")
    return f"{root}/{build}/{filename}"


def remeta_export_manifest(config, build):
    return f"results/remeta/export/{analysis_output_name(config)}.{build}.remeta_manifest.tsv"


def remeta_targets(checkpoints, wildcards, config):
    if not remeta_enabled(config):
        return []
    if not phase2_enabled(config):
        workflow_error("remeta.enabled requires phase2_regenie.enabled")
    build = inferred_build(checkpoints)
    if build not in {"GRCh37", "GRCh38"}:
        workflow_error(f"remeta.enabled requires GRCh37 or GRCh38; inferred build is {build}")
    return [remeta_export_manifest(config, build)]


def remeta_manifest_inputs(checkpoints, config, build):
    if not remeta_enabled(config):
        return []
    paths = []
    for row in phase2_active_group_rows(checkpoints):
        group = row["group"]
        for chrom in range(1, 23):
            prefix = f"results/remeta/export/{build}/ld/{group}/chr{chrom}"
            paths.extend(
                [
                    f"{prefix}.remeta.gene.ld",
                    f"{prefix}.remeta.buffer.ld",
                    f"{prefix}.remeta.ld.idx.gz",
                ]
            )
    for trait in [row["trait"] for row in phase2_active_group_rows(checkpoints)]:
        paths.append(f"results/remeta/export/{build}/htp/{trait}.PAN.regenie.gz")
    return paths


def remeta_group_htp_inputs(config, group, build):
    return [
        f"results/remeta/export/{build}/htp/{trait}.PAN.regenie.gz"
        for trait in phase2_group_traits(config, group)
    ]


def remeta_group_ld_inputs(group, build):
    paths = []
    for chrom in range(1, 23):
        prefix = f"results/remeta/export/{build}/ld/{group}/chr{chrom}"
        paths.extend(
            [
                f"{prefix}.remeta.gene.ld",
                f"{prefix}.remeta.buffer.ld",
                f"{prefix}.remeta.ld.idx.gz",
            ]
        )
    return paths


def remeta_target_summaries(checkpoints, config, build):
    return [
        f"results/remeta/work/{build}/groups/{row['group']}/{row['group']}.target.summary.tsv"
        for row in phase2_active_group_rows(checkpoints)
    ]


def remeta_trait_summaries(config):
    return [
        f"results/qc/phase2_regenie/groups/{group}/{group}.trait_summary.tsv"
        for group in phase2_group_ids(config)
    ]


def remeta_validations(checkpoints, config, build):
    return [
        f"results/remeta/work/{build}/groups/{row['group']}/{row['group']}.validation.ok"
        for row in phase2_active_group_rows(checkpoints)
    ]


def remeta_validation_for_trait(checkpoints, wildcards, config):
    if not remeta_enabled(config):
        return []
    status = phase2_trait_status(checkpoints, wildcards.trait)
    if status["skipped"] == "True":
        return []
    group = phase2_trait_group(config, wildcards.trait)
    return [f"results/remeta/work/{wildcards.build}/groups/{group}/{group}.validation.ok"]


# Keep POP-MaD outputs in the production ancestry directory.
def popmad_output_dir(config):
    return "results/qc/ancestry/production"


# Build paths inside the active POP-MaD output directory.
def popmad_path(config, name):
    return f"{popmad_output_dir(config)}/{name}"


# Resolve the package-backed POP-MaD ancestry assignment file.
def ancestry_file(config):
    return popmad_path(config, "popmad_assignments.tsv")


# Resolve the projected study PC file used by POP-MaD.
def ancestry_study_pcs_file(config, _wildcards=None):
    return "results/qc/ancestry/reference/study_projected_pcs.tsv"


# Resolve the reference PC file used by POP-MaD.
def ancestry_reference_pcs_file(config, _wildcards=None):
    return "results/qc/ancestry/reference/reference_pcs.tsv"


# Resolve the PC table that includes assigned ancestry labels.
def popmad_within_file(config):
    return f"{popmad_output_dir(config)}/popmad_projected_pcs_with_ancestry.tsv"


# Resolve the package-backed within-ancestry PC file for GWAS covariates.
def pcs_file(config):
    return "results/qc/ancestry/within_ancestry_pcs.tsv"
