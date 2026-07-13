# Snakemake reads the editable user config here; validation.smk later writes
# immutable run snapshots after config merging and reference-package resolution.
configfile: "config/config.yaml"

from glob import glob
from workflow.snake_helpers import (
    active_admixture_stratum_outputs as active_admixture_stratum_outputs_for_context,
    active_unrelated_keep_files as active_unrelated_keep_files_for_context,
    active_within_ancestry_eigenvecs as active_within_ancestry_eigenvecs_for_context,
    analysis_output_name as analysis_output_name_for_config,
    admixture_targets as admixture_targets_for_config,
    admixture_validation_inputs as admixture_validation_inputs_for_config,
    ancestry_file as ancestry_file_for_config,
    ancestry_reference_pcs_file as ancestry_reference_pcs_file_for_config,
    ancestry_reference_targets as ancestry_reference_targets_for_config,
    ancestry_reference_validation_inputs as ancestry_reference_validation_inputs_for_config,
    ancestry_study_pcs_file as ancestry_study_pcs_file_for_config,
    pcs_file as pcs_file_for_config,
    phase2_group_stage1_stats as phase2_group_stage1_stats_for_context,
    phase2_report_targets as phase2_report_targets_for_context,
    phase2_trait_group as phase2_trait_group_for_config,
    phase2_trait_regenie_done as phase2_trait_regenie_done_for_context,
    phase2_trait_stage1_summaries as phase2_trait_stage1_summaries_for_context,
    remeta_manifest_inputs as remeta_manifest_inputs_for_config,
    remeta_group_htp_inputs as remeta_group_htp_inputs_for_config,
    remeta_group_index_inputs as remeta_group_index_inputs_for_config,
    remeta_resource_file as remeta_resource_file_for_config,
    remeta_target_summaries as remeta_target_summaries_for_config,
    remeta_targets as remeta_targets_for_context,
    remeta_trait_summaries as remeta_trait_summaries_for_config,
    remeta_validations as remeta_validations_for_config,
    popmad_path,
    popmad_within_file as popmad_within_file_for_config,
    reference_package_inputs as reference_package_inputs_for_config,
    report_targets as report_targets_for_context,
    trait_ids as trait_ids_for_config,
)


# EFFECTIVE_CONFIG captures exactly what Snakemake saw at parse time.
# RUN_CONFIG is the build-resolved version consumed by all downstream scripts.
EFFECTIVE_CONFIG = "results/config/effective_config.yaml"
RUN_CONFIG = "results/config/resolved_config.yaml"

wildcard_constraints:
    build="GRCh[0-9]+",

CONFIG_SOURCE_FILES = sorted(str(path) for path in workflow.configfiles)
if not CONFIG_SOURCE_FILES:
    CONFIG_SOURCE_FILES = ["config/config.yaml"]

# Track manifest-listed package files as inputs so package edits invalidate the
# resolved config instead of silently reusing an older reference mapping.
REFERENCE_PACKAGE_INPUTS = reference_package_inputs_for_config(config)

PIPELINE_CODE = sorted(
    [
        "Snakefile",
        "workflow/snake_helpers.py",
        *glob("workflow/rules/*.smk"),
        *glob("scripts/*.R"),
        *glob("scripts/*.sh"),
        *glob("scripts/*.py"),
        *glob("scripts/lib/*.R"),
        *glob("scripts/lib/*.sh"),
        *glob("envs/*.yaml"),
    ]
)

TRAITS = trait_ids_for_config(config)
ANCESTRIES = config["analysis"]["ancestries"]
ANCESTRY_REFERENCE_ENABLED = bool(config.get("ancestry_reference", {}).get("enabled", False))
ANALYSIS_OUTPUT_NAME = analysis_output_name_for_config(config)

# These lambdas delay parts of target expansion until checkpoint outputs exist.
ancestry_reference_targets = lambda: ancestry_reference_targets_for_config(config)
ancestry_reference_validation_inputs = lambda wildcards: ancestry_reference_validation_inputs_for_config(config, wildcards)
admixture_targets = lambda: admixture_targets_for_config(config)
admixture_validation_inputs = lambda wildcards: admixture_validation_inputs_for_config(config, wildcards)
report_targets = lambda wildcards: report_targets_for_context(checkpoints, TRAITS, wildcards, config)
phase2_report_targets = lambda wildcards: phase2_report_targets_for_context(checkpoints, TRAITS, wildcards, config)
phase2_group_stage1_stats = lambda wildcards: phase2_group_stage1_stats_for_context(checkpoints, wildcards, config)
phase2_trait_stage1_summaries = lambda wildcards: phase2_trait_stage1_summaries_for_context(checkpoints, wildcards, config)
phase2_trait_regenie_done = lambda wildcards: phase2_trait_regenie_done_for_context(wildcards, config)
phase2_trait_group = lambda wildcards: phase2_trait_group_for_config(config, wildcards.trait)
remeta_targets = lambda wildcards: remeta_targets_for_context(checkpoints, wildcards, config)
remeta_manifest_inputs = lambda wildcards: remeta_manifest_inputs_for_config(config, wildcards.build)
remeta_group_htp_inputs = lambda wildcards: remeta_group_htp_inputs_for_config(config, wildcards.group, wildcards.build)
remeta_group_index_inputs = lambda wildcards: remeta_group_index_inputs_for_config(wildcards.group, wildcards.build)
remeta_target_summaries = lambda wildcards: remeta_target_summaries_for_config(config, wildcards.build)
remeta_trait_summaries = lambda wildcards: remeta_trait_summaries_for_config(config)
remeta_validations = lambda wildcards: remeta_validations_for_config(config, wildcards.build)
remeta_resource_file = lambda build, filename: remeta_resource_file_for_config(config, build, filename)
active_admixture_stratum_outputs = lambda wildcards, filename: active_admixture_stratum_outputs_for_context(checkpoints, wildcards, filename)
active_unrelated_keep_files = lambda wildcards: active_unrelated_keep_files_for_context(checkpoints, wildcards)
active_within_ancestry_eigenvecs = lambda wildcards: active_within_ancestry_eigenvecs_for_context(checkpoints, wildcards)
ancestry_file = lambda: ancestry_file_for_config(config)
popmad_assignments_file = lambda: popmad_path(config, "popmad_assignments.tsv")
popmad_study_pcs_file = lambda: popmad_path(config, "study_pcs.tsv")
popmad_distances_file = lambda: popmad_path(config, "popmad_distances.tsv")
popmad_reference_outliers_file = lambda: popmad_path(config, "reference_outliers.tsv")
popmad_model_summary_file = lambda: popmad_path(config, "population_model_summary.tsv")
popmad_excluded_file = lambda: popmad_path(config, "popmad_excluded.tsv")
popmad_counts_file = lambda: popmad_path(config, "popmad_population_counts.tsv")
ancestry_study_pcs_file = lambda wildcards=None: ancestry_study_pcs_file_for_config(config, wildcards)
ancestry_reference_pcs_file = lambda wildcards=None: ancestry_reference_pcs_file_for_config(config, wildcards)
popmad_within_file = lambda: popmad_within_file_for_config(config)
pcs_file = lambda: pcs_file_for_config(config)


rule all:
    input:
        "results/qc/input_validation/validation.ok",
        # Report filenames include the inferred genome build, so they are
        # expanded after the genome-build checkpoint completes.
        report_targets,
        phase2_report_targets,
        remeta_targets,
        ancestry_reference_targets(),
        admixture_targets(),
        "results/manifests/run_manifest.tsv",


include: "workflow/rules/validation.smk"
include: "workflow/rules/ancestry_reference.smk"
include: "workflow/rules/ancestry.smk"
include: "workflow/rules/admixture.smk"
include: "workflow/rules/sample_prep.smk"
include: "workflow/rules/gwas.smk"
include: "workflow/rules/phase2_regenie.smk"
include: "workflow/rules/remeta.smk"
include: "workflow/rules/reporting.smk"
