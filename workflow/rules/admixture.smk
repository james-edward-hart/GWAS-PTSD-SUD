ADMIXTURE_DIR = "results/qc/admixture"
ADMIXTURE_RAW_DIR = f"{ADMIXTURE_DIR}/raw"
ADMIXTURE_BY_ANCESTRY_DIR = f"{ADMIXTURE_DIR}/by_ancestry"
ADMIXTURE_STRATUM_DIR = f"{ADMIXTURE_BY_ANCESTRY_DIR}/{{ancestry}}"
ADMIXTURE_STRATUM_RAW_DIR = f"{ADMIXTURE_RAW_DIR}/{{ancestry}}"
ADMIXTURE_K = int(config.get("admixture", {}).get("k", 5))

ADMIXTURE_REFERENCE_QC_PREFIX = f"{ADMIXTURE_RAW_DIR}/reference_qc"
# TEMPORARY WORKAROUND: remove this artifact and its producing rule as soon as
# the repaired GRCh38 reference package replaces the known broken package.
ADMIXTURE_TEMPORARY_REFERENCE_REMOVE = (
    f"{ADMIXTURE_RAW_DIR}/temporary_broken_grch38_reference.remove.tsv"
)
ADMIXTURE_STUDY_QC_PREFIX = f"{ADMIXTURE_STRATUM_RAW_DIR}/study_qc"
ADMIXTURE_REFERENCE_SHARED_PREFIX = f"{ADMIXTURE_STRATUM_RAW_DIR}/reference_shared"
ADMIXTURE_STUDY_SHARED_PREFIX = f"{ADMIXTURE_STRATUM_RAW_DIR}/study_shared"
ADMIXTURE_VARIANT_MAPPING = f"{ADMIXTURE_STRATUM_RAW_DIR}/variant_harmonization.tsv"
ADMIXTURE_REFERENCE_NATIVE_VARIANTS = f"{ADMIXTURE_STRATUM_RAW_DIR}/reference_native_variants.txt"
ADMIXTURE_STUDY_NATIVE_VARIANTS = f"{ADMIXTURE_STRATUM_RAW_DIR}/study_native_variants.txt"
ADMIXTURE_REFERENCE_UPDATE_NAMES = f"{ADMIXTURE_STRATUM_RAW_DIR}/reference_update_names.tsv"
ADMIXTURE_STUDY_UPDATE_NAMES = f"{ADMIXTURE_STRATUM_RAW_DIR}/study_update_names.tsv"
ADMIXTURE_HARMONIZED_VALIDATION = f"{ADMIXTURE_STRATUM_RAW_DIR}/harmonized_variants.ok"
ADMIXTURE_REFERENCE_PRUNED_PREFIX = f"{ADMIXTURE_STRATUM_RAW_DIR}/reference_pruned"
ADMIXTURE_STUDY_PRUNED_PREFIX = f"{ADMIXTURE_STRATUM_RAW_DIR}/study_pruned"
ADMIXTURE_LD_PRUNE_PREFIX = f"{ADMIXTURE_STRATUM_RAW_DIR}/ld_prune/admixture_ld_prune"
ADMIXTURE_MERGED_PREFIX = f"{ADMIXTURE_STRATUM_RAW_DIR}/merged"

admixture_popmad_input = lambda wildcards: popmad_assignments_file()
active_admixture_outputs = lambda wildcards, filename: active_admixture_stratum_outputs(
    wildcards, filename
)


# TEMPORARY WORKAROUND — DEPRECATE AFTER THE REFERENCE PACKAGE IS REBUILT.
# The broken GRCh38 package left 629 related samples in its genotype files even
# though its metadata had already removed them. This rule creates a narrowly
# validated removal list; it must not become a general metadata-mismatch filter.
rule write_temporary_admixture_reference_remove:
    input:
        config=RUN_CONFIG,
        ok="results/qc/input_validation/validation.ok",
    output:
        exclusions=ADMIXTURE_TEMPORARY_REFERENCE_REMOVE,
    log:
        "results/logs/admixture/temporary_reference_remove.log",
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/admixture_qc.R temporary-reference-remove \
          --config {input.config} \
          --out {output.exclusions} \
          > {log} 2>&1
        """


rule convert_reference_for_admixture:
    input:
        config=RUN_CONFIG,
        ok="results/qc/input_validation/validation.ok",
        exclusions=ADMIXTURE_TEMPORARY_REFERENCE_REMOVE,
    output:
        pgen=f"{ADMIXTURE_REFERENCE_QC_PREFIX}.pgen",
        pvar=f"{ADMIXTURE_REFERENCE_QC_PREFIX}.pvar",
        psam=f"{ADMIXTURE_REFERENCE_QC_PREFIX}.psam",
    log:
        "results/logs/admixture/convert_reference.log",
    threads:
        config["runtime"]["threads_small"]
    resources:
        mem_mb=config["runtime"]["mem_mb_small"],
        runtime=config["runtime"]["time_min_small"],
    conda:
        "../../envs/gwas.yaml",
    params:
        out_prefix=lambda wildcards, output: str(output.pgen)[:-5],
    shell:
        """
        Rscript scripts/admixture_qc.R convert-reference \
          --config {input.config} \
          --remove {input.exclusions} \
          --out-prefix {params.out_prefix} \
          --threads {threads} \
          > {log} 2>&1
        """


rule convert_study_for_admixture:
    input:
        config=RUN_CONFIG,
        ok="results/qc/input_validation/validation.ok",
        keep="results/qc/strata/{ancestry}.keep.tsv",
    output:
        pgen=f"{ADMIXTURE_STUDY_QC_PREFIX}.pgen",
        pvar=f"{ADMIXTURE_STUDY_QC_PREFIX}.pvar",
        psam=f"{ADMIXTURE_STUDY_QC_PREFIX}.psam",
    log:
        "results/logs/admixture/convert_study.{ancestry}.log",
    threads:
        config["runtime"]["threads_small"]
    resources:
        mem_mb=config["runtime"]["mem_mb_small"],
        runtime=config["runtime"]["time_min_small"],
    conda:
        "../../envs/gwas.yaml",
    params:
        out_prefix=lambda wildcards, output: str(output.pgen)[:-5],
    shell:
        """
        Rscript scripts/admixture_qc.R convert-study \
          --config {input.config} \
          --keep {input.keep} \
          --out-prefix {params.out_prefix} \
          --threads {threads} \
          > {log} 2>&1
        """


rule write_admixture_shared_variants:
    input:
        config=RUN_CONFIG,
        reference=f"{ADMIXTURE_REFERENCE_QC_PREFIX}.pvar",
        study=f"{ADMIXTURE_STUDY_QC_PREFIX}.pvar",
    output:
        variants=f"{ADMIXTURE_STRATUM_RAW_DIR}/shared_variants.txt",
        mismatches=f"{ADMIXTURE_STRATUM_RAW_DIR}/shared_variant_mismatches.tsv",
        mapping=ADMIXTURE_VARIANT_MAPPING,
        reference_extract=ADMIXTURE_REFERENCE_NATIVE_VARIANTS,
        study_extract=ADMIXTURE_STUDY_NATIVE_VARIANTS,
        reference_update=ADMIXTURE_REFERENCE_UPDATE_NAMES,
        study_update=ADMIXTURE_STUDY_UPDATE_NAMES,
    log:
        "results/logs/admixture/write_shared_variants.{ancestry}.log",
    conda:
        "../../envs/gwas.yaml",
    params:
        reference_prefix=lambda wildcards, input: str(input.reference)[:-5],
        study_prefix=lambda wildcards, input: str(input.study)[:-5],
    shell:
        """
        Rscript scripts/admixture_qc.R shared-variants \
          --config {input.config} \
          --reference-prefix {params.reference_prefix} \
          --study-prefix {params.study_prefix} \
          --out {output.variants} \
          --mismatch-report {output.mismatches} \
          --mapping {output.mapping} \
          --reference-extract {output.reference_extract} \
          --study-extract {output.study_extract} \
          --reference-update {output.reference_update} \
          --study-update {output.study_update} \
          > {log} 2>&1
        """


rule extract_shared_admixture_reference:
    input:
        config=RUN_CONFIG,
        pgen=f"{ADMIXTURE_REFERENCE_QC_PREFIX}.pgen",
        pvar=f"{ADMIXTURE_REFERENCE_QC_PREFIX}.pvar",
        psam=f"{ADMIXTURE_REFERENCE_QC_PREFIX}.psam",
        variants=ADMIXTURE_REFERENCE_NATIVE_VARIANTS,
        update_names=ADMIXTURE_REFERENCE_UPDATE_NAMES,
    output:
        pgen=f"{ADMIXTURE_REFERENCE_SHARED_PREFIX}.pgen",
        pvar=f"{ADMIXTURE_REFERENCE_SHARED_PREFIX}.pvar",
        psam=f"{ADMIXTURE_REFERENCE_SHARED_PREFIX}.psam",
    log:
        "results/logs/admixture/extract_shared_reference.{ancestry}.log",
    threads:
        config["runtime"]["threads_small"]
    resources:
        mem_mb=config["runtime"]["mem_mb_small"],
        runtime=config["runtime"]["time_min_small"],
    conda:
        "../../envs/gwas.yaml",
    params:
        input_prefix=lambda wildcards, input: str(input.pgen)[:-5],
        out_prefix=lambda wildcards, output: str(output.pgen)[:-5],
    shell:
        """
        Rscript scripts/admixture_qc.R extract-variants \
          --config {input.config} \
          --input-prefix {params.input_prefix} \
          --variants {input.variants} \
          --update-name {input.update_names} \
          --out-prefix {params.out_prefix} \
          --threads {threads} \
          > {log} 2>&1
        """


rule extract_shared_admixture_study:
    input:
        config=RUN_CONFIG,
        pgen=f"{ADMIXTURE_STUDY_QC_PREFIX}.pgen",
        pvar=f"{ADMIXTURE_STUDY_QC_PREFIX}.pvar",
        psam=f"{ADMIXTURE_STUDY_QC_PREFIX}.psam",
        variants=ADMIXTURE_STUDY_NATIVE_VARIANTS,
        update_names=ADMIXTURE_STUDY_UPDATE_NAMES,
    output:
        pgen=f"{ADMIXTURE_STUDY_SHARED_PREFIX}.pgen",
        pvar=f"{ADMIXTURE_STUDY_SHARED_PREFIX}.pvar",
        psam=f"{ADMIXTURE_STUDY_SHARED_PREFIX}.psam",
    log:
        "results/logs/admixture/extract_shared_study.{ancestry}.log",
    threads:
        config["runtime"]["threads_small"]
    resources:
        mem_mb=config["runtime"]["mem_mb_small"],
        runtime=config["runtime"]["time_min_small"],
    conda:
        "../../envs/gwas.yaml",
    params:
        input_prefix=lambda wildcards, input: str(input.pgen)[:-5],
        out_prefix=lambda wildcards, output: str(output.pgen)[:-5],
    shell:
        """
        Rscript scripts/admixture_qc.R extract-variants \
          --config {input.config} \
          --input-prefix {params.input_prefix} \
          --variants {input.variants} \
          --update-name {input.update_names} \
          --out-prefix {params.out_prefix} \
          --threads {threads} \
          > {log} 2>&1
        """


rule validate_harmonized_admixture_variants:
    input:
        config=RUN_CONFIG,
        reference_pgen=f"{ADMIXTURE_REFERENCE_SHARED_PREFIX}.pgen",
        reference_pvar=f"{ADMIXTURE_REFERENCE_SHARED_PREFIX}.pvar",
        reference_psam=f"{ADMIXTURE_REFERENCE_SHARED_PREFIX}.psam",
        study_pgen=f"{ADMIXTURE_STUDY_SHARED_PREFIX}.pgen",
        study_pvar=f"{ADMIXTURE_STUDY_SHARED_PREFIX}.pvar",
        study_psam=f"{ADMIXTURE_STUDY_SHARED_PREFIX}.psam",
    output:
        ok=ADMIXTURE_HARMONIZED_VALIDATION,
    log:
        "results/logs/admixture/validate_harmonized_variants.{ancestry}.log",
    conda:
        "../../envs/gwas.yaml",
    params:
        reference_prefix=lambda wildcards, input: str(input.reference_pgen)[:-5],
        study_prefix=lambda wildcards, input: str(input.study_pgen)[:-5],
    shell:
        """
        Rscript scripts/admixture_qc.R validate-harmonized \
          --config {input.config} \
          --reference-prefix {params.reference_prefix} \
          --study-prefix {params.study_prefix} \
          --out {output.ok} \
          > {log} 2>&1
        """


rule ld_prune_admixture_markers:
    input:
        config=RUN_CONFIG,
        pgen=f"{ADMIXTURE_REFERENCE_SHARED_PREFIX}.pgen",
        pvar=f"{ADMIXTURE_REFERENCE_SHARED_PREFIX}.pvar",
        psam=f"{ADMIXTURE_REFERENCE_SHARED_PREFIX}.psam",
        harmonized=ADMIXTURE_HARMONIZED_VALIDATION,
    output:
        prune_in=f"{ADMIXTURE_LD_PRUNE_PREFIX}.prune.in",
        prune_out=f"{ADMIXTURE_LD_PRUNE_PREFIX}.prune.out",
        excluded=f"{ADMIXTURE_LD_PRUNE_PREFIX}.excluded_region_variants.txt",
    log:
        "results/logs/admixture/ld_prune.{ancestry}.log",
    threads:
        config["runtime"]["threads_small"]
    resources:
        mem_mb=config["runtime"]["mem_mb_small"],
        runtime=config["runtime"]["time_min_small"],
    conda:
        "../../envs/gwas.yaml",
    params:
        pfile_prefix=lambda wildcards, input: str(input.pgen)[:-5],
        out_prefix=lambda wildcards, output: str(output.prune_in).replace(".prune.in", ""),
    shell:
        """
        Rscript scripts/admixture_qc.R ld-prune \
          --config {input.config} \
          --pfile-prefix {params.pfile_prefix} \
          --out-prefix {params.out_prefix} \
          --prune-in {output.prune_in} \
          --excluded-regions {output.excluded} \
          --threads {threads} \
          > {log} 2>&1
        """


rule extract_pruned_admixture_reference:
    input:
        config=RUN_CONFIG,
        pgen=f"{ADMIXTURE_REFERENCE_SHARED_PREFIX}.pgen",
        pvar=f"{ADMIXTURE_REFERENCE_SHARED_PREFIX}.pvar",
        psam=f"{ADMIXTURE_REFERENCE_SHARED_PREFIX}.psam",
        variants=f"{ADMIXTURE_LD_PRUNE_PREFIX}.prune.in",
    output:
        pgen=f"{ADMIXTURE_REFERENCE_PRUNED_PREFIX}.pgen",
        pvar=f"{ADMIXTURE_REFERENCE_PRUNED_PREFIX}.pvar",
        psam=f"{ADMIXTURE_REFERENCE_PRUNED_PREFIX}.psam",
    log:
        "results/logs/admixture/extract_pruned_reference.{ancestry}.log",
    threads:
        config["runtime"]["threads_small"]
    resources:
        mem_mb=config["runtime"]["mem_mb_small"],
        runtime=config["runtime"]["time_min_small"],
    conda:
        "../../envs/gwas.yaml",
    params:
        input_prefix=lambda wildcards, input: str(input.pgen)[:-5],
        out_prefix=lambda wildcards, output: str(output.pgen)[:-5],
    shell:
        """
        Rscript scripts/admixture_qc.R extract-variants \
          --config {input.config} \
          --input-prefix {params.input_prefix} \
          --variants {input.variants} \
          --out-prefix {params.out_prefix} \
          --threads {threads} \
          > {log} 2>&1
        """


rule extract_pruned_admixture_study:
    input:
        config=RUN_CONFIG,
        pgen=f"{ADMIXTURE_STUDY_SHARED_PREFIX}.pgen",
        pvar=f"{ADMIXTURE_STUDY_SHARED_PREFIX}.pvar",
        psam=f"{ADMIXTURE_STUDY_SHARED_PREFIX}.psam",
        variants=f"{ADMIXTURE_LD_PRUNE_PREFIX}.prune.in",
    output:
        pgen=f"{ADMIXTURE_STUDY_PRUNED_PREFIX}.pgen",
        pvar=f"{ADMIXTURE_STUDY_PRUNED_PREFIX}.pvar",
        psam=f"{ADMIXTURE_STUDY_PRUNED_PREFIX}.psam",
    log:
        "results/logs/admixture/extract_pruned_study.{ancestry}.log",
    threads:
        config["runtime"]["threads_small"]
    resources:
        mem_mb=config["runtime"]["mem_mb_small"],
        runtime=config["runtime"]["time_min_small"],
    conda:
        "../../envs/gwas.yaml",
    params:
        input_prefix=lambda wildcards, input: str(input.pgen)[:-5],
        out_prefix=lambda wildcards, output: str(output.pgen)[:-5],
    shell:
        """
        Rscript scripts/admixture_qc.R extract-variants \
          --config {input.config} \
          --input-prefix {params.input_prefix} \
          --variants {input.variants} \
          --out-prefix {params.out_prefix} \
          --threads {threads} \
          > {log} 2>&1
        """


rule merge_admixture_genotypes:
    input:
        config=RUN_CONFIG,
        reference_pgen=f"{ADMIXTURE_REFERENCE_PRUNED_PREFIX}.pgen",
        reference_pvar=f"{ADMIXTURE_REFERENCE_PRUNED_PREFIX}.pvar",
        reference_psam=f"{ADMIXTURE_REFERENCE_PRUNED_PREFIX}.psam",
        study_pgen=f"{ADMIXTURE_STUDY_PRUNED_PREFIX}.pgen",
        study_pvar=f"{ADMIXTURE_STUDY_PRUNED_PREFIX}.pvar",
        study_psam=f"{ADMIXTURE_STUDY_PRUNED_PREFIX}.psam",
    output:
        merged_pgen=f"{ADMIXTURE_MERGED_PREFIX}_pmerge.pgen",
        merged_pvar=f"{ADMIXTURE_MERGED_PREFIX}_pmerge.pvar",
        merged_psam=f"{ADMIXTURE_MERGED_PREFIX}_pmerge.psam",
        bed=f"{ADMIXTURE_MERGED_PREFIX}.bed",
        bim=f"{ADMIXTURE_MERGED_PREFIX}.bim",
        fam=f"{ADMIXTURE_MERGED_PREFIX}.fam",
    log:
        "results/logs/admixture/merge_genotypes.{ancestry}.log",
    threads:
        config["runtime"]["threads_small"]
    resources:
        mem_mb=config["runtime"]["mem_mb_small"],
        runtime=config["runtime"]["time_min_small"],
    conda:
        "../../envs/gwas.yaml",
    params:
        reference_prefix=lambda wildcards, input: str(input.reference_pgen)[:-5],
        study_prefix=lambda wildcards, input: str(input.study_pgen)[:-5],
        out_prefix=lambda wildcards, output: str(output.bed)[:-4],
    shell:
        """
        Rscript scripts/admixture_qc.R merge \
          --config {input.config} \
          --reference-prefix {params.reference_prefix} \
          --study-prefix {params.study_prefix} \
          --out-prefix {params.out_prefix} \
          --threads {threads} \
          > {log} 2>&1
        """


rule write_admixture_pop_file:
    input:
        config=RUN_CONFIG,
        fam=f"{ADMIXTURE_MERGED_PREFIX}.fam",
        reference_psam=f"{ADMIXTURE_REFERENCE_PRUNED_PREFIX}.psam",
        study_psam=f"{ADMIXTURE_STUDY_PRUNED_PREFIX}.psam",
    output:
        pop_file=f"{ADMIXTURE_MERGED_PREFIX}.pop",
        samples=f"{ADMIXTURE_STRATUM_RAW_DIR}/sample_populations.tsv",
    log:
        "results/logs/admixture/write_pop.{ancestry}.log",
    conda:
        "../../envs/gwas.yaml",
    params:
        reference_prefix=lambda wildcards, input: str(input.reference_psam)[:-5],
        study_prefix=lambda wildcards, input: str(input.study_psam)[:-5],
    shell:
        """
        Rscript scripts/admixture_qc.R write-pop \
          --config {input.config} \
          --fam {input.fam} \
          --reference-prefix {params.reference_prefix} \
          --study-prefix {params.study_prefix} \
          --pop-out {output.pop_file} \
          --sample-populations {output.samples} \
          > {log} 2>&1
        """


rule run_supervised_admixture:
    input:
        config=RUN_CONFIG,
        bed=f"{ADMIXTURE_MERGED_PREFIX}.bed",
        bim=f"{ADMIXTURE_MERGED_PREFIX}.bim",
        fam=f"{ADMIXTURE_MERGED_PREFIX}.fam",
        pop_file=f"{ADMIXTURE_MERGED_PREFIX}.pop",
    output:
        q=f"{ADMIXTURE_MERGED_PREFIX}.{ADMIXTURE_K}.Q",
        p=f"{ADMIXTURE_MERGED_PREFIX}.{ADMIXTURE_K}.P",
    log:
        "results/logs/admixture/run_supervised_admixture.{ancestry}.log",
    threads:
        config["runtime"].get("threads_admixture", config["runtime"]["threads_small"])
    resources:
        mem_mb=config["runtime"].get("mem_mb_admixture", config["runtime"]["mem_mb_small"]),
        runtime=config["runtime"].get("time_min_admixture", config["runtime"]["time_min_small"]),
    conda:
        "../../envs/gwas.yaml",
    params:
        k=ADMIXTURE_K,
    shell:
        """
        Rscript scripts/admixture_qc.R run-admixture \
          --config {input.config} \
          --bed {input.bed} \
          --k {params.k} \
          --threads {threads} \
          > {log} 2>&1
        """


rule write_admixture_stratum_report:
    input:
        config=RUN_CONFIG,
        q=f"{ADMIXTURE_MERGED_PREFIX}.{ADMIXTURE_K}.Q",
        p=f"{ADMIXTURE_MERGED_PREFIX}.{ADMIXTURE_K}.P",
        fam=f"{ADMIXTURE_MERGED_PREFIX}.fam",
        bim=f"{ADMIXTURE_MERGED_PREFIX}.bim",
        pop_file=f"{ADMIXTURE_MERGED_PREFIX}.pop",
        sample_populations=f"{ADMIXTURE_STRATUM_RAW_DIR}/sample_populations.tsv",
        popmad=admixture_popmad_input,
    output:
        study=f"{ADMIXTURE_STRATUM_DIR}/study_ancestry_proportions.tsv",
        reference=f"{ADMIXTURE_STRATUM_DIR}/reference_ancestry_proportions.tsv",
        comparison=f"{ADMIXTURE_STRATUM_DIR}/popmad_admixture_comparison.tsv",
        summary=f"{ADMIXTURE_STRATUM_DIR}/admixture_run_summary.tsv",
        report=f"{ADMIXTURE_STRATUM_DIR}/admixture_report.md",
    log:
        "results/logs/admixture/write_report.{ancestry}.log",
    conda:
        "../../envs/gwas.yaml",
    params:
        popmad_arg=lambda wildcards, input: f"--popmad {input.popmad}" if input.popmad else "",
    shell:
        """
        Rscript scripts/admixture_qc.R parse-report \
          --config {input.config} \
          --analysis-ancestry {wildcards.ancestry} \
          --q {input.q} \
          --p {input.p} \
          --fam {input.fam} \
          --bim {input.bim} \
          --pop {input.pop_file} \
          --sample-populations {input.sample_populations} \
          {params.popmad_arg} \
          --study-out {output.study} \
          --reference-out {output.reference} \
          --comparison-out {output.comparison} \
          --summary-out {output.summary} \
          --report-out {output.report} \
          > {log} 2>&1
        """


rule write_admixture_report:
    input:
        config=RUN_CONFIG,
        studies=lambda wildcards: active_admixture_outputs(wildcards, "study_ancestry_proportions.tsv"),
        references=lambda wildcards: active_admixture_outputs(wildcards, "reference_ancestry_proportions.tsv"),
        comparisons=lambda wildcards: active_admixture_outputs(wildcards, "popmad_admixture_comparison.tsv"),
        summaries=lambda wildcards: active_admixture_outputs(wildcards, "admixture_run_summary.tsv"),
        reports=lambda wildcards: active_admixture_outputs(wildcards, "admixture_report.md"),
    output:
        study=f"{ADMIXTURE_DIR}/study_ancestry_proportions.tsv",
        reference=f"{ADMIXTURE_DIR}/reference_ancestry_proportions.tsv",
        comparison=f"{ADMIXTURE_DIR}/popmad_admixture_comparison.tsv",
        summary=f"{ADMIXTURE_DIR}/admixture_run_summary.tsv",
        report=f"{ADMIXTURE_DIR}/admixture_report.md",
    log:
        "results/logs/admixture/write_report.log",
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/admixture_qc.R combine-reports \
          --config {input.config} \
          --study {input.studies} \
          --reference {input.references} \
          --comparison {input.comparisons} \
          --summary {input.summaries} \
          --report {input.reports} \
          --study-out {output.study} \
          --reference-out {output.reference} \
          --comparison-out {output.comparison} \
          --summary-out {output.summary} \
          --report-out {output.report} \
          > {log} 2>&1
        """
