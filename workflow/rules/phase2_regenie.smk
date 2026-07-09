PHASE2_DIR = "results/qc/phase2_regenie"
PHASE2_PAN_PREFIX = f"{PHASE2_DIR}/pan_sample_qc/pan"
PHASE2_GLOBAL_PCA_QC_PREFIX = f"{PHASE2_DIR}/global_pca/global_pca_qc"
PHASE2_GLOBAL_PCA_PRUNE_PREFIX = f"{PHASE2_DIR}/global_pca/global_pca_ld_prune"
PHASE2_GLOBAL_PCA_PREFIX = f"{PHASE2_DIR}/global_pca/global_pca"
PHASE2_GLOBAL_PCA_SCORE_PREFIX = f"{PHASE2_DIR}/global_pca/pan_global_pcs"
PHASE2_STEP1_QC_PREFIX = f"{PHASE2_DIR}/step1/step1_qc"
PHASE2_STEP1_PRUNE_PREFIX = f"{PHASE2_DIR}/step1/step1_ld_prune"
PHASE2_REGENIE_TOOL = f"{PHASE2_DIR}/regenie_tool.tsv"


rule record_phase2_regenie_tool:
    input:
        config=RUN_CONFIG,
    output:
        tool=PHASE2_REGENIE_TOOL,
    log:
        "results/logs/phase2_regenie/record_regenie_tool.log",
    conda:
        "../../envs/regenie.yaml",
    params:
        tool=lambda wildcards: config.get("tools", {}).get("regenie", "regenie"),
    shell:
        """
        python scripts/record_regenie_tool.py \
          --tool {params.tool:q} \
          --out {output.tool:q} \
          > {log:q} 2>&1
        """


rule write_phase2_regenie_groups:
    input:
        config=RUN_CONFIG,
        ok="results/qc/input_validation/validation.ok",
    output:
        groups=f"{PHASE2_DIR}/groups/trait_groups.tsv",
    log:
        "results/logs/phase2_regenie/write_groups.log",
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/phase2_regenie.R write-groups \
          --config {input.config} \
          --out {output.groups} \
          > {log} 2>&1
        """


rule prepare_phase2_pan_genotypes:
    input:
        config=RUN_CONFIG,
        ok="results/qc/input_validation/validation.ok",
        sex_keep="results/qc/sex/sex_checked.keep.tsv",
        assignments=popmad_assignments_file(),
        excluded=popmad_excluded_file(),
    output:
        pgen=f"{PHASE2_PAN_PREFIX}.pgen",
        pvar=f"{PHASE2_PAN_PREFIX}.pvar",
        psam=f"{PHASE2_PAN_PREFIX}.psam",
        keep=f"{PHASE2_DIR}/pan.keep.tsv",
        ancestry=f"{PHASE2_DIR}/pan_ancestry_report.tsv",
        summary=f"{PHASE2_DIR}/pan_sample_summary.tsv",
    log:
        "results/logs/phase2_regenie/prepare_pan_genotypes.log",
    threads:
        config["runtime"].get("threads_phase2_pca", config["runtime"]["threads_gwas"])
    resources:
        mem_mb=config["runtime"].get("mem_mb_phase2_pca", config["runtime"]["mem_mb_gwas"]),
        runtime=config["runtime"].get("time_min_phase2_pca", config["runtime"]["time_min_gwas"]),
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/phase2_regenie.R prepare-pan-genotypes \
          --config {input.config} \
          --sex-keep {input.sex_keep} \
          --assignments {input.assignments} \
          --excluded {input.excluded} \
          --out-prefix {PHASE2_PAN_PREFIX} \
          --keep-out {output.keep} \
          --ancestry-out {output.ancestry} \
          --summary-out {output.summary} \
          --threads {threads} \
          > {log} 2>&1
        """


rule prepare_phase2_global_pca_markers:
    input:
        config=RUN_CONFIG,
        pgen=f"{PHASE2_PAN_PREFIX}.pgen",
        pvar=f"{PHASE2_PAN_PREFIX}.pvar",
        psam=f"{PHASE2_PAN_PREFIX}.psam",
        unrelated="results/qc/relatedness/unrelated.king.cutoff.in.id",
    output:
        pgen=f"{PHASE2_GLOBAL_PCA_QC_PREFIX}.pgen",
        pvar=f"{PHASE2_GLOBAL_PCA_QC_PREFIX}.pvar",
        psam=f"{PHASE2_GLOBAL_PCA_QC_PREFIX}.psam",
        prune_in=f"{PHASE2_GLOBAL_PCA_PRUNE_PREFIX}.prune.in",
        prune_out=f"{PHASE2_GLOBAL_PCA_PRUNE_PREFIX}.prune.out",
        excluded=f"{PHASE2_GLOBAL_PCA_PRUNE_PREFIX}.excluded_region_variants.txt",
    log:
        "results/logs/phase2_regenie/global_pca_markers.log",
    threads:
        config["runtime"].get("threads_phase2_pca", config["runtime"]["threads_gwas"])
    resources:
        mem_mb=config["runtime"].get("mem_mb_phase2_pca", config["runtime"]["mem_mb_gwas"]),
        runtime=config["runtime"].get("time_min_phase2_pca", config["runtime"]["time_min_gwas"]),
    conda:
        "../../envs/gwas.yaml",
    params:
        pfile_prefix=lambda wildcards, input: str(input.pgen)[:-5],
    shell:
        """
        Rscript scripts/phase2_regenie.R prepare-marker-set \
          --config {input.config} \
          --branch global_pca \
          --pfile-prefix {params.pfile_prefix} \
          --keep {input.unrelated} \
          --out-prefix {PHASE2_GLOBAL_PCA_QC_PREFIX} \
          --prune-prefix {PHASE2_GLOBAL_PCA_PRUNE_PREFIX} \
          --prune-in {output.prune_in} \
          --excluded-regions {output.excluded} \
          --threads {threads} \
          > {log} 2>&1
        """


rule prepare_phase2_step1_markers:
    input:
        config=RUN_CONFIG,
        pgen=f"{PHASE2_PAN_PREFIX}.pgen",
        pvar=f"{PHASE2_PAN_PREFIX}.pvar",
        psam=f"{PHASE2_PAN_PREFIX}.psam",
    output:
        pgen=f"{PHASE2_STEP1_QC_PREFIX}.pgen",
        pvar=f"{PHASE2_STEP1_QC_PREFIX}.pvar",
        psam=f"{PHASE2_STEP1_QC_PREFIX}.psam",
        prune_in=f"{PHASE2_STEP1_PRUNE_PREFIX}.prune.in",
        prune_out=f"{PHASE2_STEP1_PRUNE_PREFIX}.prune.out",
        excluded=f"{PHASE2_STEP1_PRUNE_PREFIX}.excluded_region_variants.txt",
    log:
        "results/logs/phase2_regenie/step1_markers.log",
    threads:
        config["runtime"].get("threads_regenie_step1", config["runtime"]["threads_gwas"])
    resources:
        mem_mb=config["runtime"].get("mem_mb_regenie_step1", config["runtime"]["mem_mb_gwas"]),
        runtime=config["runtime"].get("time_min_regenie_step1", config["runtime"]["time_min_gwas"]),
    conda:
        "../../envs/gwas.yaml",
    params:
        pfile_prefix=lambda wildcards, input: str(input.pgen)[:-5],
    shell:
        """
        Rscript scripts/phase2_regenie.R prepare-marker-set \
          --config {input.config} \
          --branch step1 \
          --pfile-prefix {params.pfile_prefix} \
          --out-prefix {PHASE2_STEP1_QC_PREFIX} \
          --prune-prefix {PHASE2_STEP1_PRUNE_PREFIX} \
          --prune-in {output.prune_in} \
          --excluded-regions {output.excluded} \
          --threads {threads} \
          > {log} 2>&1
        """


rule fit_phase2_global_pca:
    input:
        config=RUN_CONFIG,
        pgen=f"{PHASE2_GLOBAL_PCA_QC_PREFIX}.pgen",
        pvar=f"{PHASE2_GLOBAL_PCA_QC_PREFIX}.pvar",
        psam=f"{PHASE2_GLOBAL_PCA_QC_PREFIX}.psam",
        variants=f"{PHASE2_GLOBAL_PCA_PRUNE_PREFIX}.prune.in",
    output:
        eigenvec=f"{PHASE2_GLOBAL_PCA_PREFIX}.eigenvec",
        eigenval=f"{PHASE2_GLOBAL_PCA_PREFIX}.eigenval",
        allele=f"{PHASE2_GLOBAL_PCA_PREFIX}.eigenvec.allele",
        acount=f"{PHASE2_GLOBAL_PCA_PREFIX}.acount",
    log:
        "results/logs/phase2_regenie/fit_global_pca.log",
    threads:
        config["runtime"].get("threads_phase2_pca", config["runtime"]["threads_gwas"])
    resources:
        mem_mb=config["runtime"].get("mem_mb_phase2_pca", config["runtime"]["mem_mb_gwas"]),
        runtime=config["runtime"].get("time_min_phase2_pca", config["runtime"]["time_min_gwas"]),
    conda:
        "../../envs/gwas.yaml",
    params:
        pfile_prefix=lambda wildcards, input: str(input.pgen)[:-5],
    shell:
        """
        Rscript scripts/phase2_regenie.R fit-global-pca \
          --config {input.config} \
          --pfile-prefix {params.pfile_prefix} \
          --variants {input.variants} \
          --out-prefix {PHASE2_GLOBAL_PCA_PREFIX} \
          --threads {threads} \
          > {log} 2>&1
        """


rule score_phase2_global_pcs:
    input:
        config=RUN_CONFIG,
        pgen=f"{PHASE2_PAN_PREFIX}.pgen",
        pvar=f"{PHASE2_PAN_PREFIX}.pvar",
        psam=f"{PHASE2_PAN_PREFIX}.psam",
        variants=f"{PHASE2_GLOBAL_PCA_PRUNE_PREFIX}.prune.in",
        weights=f"{PHASE2_GLOBAL_PCA_PREFIX}.eigenvec.allele",
        frequencies=f"{PHASE2_GLOBAL_PCA_PREFIX}.acount",
    output:
        sscore=f"{PHASE2_GLOBAL_PCA_SCORE_PREFIX}.sscore",
    log:
        "results/logs/phase2_regenie/score_global_pcs.log",
    threads:
        config["runtime"].get("threads_phase2_pca", config["runtime"]["threads_gwas"])
    resources:
        mem_mb=config["runtime"].get("mem_mb_phase2_pca", config["runtime"]["mem_mb_gwas"]),
        runtime=config["runtime"].get("time_min_phase2_pca", config["runtime"]["time_min_gwas"]),
    conda:
        "../../envs/gwas.yaml",
    params:
        pfile_prefix=lambda wildcards, input: str(input.pgen)[:-5],
    shell:
        """
        Rscript scripts/phase2_regenie.R score-global-pcs \
          --config {input.config} \
          --pfile-prefix {params.pfile_prefix} \
          --variants {input.variants} \
          --weights {input.weights} \
          --frequencies {input.frequencies} \
          --out-prefix {PHASE2_GLOBAL_PCA_SCORE_PREFIX} \
          --threads {threads} \
          > {log} 2>&1
        """


rule write_phase2_global_pcs:
    input:
        config=RUN_CONFIG,
        sscore=f"{PHASE2_GLOBAL_PCA_SCORE_PREFIX}.sscore",
    output:
        pcs=f"{PHASE2_DIR}/global_pcs.tsv",
    log:
        "results/logs/phase2_regenie/write_global_pcs.log",
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/phase2_regenie.R write-global-pcs \
          --config {input.config} \
          --sscore {input.sscore} \
          --out {output.pcs} \
          > {log} 2>&1
        """


rule build_phase2_regenie_group_inputs:
    input:
        config=RUN_CONFIG,
        groups=f"{PHASE2_DIR}/groups/trait_groups.tsv",
        keep=f"{PHASE2_DIR}/pan.keep.tsv",
        pcs=f"{PHASE2_DIR}/global_pcs.tsv",
    output:
        pheno=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.pheno.tsv",
        covar=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.covar.tsv",
        summary=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.trait_summary.tsv",
        traits=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.analysis_traits.txt",
        covars=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.covariates.txt",
        keep=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.keep.txt",
    log:
        "results/logs/phase2_regenie/build_group_inputs.{group}.log",
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/phase2_regenie.R build-group-inputs \
          --config {input.config} \
          --group {wildcards.group} \
          --keep {input.keep} \
          --pcs {input.pcs} \
          --pheno-out {output.pheno} \
          --covar-out {output.covar} \
          --summary-out {output.summary} \
          --trait-list-out {output.traits} \
          --covar-list-out {output.covars} \
          --keep-plink-out {output.keep} \
          > {log} 2>&1
        """


rule build_phase2_stage1_union:
    input:
        config=RUN_CONFIG,
        stats=phase2_group_stage1_stats,
    output:
        extract=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.stage1_union.snplist",
        summary=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.stage1_union_summary.tsv",
    log:
        "results/logs/phase2_regenie/stage1_union.{group}.log",
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/phase2_regenie.R stage1-pass-union \
          --config {input.config} \
          --stage1-stats {input.stats} \
          --out {output.extract} \
          --summary-out {output.summary} \
          > {log} 2>&1
        """


rule prepare_phase2_regenie_assoc_genotypes:
    input:
        config=RUN_CONFIG,
        pgen=f"{PHASE2_PAN_PREFIX}.pgen",
        pvar=f"{PHASE2_PAN_PREFIX}.pvar",
        psam=f"{PHASE2_PAN_PREFIX}.psam",
        extract=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.stage1_union.snplist",
    output:
        pgen=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.assoc_qc.pgen",
        pvar=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.assoc_qc.pvar",
        psam=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.assoc_qc.psam",
    log:
        "results/logs/phase2_regenie/assoc_genotypes.{group}.log",
    threads:
        config["runtime"].get("threads_regenie_step2", config["runtime"]["threads_gwas"])
    resources:
        mem_mb=config["runtime"].get("mem_mb_regenie_step2", config["runtime"]["mem_mb_gwas"]),
        runtime=config["runtime"].get("time_min_regenie_step2", config["runtime"]["time_min_gwas"]),
    conda:
        "../../envs/gwas.yaml",
    params:
        pfile_prefix=lambda wildcards, input: str(input.pgen)[:-5],
        out_prefix=lambda wildcards, output: str(output.pgen)[:-5],
    shell:
        """
        Rscript scripts/phase2_regenie.R prepare-assoc-genotypes \
          --config {input.config} \
          --pfile-prefix {params.pfile_prefix} \
          --extract {input.extract} \
          --out-prefix {params.out_prefix} \
          --threads {threads} \
          > {log} 2>&1
        """


rule filter_phase2_regenie_step1_variants:
    input:
        config=RUN_CONFIG,
        pgen=f"{PHASE2_STEP1_QC_PREFIX}.pgen",
        pvar=f"{PHASE2_STEP1_QC_PREFIX}.pvar",
        psam=f"{PHASE2_STEP1_QC_PREFIX}.psam",
        variants=f"{PHASE2_STEP1_PRUNE_PREFIX}.prune.in",
        keep=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.keep.txt",
        traits=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.analysis_traits.txt",
    output:
        variants=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.step1.snplist",
        summary=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.step1.variant_qc.summary.tsv",
        excluded=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.step1.variant_qc.excluded.tsv",
    log:
        "results/logs/phase2_regenie/filter_step1_variants.{group}.log",
    threads:
        config["runtime"].get("threads_regenie_step1", config["runtime"]["threads_gwas"])
    resources:
        mem_mb=config["runtime"].get("mem_mb_regenie_step1", config["runtime"]["mem_mb_gwas"]),
        runtime=config["runtime"].get("time_min_regenie_step1", config["runtime"]["time_min_gwas"]),
    conda:
        "../../envs/gwas.yaml",
    params:
        pfile_prefix=lambda wildcards, input: str(input.pgen)[:-5],
    shell:
        """
        Rscript scripts/phase2_regenie.R filter-step1-variants \
          --config {input.config} \
          --pfile-prefix {params.pfile_prefix} \
          --extract {input.variants} \
          --keep {input.keep} \
          --trait-list {input.traits} \
          --out {output.variants} \
          --summary-out {output.summary} \
          --excluded-out {output.excluded} \
          --threads {threads} \
          > {log} 2>&1
        """


rule write_phase2_regenie_step1_command:
    input:
        config=RUN_CONFIG,
        pgen=f"{PHASE2_STEP1_QC_PREFIX}.pgen",
        pvar=f"{PHASE2_STEP1_QC_PREFIX}.pvar",
        psam=f"{PHASE2_STEP1_QC_PREFIX}.psam",
        variants=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.step1.snplist",
        pheno=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.pheno.tsv",
        covar=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.covar.tsv",
        keep=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.keep.txt",
        traits=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.analysis_traits.txt",
    output:
        command=f"results/gwas/PAN/regenie/groups/{{group}}/{{group}}.step1.command.sh",
    log:
        "results/logs/phase2_regenie/write_regenie_step1_command.{group}.log",
    threads:
        config["runtime"].get("threads_regenie_step1", config["runtime"]["threads_gwas"])
    resources:
        mem_mb=config["runtime"].get("mem_mb_regenie_step1", config["runtime"]["mem_mb_gwas"]),
        runtime=config["runtime"].get("time_min_regenie_step1", config["runtime"]["time_min_gwas"]),
    conda:
        "../../envs/gwas.yaml",
    params:
        pfile_prefix=lambda wildcards, input: str(input.pgen)[:-5],
        out_prefix=lambda wildcards: f"results/gwas/PAN/regenie/groups/{wildcards.group}/{wildcards.group}.step1",
        pred=lambda wildcards: f"results/gwas/PAN/regenie/groups/{wildcards.group}/{wildcards.group}.step1_pred.list",
    shell:
        """
        Rscript scripts/phase2_regenie.R write-step1-command \
          --config {input.config} \
          --group {wildcards.group} \
          --pfile-prefix {params.pfile_prefix} \
          --extract {input.variants} \
          --pheno {input.pheno} \
          --covar {input.covar} \
          --keep {input.keep} \
          --trait-list {input.traits} \
          --pred-list {params.pred} \
          --out-prefix {params.out_prefix} \
          --script-out {output.command} \
          --threads {threads} \
          > {log} 2>&1
        """


rule run_phase2_regenie_step1:
    input:
        command=f"results/gwas/PAN/regenie/groups/{{group}}/{{group}}.step1.command.sh",
        tool=PHASE2_REGENIE_TOOL,
    output:
        pred=f"results/gwas/PAN/regenie/groups/{{group}}/{{group}}.step1_pred.list",
        done=f"results/gwas/PAN/regenie/groups/{{group}}/{{group}}.step1.done",
    log:
        "results/logs/phase2_regenie/regenie_step1.{group}.log",
    threads:
        config["runtime"].get("threads_regenie_step1", config["runtime"]["threads_gwas"])
    resources:
        mem_mb=config["runtime"].get("mem_mb_regenie_step1", config["runtime"]["mem_mb_gwas"]),
        runtime=config["runtime"].get("time_min_regenie_step1", config["runtime"]["time_min_gwas"]),
    conda:
        "../../envs/regenie.yaml",
    shell:
        """
        bash {input.command:q} > {log:q} 2>&1
        """


rule write_phase2_regenie_step2_command:
    input:
        config=RUN_CONFIG,
        pgen=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.assoc_qc.pgen",
        pvar=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.assoc_qc.pvar",
        psam=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.assoc_qc.psam",
        pred=f"results/gwas/PAN/regenie/groups/{{group}}/{{group}}.step1_pred.list",
        pheno=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.pheno.tsv",
        covar=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.covar.tsv",
        traits=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.analysis_traits.txt",
        build="results/qc/genome_build/genome_build.txt",
    output:
        command=f"results/gwas/PAN/regenie/groups/{{group}}/{{group}}.{{build}}.step2.command.sh",
    log:
        "results/logs/phase2_regenie/write_regenie_step2_command.{group}.{build}.log",
    threads:
        config["runtime"].get("threads_regenie_step2", config["runtime"]["threads_gwas"])
    resources:
        mem_mb=config["runtime"].get("mem_mb_regenie_step2", config["runtime"]["mem_mb_gwas"]),
        runtime=config["runtime"].get("time_min_regenie_step2", config["runtime"]["time_min_gwas"]),
    conda:
        "../../envs/gwas.yaml",
    params:
        pfile_prefix=lambda wildcards, input: str(input.pgen)[:-5],
        out_prefix=lambda wildcards: f"results/gwas/PAN/regenie/groups/{wildcards.group}/{wildcards.group}.{wildcards.build}",
        done=lambda wildcards: f"results/gwas/PAN/regenie/groups/{wildcards.group}/{wildcards.group}.{wildcards.build}.step2.done",
    shell:
        """
        Rscript scripts/phase2_regenie.R write-step2-command \
          --config {input.config} \
          --group {wildcards.group} \
          --pfile-prefix {params.pfile_prefix} \
          --pheno {input.pheno} \
          --covar {input.covar} \
          --pred-list {input.pred} \
          --trait-list {input.traits} \
          --out-prefix {params.out_prefix} \
          --done {params.done} \
          --script-out {output.command} \
          --threads {threads} \
          > {log} 2>&1
        """


rule run_phase2_regenie_step2:
    input:
        command=f"results/gwas/PAN/regenie/groups/{{group}}/{{group}}.{{build}}.step2.command.sh",
        tool=PHASE2_REGENIE_TOOL,
    output:
        done=f"results/gwas/PAN/regenie/groups/{{group}}/{{group}}.{{build}}.step2.done",
    log:
        "results/logs/phase2_regenie/regenie_step2.{group}.{build}.log",
    threads:
        config["runtime"].get("threads_regenie_step2", config["runtime"]["threads_gwas"])
    resources:
        mem_mb=config["runtime"].get("mem_mb_regenie_step2", config["runtime"]["mem_mb_gwas"]),
        runtime=config["runtime"].get("time_min_regenie_step2", config["runtime"]["time_min_gwas"]),
    conda:
        "../../envs/regenie.yaml",
    shell:
        """
        bash {input.command:q} > {log:q} 2>&1
        """


rule stage_phase2_regenie_trait:
    input:
        config=RUN_CONFIG,
        done=phase2_trait_regenie_done,
        group_summary=lambda wildcards: f"{PHASE2_DIR}/groups/{phase2_trait_group(wildcards)}/{phase2_trait_group(wildcards)}.trait_summary.tsv",
    output:
        stats="results/gwas/{trait}/PAN/{trait}.PAN.{build}.regenie",
        summary="results/gwas/{trait}/PAN/{trait}.PAN.{build}.phase2_summary.tsv",
    log:
        "results/logs/phase2_regenie/stage_trait.{trait}.{build}.log",
    conda:
        "../../envs/gwas.yaml",
    params:
        group=lambda wildcards: phase2_trait_group(wildcards),
        raw_prefix=lambda wildcards: f"results/gwas/PAN/regenie/groups/{phase2_trait_group(wildcards)}/{phase2_trait_group(wildcards)}.{wildcards.build}",
    shell:
        """
        Rscript scripts/phase2_regenie.R stage-trait-output \
          --config {input.config} \
          --trait {wildcards.trait} \
          --group-summary {input.group_summary} \
          --raw-prefix {params.raw_prefix} \
          --out-stats {output.stats} \
          --out-summary {output.summary} \
          > {log} 2>&1
        """
