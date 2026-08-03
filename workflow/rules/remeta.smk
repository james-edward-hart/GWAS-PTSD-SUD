REMETA_DIR = "results/remeta"
REMETA_TOOL = f"{REMETA_DIR}/work/remeta_tool.tsv"


rule record_remeta_tool:
    input:
        config=RUN_CONFIG,
    output:
        tool=REMETA_TOOL,
    log:
        "results/logs/remeta/record_tool.log",
    conda:
        "../../envs/remeta.yaml",
    params:
        tool=lambda wildcards: config.get("tools", {}).get("remeta", "remeta"),
    shell:
        """
        python scripts/record_remeta_tool.py \
          --tool {params.tool:q} \
          --out {output.tool:q} \
          > {log:q} 2>&1
        """


rule prepare_remeta_target_genotypes:
    input:
        config=RUN_CONFIG,
        build="results/qc/genome_build/genome_build.txt",
        pgen=f"{PHASE2_PAN_PREFIX}.pgen",
        pvar=f"{PHASE2_PAN_PREFIX}.pvar",
        psam=f"{PHASE2_PAN_PREFIX}.psam",
        keep=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.keep.txt",
        regions=lambda wildcards: remeta_resource_file(wildcards.build, "target_regions.bed"),
        provenance=lambda wildcards: remeta_resource_file(wildcards.build, "provenance.tsv"),
    output:
        pgen=f"{REMETA_DIR}/work/{{build}}/groups/{{group}}/{{group}}.target.pgen",
        pvar=f"{REMETA_DIR}/work/{{build}}/groups/{{group}}/{{group}}.target.pvar",
        psam=f"{REMETA_DIR}/work/{{build}}/groups/{{group}}/{{group}}.target.psam",
        summary=f"{REMETA_DIR}/work/{{build}}/groups/{{group}}/{{group}}.target.summary.tsv",
    log:
        "results/logs/remeta/prepare_target.{group}.{build}.log",
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
        Rscript scripts/remeta_cohort.R prepare-target \
          --config {input.config:q} \
          --pfile-prefix {params.pfile_prefix:q} \
          --keep {input.keep:q} \
          --regions {input.regions:q} \
          --out-prefix {params.out_prefix:q} \
          --summary-out {output.summary:q} \
          --threads {threads} \
          > {log:q} 2>&1
        """


rule write_remeta_regenie_step2_command:
    input:
        config=RUN_CONFIG,
        pgen=f"{REMETA_DIR}/work/{{build}}/groups/{{group}}/{{group}}.target.pgen",
        pvar=f"{REMETA_DIR}/work/{{build}}/groups/{{group}}/{{group}}.target.pvar",
        psam=f"{REMETA_DIR}/work/{{build}}/groups/{{group}}/{{group}}.target.psam",
        keep=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.keep.txt",
        pred="results/gwas/PAN/regenie/groups/{group}/{group}.step1_pred.list",
        loco="results/gwas/PAN/regenie/groups/{group}/{group}.step1_1.loco",
        pheno=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.pheno.tsv",
        covar=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.covar.tsv",
        traits=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.analysis_traits.txt",
        covars=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.covariates.txt",
        summary=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.trait_summary.tsv",
    output:
        command=f"{REMETA_DIR}/work/{{build}}/groups/{{group}}/{{group}}.rare_{{trait}}.step2.command.sh",
    log:
        "results/logs/remeta/write_regenie_step2.{group}.{build}.{trait}.log",
    threads:
        config["runtime"].get("threads_regenie_step2", config["runtime"]["threads_gwas"])
    resources:
        mem_mb=config["runtime"].get("mem_mb_regenie_step2", config["runtime"]["mem_mb_gwas"]),
        runtime=config["runtime"].get("time_min_regenie_step2", config["runtime"]["time_min_gwas"]),
    conda:
        "../../envs/gwas.yaml",
    params:
        pfile_prefix=lambda wildcards, input: str(input.pgen)[:-5],
        out_prefix=lambda wildcards: f"{REMETA_DIR}/work/{wildcards.build}/groups/{wildcards.group}/{wildcards.group}.rare",
        done=lambda wildcards: f"{REMETA_DIR}/work/{wildcards.build}/groups/{wildcards.group}/{wildcards.group}.rare_{wildcards.trait}.step2.done",
    shell:
        """
        Rscript scripts/remeta_cohort.R write-step2-command \
          --config {input.config:q} \
          --trait {wildcards.trait:q} \
          --group-summary {input.summary:q} \
          --pfile-prefix {params.pfile_prefix:q} \
          --keep {input.keep:q} \
          --pheno {input.pheno:q} \
          --covar {input.covar:q} \
          --pred-list {input.pred:q} \
          --trait-list {input.traits:q} \
          --covar-list {input.covars:q} \
          --out-prefix {params.out_prefix:q} \
          --done {params.done:q} \
          --script-out {output.command:q} \
          --threads {threads} \
          > {log:q} 2>&1
        """


rule run_remeta_regenie_step2:
    input:
        command=f"{REMETA_DIR}/work/{{build}}/groups/{{group}}/{{group}}.rare_{{trait}}.step2.command.sh",
        tool=PHASE2_REGENIE_TOOL,
    output:
        stats=f"{REMETA_DIR}/work/{{build}}/groups/{{group}}/{{group}}.rare_{{trait}}.regenie.gz",
        ids=f"{REMETA_DIR}/work/{{build}}/groups/{{group}}/{{group}}.rare_{{trait}}.regenie.ids",
        done=f"{REMETA_DIR}/work/{{build}}/groups/{{group}}/{{group}}.rare_{{trait}}.step2.done",
    log:
        "results/logs/remeta/regenie_step2.{group}.{build}.{trait}.log",
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


rule stage_remeta_regenie_trait:
    input:
        stats=lambda wildcards: f"{REMETA_DIR}/work/{wildcards.build}/groups/{phase2_trait_group(wildcards)}/{phase2_trait_group(wildcards)}.rare_{wildcards.trait}.regenie.gz",
        ids=lambda wildcards: f"{REMETA_DIR}/work/{wildcards.build}/groups/{phase2_trait_group(wildcards)}/{phase2_trait_group(wildcards)}.rare_{wildcards.trait}.regenie.ids",
        keep=lambda wildcards: f"{PHASE2_DIR}/groups/{phase2_trait_group(wildcards)}/{phase2_trait_group(wildcards)}.keep.txt",
        summary=lambda wildcards: f"{PHASE2_DIR}/groups/{phase2_trait_group(wildcards)}/{phase2_trait_group(wildcards)}.trait_summary.tsv",
    output:
        htp=f"{REMETA_DIR}/export/{{build}}/htp/{{trait}}.PAN.regenie.gz",
    log:
        "results/logs/remeta/stage_trait.{trait}.{build}.log",
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/remeta_cohort.R stage-trait \
          --trait {wildcards.trait:q} \
          --group-summary {input.summary:q} \
          --raw-stats {input.stats:q} \
          --sample-ids {input.ids:q} \
          --keep {input.keep:q} \
          --out {output.htp:q} \
          > {log:q} 2>&1
        """


rule compute_remeta_marginal_ld:
    input:
        tool=REMETA_TOOL,
        pgen=f"{REMETA_DIR}/work/{{build}}/groups/{{group}}/{{group}}.target.pgen",
        pvar=f"{REMETA_DIR}/work/{{build}}/groups/{{group}}/{{group}}.target.pvar",
        psam=f"{REMETA_DIR}/work/{{build}}/groups/{{group}}/{{group}}.target.psam",
        genes=lambda wildcards: remeta_resource_file(wildcards.build, "gene_list.tsv"),
    output:
        gene_ld=f"{REMETA_DIR}/export/{{build}}/ld/{{group}}/chr{{chrom}}.remeta.gene.ld",
        buffer_ld=f"{REMETA_DIR}/export/{{build}}/ld/{{group}}/chr{{chrom}}.remeta.buffer.ld",
        index=f"{REMETA_DIR}/export/{{build}}/ld/{{group}}/chr{{chrom}}.remeta.ld.idx.gz",
    log:
        "results/logs/remeta/compute_ld.{group}.{build}.chr{chrom}.log",
    threads:
        config["runtime"].get("threads_remeta_ld", 4)
    resources:
        mem_mb=config["runtime"].get("mem_mb_remeta_ld", 16000),
        runtime=config["runtime"].get("time_min_remeta_ld", 240),
        # Global profile capacity limits simultaneous readers of the shared PGEN.
        remeta_ld_jobs=1,
    conda:
        "../../envs/remeta.yaml",
    params:
        target=lambda wildcards, input: str(input.pgen)[:-5],
        out=lambda wildcards, output: str(output.gene_ld)[:-len(".remeta.gene.ld")],
        remeta=lambda wildcards: config.get("tools", {}).get("remeta", "remeta"),
        target_r2=lambda wildcards: config.get("remeta", {}).get("target_r2", 0.0001),
        dosage=lambda wildcards: "--use-dosages"
        if str(config.get("remeta", {}).get("genotype_mode", "")).lower() == "dosage"
        else "",
    shell:
        """
        {params.remeta:q} compute-ref-ld \
          --target-pfile {params.target:q} \
          --gene-list {input.genes:q} \
          --chr {wildcards.chrom:q} \
          --out {params.out:q} \
          --target-r2 {params.target_r2} \
          --skip-buffer \
          --threads {threads} \
          {params.dosage} \
          > {log:q} 2>&1
        """


rule validate_remeta_group:
    input:
        pgen=f"{REMETA_DIR}/work/{{build}}/groups/{{group}}/{{group}}.target.pgen",
        pvar=f"{REMETA_DIR}/work/{{build}}/groups/{{group}}/{{group}}.target.pvar",
        psam=f"{REMETA_DIR}/work/{{build}}/groups/{{group}}/{{group}}.target.psam",
        keep=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.keep.txt",
        ordinary_ids=lambda wildcards: f"results/gwas/PAN/regenie/groups/{wildcards.group}/{wildcards.group}.{wildcards.build}_{phase2_group_trait(wildcards)}.regenie.ids",
        rare_ids=lambda wildcards: f"{REMETA_DIR}/work/{wildcards.build}/groups/{wildcards.group}/{wildcards.group}.rare_{phase2_group_trait(wildcards)}.regenie.ids",
        summary=f"{PHASE2_DIR}/groups/{{group}}/{{group}}.trait_summary.tsv",
        genes=lambda wildcards: remeta_resource_file(wildcards.build, "gene_list.tsv"),
        htp=remeta_group_htp_inputs,
        ld=remeta_group_ld_inputs,
    output:
        ok=f"{REMETA_DIR}/work/{{build}}/groups/{{group}}/{{group}}.validation.ok",
    log:
        "results/logs/remeta/validate_group.{group}.{build}.log",
    resources:
        # Validation scans every BGZF LD component and can approach LD-job cost.
        mem_mb=config["runtime"].get("mem_mb_remeta_ld", 16000),
        runtime=config["runtime"].get("time_min_remeta_ld", 240),
    conda:
        "../../envs/gwas.yaml",
    params:
        target=lambda wildcards, input: str(input.pgen)[:-5],
    shell:
        """
        Rscript scripts/remeta_cohort.R validate-group \
          --target-prefix {params.target:q} \
          --keep {input.keep:q} \
          --ordinary-sample-ids {input.ordinary_ids:q} \
          --sample-ids {input.rare_ids:q} \
          --group-summary {input.summary:q} \
          --gene-list {input.genes:q} \
          --htp {input.htp:q} \
          --ld {input.ld:q} \
          --out {output.ok:q} \
          > {log:q} 2>&1
        """


rule write_remeta_cohort_manifest:
    input:
        config=RUN_CONFIG,
        build="results/qc/genome_build/genome_build.txt",
        genes=lambda wildcards: remeta_resource_file(wildcards.build, "gene_list.tsv"),
        provenance=lambda wildcards: remeta_resource_file(wildcards.build, "provenance.tsv"),
        status=phase2_status_file,
        target_summary=remeta_target_summaries,
        trait_summary=remeta_trait_summaries,
        validation=remeta_validations,
        tools=[REMETA_TOOL, PHASE2_REGENIE_TOOL],
        artifacts=remeta_manifest_inputs,
    output:
        manifest=f"{REMETA_DIR}/export/{ANALYSIS_OUTPUT_NAME}.{{build}}.remeta_manifest.tsv",
    log:
        "results/logs/remeta/write_manifest.{build}.log",
    conda:
        "../../envs/gwas.yaml",
    params:
        target_summary_args=lambda wildcards, input: optional_path_args("--target-summary", input.target_summary),
        trait_summary_args=lambda wildcards, input: optional_path_args("--trait-summary", input.trait_summary),
        validation_args=lambda wildcards, input: optional_path_args("--validation", input.validation),
        tool_args=lambda wildcards, input: optional_path_args("--tool", input.tools),
        artifact_args=lambda wildcards, input: optional_path_args("--artifact", input.artifacts),
    shell:
        """
        Rscript scripts/remeta_cohort.R write-manifest \
          --config {input.config:q} \
          --build {wildcards.build:q} \
          --gene-list {input.genes:q} \
          --provenance {input.provenance:q} \
          --group-status {input.status:q} \
          {params.target_summary_args} \
          {params.trait_summary_args} \
          {params.validation_args} \
          {params.tool_args} \
          {params.artifact_args} \
          --out {output.manifest:q} \
          > {log:q} 2>&1
        """
