rule write_control_hwe_variants:
    input:
        config=RUN_CONFIG,
        pheno="results/qc/traits/{trait}.pheno.tsv",
        keep="results/qc/strata/{ancestry}.unrelated.keep.tsv",
        build="results/qc/genome_build/genome_build.txt",
    output:
        snplist="results/qc/gwas_hwe/{trait}/{ancestry}/{trait}.{ancestry}.{build}.control_hwe.snplist",
        plink_log="results/qc/gwas_hwe/{trait}/{ancestry}/{trait}.{ancestry}.{build}.control_hwe.log",
    log:
        "results/logs/gwas/{trait}.{ancestry}.{build}.control_hwe.log",
    threads:
        config["runtime"]["threads_gwas"]
    resources:
        mem_mb=config["runtime"]["mem_mb_gwas"],
        runtime=config["runtime"]["time_min_gwas"],
    conda:
        "../../envs/gwas.yaml",
    params:
        prefix=lambda wildcards: f"results/qc/gwas_hwe/{wildcards.trait}/{wildcards.ancestry}/{wildcards.trait}.{wildcards.ancestry}.{wildcards.build}.control_hwe",
    shell:
        """
        Rscript scripts/write_control_hwe_variants.R \
          --config {input.config} \
          --pheno {input.pheno} \
          --keep {input.keep} \
          --plink-prefix {params.prefix} \
          --out {output.snplist} \
          --threads {threads} \
          > {log} 2>&1
        """


rule run_plink2_gwas:
    input:
        config=RUN_CONFIG,
        pheno="results/qc/traits/{trait}.pheno.tsv",
        covar="results/qc/traits/{trait}.covar.tsv",
        keep="results/qc/strata/{ancestry}.unrelated.keep.tsv",
        hwe_snplist="results/qc/gwas_hwe/{trait}/{ancestry}/{trait}.{ancestry}.{build}.control_hwe.snplist",
        build="results/qc/genome_build/genome_build.txt",
    output:
        stats="results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.plink2.glm.tsv",
        plink_log="results/gwas/{trait}/{ancestry}/plink2_raw/{trait}.{ancestry}.{build}.log",
    log:
        "results/logs/gwas/{trait}.{ancestry}.{build}.plink2.log",
    threads:
        config["runtime"]["threads_gwas"]
    resources:
        mem_mb=config["runtime"]["mem_mb_gwas"],
        runtime=config["runtime"]["time_min_gwas"],
    conda:
        "../../envs/gwas.yaml",
    params:
        prefix=lambda wildcards: f"results/gwas/{wildcards.trait}/{wildcards.ancestry}/plink2_raw/{wildcards.trait}.{wildcards.ancestry}.{wildcards.build}",
    shell:
        """
        Rscript scripts/run_plink2_gwas.R \
          --config {input.config} \
          --trait {wildcards.trait} \
          --ancestry {wildcards.ancestry} \
          --build {wildcards.build} \
          --pheno {input.pheno} \
          --covar {input.covar} \
          --keep {input.keep} \
          --hwe-snplist {input.hwe_snplist} \
          --plink-prefix {params.prefix} \
          --out {output.stats} \
          --threads {threads} \
          > {log} 2>&1
        """


rule summarize_gwas_run:
    input:
        config=RUN_CONFIG,
        stats="results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.plink2.glm.tsv",
        plink_log="results/gwas/{trait}/{ancestry}/plink2_raw/{trait}.{ancestry}.{build}.log",
        hwe_snplist="results/qc/gwas_hwe/{trait}/{ancestry}/{trait}.{ancestry}.{build}.control_hwe.snplist",
        hwe_log="results/qc/gwas_hwe/{trait}/{ancestry}/{trait}.{ancestry}.{build}.control_hwe.log",
        pheno="results/qc/traits/{trait}.pheno.tsv",
        covar="results/qc/traits/{trait}.covar.tsv",
        keep="results/qc/strata/{ancestry}.unrelated.keep.tsv",
    output:
        summary="results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.gwas_filter_summary.tsv",
    log:
        "results/logs/gwas/{trait}.{ancestry}.{build}.summary.log",
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/summarize_gwas_run.R \
          --config {input.config} \
          --stats {input.stats} \
          --plink-log {input.plink_log} \
          --hwe-snplist {input.hwe_snplist} \
          --hwe-log {input.hwe_log} \
          --pheno {input.pheno} \
          --covar {input.covar} \
          --keep {input.keep} \
          --out {output.summary} \
          > {log} 2>&1
        """
