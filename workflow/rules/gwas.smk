rule run_plink2_gwas:
    input:
        config=RUN_CONFIG,
        pheno="results/qc/traits/{trait}.pheno.tsv",
        covar="results/qc/traits/{trait}.covar.tsv",
        keep="results/qc/strata/{ancestry}.unrelated.keep.tsv",
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
          --plink-prefix {params.prefix} \
          --out {output.stats} \
          --threads {threads} \
          > {log} 2>&1
        """
