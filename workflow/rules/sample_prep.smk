rule run_sex_check:
    input:
        config=RUN_CONFIG,
        ok="results/qc/input_validation/validation.ok",
    output:
        sexcheck="results/qc/sex/sexcheck.tsv",
        drop="results/qc/sex/sex_mismatches.remove.tsv",
        keep="results/qc/sex/sex_checked.keep.tsv",
        summary="results/qc/sex/sex_check_summary.tsv",
    log:
        "results/logs/sample_prep/sex_check.log",
    threads:
        config["runtime"]["threads_small"]
    resources:
        mem_mb=config["runtime"]["mem_mb_small"],
        runtime=config["runtime"]["time_min_small"],
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/sex_check.R \
          --config {input.config} \
          --plink-out-prefix results/qc/sex/plink_sex_check \
          --sexcheck-out {output.sexcheck} \
          --remove-out {output.drop} \
          --keep-out {output.keep} \
          --summary-out {output.summary} \
          --threads {threads} \
          > {log} 2>&1
        """


rule make_strata_files:
    input:
        config=RUN_CONFIG,
        ok="results/qc/input_validation/validation.ok",
        ancestry=ancestry_file(),
    output:
        counts="results/qc/strata/strata_counts.tsv",
        keep=expand("results/qc/strata/{ancestry}.keep.tsv", ancestry=ANCESTRIES),
    log:
        "results/logs/sample_prep/make_strata_files.log",
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/make_strata_files.R \
          --config {input.config} \
          --ancestry-file {input.ancestry} \
          --outdir results/qc/strata \
          > {log} 2>&1
        """


rule build_trait_files:
    input:
        config=RUN_CONFIG,
        ok="results/qc/input_validation/validation.ok",
        pcs=pcs_file(),
        # All final keep files are inputs so PC completeness is checked only
        # for samples that can actually enter one of the ancestry GWAS jobs.
        keep=expand("results/qc/strata/{ancestry}.unrelated.keep.tsv", ancestry=ANCESTRIES),
    output:
        pheno="results/qc/traits/{trait}.pheno.tsv",
        covar="results/qc/traits/{trait}.covar.tsv",
    log:
        "results/logs/sample_prep/build_trait_files.{trait}.log",
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/build_trait_files.R \
          --config {input.config} \
          --trait {wildcards.trait} \
          --pcs-file {input.pcs} \
          --keep {input.keep} \
          --pheno-out {output.pheno} \
          --covar-out {output.covar} \
          > {log} 2>&1
        """


rule make_unrelated_keep:
    input:
        config=RUN_CONFIG,
        ok="results/qc/input_validation/validation.ok",
        pgen="results/qc/relatedness/relatedness_qc.pgen",
        pvar="results/qc/relatedness/relatedness_qc.pvar",
        psam="results/qc/relatedness/relatedness_qc.psam",
        markers="results/qc/relatedness/relatedness_ld_prune.prune.in",
    output:
        keep="results/qc/relatedness/unrelated.king.cutoff.in.id",
    log:
        "results/logs/sample_prep/make_unrelated_keep.log",
    threads:
        config["runtime"]["threads_small"]
    resources:
        mem_mb=config["runtime"]["mem_mb_small"],
        runtime=config["runtime"]["time_min_small"],
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/make_unrelated_keep.R \
          --config {input.config} \
          --out {output.keep} \
          --pfile-prefix results/qc/relatedness/relatedness_qc \
          --extract {input.markers} \
          --plink-out-prefix results/qc/relatedness/unrelated \
          --threads {threads} \
          > {log} 2>&1
        """


rule prepare_relatedness_markers:
    input:
        config=RUN_CONFIG,
        ok="results/qc/input_validation/validation.ok",
        sex_keep="results/qc/sex/sex_checked.keep.tsv",
    output:
        pgen="results/qc/relatedness/relatedness_qc.pgen",
        pvar="results/qc/relatedness/relatedness_qc.pvar",
        psam="results/qc/relatedness/relatedness_qc.psam",
        prune_in="results/qc/relatedness/relatedness_ld_prune.prune.in",
        prune_out="results/qc/relatedness/relatedness_ld_prune.prune.out",
        summary="results/qc/relatedness/relatedness_summary.tsv",
    log:
        "results/logs/sample_prep/prepare_relatedness_markers.log",
    threads:
        config["runtime"]["threads_small"]
    resources:
        mem_mb=config["runtime"]["mem_mb_small"],
        runtime=config["runtime"]["time_min_small"],
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/prepare_relatedness_markers.R \
          --config {input.config} \
          --keep {input.sex_keep} \
          --out-prefix results/qc/relatedness/relatedness_qc \
          --prune-prefix results/qc/relatedness/relatedness_ld_prune \
          --summary {output.summary} \
          --threads {threads} \
          > {log} 2>&1
        """


rule intersect_stratum_unrelated:
    input:
        stratum="results/qc/strata/{ancestry}.keep.tsv",
        unrelated="results/qc/relatedness/unrelated.king.cutoff.in.id",
        # sex_check.action controls whether this is all samples or excludes
        # genetic-sex mismatches; downstream rules just consume the keep list.
        sex_keep="results/qc/sex/sex_checked.keep.tsv",
    output:
        keep="results/qc/strata/{ancestry}.unrelated.keep.tsv",
    log:
        "results/logs/sample_prep/intersect_stratum_unrelated.{ancestry}.log",
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/intersect_keep_files.R \
          --inputs {input.stratum} {input.unrelated} {input.sex_keep} \
          --out {output.keep} \
          > {log} 2>&1
        """
