ANCESTRY_REF_DIR = "results/qc/ancestry/reference"
ANCESTRY_WITHIN_DIR = "results/qc/ancestry/within"

REFERENCE_QC_PREFIX = f"{ANCESTRY_REF_DIR}/reference_qc"
STUDY_QC_PREFIX = f"{ANCESTRY_REF_DIR}/study_qc"
REFERENCE_SHARED_PREFIX = f"{ANCESTRY_REF_DIR}/reference_shared"
STUDY_SHARED_PREFIX = f"{ANCESTRY_REF_DIR}/study_shared"
LD_PRUNE_PREFIX = f"{ANCESTRY_REF_DIR}/ld_prune/ancestry_ld_prune"
REFERENCE_PCA_PREFIX = f"{ANCESTRY_REF_DIR}/reference_pca/reference"
REFERENCE_SCORE_PREFIX = f"{ANCESTRY_REF_DIR}/reference_pca/reference_projected"
STUDY_SCORE_PREFIX = f"{ANCESTRY_REF_DIR}/reference_pca/study_projected"

rule convert_reference_for_ancestry:
    input:
        config=RUN_CONFIG,
        ok="results/qc/input_validation/validation.ok",
    output:
        pgen=f"{REFERENCE_QC_PREFIX}.pgen",
        pvar=f"{REFERENCE_QC_PREFIX}.pvar",
        psam=f"{REFERENCE_QC_PREFIX}.psam",
    log:
        "results/logs/ancestry/reference/convert_reference.log",
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
        Rscript scripts/ancestry_reference.R convert-reference \
          --config {input.config} \
          --out-prefix {params.out_prefix} \
          --threads {threads} \
          > {log} 2>&1
        """


rule convert_study_for_ancestry:
    input:
        config=RUN_CONFIG,
        ok="results/qc/input_validation/validation.ok",
    output:
        pgen=f"{STUDY_QC_PREFIX}.pgen",
        pvar=f"{STUDY_QC_PREFIX}.pvar",
        psam=f"{STUDY_QC_PREFIX}.psam",
    log:
        "results/logs/ancestry/reference/convert_study.log",
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
        Rscript scripts/ancestry_reference.R convert-study \
          --config {input.config} \
          --out-prefix {params.out_prefix} \
          --threads {threads} \
          > {log} 2>&1
        """


rule write_ancestry_shared_variants:
    input:
        config=RUN_CONFIG,
        reference=f"{REFERENCE_QC_PREFIX}.pvar",
        study=f"{STUDY_QC_PREFIX}.pvar",
    output:
        variants=f"{ANCESTRY_REF_DIR}/shared_variants.txt",
        mismatches=f"{ANCESTRY_REF_DIR}/shared_variant_mismatches.tsv",
    log:
        "results/logs/ancestry/reference/write_shared_variants.log",
    conda:
        "../../envs/gwas.yaml",
    params:
        reference_prefix=lambda wildcards, input: str(input.reference)[:-5],
        study_prefix=lambda wildcards, input: str(input.study)[:-5],
    shell:
        """
        Rscript scripts/ancestry_reference.R shared-variants \
          --config {input.config} \
          --reference-prefix {params.reference_prefix} \
          --study-prefix {params.study_prefix} \
          --out {output.variants} \
          --mismatch-report {output.mismatches} \
          > {log} 2>&1
        """


rule extract_shared_reference:
    input:
        config=RUN_CONFIG,
        pgen=f"{REFERENCE_QC_PREFIX}.pgen",
        pvar=f"{REFERENCE_QC_PREFIX}.pvar",
        psam=f"{REFERENCE_QC_PREFIX}.psam",
        variants=f"{ANCESTRY_REF_DIR}/shared_variants.txt",
    output:
        pgen=f"{REFERENCE_SHARED_PREFIX}.pgen",
        pvar=f"{REFERENCE_SHARED_PREFIX}.pvar",
        psam=f"{REFERENCE_SHARED_PREFIX}.psam",
    log:
        "results/logs/ancestry/reference/extract_shared_reference.log",
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
        Rscript scripts/ancestry_reference.R extract-shared \
          --config {input.config} \
          --input-prefix {params.input_prefix} \
          --variants {input.variants} \
          --out-prefix {params.out_prefix} \
          --threads {threads} \
          > {log} 2>&1
        """


rule extract_shared_study:
    input:
        config=RUN_CONFIG,
        pgen=f"{STUDY_QC_PREFIX}.pgen",
        pvar=f"{STUDY_QC_PREFIX}.pvar",
        psam=f"{STUDY_QC_PREFIX}.psam",
        variants=f"{ANCESTRY_REF_DIR}/shared_variants.txt",
    output:
        pgen=f"{STUDY_SHARED_PREFIX}.pgen",
        pvar=f"{STUDY_SHARED_PREFIX}.pvar",
        psam=f"{STUDY_SHARED_PREFIX}.psam",
    log:
        "results/logs/ancestry/reference/extract_shared_study.log",
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
        Rscript scripts/ancestry_reference.R extract-shared \
          --config {input.config} \
          --input-prefix {params.input_prefix} \
          --variants {input.variants} \
          --out-prefix {params.out_prefix} \
          --threads {threads} \
          > {log} 2>&1
        """


rule ld_prune_ancestry_markers:
    input:
        config=RUN_CONFIG,
        pgen=f"{REFERENCE_SHARED_PREFIX}.pgen",
        pvar=f"{REFERENCE_SHARED_PREFIX}.pvar",
        psam=f"{REFERENCE_SHARED_PREFIX}.psam",
        shared=f"{ANCESTRY_REF_DIR}/shared_variants.txt",
    output:
        prune_in=f"{LD_PRUNE_PREFIX}.prune.in",
        prune_out=f"{LD_PRUNE_PREFIX}.prune.out",
        excluded=f"{LD_PRUNE_PREFIX}.excluded_region_variants.txt",
    log:
        "results/logs/ancestry/reference/ld_prune.log",
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
        Rscript scripts/ancestry_reference.R ld-prune \
          --config {input.config} \
          --pfile-prefix {params.pfile_prefix} \
          --shared-variants {input.shared} \
          --out-prefix {params.out_prefix} \
          --prune-in {output.prune_in} \
          --excluded-regions {output.excluded} \
          --threads {threads} \
          > {log} 2>&1
        """


rule fit_reference_pca:
    input:
        config=RUN_CONFIG,
        pgen=f"{REFERENCE_SHARED_PREFIX}.pgen",
        pvar=f"{REFERENCE_SHARED_PREFIX}.pvar",
        psam=f"{REFERENCE_SHARED_PREFIX}.psam",
        variants=f"{LD_PRUNE_PREFIX}.prune.in",
    output:
        eigenvec=f"{REFERENCE_PCA_PREFIX}.eigenvec",
        eigenval=f"{REFERENCE_PCA_PREFIX}.eigenval",
        allele=f"{REFERENCE_PCA_PREFIX}.eigenvec.allele",
        acount=f"{REFERENCE_PCA_PREFIX}.acount",
    log:
        "results/logs/ancestry/reference/fit_reference_pca.log",
    threads:
        config["runtime"]["threads_small"]
    resources:
        mem_mb=config["runtime"]["mem_mb_small"],
        runtime=config["runtime"]["time_min_small"],
    conda:
        "../../envs/gwas.yaml",
    params:
        pfile_prefix=lambda wildcards, input: str(input.pgen)[:-5],
        out_prefix=lambda wildcards, output: str(output.eigenvec).replace(".eigenvec", ""),
    shell:
        """
        Rscript scripts/ancestry_reference.R fit-reference-pca \
          --config {input.config} \
          --pfile-prefix {params.pfile_prefix} \
          --variants {input.variants} \
          --out-prefix {params.out_prefix} \
          --threads {threads} \
          > {log} 2>&1
        """


rule score_reference_pcs:
    input:
        config=RUN_CONFIG,
        pgen=f"{REFERENCE_SHARED_PREFIX}.pgen",
        pvar=f"{REFERENCE_SHARED_PREFIX}.pvar",
        psam=f"{REFERENCE_SHARED_PREFIX}.psam",
        variants=f"{LD_PRUNE_PREFIX}.prune.in",
        weights=f"{REFERENCE_PCA_PREFIX}.eigenvec.allele",
        frequencies=f"{REFERENCE_PCA_PREFIX}.acount",
    output:
        sscore=f"{REFERENCE_SCORE_PREFIX}.sscore",
    log:
        "results/logs/ancestry/reference/score_reference_pcs.log",
    threads:
        config["runtime"]["threads_small"]
    resources:
        mem_mb=config["runtime"]["mem_mb_small"],
        runtime=config["runtime"]["time_min_small"],
    conda:
        "../../envs/gwas.yaml",
    params:
        pfile_prefix=lambda wildcards, input: str(input.pgen)[:-5],
        out_prefix=lambda wildcards, output: str(output.sscore).replace(".sscore", ""),
    shell:
        """
        Rscript scripts/ancestry_reference.R score-pcs \
          --config {input.config} \
          --pfile-prefix {params.pfile_prefix} \
          --variants {input.variants} \
          --weights {input.weights} \
          --frequencies {input.frequencies} \
          --out-prefix {params.out_prefix} \
          --threads {threads} \
          > {log} 2>&1
        """


rule score_study_pcs:
    input:
        config=RUN_CONFIG,
        pgen=f"{STUDY_SHARED_PREFIX}.pgen",
        pvar=f"{STUDY_SHARED_PREFIX}.pvar",
        psam=f"{STUDY_SHARED_PREFIX}.psam",
        variants=f"{LD_PRUNE_PREFIX}.prune.in",
        weights=f"{REFERENCE_PCA_PREFIX}.eigenvec.allele",
        frequencies=f"{REFERENCE_PCA_PREFIX}.acount",
    output:
        sscore=f"{STUDY_SCORE_PREFIX}.sscore",
    log:
        "results/logs/ancestry/reference/score_study_pcs.log",
    threads:
        config["runtime"]["threads_small"]
    resources:
        mem_mb=config["runtime"]["mem_mb_small"],
        runtime=config["runtime"]["time_min_small"],
    conda:
        "../../envs/gwas.yaml",
    params:
        pfile_prefix=lambda wildcards, input: str(input.pgen)[:-5],
        out_prefix=lambda wildcards, output: str(output.sscore).replace(".sscore", ""),
    shell:
        """
        Rscript scripts/ancestry_reference.R score-pcs \
          --config {input.config} \
          --pfile-prefix {params.pfile_prefix} \
          --variants {input.variants} \
          --weights {input.weights} \
          --frequencies {input.frequencies} \
          --out-prefix {params.out_prefix} \
          --threads {threads} \
          > {log} 2>&1
        """


rule write_reference_pcs:
    input:
        config=RUN_CONFIG,
        sscore=f"{REFERENCE_SCORE_PREFIX}.sscore",
    output:
        pcs=f"{ANCESTRY_REF_DIR}/reference_pcs.tsv",
    log:
        "results/logs/ancestry/reference/write_reference_pcs.log",
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/ancestry_reference.R write-reference-pcs \
          --config {input.config} \
          --sscore {input.sscore} \
          --out {output.pcs} \
          > {log} 2>&1
        """


rule write_study_projected_pcs:
    input:
        config=RUN_CONFIG,
        sscore=f"{STUDY_SCORE_PREFIX}.sscore",
    output:
        pcs=f"{ANCESTRY_REF_DIR}/study_projected_pcs.tsv",
    log:
        "results/logs/ancestry/reference/write_study_projected_pcs.log",
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/ancestry_reference.R write-study-pcs \
          --config {input.config} \
          --sscore {input.sscore} \
          --out {output.pcs} \
          > {log} 2>&1
        """


rule write_ancestry_reference_report:
    input:
        config=RUN_CONFIG,
        variants=f"{ANCESTRY_REF_DIR}/shared_variants.txt",
        prune_in=f"{LD_PRUNE_PREFIX}.prune.in",
        mismatches=f"{ANCESTRY_REF_DIR}/shared_variant_mismatches.tsv",
        projection_validation=f"{ANCESTRY_REF_DIR}/reference_projection_validation.tsv",
        reference_pcs=f"{ANCESTRY_REF_DIR}/reference_pcs.tsv",
        study_pcs=f"{ANCESTRY_REF_DIR}/study_projected_pcs.tsv",
    output:
        report=f"{ANCESTRY_REF_DIR}/reference_prep_report.md",
    log:
        "results/logs/ancestry/reference/write_reference_report.log",
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/ancestry_reference.R write-report \
          --config {input.config} \
          --shared-variants {input.variants} \
          --prune-in {input.prune_in} \
          --mismatch-report {input.mismatches} \
          --projection-validation {input.projection_validation} \
          --reference-pcs {input.reference_pcs} \
          --study-pcs {input.study_pcs} \
          --out {output.report} \
          > {log} 2>&1
        """


rule validate_reference_projection:
    input:
        config=RUN_CONFIG,
        eigenvec=f"{REFERENCE_PCA_PREFIX}.eigenvec",
        projected=f"{ANCESTRY_REF_DIR}/reference_pcs.tsv",
    output:
        validation=f"{ANCESTRY_REF_DIR}/reference_projection_validation.tsv",
    log:
        "results/logs/ancestry/reference/validate_reference_projection.log",
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/ancestry_reference.R validate-projection \
          --config {input.config} \
          --eigenvec {input.eigenvec} \
          --projected-pcs {input.projected} \
          --out {output.validation} \
          > {log} 2>&1
        """


rule run_within_ancestry_pca:
    input:
        config=RUN_CONFIG,
        pgen=f"{STUDY_SHARED_PREFIX}.pgen",
        pvar=f"{STUDY_SHARED_PREFIX}.pvar",
        psam=f"{STUDY_SHARED_PREFIX}.psam",
        variants=f"{LD_PRUNE_PREFIX}.prune.in",
        keep="results/qc/strata/{ancestry}.unrelated.keep.tsv",
    output:
        eigenvec=f"{ANCESTRY_WITHIN_DIR}/{{ancestry}}.eigenvec",
        eigenval=f"{ANCESTRY_WITHIN_DIR}/{{ancestry}}.eigenval",
    log:
        "results/logs/ancestry/within_pca.{ancestry}.log",
    threads:
        config["runtime"]["threads_small"]
    resources:
        mem_mb=config["runtime"]["mem_mb_small"],
        runtime=config["runtime"]["time_min_small"],
    conda:
        "../../envs/gwas.yaml",
    params:
        prefix=lambda wildcards: f"{ANCESTRY_WITHIN_DIR}/{wildcards.ancestry}",
        pfile_prefix=lambda wildcards, input: str(input.pgen)[:-5],
    shell:
        """
        Rscript scripts/ancestry_reference.R within-ancestry-pca \
          --config {input.config} \
          --pfile-prefix {params.pfile_prefix} \
          --keep {input.keep} \
          --variants {input.variants} \
          --out-prefix {params.prefix} \
          --threads {threads} \
          > {log} 2>&1
        """


rule combine_within_ancestry_pcs:
    input:
        config=RUN_CONFIG,
        eigenvecs=active_within_ancestry_eigenvecs,
    output:
        pcs="results/qc/ancestry/within_ancestry_pcs.tsv",
    log:
        "results/logs/ancestry/combine_within_ancestry_pcs.log",
    conda:
        "../../envs/gwas.yaml",
    params:
        eigenvecs=lambda wildcards, input: " ".join(
            f"--eigenvec {str(path).split('/')[-1].removesuffix('.eigenvec')}:{path}"
            for path in input.eigenvecs
        ),
    shell:
        """
        Rscript scripts/ancestry_reference.R combine-within-pcs \
          --config {input.config} \
          {params.eigenvecs} \
          --out {output.pcs} \
          > {log} 2>&1
        """
