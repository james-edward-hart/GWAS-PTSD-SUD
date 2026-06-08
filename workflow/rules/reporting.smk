REFERENCE_PREP_REPORT = "results/qc/ancestry/reference/reference_prep_report.md"

# Per-trait reports should wait for the reference-prep report when that branch
# is enabled; otherwise they record "NA" for that link.
reference_prep_report_input = lambda wildcards: [REFERENCE_PREP_REPORT] \
    if ANCESTRY_REFERENCE_ENABLED \
    else []

ancestry_counts_input = lambda wildcards: popmad_counts_file()
POPMAD_PROJECTION_PLOT = "results/plots/ancestry/popmad_reference_study_pcs.png"


rule plot_gwas:
    input:
        stats="results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.plink2.glm.tsv",
    output:
        qq="results/plots/{trait}/{ancestry}/{trait}.{ancestry}.{build}.qq.png",
        manhattan="results/plots/{trait}/{ancestry}/{trait}.{ancestry}.{build}.manhattan.png",
        manhattan_pdf="results/plots/{trait}/{ancestry}/{trait}.{ancestry}.{build}.manhattan.pdf",
    log:
        "results/logs/reporting/plot_gwas.{trait}.{ancestry}.{build}.log",
    conda:
        "../../envs/reporting.yaml",
    shell:
        """
        Rscript scripts/plot_gwas.R \
          --stats {input.stats} \
          --qq {output.qq} \
          --manhattan {output.manhattan} \
          --manhattan-pdf {output.manhattan_pdf} \
          > {log} 2>&1
        """


rule plot_popmad_projection:
    input:
        reference_pcs=ancestry_reference_pcs_file(),
        study_pcs=popmad_study_pcs_file(),
        assignments=popmad_assignments_file(),
        excluded=popmad_excluded_file(),
    output:
        plot=POPMAD_PROJECTION_PLOT,
    log:
        "results/logs/reporting/plot_popmad_projection.log",
    conda:
        "../../envs/reporting.yaml",
    shell:
        """
        Rscript scripts/plot_popmad_projection.R \
          --reference-pcs {input.reference_pcs} \
          --study-pcs {input.study_pcs} \
          --assignments {input.assignments} \
          --excluded {input.excluded} \
          --out {output.plot} \
          > {log} 2>&1
        """


rule make_report:
    input:
        config=RUN_CONFIG,
        stats="results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.plink2.glm.tsv",
        gwas_summary="results/gwas/{trait}/{ancestry}/{trait}.{ancestry}.{build}.gwas_filter_summary.tsv",
        plink_log="results/gwas/{trait}/{ancestry}/plink2_raw/{trait}.{ancestry}.{build}.log",
        qq="results/plots/{trait}/{ancestry}/{trait}.{ancestry}.{build}.qq.png",
        manhattan="results/plots/{trait}/{ancestry}/{trait}.{ancestry}.{build}.manhattan.png",
        manhattan_pdf="results/plots/{trait}/{ancestry}/{trait}.{ancestry}.{build}.manhattan.pdf",
        popmad_plot=POPMAD_PROJECTION_PLOT,
        strata_counts="results/qc/strata/strata_counts.tsv",
        pheno="results/qc/traits/{trait}.pheno.tsv",
        covar="results/qc/traits/{trait}.covar.tsv",
        keep="results/qc/strata/{ancestry}.unrelated.keep.tsv",
        relatedness="results/qc/relatedness/relatedness_summary.tsv",
        sex_check="results/qc/sex/sex_check_summary.tsv",
        genome_details="results/qc/genome_build/genome_build_marker_matches.tsv",
        ancestry_counts=ancestry_counts_input,
        admixture_summary="results/qc/admixture/admixture_run_summary.tsv",
        admixture_study="results/qc/admixture/study_ancestry_proportions.tsv",
        admixture_comparison="results/qc/admixture/popmad_admixture_comparison.tsv",
        admixture_report="results/qc/admixture/admixture_report.md",
        software=lambda wildcards: config["resources"]["software_manifest"],
        reference=lambda wildcards: config["resources"]["reference_manifest"],
        reference_prep_report=reference_prep_report_input,
    output:
        report="results/reports/{trait}/{trait}.{ancestry}.{build}.report.md",
    log:
        "results/logs/reporting/make_report.{trait}.{ancestry}.{build}.log",
    conda:
        "../../envs/gwas.yaml",
    params:
        reference_prep_report=lambda wildcards, input: input.reference_prep_report[0]
        if input.reference_prep_report
        else "NA",
    shell:
        """
        Rscript scripts/make_report.R \
          --config {input.config} \
          --trait {wildcards.trait} \
          --ancestry {wildcards.ancestry} \
          --build {wildcards.build} \
          --stats {input.stats} \
          --gwas-summary {input.gwas_summary} \
          --plink-log {input.plink_log} \
          --qq {input.qq} \
          --manhattan {input.manhattan} \
          --manhattan-pdf {input.manhattan_pdf} \
          --popmad-plot {input.popmad_plot} \
          --strata-counts {input.strata_counts} \
          --pheno {input.pheno} \
          --covar {input.covar} \
          --keep {input.keep} \
          --relatedness-summary {input.relatedness} \
          --sex-check-summary {input.sex_check} \
          --genome-build-details {input.genome_details} \
          --ancestry-counts {input.ancestry_counts} \
          --admixture-summary {input.admixture_summary} \
          --admixture-study {input.admixture_study} \
          --admixture-comparison {input.admixture_comparison} \
          --admixture-report {input.admixture_report} \
          --software {input.software} \
          --reference {input.reference} \
          --reference-prep-report {params.reference_prep_report} \
          --out {output.report} \
          > {log} 2>&1
        """


rule write_run_manifest:
    input:
        config=RUN_CONFIG,
        reports=report_targets,
        build="results/qc/genome_build/genome_build.txt",
    output:
        manifest="results/manifests/run_manifest.tsv",
    log:
        "results/logs/reporting/write_run_manifest.log",
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/write_run_manifest.R \
          --config {input.config} \
          --genome-build-file {input.build} \
          --out {output.manifest} \
          > {log} 2>&1
        """
