rule infer_popmad_ancestry:
    input:
        ok="results/qc/input_validation/validation.ok",
        study=ancestry_study_pcs_file,
        reference=ancestry_reference_pcs_file,
    output:
        assignments=popmad_assignments_file(),
        pcs=popmad_study_pcs_file(),
        distances=popmad_distances_file(),
        reference_outliers=popmad_reference_outliers_file(),
        model_summary=popmad_model_summary_file(),
        within=popmad_within_file(),
        excluded=popmad_excluded_file(),
        counts=popmad_counts_file(),
    log:
        "results/logs/ancestry/infer_popmad_ancestry.log",
    params:
        ancestries=",".join(ANCESTRIES),
        pcs=config["popmad"]["pcs"],
        min_confidence=config["popmad"]["min_confidence"],
        outlier_sd=config["popmad"]["reference_outlier_sd"],
        min_reference_n=config["popmad"].get("min_reference_population_n", 20),
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/infer_popmad_ancestry.R \
          --study-pcs {input.study} \
          --reference-pcs {input.reference} \
          --ancestries {params.ancestries} \
          --pcs {params.pcs} \
          --assignments {output.assignments} \
          --study-pcs-out {output.pcs} \
          --distances {output.distances} \
          --reference-outliers {output.reference_outliers} \
          --model-summary {output.model_summary} \
          --within-pcs {output.within} \
          --excluded {output.excluded} \
          --counts {output.counts} \
          --min-confidence {params.min_confidence} \
          --outlier-sd {params.outlier_sd} \
          --min-reference-n {params.min_reference_n} \
          > {log} 2>&1
        """
