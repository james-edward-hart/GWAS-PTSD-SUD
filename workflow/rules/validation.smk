rule write_effective_config:
    input:
        # Treat workflow code and environment YAMLs as provenance inputs.
        # Changing code should refresh the config fingerprint in results/config.
        code=PIPELINE_CODE,
        sources=CONFIG_SOURCE_FILES,
    output:
        config=EFFECTIVE_CONFIG,
    log:
        "results/logs/validation/write_effective_config.log",
    conda:
        "../../envs/gwas.yaml",
    script:
        "../../scripts/write_effective_config.R"


checkpoint infer_genome_build:
    input:
        config=EFFECTIVE_CONFIG,
        markers=lambda wildcards: config["genome_build"]["marker_file"],
    output:
        build="results/qc/genome_build/genome_build.txt",
        details="results/qc/genome_build/genome_build_marker_matches.tsv",
    log:
        "results/logs/validation/infer_genome_build.log",
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/infer_genome_build.R \
          --config {input.config} \
          --out {output.build} \
          --details {output.details} \
          > {log} 2>&1
        """


rule resolve_reference_package:
    input:
        config=EFFECTIVE_CONFIG,
        build="results/qc/genome_build/genome_build.txt",
        # This is the full package file list from file_manifest.tsv, not just
        # the manifest files, so edits to panel data trigger resolution again.
        package=REFERENCE_PACKAGE_INPUTS,
    output:
        config=RUN_CONFIG,
    log:
        "results/logs/validation/resolve_reference_package.log",
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/resolve_reference_package.R \
          --config {input.config} \
          --genome-build-file {input.build} \
          --out {output.config} \
          > {log} 2>&1
        """


rule validate_config:
    input:
        config=RUN_CONFIG,
        traits=lambda wildcards: config["inputs"]["trait_registry"],
        samples=lambda wildcards: config["inputs"]["sample_manifest"],
        input_manifest=lambda wildcards: [config["resources"]["input_manifest"]]
        if config["resources"].get("input_manifest")
        else [],
        software=lambda wildcards: config["resources"]["software_manifest"],
        reference=lambda wildcards: config["resources"]["reference_manifest"],
        markers=lambda wildcards: config["genome_build"]["marker_file"],
        build="results/qc/genome_build/genome_build.txt",
        # Branch inputs are injected only when enabled so Snakemake tracks the
        # concrete reference files selected by package resolution.
        ancestry_reference=ancestry_reference_validation_inputs,
        admixture=admixture_validation_inputs,
    output:
        ok="results/qc/input_validation/validation.ok",
    log:
        "results/logs/validation/validate_config.log",
    conda:
        "../../envs/gwas.yaml",
    shell:
        """
        Rscript scripts/validate_config.R \
          --config {input.config} \
          --genome-build-file {input.build} \
          --out {output.ok} \
          > {log} 2>&1
        """
