#!/usr/bin/env python3

"""Focused tests for trait-specific Phase 2 DAG expansion."""

from pathlib import Path
from types import SimpleNamespace
import os
import shutil
import subprocess
import tempfile
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from workflow.snake_helpers import (
    phase2_trait_groups,
    phase2_trait_regenie_outputs,
    remeta_manifest_inputs,
)


def write(path, text):
    path.write_text(text)
    return str(path)


def config_for(samples, traits):
    return {
        "inputs": {"sample_manifest": samples, "trait_registry": traits},
        "phase2_regenie": {
            "enabled": True,
            "global_pcs": 2,
            "default_covariates": ["age", "sex", "PC1", "PC2"],
        },
        "remeta": {"enabled": True},
    }


def checkpoints_for(status):
    result = SimpleNamespace(output=SimpleNamespace(status=status))
    checkpoint = SimpleNamespace(get=lambda: result)
    return SimpleNamespace(select_phase2_active_groups=checkpoint)


with tempfile.TemporaryDirectory(prefix="phase2-dag-") as tmpdir:
    root = Path(tmpdir)
    samples = write(
        root / "samples.tsv",
        "FID\tIID\tage\tsex\tA\tB\tQ\n"
        "F1\tI1\t40\t1\t1\t0\t1.2\n"
        "F2\tI2\t50\t2\t0\tNA\t2.4\n",
    )
    traits = write(
        root / "traits.tsv",
        "trait_id\tphenotype_column\tcase_value\tcontrol_value\tmissing_values\tcovariates\n"
        "A\tA\t1\t0\tNA\t\n"
        "B\tB\t1\t0\tNA\t\n"
        "Q\tQ\t\t\tNA\t\n",
    )
    config = config_for(samples, traits)
    groups = phase2_trait_groups(config)
    observed = {row["traits"][0]: row["group"] for row in groups}
    assert observed == {"A": "bt__A", "B": "bt__B", "Q": "qt__Q"}
    assert all(len(row["traits"]) == 1 for row in groups)

    reordered = write(
        root / "traits_reordered.tsv",
        Path(traits).read_text().splitlines()[0] + "\n"
        + "\n".join(reversed(Path(traits).read_text().splitlines()[1:]))
        + "\n",
    )
    reordered_groups = phase2_trait_groups(config_for(samples, reordered))
    reordered_map = {row["traits"][0]: row["group"] for row in reordered_groups}
    assert reordered_map == observed

    status = write(
        root / "status.tsv",
        "group\ttrait\ttrait_type\tusable_n\tmodel_sample_count\tkeep_count\tmodel_keep_sha256\tskipped\tremeta_eligible\n"
        "bt__A\tA\tbt\t2\t2\t2\thash-a\tFalse\tTrue\n"
        "bt__B\tB\tbt\t1\t0\t0\thash-b\tTrue\tFalse\n"
        "qt__Q\tQ\tqt\t2\t2\t2\thash-q\tFalse\tTrue\n",
    )
    checkpoints = checkpoints_for(status)
    active_wc = SimpleNamespace(trait="A", build="GRCh38")
    skipped_wc = SimpleNamespace(trait="B", build="GRCh38")
    active_outputs = phase2_trait_regenie_outputs(checkpoints, active_wc, config)
    assert len(active_outputs) == 3 and all("bt__A" in path for path in active_outputs)
    assert phase2_trait_regenie_outputs(checkpoints, skipped_wc, config) == []

    artifacts = remeta_manifest_inputs(checkpoints, config, "GRCh38")
    assert len(artifacts) == 2 * (22 * 3 + 1)
    assert not any("bt__B" in path or "/B.PAN" in path for path in artifacts)

    all_skipped = write(
        root / "all_skipped.tsv",
        "group\ttrait\ttrait_type\tusable_n\tmodel_sample_count\tkeep_count\tmodel_keep_sha256\tskipped\tremeta_eligible\n"
        "bt__A\tA\tbt\t2\t0\t0\thash-a\tTrue\tFalse\n"
        "bt__B\tB\tbt\t1\t0\t0\thash-b\tTrue\tFalse\n"
        "qt__Q\tQ\tqt\t2\t0\t0\thash-q\tTrue\tFalse\n",
    )
    skipped_checkpoints = checkpoints_for(all_skipped)
    assert phase2_trait_regenie_outputs(skipped_checkpoints, active_wc, config) == []
    assert remeta_manifest_inputs(skipped_checkpoints, config, "GRCh38") == []

    duplicate_traits = write(
        root / "duplicate_traits.tsv",
        Path(traits).read_text() + "A\tA\t1\t0\tNA\t\n",
    )
    try:
        phase2_trait_groups(config_for(samples, duplicate_traits))
    except Exception:
        pass
    else:
        raise AssertionError("duplicate trait IDs were accepted")

    unsafe_traits = write(
        root / "unsafe_traits.tsv",
        Path(traits).read_text().splitlines()[0] + "\nbad/trait\tA\t1\t0\tNA\t\n",
    )
    try:
        phase2_trait_groups(config_for(samples, unsafe_traits))
    except Exception:
        pass
    else:
        raise AssertionError("path-unsafe trait ID was accepted")

    whitespace_traits = write(
        root / "whitespace_traits.tsv",
        Path(traits).read_text().splitlines()[0] + "\n A\tA\t1\t0\tNA\t\n",
    )
    try:
        phase2_trait_groups(config_for(samples, whitespace_traits))
    except Exception:
        pass
    else:
        raise AssertionError("trait ID with surrounding whitespace was accepted")

    # When run inside the Snakemake environment, exercise real checkpoint
    # reevaluation and dynamic ReMeta expansion in addition to the helper mocks.
    snakemake = shutil.which("snakemake")
    if snakemake:
        repo = Path(__file__).resolve().parents[1]

        def run_checkpoint_case(name, status_rows, expected_active):
            case = root / name
            (case / "fixture").mkdir(parents=True)
            write(
                case / "fixture" / "status.tsv",
                "group\ttrait\ttrait_type\tusable_n\tmodel_sample_count\tkeep_count\tmodel_keep_sha256\tskipped\tremeta_eligible\n"
                + "\n".join(status_rows)
                + "\n",
            )
            snakefile = case / "Snakefile"
            snakefile.write_text(
                f'''from pathlib import Path
import shutil
import sys
sys.path.insert(0, {str(repo)!r})
from workflow.snake_helpers import remeta_manifest_inputs

def dynamic_targets(wildcards):
    return remeta_manifest_inputs(checkpoints, {{"remeta": {{"enabled": True}}}}, "GRCh38")

rule all:
    input:
        dynamic_targets

checkpoint select_phase2_active_groups:
    input:
        source="fixture/status.tsv"
    output:
        status="results/qc/phase2_regenie/trait_group_status.tsv"
    run:
        Path(output.status).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(input.source, output.status)

rule htp:
    output:
        "results/remeta/export/GRCh38/htp/{{trait}}.PAN.regenie.gz"
    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        Path(output[0]).touch()

rule ld:
    output:
        gene="results/remeta/export/GRCh38/ld/{{group}}/chr{{chrom}}.remeta.gene.ld",
        buffer="results/remeta/export/GRCh38/ld/{{group}}/chr{{chrom}}.remeta.buffer.ld",
        idx="results/remeta/export/GRCh38/ld/{{group}}/chr{{chrom}}.remeta.ld.idx.gz"
    wildcard_constraints:
        chrom="[0-9]+"
    run:
        for path in output:
            Path(path).parent.mkdir(parents=True, exist_ok=True)
            Path(path).touch()
'''
            )
            environment = dict(os.environ)
            environment["XDG_CACHE_HOME"] = str(case / ".cache")
            completed = subprocess.run(
                [
                    snakemake,
                    "--snakefile", str(snakefile),
                    "--directory", str(case),
                    "--cores", "1",
                    "--force-use-threads",
                    "all",
                    # Avoid the user-level source cache so this fixture remains
                    # self-contained in restricted development environments.
                    "--shared-fs-usage",
                    "input-output", "persistence", "software-deployment",
                    "software-deployment-cache", "sources", "storage-local-copies",
                ],
                text=True,
                capture_output=True,
                env=environment,
            )
            if completed.returncode:
                raise AssertionError(
                    "checkpoint integration failed:\n" + completed.stdout + "\n" + completed.stderr
                )
            assert (case / "results/qc/phase2_regenie/trait_group_status.tsv").exists()
            htp = list((case / "results/remeta/export/GRCh38/htp").glob("*.gz"))
            ld = list((case / "results/remeta/export/GRCh38/ld").glob("*/chr*.remeta.*"))
            assert len(htp) == expected_active
            assert len(ld) == expected_active * 22 * 3
            assert not any("TUD" in str(path) for path in htp + ld)

        active_rows = [
            "bt__ANY\tANY\tbt\t625\t625\t625\th1\tFalse\tTrue",
            "bt__AUD\tAUD\tbt\t885\t885\t885\th2\tFalse\tTrue",
            "bt__CUD\tCUD\tbt\t768\t768\t768\th3\tFalse\tTrue",
            "bt__OUD\tOUD\tbt\t1112\t1112\t1112\th4\tFalse\tTrue",
            "bt__TUD\tTUD\tbt\t416\t0\t0\th5\tTrue\tFalse",
        ]
        run_checkpoint_case("checkpoint_mixed", active_rows, 4)
        all_skipped_rows = [
            "bt__ANY\tANY\tbt\t625\t0\t0\th1\tTrue\tFalse",
            "bt__AUD\tAUD\tbt\t885\t0\t0\th2\tTrue\tFalse",
            "bt__CUD\tCUD\tbt\t768\t0\t0\th3\tTrue\tFalse",
            "bt__OUD\tOUD\tbt\t1112\t0\t0\th4\tTrue\tFalse",
            "bt__TUD\tTUD\tbt\t416\t0\t0\th5\tTrue\tFalse",
        ]
        run_checkpoint_case(
            "checkpoint_all_skipped",
            all_skipped_rows,
            0,
        )

print("Phase 2 DAG helper tests passed")
