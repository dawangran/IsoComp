from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from typer.testing import CliRunner

from isocomp.cli import app


def test_cli_pipeline_smoke(tmp_path: Path, bed12_path: Path, synthetic_bam: Path) -> None:
    out_prefix = tmp_path / "sample.isocomp"

    result = CliRunner().invoke(
        app,
        [
            "--bam",
            str(synthetic_bam),
            "--annotation",
            str(bed12_path),
            "--out",
            str(out_prefix),
            "--bin-num",
            "10",
            "--min-mapq",
            "20",
            "--tss-tol",
            "100",
            "--tes-tol",
            "100",
            "--min-overlap",
            "50",
            "--junction-tol",
            "5",
            "--strandness",
            "unstranded",
            "--threads",
            "1",
        ],
    )
    assert result.exit_code == 0, result.output

    read_assignment = Path(f"{out_prefix}.read_assignment.tsv")
    transcript_metrics = Path(f"{out_prefix}.transcript_metrics.tsv")
    sample_summary = Path(f"{out_prefix}.sample_summary.tsv")
    body_coverage = Path(f"{out_prefix}.transcript_body_coverage.tsv")
    assignment_stats = Path(f"{out_prefix}.assignment_stats.json")
    plots_dir = Path(f"{out_prefix}.plots")

    for path in [read_assignment, transcript_metrics, sample_summary, body_coverage, assignment_stats]:
        assert path.exists()
    for plot_name in [
        "transcript_body_coverage.png",
        "transcript_body_heatmap.png",
        "read_body_heatmap.png",
        "read_coverage_fraction.png",
        "dist_to_5p.png",
        "dist_to_3p.png",
        "full_length_fraction.png",
    ]:
        assert (plots_dir / plot_name).exists()
        assert (plots_dir / plot_name).with_suffix(".pdf").exists()

    reads = pd.read_csv(read_assignment, sep="\t")
    summary = pd.read_csv(sample_summary, sep="\t").iloc[0].to_dict()
    coverage = pd.read_csv(body_coverage, sep="\t")
    stats = json.loads(assignment_stats.read_text(encoding="utf-8"))

    statuses = dict(zip(reads["read_id"], reads["assignment_status"]))
    assert statuses["full_pos"] == "unique"
    assert statuses["trunc_5p"] == "unique"
    assert statuses["trunc_3p"] == "unique"
    assert statuses["ambiguous"] == "ambiguous"
    assert statuses["low_conf"] == "low_confidence"
    assert statuses["unassigned"] == "unassigned"
    assert statuses["full_neg"] == "unique"

    negative = reads.loc[reads["read_id"] == "full_neg"].iloc[0]
    assert negative["dist_to_5p"] == 0
    assert negative["dist_to_3p"] == 0

    assert summary["total_reads"] == 7
    assert summary["unique_assigned_reads"] == 4
    assert summary["ambiguous_reads"] == 1
    assert summary["low_confidence_reads"] == 1
    assert summary["unassigned_reads"] == 1
    assert len(coverage) == 10
    assert list(coverage.columns) == [
        "bin",
        "coverage",
        "mean_normalized_coverage",
        "max_normalized_coverage",
    ]
    assert round(float(coverage["mean_normalized_coverage"].mean()), 6) == 1.0
    assert round(float(coverage["max_normalized_coverage"].max()), 6) == 1.0
    assert stats["parameters"]["bin_num"] == 10
