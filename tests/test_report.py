from __future__ import annotations

import base64
from pathlib import Path

from isocomp.report import write_html_report, write_pdf_report


PNG_1X1 = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+/p9sAAAAASUVORK5CYII="
)


def _sample_row() -> dict[str, object]:
    return {
        "sample": "sample <A>",
        "total_reads": 10,
        "usable_reads": 8,
        "assigned_reads": 6,
        "unique_assigned_reads": 5,
        "ambiguous_reads": 1,
        "low_confidence_reads": 1,
        "unassigned_reads": 1,
        "terminal_evaluable_unique_reads": 5,
        "median_transcript_coverage_fraction": 0.75,
        "mean_transcript_coverage_fraction": 0.7,
        "5p_complete_fraction": 0.6,
        "3p_complete_fraction": 0.8,
        "full_length_like_fraction": 0.5,
        "median_dist_to_5p": 20,
        "median_dist_to_3p": 0.0,
    }


def _stats_payload() -> dict[str, object]:
    return {
        "version": "0.1.0",
        "parameters": {
            "bam": "reads.bam",
            "annotation": "tx.bed12",
            "bin_num": 100,
            "resolved_annotation_format": "bed12",
            "resolved_strandness": "unstranded",
        },
    }


def test_write_html_report_embeds_plots_and_escapes_text(tmp_path: Path) -> None:
    plots_dir = tmp_path / "sample.plots"
    plots_dir.mkdir()
    (plots_dir / "transcript_body_coverage.png").write_bytes(PNG_1X1)

    output_path = tmp_path / "sample.report.html"
    write_html_report(
        output_path,
        sample_row=_sample_row(),
        stats_payload=_stats_payload(),
        plots_dir=plots_dir,
        output_paths={"read_assignment": tmp_path / "sample.read_assignment.tsv"},
    )

    html = output_path.read_text(encoding="utf-8")
    assert "IsoComp QC Report" in html
    assert "sample &lt;A&gt;" in html
    assert "data:image/png;base64," in html
    assert "Full-length-like" in html
    assert "sample.read_assignment.tsv" in html
    assert "<td>Median dist to 3&#x27;</td><td>0</td>" in html


def test_write_pdf_report_creates_valid_pdf_with_missing_plots(tmp_path: Path) -> None:
    output_path = tmp_path / "sample.report.pdf"

    write_pdf_report(
        output_path,
        sample_row=_sample_row(),
        stats_payload=_stats_payload(),
        plots_dir=tmp_path / "missing.plots",
    )

    assert output_path.read_bytes().startswith(b"%PDF")
    assert output_path.stat().st_size > 1000
