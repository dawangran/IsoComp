"""Static HTML and PDF report generation for IsoComp runs."""

from __future__ import annotations

import base64
import html
import math
import os
import tempfile
from collections.abc import Mapping, Sequence
from datetime import datetime
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any

_MPLCONFIGDIR = Path(tempfile.gettempdir()) / "isocomp-matplotlib"
_MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPLCONFIGDIR))
_XDG_CACHE_HOME = Path(tempfile.gettempdir()) / "isocomp-cache"
_XDG_CACHE_HOME.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("XDG_CACHE_HOME", str(_XDG_CACHE_HOME))

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


PLOT_SPECS = [
    (
        "transcript_body_coverage.png",
        "Transcript Body Coverage",
        "Aggregate unique assigned-read coverage across normalized transcript coordinates.",
    ),
    (
        "transcript_body_heatmap.png",
        "Transcript Body Heatmap",
        "Mean-normalized transcript-by-bin coverage for the most-covered transcripts.",
    ),
    (
        "transcript_body_heatmap_full.png",
        "Full Transcript Body Heatmap",
        "Mean-normalized transcript-by-bin coverage without the display row cap.",
    ),
    (
        "read_body_heatmap.png",
        "Read Body Heatmap",
        "Unique reads projected into assigned transcript coordinates.",
    ),
    (
        "read_coverage_fraction.png",
        "Read Coverage Fraction",
        "Distribution of assigned transcript fractions covered by unique reads.",
    ),
    (
        "dist_to_5p.png",
        "Distance To 5' End",
        "Projected distance from read start to the annotated transcript 5' end.",
    ),
    (
        "dist_to_3p.png",
        "Distance To 3' End",
        "Projected distance from read end to the annotated transcript 3' end.",
    ),
    (
        "full_length_fraction.png",
        "Terminal Completeness",
        "5' complete, 3' complete, and full-length-like fractions.",
    ),
]

KPI_SPECS = [
    ("Unique reads", "unique_assigned_reads", "Unique transcript assignments"),
    ("Full-length-like", "full_length_like_fraction", "Terminal-evaluable unique reads"),
    ("5' complete", "5p_complete_fraction", "Terminal-evaluable unique reads"),
    ("3' complete", "3p_complete_fraction", "Terminal-evaluable unique reads"),
    ("Median coverage", "median_transcript_coverage_fraction", "Assigned transcript fraction"),
    ("Median 5' distance", "median_dist_to_5p", "bp"),
    ("Median 3' distance", "median_dist_to_3p", "bp"),
]

SUMMARY_KEYS = [
    "total_reads",
    "mapped_reads",
    "primary_reads",
    "usable_reads",
    "assigned_reads",
    "unique_assigned_reads",
    "ambiguous_reads",
    "low_confidence_reads",
    "unassigned_reads",
    "terminal_evaluable_unique_reads",
    "median_read_aligned_length",
    "median_transcript_coverage_fraction",
    "mean_transcript_coverage_fraction",
    "5p_complete_fraction",
    "3p_complete_fraction",
    "full_length_like_fraction",
    "median_dist_to_5p",
    "median_dist_to_3p",
    "all_assigned_median_transcript_coverage_fraction",
    "all_assigned_5p_complete_fraction",
    "all_assigned_3p_complete_fraction",
    "all_assigned_full_length_like_fraction",
]

PARAMETER_KEYS = [
    "bam",
    "annotation",
    "annotation_format",
    "resolved_annotation_format",
    "out",
    "bin_num",
    "min_mapq",
    "tss_tol",
    "tes_tol",
    "min_overlap",
    "junction_tol",
    "unique_threshold",
    "margin_threshold",
    "full_length_coverage",
    "min_terminal_anchor",
    "min_unspliced_coverage_for_unique",
    "strandness",
    "resolved_strandness",
    "threads",
]


def write_reports(
    *,
    html_path: Path,
    pdf_path: Path,
    sample_row: Mapping[str, Any],
    stats_payload: Mapping[str, Any],
    plots_dir: Path,
    output_paths: Mapping[str, Path],
) -> None:
    write_html_report(
        html_path,
        sample_row=sample_row,
        stats_payload=stats_payload,
        plots_dir=plots_dir,
        output_paths=output_paths,
    )
    write_pdf_report(
        pdf_path,
        sample_row=sample_row,
        stats_payload=stats_payload,
        plots_dir=plots_dir,
    )


def write_html_report(
    path: Path,
    *,
    sample_row: Mapping[str, Any],
    stats_payload: Mapping[str, Any],
    plots_dir: Path,
    output_paths: Mapping[str, Path] | None = None,
) -> None:
    sample_name = _string_value(sample_row.get("sample", "IsoComp sample"))
    generated_at = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")
    version = _string_value(stats_payload.get("version", "unknown"))
    parameters = _mapping_value(stats_payload.get("parameters"))

    html_text = "\n".join(
        [
            "<!doctype html>",
            '<html lang="en">',
            "<head>",
            '<meta charset="utf-8">',
            '<meta name="viewport" content="width=device-width, initial-scale=1">',
            f"<title>IsoComp QC Report - {_escape(sample_name)}</title>",
            "<style>",
            _html_css(),
            "</style>",
            "</head>",
            "<body>",
            '<main class="page">',
            '<section class="hero">',
            "<div>",
            "<p>IsoComp</p>",
            "<h1>IsoComp QC Report</h1>",
            f"<h2>{_escape(sample_name)}</h2>",
            "</div>",
            '<dl class="run-meta">',
            f"<div><dt>Version</dt><dd>{_escape(version)}</dd></div>",
            f"<div><dt>Generated</dt><dd>{_escape(generated_at)}</dd></div>",
            f"<div><dt>Annotation</dt><dd>{_escape(_string_value(parameters.get('resolved_annotation_format', 'unknown')))}</dd></div>",
            f"<div><dt>Strandness</dt><dd>{_escape(_string_value(parameters.get('resolved_strandness', 'unknown')))}</dd></div>",
            "</dl>",
            "</section>",
            _kpi_section(sample_row),
            _summary_section(sample_row),
            _plots_section(plots_dir),
            _parameters_section(parameters),
            _output_files_section(output_paths or {}),
            "</main>",
            "</body>",
            "</html>",
        ]
    )
    _atomic_write_text(path, html_text)


def write_pdf_report(
    path: Path,
    *,
    sample_row: Mapping[str, Any],
    stats_payload: Mapping[str, Any],
    plots_dir: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("wb", dir=path.parent, delete=False) as handle:
        tmp_name = handle.name
    try:
        with PdfPages(tmp_name) as pdf:
            _write_pdf_summary_page(pdf, sample_row, stats_payload)
            _write_pdf_plot_pages(pdf, plots_dir)
        os.replace(tmp_name, path)
    except Exception:
        try:
            os.unlink(tmp_name)
        except FileNotFoundError:
            pass
        raise


def _write_pdf_summary_page(
    pdf: PdfPages,
    sample_row: Mapping[str, Any],
    stats_payload: Mapping[str, Any],
) -> None:
    sample_name = _string_value(sample_row.get("sample", "IsoComp sample"))
    version = _string_value(stats_payload.get("version", "unknown"))
    parameters = _mapping_value(stats_payload.get("parameters"))

    fig = plt.figure(figsize=(8.27, 11.69))
    fig.patch.set_facecolor("white")
    fig.text(0.08, 0.95, "IsoComp QC Report", fontsize=22, weight="bold")
    fig.text(0.08, 0.92, sample_name, fontsize=13, color="#333333")
    fig.text(
        0.08,
        0.89,
        (
            f"Version {version} | annotation "
            f"{_string_value(parameters.get('resolved_annotation_format', 'unknown'))} | "
            f"strandness {_string_value(parameters.get('resolved_strandness', 'unknown'))}"
        ),
        fontsize=9,
        color="#555555",
    )
    fig.text(
        0.08,
        0.86,
        "Read-centric transcript completeness summary from unique assignments.",
        fontsize=10,
        color="#333333",
    )

    _draw_pdf_table(
        fig,
        [0.08, 0.58, 0.84, 0.22],
        [
            [label, _format_field_value(key, sample_row.get(key)), description]
            for label, key, description in KPI_SPECS
        ],
        ["Metric", "Value", "Definition"],
        font_size=8,
    )
    _draw_pdf_table(
        fig,
        [0.08, 0.18, 0.84, 0.32],
        [
            [_label_from_key(key), _format_field_value(key, sample_row.get(key))]
            for key in SUMMARY_KEYS
            if key in sample_row
        ],
        ["Sample field", "Value"],
        font_size=7,
    )
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def _write_pdf_plot_pages(pdf: PdfPages, plots_dir: Path) -> None:
    for page_index in range(0, len(PLOT_SPECS), 4):
        page_specs = PLOT_SPECS[page_index : page_index + 4]
        fig, axes = plt.subplots(2, 2, figsize=(8.27, 11.69))
        fig.patch.set_facecolor("white")
        fig.suptitle("IsoComp Plots", fontsize=16, weight="bold", y=0.98)
        for ax, (filename, title, caption) in zip(axes.ravel(), page_specs):
            plot_path = plots_dir / filename
            ax.set_title(title, fontsize=10, pad=8)
            ax.set_axis_off()
            if plot_path.exists():
                image = mpimg.imread(plot_path)
                ax.imshow(image, aspect="auto")
            else:
                ax.text(
                    0.5,
                    0.56,
                    "Plot not generated",
                    ha="center",
                    va="center",
                    fontsize=10,
                    color="#555555",
                    transform=ax.transAxes,
                )
            ax.text(
                0.0,
                -0.05,
                caption,
                ha="left",
                va="top",
                fontsize=7,
                color="#333333",
                wrap=True,
                transform=ax.transAxes,
            )
        for ax in axes.ravel()[len(page_specs) :]:
            ax.set_axis_off()
        fig.tight_layout(rect=(0.04, 0.04, 0.96, 0.95))
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)


def _draw_pdf_table(
    fig: plt.Figure,
    bbox: list[float],
    rows: Sequence[Sequence[str]],
    columns: Sequence[str],
    *,
    font_size: int,
) -> None:
    ax = fig.add_axes(bbox)
    ax.set_axis_off()
    table = ax.table(
        cellText=rows,
        colLabels=columns,
        loc="center",
        cellLoc="left",
        colLoc="left",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(font_size)
    table.scale(1, 1.25)
    for (row_index, _), cell in table.get_celld().items():
        cell.set_edgecolor("#D8D8D8")
        cell.set_linewidth(0.4)
        if row_index == 0:
            cell.set_facecolor("#F0F3F6")
            cell.set_text_props(weight="bold", color="#222222")
        else:
            cell.set_facecolor("white")


def _kpi_section(sample_row: Mapping[str, Any]) -> str:
    cards = []
    for label, key, description in KPI_SPECS:
        cards.append(
            "\n".join(
                [
                    '<article class="kpi-card">',
                    f"<p>{_escape(description)}</p>",
                    f"<strong>{_escape(_format_field_value(key, sample_row.get(key)))}</strong>",
                    f"<span>{_escape(label)}</span>",
                    "</article>",
                ]
            )
        )
    return f'<section class="kpi-grid">{"".join(cards)}</section>'


def _summary_section(sample_row: Mapping[str, Any]) -> str:
    rows = [
        (key, sample_row.get(key))
        for key in SUMMARY_KEYS
        if key in sample_row
    ]
    return "\n".join(
        [
            '<section class="panel">',
            "<h2>Sample Summary</h2>",
            _html_table(
                [("Metric", "Value")]
                + [(_label_from_key(key), _format_field_value(key, value)) for key, value in rows]
            ),
            "</section>",
        ]
    )


def _plots_section(plots_dir: Path) -> str:
    figures = []
    for filename, title, caption in PLOT_SPECS:
        plot_path = plots_dir / filename
        if plot_path.exists():
            image_html = (
                f'<img src="{_escape(_image_data_uri(plot_path))}" '
                f'alt="{_escape(title)}">'
            )
        else:
            image_html = '<div class="missing-plot">Plot not generated</div>'
        figures.append(
            "\n".join(
                [
                    "<figure>",
                    image_html,
                    f"<figcaption><strong>{_escape(title)}</strong>{_escape(caption)}</figcaption>",
                    "</figure>",
                ]
            )
        )
    return "\n".join(
        [
            '<section class="panel">',
            "<h2>Plots</h2>",
            '<div class="plot-grid">',
            "".join(figures),
            "</div>",
            "</section>",
        ]
    )


def _parameters_section(parameters: Mapping[str, Any]) -> str:
    rows = [
        (key, parameters.get(key))
        for key in PARAMETER_KEYS
        if key in parameters
    ]
    return "\n".join(
        [
            '<section class="panel">',
            "<h2>Run Parameters</h2>",
            _html_table(
                [("Parameter", "Value")]
                + [(_label_from_key(key), _format_field_value(key, value)) for key, value in rows]
            ),
            "</section>",
        ]
    )


def _output_files_section(output_paths: Mapping[str, Path]) -> str:
    if not output_paths:
        return ""
    rows = [
        (_label_from_key(key), str(path))
        for key, path in output_paths.items()
        if key != "plots_dir"
    ]
    if "plots_dir" in output_paths:
        rows.append(("Plots directory", str(output_paths["plots_dir"])))
    return "\n".join(
        [
            '<section class="panel">',
            "<h2>Output Files</h2>",
            _html_table([("Output", "Path")] + rows),
            "</section>",
        ]
    )


def _html_table(rows: Sequence[Sequence[str]]) -> str:
    if not rows:
        return ""
    header = rows[0]
    body = rows[1:]
    header_html = "".join(f"<th>{_escape(str(item))}</th>" for item in header)
    body_html = "\n".join(
        "<tr>" + "".join(f"<td>{_escape(str(item))}</td>" for item in row) + "</tr>"
        for row in body
    )
    return "\n".join(
        [
            "<table>",
            f"<thead><tr>{header_html}</tr></thead>",
            f"<tbody>{body_html}</tbody>",
            "</table>",
        ]
    )


def _html_css() -> str:
    return """
:root {
  color-scheme: light;
  --ink: #1f2933;
  --muted: #596773;
  --line: #d8dee6;
  --surface: #ffffff;
  --band: #f4f7fa;
  --accent: #1f77b4;
  --accent-2: #009e73;
}
* { box-sizing: border-box; }
body {
  margin: 0;
  background: var(--band);
  color: var(--ink);
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  font-size: 14px;
  line-height: 1.45;
}
.page {
  max-width: 1180px;
  margin: 0 auto;
  padding: 32px 24px 48px;
}
.hero {
  display: flex;
  gap: 24px;
  align-items: flex-end;
  justify-content: space-between;
  padding: 28px 0 22px;
  border-bottom: 1px solid var(--line);
}
.hero p {
  margin: 0 0 8px;
  color: var(--accent);
  font-size: 13px;
  font-weight: 700;
  letter-spacing: 0.04em;
  text-transform: uppercase;
}
h1, h2 {
  margin: 0;
  line-height: 1.15;
}
h1 { font-size: 34px; }
.hero h2 {
  margin-top: 8px;
  color: var(--muted);
  font-size: 17px;
  font-weight: 500;
}
.run-meta {
  display: grid;
  grid-template-columns: repeat(2, minmax(140px, 1fr));
  gap: 12px 20px;
  min-width: min(420px, 100%);
  margin: 0;
}
.run-meta dt {
  color: var(--muted);
  font-size: 12px;
}
.run-meta dd {
  margin: 2px 0 0;
  font-weight: 650;
  overflow-wrap: anywhere;
}
.kpi-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
  gap: 12px;
  margin: 22px 0;
}
.kpi-card {
  min-height: 112px;
  padding: 16px;
  background: var(--surface);
  border: 1px solid var(--line);
  border-radius: 8px;
}
.kpi-card p {
  min-height: 34px;
  margin: 0;
  color: var(--muted);
  font-size: 12px;
}
.kpi-card strong {
  display: block;
  margin-top: 12px;
  font-size: 27px;
  line-height: 1;
}
.kpi-card span {
  display: block;
  margin-top: 7px;
  color: var(--muted);
  font-size: 13px;
}
.panel {
  margin-top: 22px;
  padding: 22px;
  background: var(--surface);
  border: 1px solid var(--line);
  border-radius: 8px;
}
.panel h2 {
  margin-bottom: 16px;
  font-size: 20px;
}
table {
  width: 100%;
  border-collapse: collapse;
  table-layout: fixed;
}
th, td {
  padding: 8px 9px;
  border-bottom: 1px solid var(--line);
  text-align: left;
  vertical-align: top;
  overflow-wrap: anywhere;
}
th {
  color: var(--muted);
  font-size: 12px;
  text-transform: uppercase;
}
.plot-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
  gap: 18px;
}
figure {
  margin: 0;
  padding: 12px;
  border: 1px solid var(--line);
  border-radius: 8px;
}
figure img {
  display: block;
  width: 100%;
  height: auto;
  background: white;
}
figcaption {
  margin-top: 10px;
  color: var(--muted);
  font-size: 12px;
}
figcaption strong {
  display: block;
  margin-bottom: 2px;
  color: var(--ink);
  font-size: 13px;
}
.missing-plot {
  display: grid;
  min-height: 220px;
  place-items: center;
  color: var(--muted);
  background: var(--band);
  border: 1px dashed var(--line);
  border-radius: 6px;
}
@media print {
  body { background: white; }
  .page { max-width: none; padding: 0; }
  .panel, .kpi-card, figure { break-inside: avoid; }
}
@media (max-width: 760px) {
  .hero { display: block; }
  .run-meta { margin-top: 18px; grid-template-columns: 1fr; }
  h1 { font-size: 28px; }
  .page { padding: 20px 14px 32px; }
}
"""


def _image_data_uri(path: Path) -> str:
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def _format_field_value(key: str, value: Any) -> str:
    return _format_value(value, as_fraction=_is_fraction_key(key))


def _format_value(value: Any, *, as_fraction: bool = False) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return "NA"
        if as_fraction and 0 <= value <= 1:
            return f"{value * 100:.1f}%"
        if abs(value) >= 1000:
            return f"{value:,.0f}"
        return f"{value:.3g}"
    return str(value)


def _is_fraction_key(key: str) -> bool:
    return (
        key.endswith("_fraction")
        or key.endswith("_pct")
        or "coverage_fraction" in key
        or key in {
            "full_length_coverage",
            "margin_threshold",
            "min_unspliced_coverage_for_unique",
            "unique_threshold",
        }
    )


def _label_from_key(key: str) -> str:
    label = key.replace("_", " ")
    label = label.replace("5p", "5'").replace("3p", "3'")
    return label[:1].upper() + label[1:]


def _mapping_value(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _string_value(value: Any) -> str:
    if value is None:
        return "NA"
    return str(value)


def _escape(value: str) -> str:
    return html.escape(value, quote=True)


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("wt", encoding="utf-8", dir=path.parent, delete=False) as handle:
        tmp_name = handle.name
        handle.write(text)
    os.replace(tmp_name, path)
