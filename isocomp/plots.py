"""Matplotlib plot generation."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

_MPLCONFIGDIR = Path(tempfile.gettempdir()) / "isocomp-matplotlib"
_MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPLCONFIGDIR))
_XDG_CACHE_HOME = Path(tempfile.gettempdir()) / "isocomp-cache"
_XDG_CACHE_HOME.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("XDG_CACHE_HOME", str(_XDG_CACHE_HOME))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .models import AssignmentResult, ReadMetrics

JOURNAL_DPI = 300
JOURNAL_RC_PARAMS = {
    "figure.dpi": 150,
    "savefig.dpi": JOURNAL_DPI,
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.unicode_minus": False,
}
COLORS = {
    "blue": "#0072B2",
    "green": "#009E73",
    "vermillion": "#D55E00",
    "axis": "#222222",
    "grid": "#D9D9D9",
}


def write_plots(
    plots_dir: Path,
    read_metrics: list[ReadMetrics],
    body_coverage: np.ndarray,
    per_transcript_coverage: dict[str, np.ndarray] | None = None,
    assignments: list[AssignmentResult] | None = None,
) -> None:
    plots_dir.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(JOURNAL_RC_PARAMS):
        _plot_body_coverage(plots_dir / "transcript_body_coverage.png", body_coverage)
        _plot_transcript_body_heatmap(
            plots_dir / "transcript_body_heatmap.png",
            per_transcript_coverage or {},
        )
        _plot_read_body_heatmap(
            plots_dir / "read_body_heatmap.png",
            assignments or [],
            bin_num=len(body_coverage),
        )
        _plot_hist(
            plots_dir / "read_coverage_fraction.png",
            [
                item.coverage_fraction
                for item in read_metrics
                if item.assignment_status == "unique" and item.coverage_fraction is not None
            ],
            "Transcript coverage fraction",
            "Read count",
            bins=30,
            range_=(0, 1),
        )
        _plot_hist(
            plots_dir / "dist_to_5p.png",
            [
                item.dist_to_5p
                for item in read_metrics
                if item.assignment_status == "unique" and item.dist_to_5p is not None
            ],
            "Distance to 5' end (bp)",
            "Read count",
            bins=30,
        )
        _plot_hist(
            plots_dir / "dist_to_3p.png",
            [
                item.dist_to_3p
                for item in read_metrics
                if item.assignment_status == "unique" and item.dist_to_3p is not None
            ],
            "Distance to 3' end (bp)",
            "Read count",
            bins=30,
        )
        _plot_full_length_fraction(plots_dir / "full_length_fraction.png", read_metrics)


def _plot_body_coverage(path: Path, coverage: np.ndarray) -> None:
    normalized = _max_normalize_coverage(coverage)
    fig, ax = plt.subplots(figsize=(3.6, 2.6))
    x = _bin_centers(len(normalized))
    if x.size:
        ax.fill_between(x, normalized, color=COLORS["blue"], alpha=0.14, linewidth=0)
        ax.plot(x, normalized, color=COLORS["blue"], linewidth=1.7)
    ax.set_xlabel("Transcript position (5' to 3', %)")
    ax.set_ylabel("Max-normalized coverage")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 1.05)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_yticks([0, 0.25, 0.50, 0.75, 1.00])
    _style_axis(ax, grid_axis="y")
    _save_figure(fig, path)


def _plot_transcript_body_heatmap(
    path: Path,
    per_transcript_coverage: dict[str, np.ndarray],
    *,
    max_rows: int = 200,
) -> None:
    covered = [
        (transcript_id, values)
        for transcript_id, values in per_transcript_coverage.items()
        if values.size and float(np.sum(values)) > 0
    ]
    covered.sort(key=lambda item: (-float(np.sum(item[1])), item[0]))
    covered = covered[:max_rows]
    bin_num = _infer_bin_num(per_transcript_coverage)

    fig, ax = plt.subplots(figsize=(4.8, max(2.8, min(7.2, 0.07 * len(covered) + 2.2))))
    if covered:
        matrix = np.vstack([_mean_normalize_coverage(values) for _, values in covered])
        image = ax.imshow(
            matrix,
            aspect="auto",
            interpolation="nearest",
            cmap="cividis",
            vmin=0,
            vmax=_robust_vmax(matrix, default=1.0),
        )
        image.set_rasterized(True)
        colorbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
        colorbar.set_label("Mean-normalized coverage")
        _style_colorbar(colorbar)
        ax.set_yticks(_heatmap_yticks(len(covered)))
        ax.set_yticklabels([covered[index][0] for index in _heatmap_yticks(len(covered))], fontsize=7)
    else:
        image = ax.imshow(
            np.zeros((1, max(1, bin_num))),
            aspect="auto",
            interpolation="nearest",
            cmap="cividis",
            vmin=0,
            vmax=1,
        )
        image.set_rasterized(True)
        ax.set_yticks([])
        _empty_panel(ax, "No uniquely assigned transcripts")

    ax.set_xlabel("Transcript position (5' to 3', %)")
    ax.set_ylabel("Transcript")
    ax.set_xticks(_bin_axis_ticks(bin_num))
    ax.set_xticklabels(_bin_axis_tick_labels(bin_num))
    _style_axis(ax)
    _save_figure(fig, path)


def _plot_read_body_heatmap(
    path: Path,
    assignments: list[AssignmentResult],
    *,
    bin_num: int,
    max_rows: int = 500,
) -> None:
    rows: list[tuple[str, str, int, int, np.ndarray]] = []
    for assignment in assignments:
        if assignment.status != "unique" or assignment.transcript is None or assignment.projection is None:
            continue
        projection = assignment.projection
        if projection.read_start_tx is None or projection.read_end_tx is None:
            continue
        row = np.zeros(bin_num, dtype=float)
        for interval in projection.intervals:
            _add_interval_to_bins(row, interval, assignment.transcript.transcript_length)
        rows.append(
            (
                assignment.read_id,
                assignment.transcript.transcript_id,
                projection.read_start_tx,
                projection.read_end_tx,
                row,
            )
        )

    rows.sort(key=lambda item: (item[2], item[3], item[1], item[0]))
    rows = rows[:max_rows]

    fig, ax = plt.subplots(figsize=(4.8, max(2.8, min(7.2, 0.035 * len(rows) + 2.2))))
    if rows:
        matrix = np.vstack([item[4] for item in rows])
        image = ax.imshow(
            matrix,
            aspect="auto",
            interpolation="nearest",
            cmap="Greys",
            vmin=0,
            vmax=1,
        )
        image.set_rasterized(True)
        colorbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
        colorbar.set_label("Covered fraction of bin")
        _style_colorbar(colorbar)
        tick_indexes = _heatmap_yticks(len(rows))
        ax.set_yticks(tick_indexes)
        ax.set_yticklabels([rows[index][0] for index in tick_indexes], fontsize=7)
    else:
        image = ax.imshow(
            np.zeros((1, max(1, bin_num))),
            aspect="auto",
            interpolation="nearest",
            cmap="Greys",
            vmin=0,
            vmax=1,
        )
        image.set_rasterized(True)
        ax.set_yticks([])
        _empty_panel(ax, "No uniquely assigned reads")

    ax.set_xlabel("Transcript position (5' to 3', %)")
    ax.set_ylabel("Read")
    ax.set_xticks(_bin_axis_ticks(bin_num))
    ax.set_xticklabels(_bin_axis_tick_labels(bin_num))
    _style_axis(ax)
    _save_figure(fig, path)


def _plot_hist(
    path: Path,
    values: list[float | int],
    xlabel: str,
    ylabel: str,
    *,
    bins: int,
    range_: tuple[float, float] | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    if values:
        ax.hist(
            values,
            bins=bins,
            range=range_,
            color=COLORS["blue"],
            edgecolor="white",
            linewidth=0.5,
        )
        _add_median_marker(ax, values)
    else:
        ax.hist([], bins=bins, range=range_, color=COLORS["blue"], edgecolor="white")
        _empty_panel(ax, "No uniquely assigned reads")
    if range_ is not None:
        ax.set_xlim(*range_)
    else:
        ax.set_xlim(left=0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    _style_axis(ax, grid_axis="y")
    _save_figure(fig, path)


def _plot_full_length_fraction(path: Path, read_metrics: list[ReadMetrics]) -> None:
    unique = [item for item in read_metrics if item.assignment_status == "unique"]
    fractions = [
        sum(1 for item in unique if item.is_5p_complete) / len(unique) if unique else 0.0,
        sum(1 for item in unique if item.is_3p_complete) / len(unique) if unique else 0.0,
        sum(1 for item in unique if item.is_full_length_like) / len(unique) if unique else 0.0,
    ]
    labels = ["5' complete", "3' complete", "Full-length\nlike"]
    colors = [COLORS["blue"], COLORS["green"], COLORS["vermillion"]]
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    x = np.arange(len(fractions))
    ax.bar(x, fractions, color=colors, edgecolor="white", linewidth=0.6, width=0.68)
    for index, value in enumerate(fractions):
        ax.text(index, min(value + 0.035, 1.02), f"{value:.2f}", ha="center", va="bottom", fontsize=7)
    if not unique:
        _empty_panel(ax, "No uniquely assigned reads", y=0.58)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of unique reads")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    _style_axis(ax, grid_axis="y")
    _save_figure(fig, path)


def _style_axis(ax: plt.Axes, *, grid_axis: str | None = None) -> None:
    for spine_name in ("top", "right"):
        ax.spines[spine_name].set_visible(False)
    for spine_name in ("left", "bottom"):
        ax.spines[spine_name].set_color(COLORS["axis"])
        ax.spines[spine_name].set_linewidth(0.8)
    ax.tick_params(axis="both", direction="out", length=3, width=0.8, colors=COLORS["axis"])
    ax.xaxis.label.set_color(COLORS["axis"])
    ax.yaxis.label.set_color(COLORS["axis"])
    if grid_axis is not None:
        ax.grid(axis=grid_axis, color=COLORS["grid"], linewidth=0.5)
        ax.set_axisbelow(True)


def _style_colorbar(colorbar: object) -> None:
    colorbar.outline.set_linewidth(0.6)
    colorbar.ax.tick_params(length=2.5, width=0.6)


def _empty_panel(ax: plt.Axes, message: str, *, y: float = 0.5) -> None:
    ax.text(
        0.5,
        y,
        message,
        ha="center",
        va="center",
        color="#555555",
        transform=ax.transAxes,
    )


def _add_median_marker(ax: plt.Axes, values: list[float | int]) -> None:
    median_value = float(np.median(values))
    ax.axvline(
        median_value,
        color=COLORS["vermillion"],
        linestyle=(0, (3, 2)),
        linewidth=1.0,
        label="Median",
    )
    ax.legend(frameon=False, handlelength=1.6, loc="upper right")


def _save_figure(fig: plt.Figure, path: Path) -> None:
    fig.tight_layout(pad=0.35)
    fig.savefig(path, dpi=JOURNAL_DPI, bbox_inches="tight", facecolor="white")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _robust_vmax(matrix: np.ndarray, *, default: float) -> float:
    positive = matrix[np.isfinite(matrix) & (matrix > 0)]
    if positive.size == 0:
        return default
    return max(default, float(np.percentile(positive, 99)))


def _mean_normalize_coverage(coverage: np.ndarray) -> np.ndarray:
    if coverage.size == 0:
        return coverage.astype(float)
    mean_value = float(np.mean(coverage))
    if mean_value <= 0:
        return np.zeros_like(coverage, dtype=float)
    return coverage / mean_value


def _max_normalize_coverage(coverage: np.ndarray) -> np.ndarray:
    if coverage.size == 0:
        return coverage.astype(float)
    max_value = float(np.max(coverage))
    if max_value <= 0:
        return np.zeros_like(coverage, dtype=float)
    return coverage / max_value


def _add_interval_to_bins(
    bins: np.ndarray,
    interval: tuple[int, int],
    transcript_length: int,
) -> None:
    if transcript_length <= 0:
        return
    start, end = interval
    if end <= start:
        return
    bin_num = len(bins)
    first_bin = max(0, min(bin_num - 1, int(start * bin_num / transcript_length)))
    last_bin = max(0, min(bin_num - 1, int((end - 1) * bin_num / transcript_length)))
    for bin_index in range(first_bin, last_bin + 1):
        bin_start = bin_index * transcript_length / bin_num
        bin_end = (bin_index + 1) * transcript_length / bin_num
        overlap = max(0.0, min(end, bin_end) - max(start, bin_start))
        bin_width = bin_end - bin_start
        if bin_width > 0:
            bins[bin_index] = max(bins[bin_index], overlap / bin_width)


def _bin_centers(bin_num: int) -> np.ndarray:
    if bin_num <= 0:
        return np.array([], dtype=float)
    return (np.arange(bin_num, dtype=float) + 0.5) * 100 / bin_num


def _bin_axis_ticks(bin_num: int) -> list[int]:
    if bin_num <= 1:
        return [0]
    return [0, bin_num // 4, bin_num // 2, (3 * bin_num) // 4, bin_num - 1]


def _bin_axis_tick_labels(bin_num: int) -> list[str]:
    if bin_num <= 1:
        return ["0"]
    return ["0", "25", "50", "75", "100"]


def _infer_bin_num(per_transcript_coverage: dict[str, np.ndarray]) -> int:
    for values in per_transcript_coverage.values():
        return int(values.size)
    return 1


def _heatmap_yticks(row_count: int) -> list[int]:
    if row_count <= 30:
        return list(range(row_count))
    return sorted(set(np.linspace(0, row_count - 1, num=12, dtype=int).tolist()))
