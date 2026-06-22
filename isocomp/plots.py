"""Matplotlib plot generation."""

from __future__ import annotations

import os
import math
import random
import tempfile
from dataclasses import dataclass, field
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
from .metrics import OnlineNumericSummary

JOURNAL_DPI = 300
DISTANCE_PLOT_PERCENTILE = 95.0
DISTANCE_PLOT_SOFT_MAX_BP = 2_000.0
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


@dataclass
class BoundedPlotValues:
    max_values: int = 100_000
    values: list[float] = field(default_factory=list)
    summary: OnlineNumericSummary = field(default_factory=OnlineNumericSummary)
    seen_count: int = 0
    rng: random.Random = field(default_factory=lambda: random.Random(0), repr=False)

    def add(self, value: float | int) -> None:
        numeric = float(value)
        self.summary.add(numeric)
        self.seen_count += 1
        if len(self.values) < self.max_values:
            self.values.append(numeric)
            return

        replacement_index = self.rng.randrange(self.seen_count)
        if replacement_index < self.max_values:
            self.values[replacement_index] = numeric


@dataclass
class PlotData:
    coverage_fractions: BoundedPlotValues = field(default_factory=BoundedPlotValues)
    dist_to_5p: BoundedPlotValues = field(default_factory=BoundedPlotValues)
    dist_to_3p: BoundedPlotValues = field(default_factory=BoundedPlotValues)
    unique_read_count: int = 0
    unique_terminal_read_count: int = 0
    unique_5p_complete_count: int = 0
    unique_3p_complete_count: int = 0
    unique_full_length_like_count: int = 0
    read_body_rows: list[tuple[str, str, float, float, np.ndarray]] = field(default_factory=list)
    read_body_row_seen_count: int = 0
    read_body_row_rng: random.Random = field(default_factory=lambda: random.Random(1), repr=False)


def write_plots(
    plots_dir: Path,
    read_metrics: list[ReadMetrics],
    body_coverage: np.ndarray,
    per_transcript_coverage: dict[str, np.ndarray] | None = None,
    assignments: list[AssignmentResult] | None = None,
    plot_data: PlotData | None = None,
    tss_tol: int = 100,
    tes_tol: int = 100,
) -> None:
    plots_dir.mkdir(parents=True, exist_ok=True)
    if plot_data is None:
        plot_data = build_plot_data(
            read_metrics,
            assignments or [],
            bin_num=len(body_coverage),
        )
    with plt.rc_context(JOURNAL_RC_PARAMS):
        _plot_body_coverage(plots_dir / "transcript_body_coverage.png", body_coverage)
        _plot_transcript_body_heatmap(
            plots_dir / "transcript_body_heatmap.png",
            per_transcript_coverage or {},
            bin_num=len(body_coverage),
        )
        _plot_transcript_body_heatmap(
            plots_dir / "transcript_body_heatmap_full.png",
            per_transcript_coverage or {},
            max_rows=None,
            bin_num=len(body_coverage),
        )
        _plot_read_body_heatmap_from_rows(
            plots_dir / "read_body_heatmap.png",
            plot_data.read_body_rows,
            bin_num=len(body_coverage),
        )
        _plot_hist(
            plots_dir / "read_coverage_fraction.png",
            plot_data.coverage_fractions,
            "Transcript coverage fraction",
            "Read count",
            bins=30,
            range_=(0, 1),
        )
        _plot_distance_hist(
            plots_dir / "dist_to_5p.png",
            plot_data.dist_to_5p,
            "Distance to 5' end (bp)",
            tolerance=tss_tol,
        )
        _plot_distance_hist(
            plots_dir / "dist_to_3p.png",
            plot_data.dist_to_3p,
            "Distance to 3' end (bp)",
            tolerance=tes_tol,
        )
        _plot_full_length_fraction_from_counts(
            plots_dir / "full_length_fraction.png",
            plot_data,
        )


def build_plot_data(
    read_metrics: list[ReadMetrics],
    assignments: list[AssignmentResult],
    *,
    bin_num: int,
    max_read_heatmap_rows: int = 500,
) -> PlotData:
    plot_data = PlotData()
    assignment_by_read = {assignment.read_id: assignment for assignment in assignments}
    for metric in read_metrics:
        update_plot_data(
            plot_data,
            metric,
            assignment_by_read.get(metric.read_id),
            bin_num=bin_num,
            max_read_heatmap_rows=max_read_heatmap_rows,
        )
    return plot_data


def update_plot_data(
    plot_data: PlotData,
    metric: ReadMetrics,
    assignment: AssignmentResult | None,
    *,
    bin_num: int,
    max_read_heatmap_rows: int = 500,
) -> None:
    if metric.assignment_status != "unique":
        return
    plot_data.unique_read_count += 1
    if metric.coverage_fraction is not None:
        plot_data.coverage_fractions.add(metric.coverage_fraction)
    if metric.dist_to_5p is not None:
        plot_data.dist_to_5p.add(metric.dist_to_5p)
    if metric.dist_to_3p is not None:
        plot_data.dist_to_3p.add(metric.dist_to_3p)
    if metric.is_5p_complete is not None and metric.is_3p_complete is not None:
        plot_data.unique_terminal_read_count += 1
    if metric.is_5p_complete is True:
        plot_data.unique_5p_complete_count += 1
    if metric.is_3p_complete is True:
        plot_data.unique_3p_complete_count += 1
    if metric.is_full_length_like:
        plot_data.unique_full_length_like_count += 1

    row = _assignment_to_read_body_row(assignment, bin_num=bin_num)
    if row is not None:
        _add_read_body_row(
            plot_data,
            row,
            max_read_heatmap_rows=max_read_heatmap_rows,
        )


def _plot_body_coverage(path: Path, coverage: np.ndarray) -> None:
    normalized = _max_normalize_coverage(coverage)
    fig, ax = plt.subplots(figsize=(4.2, 3.0))
    x = _bin_centers(len(normalized))
    if x.size:
        ax.plot(
            x,
            normalized,
            color=COLORS["blue"],
            linewidth=2.0,
            solid_capstyle="round",
            solid_joinstyle="round",
        )
    ax.set_xlabel("Transcript position (5' to 3', %)")
    ax.set_ylabel("Relative coverage (max = 1)")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 1.05)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xticklabels(["5' end", "25", "50", "75", "3' end"])
    ax.set_yticks([0, 0.25, 0.50, 0.75, 1.00])
    _style_axis(ax, grid_axis="y")
    _save_figure(fig, path)


def _plot_transcript_body_heatmap(
    path: Path,
    per_transcript_coverage: dict[str, np.ndarray],
    *,
    max_rows: int | None = 200,
    bin_num: int | None = None,
) -> None:
    covered = _sorted_transcript_body_heatmap_rows(
        per_transcript_coverage,
        max_rows=max_rows,
    )
    bin_num = bin_num or _infer_bin_num(per_transcript_coverage)

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
        ax.set_yticks([])
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
    ax.set_ylabel("Transcripts sorted by coverage")
    ax.set_xticks(_bin_axis_ticks(bin_num))
    ax.set_xticklabels(_terminal_axis_tick_labels(bin_num))
    _style_axis(ax)
    _save_figure(fig, path)


def _sorted_transcript_body_heatmap_rows(
    per_transcript_coverage: dict[str, np.ndarray],
    *,
    max_rows: int | None,
) -> list[tuple[str, np.ndarray]]:
    covered = [
        (transcript_id, values)
        for transcript_id, values in per_transcript_coverage.items()
        if values.size and float(np.sum(values)) > 0
    ]
    covered.sort(key=lambda item: (-float(np.sum(item[1])), item[0]))
    if max_rows is None:
        return covered
    return covered[:max_rows]


def _assignment_to_read_body_row(
    assignment: AssignmentResult | None,
    *,
    bin_num: int,
) -> tuple[str, str, float, float, np.ndarray] | None:
    if assignment is None:
        return None
    if assignment.status != "unique" or assignment.transcript is None or assignment.projection is None:
        return None
    projection = assignment.projection
    if projection.read_start_tx is None or projection.read_end_tx is None:
        return None
    transcript_length = assignment.transcript.transcript_length
    if transcript_length <= 0:
        return None
    row = np.zeros(bin_num, dtype=float)
    for interval in projection.intervals:
        _add_interval_to_bins(row, interval, transcript_length)
    start_fraction = projection.read_start_tx / transcript_length
    dist_3p_fraction = (transcript_length - projection.read_end_tx) / transcript_length
    return (
        assignment.read_id,
        assignment.transcript.transcript_id,
        start_fraction,
        dist_3p_fraction,
        row,
    )


def _plot_read_body_heatmap_from_rows(
    path: Path,
    rows: list[tuple[str, str, float, float, np.ndarray]],
    *,
    bin_num: int,
    max_rows: int = 500,
) -> None:
    rows.sort(key=lambda item: (item[2], item[3], item[1], item[0]))
    rows = _evenly_downsample_rows(rows, max_rows)

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
    ax.set_ylabel("Unique reads sorted by 5'/3' distance")
    ax.set_xticks(_bin_axis_ticks(bin_num))
    ax.set_xticklabels(_terminal_axis_tick_labels(bin_num))
    _style_axis(ax)
    _save_figure(fig, path)


def _plot_hist(
    path: Path,
    values: BoundedPlotValues,
    xlabel: str,
    ylabel: str,
    *,
    bins: int,
    range_: tuple[float, float] | None = None,
    clip_to_range: bool = False,
    reference_value: float | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    sampled_values = np.asarray(values.values, dtype=float)
    has_upper_overflow = bool(
        clip_to_range
        and range_ is not None
        and sampled_values.size
        and np.any(sampled_values > range_[1])
    )
    if clip_to_range and range_ is not None:
        sampled_values = np.clip(sampled_values, range_[0], range_[1])
    if sampled_values.size:
        weights = None
        if values.seen_count > len(sampled_values):
            scale = values.seen_count / len(sampled_values)
            weights = np.full(len(sampled_values), scale, dtype=float)
        ax.hist(
            sampled_values,
            bins=bins,
            range=range_,
            weights=weights,
            color=COLORS["blue"],
            edgecolor="white",
            linewidth=0.5,
        )
        if reference_value is not None:
            ax.axvline(
                reference_value,
                color=COLORS["green"],
                linestyle=(0, (2, 2)),
                linewidth=1.0,
                label=f"Tolerance ({reference_value:g} bp)",
            )
        _add_median_marker(
            ax,
            values.summary.median(),
            display_upper=range_[1] if clip_to_range and range_ is not None else None,
        )
    else:
        ax.hist([], bins=bins, range=range_, color=COLORS["blue"], edgecolor="white")
        _empty_panel(ax, "No uniquely assigned reads")
    if range_ is not None:
        ax.set_xlim(*range_)
    else:
        ax.set_xlim(left=0)
    if has_upper_overflow and range_ is not None:
        xlabel = f"{xlabel}; >{range_[1]:g} pooled"
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    _style_axis(ax, grid_axis="y")
    _save_figure(fig, path)


def _plot_distance_hist(
    path: Path,
    values: BoundedPlotValues,
    xlabel: str,
    *,
    tolerance: int,
) -> None:
    display_upper = _distance_display_upper(values.values, tolerance=tolerance)
    _plot_hist(
        path,
        values,
        xlabel,
        "Read count",
        bins=40,
        range_=(0.0, display_upper),
        clip_to_range=True,
        reference_value=float(tolerance),
    )


def _distance_display_upper(values: list[float], *, tolerance: int) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite) & (finite >= 0)]
    percentile = (
        float(np.percentile(finite, DISTANCE_PLOT_PERCENTILE))
        if finite.size
        else 0.0
    )
    target = max(
        1.0,
        2.0 * tolerance,
        min(percentile, DISTANCE_PLOT_SOFT_MAX_BP),
    )
    return _nice_axis_upper(target)


def _nice_axis_upper(value: float) -> float:
    magnitude = 10 ** math.floor(math.log10(value))
    normalized = value / magnitude
    for step in (1.0, 2.0, 5.0, 10.0):
        if normalized <= step:
            return step * magnitude
    return 10.0 * magnitude


def _plot_full_length_fraction(path: Path, read_metrics: list[ReadMetrics]) -> None:
    unique = [
        item
        for item in read_metrics
        if item.assignment_status == "unique"
        and item.is_5p_complete is not None
        and item.is_3p_complete is not None
    ]
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
        _empty_panel(ax, "No terminal-evaluable unique reads", y=0.58)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of evaluable unique reads")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    _style_axis(ax, grid_axis="y")
    _save_figure(fig, path)


def _plot_full_length_fraction_from_counts(path: Path, plot_data: PlotData) -> None:
    unique_count = plot_data.unique_terminal_read_count
    fractions = [
        safe_fraction(plot_data.unique_5p_complete_count, unique_count),
        safe_fraction(plot_data.unique_3p_complete_count, unique_count),
        safe_fraction(plot_data.unique_full_length_like_count, unique_count),
    ]
    labels = ["5' complete", "3' complete", "Full-length\nlike"]
    colors = [COLORS["blue"], COLORS["green"], COLORS["vermillion"]]
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    x = np.arange(len(fractions))
    ax.bar(x, fractions, color=colors, edgecolor="white", linewidth=0.6, width=0.68)
    for index, value in enumerate(fractions):
        ax.text(index, min(value + 0.035, 1.02), f"{value:.2f}", ha="center", va="bottom", fontsize=7)
    if unique_count == 0:
        _empty_panel(ax, "No terminal-evaluable unique reads", y=0.58)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of evaluable unique reads")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    _style_axis(ax, grid_axis="y")
    _save_figure(fig, path)


def safe_fraction(numerator: int, denominator: int) -> float:
    if denominator == 0:
        return 0.0
    return numerator / denominator


def _add_read_body_row(
    plot_data: PlotData,
    row: tuple[str, str, float, float, np.ndarray],
    *,
    max_read_heatmap_rows: int,
) -> None:
    plot_data.read_body_row_seen_count += 1
    if len(plot_data.read_body_rows) < max_read_heatmap_rows:
        plot_data.read_body_rows.append(row)
        return

    replacement_index = plot_data.read_body_row_rng.randrange(
        plot_data.read_body_row_seen_count
    )
    if replacement_index < max_read_heatmap_rows:
        plot_data.read_body_rows[replacement_index] = row


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


def _add_median_marker(
    ax: plt.Axes,
    median_value: float,
    *,
    display_upper: float | None = None,
) -> None:
    if math.isnan(median_value):
        return
    marker_value = median_value
    label = "Median"
    if display_upper is not None and median_value > display_upper:
        marker_value = display_upper
        label = f"Median > {display_upper:g} bp"
    ax.axvline(
        marker_value,
        color=COLORS["vermillion"],
        linestyle=(0, (3, 2)),
        linewidth=1.0,
        label=label,
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


def _terminal_axis_tick_labels(bin_num: int) -> list[str]:
    if bin_num <= 1:
        return ["5' end"]
    return ["5' end", "25", "50", "75", "3' end"]


def _evenly_downsample_rows(
    rows: list[tuple[str, str, float, float, np.ndarray]],
    max_rows: int,
) -> list[tuple[str, str, float, float, np.ndarray]]:
    if max_rows <= 0:
        return []
    if len(rows) <= max_rows:
        return rows
    indexes = np.linspace(0, len(rows) - 1, num=max_rows, dtype=int)
    return [rows[int(index)] for index in indexes]


def _infer_bin_num(per_transcript_coverage: dict[str, np.ndarray]) -> int:
    for values in per_transcript_coverage.values():
        return int(values.size)
    return 1


def _heatmap_yticks(row_count: int) -> list[int]:
    if row_count <= 30:
        return list(range(row_count))
    return sorted(set(np.linspace(0, row_count - 1, num=12, dtype=int).tolist()))
