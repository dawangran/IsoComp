"""Matplotlib plot generation."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

_MPLCONFIGDIR = Path(tempfile.gettempdir()) / "isocomp-matplotlib"
_MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPLCONFIGDIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .models import AssignmentResult, ReadMetrics


def write_plots(
    plots_dir: Path,
    read_metrics: list[ReadMetrics],
    body_coverage: np.ndarray,
    per_transcript_coverage: dict[str, np.ndarray] | None = None,
    assignments: list[AssignmentResult] | None = None,
) -> None:
    plots_dir.mkdir(parents=True, exist_ok=True)
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
        [item.coverage_fraction for item in read_metrics if item.assignment_status == "unique" and item.coverage_fraction is not None],
        "Transcript coverage fraction",
        "Read count",
        bins=30,
        range_=(0, 1),
    )
    _plot_hist(
        plots_dir / "dist_to_5p.png",
        [item.dist_to_5p for item in read_metrics if item.assignment_status == "unique" and item.dist_to_5p is not None],
        "Distance to 5' end (bp)",
        "Read count",
        bins=30,
    )
    _plot_hist(
        plots_dir / "dist_to_3p.png",
        [item.dist_to_3p for item in read_metrics if item.assignment_status == "unique" and item.dist_to_3p is not None],
        "Distance to 3' end (bp)",
        "Read count",
        bins=30,
    )
    _plot_full_length_fraction(plots_dir / "full_length_fraction.png", read_metrics)


def _plot_body_coverage(path: Path, coverage: np.ndarray) -> None:
    normalized = _max_normalize_coverage(coverage)
    fig, ax = plt.subplots(figsize=(7, 4))
    x = _bin_centers(len(normalized))
    ax.plot(x, normalized, color="#006d77", linewidth=2)
    ax.set_xlabel("Transcript body percentile (5'->3')")
    ax.set_ylabel("Coverage")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


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

    fig, ax = plt.subplots(figsize=(8, max(3.5, min(12, 0.16 * len(covered) + 2.5))))
    if covered:
        matrix = np.vstack([_mean_normalize_coverage(values) for _, values in covered])
        image = ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap="viridis", vmin=0)
        colorbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
        colorbar.set_label("Mean-normalized coverage")
        ax.set_yticks(_heatmap_yticks(len(covered)))
        ax.set_yticklabels([covered[index][0] for index in _heatmap_yticks(len(covered))], fontsize=7)
    else:
        ax.imshow(np.zeros((1, 1)), aspect="auto", interpolation="nearest", cmap="viridis", vmin=0, vmax=1)
        ax.set_yticks([])
        ax.text(0.5, 0.5, "No uniquely assigned transcripts", ha="center", va="center", transform=ax.transAxes)

    ax.set_xlabel("Normalized transcript position (%)")
    ax.set_ylabel("Transcript")
    ax.set_xticks(_heatmap_xticks(per_transcript_coverage))
    ax.set_xticklabels(_heatmap_xtick_labels(per_transcript_coverage))
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


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

    fig, ax = plt.subplots(figsize=(8, max(3.5, min(14, 0.08 * len(rows) + 2.5))))
    if rows:
        matrix = np.vstack([item[4] for item in rows])
        image = ax.imshow(
            matrix,
            aspect="auto",
            interpolation="nearest",
            cmap="magma",
            vmin=0,
            vmax=1,
        )
        colorbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
        colorbar.set_label("Covered fraction of bin")
        tick_indexes = _heatmap_yticks(len(rows))
        ax.set_yticks(tick_indexes)
        ax.set_yticklabels([rows[index][0] for index in tick_indexes], fontsize=7)
    else:
        ax.imshow(np.zeros((1, max(1, bin_num))), aspect="auto", interpolation="nearest", cmap="magma", vmin=0, vmax=1)
        ax.set_yticks([])
        ax.text(0.5, 0.5, "No uniquely assigned reads", ha="center", va="center", transform=ax.transAxes)

    ax.set_xlabel("Assigned transcript percentile (5'->3')")
    ax.set_ylabel("Read")
    ax.set_xticks(_bin_axis_ticks(bin_num))
    ax.set_xticklabels(_bin_axis_tick_labels(bin_num))
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_hist(
    path: Path,
    values: list[float | int],
    xlabel: str,
    ylabel: str,
    *,
    bins: int,
    range_: tuple[float, float] | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(7, 4))
    if values:
        ax.hist(values, bins=bins, range=range_, color="#457b9d", edgecolor="white")
    else:
        ax.hist([], bins=bins, range=range_, color="#457b9d", edgecolor="white")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_full_length_fraction(path: Path, read_metrics: list[ReadMetrics]) -> None:
    unique = [item for item in read_metrics if item.assignment_status == "unique"]
    fraction = (
        sum(1 for item in unique if item.is_full_length_like) / len(unique)
        if unique
        else 0.0
    )
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(["sample"], [fraction], color="#2a9d8f")
    ax.set_ylim(0, 1)
    ax.set_ylabel("Full-length-like fraction")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


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


def _heatmap_xticks(per_transcript_coverage: dict[str, np.ndarray]) -> list[int]:
    bin_num = _infer_bin_num(per_transcript_coverage)
    return _bin_axis_ticks(bin_num)


def _heatmap_xtick_labels(per_transcript_coverage: dict[str, np.ndarray]) -> list[str]:
    bin_num = _infer_bin_num(per_transcript_coverage)
    return _bin_axis_tick_labels(bin_num)


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
