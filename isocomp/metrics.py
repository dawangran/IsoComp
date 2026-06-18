"""Completeness metrics and aggregate summaries."""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import asdict, dataclass, field

import numpy as np

from .models import AssignmentResult, ReadAlignment, ReadMetrics, RunSummary, Transcript
from .utils import bool_to_int, safe_divide

READ_ASSIGNMENT_COLUMNS = [
    "read_id",
    "chrom",
    "gene_id",
    "transcript_id",
    "assignment_status",
    "assignment_score",
    "second_best_transcript",
    "second_best_score",
    "transcript_length",
    "read_aligned_length",
    "read_start_tx",
    "read_end_tx",
    "start_pct",
    "end_pct",
    "coverage_fraction",
    "dist_to_5p",
    "dist_to_3p",
    "is_5p_complete",
    "is_3p_complete",
    "is_full_length_like",
    "junction_match_count",
    "junction_precision",
    "junction_recall",
    "exon_overlap_score",
    "mapq",
    "softclip_5p",
    "softclip_3p",
]


@dataclass
class OnlineNumericSummary:
    count: int = 0
    total: float = 0.0
    max_exact_values: int = 100_000
    exact_values: list[float] = field(default_factory=list)
    median_estimator: "P2MedianEstimator" = field(default_factory=lambda: P2MedianEstimator())

    def add(self, value: float | int) -> None:
        numeric = float(value)
        self.count += 1
        self.total += numeric
        if len(self.exact_values) < self.max_exact_values:
            self.exact_values.append(numeric)
        elif self.exact_values:
            self.exact_values.clear()
        self.median_estimator.add(numeric)

    def mean(self) -> float:
        return safe_divide(self.total, self.count) if self.count else math.nan

    def median(self) -> float:
        if self.exact_values:
            return _median(self.exact_values)
        return self.median_estimator.median()


@dataclass
class P2MedianEstimator:
    """Constant-memory P² estimator for the 0.5 quantile.

    The first five observations are exact. After that, the estimate is updated
    online with Jain and Chlamtac's marker algorithm.
    """

    initial_values: list[float] = field(default_factory=list)
    marker_positions: list[int] = field(default_factory=list)
    desired_positions: list[float] = field(default_factory=list)
    marker_heights: list[float] = field(default_factory=list)

    def add(self, value: float) -> None:
        if len(self.initial_values) < 5:
            self.initial_values.append(value)
            if len(self.initial_values) == 5:
                self.initial_values.sort()
                self.marker_heights = list(self.initial_values)
                self.marker_positions = [1, 2, 3, 4, 5]
                self.desired_positions = [1.0, 2.0, 3.0, 4.0, 5.0]
            return

        heights = self.marker_heights
        if value < heights[0]:
            heights[0] = value
            bucket = 0
        elif value >= heights[4]:
            heights[4] = value
            bucket = 3
        else:
            bucket = 0
            for index in range(4):
                if heights[index] <= value < heights[index + 1]:
                    bucket = index
                    break

        for index in range(bucket + 1, 5):
            self.marker_positions[index] += 1
        increments = [0.0, 0.25, 0.5, 0.75, 1.0]
        for index, increment in enumerate(increments):
            self.desired_positions[index] += increment

        for index in range(1, 4):
            delta = self.desired_positions[index] - self.marker_positions[index]
            if (
                delta >= 1
                and self.marker_positions[index + 1] - self.marker_positions[index] > 1
            ):
                self._adjust_marker(index, 1)
            elif (
                delta <= -1
                and self.marker_positions[index - 1] - self.marker_positions[index] < -1
            ):
                self._adjust_marker(index, -1)

    def median(self) -> float:
        if not self.initial_values:
            return math.nan
        if len(self.initial_values) < 5:
            values = sorted(self.initial_values)
            middle = len(values) // 2
            if len(values) % 2:
                return values[middle]
            return (values[middle - 1] + values[middle]) / 2
        return self.marker_heights[2]

    def _adjust_marker(self, index: int, direction: int) -> None:
        positions = self.marker_positions
        heights = self.marker_heights
        lower = index - 1
        upper = index + 1
        denominator = positions[upper] - positions[lower]
        if denominator == 0:
            return

        parabolic = heights[index] + direction / denominator * (
            (positions[index] - positions[lower] + direction)
            * (heights[upper] - heights[index])
            / (positions[upper] - positions[index])
            + (positions[upper] - positions[index] - direction)
            * (heights[index] - heights[lower])
            / (positions[index] - positions[lower])
        )
        if heights[lower] < parabolic < heights[upper]:
            heights[index] = parabolic
        else:
            neighbor = index + direction
            heights[index] = heights[index] + direction * (
                heights[neighbor] - heights[index]
            ) / (positions[neighbor] - positions[index])
        positions[index] += direction


@dataclass
class SampleMetricAccumulator:
    aligned_lengths: OnlineNumericSummary = field(default_factory=OnlineNumericSummary)
    unique_coverage_values: OnlineNumericSummary = field(default_factory=OnlineNumericSummary)
    unique_dist_5p: OnlineNumericSummary = field(default_factory=OnlineNumericSummary)
    unique_dist_3p: OnlineNumericSummary = field(default_factory=OnlineNumericSummary)
    unique_5p_complete_count: int = 0
    unique_3p_complete_count: int = 0
    unique_full_length_like_count: int = 0
    unique_read_count: int = 0
    all_assigned_coverage_values: OnlineNumericSummary = field(default_factory=OnlineNumericSummary)
    all_assigned_dist_5p: OnlineNumericSummary = field(default_factory=OnlineNumericSummary)
    all_assigned_dist_3p: OnlineNumericSummary = field(default_factory=OnlineNumericSummary)
    all_assigned_5p_complete_count: int = 0
    all_assigned_3p_complete_count: int = 0
    all_assigned_full_length_like_count: int = 0
    all_assigned_read_count: int = 0


@dataclass
class TranscriptMetricAccumulator:
    assigned_read_count: int = 0
    unique_read_count: int = 0
    coverage_values: OnlineNumericSummary = field(default_factory=OnlineNumericSummary)
    dist_5p: OnlineNumericSummary = field(default_factory=OnlineNumericSummary)
    dist_3p: OnlineNumericSummary = field(default_factory=OnlineNumericSummary)
    full_length_like_count: int = 0


def build_read_metrics(
    read: ReadAlignment,
    assignment: AssignmentResult,
    *,
    tss_tol: int,
    tes_tol: int,
) -> ReadMetrics:
    transcript = assignment.transcript
    projection = assignment.projection
    transcript_length = transcript.transcript_length if transcript else None
    read_start_tx = projection.read_start_tx if projection else None
    read_end_tx = projection.read_end_tx if projection else None
    coverage_fraction = projection.coverage_fraction if projection else None
    dist_to_5p = projection.dist_to_5p if projection else None
    dist_to_3p = projection.dist_to_3p if projection else None
    start_pct = safe_divide(read_start_tx, transcript_length) if read_start_tx is not None and transcript_length else None
    end_pct = safe_divide(read_end_tx, transcript_length) if read_end_tx is not None and transcript_length else None
    is_5p_complete = dist_to_5p <= tss_tol if dist_to_5p is not None else None
    is_3p_complete = dist_to_3p <= tes_tol if dist_to_3p is not None else None
    softclip_5p, softclip_3p = _terminal_softclips(read, transcript)
    is_full_length_like = bool(
        assignment.status == "unique"
        and coverage_fraction is not None
        and coverage_fraction >= 0.8
        and is_5p_complete
        and is_3p_complete
    )

    return ReadMetrics(
        read_id=read.read_id,
        chrom=read.chrom,
        gene_id=transcript.gene_id if transcript else None,
        transcript_id=transcript.transcript_id if transcript else None,
        assignment_status=assignment.status,
        assignment_score=assignment.score,
        second_best_transcript=assignment.second_best_transcript,
        second_best_score=assignment.second_best_score,
        transcript_length=transcript_length,
        read_aligned_length=read.aligned_length,
        read_start_tx=read_start_tx,
        read_end_tx=read_end_tx,
        start_pct=start_pct,
        end_pct=end_pct,
        coverage_fraction=coverage_fraction,
        dist_to_5p=dist_to_5p,
        dist_to_3p=dist_to_3p,
        is_5p_complete=is_5p_complete,
        is_3p_complete=is_3p_complete,
        is_full_length_like=is_full_length_like,
        junction_match_count=assignment.junction_match_count,
        junction_precision=assignment.junction_precision,
        junction_recall=assignment.junction_recall,
        exon_overlap_score=assignment.exon_overlap_score,
        mapq=read.mapq,
        softclip_5p=softclip_5p,
        softclip_3p=softclip_3p,
    )


def read_metrics_to_row(metrics: ReadMetrics) -> dict[str, object]:
    row = asdict(metrics)
    row["is_5p_complete"] = bool_to_int(metrics.is_5p_complete)
    row["is_3p_complete"] = bool_to_int(metrics.is_3p_complete)
    row["is_full_length_like"] = bool_to_int(metrics.is_full_length_like)
    return row


def summarize_sample(read_metrics: list[ReadMetrics], summary: RunSummary, sample_name: str) -> dict[str, object]:
    unique_reads = [item for item in read_metrics if item.assignment_status == "unique"]
    assigned_reads = [
        item
        for item in read_metrics
        if item.assignment_status in {"unique", "ambiguous"}
    ]
    coverage_values = [
        item.coverage_fraction
        for item in unique_reads
        if item.coverage_fraction is not None
    ]
    all_assigned_coverage_values = [
        item.coverage_fraction
        for item in assigned_reads
        if item.coverage_fraction is not None
    ]
    dist_5p = [item.dist_to_5p for item in unique_reads if item.dist_to_5p is not None]
    dist_3p = [item.dist_to_3p for item in unique_reads if item.dist_to_3p is not None]
    all_assigned_dist_5p = [
        item.dist_to_5p for item in assigned_reads if item.dist_to_5p is not None
    ]
    all_assigned_dist_3p = [
        item.dist_to_3p for item in assigned_reads if item.dist_to_3p is not None
    ]
    aligned_lengths = [item.read_aligned_length for item in read_metrics]

    return {
        "sample": sample_name,
        "total_reads": summary.total_reads,
        "mapped_reads": summary.mapped_reads,
        "primary_reads": summary.primary_reads,
        "assigned_reads": summary.assigned_reads,
        "unique_assigned_reads": summary.unique_assigned_reads,
        "ambiguous_reads": summary.ambiguous_reads,
        "gene_only_reads": summary.gene_only_reads,
        "low_confidence_reads": summary.low_confidence_reads,
        "unassigned_reads": summary.unassigned_reads,
        "median_read_aligned_length": _median(aligned_lengths),
        "median_transcript_coverage_fraction": _median(coverage_values),
        "mean_transcript_coverage_fraction": _mean(coverage_values),
        "5p_complete_fraction": safe_divide(
            sum(1 for item in unique_reads if item.is_5p_complete),
            len(unique_reads),
        ),
        "3p_complete_fraction": safe_divide(
            sum(1 for item in unique_reads if item.is_3p_complete),
            len(unique_reads),
        ),
        "full_length_like_fraction": safe_divide(
            sum(1 for item in unique_reads if item.is_full_length_like),
            len(unique_reads),
        ),
        "median_dist_to_5p": _median(dist_5p),
        "median_dist_to_3p": _median(dist_3p),
        "all_assigned_median_transcript_coverage_fraction": _median(
            all_assigned_coverage_values
        ),
        "all_assigned_mean_transcript_coverage_fraction": _mean(
            all_assigned_coverage_values
        ),
        "all_assigned_5p_complete_fraction": safe_divide(
            sum(1 for item in assigned_reads if item.is_5p_complete is True),
            len(assigned_reads),
        ),
        "all_assigned_3p_complete_fraction": safe_divide(
            sum(1 for item in assigned_reads if item.is_3p_complete is True),
            len(assigned_reads),
        ),
        "all_assigned_full_length_like_fraction": safe_divide(
            sum(1 for item in assigned_reads if _is_terminal_full_length_like(item)),
            len(assigned_reads),
        ),
        "all_assigned_median_dist_to_5p": _median(all_assigned_dist_5p),
        "all_assigned_median_dist_to_3p": _median(all_assigned_dist_3p),
    }


def update_sample_accumulator(
    accumulator: SampleMetricAccumulator,
    metric: ReadMetrics,
) -> None:
    accumulator.aligned_lengths.add(metric.read_aligned_length)
    if metric.assignment_status == "unique":
        accumulator.unique_read_count += 1
        if metric.coverage_fraction is not None:
            accumulator.unique_coverage_values.add(metric.coverage_fraction)
        if metric.dist_to_5p is not None:
            accumulator.unique_dist_5p.add(metric.dist_to_5p)
        if metric.dist_to_3p is not None:
            accumulator.unique_dist_3p.add(metric.dist_to_3p)
        if metric.is_5p_complete is True:
            accumulator.unique_5p_complete_count += 1
        if metric.is_3p_complete is True:
            accumulator.unique_3p_complete_count += 1
        if metric.is_full_length_like:
            accumulator.unique_full_length_like_count += 1

    if metric.assignment_status in {"unique", "ambiguous"}:
        accumulator.all_assigned_read_count += 1
        if metric.coverage_fraction is not None:
            accumulator.all_assigned_coverage_values.add(metric.coverage_fraction)
        if metric.dist_to_5p is not None:
            accumulator.all_assigned_dist_5p.add(metric.dist_to_5p)
        if metric.dist_to_3p is not None:
            accumulator.all_assigned_dist_3p.add(metric.dist_to_3p)
        if metric.is_5p_complete is True:
            accumulator.all_assigned_5p_complete_count += 1
        if metric.is_3p_complete is True:
            accumulator.all_assigned_3p_complete_count += 1
        if _is_terminal_full_length_like(metric):
            accumulator.all_assigned_full_length_like_count += 1


def summarize_sample_accumulator(
    accumulator: SampleMetricAccumulator,
    summary: RunSummary,
    sample_name: str,
) -> dict[str, object]:
    return {
        "sample": sample_name,
        "total_reads": summary.total_reads,
        "mapped_reads": summary.mapped_reads,
        "primary_reads": summary.primary_reads,
        "assigned_reads": summary.assigned_reads,
        "unique_assigned_reads": summary.unique_assigned_reads,
        "ambiguous_reads": summary.ambiguous_reads,
        "gene_only_reads": summary.gene_only_reads,
        "low_confidence_reads": summary.low_confidence_reads,
        "unassigned_reads": summary.unassigned_reads,
        "median_read_aligned_length": accumulator.aligned_lengths.median(),
        "median_transcript_coverage_fraction": (
            accumulator.unique_coverage_values.median()
        ),
        "mean_transcript_coverage_fraction": (
            accumulator.unique_coverage_values.mean()
        ),
        "5p_complete_fraction": safe_divide(
            accumulator.unique_5p_complete_count,
            accumulator.unique_read_count,
        ),
        "3p_complete_fraction": safe_divide(
            accumulator.unique_3p_complete_count,
            accumulator.unique_read_count,
        ),
        "full_length_like_fraction": safe_divide(
            accumulator.unique_full_length_like_count,
            accumulator.unique_read_count,
        ),
        "median_dist_to_5p": accumulator.unique_dist_5p.median(),
        "median_dist_to_3p": accumulator.unique_dist_3p.median(),
        "all_assigned_median_transcript_coverage_fraction": (
            accumulator.all_assigned_coverage_values.median()
        ),
        "all_assigned_mean_transcript_coverage_fraction": (
            accumulator.all_assigned_coverage_values.mean()
        ),
        "all_assigned_5p_complete_fraction": safe_divide(
            accumulator.all_assigned_5p_complete_count,
            accumulator.all_assigned_read_count,
        ),
        "all_assigned_3p_complete_fraction": safe_divide(
            accumulator.all_assigned_3p_complete_count,
            accumulator.all_assigned_read_count,
        ),
        "all_assigned_full_length_like_fraction": safe_divide(
            accumulator.all_assigned_full_length_like_count,
            accumulator.all_assigned_read_count,
        ),
        "all_assigned_median_dist_to_5p": accumulator.all_assigned_dist_5p.median(),
        "all_assigned_median_dist_to_3p": accumulator.all_assigned_dist_3p.median(),
    }


def summarize_transcripts(
    transcripts: dict[str, Transcript],
    read_metrics: list[ReadMetrics],
    transcript_bin_coverage: dict[str, np.ndarray],
) -> list[dict[str, object]]:
    by_transcript: dict[str, list[ReadMetrics]] = defaultdict(list)
    assigned_counts: dict[str, int] = defaultdict(int)
    for metric in read_metrics:
        if metric.assignment_status in {"unique", "ambiguous"} and metric.transcript_id:
            assigned_counts[metric.transcript_id] += 1
        if metric.assignment_status == "unique" and metric.transcript_id:
            by_transcript[metric.transcript_id].append(metric)

    rows: list[dict[str, object]] = []
    for transcript_id in sorted(transcripts):
        transcript = transcripts[transcript_id]
        unique_items = by_transcript.get(transcript_id, [])
        coverage = [item.coverage_fraction for item in unique_items if item.coverage_fraction is not None]
        dist_5p = [item.dist_to_5p for item in unique_items if item.dist_to_5p is not None]
        dist_3p = [item.dist_to_3p for item in unique_items if item.dist_to_3p is not None]
        full_length_like_count = sum(1 for item in unique_items if item.is_full_length_like)
        bins = transcript_bin_coverage.get(transcript_id)
        bin_cv = _coefficient_of_variation(bins) if bins is not None else math.nan

        rows.append(
            {
                "transcript_id": transcript.transcript_id,
                "gene_id": transcript.gene_id or "",
                "transcript_length": transcript.transcript_length,
                "assigned_read_count": assigned_counts.get(transcript_id, 0),
                "unique_read_count": len(unique_items),
                "mean_coverage_fraction": _mean(coverage),
                "median_coverage_fraction": _median(coverage),
                "full_length_like_count": full_length_like_count,
                "full_length_like_fraction": safe_divide(full_length_like_count, len(unique_items)),
                "mean_dist_to_5p": _mean(dist_5p),
                "median_dist_to_5p": _median(dist_5p),
                "mean_dist_to_3p": _mean(dist_3p),
                "median_dist_to_3p": _median(dist_3p),
                "bin_coverage_cv": bin_cv,
                "body_uniformity_score": 1 / (1 + bin_cv) if not math.isnan(bin_cv) else math.nan,
            }
        )
    return rows


def init_transcript_accumulators(
    transcripts: dict[str, Transcript],
) -> dict[str, TranscriptMetricAccumulator]:
    return {
        transcript_id: TranscriptMetricAccumulator()
        for transcript_id in transcripts
    }


def update_transcript_accumulators(
    accumulators: dict[str, TranscriptMetricAccumulator],
    metric: ReadMetrics,
) -> None:
    if metric.transcript_id is None:
        return
    accumulator = accumulators.get(metric.transcript_id)
    if accumulator is None:
        return
    if metric.assignment_status in {"unique", "ambiguous"}:
        accumulator.assigned_read_count += 1
    if metric.assignment_status != "unique":
        return
    accumulator.unique_read_count += 1
    if metric.coverage_fraction is not None:
        accumulator.coverage_values.add(metric.coverage_fraction)
    if metric.dist_to_5p is not None:
        accumulator.dist_5p.add(metric.dist_to_5p)
    if metric.dist_to_3p is not None:
        accumulator.dist_3p.add(metric.dist_to_3p)
    if metric.is_full_length_like:
        accumulator.full_length_like_count += 1


def summarize_transcript_accumulators(
    transcripts: dict[str, Transcript],
    accumulators: dict[str, TranscriptMetricAccumulator],
    transcript_bin_coverage: dict[str, np.ndarray],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for transcript_id in sorted(transcripts):
        transcript = transcripts[transcript_id]
        accumulator = accumulators.get(transcript_id, TranscriptMetricAccumulator())
        bins = transcript_bin_coverage.get(transcript_id)
        bin_cv = _coefficient_of_variation(bins) if bins is not None else math.nan

        rows.append(
            {
                "transcript_id": transcript.transcript_id,
                "gene_id": transcript.gene_id or "",
                "transcript_length": transcript.transcript_length,
                "assigned_read_count": accumulator.assigned_read_count,
                "unique_read_count": accumulator.unique_read_count,
                "mean_coverage_fraction": accumulator.coverage_values.mean(),
                "median_coverage_fraction": accumulator.coverage_values.median(),
                "full_length_like_count": accumulator.full_length_like_count,
                "full_length_like_fraction": safe_divide(
                    accumulator.full_length_like_count,
                    accumulator.unique_read_count,
                ),
                "mean_dist_to_5p": accumulator.dist_5p.mean(),
                "median_dist_to_5p": accumulator.dist_5p.median(),
                "mean_dist_to_3p": accumulator.dist_3p.mean(),
                "median_dist_to_3p": accumulator.dist_3p.median(),
                "bin_coverage_cv": bin_cv,
                "body_uniformity_score": 1 / (1 + bin_cv) if not math.isnan(bin_cv) else math.nan,
            }
        )
    return rows


def compute_transcript_body_coverage(
    transcripts: dict[str, Transcript],
    assignments: list[AssignmentResult],
    *,
    bin_num: int,
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    aggregate = np.zeros(bin_num, dtype=float)
    per_transcript = {
        transcript_id: np.zeros(bin_num, dtype=float)
        for transcript_id in transcripts
    }

    for assignment in assignments:
        add_assignment_to_body_coverage(aggregate, per_transcript, assignment)

    return aggregate, per_transcript


def add_assignment_to_body_coverage(
    aggregate: np.ndarray,
    per_transcript: dict[str, np.ndarray],
    assignment: AssignmentResult,
) -> None:
    if assignment.status != "unique" or assignment.transcript is None or assignment.projection is None:
        return
    transcript = assignment.transcript
    tx_bins = per_transcript[transcript.transcript_id]
    for interval in assignment.projection.intervals:
        _add_interval_to_bins(tx_bins, interval, transcript.transcript_length)
        _add_interval_to_bins(aggregate, interval, transcript.transcript_length)


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
            bins[bin_index] += overlap / bin_width


def update_summary_for_metric(summary: RunSummary, metric: ReadMetrics) -> None:
    if metric.assignment_status in {"unique", "ambiguous"}:
        summary.assigned_reads += 1
    if metric.assignment_status == "unique":
        summary.unique_assigned_reads += 1
    elif metric.assignment_status == "ambiguous":
        summary.ambiguous_reads += 1
    elif metric.assignment_status == "low_confidence":
        summary.low_confidence_reads += 1
    elif metric.assignment_status == "unassigned":
        summary.unassigned_reads += 1


def _mean(values: list[float | int]) -> float:
    if not values:
        return math.nan
    return float(np.mean(values))


def _median(values: list[float | int]) -> float:
    if not values:
        return math.nan
    return float(np.median(values))


def _coefficient_of_variation(values: np.ndarray | None) -> float:
    if values is None or values.size == 0:
        return math.nan
    mean_value = float(np.mean(values))
    if mean_value == 0:
        return math.nan
    return float(np.std(values) / mean_value)


def _is_terminal_full_length_like(metric: ReadMetrics) -> bool:
    return bool(
        metric.coverage_fraction is not None
        and metric.coverage_fraction >= 0.8
        and metric.is_5p_complete is True
        and metric.is_3p_complete is True
    )


def _terminal_softclips(
    read: ReadAlignment,
    transcript: Transcript | None,
) -> tuple[int, int]:
    if transcript is None or transcript.strand == ".":
        return read.softclip_5p, read.softclip_3p
    if transcript.strand == "-":
        return read.softclip_right, read.softclip_left
    return read.softclip_left, read.softclip_right
