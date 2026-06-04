"""Small interval and numeric helpers."""

from __future__ import annotations

from collections.abc import Iterable

from .models import Interval


def safe_divide(numerator: float, denominator: float) -> float:
    if denominator == 0:
        return 0.0
    return numerator / denominator


def overlap_interval(a: Interval, b: Interval) -> Interval | None:
    start = max(a[0], b[0])
    end = min(a[1], b[1])
    if end <= start:
        return None
    return start, end


def overlap_length(a: Interval, b: Interval) -> int:
    overlap = overlap_interval(a, b)
    if overlap is None:
        return 0
    return overlap[1] - overlap[0]


def total_overlap_length(query_intervals: Iterable[Interval], target_intervals: Iterable[Interval]) -> int:
    query = merge_intervals(query_intervals)
    target = merge_intervals(target_intervals)
    i = 0
    j = 0
    total = 0
    while i < len(query) and j < len(target):
        q = query[i]
        t = target[j]
        total += overlap_length(q, t)
        if q[1] <= t[1]:
            i += 1
        else:
            j += 1
    return total


def merge_intervals(intervals: Iterable[Interval]) -> list[Interval]:
    sorted_intervals = sorted((start, end) for start, end in intervals if end > start)
    if not sorted_intervals:
        return []

    merged = [sorted_intervals[0]]
    for start, end in sorted_intervals[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return merged


def harmonic_mean(a: float, b: float) -> float:
    if a <= 0 or b <= 0:
        return 0.0
    return 2 * a * b / (a + b)


def bool_to_int(value: bool | None) -> int | str:
    if value is None:
        return ""
    return 1 if value else 0

