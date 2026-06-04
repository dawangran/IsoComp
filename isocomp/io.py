"""Output writing helpers with atomic replacement."""

from __future__ import annotations

import json
import os
from collections.abc import Iterable
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any

import pandas as pd
import numpy as np

from .metrics import READ_ASSIGNMENT_COLUMNS
from .models import ReadMetrics


class OutputError(ValueError):
    """Raised when output paths are invalid or unsafe to overwrite."""


def output_paths(out_prefix: str | Path) -> dict[str, Path]:
    prefix = Path(out_prefix)
    return {
        "read_assignment": prefix.with_suffix(prefix.suffix + ".read_assignment.tsv")
        if prefix.suffix
        else Path(f"{prefix}.read_assignment.tsv"),
        "transcript_metrics": prefix.with_suffix(prefix.suffix + ".transcript_metrics.tsv")
        if prefix.suffix
        else Path(f"{prefix}.transcript_metrics.tsv"),
        "sample_summary": prefix.with_suffix(prefix.suffix + ".sample_summary.tsv")
        if prefix.suffix
        else Path(f"{prefix}.sample_summary.tsv"),
        "transcript_body_coverage": prefix.with_suffix(prefix.suffix + ".transcript_body_coverage.tsv")
        if prefix.suffix
        else Path(f"{prefix}.transcript_body_coverage.tsv"),
        "assignment_stats": prefix.with_suffix(prefix.suffix + ".assignment_stats.json")
        if prefix.suffix
        else Path(f"{prefix}.assignment_stats.json"),
        "plots_dir": Path(f"{prefix}.plots"),
    }


def ensure_outputs_available(paths: dict[str, Path], *, force: bool) -> None:
    existing = [path for key, path in paths.items() if key != "plots_dir" and path.exists()]
    plots_dir = paths["plots_dir"]
    if plots_dir.exists() and any(plots_dir.iterdir()):
        existing.append(plots_dir)
    if existing and not force:
        formatted = ", ".join(str(path) for path in existing)
        raise OutputError(
            "Output path(s) already exist. Use --force to overwrite intentionally: "
            f"{formatted}"
        )
    for key, path in paths.items():
        if key == "plots_dir":
            path.mkdir(parents=True, exist_ok=True)
        else:
            path.parent.mkdir(parents=True, exist_ok=True)


def write_read_assignment(path: Path, read_metrics: Iterable[ReadMetrics]) -> None:
    from .metrics import read_metrics_to_row

    rows = [read_metrics_to_row(item) for item in read_metrics]
    frame = pd.DataFrame(rows, columns=READ_ASSIGNMENT_COLUMNS)
    write_dataframe(path, frame)


def write_dataframe(path: Path, frame: pd.DataFrame) -> None:
    _atomic_write_text(path, lambda handle: frame.to_csv(handle, sep="\t", index=False))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    def writer(handle: Any) -> None:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")

    _atomic_write_text(path, writer)


def write_transcript_body_coverage(path: Path, coverage: Any) -> None:
    mean_normalized = _mean_normalize_coverage(coverage)
    max_normalized = _max_normalize_coverage(coverage)
    frame = pd.DataFrame(
        {
            "bin": list(range(1, len(coverage) + 1)),
            "coverage": coverage,
            "mean_normalized_coverage": mean_normalized,
            "max_normalized_coverage": max_normalized,
        }
    )
    write_dataframe(path, frame)


def _atomic_write_text(path: Path, writer: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("wt", encoding="utf-8", dir=path.parent, delete=False) as handle:
        tmp_name = handle.name
        writer(handle)
    os.replace(tmp_name, path)


def _mean_normalize_coverage(coverage: Any) -> Any:
    values = np.asarray(coverage, dtype=float)
    mean_value = float(np.mean(values)) if values.size else 0.0
    if mean_value <= 0:
        return np.zeros_like(values, dtype=float)
    return values / mean_value


def _max_normalize_coverage(coverage: Any) -> Any:
    values = np.asarray(coverage, dtype=float)
    max_value = float(np.max(values)) if values.size else 0.0
    if max_value <= 0:
        return np.zeros_like(values, dtype=float)
    return values / max_value
