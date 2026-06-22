from __future__ import annotations

import numpy as np

from isocomp.plots import (
    BoundedPlotValues,
    PlotData,
    _add_read_body_row,
    _distance_display_upper,
    _plot_transcript_body_heatmap,
)


def test_bounded_plot_values_keep_memory_bounded_and_summary_online() -> None:
    values = BoundedPlotValues(max_values=5)

    for value in range(20):
        values.add(value)

    assert values.seen_count == 20
    assert len(values.values) == 5
    assert values.summary.mean() == 9.5
    assert values.summary.median() == 9.5


def test_read_body_rows_keep_memory_bounded() -> None:
    plot_data = PlotData()

    for index in range(20):
        _add_read_body_row(
            plot_data,
            (f"read{index}", "tx", index / 20, 0.0, np.ones(4)),
            max_read_heatmap_rows=5,
        )

    assert plot_data.read_body_row_seen_count == 20
    assert len(plot_data.read_body_rows) == 5


def test_distance_display_upper_ignores_extreme_tail() -> None:
    values = list(range(100)) + [1_000_000]

    assert _distance_display_upper(values, tolerance=100) == 200


def test_distance_display_upper_keeps_tolerance_visible() -> None:
    assert _distance_display_upper([0, 10, 20], tolerance=1_500) >= 3_000


def test_empty_transcript_body_heatmap_uses_requested_bin_count(tmp_path) -> None:
    path = tmp_path / "empty_heatmap.png"

    _plot_transcript_body_heatmap(path, {}, bin_num=10)

    assert path.exists()
    assert path.with_suffix(".pdf").exists()
