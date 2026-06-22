from __future__ import annotations

import numpy as np

from isocomp.plots import (
    BoundedPlotValues,
    PlotData,
    _add_read_body_row,
    _distance_display_upper,
    _plot_transcript_body_heatmap,
    _sorted_transcript_body_heatmap_rows,
    write_plots,
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


def test_write_plots_creates_full_transcript_body_heatmap(tmp_path) -> None:
    write_plots(
        tmp_path,
        read_metrics=[],
        body_coverage=np.ones(2),
        per_transcript_coverage={"tx": np.ones(2)},
    )

    assert (tmp_path / "transcript_body_heatmap_full.png").exists()
    assert (tmp_path / "transcript_body_heatmap_full.pdf").exists()


def test_transcript_body_heatmap_rows_are_sorted_by_coverage() -> None:
    rows = _sorted_transcript_body_heatmap_rows(
        {
            "tx_low": np.array([0.5, 0.5]),
            "tx_high": np.array([2.0, 3.0]),
            "tx_mid_b": np.array([1.0, 1.0]),
            "tx_mid_a": np.array([1.0, 1.0]),
            "tx_empty": np.array([0.0, 0.0]),
        },
        max_rows=3,
    )

    assert [transcript_id for transcript_id, _ in rows] == [
        "tx_high",
        "tx_mid_a",
        "tx_mid_b",
    ]


def test_transcript_body_heatmap_rows_can_include_all_transcripts() -> None:
    rows = _sorted_transcript_body_heatmap_rows(
        {
            "tx_low": np.array([0.5, 0.5]),
            "tx_high": np.array([2.0, 3.0]),
            "tx_mid": np.array([1.0, 1.0]),
            "tx_empty": np.array([0.0, 0.0]),
        },
        max_rows=None,
    )

    assert [transcript_id for transcript_id, _ in rows] == [
        "tx_high",
        "tx_mid",
        "tx_low",
    ]
