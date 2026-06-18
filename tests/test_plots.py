from __future__ import annotations

import numpy as np

from isocomp.plots import BoundedPlotValues, PlotData, _add_read_body_row


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
