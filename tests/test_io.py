from __future__ import annotations

import json
import math

from isocomp.io import write_json


def test_write_json_converts_nan_to_null(tmp_path) -> None:
    path = tmp_path / "stats.json"

    write_json(path, {"value": math.nan, "nested": {"ok": 1.0}})

    assert json.loads(path.read_text(encoding="utf-8")) == {
        "value": None,
        "nested": {"ok": 1.0},
    }
