from __future__ import annotations

from isocomp.annotation import parse_bed12_line
from isocomp.projection import project_blocks_to_transcript


def test_positive_projection_full_length() -> None:
    transcript = parse_bed12_line("chr1\t100\t400\tT1\t0\t+\t100\t400\t0\t2\t100,100\t0,200")
    projection = project_blocks_to_transcript([(100, 200), (300, 400)], transcript)

    assert projection.intervals == [(0, 200)]
    assert projection.coverage_fraction == 1.0
    assert projection.dist_to_5p == 0
    assert projection.dist_to_3p == 0


def test_positive_projection_truncated_ends() -> None:
    transcript = parse_bed12_line("chr1\t100\t400\tT1\t0\t+\t100\t400\t0\t2\t100,100\t0,200")
    projection = project_blocks_to_transcript([(150, 200), (300, 350)], transcript)

    assert projection.intervals == [(50, 150)]
    assert projection.coverage_fraction == 0.5
    assert projection.dist_to_5p == 50
    assert projection.dist_to_3p == 50


def test_negative_projection_is_oriented_5p_to_3p() -> None:
    transcript = parse_bed12_line("chr2\t100\t400\tTneg\t0\t-\t100\t400\t0\t2\t100,100\t0,200")
    projection = project_blocks_to_transcript([(100, 200), (300, 400)], transcript)

    assert projection.intervals == [(0, 200)]
    assert projection.coverage_fraction == 1.0
    assert projection.dist_to_5p == 0
    assert projection.dist_to_3p == 0

