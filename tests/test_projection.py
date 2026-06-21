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


def test_negative_projection_distinguishes_5p_and_3p_truncation() -> None:
    transcript = parse_bed12_line("chr2\t100\t400\tTneg\t0\t-\t100\t400\t0\t2\t100,100\t0,200")

    five_prime_projection = project_blocks_to_transcript([(350, 400)], transcript)
    three_prime_projection = project_blocks_to_transcript([(100, 150)], transcript)

    assert five_prime_projection.intervals == [(0, 50)]
    assert five_prime_projection.dist_to_5p == 0
    assert five_prime_projection.dist_to_3p == 150
    assert three_prime_projection.intervals == [(150, 200)]
    assert three_prime_projection.dist_to_5p == 150
    assert three_prime_projection.dist_to_3p == 0


def test_unknown_strand_projection_does_not_report_terminal_coordinates() -> None:
    transcript = parse_bed12_line(
        "chr3\t100\t400\tTunknown\t0\t.\t100\t400\t0\t2\t100,100\t0,200"
    )

    projection = project_blocks_to_transcript([(100, 150)], transcript)

    assert projection.coverage_fraction == 0.25
    assert projection.read_start_tx is None
    assert projection.read_end_tx is None
    assert projection.dist_to_5p is None
    assert projection.dist_to_3p is None
