from __future__ import annotations

import pytest

from isocomp.annotation import AnnotationError, parse_bed12_line, read_bed12


def test_bed12_positive_transcript_coordinates() -> None:
    transcript = parse_bed12_line("chr1\t100\t400\tT1\t0\t+\t100\t400\t0\t2\t100,100\t0,200")

    assert transcript.exons == [(100, 200), (300, 400)]
    assert transcript.transcript_length == 200
    assert transcript.junctions == [(200, 300)]
    assert [(exon.genomic_start, exon.genomic_end, exon.tx_start, exon.tx_end) for exon in transcript.tx_exons] == [
        (100, 200, 0, 100),
        (300, 400, 100, 200),
    ]


def test_bed12_negative_transcript_coordinates_reverse_exon_order() -> None:
    transcript = parse_bed12_line("chr2\t100\t400\tTneg\t0\t-\t100\t400\t0\t2\t100,100\t0,200")

    assert transcript.exons == [(100, 200), (300, 400)]
    assert transcript.junctions == [(200, 300)]
    assert [(exon.genomic_start, exon.genomic_end, exon.tx_start, exon.tx_end) for exon in transcript.tx_exons] == [
        (300, 400, 0, 100),
        (100, 200, 100, 200),
    ]


def test_read_bed12_rejects_malformed_line(tmp_path) -> None:
    path = tmp_path / "bad.bed12"
    path.write_text("chr1\t100\t200\n", encoding="utf-8")

    with pytest.raises(AnnotationError, match="expected at least 12"):
        read_bed12(path)


def test_bed12_rejects_overlapping_blocks() -> None:
    with pytest.raises(AnnotationError, match="Overlapping exons"):
        parse_bed12_line("chr1\t100\t250\tToverlap\t0\t+\t100\t250\t0\t2\t100,100\t0,50")
