from __future__ import annotations

import pysam

from isocomp.bam_parser import cigar_to_blocks_and_junctions, parse_alignment_record, softclips_from_cigar


def test_cigar_block_and_junction_extraction() -> None:
    blocks, junctions = cigar_to_blocks_and_junctions(
        100,
        [(0, 100), (3, 400), (0, 150), (3, 300), (0, 200)],
    )

    assert blocks == [(100, 200), (600, 750), (1050, 1250)]
    assert junctions == [(200, 600), (750, 1050)]


def test_cigar_deletion_does_not_count_as_covered_block() -> None:
    blocks, junctions = cigar_to_blocks_and_junctions(
        100,
        [(0, 50), (2, 5), (0, 50)],
    )

    assert blocks == [(100, 150), (155, 205)]
    assert junctions == []


def test_softclip_extraction() -> None:
    assert softclips_from_cigar([(4, 5), (0, 90), (4, 10)]) == (5, 10)
    assert softclips_from_cigar(
        [(4, 5), (0, 90), (4, 10)],
        is_reverse=True,
    ) == (10, 5)


def test_parse_alignment_record(synthetic_bam) -> None:
    with pysam.AlignmentFile(synthetic_bam, "rb") as bam:
        read = next(bam.fetch(until_eof=True))
        parsed = parse_alignment_record(read)

    assert parsed.read_id == "full_pos"
    assert parsed.chrom == "chr1"
    assert parsed.blocks == [(100, 200), (300, 400)]
    assert parsed.junctions == [(200, 300)]
    assert parsed.aligned_length == 200
