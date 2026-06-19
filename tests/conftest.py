from __future__ import annotations

from pathlib import Path

import pysam
import pytest


@pytest.fixture
def bed12_path(tmp_path: Path) -> Path:
    path = tmp_path / "transcripts.bed12"
    path.write_text(
        "\n".join(
            [
                "chr1\t100\t400\tTpos\t0\t+\t100\t400\t0\t2\t100,100\t0,200",
                "chr1\t500\t600\tTamb1\t0\t+\t500\t600\t0\t1\t100\t0",
                "chr1\t500\t600\tTamb2\t0\t+\t500\t600\t0\t1\t100\t0",
                "chr2\t100\t400\tTneg\t0\t-\t100\t400\t0\t2\t100,100\t0,200",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return path


@pytest.fixture
def gtf_path(tmp_path: Path) -> Path:
    path = tmp_path / "transcripts.gtf"
    path.write_text(
        "\n".join(
            [
                'chr1\tIsoComp\texon\t101\t200\t.\t+\t.\tgene_id "Gpos"; transcript_id "Tpos";',
                'chr1\tIsoComp\texon\t301\t400\t.\t+\t.\tgene_id "Gpos"; transcript_id "Tpos";',
                'chr1\tIsoComp\texon\t501\t600\t.\t+\t.\tgene_id "Gamb"; transcript_id "Tamb1";',
                'chr1\tIsoComp\texon\t501\t600\t.\t+\t.\tgene_id "Gamb"; transcript_id "Tamb2";',
                'chr2\tIsoComp\texon\t101\t200\t.\t-\t.\tgene_id "Gneg"; transcript_id "Tneg";',
                'chr2\tIsoComp\texon\t301\t400\t.\t-\t.\tgene_id "Gneg"; transcript_id "Tneg";',
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return path


def write_bam(path: Path, records: list[dict[str, object]]) -> Path:
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [
            {"SN": "chr1", "LN": 2000},
            {"SN": "chr2", "LN": 2000},
        ],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for record in records:
            segment = pysam.AlignedSegment()
            segment.query_name = str(record["name"])
            cigar = record["cigar"]
            segment.query_sequence = "A" * _query_length(cigar)
            segment.flag = int(record.get("flag", 0))
            segment.reference_id = bam.get_tid(str(record["chrom"]))
            segment.reference_start = int(record["start"])
            segment.mapping_quality = int(record.get("mapq", 60))
            segment.cigar = cigar
            segment.query_qualities = pysam.qualitystring_to_array("I" * len(segment.query_sequence))
            bam.write(segment)
    return path


def _query_length(cigar: list[tuple[int, int]]) -> int:
    query_consuming = {0, 1, 4, 7, 8}
    return sum(length for op, length in cigar if op in query_consuming)


@pytest.fixture
def synthetic_bam(tmp_path: Path) -> Path:
    return write_bam(
        tmp_path / "reads.bam",
        [
            {"name": "full_pos", "chrom": "chr1", "start": 100, "cigar": [(0, 100), (3, 100), (0, 100)]},
            {"name": "trunc_5p", "chrom": "chr1", "start": 150, "cigar": [(0, 50), (3, 100), (0, 100)]},
            {"name": "trunc_3p", "chrom": "chr1", "start": 100, "cigar": [(0, 100), (3, 100), (0, 50)]},
            {"name": "ambiguous", "chrom": "chr1", "start": 500, "cigar": [(0, 100)]},
            {"name": "low_conf", "chrom": "chr1", "start": 100, "cigar": [(0, 200)]},
            {"name": "unassigned", "chrom": "chr1", "start": 1200, "cigar": [(0, 100)]},
            {"name": "full_neg", "chrom": "chr2", "start": 100, "cigar": [(0, 100), (3, 100), (0, 100)], "flag": 16},
        ],
    )
