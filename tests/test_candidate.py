from __future__ import annotations

from isocomp.annotation import parse_bed12_line, read_bed12
from isocomp.candidate import CandidateIndex
from isocomp.models import ReadAlignment


def test_candidate_lookup_filters_by_exonic_overlap(bed12_path) -> None:
    transcripts = read_bed12(bed12_path)
    index = CandidateIndex(transcripts)
    read = ReadAlignment(
        read_id="r1",
        chrom="chr1",
        genomic_start=100,
        genomic_end=400,
        blocks=[(100, 200), (300, 400)],
        junctions=[(200, 300)],
        aligned_length=200,
        mapq=60,
        cigar="100M100N100M",
        is_reverse=False,
    )

    hits = index.query(read, min_overlap=50, strandness="unstranded")

    assert [hit.transcript.transcript_id for hit in hits] == ["Tpos"]
    assert hits[0].exonic_overlap == 200


def test_candidate_lookup_honors_forward_strandness(bed12_path) -> None:
    transcripts = read_bed12(bed12_path)
    index = CandidateIndex(transcripts)
    reverse_read = ReadAlignment(
        read_id="r2",
        chrom="chr2",
        genomic_start=100,
        genomic_end=400,
        blocks=[(100, 200), (300, 400)],
        junctions=[(200, 300)],
        aligned_length=200,
        mapq=60,
        cigar="100M100N100M",
        is_reverse=True,
    )

    hits = index.query(reverse_read, min_overlap=50, strandness="forward")

    assert [hit.transcript.transcript_id for hit in hits] == ["Tneg"]


def test_candidate_lookup_auto_uses_clear_orientation() -> None:
    plus_transcript = parse_bed12_line(
        "chr1\t100\t200\tTplus\t0\t+\t100\t200\t0\t1\t100\t0"
    )
    minus_transcript = parse_bed12_line(
        "chr1\t300\t400\tTminus\t0\t-\t300\t400\t0\t1\t100\t0"
    )
    index = CandidateIndex(
        {
            plus_transcript.transcript_id: plus_transcript,
            minus_transcript.transcript_id: minus_transcript,
        }
    )
    reverse_read = ReadAlignment(
        read_id="r3",
        chrom="chr1",
        genomic_start=100,
        genomic_end=200,
        blocks=[(100, 200)],
        junctions=[],
        aligned_length=100,
        mapq=60,
        cigar="100M",
        is_reverse=True,
    )

    assert index.query(reverse_read, min_overlap=50, strandness="forward") == []
    reverse_hits = index.query(reverse_read, min_overlap=50, strandness="reverse")
    auto_hits = index.query(reverse_read, min_overlap=50, strandness="auto")

    assert [hit.transcript.transcript_id for hit in reverse_hits] == ["Tplus"]
    assert [hit.transcript.transcript_id for hit in auto_hits] == ["Tplus"]
