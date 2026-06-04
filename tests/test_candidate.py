from __future__ import annotations

from isocomp.annotation import read_bed12
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

