from __future__ import annotations

import numpy as np

from isocomp.annotation import parse_bed12_line, read_bed12
from isocomp.assigner import assign_read
from isocomp.candidate import CandidateIndex
from isocomp.metrics import build_read_metrics, compute_transcript_body_coverage
from isocomp.models import AssignmentResult, ReadAlignment
from isocomp.projection import project_blocks_to_transcript


def test_full_length_like_classification(bed12_path) -> None:
    transcripts = read_bed12(bed12_path)
    read = ReadAlignment(
        read_id="full_pos",
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
    assignment = assign_read(read, CandidateIndex(transcripts).query(read, min_overlap=50))

    metrics = build_read_metrics(read, assignment, tss_tol=100, tes_tol=100)

    assert metrics.assignment_status == "unique"
    assert metrics.is_5p_complete is True
    assert metrics.is_3p_complete is True
    assert metrics.is_full_length_like is True


def test_ambiguous_reads_do_not_contribute_to_body_coverage(bed12_path) -> None:
    transcripts = read_bed12(bed12_path)
    index = CandidateIndex(transcripts)
    unique_read = ReadAlignment(
        read_id="full_pos",
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
    ambiguous_read = ReadAlignment(
        read_id="ambiguous",
        chrom="chr1",
        genomic_start=500,
        genomic_end=600,
        blocks=[(500, 600)],
        junctions=[],
        aligned_length=100,
        mapq=60,
        cigar="100M",
        is_reverse=False,
    )
    unique_assignment = assign_read(unique_read, index.query(unique_read, min_overlap=50))
    ambiguous_assignment = assign_read(ambiguous_read, index.query(ambiguous_read, min_overlap=50))

    aggregate, per_transcript = compute_transcript_body_coverage(
        transcripts,
        [unique_assignment, ambiguous_assignment],
        bin_num=10,
    )

    assert unique_assignment.status == "unique"
    assert ambiguous_assignment.status == "ambiguous"
    assert np.allclose(per_transcript["Tpos"], np.ones(10))
    assert np.allclose(per_transcript["Tamb1"], np.zeros(10))
    assert np.allclose(per_transcript["Tamb2"], np.zeros(10))
    assert np.allclose(aggregate, np.ones(10))


def test_negative_body_coverage_bins_are_oriented_5p_to_3p() -> None:
    transcript = parse_bed12_line("chr2\t100\t400\tTneg\t0\t-\t100\t400\t0\t2\t100,100\t0,200")
    five_prime_assignment = AssignmentResult(
        read_id="neg_5p",
        status="unique",
        transcript=transcript,
        score=1.0,
        second_best_transcript=None,
        second_best_score=None,
        exon_overlap_score=1.0,
        junction_match_count=0,
        junction_precision=1.0,
        junction_recall=1.0,
        projection=project_blocks_to_transcript([(350, 400)], transcript),
    )
    three_prime_assignment = AssignmentResult(
        read_id="neg_3p",
        status="unique",
        transcript=transcript,
        score=1.0,
        second_best_transcript=None,
        second_best_score=None,
        exon_overlap_score=1.0,
        junction_match_count=0,
        junction_precision=1.0,
        junction_recall=1.0,
        projection=project_blocks_to_transcript([(100, 150)], transcript),
    )

    _, five_prime_coverage = compute_transcript_body_coverage(
        {"Tneg": transcript},
        [five_prime_assignment],
        bin_num=4,
    )
    _, three_prime_coverage = compute_transcript_body_coverage(
        {"Tneg": transcript},
        [three_prime_assignment],
        bin_num=4,
    )

    assert np.allclose(five_prime_coverage["Tneg"], [1.0, 0.0, 0.0, 0.0])
    assert np.allclose(three_prime_coverage["Tneg"], [0.0, 0.0, 0.0, 1.0])
