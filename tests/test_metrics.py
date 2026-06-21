from __future__ import annotations

import numpy as np

from isocomp.annotation import parse_bed12_line, read_bed12
from isocomp.assigner import assign_read
from isocomp.candidate import CandidateHit, CandidateIndex
from isocomp.metrics import (
    OnlineNumericSummary,
    SampleMetricAccumulator,
    TranscriptMetricAccumulator,
    build_read_metrics,
    compute_transcript_body_coverage,
    summarize_sample,
    summarize_sample_accumulator,
    summarize_transcripts,
    update_summary_for_metric,
    update_sample_accumulator,
)
from isocomp.models import AssignmentResult, ReadAlignment, RunSummary
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


def test_online_numeric_summary_is_exact_for_small_samples() -> None:
    summary = OnlineNumericSummary()
    for value in [5, 1, 3, 2]:
        summary.add(value)

    assert summary.mean() == 2.75
    assert summary.median() == 2.5


def test_online_numeric_summary_does_not_repopulate_exact_values_after_limit() -> None:
    summary = OnlineNumericSummary(max_exact_values=5)
    for value in [0, 0, 0, 1000, 1000, 1000, 1000]:
        summary.add(value)

    assert summary.exact_values == []
    assert summary.median() == summary.median_estimator.median()


def test_transcript_numeric_summaries_use_constant_memory() -> None:
    accumulator = TranscriptMetricAccumulator()
    for value in range(20):
        accumulator.coverage_values.add(value / 20)
        accumulator.dist_5p.add(value)
        accumulator.dist_3p.add(value)

    for summary in [
        accumulator.coverage_values,
        accumulator.dist_5p,
        accumulator.dist_3p,
    ]:
        assert summary.max_exact_values == 0
        assert summary.exact_values == []
        assert len(summary.median_estimator.initial_values) == 5


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
    assert "Tamb1" not in per_transcript
    assert "Tamb2" not in per_transcript
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


def test_low_confidence_read_is_not_counted_as_assigned() -> None:
    transcript = parse_bed12_line("chr1\t100\t400\tTpos\t0\t+\t100\t400\t0\t2\t100,100\t0,200")
    read = ReadAlignment(
        read_id="low_conf",
        chrom="chr1",
        genomic_start=100,
        genomic_end=300,
        blocks=[(100, 300)],
        junctions=[],
        aligned_length=200,
        mapq=60,
        cigar="200M",
        is_reverse=False,
    )
    assignment = assign_read(read, [CandidateHit(transcript, 100)])
    metric = build_read_metrics(read, assignment, tss_tol=100, tes_tol=100)
    summary = RunSummary()
    update_summary_for_metric(summary, metric)
    transcript_rows = summarize_transcripts(
        {transcript.transcript_id: transcript},
        [metric],
        {transcript.transcript_id: np.zeros(10)},
    )
    sample_row = summarize_sample([metric], summary, sample_name="sample")

    assert assignment.status == "low_confidence"
    assert assignment.transcript is None
    assert metric.transcript_id is None
    assert summary.assigned_reads == 0
    assert summary.low_confidence_reads == 1
    assert transcript_rows[0]["assigned_read_count"] == 0
    assert sample_row["assigned_reads"] == 0
    assert sample_row["low_confidence_reads"] == 1


def test_negative_strand_metrics_report_transcript_terminal_softclips() -> None:
    transcript = parse_bed12_line("chr2\t100\t400\tTneg\t0\t-\t100\t400\t0\t2\t100,100\t0,200")
    read = ReadAlignment(
        read_id="neg_softclip",
        chrom="chr2",
        genomic_start=100,
        genomic_end=400,
        blocks=[(100, 200), (300, 400)],
        junctions=[(200, 300)],
        aligned_length=200,
        mapq=60,
        cigar="5S100M100N100M10S",
        is_reverse=True,
        softclip_5p=10,
        softclip_3p=5,
        softclip_left=5,
        softclip_right=10,
    )
    assignment = AssignmentResult(
        read_id=read.read_id,
        status="unique",
        transcript=transcript,
        score=1.0,
        second_best_transcript=None,
        second_best_score=None,
        exon_overlap_score=1.0,
        junction_match_count=1,
        junction_precision=1.0,
        junction_recall=1.0,
        projection=project_blocks_to_transcript(read.blocks, transcript),
    )

    metric = build_read_metrics(read, assignment, tss_tol=100, tes_tol=100)

    assert metric.softclip_5p == 10
    assert metric.softclip_3p == 5


def test_unknown_strand_read_is_excluded_from_terminal_metrics_and_body_coverage() -> None:
    transcript = parse_bed12_line(
        "chr3\t100\t200\tTunknown\t0\t.\t100\t200\t0\t1\t100\t0"
    )
    read = ReadAlignment(
        read_id="unknown",
        chrom="chr3",
        genomic_start=100,
        genomic_end=200,
        blocks=[(100, 200)],
        junctions=[],
        aligned_length=100,
        mapq=60,
        cigar="100M",
        is_reverse=False,
    )
    assignment = assign_read(read, [CandidateHit(transcript, 100)])
    metric = build_read_metrics(read, assignment, tss_tol=100, tes_tol=100)
    aggregate, per_transcript = compute_transcript_body_coverage(
        {transcript.transcript_id: transcript},
        [assignment],
        bin_num=10,
    )

    assert assignment.status == "unique"
    assert metric.coverage_fraction == 1.0
    assert metric.is_5p_complete is None
    assert metric.is_3p_complete is None
    assert metric.is_full_length_like is False
    assert not np.any(aggregate)
    assert per_transcript == {}


def test_unknown_strand_read_does_not_dilute_terminal_fractions() -> None:
    known = parse_bed12_line(
        "chr1\t100\t200\tTknown\t0\t+\t100\t200\t0\t1\t100\t0"
    )
    unknown = parse_bed12_line(
        "chr2\t100\t200\tTunknown\t0\t.\t100\t200\t0\t1\t100\t0"
    )
    accumulator = SampleMetricAccumulator()
    for transcript in [known, unknown]:
        read = ReadAlignment(
            read_id=transcript.transcript_id,
            chrom=transcript.chrom,
            genomic_start=100,
            genomic_end=200,
            blocks=[(100, 200)],
            junctions=[],
            aligned_length=100,
            mapq=60,
            cigar="100M",
            is_reverse=False,
        )
        assignment = assign_read(read, [CandidateHit(transcript, 100)])
        metric = build_read_metrics(read, assignment, tss_tol=100, tes_tol=100)
        update_sample_accumulator(accumulator, metric)

    row = summarize_sample_accumulator(
        accumulator,
        RunSummary(unique_assigned_reads=2, assigned_reads=2),
        sample_name="sample",
    )

    assert row["terminal_evaluable_unique_reads"] == 1
    assert row["5p_complete_fraction"] == 1.0
    assert row["3p_complete_fraction"] == 1.0
    assert row["full_length_like_fraction"] == 1.0


def test_full_length_coverage_threshold_is_configurable(bed12_path) -> None:
    transcript = read_bed12(bed12_path)["Tpos"]
    read = ReadAlignment(
        read_id="truncated_but_terminal",
        chrom="chr1",
        genomic_start=150,
        genomic_end=400,
        blocks=[(150, 200), (300, 400)],
        junctions=[(200, 300)],
        aligned_length=150,
        mapq=60,
        cigar="50M100N100M",
        is_reverse=False,
    )
    assignment = assign_read(read, [CandidateHit(transcript, 150)])

    strict = build_read_metrics(
        read,
        assignment,
        tss_tol=100,
        tes_tol=100,
        full_length_coverage=0.8,
    )
    relaxed = build_read_metrics(
        read,
        assignment,
        tss_tol=100,
        tes_tol=100,
        full_length_coverage=0.7,
    )

    assert strict.coverage_fraction == 0.75
    assert strict.is_full_length_like is False
    assert relaxed.is_full_length_like is True


def test_terminal_anchor_rejects_tiny_spurious_terminal_overlap() -> None:
    transcript = parse_bed12_line(
        "chr1\t100\t1100\tTsingle\t0\t+\t100\t1100\t0\t1\t1000\t0"
    )
    read = ReadAlignment(
        read_id="tiny_terminal_anchor",
        chrom="chr1",
        genomic_start=100,
        genomic_end=800,
        blocks=[(100, 101), (600, 800)],
        junctions=[],
        aligned_length=201,
        mapq=60,
        cigar="1M499D200M",
        is_reverse=False,
    )
    assignment = assign_read(read, [CandidateHit(transcript, 201)])

    guarded = build_read_metrics(
        read,
        assignment,
        tss_tol=100,
        tes_tol=100,
        min_terminal_anchor=10,
    )
    legacy = build_read_metrics(
        read,
        assignment,
        tss_tol=100,
        tes_tol=100,
        min_terminal_anchor=0,
    )

    assert guarded.dist_to_5p == 0
    assert guarded.terminal_anchor_5p == 1
    assert guarded.is_5p_complete is False
    assert legacy.is_5p_complete is True
