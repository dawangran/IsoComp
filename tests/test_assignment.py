from __future__ import annotations

from isocomp.annotation import read_bed12
from isocomp.assigner import assign_read, score_candidate, score_junctions
from isocomp.candidate import CandidateHit, CandidateIndex
from isocomp.models import ReadAlignment


def test_full_length_read_assigns_unique(bed12_path) -> None:
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

    assignment = assign_read(read, [CandidateHit(transcripts["Tpos"], 200)])

    assert assignment.status == "unique"
    assert assignment.score == 1.0
    assert assignment.transcript is transcripts["Tpos"]


def test_ambiguous_identical_single_exon_transcripts(bed12_path) -> None:
    transcripts = read_bed12(bed12_path)
    read = ReadAlignment(
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
    hits = CandidateIndex(transcripts).query(read, min_overlap=50)

    assignment = assign_read(read, hits)

    assert assignment.status == "ambiguous"
    assert assignment.score == 1.0
    assert assignment.second_best_score == 1.0


def test_junction_scoring_requires_both_splice_sites(bed12_path) -> None:
    transcript = read_bed12(bed12_path)["Tpos"]
    read = ReadAlignment(
        read_id="bad_junction",
        chrom="chr1",
        genomic_start=100,
        genomic_end=400,
        blocks=[(100, 200), (301, 400)],
        junctions=[(200, 301)],
        aligned_length=199,
        mapq=60,
        cigar="100M101N99M",
        is_reverse=False,
    )

    matched, precision, recall = score_junctions(read, transcript, junction_tol=0)

    assert matched == 0
    assert precision == 0.0
    assert recall == 0.0


def test_low_confidence_when_exon_chain_is_poor(bed12_path) -> None:
    transcript = read_bed12(bed12_path)["Tpos"]
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
    score = score_candidate(read, transcript)
    assignment = assign_read(read, [CandidateHit(transcript, 100)])

    assert score.final_score < 0.8
    assert assignment.status == "low_confidence"


def test_assignment_score_does_not_penalize_truncated_but_compatible_read(bed12_path) -> None:
    transcript = read_bed12(bed12_path)["Tpos"]
    read = ReadAlignment(
        read_id="trunc_5p",
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

    score = score_candidate(read, transcript)
    assignment = assign_read(read, [CandidateHit(transcript, 150)])

    assert score.coverage_fraction == 0.75
    assert score.final_score == 1.0
    assert assignment.status == "unique"
