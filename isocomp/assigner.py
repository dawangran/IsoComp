"""Read-to-transcript compatibility scoring and assignment."""

from __future__ import annotations

from dataclasses import dataclass

from .candidate import CandidateHit
from .models import AssignmentResult, ProjectionResult, ReadAlignment, Transcript
from .projection import project_blocks_to_transcript
from .utils import harmonic_mean, safe_divide, total_overlap_length


@dataclass(frozen=True)
class CandidateScore:
    transcript: Transcript
    final_score: float
    exon_overlap_score: float
    junction_match_count: int
    junction_precision: float
    junction_recall: float
    junction_score: float
    coverage_fraction: float
    projection: ProjectionResult


def assign_read(
    read: ReadAlignment,
    candidates: list[CandidateHit],
    *,
    junction_tol: int = 5,
    unique_threshold: float = 0.8,
    margin_threshold: float = 0.1,
) -> AssignmentResult:
    if not candidates:
        return AssignmentResult(
            read_id=read.read_id,
            status="unassigned",
            transcript=None,
            score=0.0,
            second_best_transcript=None,
            second_best_score=None,
            exon_overlap_score=0.0,
            junction_match_count=0,
            junction_precision=0.0,
            junction_recall=0.0,
            projection=None,
        )

    scored = [
        score_candidate(read, hit.transcript, junction_tol=junction_tol)
        for hit in candidates
    ]
    scored.sort(
        key=lambda item: (
            -item.final_score,
            -item.exon_overlap_score,
            -item.coverage_fraction,
            item.transcript.transcript_id,
        )
    )

    top = scored[0]
    second = scored[1] if len(scored) > 1 else None
    second_score = second.final_score if second else None
    second_id = second.transcript.transcript_id if second else None
    margin = top.final_score - (second_score if second_score is not None else 0.0)

    if top.final_score >= unique_threshold and margin >= margin_threshold:
        status = "unique"
    elif top.final_score >= unique_threshold:
        status = "ambiguous"
    else:
        status = "low_confidence"

    return AssignmentResult(
        read_id=read.read_id,
        status=status,
        transcript=top.transcript,
        score=top.final_score,
        second_best_transcript=second_id,
        second_best_score=second_score,
        exon_overlap_score=top.exon_overlap_score,
        junction_match_count=top.junction_match_count,
        junction_precision=top.junction_precision,
        junction_recall=top.junction_recall,
        projection=top.projection,
    )


def score_candidate(
    read: ReadAlignment,
    transcript: Transcript,
    *,
    junction_tol: int = 5,
) -> CandidateScore:
    exonic_overlap = total_overlap_length(read.blocks, transcript.exons)
    exon_overlap_score = safe_divide(exonic_overlap, read.aligned_length)
    projection = project_blocks_to_transcript(read.blocks, transcript)
    (
        junction_match_count,
        junction_precision,
        junction_recall,
    ) = score_junctions(read, transcript, junction_tol=junction_tol)
    junction_score = harmonic_mean(junction_precision, junction_recall)
    final_score = (
        0.50 * exon_overlap_score
        + 0.30 * junction_score
        + 0.20 * projection.coverage_fraction
    )
    return CandidateScore(
        transcript=transcript,
        final_score=final_score,
        exon_overlap_score=exon_overlap_score,
        junction_match_count=junction_match_count,
        junction_precision=junction_precision,
        junction_recall=junction_recall,
        junction_score=junction_score,
        coverage_fraction=projection.coverage_fraction,
        projection=projection,
    )


def score_junctions(
    read: ReadAlignment,
    transcript: Transcript,
    *,
    junction_tol: int = 5,
) -> tuple[int, float, float]:
    if not read.junctions:
        transcript_junctions_in_span = _transcript_junctions_in_read_span(read, transcript)
        recall = 0.0 if transcript_junctions_in_span else 1.0
        return 0, 1.0, recall

    matched = 0
    used_transcript_junctions: set[int] = set()
    for read_junction in read.junctions:
        match_index = _find_matching_junction(read_junction, transcript.junctions, junction_tol)
        if match_index is not None:
            matched += 1
            used_transcript_junctions.add(match_index)

    transcript_junctions_in_span = _transcript_junctions_in_read_span(read, transcript)
    precision = safe_divide(matched, len(read.junctions))
    recall = safe_divide(
        len(used_transcript_junctions),
        len(transcript_junctions_in_span),
    )
    return matched, precision, recall


def _find_matching_junction(
    read_junction: tuple[int, int],
    transcript_junctions: list[tuple[int, int]],
    junction_tol: int,
) -> int | None:
    for index, transcript_junction in enumerate(transcript_junctions):
        if (
            abs(read_junction[0] - transcript_junction[0]) <= junction_tol
            and abs(read_junction[1] - transcript_junction[1]) <= junction_tol
        ):
            return index
    return None


def _transcript_junctions_in_read_span(
    read: ReadAlignment,
    transcript: Transcript,
) -> list[tuple[int, int]]:
    read_start = read.genomic_start
    read_end = read.genomic_end
    return [
        junction
        for junction in transcript.junctions
        if junction[0] >= read_start and junction[1] <= read_end
    ]
