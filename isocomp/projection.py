"""Projection of genomic read blocks into transcript coordinates."""

from __future__ import annotations

from .models import ProjectionResult, Transcript
from .utils import merge_intervals, overlap_interval


def project_blocks_to_transcript(
    blocks: list[tuple[int, int]],
    transcript: Transcript,
) -> ProjectionResult:
    projected: list[tuple[int, int]] = []

    for block in blocks:
        for exon in transcript.tx_exons:
            overlap = overlap_interval(block, (exon.genomic_start, exon.genomic_end))
            if overlap is None:
                continue
            projected.append(_genomic_overlap_to_tx(overlap, exon, transcript.strand))

    merged = merge_intervals(projected)
    covered_bases = sum(end - start for start, end in merged)
    if merged and transcript.strand != ".":
        read_start_tx = min(start for start, _ in merged)
        read_end_tx = max(end for _, end in merged)
        dist_to_5p = read_start_tx
        dist_to_3p = transcript.transcript_length - read_end_tx
    else:
        read_start_tx = None
        read_end_tx = None
        dist_to_5p = None
        dist_to_3p = None

    coverage_fraction = covered_bases / transcript.transcript_length if transcript.transcript_length else 0.0
    return ProjectionResult(
        transcript_id=transcript.transcript_id,
        intervals=merged,
        covered_bases=covered_bases,
        read_start_tx=read_start_tx,
        read_end_tx=read_end_tx,
        coverage_fraction=coverage_fraction,
        dist_to_5p=dist_to_5p,
        dist_to_3p=dist_to_3p,
    )


def _genomic_overlap_to_tx(
    overlap: tuple[int, int],
    exon: object,
    strand: str,
) -> tuple[int, int]:
    start, end = overlap
    genomic_start = getattr(exon, "genomic_start")
    genomic_end = getattr(exon, "genomic_end")
    tx_start = getattr(exon, "tx_start")
    if strand == "-":
        return tx_start + (genomic_end - end), tx_start + (genomic_end - start)
    return tx_start + (start - genomic_start), tx_start + (end - genomic_start)
