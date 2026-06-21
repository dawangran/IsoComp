"""Core dataclasses shared by the IsoComp pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

Interval = tuple[int, int]
AssignmentStatus = Literal["unique", "ambiguous", "low_confidence", "unassigned"]


@dataclass(frozen=True)
class TranscriptExon:
    """One genomic exon and its transcript-coordinate interval."""

    genomic_start: int
    genomic_end: int
    tx_start: int
    tx_end: int


@dataclass(frozen=True)
class Transcript:
    transcript_id: str
    gene_id: str | None
    chrom: str
    strand: str
    exons: list[Interval]
    transcript_length: int = field(init=False)
    junctions: list[Interval] = field(init=False)
    tx_exons: list[TranscriptExon] = field(init=False)
    span: Interval = field(init=False)

    def __post_init__(self) -> None:
        if self.strand not in {"+", "-", "."}:
            raise ValueError(f"Invalid transcript strand for {self.transcript_id}: {self.strand}")
        if not self.exons:
            raise ValueError(f"Transcript {self.transcript_id} has no exons")

        exons = sorted(self.exons)
        for start, end in exons:
            if start < 0 or end <= start:
                raise ValueError(f"Invalid exon for {self.transcript_id}: {start}-{end}")
        for previous, current in zip(exons, exons[1:]):
            if current[0] < previous[1]:
                raise ValueError(
                    f"Overlapping exons for {self.transcript_id}: "
                    f"{previous[0]}-{previous[1]} and {current[0]}-{current[1]}"
                )

        transcript_length = sum(end - start for start, end in exons)
        junctions = [
            (left_end, right_start)
            for (_, left_end), (right_start, _) in zip(exons, exons[1:])
            if right_start > left_end
        ]

        ordered_exons = exons if self.strand != "-" else list(reversed(exons))
        tx_exons: list[TranscriptExon] = []
        tx_pos = 0
        for genomic_start, genomic_end in ordered_exons:
            exon_len = genomic_end - genomic_start
            tx_exons.append(
                TranscriptExon(
                    genomic_start=genomic_start,
                    genomic_end=genomic_end,
                    tx_start=tx_pos,
                    tx_end=tx_pos + exon_len,
                )
            )
            tx_pos += exon_len

        object.__setattr__(self, "exons", exons)
        object.__setattr__(self, "transcript_length", transcript_length)
        object.__setattr__(self, "junctions", junctions)
        object.__setattr__(self, "tx_exons", tx_exons)
        object.__setattr__(self, "span", (exons[0][0], exons[-1][1]))


@dataclass(frozen=True)
class ReadAlignment:
    read_id: str
    chrom: str
    genomic_start: int
    genomic_end: int
    blocks: list[Interval]
    junctions: list[Interval]
    aligned_length: int
    mapq: int
    cigar: str
    is_reverse: bool
    softclip_5p: int = 0
    softclip_3p: int = 0
    softclip_left: int = 0
    softclip_right: int = 0


@dataclass(frozen=True)
class ProjectionResult:
    transcript_id: str
    intervals: list[Interval]
    covered_bases: int
    read_start_tx: int | None
    read_end_tx: int | None
    coverage_fraction: float
    dist_to_5p: int | None
    dist_to_3p: int | None


@dataclass(frozen=True)
class AssignmentResult:
    read_id: str
    status: AssignmentStatus
    transcript: Transcript | None
    score: float
    second_best_transcript: str | None
    second_best_score: float | None
    exon_overlap_score: float
    junction_match_count: int
    junction_precision: float
    junction_recall: float
    projection: ProjectionResult | None


@dataclass(frozen=True)
class ReadMetrics:
    read_id: str
    chrom: str
    gene_id: str | None
    transcript_id: str | None
    assignment_status: AssignmentStatus
    assignment_score: float
    second_best_transcript: str | None
    second_best_score: float | None
    transcript_length: int | None
    read_aligned_length: int
    read_start_tx: int | None
    read_end_tx: int | None
    start_pct: float | None
    end_pct: float | None
    coverage_fraction: float | None
    dist_to_5p: int | None
    dist_to_3p: int | None
    terminal_anchor_5p: int | None
    terminal_anchor_3p: int | None
    is_5p_complete: bool | None
    is_3p_complete: bool | None
    is_full_length_like: bool
    junction_match_count: int
    junction_precision: float
    junction_recall: float
    exon_overlap_score: float
    mapq: int
    softclip_5p: int
    softclip_3p: int


@dataclass
class AlignmentFilterStats:
    total_reads: int = 0
    mapped_reads: int = 0
    primary_reads: int = 0
    duplicate_reads: int = 0
    low_mapq_reads: int = 0
    empty_alignment_reads: int = 0
    usable_reads: int = 0


@dataclass
class RunSummary:
    total_reads: int = 0
    mapped_reads: int = 0
    primary_reads: int = 0
    duplicate_reads: int = 0
    low_mapq_reads: int = 0
    empty_alignment_reads: int = 0
    usable_reads: int = 0
    assigned_reads: int = 0
    unique_assigned_reads: int = 0
    ambiguous_reads: int = 0
    gene_only_reads: int = 0
    low_confidence_reads: int = 0
    unassigned_reads: int = 0
