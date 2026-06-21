"""Candidate transcript interval lookup."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

try:
    from intervaltree import IntervalTree
except ModuleNotFoundError:
    @dataclass(frozen=True)
    class _FallbackInterval:
        begin: int
        end: int
        data: object

        def __lt__(self, other: object) -> bool:
            if not isinstance(other, _FallbackInterval):
                return NotImplemented
            return (self.begin, self.end, repr(self.data)) < (
                other.begin,
                other.end,
                repr(other.data),
            )

    class IntervalTree:  # type: ignore[no-redef]
        """Small fallback used when intervaltree is unavailable.

        The real dependency is preferred for production-size annotations. This
        linear fallback keeps tiny synthetic tests and examples runnable in
        minimal environments.
        """

        def __init__(self) -> None:
            self._intervals: list[_FallbackInterval] = []

        def addi(self, start: int, end: int, data: object) -> None:
            self._intervals.append(_FallbackInterval(start, end, data))

        def overlap(self, start: int, end: int) -> list[_FallbackInterval]:
            return [
                interval
                for interval in self._intervals
                if interval.begin < end and start < interval.end
            ]

from .models import ReadAlignment, Transcript
from .utils import total_overlap_length


@dataclass(frozen=True)
class CandidateHit:
    transcript: Transcript
    exonic_overlap: int


@dataclass(frozen=True)
class StrandInferenceResult:
    resolved_strandness: str
    scanned_reads: int
    informative_reads: int
    forward_votes: int
    reverse_votes: int
    dominant_fraction: float


class CandidateIndex:
    def __init__(self, transcripts: dict[str, Transcript]) -> None:
        self.transcripts = transcripts
        self.by_chrom: dict[str, IntervalTree] = {}
        for transcript in transcripts.values():
            start, end = transcript.span
            self.by_chrom.setdefault(transcript.chrom, IntervalTree()).addi(start, end, transcript)

    def query(
        self,
        read: ReadAlignment,
        *,
        min_overlap: int = 50,
        strandness: str = "unstranded",
    ) -> list[CandidateHit]:
        tree = self.by_chrom.get(read.chrom)
        if tree is None:
            return []
        if strandness not in {"unstranded", "forward", "reverse", "auto"}:
            raise ValueError(
                "strandness must be one of unstranded, forward, reverse, auto; "
                f"got {strandness!r}"
            )

        raw_hits: list[CandidateHit] = []
        for interval in sorted(tree.overlap(read.genomic_start, read.genomic_end)):
            transcript = interval.data
            exonic_overlap = total_overlap_length(read.blocks, transcript.exons)
            if exonic_overlap >= min_overlap:
                raw_hits.append(CandidateHit(transcript=transcript, exonic_overlap=exonic_overlap))

        if strandness == "auto":
            hits = _auto_strand_hits(read, raw_hits)
        else:
            hits = [
                hit
                for hit in raw_hits
                if _strand_compatible(read, hit.transcript, strandness)
            ]

        hits.sort(key=lambda hit: (-hit.exonic_overlap, hit.transcript.transcript_id))
        return hits


def infer_library_strandness(
    reads: Iterable[ReadAlignment],
    candidate_index: CandidateIndex,
    *,
    min_overlap: int = 50,
    max_reads: int = 100_000,
    min_informative_reads: int = 100,
    min_dominant_fraction: float = 0.8,
) -> StrandInferenceResult:
    if max_reads < 1:
        raise ValueError("max_reads must be >= 1")
    if min_informative_reads < 1:
        raise ValueError("min_informative_reads must be >= 1")
    if not 0.5 <= min_dominant_fraction <= 1:
        raise ValueError("min_dominant_fraction must be between 0.5 and 1")

    scanned_reads = 0
    forward_votes = 0
    reverse_votes = 0
    for read in reads:
        scanned_reads += 1
        hits = candidate_index.query(
            read,
            min_overlap=min_overlap,
            strandness="unstranded",
        )
        candidate_strands = {
            hit.transcript.strand
            for hit in hits
            if hit.transcript.strand in {"+", "-"}
        }
        if len(candidate_strands) == 1:
            transcript_strand = next(iter(candidate_strands))
            read_strand = "-" if read.is_reverse else "+"
            if read_strand == transcript_strand:
                forward_votes += 1
            else:
                reverse_votes += 1
        if scanned_reads >= max_reads:
            break

    informative_reads = forward_votes + reverse_votes
    dominant_votes = max(forward_votes, reverse_votes)
    dominant_fraction = (
        dominant_votes / informative_reads if informative_reads else 0.0
    )
    resolved = "unstranded"
    if (
        informative_reads >= min_informative_reads
        and dominant_fraction >= min_dominant_fraction
    ):
        resolved = "forward" if forward_votes > reverse_votes else "reverse"

    return StrandInferenceResult(
        resolved_strandness=resolved,
        scanned_reads=scanned_reads,
        informative_reads=informative_reads,
        forward_votes=forward_votes,
        reverse_votes=reverse_votes,
        dominant_fraction=dominant_fraction,
    )


def _auto_strand_hits(read: ReadAlignment, hits: list[CandidateHit]) -> list[CandidateHit]:
    forward_hits = [
        hit
        for hit in hits
        if hit.transcript.strand != "." and _strand_compatible(read, hit.transcript, "forward")
    ]
    reverse_hits = [
        hit
        for hit in hits
        if hit.transcript.strand != "." and _strand_compatible(read, hit.transcript, "reverse")
    ]
    unknown_hits = [hit for hit in hits if hit.transcript.strand == "."]

    if forward_hits and not reverse_hits:
        return forward_hits + unknown_hits
    if reverse_hits and not forward_hits:
        return reverse_hits + unknown_hits
    return hits


def _strand_compatible(read: ReadAlignment, transcript: Transcript, strandness: str) -> bool:
    if strandness == "unstranded" or transcript.strand == ".":
        return True
    if strandness == "forward":
        read_strand = "-" if read.is_reverse else "+"
    elif strandness == "reverse":
        read_strand = "+" if read.is_reverse else "-"
    else:
        raise ValueError(
            "strandness must be one of unstranded, forward, reverse, auto; "
            f"got {strandness!r}"
        )
    return read_strand == transcript.strand
