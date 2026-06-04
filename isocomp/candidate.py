"""Candidate transcript interval lookup."""

from __future__ import annotations

from dataclasses import dataclass

from intervaltree import IntervalTree

from .models import ReadAlignment, Transcript
from .utils import total_overlap_length


@dataclass(frozen=True)
class CandidateHit:
    transcript: Transcript
    exonic_overlap: int


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

        hits: list[CandidateHit] = []
        for interval in sorted(tree.overlap(read.genomic_start, read.genomic_end)):
            transcript = interval.data
            if not _strand_compatible(read, transcript, strandness):
                continue
            exonic_overlap = total_overlap_length(read.blocks, transcript.exons)
            if exonic_overlap >= min_overlap:
                hits.append(CandidateHit(transcript=transcript, exonic_overlap=exonic_overlap))

        hits.sort(key=lambda hit: (-hit.exonic_overlap, hit.transcript.transcript_id))
        return hits


def _strand_compatible(read: ReadAlignment, transcript: Transcript, strandness: str) -> bool:
    if strandness == "unstranded" or transcript.strand == ".":
        return True
    if strandness == "forward":
        read_strand = "-" if read.is_reverse else "+"
    elif strandness == "reverse":
        read_strand = "+" if read.is_reverse else "-"
    elif strandness == "auto":
        return True
    else:
        raise ValueError(
            "strandness must be one of unstranded, forward, reverse, auto; "
            f"got {strandness!r}"
        )
    return read_strand == transcript.strand

