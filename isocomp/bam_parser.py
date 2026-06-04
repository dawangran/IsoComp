"""Streaming BAM parsing and CIGAR extraction."""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import pysam

from .models import AlignmentFilterStats, ReadAlignment
from .utils import merge_intervals

CMATCH = 0
CINS = 1
CDEL = 2
CREF_SKIP = 3
CSOFT_CLIP = 4
CHARD_CLIP = 5
CEQUAL = 7
CDIFF = 8

REFERENCE_FOOTPRINT_OPS = {CMATCH, CDEL, CEQUAL, CDIFF}
REFERENCE_GAP_OPS = {CREF_SKIP}
REFERENCE_CONSUMING_OPS = {CMATCH, CDEL, CREF_SKIP, CEQUAL, CDIFF}


class BamParserError(ValueError):
    """Raised when BAM input cannot be read or interpreted."""


def iter_read_alignments(
    bam_path: str | Path,
    *,
    min_mapq: int = 20,
    threads: int = 1,
    skip_supplementary: bool = True,
) -> tuple[Iterator[ReadAlignment], AlignmentFilterStats]:
    stats = AlignmentFilterStats()
    iterator = _alignment_iterator(
        Path(bam_path),
        min_mapq=min_mapq,
        threads=threads,
        skip_supplementary=skip_supplementary,
        stats=stats,
    )
    return iterator, stats


def _alignment_iterator(
    bam_path: Path,
    *,
    min_mapq: int,
    threads: int,
    skip_supplementary: bool,
    stats: AlignmentFilterStats,
) -> Iterator[ReadAlignment]:
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM file does not exist: {bam_path}")
    if threads < 1:
        raise BamParserError(f"threads must be >= 1, got {threads}")
    if min_mapq < 0:
        raise BamParserError(f"min_mapq must be >= 0, got {min_mapq}")

    try:
        with pysam.AlignmentFile(str(bam_path), "rb", threads=threads) as bam:
            for read in bam.fetch(until_eof=True):
                stats.total_reads += 1
                if read.is_unmapped:
                    continue
                stats.mapped_reads += 1
                if read.is_secondary or (skip_supplementary and read.is_supplementary):
                    continue
                stats.primary_reads += 1
                if read.is_duplicate:
                    stats.duplicate_reads += 1
                    continue
                if read.mapping_quality < min_mapq:
                    stats.low_mapq_reads += 1
                    continue

                parsed = parse_alignment_record(read)
                if parsed.aligned_length <= 0:
                    stats.empty_alignment_reads += 1
                    continue
                stats.usable_reads += 1
                yield parsed
    except OSError as exc:
        raise BamParserError(f"Could not read BAM file {bam_path}: {exc}") from exc


def parse_alignment_record(read: pysam.AlignedSegment) -> ReadAlignment:
    if read.reference_name is None or read.reference_start is None:
        raise BamParserError(f"Read {read.query_name!r} has no reference name/start")
    if read.cigartuples is None:
        raise BamParserError(f"Read {read.query_name!r} has no CIGAR")

    blocks, junctions = cigar_to_blocks_and_junctions(read.reference_start, read.cigartuples)
    merged_blocks = merge_intervals(blocks)
    if not merged_blocks:
        genomic_start = read.reference_start
        genomic_end = read.reference_start
    else:
        genomic_start = merged_blocks[0][0]
        genomic_end = merged_blocks[-1][1]

    softclip_5p, softclip_3p = softclips_from_cigar(read.cigartuples)
    return ReadAlignment(
        read_id=read.query_name,
        chrom=read.reference_name,
        genomic_start=genomic_start,
        genomic_end=genomic_end,
        blocks=merged_blocks,
        junctions=junctions,
        aligned_length=sum(end - start for start, end in merged_blocks),
        mapq=read.mapping_quality,
        cigar=read.cigarstring or "",
        is_reverse=read.is_reverse,
        softclip_5p=softclip_5p,
        softclip_3p=softclip_3p,
    )


def cigar_to_blocks_and_junctions(
    reference_start: int,
    cigartuples: list[tuple[int, int]],
) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    ref_pos = reference_start
    current_block_start: int | None = None
    blocks: list[tuple[int, int]] = []
    junctions: list[tuple[int, int]] = []

    for op, length in cigartuples:
        if length < 0:
            raise BamParserError(f"Invalid negative CIGAR length: {length}")

        if op in REFERENCE_FOOTPRINT_OPS:
            if current_block_start is None:
                current_block_start = ref_pos
            ref_pos += length
        elif op in REFERENCE_GAP_OPS:
            if current_block_start is not None and ref_pos > current_block_start:
                blocks.append((current_block_start, ref_pos))
            junction_start = ref_pos
            ref_pos += length
            junctions.append((junction_start, ref_pos))
            current_block_start = None
        elif op in REFERENCE_CONSUMING_OPS:
            ref_pos += length
        else:
            continue

    if current_block_start is not None and ref_pos > current_block_start:
        blocks.append((current_block_start, ref_pos))

    return blocks, junctions


def softclips_from_cigar(cigartuples: list[tuple[int, int]]) -> tuple[int, int]:
    left = cigartuples[0][1] if cigartuples and cigartuples[0][0] == CSOFT_CLIP else 0
    right = cigartuples[-1][1] if cigartuples and cigartuples[-1][0] == CSOFT_CLIP else 0
    return left, right

