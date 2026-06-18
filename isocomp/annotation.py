"""BED12 annotation parsing."""

from __future__ import annotations

from pathlib import Path

from .models import Transcript


class AnnotationError(ValueError):
    """Raised when an annotation file cannot be parsed safely."""


def read_bed12(path: str | Path) -> dict[str, Transcript]:
    annotation_path = Path(path)
    if not annotation_path.exists():
        raise FileNotFoundError(f"Annotation file does not exist: {annotation_path}")
    if not annotation_path.is_file():
        raise AnnotationError(f"Annotation path is not a regular file: {annotation_path}")

    transcripts: dict[str, Transcript] = {}
    with annotation_path.open("rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            transcript = parse_bed12_line(stripped, line_number)
            if transcript.transcript_id in transcripts:
                raise AnnotationError(
                    f"Duplicate transcript_id {transcript.transcript_id!r} at line {line_number}"
                )
            transcripts[transcript.transcript_id] = transcript

    if not transcripts:
        raise AnnotationError(f"Annotation contains no BED12 transcript records: {annotation_path}")
    return transcripts


def parse_bed12_line(line: str, line_number: int = 1) -> Transcript:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 12:
        raise AnnotationError(
            f"BED12 line {line_number} has {len(fields)} columns; expected at least 12"
        )

    chrom = fields[0]
    name = fields[3]
    strand = fields[5]
    if not chrom:
        raise AnnotationError(f"BED12 line {line_number} has an empty chromosome")
    if not name:
        raise AnnotationError(f"BED12 line {line_number} has an empty transcript name")
    if strand not in {"+", "-", "."}:
        raise AnnotationError(f"BED12 line {line_number} has invalid strand: {strand!r}")

    try:
        chrom_start = int(fields[1])
        chrom_end = int(fields[2])
        block_count = int(fields[9])
        block_sizes = _parse_int_list(fields[10])
        block_starts = _parse_int_list(fields[11])
    except ValueError as exc:
        raise AnnotationError(f"BED12 line {line_number} has non-integer coordinate fields") from exc

    if chrom_start < 0 or chrom_end <= chrom_start:
        raise AnnotationError(
            f"BED12 line {line_number} has invalid chromStart/chromEnd: {chrom_start}-{chrom_end}"
        )
    if block_count <= 0:
        raise AnnotationError(f"BED12 line {line_number} has non-positive blockCount")
    if block_count != len(block_sizes) or block_count != len(block_starts):
        raise AnnotationError(
            f"BED12 line {line_number} blockCount does not match blockSizes/blockStarts"
        )

    exons: list[tuple[int, int]] = []
    for size, relative_start in zip(block_sizes, block_starts):
        if size <= 0 or relative_start < 0:
            raise AnnotationError(f"BED12 line {line_number} has invalid block size/start")
        exon_start = chrom_start + relative_start
        exon_end = exon_start + size
        if exon_end > chrom_end:
            raise AnnotationError(f"BED12 line {line_number} exon extends beyond chromEnd")
        exons.append((exon_start, exon_end))

    try:
        return Transcript(
            transcript_id=name,
            gene_id=None,
            chrom=chrom,
            strand=strand,
            exons=exons,
        )
    except ValueError as exc:
        raise AnnotationError(f"BED12 line {line_number} has invalid exon structure: {exc}") from exc


def _parse_int_list(value: str) -> list[int]:
    return [int(item) for item in value.rstrip(",").split(",") if item != ""]
