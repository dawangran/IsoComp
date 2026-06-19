"""BED12 and GTF annotation parsing."""

from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO

from .models import Transcript

SUPPORTED_ANNOTATION_FORMATS = {"auto", "bed12", "gtf"}


class AnnotationError(ValueError):
    """Raised when an annotation file cannot be parsed safely."""


@dataclass
class _GtfTranscriptRecord:
    transcript_id: str
    gene_id: str | None
    chrom: str
    strand: str
    exons: list[tuple[int, int]] = field(default_factory=list)


def read_annotation(
    path: str | Path,
    annotation_format: str = "auto",
) -> dict[str, Transcript]:
    normalized_format = annotation_format.lower()
    if normalized_format not in SUPPORTED_ANNOTATION_FORMATS:
        raise AnnotationError(
            "annotation format must be one of auto, bed12, gtf; "
            f"got {annotation_format!r}"
        )

    annotation_path = _ensure_annotation_path(path)
    if normalized_format == "auto":
        normalized_format = infer_annotation_format(annotation_path)
    if normalized_format == "bed12":
        return read_bed12(annotation_path)
    if normalized_format == "gtf":
        return read_gtf(annotation_path)
    raise AnnotationError(f"Unsupported annotation format: {annotation_format!r}")


def infer_annotation_format(path: str | Path) -> str:
    annotation_path = Path(path)
    suffixes = {suffix.lower() for suffix in annotation_path.suffixes}
    if ".gtf" in suffixes:
        return "gtf"
    if ".bed" in suffixes or ".bed12" in suffixes:
        return "bed12"

    with _open_text(annotation_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split("\t")
            if len(fields) >= 12 and _is_int(fields[1]) and _is_int(fields[2]):
                return "bed12"
            if len(fields) >= 9 and _is_int(fields[3]) and _is_int(fields[4]):
                return "gtf"
            raise AnnotationError(
                f"Could not infer annotation format from line {line_number}: "
                "expected BED12 or GTF-like columns"
            )

    raise AnnotationError(f"Annotation contains no records: {annotation_path}")


def read_bed12(path: str | Path) -> dict[str, Transcript]:
    annotation_path = _ensure_annotation_path(path)

    transcripts: dict[str, Transcript] = {}
    with _open_text(annotation_path) as handle:
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


def read_gtf(path: str | Path) -> dict[str, Transcript]:
    annotation_path = _ensure_annotation_path(path)
    records: dict[str, _GtfTranscriptRecord] = {}

    with _open_text(annotation_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split("\t")
            if len(fields) < 9:
                raise AnnotationError(
                    f"GTF line {line_number} has {len(fields)} columns; expected 9"
                )
            if fields[2] != "exon":
                continue

            transcript_id, gene_id = _parse_gtf_ids(fields[8], line_number)
            chrom = fields[0]
            strand = fields[6]
            if not chrom:
                raise AnnotationError(f"GTF line {line_number} has an empty chromosome")
            if strand not in {"+", "-", "."}:
                raise AnnotationError(f"GTF line {line_number} has invalid strand: {strand!r}")
            try:
                gtf_start = int(fields[3])
                gtf_end = int(fields[4])
            except ValueError as exc:
                raise AnnotationError(
                    f"GTF line {line_number} has non-integer exon coordinates"
                ) from exc
            if gtf_start < 1 or gtf_end < gtf_start:
                raise AnnotationError(
                    f"GTF line {line_number} has invalid start/end: {gtf_start}-{gtf_end}"
                )

            record = records.get(transcript_id)
            if record is None:
                record = _GtfTranscriptRecord(
                    transcript_id=transcript_id,
                    gene_id=gene_id,
                    chrom=chrom,
                    strand=strand,
                )
                records[transcript_id] = record
            else:
                _validate_gtf_record_consistency(record, gene_id, chrom, strand, line_number)

            record.exons.append((gtf_start - 1, gtf_end))

    if not records:
        raise AnnotationError(f"Annotation contains no GTF exon records: {annotation_path}")

    transcripts: dict[str, Transcript] = {}
    for transcript_id in sorted(records):
        record = records[transcript_id]
        try:
            transcripts[transcript_id] = Transcript(
                transcript_id=record.transcript_id,
                gene_id=record.gene_id,
                chrom=record.chrom,
                strand=record.strand,
                exons=record.exons,
            )
        except ValueError as exc:
            raise AnnotationError(
                f"GTF transcript {transcript_id!r} has invalid exon structure: {exc}"
            ) from exc
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


def _parse_gtf_ids(attributes_text: str, line_number: int) -> tuple[str, str | None]:
    attributes = _parse_gtf_attributes(attributes_text)
    transcript_id = attributes.get("transcript_id")
    if not transcript_id:
        raise AnnotationError(f"GTF line {line_number} exon is missing transcript_id")
    return transcript_id, attributes.get("gene_id")


def _parse_gtf_attributes(attributes_text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in attributes_text.strip().rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item and (" " not in item or item.index("=") < item.index(" ")):
            key, raw_value = item.split("=", 1)
        else:
            parts = item.split(None, 1)
            if len(parts) != 2:
                continue
            key, raw_value = parts
        key = key.strip()
        value = raw_value.strip().strip('"')
        if key:
            attributes[key] = value
    return attributes


def _validate_gtf_record_consistency(
    record: _GtfTranscriptRecord,
    gene_id: str | None,
    chrom: str,
    strand: str,
    line_number: int,
) -> None:
    if chrom != record.chrom:
        raise AnnotationError(
            f"GTF line {line_number} transcript {record.transcript_id!r} "
            f"appears on multiple chromosomes: {record.chrom!r} and {chrom!r}"
        )
    if strand != record.strand:
        raise AnnotationError(
            f"GTF line {line_number} transcript {record.transcript_id!r} "
            f"has inconsistent strand: {record.strand!r} and {strand!r}"
        )
    if gene_id and record.gene_id and gene_id != record.gene_id:
        raise AnnotationError(
            f"GTF line {line_number} transcript {record.transcript_id!r} "
            f"has inconsistent gene_id: {record.gene_id!r} and {gene_id!r}"
        )
    if record.gene_id is None and gene_id is not None:
        record.gene_id = gene_id


def _ensure_annotation_path(path: str | Path) -> Path:
    annotation_path = Path(path)
    if not annotation_path.exists():
        raise FileNotFoundError(f"Annotation file does not exist: {annotation_path}")
    if not annotation_path.is_file():
        raise AnnotationError(f"Annotation path is not a regular file: {annotation_path}")
    return annotation_path


def _open_text(path: Path) -> TextIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def _is_int(value: str) -> bool:
    try:
        int(value)
    except ValueError:
        return False
    return True
