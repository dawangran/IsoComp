#!/usr/bin/env python3
"""Generate synthetic BED12/BAM data for testing IsoComp and RSeQC.

The generated transcriptome contains shared exons, alternative junctions,
terminally similar isoforms, negative-strand transcripts, ambiguous single-exon
transcripts, low-confidence reads, and unassigned reads.

Outputs are intended for manual tool testing:

  isocomp --bam synthetic.sorted.bam --annotation synthetic.complete.bed12 ...
  geneBody_coverage.py -i synthetic.sorted.bam -r synthetic.complete.bed12 ...
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path

CMATCH = 0
CREF_SKIP = 3


@dataclass(frozen=True)
class Transcript:
    transcript_id: str
    chrom: str
    strand: str
    exons: tuple[tuple[int, int], ...]
    expressed: bool = False
    transcript_length: int = field(init=False)
    tx_exons: tuple[tuple[int, int, int, int], ...] = field(init=False)
    span: tuple[int, int] = field(init=False)

    def __post_init__(self) -> None:
        exons = tuple(sorted(self.exons))
        tx_order = exons if self.strand != "-" else tuple(reversed(exons))
        tx_pos = 0
        tx_exons = []
        for genomic_start, genomic_end in tx_order:
            exon_len = genomic_end - genomic_start
            tx_exons.append((genomic_start, genomic_end, tx_pos, tx_pos + exon_len))
            tx_pos += exon_len
        object.__setattr__(self, "exons", exons)
        object.__setattr__(self, "transcript_length", tx_pos)
        object.__setattr__(self, "tx_exons", tuple(tx_exons))
        object.__setattr__(self, "span", (exons[0][0], exons[-1][1]))


@dataclass(frozen=True)
class SyntheticRead:
    read_id: str
    chrom: str
    start: int
    cigar: tuple[tuple[int, int], ...]
    is_reverse: bool
    mapq: int
    scenario: str
    truth_transcript_id: str
    truth_status: str
    truth_is_full_length_like: bool
    truth_tx_start: int | None
    truth_tx_end: int | None
    truth_coverage_fraction: float | None
    truth_dist_to_5p: int | None
    truth_dist_to_3p: int | None


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    transcripts = build_transcripts()
    reads = build_reads(transcripts, replicates=args.replicates, mapq=args.mapq)

    complete_bed = args.out_dir / "synthetic.complete.bed12"
    expressed_bed = args.out_dir / "synthetic.expressed_only.bed12"
    truth_tsv = args.out_dir / "synthetic.truth.tsv"
    unsorted_bam = args.out_dir / "synthetic.unsorted.bam"
    sorted_bam = args.out_dir / "synthetic.sorted.bam"

    write_bed12(complete_bed, transcripts)
    write_bed12(
        expressed_bed,
        {name: tx for name, tx in transcripts.items() if tx.expressed},
    )
    write_truth_tsv(truth_tsv, reads)
    write_bam_and_index(
        unsorted_bam=unsorted_bam,
        sorted_bam=sorted_bam,
        reads=reads,
        reference_length=args.reference_length,
    )
    write_run_commands(
        args.out_dir / "run_commands.sh",
        sorted_bam=sorted_bam,
        complete_bed=complete_bed,
        expressed_bed=expressed_bed,
    )

    print(f"Wrote synthetic data to {args.out_dir}")
    print(f"Complete annotation: {complete_bed}")
    print(f"Expressed-only annotation: {expressed_bed}")
    print(f"Sorted BAM: {sorted_bam}")
    print(f"Truth labels: {truth_tsv}")
    print(f"Commands: {args.out_dir / 'run_commands.sh'}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate synthetic BED12/BAM files for IsoComp and RSeQC testing."
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("benchmark_runs/manual_synthetic"),
        help="Output directory.",
    )
    parser.add_argument(
        "--replicates",
        type=int,
        default=10,
        help="Number of read replicates per scenario.",
    )
    parser.add_argument(
        "--mapq",
        type=int,
        default=60,
        help="Mapping quality assigned to synthetic reads.",
    )
    parser.add_argument(
        "--reference-length",
        type=int,
        default=12000,
        help="Length of synthetic chromosome chrSyn.",
    )
    return parser.parse_args()


def build_transcripts() -> dict[str, Transcript]:
    items = [
        Transcript(
            "TX_POS_MAIN",
            "chrSyn",
            "+",
            ((1000, 1100), (1300, 1400), (1600, 1700)),
            expressed=True,
        ),
        Transcript(
            "TX_POS_ALT_JUNC",
            "chrSyn",
            "+",
            ((1000, 1100), (1350, 1450), (1600, 1700)),
            expressed=True,
        ),
        Transcript(
            "TX_POS_SKIP_DECOY",
            "chrSyn",
            "+",
            ((1000, 1100), (1600, 1700)),
            expressed=False,
        ),
        Transcript(
            "TX_TERM_LONG_DECOY",
            "chrSyn",
            "+",
            ((2500, 2600), (2800, 2900)),
            expressed=False,
        ),
        Transcript(
            "TX_TERM_SHORT",
            "chrSyn",
            "+",
            ((2540, 2600), (2800, 2900)),
            expressed=True,
        ),
        Transcript(
            "TX_NEG_MAIN",
            "chrSyn",
            "-",
            ((4000, 4100), (4300, 4400), (4600, 4700)),
            expressed=True,
        ),
        Transcript(
            "TX_AMB_A",
            "chrSyn",
            "+",
            ((6000, 6200),),
            expressed=True,
        ),
        Transcript(
            "TX_AMB_B_DECOY",
            "chrSyn",
            "+",
            ((6000, 6200),),
            expressed=False,
        ),
        Transcript(
            "TX_LONG_MAIN",
            "chrSyn",
            "+",
            ((7000, 7300), (7500, 7800), (8100, 8500)),
            expressed=True,
        ),
        Transcript(
            "TX_FOUR_EXON",
            "chrSyn",
            "+",
            ((9000, 9100), (9200, 9300), (9400, 9500), (9600, 9700)),
            expressed=False,
        ),
    ]
    return {item.transcript_id: item for item in items}


def build_reads(
    transcripts: dict[str, Transcript],
    *,
    replicates: int,
    mapq: int,
) -> list[SyntheticRead]:
    if replicates < 1:
        raise ValueError("--replicates must be >= 1")

    reads: list[SyntheticRead] = []
    read_cases = [
        ("pos_full", "TX_POS_MAIN", (0, 300), "unique", True),
        ("pos_5p_trunc", "TX_POS_MAIN", (120, 300), "unique", False),
        ("pos_3p_trunc", "TX_POS_MAIN", (0, 180), "unique", False),
        ("pos_internal_partial", "TX_POS_MAIN", (60, 240), "unique", False),
        ("pos_low_cov_terminal_complete", "TX_POS_MAIN", (40, 265), "unique", False),
        ("pos_alt_junction_full", "TX_POS_ALT_JUNC", (0, 300), "unique", True),
        ("neg_full", "TX_NEG_MAIN", (0, 300), "unique", True),
        ("neg_5p_trunc", "TX_NEG_MAIN", (120, 300), "unique", False),
        ("long_5p_near_threshold", "TX_LONG_MAIN", (120, 1000), "unique", False),
        ("terminal_ambiguous", "TX_TERM_SHORT", (0, 160), "ambiguous", False),
        ("single_exon_ambiguous", "TX_AMB_A", (0, 200), "ambiguous", False),
    ]
    for scenario, tx_id, tx_interval, status, is_full_length_like in read_cases:
        transcript = transcripts[tx_id]
        for index in range(replicates):
            reads.append(
                read_from_transcript_interval(
                    read_id=f"{scenario}_{index:03d}",
                    transcript=transcript,
                    tx_interval=tx_interval,
                    truth_status=status,
                    truth_is_full_length_like=is_full_length_like,
                    scenario=scenario,
                    mapq=mapq,
                )
            )

    for index in range(replicates):
        reads.append(
            intronic_low_confidence_read(
                read_id=f"intronic_low_conf_{index:03d}",
                transcript=transcripts["TX_POS_MAIN"],
                mapq=mapq,
            )
        )
        reads.append(
            junction_conflict_low_confidence_read(
                read_id=f"junction_conflict_low_conf_{index:03d}",
                transcript=transcripts["TX_FOUR_EXON"],
                mapq=mapq,
            )
        )
        reads.append(unassigned_read(read_id=f"unassigned_{index:03d}", mapq=mapq))

    return reads


def read_from_transcript_interval(
    *,
    read_id: str,
    transcript: Transcript,
    tx_interval: tuple[int, int],
    truth_status: str,
    truth_is_full_length_like: bool,
    scenario: str,
    mapq: int,
) -> SyntheticRead:
    blocks = tx_interval_to_genomic_blocks(transcript, tx_interval)
    start, cigar = cigar_from_blocks(blocks)
    tx_start, tx_end = tx_interval
    coverage_fraction = (tx_end - tx_start) / transcript.transcript_length
    return SyntheticRead(
        read_id=read_id,
        chrom=transcript.chrom,
        start=start,
        cigar=tuple(cigar),
        is_reverse=transcript.strand == "-",
        mapq=mapq,
        scenario=scenario,
        truth_transcript_id=transcript.transcript_id,
        truth_status=truth_status,
        truth_is_full_length_like=truth_is_full_length_like,
        truth_tx_start=tx_start,
        truth_tx_end=tx_end,
        truth_coverage_fraction=coverage_fraction,
        truth_dist_to_5p=tx_start,
        truth_dist_to_3p=transcript.transcript_length - tx_end,
    )


def intronic_low_confidence_read(
    *,
    read_id: str,
    transcript: Transcript,
    mapq: int,
) -> SyntheticRead:
    start, end = transcript.span
    return SyntheticRead(
        read_id=read_id,
        chrom=transcript.chrom,
        start=start,
        cigar=((CMATCH, end - start),),
        is_reverse=False,
        mapq=mapq,
        scenario="intronic_low_conf",
        truth_transcript_id=transcript.transcript_id,
        truth_status="low_confidence",
        truth_is_full_length_like=False,
        truth_tx_start=None,
        truth_tx_end=None,
        truth_coverage_fraction=None,
        truth_dist_to_5p=None,
        truth_dist_to_3p=None,
    )


def junction_conflict_low_confidence_read(
    *,
    read_id: str,
    transcript: Transcript,
    mapq: int,
) -> SyntheticRead:
    blocks = [(9000, 9100), (9200, 9300), (9410, 9510), (9610, 9710)]
    start, cigar = cigar_from_blocks(blocks)
    return SyntheticRead(
        read_id=read_id,
        chrom=transcript.chrom,
        start=start,
        cigar=tuple(cigar),
        is_reverse=False,
        mapq=mapq,
        scenario="junction_conflict_low_conf",
        truth_transcript_id=transcript.transcript_id,
        truth_status="low_confidence",
        truth_is_full_length_like=False,
        truth_tx_start=None,
        truth_tx_end=None,
        truth_coverage_fraction=None,
        truth_dist_to_5p=None,
        truth_dist_to_3p=None,
    )


def unassigned_read(*, read_id: str, mapq: int) -> SyntheticRead:
    return SyntheticRead(
        read_id=read_id,
        chrom="chrSyn",
        start=11000,
        cigar=((CMATCH, 150),),
        is_reverse=False,
        mapq=mapq,
        scenario="unassigned",
        truth_transcript_id="",
        truth_status="unassigned",
        truth_is_full_length_like=False,
        truth_tx_start=None,
        truth_tx_end=None,
        truth_coverage_fraction=None,
        truth_dist_to_5p=None,
        truth_dist_to_3p=None,
    )


def tx_interval_to_genomic_blocks(
    transcript: Transcript,
    tx_interval: tuple[int, int],
) -> list[tuple[int, int]]:
    tx_start, tx_end = tx_interval
    if tx_start < 0 or tx_end > transcript.transcript_length or tx_end <= tx_start:
        raise ValueError(
            f"Invalid transcript interval for {transcript.transcript_id}: {tx_interval}"
        )

    blocks: list[tuple[int, int]] = []
    for genomic_start, genomic_end, exon_tx_start, exon_tx_end in transcript.tx_exons:
        overlap_start = max(tx_start, exon_tx_start)
        overlap_end = min(tx_end, exon_tx_end)
        if overlap_end <= overlap_start:
            continue
        if transcript.strand == "-":
            block_start = genomic_end - (overlap_end - exon_tx_start)
            block_end = genomic_end - (overlap_start - exon_tx_start)
        else:
            block_start = genomic_start + (overlap_start - exon_tx_start)
            block_end = genomic_start + (overlap_end - exon_tx_start)
        blocks.append((block_start, block_end))

    return sorted(blocks)


def cigar_from_blocks(blocks: list[tuple[int, int]]) -> tuple[int, list[tuple[int, int]]]:
    if not blocks:
        raise ValueError("Cannot make CIGAR from empty blocks")
    blocks = sorted(blocks)
    start = blocks[0][0]
    cigar: list[tuple[int, int]] = []
    previous_end: int | None = None
    for block_start, block_end in blocks:
        if previous_end is not None:
            gap = block_start - previous_end
            if gap > 0:
                cigar.append((CREF_SKIP, gap))
        cigar.append((CMATCH, block_end - block_start))
        previous_end = block_end
    return start, cigar


def query_length(cigar: tuple[tuple[int, int], ...]) -> int:
    return sum(length for op, length in cigar if op == CMATCH)


def cigar_to_string(cigar: tuple[tuple[int, int], ...]) -> str:
    op_names = {CMATCH: "M", CREF_SKIP: "N"}
    return "".join(f"{length}{op_names[op]}" for op, length in cigar)


def write_bed12(path: Path, transcripts: dict[str, Transcript]) -> None:
    rows = []
    for transcript in transcripts.values():
        chrom_start, chrom_end = transcript.span
        block_sizes = [end - start for start, end in transcript.exons]
        block_starts = [start - chrom_start for start, _ in transcript.exons]
        rows.append(
            "\t".join(
                [
                    transcript.chrom,
                    str(chrom_start),
                    str(chrom_end),
                    transcript.transcript_id,
                    "0",
                    transcript.strand,
                    str(chrom_start),
                    str(chrom_end),
                    "0",
                    str(len(transcript.exons)),
                    ",".join(str(value) for value in block_sizes),
                    ",".join(str(value) for value in block_starts),
                ]
            )
        )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def write_truth_tsv(path: Path, reads: list[SyntheticRead]) -> None:
    columns = [
        "read_id",
        "scenario",
        "truth_transcript_id",
        "truth_status",
        "truth_is_full_length_like",
        "truth_tx_start",
        "truth_tx_end",
        "truth_coverage_fraction",
        "truth_dist_to_5p",
        "truth_dist_to_3p",
        "cigar",
        "mapq",
    ]
    lines = ["\t".join(columns)]
    for read in reads:
        values = [
            read.read_id,
            read.scenario,
            read.truth_transcript_id,
            read.truth_status,
            int(read.truth_is_full_length_like),
            nullable(read.truth_tx_start),
            nullable(read.truth_tx_end),
            nullable(read.truth_coverage_fraction),
            nullable(read.truth_dist_to_5p),
            nullable(read.truth_dist_to_3p),
            cigar_to_string(read.cigar),
            read.mapq,
        ]
        lines.append("\t".join(str(value) for value in values))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def nullable(value: object | None) -> str:
    return "" if value is None else str(value)


def write_bam_and_index(
    *,
    unsorted_bam: Path,
    sorted_bam: Path,
    reads: list[SyntheticRead],
    reference_length: int,
) -> None:
    try:
        import pysam
    except ImportError as exc:
        raise SystemExit(
            "pysam is required to write BAM files. Install it with `pip install pysam`, "
            "or run this script inside an RSeQC environment."
        ) from exc

    header = {
        "HD": {"VN": "1.0"},
        "SQ": [{"SN": "chrSyn", "LN": reference_length}],
    }
    with pysam.AlignmentFile(unsorted_bam, "wb", header=header) as bam:
        for read in reads:
            segment = pysam.AlignedSegment()
            segment.query_name = read.read_id
            segment.query_sequence = "A" * query_length(read.cigar)
            segment.flag = 16 if read.is_reverse else 0
            segment.reference_id = bam.get_tid(read.chrom)
            segment.reference_start = read.start
            segment.mapping_quality = read.mapq
            segment.cigar = list(read.cigar)
            segment.query_qualities = pysam.qualitystring_to_array(
                "I" * len(segment.query_sequence)
            )
            bam.write(segment)

    pysam.sort("-o", str(sorted_bam), str(unsorted_bam))
    pysam.index(str(sorted_bam))


def write_run_commands(
    path: Path,
    *,
    sorted_bam: Path,
    complete_bed: Path,
    expressed_bed: Path,
) -> None:
    repo_root = Path.cwd().resolve()
    output_dir = path.parent
    sorted_bam = sorted_bam.resolve()
    complete_bed = complete_bed.resolve()
    expressed_bed = expressed_bed.resolve()
    output_dir = output_dir.resolve()
    text = f"""#!/usr/bin/env bash
set -euo pipefail

# Override these if the commands are not on PATH.
# Example:
#   ISOCOMP_CMD="benchmark_runs/rseqc_venv/bin/python -m isocomp.cli" \\
#   RSEQC_CMD="benchmark_runs/rseqc_venv/bin/geneBody_coverage.py" \\
#   bash {path.name}
ISOCOMP_CMD="${{ISOCOMP_CMD:-isocomp}}"
RSEQC_CMD="${{RSEQC_CMD:-geneBody_coverage.py}}"
export PYTHONPATH="{repo_root}:${{PYTHONPATH:-}}"

# IsoComp on the complete annotation.
$ISOCOMP_CMD \\
  --bam {sorted_bam} \\
  --annotation {complete_bed} \\
  --out {output_dir / 'isocomp.complete'} \\
  --bin-num 100 \\
  --min-mapq 20 \\
  --tss-tol 100 \\
  --tes-tol 100 \\
  --force

# RSeQC gene-body coverage on the complete annotation.
$RSEQC_CMD \\
  -i {sorted_bam} \\
  -r {complete_bed} \\
  -o {output_dir / 'rseqc.complete'} \\
  -f png

# RSeQC on an expressed-only annotation, useful as a matched-annotation control.
$RSEQC_CMD \\
  -i {sorted_bam} \\
  -r {expressed_bed} \\
  -o {output_dir / 'rseqc.expressed_only'} \\
  -f png
"""
    path.write_text(text, encoding="utf-8")
    path.chmod(0o755)


if __name__ == "__main__":
    main()
