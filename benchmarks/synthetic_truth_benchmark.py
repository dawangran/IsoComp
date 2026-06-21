#!/usr/bin/env python3
"""Synthetic truth benchmark and parameter sensitivity analysis for IsoComp.

The benchmark builds a small transcriptome with shared exons, alternative
junctions, terminal ambiguity, negative-strand transcripts, low-confidence
intronic reads, and unassigned reads. It then runs IsoComp's core assignment
logic against known per-read truth labels and sweeps key thresholds.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from isocomp.annotation import read_bed12
from isocomp.assigner import assign_read
from isocomp.candidate import CandidateIndex
from isocomp.metrics import build_read_metrics
from isocomp.models import AssignmentStatus, ReadAlignment, ReadMetrics, Transcript
from isocomp.projection import project_blocks_to_transcript

CMATCH = 0
CREF_SKIP = 3


@dataclass(frozen=True)
class TranscriptSpec:
    transcript_id: str
    chrom: str
    strand: str
    exons: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class ReadSpec:
    read_id: str
    chrom: str
    start: int
    cigar: tuple[tuple[int, int], ...]
    is_reverse: bool
    mapq: int
    truth_transcript_id: str
    truth_status: AssignmentStatus
    truth_is_full_length_like: bool
    truth_is_5p_complete: bool | None
    truth_is_3p_complete: bool | None
    truth_coverage_fraction: float | None
    truth_dist_to_5p: int | None
    truth_dist_to_3p: int | None
    scenario: str


@dataclass(frozen=True)
class ParameterSet:
    parameter_set: str
    min_overlap: int
    junction_tol: int
    unique_threshold: float
    margin_threshold: float
    min_unspliced_coverage_for_unique: float
    tss_tol: int
    tes_tol: int
    coverage_threshold: float


def main() -> None:
    args = parse_args()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    configure_plot_cache(out_dir)

    transcripts = build_transcripts()
    read_specs = build_read_specs(transcripts, replicates=args.replicates)
    bed_path = out_dir / "synthetic_truth.bed12"
    truth_path = out_dir / "synthetic_truth.tsv"

    write_bed12(bed_path, transcripts)
    if args.write_bam:
        write_bam(out_dir / "synthetic_truth.bam", read_specs)
    write_truth_table(truth_path, read_specs)

    parsed_transcripts = read_bed12(bed_path)
    reads = load_reads_from_specs(read_specs, min_mapq=args.min_mapq)
    truth = pd.read_csv(truth_path, sep="\t")
    parameter_sets = build_parameter_grid(
        args.min_overlaps,
        args.junction_tols,
        args.unique_thresholds,
        args.margin_thresholds,
        args.min_unspliced_coverages,
        args.terminal_tolerances,
        args.coverage_thresholds,
    )

    per_read_frames: list[pd.DataFrame] = []
    summary_rows: list[dict[str, Any]] = []
    for parameters in parameter_sets:
        per_read = run_parameter_set(reads, parsed_transcripts, truth, parameters)
        metrics = summarize_predictions(per_read, parameters)
        per_read_frames.append(per_read)
        summary_rows.append(metrics)

    sensitivity_per_read = pd.concat(per_read_frames, ignore_index=True)
    sensitivity_summary = pd.DataFrame(summary_rows)
    sensitivity_per_read.to_csv(out_dir / "sensitivity_per_read.tsv", sep="\t", index=False)
    sensitivity_summary.to_csv(out_dir / "sensitivity_summary.tsv", sep="\t", index=False)

    default_per_read = sensitivity_per_read.loc[
        sensitivity_per_read["parameter_set"] == "default"
    ].copy()
    default_metrics = sensitivity_summary.loc[
        sensitivity_summary["parameter_set"] == "default"
    ].copy()
    default_per_read.to_csv(out_dir / "default.per_read.tsv", sep="\t", index=False)
    default_metrics.to_csv(out_dir / "default.metrics.tsv", sep="\t", index=False)

    coverage_profiles, coverage_metrics = compute_coverage_profiles(
        reads,
        parsed_transcripts,
        truth,
        default_per_read,
        bin_num=args.coverage_bin_num,
    )
    coverage_profiles.to_csv(out_dir / "rseqc_style_coverage_profiles.tsv", sep="\t", index=False)
    coverage_metrics.to_csv(out_dir / "rseqc_style_coverage_metrics.tsv", sep="\t", index=False)

    if not args.no_plots:
        plot_status_confusion(default_per_read, out_dir / "synthetic_status_confusion.png")
        plot_sensitivity_f1(sensitivity_summary, out_dir / "sensitivity_full_length_f1.png")
        plot_coverage_profiles(
            coverage_profiles,
            out_dir / "rseqc_style_coverage_profiles.png",
            normalization=args.profile_normalization,
        )
        write_article_panel_outputs(
            default_metrics=default_metrics,
            sensitivity_summary=sensitivity_summary,
            coverage_profiles=coverage_profiles,
            coverage_metrics=coverage_metrics,
            out_dir=out_dir,
            normalization=args.profile_normalization,
        )
    else:
        write_article_panel_outputs(
            default_metrics=default_metrics,
            sensitivity_summary=sensitivity_summary,
            coverage_profiles=coverage_profiles,
            coverage_metrics=coverage_metrics,
            out_dir=out_dir,
            normalization=args.profile_normalization,
            write_plots=False,
        )

    print(f"Wrote synthetic benchmark to {out_dir}")
    if not default_metrics.empty:
        row = default_metrics.iloc[0]
        print(
            "Default metrics: "
            f"status_accuracy={row['status_accuracy']:.3f}, "
            f"unique_transcript_accuracy={row['unique_transcript_accuracy']:.3f}, "
            f"full_length_f1={row['full_length_f1']:.3f}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate an IsoComp synthetic truth benchmark and threshold sensitivity grid."
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("benchmark_runs/synthetic_truth"),
        help="Output directory for synthetic data and benchmark tables.",
    )
    parser.add_argument(
        "--replicates",
        type=int,
        default=20,
        help="Number of read replicates per scenario.",
    )
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=20,
        help="Minimum MAPQ used when loading synthetic reads.",
    )
    parser.add_argument(
        "--write-bam",
        action="store_true",
        help="Also write synthetic_truth.bam. Requires pysam in the active environment.",
    )
    parser.add_argument(
        "--min-overlaps",
        type=int_list,
        default=[20, 50, 100],
        help="Comma-separated candidate minimum exonic overlaps to sweep.",
    )
    parser.add_argument(
        "--junction-tols",
        type=int_list,
        default=[0, 5, 10],
        help="Comma-separated splice-junction tolerances to sweep.",
    )
    parser.add_argument(
        "--unique-thresholds",
        type=float_list,
        default=[0.7, 0.8, 0.9],
        help="Comma-separated top-score thresholds to sweep.",
    )
    parser.add_argument(
        "--margin-thresholds",
        type=float_list,
        default=[0.05, 0.1, 0.2],
        help="Comma-separated top-minus-second score margins to sweep.",
    )
    parser.add_argument(
        "--min-unspliced-coverages",
        type=float_list,
        default=[0.0, 0.2, 0.4],
        help=(
            "Comma-separated minimum transcript coverage values for a "
            "junctionless read to be unique on a spliced transcript."
        ),
    )
    parser.add_argument(
        "--terminal-tolerances",
        type=int_list,
        default=[50, 100, 200],
        help="Comma-separated 5'/3' terminal tolerances to sweep.",
    )
    parser.add_argument(
        "--coverage-thresholds",
        type=float_list,
        default=[0.7, 0.8, 0.9],
        help="Comma-separated coverage-fraction thresholds for full-length-like calls.",
    )
    parser.add_argument(
        "--coverage-bin-num",
        type=int,
        default=100,
        help="Number of transcript-body bins for IsoComp and RSeQC-style coverage profiles.",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip PNG plot generation.",
    )
    parser.add_argument(
        "--profile-normalization",
        choices=["max", "mean"],
        default="max",
        help="Normalization used for manuscript coverage profile plots.",
    )
    return parser.parse_args()


def int_list(value: str) -> list[int]:
    return [int(item) for item in value.split(",") if item]


def float_list(value: str) -> list[float]:
    return [float(item) for item in value.split(",") if item]


def configure_plot_cache(out_dir: Path) -> None:
    cache_dir = out_dir / ".plot-cache"
    matplotlib_cache = cache_dir / "matplotlib"
    xdg_cache = cache_dir / "xdg"
    matplotlib_cache.mkdir(parents=True, exist_ok=True)
    xdg_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_cache))
    os.environ.setdefault("XDG_CACHE_HOME", str(xdg_cache))


def build_transcripts() -> dict[str, Transcript]:
    specs = [
        TranscriptSpec(
            "TX_POS_MAIN",
            "chrSyn",
            "+",
            ((1000, 1100), (1300, 1400), (1600, 1700)),
        ),
        TranscriptSpec(
            "TX_POS_ALT_JUNC",
            "chrSyn",
            "+",
            ((1000, 1100), (1350, 1450), (1600, 1700)),
        ),
        TranscriptSpec(
            "TX_POS_SKIP",
            "chrSyn",
            "+",
            ((1000, 1100), (1600, 1700)),
        ),
        TranscriptSpec(
            "TX_TERM_LONG",
            "chrSyn",
            "+",
            ((2500, 2600), (2800, 2900)),
        ),
        TranscriptSpec(
            "TX_TERM_SHORT",
            "chrSyn",
            "+",
            ((2540, 2600), (2800, 2900)),
        ),
        TranscriptSpec(
            "TX_NEG_MAIN",
            "chrSyn",
            "-",
            ((4000, 4100), (4300, 4400), (4600, 4700)),
        ),
        TranscriptSpec(
            "TX_LONG_MAIN",
            "chrSyn",
            "+",
            ((7000, 7300), (7500, 7800), (8100, 8500)),
        ),
        TranscriptSpec(
            "TX_FOUR_EXON",
            "chrSyn",
            "+",
            ((9000, 9100), (9200, 9300), (9400, 9500), (9600, 9700)),
        ),
        TranscriptSpec(
            "TX_AMB_A",
            "chrSyn",
            "+",
            ((6000, 6200),),
        ),
        TranscriptSpec(
            "TX_AMB_B",
            "chrSyn",
            "+",
            ((6000, 6200),),
        ),
    ]
    return {
        spec.transcript_id: Transcript(
            transcript_id=spec.transcript_id,
            gene_id=None,
            chrom=spec.chrom,
            strand=spec.strand,
            exons=list(spec.exons),
        )
        for spec in specs
    }


def build_redundant_transcripts(
    transcripts: dict[str, Transcript],
) -> dict[str, Transcript]:
    specs = [
        TranscriptSpec(
            "TX_POS_MAIN_REDUNDANT_DUP",
            "chrSyn",
            "+",
            ((1000, 1100), (1300, 1400), (1600, 1700)),
        ),
        TranscriptSpec(
            "TX_POS_MAIN_REDUNDANT_5P_SHORT",
            "chrSyn",
            "+",
            ((1040, 1100), (1300, 1400), (1600, 1700)),
        ),
        TranscriptSpec(
            "TX_POS_SHARED_EXON_DECOY",
            "chrSyn",
            "+",
            ((1000, 1100), (1300, 1400)),
        ),
        TranscriptSpec(
            "TX_NEG_MAIN_REDUNDANT_DUP",
            "chrSyn",
            "-",
            ((4000, 4100), (4300, 4400), (4600, 4700)),
        ),
        TranscriptSpec(
            "TX_LONG_INTERNAL_DECOY",
            "chrSyn",
            "+",
            ((7500, 7800), (8100, 8500)),
        ),
    ]
    redundant = dict(transcripts)
    for spec in specs:
        redundant[spec.transcript_id] = Transcript(
            transcript_id=spec.transcript_id,
            gene_id=None,
            chrom=spec.chrom,
            strand=spec.strand,
            exons=list(spec.exons),
        )
    return redundant


def build_read_specs(transcripts: dict[str, Transcript], *, replicates: int) -> list[ReadSpec]:
    if replicates < 1:
        raise ValueError("--replicates must be >= 1")

    reads: list[ReadSpec] = []
    scenarios = [
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

    for scenario, transcript_id, tx_interval, status, full_length in scenarios:
        transcript = transcripts[transcript_id]
        for index in range(replicates):
            reads.append(
                read_from_transcript_interval(
                    read_id=f"{scenario}_{index:03d}",
                    transcript=transcript,
                    tx_interval=tx_interval,
                    truth_status=status,
                    truth_is_full_length_like=full_length,
                    scenario=scenario,
                )
            )

    for index in range(replicates):
        reads.append(
            intronic_low_confidence_read(
                read_id=f"intronic_low_conf_{index:03d}",
                transcript=transcripts["TX_POS_MAIN"],
            )
        )
        reads.append(
            junction_conflict_low_confidence_read(
                read_id=f"junction_conflict_low_conf_{index:03d}",
                transcript=transcripts["TX_FOUR_EXON"],
            )
        )
        reads.append(
            unassigned_read(
                read_id=f"unassigned_{index:03d}",
            )
        )
    return reads


def read_from_transcript_interval(
    *,
    read_id: str,
    transcript: Transcript,
    tx_interval: tuple[int, int],
    truth_status: str,
    truth_is_full_length_like: bool,
    scenario: str,
) -> ReadSpec:
    blocks = tx_interval_to_genomic_blocks(transcript, tx_interval)
    start, cigar = cigar_from_blocks(blocks)
    coverage_fraction = (tx_interval[1] - tx_interval[0]) / transcript.transcript_length
    dist_to_5p = tx_interval[0]
    dist_to_3p = transcript.transcript_length - tx_interval[1]
    return ReadSpec(
        read_id=read_id,
        chrom=transcript.chrom,
        start=start,
        cigar=tuple(cigar),
        is_reverse=transcript.strand == "-",
        mapq=60,
        truth_transcript_id=transcript.transcript_id,
        truth_status=truth_status,  # type: ignore[arg-type]
        truth_is_full_length_like=truth_is_full_length_like,
        truth_is_5p_complete=dist_to_5p <= 100,
        truth_is_3p_complete=dist_to_3p <= 100,
        truth_coverage_fraction=coverage_fraction,
        truth_dist_to_5p=dist_to_5p,
        truth_dist_to_3p=dist_to_3p,
        scenario=scenario,
    )


def intronic_low_confidence_read(*, read_id: str, transcript: Transcript) -> ReadSpec:
    block = (transcript.span[0], transcript.span[1])
    return ReadSpec(
        read_id=read_id,
        chrom=transcript.chrom,
        start=block[0],
        cigar=((CMATCH, block[1] - block[0]),),
        is_reverse=False,
        mapq=60,
        truth_transcript_id=transcript.transcript_id,
        truth_status="low_confidence",
        truth_is_full_length_like=False,
        truth_is_5p_complete=None,
        truth_is_3p_complete=None,
        truth_coverage_fraction=None,
        truth_dist_to_5p=None,
        truth_dist_to_3p=None,
        scenario="intronic_low_conf",
    )


def junction_conflict_low_confidence_read(*, read_id: str, transcript: Transcript) -> ReadSpec:
    blocks = [(9000, 9100), (9200, 9300), (9410, 9510), (9610, 9710)]
    start, cigar = cigar_from_blocks(blocks)
    return ReadSpec(
        read_id=read_id,
        chrom=transcript.chrom,
        start=start,
        cigar=tuple(cigar),
        is_reverse=False,
        mapq=60,
        truth_transcript_id=transcript.transcript_id,
        truth_status="low_confidence",
        truth_is_full_length_like=False,
        truth_is_5p_complete=None,
        truth_is_3p_complete=None,
        truth_coverage_fraction=None,
        truth_dist_to_5p=None,
        truth_dist_to_3p=None,
        scenario="junction_conflict_low_conf",
    )


def unassigned_read(*, read_id: str) -> ReadSpec:
    return ReadSpec(
        read_id=read_id,
        chrom="chrSyn",
        start=11000,
        cigar=((CMATCH, 150),),
        is_reverse=False,
        mapq=60,
        truth_transcript_id="",
        truth_status="unassigned",
        truth_is_full_length_like=False,
        truth_is_5p_complete=None,
        truth_is_3p_complete=None,
        truth_coverage_fraction=None,
        truth_dist_to_5p=None,
        truth_dist_to_3p=None,
        scenario="unassigned",
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
    for exon in transcript.tx_exons:
        overlap_start = max(tx_start, exon.tx_start)
        overlap_end = min(tx_end, exon.tx_end)
        if overlap_end <= overlap_start:
            continue
        if transcript.strand == "-":
            genomic_start = exon.genomic_end - (overlap_end - exon.tx_start)
            genomic_end = exon.genomic_end - (overlap_start - exon.tx_start)
        else:
            genomic_start = exon.genomic_start + (overlap_start - exon.tx_start)
            genomic_end = exon.genomic_start + (overlap_end - exon.tx_start)
        blocks.append((genomic_start, genomic_end))

    return sorted(blocks)


def cigar_from_blocks(blocks: list[tuple[int, int]]) -> tuple[int, list[tuple[int, int]]]:
    if not blocks:
        raise ValueError("Cannot create a CIGAR from empty blocks")
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


def write_bed12(path: Path, transcripts: dict[str, Transcript]) -> None:
    rows = []
    for transcript in sorted(transcripts.values(), key=lambda item: item.transcript_id):
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
                    ",".join(str(item) for item in block_sizes),
                    ",".join(str(item) for item in block_starts),
                ]
            )
        )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def write_bam(path: Path, reads: list[ReadSpec]) -> None:
    import pysam

    header = {
        "HD": {"VN": "1.0"},
        "SQ": [{"SN": "chrSyn", "LN": 12000}],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
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


def query_length(cigar: tuple[tuple[int, int], ...]) -> int:
    query_consuming = {0, 1, 4, 7, 8}
    return sum(length for op, length in cigar if op in query_consuming)


def write_truth_table(path: Path, reads: list[ReadSpec]) -> None:
    pd.DataFrame(
        [
            {
                "read_id": read.read_id,
                "truth_transcript_id": read.truth_transcript_id,
                "truth_status": read.truth_status,
                "truth_is_full_length_like": int(read.truth_is_full_length_like),
                "truth_is_5p_complete": nullable_bool(read.truth_is_5p_complete),
                "truth_is_3p_complete": nullable_bool(read.truth_is_3p_complete),
                "truth_coverage_fraction": read.truth_coverage_fraction,
                "truth_dist_to_5p": read.truth_dist_to_5p,
                "truth_dist_to_3p": read.truth_dist_to_3p,
                "scenario": read.scenario,
            }
            for read in reads
        ]
    ).to_csv(path, sep="\t", index=False)


def nullable_bool(value: bool | None) -> int | None:
    if value is None:
        return None
    return int(value)


def load_reads_from_specs(reads: list[ReadSpec], *, min_mapq: int) -> list[ReadAlignment]:
    parsed: list[ReadAlignment] = []
    for read in reads:
        if read.mapq < min_mapq:
            continue
        blocks, junctions = blocks_and_junctions_from_cigar(read.start, read.cigar)
        parsed.append(
            ReadAlignment(
                read_id=read.read_id,
                chrom=read.chrom,
                genomic_start=blocks[0][0],
                genomic_end=blocks[-1][1],
                blocks=blocks,
                junctions=junctions,
                aligned_length=sum(end - start for start, end in blocks),
                mapq=read.mapq,
                cigar=cigar_to_string(read.cigar),
                is_reverse=read.is_reverse,
                softclip_5p=0,
                softclip_3p=0,
            )
        )
    return parsed


def blocks_and_junctions_from_cigar(
    reference_start: int,
    cigar: tuple[tuple[int, int], ...],
) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    ref_pos = reference_start
    current_block_start: int | None = None
    blocks: list[tuple[int, int]] = []
    junctions: list[tuple[int, int]] = []
    for op, length in cigar:
        if op == CMATCH:
            if current_block_start is None:
                current_block_start = ref_pos
            ref_pos += length
        elif op == CREF_SKIP:
            if current_block_start is not None and ref_pos > current_block_start:
                blocks.append((current_block_start, ref_pos))
            junction_start = ref_pos
            ref_pos += length
            junctions.append((junction_start, ref_pos))
            current_block_start = None
        else:
            raise ValueError(f"Unsupported synthetic CIGAR op: {op}")
    if current_block_start is not None and ref_pos > current_block_start:
        blocks.append((current_block_start, ref_pos))
    if not blocks:
        raise ValueError("Synthetic read produced no aligned blocks")
    return blocks, junctions


def cigar_to_string(cigar: tuple[tuple[int, int], ...]) -> str:
    op_names = {CMATCH: "M", CREF_SKIP: "N"}
    return "".join(f"{length}{op_names[op]}" for op, length in cigar)


def build_parameter_grid(
    min_overlaps: list[int],
    junction_tols: list[int],
    unique_thresholds: list[float],
    margin_thresholds: list[float],
    min_unspliced_coverages: list[float],
    terminal_tolerances: list[int],
    coverage_thresholds: list[float],
) -> list[ParameterSet]:
    parameter_sets: list[ParameterSet] = [
        ParameterSet("default", 50, 5, 0.8, 0.1, 0.2, 100, 100, 0.8)
    ]
    seen = {signature(parameter_sets[0])}
    for min_overlap in min_overlaps:
        for junction_tol in junction_tols:
            for unique_threshold in unique_thresholds:
                for margin_threshold in margin_thresholds:
                    for min_unspliced_coverage in min_unspliced_coverages:
                        for terminal_tolerance in terminal_tolerances:
                            for coverage_threshold in coverage_thresholds:
                                item = ParameterSet(
                                    parameter_set=(
                                        f"mo{min_overlap}_jt{junction_tol}_"
                                        f"ut{unique_threshold:g}_mt{margin_threshold:g}_"
                                        f"uc{min_unspliced_coverage:g}_"
                                        f"tt{terminal_tolerance}_ct{coverage_threshold:g}"
                                    ),
                                    min_overlap=min_overlap,
                                    junction_tol=junction_tol,
                                    unique_threshold=unique_threshold,
                                    margin_threshold=margin_threshold,
                                    min_unspliced_coverage_for_unique=min_unspliced_coverage,
                                    tss_tol=terminal_tolerance,
                                    tes_tol=terminal_tolerance,
                                    coverage_threshold=coverage_threshold,
                                )
                                key = signature(item)
                                if key in seen:
                                    continue
                                seen.add(key)
                                parameter_sets.append(item)
    return parameter_sets


def signature(parameters: ParameterSet) -> tuple[int, int, float, float, float, int, int, float]:
    return (
        parameters.min_overlap,
        parameters.junction_tol,
        parameters.unique_threshold,
        parameters.margin_threshold,
        parameters.min_unspliced_coverage_for_unique,
        parameters.tss_tol,
        parameters.tes_tol,
        parameters.coverage_threshold,
    )


def run_parameter_set(
    reads: list[ReadAlignment],
    transcripts: dict[str, Transcript],
    truth: pd.DataFrame,
    parameters: ParameterSet,
) -> pd.DataFrame:
    candidate_index = CandidateIndex(transcripts)
    rows: list[dict[str, Any]] = []
    truth_by_read = truth.set_index("read_id")

    for read in reads:
        candidates = candidate_index.query(
            read,
            min_overlap=parameters.min_overlap,
            strandness="unstranded",
        )
        assignment = assign_read(
            read,
            candidates,
            junction_tol=parameters.junction_tol,
            unique_threshold=parameters.unique_threshold,
            margin_threshold=parameters.margin_threshold,
            min_unspliced_coverage_for_unique=parameters.min_unspliced_coverage_for_unique,
        )
        metric = build_read_metrics(
            read,
            assignment,
            tss_tol=parameters.tss_tol,
            tes_tol=parameters.tes_tol,
            full_length_coverage=parameters.coverage_threshold,
            min_terminal_anchor=10,
        )
        truth_row = truth_by_read.loc[read.read_id]
        rows.append(prediction_row(metric, truth_row, parameters))
    return pd.DataFrame(rows)


def prediction_row(
    metric: ReadMetrics,
    truth_row: pd.Series,
    parameters: ParameterSet,
) -> dict[str, Any]:
    pred_is_5p_complete = metric.is_5p_complete
    pred_is_3p_complete = metric.is_3p_complete
    pred_full_length_like = bool(
        metric.assignment_status == "unique"
        and metric.coverage_fraction is not None
        and metric.coverage_fraction >= parameters.coverage_threshold
        and pred_is_5p_complete
        and pred_is_3p_complete
    )
    second_score = metric.second_best_score
    score_margin = (
        metric.assignment_score - second_score
        if second_score is not None and not math.isnan(second_score)
        else metric.assignment_score
    )
    return {
        "parameter_set": parameters.parameter_set,
        "min_overlap": parameters.min_overlap,
        "junction_tol": parameters.junction_tol,
        "unique_threshold": parameters.unique_threshold,
        "margin_threshold": parameters.margin_threshold,
        "min_unspliced_coverage_for_unique": parameters.min_unspliced_coverage_for_unique,
        "tss_tol": parameters.tss_tol,
        "tes_tol": parameters.tes_tol,
        "coverage_threshold": parameters.coverage_threshold,
        "read_id": metric.read_id,
        "scenario": truth_row["scenario"],
        "truth_transcript_id": truth_row["truth_transcript_id"],
        "pred_transcript_id": metric.transcript_id or "",
        "truth_status": truth_row["truth_status"],
        "pred_status": metric.assignment_status,
        "status_correct": int(metric.assignment_status == truth_row["truth_status"]),
        "truth_is_full_length_like": int(truth_row["truth_is_full_length_like"]),
        "pred_is_full_length_like": int(pred_full_length_like),
        "truth_dist_to_5p": truth_row["truth_dist_to_5p"],
        "pred_dist_to_5p": metric.dist_to_5p,
        "truth_dist_to_3p": truth_row["truth_dist_to_3p"],
        "pred_dist_to_3p": metric.dist_to_3p,
        "truth_coverage_fraction": truth_row["truth_coverage_fraction"],
        "pred_coverage_fraction": metric.coverage_fraction,
        "assignment_score": metric.assignment_score,
        "second_best_score": second_score,
        "score_margin": score_margin,
        "junction_match_count": metric.junction_match_count,
        "junction_precision": metric.junction_precision,
        "junction_recall": metric.junction_recall,
        "exon_overlap_score": metric.exon_overlap_score,
    }


def summarize_predictions(per_read: pd.DataFrame, parameters: ParameterSet) -> dict[str, Any]:
    truth_full = per_read["truth_is_full_length_like"].astype(bool)
    pred_full = per_read["pred_is_full_length_like"].astype(bool)
    tp = int((truth_full & pred_full).sum())
    fp = int((~truth_full & pred_full).sum())
    fn = int((truth_full & ~pred_full).sum())
    precision = safe_divide(tp, tp + fp)
    recall = safe_divide(tp, tp + fn)
    f1 = safe_divide(2 * precision * recall, precision + recall)

    truth_unique = per_read["truth_status"] == "unique"
    unique_transcript_accuracy = safe_divide(
        int(
            (
                truth_unique
                & (per_read["pred_status"] == "unique")
                & (per_read["pred_transcript_id"] == per_read["truth_transcript_id"])
            ).sum()
        ),
        int(truth_unique.sum()),
    )
    distance_rows = per_read.loc[
        per_read["truth_dist_to_5p"].notna()
        & per_read["pred_dist_to_5p"].notna()
        & per_read["truth_dist_to_3p"].notna()
        & per_read["pred_dist_to_3p"].notna()
    ].copy()
    if distance_rows.empty:
        dist_5p_mae = math.nan
        dist_3p_mae = math.nan
        coverage_mae = math.nan
    else:
        dist_5p_mae = float(
            (distance_rows["truth_dist_to_5p"] - distance_rows["pred_dist_to_5p"]).abs().mean()
        )
        dist_3p_mae = float(
            (distance_rows["truth_dist_to_3p"] - distance_rows["pred_dist_to_3p"]).abs().mean()
        )
        coverage_mae = float(
            (
                distance_rows["truth_coverage_fraction"]
                - distance_rows["pred_coverage_fraction"]
            )
            .abs()
            .mean()
        )

    status_counts = Counter(per_read["pred_status"])
    truth_counts = Counter(per_read["truth_status"])
    row = {
        "parameter_set": parameters.parameter_set,
        "min_overlap": parameters.min_overlap,
        "junction_tol": parameters.junction_tol,
        "unique_threshold": parameters.unique_threshold,
        "margin_threshold": parameters.margin_threshold,
        "min_unspliced_coverage_for_unique": parameters.min_unspliced_coverage_for_unique,
        "tss_tol": parameters.tss_tol,
        "tes_tol": parameters.tes_tol,
        "coverage_threshold": parameters.coverage_threshold,
        "read_count": len(per_read),
        "status_accuracy": float(per_read["status_correct"].mean()),
        "unique_transcript_accuracy": unique_transcript_accuracy,
        "full_length_precision": precision,
        "full_length_recall": recall,
        "full_length_f1": f1,
        "dist_to_5p_mae": dist_5p_mae,
        "dist_to_3p_mae": dist_3p_mae,
        "coverage_fraction_mae": coverage_mae,
        "pred_unique_reads": status_counts["unique"],
        "pred_ambiguous_reads": status_counts["ambiguous"],
        "pred_low_confidence_reads": status_counts["low_confidence"],
        "pred_unassigned_reads": status_counts["unassigned"],
        "truth_unique_reads": truth_counts["unique"],
        "truth_ambiguous_reads": truth_counts["ambiguous"],
        "truth_low_confidence_reads": truth_counts["low_confidence"],
        "truth_unassigned_reads": truth_counts["unassigned"],
        "full_length_tp": tp,
        "full_length_fp": fp,
        "full_length_fn": fn,
    }
    for status in ["unique", "ambiguous", "low_confidence", "unassigned"]:
        mask = per_read["truth_status"] == status
        row[f"{status}_recall"] = safe_divide(
            int((mask & (per_read["pred_status"] == status)).sum()),
            int(mask.sum()),
        )
    return row


def compute_coverage_profiles(
    reads: list[ReadAlignment],
    transcripts: dict[str, Transcript],
    truth: pd.DataFrame,
    default_per_read: pd.DataFrame,
    *,
    bin_num: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if bin_num < 1:
        raise ValueError("--coverage-bin-num must be >= 1")

    read_by_id = {read.read_id: read for read in reads}
    truth_profiles = {
        "truth_read_centric": [0.0] * bin_num,
        "isocomp_unique_assignment": [0.0] * bin_num,
        "rseqc_complete_annotation": [0.0] * bin_num,
        "rseqc_truth_expressed_annotation": [0.0] * bin_num,
        "rseqc_redundant_annotation": [0.0] * bin_num,
    }
    projection_counts = {name: 0 for name in truth_profiles}
    redundant_transcripts = build_redundant_transcripts(transcripts)

    truth_unique = truth.loc[truth["truth_status"] == "unique"]
    for row in truth_unique.itertuples(index=False):
        read = read_by_id.get(row.read_id)
        transcript = transcripts.get(row.truth_transcript_id)
        if read is None or transcript is None:
            continue
        if add_projected_read_to_profile(
            truth_profiles["truth_read_centric"],
            read,
            transcript,
        ):
            projection_counts["truth_read_centric"] += 1

    for row in default_per_read.itertuples(index=False):
        if row.pred_status != "unique":
            continue
        read = read_by_id.get(row.read_id)
        transcript = transcripts.get(row.pred_transcript_id)
        if read is None or transcript is None:
            continue
        if add_projected_read_to_profile(
            truth_profiles["isocomp_unique_assignment"],
            read,
            transcript,
        ):
            projection_counts["isocomp_unique_assignment"] += 1

    truth_expressed_transcripts = {
        row.truth_transcript_id
        for row in truth.itertuples(index=False)
        if row.truth_status in {"unique", "ambiguous"}
        and isinstance(row.truth_transcript_id, str)
        and row.truth_transcript_id
    }
    truth_expressed_models = [
        transcript
        for transcript_id, transcript in transcripts.items()
        if transcript_id in truth_expressed_transcripts
    ]

    for read in reads:
        for transcript in transcripts.values():
            if add_projected_read_to_profile(
                truth_profiles["rseqc_complete_annotation"],
                read,
                transcript,
            ):
                projection_counts["rseqc_complete_annotation"] += 1
        for transcript in truth_expressed_models:
            if add_projected_read_to_profile(
                truth_profiles["rseqc_truth_expressed_annotation"],
                read,
                transcript,
            ):
                projection_counts["rseqc_truth_expressed_annotation"] += 1
        for transcript in redundant_transcripts.values():
            if add_projected_read_to_profile(
                truth_profiles["rseqc_redundant_annotation"],
                read,
                transcript,
            ):
                projection_counts["rseqc_redundant_annotation"] += 1

    profile_rows: list[dict[str, Any]] = []
    mean_normalized = {
        name: mean_normalize(values)
        for name, values in truth_profiles.items()
    }
    max_normalized = {
        name: max_normalize(values)
        for name, values in truth_profiles.items()
    }
    for index in range(bin_num):
        row: dict[str, Any] = {"bin": index + 1}
        for name, values in truth_profiles.items():
            row[f"{name}_raw"] = values[index]
            row[f"{name}_mean_normalized"] = mean_normalized[name][index]
            row[f"{name}_max_normalized"] = max_normalized[name][index]
        profile_rows.append(row)

    truth_mean_normalized = mean_normalized["truth_read_centric"]
    truth_max_normalized = max_normalized["truth_read_centric"]
    truth_raw = truth_profiles["truth_read_centric"]
    metrics_rows = []
    for name, values in truth_profiles.items():
        values_mean_normalized = mean_normalized[name]
        values_max_normalized = max_normalized[name]
        metrics_rows.append(
            {
                "profile": name,
                "projection_events": projection_counts[name],
                "total_coverage": sum(values),
                "total_coverage_ratio_vs_truth": safe_divide(sum(values), sum(truth_raw)),
                "mae_mean_normalized_vs_truth": mean_absolute_error(
                    values_mean_normalized,
                    truth_mean_normalized,
                ),
                "rmse_mean_normalized_vs_truth": root_mean_squared_error(
                    values_mean_normalized,
                    truth_mean_normalized,
                ),
                "pearson_mean_normalized_vs_truth": pearson_correlation(
                    values_mean_normalized,
                    truth_mean_normalized,
                ),
                "mae_max_normalized_vs_truth": mean_absolute_error(
                    values_max_normalized,
                    truth_max_normalized,
                ),
                "rmse_max_normalized_vs_truth": root_mean_squared_error(
                    values_max_normalized,
                    truth_max_normalized,
                ),
                "pearson_max_normalized_vs_truth": pearson_correlation(
                    values_max_normalized,
                    truth_max_normalized,
                ),
                "first_decile_mean_normalized": decile_mean(
                    values_mean_normalized,
                    0,
                ),
                "last_decile_mean_normalized": decile_mean(
                    values_mean_normalized,
                    9,
                ),
                "first_to_last_decile_ratio": safe_divide(
                    decile_mean(values_mean_normalized, 0),
                    decile_mean(values_mean_normalized, 9),
                ),
                "first_decile_max_normalized": decile_mean(values_max_normalized, 0),
                "last_decile_max_normalized": decile_mean(values_max_normalized, 9),
                "first_to_last_decile_ratio_max_normalized": safe_divide(
                    decile_mean(values_max_normalized, 0),
                    decile_mean(values_max_normalized, 9),
                ),
                "transcript_models_scanned": (
                    len(transcripts)
                    if name == "rseqc_complete_annotation"
                    else len(truth_expressed_models)
                    if name == "rseqc_truth_expressed_annotation"
                    else len(redundant_transcripts)
                    if name == "rseqc_redundant_annotation"
                    else math.nan
                ),
            }
        )

    return pd.DataFrame(profile_rows), pd.DataFrame(metrics_rows)


def add_projected_read_to_profile(
    bins: list[float],
    read: ReadAlignment,
    transcript: Transcript,
) -> bool:
    projection = project_blocks_to_transcript(read.blocks, transcript)
    if not projection.intervals:
        return False
    for interval in projection.intervals:
        add_interval_to_bins(bins, interval, transcript.transcript_length)
    return True


def add_interval_to_bins(
    bins: list[float],
    interval: tuple[int, int],
    transcript_length: int,
) -> None:
    if transcript_length <= 0:
        return
    start, end = interval
    if end <= start:
        return
    bin_num = len(bins)
    first_bin = max(0, min(bin_num - 1, int(start * bin_num / transcript_length)))
    last_bin = max(0, min(bin_num - 1, int((end - 1) * bin_num / transcript_length)))
    for bin_index in range(first_bin, last_bin + 1):
        bin_start = bin_index * transcript_length / bin_num
        bin_end = (bin_index + 1) * transcript_length / bin_num
        overlap = max(0.0, min(end, bin_end) - max(start, bin_start))
        bin_width = bin_end - bin_start
        if bin_width > 0:
            bins[bin_index] += overlap / bin_width


def mean_normalize(values: list[float]) -> list[float]:
    mean_value = sum(values) / len(values) if values else 0.0
    if mean_value <= 0:
        return [0.0 for _ in values]
    return [value / mean_value for value in values]


def max_normalize(values: list[float]) -> list[float]:
    max_value = max(values) if values else 0.0
    if max_value <= 0:
        return [0.0 for _ in values]
    return [value / max_value for value in values]


def mean_absolute_error(values: list[float], truth_values: list[float]) -> float:
    if not values:
        return math.nan
    return sum(abs(value - truth) for value, truth in zip(values, truth_values)) / len(values)


def root_mean_squared_error(values: list[float], truth_values: list[float]) -> float:
    if not values:
        return math.nan
    return math.sqrt(
        sum((value - truth) ** 2 for value, truth in zip(values, truth_values))
        / len(values)
    )


def pearson_correlation(values: list[float], truth_values: list[float]) -> float:
    if not values:
        return math.nan
    mean_value = sum(values) / len(values)
    mean_truth = sum(truth_values) / len(truth_values)
    numerator = sum(
        (value - mean_value) * (truth - mean_truth)
        for value, truth in zip(values, truth_values)
    )
    value_var = sum((value - mean_value) ** 2 for value in values)
    truth_var = sum((truth - mean_truth) ** 2 for truth in truth_values)
    denominator = math.sqrt(value_var * truth_var)
    return safe_divide(numerator, denominator)


def decile_mean(values: list[float], decile_index: int) -> float:
    if not values:
        return math.nan
    start = int(decile_index * len(values) / 10)
    end = int((decile_index + 1) * len(values) / 10)
    subset = values[start:end]
    if not subset:
        return math.nan
    return sum(subset) / len(subset)


def safe_divide(numerator: float, denominator: float) -> float:
    if denominator == 0:
        return math.nan
    return float(numerator / denominator)


def plot_status_confusion(per_read: pd.DataFrame, path: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipping status confusion plot")
        return

    statuses = ["unique", "ambiguous", "low_confidence", "unassigned"]
    matrix = pd.crosstab(
        per_read["truth_status"],
        per_read["pred_status"],
    ).reindex(index=statuses, columns=statuses, fill_value=0)

    fig, ax = plt.subplots(figsize=(6.8, 5.8))
    image = ax.imshow(matrix.values, cmap="Blues")
    ax.set_xticks(range(len(statuses)), statuses, rotation=35, ha="right")
    ax.set_yticks(range(len(statuses)), statuses)
    ax.set_xlabel("Predicted status")
    ax.set_ylabel("Truth status")
    ax.set_title("Synthetic status confusion matrix")
    for row_index, row in enumerate(matrix.values):
        for col_index, value in enumerate(row):
            ax.text(col_index, row_index, str(value), ha="center", va="center")
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def plot_sensitivity_f1(summary: pd.DataFrame, path: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipping sensitivity F1 plot")
        return

    data = summary.loc[summary["parameter_set"] != "default"].copy()
    if data.empty:
        return
    grouped = (
        data.groupby(["tss_tol", "coverage_threshold"], as_index=False)
        .agg(full_length_f1=("full_length_f1", "mean"))
        .sort_values(["coverage_threshold", "tss_tol"])
    )

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for coverage_threshold, subset in grouped.groupby("coverage_threshold"):
        ax.plot(
            subset["tss_tol"],
            subset["full_length_f1"],
            marker="o",
            label=f"coverage >= {coverage_threshold:g}",
        )
    ax.set_xlabel("Terminal tolerance (bp)")
    ax.set_ylabel("Mean full-length-like F1")
    ax.set_title("Synthetic threshold sensitivity")
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def plot_coverage_profiles(
    profiles: pd.DataFrame,
    path: Path,
    *,
    normalization: str = "max",
) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipping RSeQC-style coverage profile plot")
        return

    profile_names = [
        "truth_read_centric",
        "isocomp_unique_assignment",
        "rseqc_truth_expressed_annotation",
        "rseqc_complete_annotation",
        "rseqc_redundant_annotation",
    ]
    labels = {
        "truth_read_centric": "Truth read-centric",
        "isocomp_unique_assignment": "IsoComp unique assignments",
        "rseqc_truth_expressed_annotation": "RSeQC-style expressed-only annotation",
        "rseqc_complete_annotation": "RSeQC-style complete annotation",
        "rseqc_redundant_annotation": "RSeQC-style redundant annotation",
    }
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for profile_name in profile_names:
        column = f"{profile_name}_{normalization}_normalized"
        ax.plot(
            profiles["bin"],
            profiles[column],
            linewidth=2,
            label=labels[profile_name],
        )
    ax.set_xlabel("Transcript body bin, 5' to 3'")
    ax.set_ylabel(f"{normalization.capitalize()}-normalized coverage")
    ax.set_title("Synthetic IsoComp vs RSeQC-style coverage profiles")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def write_article_panel_outputs(
    *,
    default_metrics: pd.DataFrame,
    sensitivity_summary: pd.DataFrame,
    coverage_profiles: pd.DataFrame,
    coverage_metrics: pd.DataFrame,
    out_dir: Path,
    normalization: str,
    write_plots: bool = True,
) -> None:
    truth_metrics = article_truth_metrics(default_metrics)
    truth_metrics.to_csv(
        out_dir / "figure1a_isocomp_truth_metrics.tsv",
        sep="\t",
        index=False,
    )

    sensitivity_panel = article_sensitivity_panel(sensitivity_summary)
    sensitivity_panel.to_csv(
        out_dir / "figure1c_parameter_sensitivity.tsv",
        sep="\t",
        index=False,
    )

    profile_panel = article_annotation_dependence_profiles(
        coverage_profiles,
        normalization=normalization,
    )
    profile_panel.to_csv(
        out_dir / "figure1b_rseqc_style_annotation_dependence_profiles.tsv",
        sep="\t",
        index=False,
    )

    metrics_panel = article_annotation_dependence_metrics(
        coverage_metrics,
        normalization=normalization,
    )
    metrics_panel.to_csv(
        out_dir / "figure1b_rseqc_style_annotation_dependence_metrics.tsv",
        sep="\t",
        index=False,
    )

    table_s2 = supplementary_table_s2(
        truth_metrics=truth_metrics,
        sensitivity_summary=sensitivity_summary,
        annotation_metrics=metrics_panel,
    )
    table_s2.to_csv(
        out_dir / "supplementary_table_s2_synthetic_summary.tsv",
        sep="\t",
        index=False,
    )

    if not write_plots:
        return
    plot_truth_metrics(
        truth_metrics,
        out_dir / "figure1a_isocomp_truth_metrics.png",
    )
    plot_sensitivity_heatmap(
        sensitivity_panel,
        out_dir / "figure1c_parameter_sensitivity_heatmap.png",
    )
    plot_article_annotation_dependence(
        profile_panel,
        out_dir / "figure1b_rseqc_style_annotation_dependence.png",
        normalization=normalization,
    )


def article_truth_metrics(default_metrics: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        ("assignment status accuracy", "status_accuracy"),
        ("unique transcript accuracy", "unique_transcript_accuracy"),
        ("full-length precision", "full_length_precision"),
        ("full-length recall", "full_length_recall"),
        ("full-length F1", "full_length_f1"),
    ]
    if default_metrics.empty:
        return pd.DataFrame(columns=["metric", "value"])
    row = default_metrics.iloc[0]
    return pd.DataFrame(
        [{"metric": label, "value": float(row[column])} for label, column in metrics]
    )


def article_sensitivity_panel(summary: pd.DataFrame) -> pd.DataFrame:
    defaults = {
        "min_overlap": 50,
        "junction_tol": 5,
        "unique_threshold": 0.8,
        "margin_threshold": 0.1,
        "min_unspliced_coverage_for_unique": 0.2,
    }
    data = summary.copy()
    for column, value in defaults.items():
        data = data.loc[data[column] == value]
    if data.empty:
        data = summary.copy()
    return (
        data.groupby(["tss_tol", "tes_tol", "coverage_threshold"], as_index=False)
        .agg(
            full_length_precision=("full_length_precision", "mean"),
            full_length_recall=("full_length_recall", "mean"),
            full_length_f1=("full_length_f1", "mean"),
        )
        .sort_values(["coverage_threshold", "tss_tol", "tes_tol"])
        .rename(columns={"tss_tol": "terminal_tolerance_bp"})
    )


def article_annotation_dependence_profiles(
    profiles: pd.DataFrame,
    *,
    normalization: str,
) -> pd.DataFrame:
    profile_names = [
        "truth_read_centric",
        "isocomp_unique_assignment",
        "rseqc_truth_expressed_annotation",
        "rseqc_complete_annotation",
        "rseqc_redundant_annotation",
    ]
    labels = {
        "truth_read_centric": "Truth read-centric",
        "isocomp_unique_assignment": "IsoComp unique assignment",
        "rseqc_truth_expressed_annotation": "RSeQC-style expressed-only annotation",
        "rseqc_complete_annotation": "RSeQC-style complete annotation",
        "rseqc_redundant_annotation": "RSeQC-style redundant annotation",
    }
    rows = []
    for profile_name in profile_names:
        column = f"{profile_name}_{normalization}_normalized"
        raw_column = f"{profile_name}_raw"
        if column not in profiles:
            continue
        for row in profiles.itertuples(index=False):
            rows.append(
                {
                    "bin": int(row.bin),
                    "profile": profile_name,
                    "label": labels[profile_name],
                    "raw_coverage": getattr(row, raw_column),
                    f"{normalization}_normalized_coverage": getattr(row, column),
                    "normalization": normalization,
                }
            )
    return pd.DataFrame(rows)


def article_annotation_dependence_metrics(
    metrics: pd.DataFrame,
    *,
    normalization: str,
) -> pd.DataFrame:
    rename_map = {
        f"mae_{normalization}_normalized_vs_truth": "mae_vs_truth",
        f"rmse_{normalization}_normalized_vs_truth": "rmse_vs_truth",
        f"pearson_{normalization}_normalized_vs_truth": "pearson_vs_truth",
        f"first_decile_{normalization}_normalized": "first_decile",
        f"last_decile_{normalization}_normalized": "last_decile",
        (
            f"first_to_last_decile_ratio_{normalization}_normalized"
            if normalization == "max"
            else "first_to_last_decile_ratio"
        ): "first_to_last_decile_ratio",
    }
    base_columns = ["profile", "transcript_models_scanned", "projection_events"]
    selected_columns = base_columns + [
        column for column in rename_map if column in metrics.columns
    ]
    return metrics.loc[:, selected_columns].rename(columns=rename_map)


def supplementary_table_s2(
    *,
    truth_metrics: pd.DataFrame,
    sensitivity_summary: pd.DataFrame,
    annotation_metrics: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for row in truth_metrics.itertuples(index=False):
        rows.append(
            {
                "section": "synthetic_truth",
                "setting": "default",
                "metric": row.metric,
                "value": row.value,
            }
        )

    if not sensitivity_summary.empty:
        rows.extend(
            [
                {
                    "section": "parameter_sensitivity",
                    "setting": "best_parameter_set",
                    "metric": "full_length_F1",
                    "value": float(sensitivity_summary["full_length_f1"].max()),
                },
                {
                    "section": "parameter_sensitivity",
                    "setting": "worst_parameter_set",
                    "metric": "full_length_F1",
                    "value": float(sensitivity_summary["full_length_f1"].min()),
                },
                {
                    "section": "parameter_sensitivity",
                    "setting": "tested_terminal_tolerance_bp",
                    "metric": "values",
                    "value": ",".join(
                        str(int(value))
                        for value in sorted(sensitivity_summary["tss_tol"].unique())
                    ),
                },
                {
                    "section": "parameter_sensitivity",
                    "setting": "tested_coverage_threshold",
                    "metric": "values",
                    "value": ",".join(
                        f"{value:g}"
                        for value in sorted(
                            sensitivity_summary["coverage_threshold"].unique()
                        )
                    ),
                },
                {
                    "section": "parameter_sensitivity",
                    "setting": "tested_min_unspliced_coverage_for_unique",
                    "metric": "values",
                    "value": ",".join(
                        f"{value:g}"
                        for value in sorted(
                            sensitivity_summary[
                                "min_unspliced_coverage_for_unique"
                            ].unique()
                        )
                    ),
                },
            ]
        )

    for row in annotation_metrics.itertuples(index=False):
        rows.extend(
            [
                {
                    "section": "RSeQC_style_annotation_dependence",
                    "setting": row.profile,
                    "metric": "transcript_models_scanned",
                    "value": row.transcript_models_scanned,
                },
                {
                    "section": "RSeQC_style_annotation_dependence",
                    "setting": row.profile,
                    "metric": "RMSE_vs_truth",
                    "value": getattr(row, "rmse_vs_truth", math.nan),
                },
                {
                    "section": "RSeQC_style_annotation_dependence",
                    "setting": row.profile,
                    "metric": "first_to_last_decile_ratio",
                    "value": getattr(row, "first_to_last_decile_ratio", math.nan),
                },
            ]
        )
    return pd.DataFrame(rows)


def plot_truth_metrics(metrics: pd.DataFrame, path: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipping truth metric plot")
        return
    if metrics.empty:
        return

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.bar(metrics["metric"], metrics["value"], color="#4C78A8")
    ax.set_ylabel("Performance")
    ax.set_ylim(0, 1.05)
    ax.set_title("IsoComp synthetic truth benchmark")
    ax.tick_params(axis="x", rotation=30)
    for label in ax.get_xticklabels():
        label.set_horizontalalignment("right")
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def plot_sensitivity_heatmap(panel: pd.DataFrame, path: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipping sensitivity heatmap")
        return
    if panel.empty:
        return

    matrix = panel.pivot_table(
        index="coverage_threshold",
        columns="terminal_tolerance_bp",
        values="full_length_f1",
    ).sort_index(ascending=False)
    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    image = ax.imshow(matrix.values, vmin=0, vmax=1, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(matrix.columns)), [int(value) for value in matrix.columns])
    ax.set_yticks(range(len(matrix.index)), [f"{value:g}" for value in matrix.index])
    ax.set_xlabel("Terminal tolerance (bp)")
    ax.set_ylabel("Coverage threshold")
    ax.set_title("Full-length-like F1 sensitivity")
    for row_index, coverage_threshold in enumerate(matrix.index):
        for col_index, terminal_tolerance in enumerate(matrix.columns):
            value = matrix.loc[coverage_threshold, terminal_tolerance]
            ax.text(col_index, row_index, f"{value:.2f}", ha="center", va="center")
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label="F1")
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def plot_article_annotation_dependence(
    panel: pd.DataFrame,
    path: Path,
    *,
    normalization: str,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipping annotation-dependence plot")
        return
    if panel.empty:
        return

    y_column = f"{normalization}_normalized_coverage"
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for label, subset in panel.groupby("label", sort=False):
        ax.plot(subset["bin"], subset[y_column], linewidth=2, label=label)
    ax.set_xlabel("Transcript body bin, 5' to 3'")
    ax.set_ylabel(f"{normalization.capitalize()}-normalized coverage")
    ax.set_title("Synthetic RSeQC-style annotation dependence")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
