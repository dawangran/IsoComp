"""Command-line entry point for IsoComp."""

from __future__ import annotations

import csv
import logging
import os
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Annotated

import numpy as np
import pandas as pd
import typer

from . import __version__
from .annotation import AnnotationError, infer_annotation_format, read_annotation
from .assigner import assign_read
from .bam_parser import BamParserError, iter_read_alignments, read_bam_reference_names
from .candidate import CandidateIndex, infer_library_strandness
from .io import (
    OutputError,
    ensure_outputs_available,
    output_paths,
    write_dataframe,
    write_json,
    write_transcript_body_coverage,
)
from .metrics import (
    READ_ASSIGNMENT_COLUMNS,
    SampleMetricAccumulator,
    add_assignment_to_body_coverage,
    build_read_metrics,
    init_transcript_accumulators,
    read_metrics_to_row,
    summarize_sample_accumulator,
    summarize_transcript_accumulators,
    update_sample_accumulator,
    update_summary_for_metric,
    update_transcript_accumulators,
)
from .models import RunSummary

app = typer.Typer(add_completion=False, help="Read-centric isoform completeness QC.")
LOGGER = logging.getLogger("isocomp")


@app.callback(invoke_without_command=True)
def run(
    ctx: typer.Context,
    bam: Annotated[Path | None, typer.Option("--bam", help="Genome-aligned long-read RNA BAM.")] = None,
    annotation: Annotated[Path | None, typer.Option("--annotation", help="BED12 or GTF transcript annotation.")] = None,
    annotation_format: Annotated[
        str,
        typer.Option("--annotation-format", help="One of auto, bed12, gtf."),
    ] = "auto",
    out: Annotated[str | None, typer.Option("--out", help="Output prefix.")] = None,
    bin_num: Annotated[int, typer.Option("--bin-num", min=1, help="Number of transcript body bins.")] = 100,
    min_mapq: Annotated[int, typer.Option("--min-mapq", min=0, help="Minimum mapping quality.")] = 20,
    tss_tol: Annotated[int, typer.Option("--tss-tol", min=0, help="5' completeness tolerance in bp.")] = 100,
    tes_tol: Annotated[int, typer.Option("--tes-tol", min=0, help="3' completeness tolerance in bp.")] = 100,
    min_overlap: Annotated[int, typer.Option("--min-overlap", min=1, help="Minimum exonic overlap in bp.")] = 50,
    junction_tol: Annotated[int, typer.Option("--junction-tol", min=0, help="Splice-site tolerance in bp.")] = 5,
    unique_threshold: Annotated[
        float,
        typer.Option("--unique-threshold", min=0.0, max=1.0, help="Minimum top assignment score."),
    ] = 0.8,
    margin_threshold: Annotated[
        float,
        typer.Option("--margin-threshold", min=0.0, max=1.0, help="Minimum top-minus-second score margin."),
    ] = 0.1,
    full_length_coverage: Annotated[
        float,
        typer.Option("--full-length-coverage", min=0.0, max=1.0, help="Minimum transcript coverage for a full-length-like read."),
    ] = 0.8,
    min_terminal_anchor: Annotated[
        int,
        typer.Option("--min-terminal-anchor", min=0, help="Minimum aligned bases supporting each terminal call."),
    ] = 10,
    min_unspliced_coverage_for_unique: Annotated[
        float,
        typer.Option(
            "--min-unspliced-coverage-for-unique",
            min=0.0,
            max=1.0,
            help="Minimum transcript coverage for a junctionless read to be unique on a spliced transcript.",
        ),
    ] = 0.2,
    strandness: Annotated[
        str,
        typer.Option("--strandness", help="One of unstranded, forward, reverse, auto."),
    ] = "unstranded",
    threads: Annotated[int, typer.Option("--threads", min=1, help="BAM decompression/read threads.")] = 1,
    log_level: Annotated[
        str,
        typer.Option("--log-level", help="One of DEBUG, INFO, WARNING, ERROR."),
    ] = "INFO",
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing outputs.")] = False,
    version: Annotated[bool, typer.Option("--version", help="Show version and exit.")] = False,
) -> None:
    if version:
        typer.echo(__version__)
        raise typer.Exit()
    if ctx.invoked_subcommand is not None:
        return
    if bam is None or annotation is None or out is None:
        typer.echo(ctx.get_help())
        raise typer.Exit(code=2)

    try:
        configure_logging(log_level)
        run_pipeline(
            bam=bam,
            annotation=annotation,
            annotation_format=annotation_format,
            out=out,
            bin_num=bin_num,
            min_mapq=min_mapq,
            tss_tol=tss_tol,
            tes_tol=tes_tol,
            min_overlap=min_overlap,
            junction_tol=junction_tol,
            unique_threshold=unique_threshold,
            margin_threshold=margin_threshold,
            full_length_coverage=full_length_coverage,
            min_terminal_anchor=min_terminal_anchor,
            min_unspliced_coverage_for_unique=min_unspliced_coverage_for_unique,
            strandness=strandness,
            threads=threads,
            force=force,
        )
    except (AnnotationError, BamParserError, OutputError, FileNotFoundError, ValueError) as exc:
        LOGGER.error("%s", exc)
        raise typer.Exit(code=1) from exc


def run_pipeline(
    *,
    bam: Path,
    annotation: Path,
    annotation_format: str,
    out: str,
    bin_num: int,
    min_mapq: int,
    tss_tol: int,
    tes_tol: int,
    min_overlap: int,
    junction_tol: int,
    unique_threshold: float,
    margin_threshold: float,
    full_length_coverage: float,
    min_terminal_anchor: int,
    min_unspliced_coverage_for_unique: float,
    strandness: str,
    threads: int,
    force: bool,
) -> None:
    annotation_format = annotation_format.lower()
    if strandness not in {"unstranded", "forward", "reverse", "auto"}:
        raise ValueError(
            "Invalid --strandness. Expected one of unstranded, forward, reverse, auto; "
            f"got {strandness!r}"
        )
    if annotation_format not in {"auto", "bed12", "gtf"}:
        raise ValueError(
            "Invalid --annotation-format. Expected one of auto, bed12, gtf; "
            f"got {annotation_format!r}"
        )
    if not 0 <= min_unspliced_coverage_for_unique <= 1:
        raise ValueError(
            "--min-unspliced-coverage-for-unique must be between 0 and 1; "
            f"got {min_unspliced_coverage_for_unique!r}"
        )

    paths = output_paths(out)
    ensure_outputs_available(paths, force=force)
    LOGGER.info("Reading annotation: %s", annotation)
    resolved_annotation_format = (
        infer_annotation_format(annotation)
        if annotation_format == "auto"
        else annotation_format
    )
    transcripts = read_annotation(annotation, resolved_annotation_format)
    LOGGER.info("Loaded %d transcripts", len(transcripts))
    annotation_references = {transcript.chrom for transcript in transcripts.values()}
    bam_references = read_bam_reference_names(bam)
    shared_references = annotation_references & bam_references
    if not shared_references:
        raise BamParserError(
            "BAM and annotation chromosome names do not overlap; "
            "check conventions such as 'chr1' versus '1'"
        )
    missing_annotation_references = annotation_references - bam_references
    if missing_annotation_references:
        LOGGER.warning(
            "%d annotation reference(s) are absent from the BAM header",
            len(missing_annotation_references),
        )
    candidate_index = CandidateIndex(transcripts)

    resolved_strandness = strandness
    strand_inference: dict[str, object] | None = None
    if strandness == "auto":
        inference_iter, _ = iter_read_alignments(
            bam,
            min_mapq=min_mapq,
            threads=threads,
        )
        try:
            inferred = infer_library_strandness(
                inference_iter,
                candidate_index,
                min_overlap=min_overlap,
            )
        finally:
            close = getattr(inference_iter, "close", None)
            if close is not None:
                close()
        resolved_strandness = inferred.resolved_strandness
        strand_inference = {
            "scanned_reads": inferred.scanned_reads,
            "informative_reads": inferred.informative_reads,
            "forward_votes": inferred.forward_votes,
            "reverse_votes": inferred.reverse_votes,
            "dominant_fraction": inferred.dominant_fraction,
        }
        LOGGER.info(
            "Resolved --strandness auto to %s using %d informative reads (dominant fraction %.3f)",
            resolved_strandness,
            inferred.informative_reads,
            inferred.dominant_fraction,
        )

    read_iter, filter_stats = iter_read_alignments(
        bam,
        min_mapq=min_mapq,
        threads=threads,
    )
    summary = RunSummary()
    sample_accumulator = SampleMetricAccumulator()
    transcript_accumulators = init_transcript_accumulators(transcripts)
    aggregate_coverage = np.zeros(bin_num, dtype=float)
    per_transcript_coverage: dict[str, np.ndarray] = {}

    from .plots import PlotData, update_plot_data, write_plots

    plot_data = PlotData()

    LOGGER.info("Streaming BAM: %s", bam)
    tmp_read_assignment: str | None = None
    try:
        with NamedTemporaryFile(
            "wt",
            encoding="utf-8",
            dir=paths["read_assignment"].parent,
            delete=False,
            newline="",
        ) as handle:
            tmp_read_assignment = handle.name
            writer = csv.DictWriter(
                handle,
                fieldnames=READ_ASSIGNMENT_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()

            for read in read_iter:
                candidates = candidate_index.query(
                    read,
                    min_overlap=min_overlap,
                    strandness=resolved_strandness,
                )
                assignment = assign_read(
                    read,
                    candidates,
                    junction_tol=junction_tol,
                    unique_threshold=unique_threshold,
                    margin_threshold=margin_threshold,
                    min_unspliced_coverage_for_unique=min_unspliced_coverage_for_unique,
                )
                metric = build_read_metrics(
                    read,
                    assignment,
                    tss_tol=tss_tol,
                    tes_tol=tes_tol,
                    full_length_coverage=full_length_coverage,
                    min_terminal_anchor=min_terminal_anchor,
                )
                update_summary_for_metric(summary, metric)
                update_sample_accumulator(
                    sample_accumulator,
                    metric,
                    full_length_coverage=full_length_coverage,
                )
                update_transcript_accumulators(transcript_accumulators, metric)
                add_assignment_to_body_coverage(
                    aggregate_coverage,
                    per_transcript_coverage,
                    assignment,
                )
                update_plot_data(
                    plot_data,
                    metric,
                    assignment,
                    bin_num=bin_num,
                )
                writer.writerow(read_metrics_to_row(metric))
        os.replace(tmp_read_assignment, paths["read_assignment"])
    except Exception:
        if tmp_read_assignment is not None:
            try:
                os.unlink(tmp_read_assignment)
            except FileNotFoundError:
                pass
        raise

    summary.total_reads = filter_stats.total_reads
    summary.mapped_reads = filter_stats.mapped_reads
    summary.primary_reads = filter_stats.primary_reads
    summary.duplicate_reads = filter_stats.duplicate_reads
    summary.low_mapq_reads = filter_stats.low_mapq_reads
    summary.empty_alignment_reads = filter_stats.empty_alignment_reads
    summary.usable_reads = filter_stats.usable_reads

    LOGGER.info("Parsed %d usable reads", summary.usable_reads)
    transcript_rows = summarize_transcript_accumulators(
        transcripts,
        transcript_accumulators,
        per_transcript_coverage,
    )
    sample_row = summarize_sample_accumulator(
        sample_accumulator,
        summary,
        sample_name=Path(out).name,
    )

    write_dataframe(paths["transcript_metrics"], pd.DataFrame(transcript_rows))
    write_dataframe(paths["sample_summary"], pd.DataFrame([sample_row]))
    write_transcript_body_coverage(paths["transcript_body_coverage"], aggregate_coverage)
    stats_payload = {
        "version": __version__,
        "parameters": {
            "bam": str(bam),
            "annotation": str(annotation),
            "annotation_format": annotation_format,
            "resolved_annotation_format": resolved_annotation_format,
            "out": out,
            "bin_num": bin_num,
            "min_mapq": min_mapq,
            "tss_tol": tss_tol,
            "tes_tol": tes_tol,
            "min_overlap": min_overlap,
            "junction_tol": junction_tol,
            "unique_threshold": unique_threshold,
            "margin_threshold": margin_threshold,
            "full_length_coverage": full_length_coverage,
            "min_terminal_anchor": min_terminal_anchor,
            "min_unspliced_coverage_for_unique": min_unspliced_coverage_for_unique,
            "strandness": strandness,
            "resolved_strandness": resolved_strandness,
            "strand_inference": strand_inference,
            "threads": threads,
        },
        "counts": sample_row,
    }
    write_json(paths["assignment_stats"], stats_payload)
    write_plots(
        paths["plots_dir"],
        [],
        aggregate_coverage,
        per_transcript_coverage=per_transcript_coverage,
        assignments=[],
        plot_data=plot_data,
        tss_tol=tss_tol,
        tes_tol=tes_tol,
    )
    from .report import write_reports

    write_reports(
        html_path=paths["html_report"],
        pdf_path=paths["pdf_report"],
        sample_row=sample_row,
        stats_payload=stats_payload,
        plots_dir=paths["plots_dir"],
        output_paths=paths,
    )
    LOGGER.info("IsoComp completed: %s", out)


def configure_logging(log_level: str) -> None:
    normalized = log_level.upper()
    if normalized not in {"DEBUG", "INFO", "WARNING", "ERROR"}:
        raise ValueError(
            "Invalid --log-level. Expected one of DEBUG, INFO, WARNING, ERROR; "
            f"got {log_level!r}"
        )
    logging.basicConfig(
        level=getattr(logging, normalized),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )


def main() -> None:
    app()


if __name__ == "__main__":
    main()
