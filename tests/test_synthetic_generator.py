from __future__ import annotations

from pathlib import Path

import pandas as pd

from isocomp.annotation import read_bed12
from isocomp.assigner import assign_read
from isocomp.bam_parser import cigar_to_blocks_and_junctions
from isocomp.candidate import CandidateIndex
from isocomp.models import ReadAlignment

from benchmarks import generate_synthetic_data as synth
from benchmarks import rseqc_annotation_dependence as rseqc


def test_all_simulation_mode_adds_realistic_and_stress_reads() -> None:
    transcripts = synth.build_transcripts()

    reads = synth.build_reads(
        transcripts,
        replicates=1,
        mapq=60,
        simulation_mode="all",
    )

    tiers = {read.simulation_tier for read in reads}
    artifacts = {read.artifact for read in reads}
    read_ids = {read.read_id for read in reads}

    assert tiers == {"clean", "realistic", "stress"}
    assert "terminal_softclip" in artifacts
    assert "internal_insertion" in artifacts
    assert "junction_shift" in artifacts
    assert "chimeric_alignment" in artifacts
    assert "novel_exon_chain" in artifacts
    assert "complex_exon_chain" in artifacts
    assert "exon_skipping" in artifacts
    assert "retained_intron" in artifacts
    assert "alternative_tss" in artifacts
    assert "alternative_tes" in artifacts
    assert "softclip_5p_full_000" in read_ids
    assert "chimeric_low_conf_000" in read_ids
    assert "five_exon_full_000" in read_ids
    assert "ten_exon_full_000" in read_ids
    assert "exon_skip_full_000" in read_ids
    assert "retained_intron_full_000" in read_ids
    assert "alt_tss_full_000" in read_ids
    assert "alt_tes_full_000" in read_ids


def test_transcriptome_includes_complex_isoform_structures() -> None:
    transcripts = synth.build_transcripts()

    assert len(transcripts["TX_FIVE_EXON_MAIN"].exons) == 5
    assert len(transcripts["TX_TEN_EXON_MAIN"].exons) == 10
    assert transcripts["TX_FIVE_EXON_SKIP"].exons == (
        (12000, 12100),
        (12200, 12300),
        (12700, 12800),
        (13000, 13100),
    )
    assert transcripts["TX_FIVE_RETAINED_INTRON"].exons == (
        (12000, 12100),
        (12200, 12550),
        (12700, 12800),
        (13000, 13100),
    )
    assert transcripts["TX_FIVE_ALT_TSS"].exons[0] == (12040, 12100)
    assert transcripts["TX_FIVE_ALT_TES"].exons[-1] == (13000, 13060)


def test_terminally_similar_complex_isoforms_are_labeled_ambiguous() -> None:
    transcripts = synth.build_transcripts()

    reads = synth.build_reads(
        transcripts,
        replicates=1,
        mapq=60,
        simulation_mode="clean",
    )
    statuses = {read.read_id: read.truth_status for read in reads}
    full_length_flags = {
        read.read_id: read.truth_is_full_length_like
        for read in reads
    }

    assert statuses["five_exon_full_000"] == "ambiguous"
    assert statuses["five_exon_internal_partial_000"] == "ambiguous"
    assert statuses["alt_tss_full_000"] == "ambiguous"
    assert statuses["alt_tes_full_000"] == "ambiguous"
    assert full_length_flags["five_exon_full_000"] is False


def test_boundary_truth_labels_match_current_assignment_rules(tmp_path: Path) -> None:
    transcripts = synth.build_transcripts()
    bed_path = tmp_path / "synthetic.bed12"
    synth.write_bed12(bed_path, transcripts)
    parsed_transcripts = read_bed12(bed_path)
    index = CandidateIndex(parsed_transcripts)
    reads = {
        read.read_id: read
        for read in synth.build_reads(
            transcripts,
            replicates=1,
            mapq=60,
            simulation_mode="all",
        )
    }

    expected = {
        "junction_shift_low_conf_000": "low_confidence",
        "chimeric_low_conf_000": "low_confidence",
    }
    for read_id, expected_status in expected.items():
        read = reads[read_id]
        assignment = assign_read(
            _to_read_alignment(read),
            index.query(_to_read_alignment(read), min_overlap=1),
        )

        assert read.truth_status == expected_status
        assert assignment.status == expected_status


def test_truth_table_records_simulation_tier_and_artifact(tmp_path: Path) -> None:
    transcripts = synth.build_transcripts()
    reads = synth.build_reads(
        transcripts,
        replicates=1,
        mapq=60,
        simulation_mode="realistic",
    )
    truth_path = tmp_path / "truth.tsv"

    synth.write_truth_tsv(truth_path, reads)

    lines = truth_path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    first_realistic = dict(zip(header, lines[-1].split("\t")))

    assert "simulation_tier" in header
    assert "artifact" in header
    assert first_realistic["simulation_tier"] == "realistic"
    assert first_realistic["artifact"] != "ideal"


def test_realistic_cigar_artifacts_preserve_reference_projection() -> None:
    transcript = synth.build_transcripts()["TX_POS_MAIN"]
    read = synth.read_from_transcript_interval(
        read_id="artifact_read",
        transcript=transcript,
        tx_interval=(0, transcript.transcript_length),
        truth_status="unique",
        truth_is_full_length_like=True,
        scenario="artifact_read",
        mapq=60,
        simulation_tier="realistic",
        artifact="terminal_softclip",
        left_softclip=12,
        right_softclip=7,
        insertion_after_first_match=5,
    )

    assert synth.cigar_to_string(read.cigar).startswith("12S")
    assert "5I" in synth.cigar_to_string(read.cigar)
    assert synth.cigar_to_string(read.cigar).endswith("7S")
    assert synth.query_length(read.cigar) == transcript.transcript_length + 24
    assert synth.project_read_to_transcript(read, transcript) == [
        (0, 100),
        (100, 200),
        (200, 300),
    ]


def test_incomplete_annotation_omits_selected_expressed_isoforms() -> None:
    transcripts = synth.build_transcripts()

    incomplete = synth.build_incomplete_transcripts(transcripts)

    assert "TX_POS_MAIN" in incomplete
    assert "TX_POS_ALT_JUNC" not in incomplete
    assert "TX_TERM_SHORT" not in incomplete


def test_rseqc_expected_profile_rows_include_incomplete_annotation() -> None:
    expected_profiles = pd.DataFrame(
        {
            "bin": [1],
            "truth_read_centric_raw": [1.0],
            "truth_read_centric_max_normalized": [1.0],
            "rseqc_style_incomplete_annotation_raw": [0.5],
            "rseqc_style_incomplete_annotation_max_normalized": [0.5],
        }
    )

    rows = rseqc.build_profile_rows(
        expected_profiles=expected_profiles,
        rseqc_prefixes={},
        normalization="max",
    )

    assert "rseqc_style_incomplete_annotation" in set(rows["profile"])


def _to_read_alignment(read: synth.SyntheticRead) -> ReadAlignment:
    blocks, junctions = cigar_to_blocks_and_junctions(read.start, list(read.cigar))
    return ReadAlignment(
        read_id=read.read_id,
        chrom=read.chrom,
        genomic_start=blocks[0][0],
        genomic_end=blocks[-1][1],
        blocks=blocks,
        junctions=junctions,
        aligned_length=sum(end - start for start, end in blocks),
        mapq=read.mapq,
        cigar=synth.cigar_to_string(read.cigar),
        is_reverse=read.is_reverse,
    )
