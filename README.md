# IsoComp

IsoComp is a read-centric isoform completeness QC tool for long-read RNA sequencing.
It assigns each genome-aligned read to its best supported transcript model, projects
the read into transcript coordinates, and reports 5' completeness, 3' completeness,
coverage fraction, and transcript-normalized body coverage.

IsoComp intentionally avoids the RSeQC-style strategy of scanning every annotated
transcript model and querying BAM coverage. All v0.1 metrics come from reads that
are assigned and projected in transcript coordinates.

## Install

```bash
python -m pip install -e ".[test]"
```

IsoComp requires Python 3.10 or newer.

## Minimal Run

```bash
isocomp \
  --bam sample.bam \
  --annotation transcripts.bed12 \
  --out sample.isocomp
```

Common HPC-friendly options:

```bash
isocomp \
  --bam sample.bam \
  --annotation transcripts.bed12 \
  --out sample.isocomp \
  --threads 4 \
  --min-mapq 20 \
  --bin-num 100 \
  --log-level INFO
```

By default, IsoComp refuses to overwrite existing output files. Use `--force` when
you intentionally want to replace a previous run.

## Inputs

- `--bam`: genome-aligned long-read RNA BAM. IsoComp streams records with
  `pysam.AlignmentFile(...).fetch(until_eof=True)`, so a BAM index is not required.
- `--annotation`: BED12 transcript annotation. Coordinates are interpreted as
  0-based, half-open.
- BED12 `name` is used as `transcript_id` in v0.1. `gene_id` is left empty until
  GTF or explicit transcript-to-gene mapping is added.

## Outputs

For `--out sample.isocomp`, IsoComp writes:

- `sample.isocomp.read_assignment.tsv`
- `sample.isocomp.transcript_metrics.tsv`
- `sample.isocomp.sample_summary.tsv`
- `sample.isocomp.transcript_body_coverage.tsv`
- `sample.isocomp.assignment_stats.json`
- `sample.isocomp.plots/`

Transcript-level metrics and transcript body coverage use uniquely assigned reads
by default. Ambiguous reads are reported in read-level and sample-level outputs but
do not contribute to default transcript-level coverage. The body coverage TSV keeps
raw `coverage`, `mean_normalized_coverage`, and `max_normalized_coverage`.

The plots directory contains:

- `transcript_body_coverage.png`: aggregate line plot using max-normalized coverage
  on a fixed 0-1 y-axis.
- `transcript_body_heatmap.png`: transcript-by-bin heatmap for transcripts with at
  least one unique read; each row is mean-normalized and the display is capped at
  the 200 most-covered transcripts.
- `read_body_heatmap.png`: read-by-bin heatmap for unique reads; each row shows
  the covered fraction of bins in that read's assigned transcript coordinates and
  the display is capped at 500 reads.

## Assignment Defaults

Candidate transcripts are looked up through a chromosome-wise interval index. Each
candidate is scored as:

```text
0.50 * exon_overlap_score + 0.30 * junction_score + 0.20 * coverage_fraction
```

Statuses:

- `unique`: top score is at least `0.8` and exceeds second-best by at least `0.1`
- `ambiguous`: top score is at least `0.8` but the margin is below `0.1`
- `low_confidence`: candidates exist but top score is below `0.8`
- `unassigned`: no candidate passes the candidate filters

## Development

```bash
python -m pip install -e ".[test]"
pytest
```

The tests use tiny synthetic BED12/BAM fixtures and do not require real sequencing
data.
