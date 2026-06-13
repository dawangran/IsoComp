# Synthetic Truth Benchmark

This directory contains manuscript-oriented benchmark helpers for IsoComp.

## Manual Synthetic Data

To generate only the synthetic BED12/BAM/truth files and test IsoComp or RSeQC
yourself:

```bash
python benchmarks/generate_synthetic_data.py --out-dir benchmark_runs/manual_synthetic
```

This writes:

- `synthetic.complete.bed12`: complete annotation with expressed and decoy
  transcripts.
- `synthetic.expressed_only.bed12`: expressed-only matched annotation control.
- `synthetic.sorted.bam` and `synthetic.sorted.bam.bai`: sorted/indexed synthetic
  alignments.
- `synthetic.truth.tsv`: expected read-level status and completeness labels.
- `run_commands.sh`: example commands for IsoComp and RSeQC.

The generator requires `pysam` to write BAM files. Running it inside an RSeQC
environment is sufficient because RSeQC installs `pysam`.

## Synthetic Truth Benchmark

Run the synthetic truth benchmark and parameter sensitivity analysis from the
repository root:

```bash
python benchmarks/synthetic_truth_benchmark.py --out-dir benchmark_runs/synthetic_truth
```

By default, the script runs the benchmark in memory against IsoComp's core
assignment/projection functions, so it does not require `pysam`. Add
`--write-bam` when you want a physical BAM fixture for CLI or BAM parser checks.

The script writes:

- `synthetic_truth.bed12`: synthetic isoform annotation.
- `synthetic_truth.bam`: synthetic genome-aligned long-read BAM when `--write-bam`
  is used.
- `synthetic_truth.tsv`: per-read truth labels.
- `default.per_read.tsv`: default-parameter per-read predictions.
- `default.metrics.tsv`: default-parameter benchmark metrics.
- `sensitivity_summary.tsv`: one row per parameter combination.
- `sensitivity_per_read.tsv`: per-read predictions for every parameter combination.
- `rseqc_style_coverage_profiles.tsv`: truth, IsoComp, and RSeQC-style
  transcript-body coverage profiles on the same synthetic reads.
- `rseqc_style_coverage_metrics.tsv`: profile error, correlation, and coverage
  inflation summaries for the RSeQC-style comparison.
- `synthetic_status_confusion.png`: default status confusion matrix.
- `sensitivity_full_length_f1.png`: full-length-like F1 across the sensitivity grid.
- `rseqc_style_coverage_profiles.png`: optional coverage profile plot.

Use these outputs to support Application Note claims about assignment accuracy,
terminal completeness accuracy, threshold robustness, and the distinction between
read-centric completeness and RSeQC-style annotation scanning.
