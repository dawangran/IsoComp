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
- `synthetic.redundant.bed12`: complete annotation plus extra overlapping
  unexpressed transcript models for the RSeQC annotation-dependence test.
- `synthetic.annotation_manifest.tsv`: transcript membership across the three
  annotation sets.
- `synthetic.sorted.bam` and `synthetic.sorted.bam.bai`: sorted/indexed synthetic
  alignments.
- `synthetic.truth.tsv`: expected read-level status and completeness labels.
- `synthetic.expected_profiles.tsv`: truth read-centric and RSeQC-style expected
  transcript-body coverage profiles with raw, mean-normalized, and
  max-normalized values.
- `synthetic.expected_profile_metrics.tsv`: profile error and decile-ratio
  metrics against the truth read-centric profile.
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
- `figure1a_isocomp_truth_metrics.tsv` and
  `figure1a_isocomp_truth_metrics.png`: manuscript-ready source data and plot
  for the IsoComp synthetic truth benchmark.
- `figure1b_rseqc_style_annotation_dependence_profiles.tsv`,
  `figure1b_rseqc_style_annotation_dependence_metrics.tsv`, and
  `figure1b_rseqc_style_annotation_dependence.png`: manuscript-ready RSeQC-style
  annotation-dependence source data and plot.
- `figure1c_parameter_sensitivity.tsv` and
  `figure1c_parameter_sensitivity_heatmap.png`: manuscript-ready threshold
  sensitivity source data and heatmap.
- `supplementary_table_s2_synthetic_summary.tsv`: long-format synthetic
  benchmark, sensitivity, and RSeQC-style comparison summary.
- `synthetic_status_confusion.png`: default status confusion matrix.
- `sensitivity_full_length_f1.png`: full-length-like F1 across the sensitivity grid.
- `rseqc_style_coverage_profiles.png`: optional coverage profile plot.

Use these outputs to support Application Note claims about assignment accuracy,
terminal completeness accuracy, threshold robustness, and the distinction between
read-centric completeness and RSeQC-style annotation scanning.

## Real RSeQC Annotation-Dependence Panel

To generate a real RSeQC version of the synthetic annotation-dependence panel,
first create synthetic BAM/BED inputs:

```bash
python benchmarks/generate_synthetic_data.py --out-dir benchmark_runs/article_synthetic_input --replicates 20
```

Then run RSeQC on the three annotation sets and parse the outputs:

```bash
python benchmarks/rseqc_annotation_dependence.py \
  --input-dir benchmark_runs/article_synthetic_input \
  --out-dir benchmark_runs/article_synthetic_rseqc \
  --rseqc-cmd geneBody_coverage.py
```

If RSeQC has already been run by `run_commands.sh`, parse the existing outputs:

```bash
python benchmarks/rseqc_annotation_dependence.py \
  --input-dir benchmark_runs/article_synthetic_input \
  --out-dir benchmark_runs/article_synthetic_input/rseqc_annotation_dependence \
  --rseqc-output-dir benchmark_runs/article_synthetic_input \
  --skip-rseqc
```

This writes:

- `rseqc_annotation_dependence_profiles.tsv`
- `rseqc_annotation_dependence_metrics.tsv`
- `figure1b_real_rseqc_annotation_dependence.png`
