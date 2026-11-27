# IsoComp

**Read-centric isoform completeness assessment for long-read RNA sequencing.**

IsoComp assesses how completely genome-aligned Nanopore and PacBio RNA reads
cover their best-supported annotated transcript isoforms. It reports transcript
assignment, transcript coverage fraction, separate 5′ and 3′ completeness, and
sample-level QC summaries, with publication-friendly plots and HTML/PDF reports.

```text
long read → candidate transcripts → isoform assignment → transcript coordinates → completeness QC
```

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [Inputs and filtering](#inputs-and-filtering)
- [Command-line options](#command-line-options)
- [Assignment and completeness](#assignment-and-completeness)
- [Outputs and interpretation](#outputs-and-interpretation)
- [Plots and reports](#plots-and-reports)
- [Practical notes](#practical-notes)
- [Support and license](#support-and-license)

## Installation

Python **3.10 or newer** is required. Install from the GitHub source repository:

```bash
git clone https://github.com/dawangran/IsoComp.git
cd IsoComp
python -m venv .venv
source .venv/bin/activate
python -m pip install .
```

The environment activation command above is for Linux/macOS shells. Installation
automatically includes `pysam`, `numpy`, `pandas`, `intervaltree`, `matplotlib`,
and `typer`.

Check the installed command:

```bash
isocomp --version
isocomp --help
```

For an editable source installation, use `python -m pip install -e .` instead.

## Quick start

With a genome-aligned long-read BAM and a BED12 transcript annotation:

```bash
isocomp \
  --bam sample.bam \
  --annotation transcripts.bed12 \
  --out sample.isocomp
```

GTF annotations are also supported:

```bash
isocomp \
  --bam sample.bam \
  --annotation transcripts.gtf \
  --annotation-format gtf \
  --out sample.isocomp
```

A run with explicit QC thresholds and four BAM read/decompression threads:

```bash
isocomp \
  --bam sample.bam \
  --annotation transcripts.bed12 \
  --out sample.isocomp \
  --threads 4 \
  --min-mapq 20 \
  --tss-tol 100 \
  --tes-tol 100 \
  --full-length-coverage 0.8 \
  --strandness unstranded \
  --bin-num 100
```

`--out` is a filename prefix. Start with `sample.isocomp.report.html` for a
run-level overview, and use the TSV files for downstream analysis. Existing
outputs are protected by default; add `--force` to overwrite them intentionally.

## Inputs and filtering

**BAM:** use genome-aligned long-read RNA alignments with splice junctions encoded
as CIGAR `N` operations. IsoComp reads the BAM sequentially, so a BAM index is
not required. The parser excludes unmapped, secondary, supplementary,
duplicate-flagged, low-MAPQ, and empty-alignment records. The default minimum MAPQ
is 20. Aligned reference bases from `M`, `=`, and `X` contribute coverage;
insertions, deletions, introns, and clipping do not count as covered bases.

**Annotation:** provide transcript models in BED12 or GTF format, using the same
reference assembly and chromosome naming convention as the BAM.

- BED12 coordinates are 0-based, half-open. The `name` field supplies
  `transcript_id`; BED12 does not supply a gene identifier, so `gene_id` is empty.
- GTF parsing uses `exon` records with `transcript_id` and `gene_id` attributes.
  Coordinates are converted from 1-based, closed to 0-based, half-open.
- For transcripts with a known strand, transcript coordinates run from 5′ to 3′,
  including negative-strand transcripts. Discontinuous projected intervals are
  supported.
- Transcripts with unknown strand (`.`) can contribute assignment and coverage
  fractions, but are excluded from directional terminal metrics and body coverage.
- `--annotation-format auto` uses the file suffix, then inspects the first data
  line when needed. Set `bed12` or `gtf` explicitly if detection is unsuitable.

Read and filtering counts describe BAM alignment records. They are not UMI-based
molecule counts.

## Command-line options

`--bam`, `--annotation`, and `--out` are required. Other options are:

| Option | Default | Purpose |
| --- | --- | --- |
| `--annotation-format` | `auto` | Annotation format: `auto`, `bed12`, or `gtf`. |
| `--bin-num` | `100` | Number of normalized transcript-body bins. |
| `--min-mapq` | `20` | Minimum alignment mapping quality. |
| `--min-overlap` | `50` | Minimum read–transcript exonic overlap, in bp. |
| `--junction-tol` | `5` | Allowed difference at each splice site, in bp. |
| `--unique-threshold` | `0.8` | Minimum top assignment score for an assigned read. |
| `--margin-threshold` | `0.1` | Minimum top-minus-second score margin for a unique assignment. |
| `--tss-tol` | `100` | Maximum distance to the annotated 5′ end for a terminal call, in bp. |
| `--tes-tol` | `100` | Maximum distance to the annotated 3′ end for a terminal call, in bp. |
| `--full-length-coverage` | `0.8` | Minimum coverage fraction for a full-length-like read. |
| `--min-terminal-anchor` | `10` | Minimum covered transcript bases supporting each terminal call. |
| `--min-unspliced-coverage-for-unique` | `0.2` | Minimum coverage for a junctionless read to be unique on a spliced transcript. |
| `--strandness` | `unstranded` | Library orientation: `unstranded`, `forward`, `reverse`, or `auto`. |
| `--threads` | `1` | BAM reading/decompression threads. |
| `--log-level` | `INFO` | Logging level: `DEBUG`, `INFO`, `WARNING`, or `ERROR`. |
| `--force` | Off | Overwrite existing outputs. |
| `--version` | — | Print the installed version and exit. |

`--threads` controls BAM I/O; it does not parallelize transcript assignment.

For strand-aware analysis, `forward` requires the BAM alignment strand to match
the transcript strand, and `reverse` requires the opposite orientation.
`unstranded` does not filter candidates by alignment strand. `auto` examines up
to the first 100,000 usable alignments and resolves a sample-level orientation
when at least 100 informative reads provide at least 80% support for one mode;
otherwise it uses `unstranded`. Requested and resolved settings are recorded in
the JSON output.

## Assignment and completeness

### Read-centric design

Traditional gene-body coverage scans annotation models and aggregates genomic
coverage at normalized positions. Overlapping isoforms and unexpressed models
can reuse the same genomic signal. IsoComp instead finds candidate transcripts
for each read through a chromosome-wise interval index, scores their
compatibility, and projects the selected alignment into transcript coordinates.
This reduces reuse of coverage across redundant annotation models.

Assignment and completeness are separate: a truncated read can still provide
strong support for a transcript. Completeness is then evaluated relative to that
transcript's annotated ends.

### Transcript assignment

The compatibility score is:

```text
assignment_score = 0.65 × exon_overlap_score + 0.35 × junction_score
```

`exon_overlap_score` is the fraction of aligned read reference bases that overlap
the transcript's exons. `junction_score` is the harmonic mean of junction
precision and recall, allowing the configured splice-site tolerance. Recall
considers transcript junctions within the read's genomic span.

| Status | Default interpretation |
| --- | --- |
| `unique` | Top score ≥ 0.8 and score margin ≥ 0.1, with the junctionless-read guard satisfied. |
| `ambiguous` | Top score ≥ 0.8, but score margin < 0.1. |
| `low_confidence` | Top score < 0.8, or an otherwise unique junctionless read covers less than 20% of a spliced transcript. |
| `unassigned` | No transcript passes the candidate filters. |

Set `--min-unspliced-coverage-for-unique 0` to disable the junctionless-read guard.
An ambiguous row reports the top transcript as a representative model; it does
not establish a unique isoform identity. No fractional assignment is performed.

Coverage fraction and terminal distances do not enter the compatibility score.
Coverage can break an exact score tie and is used by the junctionless-read guard.

### Completeness metrics

For a transcript of length `L`, projected intervals define:

```text
coverage_fraction = covered transcript bases / L
dist_to_5p        = read_start_tx
dist_to_3p        = L - read_end_tx
start_pct        = read_start_tx / L
end_pct          = read_end_tx / L
```

Covered bases are the union of projected intervals, rather than the entire span
between the first and last projected base. Despite their names, `start_pct` and
`end_pct` are fractions between 0 and 1.

With the default thresholds, a 5′ or 3′ completeness call requires both a distance
of at most 100 bp and at least 10 covered bases in the corresponding terminal
window. Each window extends 110 bp from the annotated end (tolerance + anchor),
limited by transcript length. For transcripts shorter than 10 bp, the required
anchor is their full length. Set `--min-terminal-anchor 0` for endpoint-only calls.

A read-level full-length-like call requires:

```text
assignment_status == "unique"
and coverage_fraction >= 0.8
and is_5p_complete
and is_3p_complete
```

For example, a uniquely assigned read covering all exons of a 1,000 bp transcript
has coverage fraction 1.0, both terminal distances 0, and is full-length-like.
A read covering transcript coordinates 300–1,000 has coverage fraction 0.7,
5′ distance 300 bp, and 3′ distance 0 bp: it is 3′ complete but not full-length-like.

## Outputs and interpretation

For `--out sample.isocomp`, the main outputs are:

| File | Contents |
| --- | --- |
| `sample.isocomp.read_assignment.tsv` | One row per usable alignment, with assignment status, selected transcript, score, projected coverage, terminal metrics, MAPQ, and clipping. |
| `sample.isocomp.transcript_metrics.tsv` | Transcript-level read counts, completeness summaries, terminal distances, and coverage-uniformity metrics. |
| `sample.isocomp.sample_summary.tsv` | Filtering and assignment counts, unique-read QC, and `all_assigned_*` summaries. |
| `sample.isocomp.transcript_body_coverage.tsv` | Aggregate transcript-normalized coverage, including raw, mean-normalized, and max-normalized values. |
| `sample.isocomp.assignment_stats.json` | Version, run parameters, resolved annotation format/strandness, and summary counts. |
| `sample.isocomp.report.html` | Self-contained QC report with embedded figures. |
| `sample.isocomp.report.pdf` | Multi-page QC report. |
| `sample.isocomp.plots/` | Individual PNG/PDF figures. |

**Read-level results:** use `assignment_status` before interpreting a transcript
identity. Unassigned and low-confidence rows do not carry an assigned transcript
or projected completeness metrics. The `is_full_length_like` flag applies only
to unique reads.

**Transcript-level results:** assigned-read counts include unique and ambiguous
assignments; coverage and completeness summaries use unique reads. Body coverage
also requires a known transcript strand. Each read contributes the fraction of
each normalized bin that its projected intervals cover.

**Sample-level results:** default coverage summaries use unique reads. The
`5p_complete_fraction`, `3p_complete_fraction`, and `full_length_like_fraction`
denominator is `terminal_evaluable_unique_reads`, excluding unknown-strand
transcripts. The `all_assigned_*` summaries additionally include ambiguous reads
using their representative top transcript, with terminal fractions divided by
`terminal_evaluable_assigned_reads`. The all-assigned full-length-like fraction
applies coverage and terminal criteria to this broader set; it is not the sum of
the unique-only read-level full-length flags. Low-confidence reads are counted
separately and do not contribute assigned-transcript evidence.

Inspect assignment counts alongside completeness fractions: a high complete-read
fraction among very few uniquely assigned reads does not describe the whole
library.

## Plots and reports

Each figure is written as a 300 dpi PNG and a PDF:

| Figure | Interpretation |
| --- | --- |
| `transcript_body_coverage` | Aggregate assigned-read coverage from 5′ to 3′, displayed with max normalization on a 0–1 y-axis. |
| `transcript_body_heatmap` | Mean-normalized transcript-by-bin coverage for up to the 200 most-covered transcripts with unique directional evidence. |
| `transcript_body_heatmap_full` | The same transcript heatmap without the 200-transcript limit. |
| `read_body_heatmap` | Projected unique-read coverage, ordered by normalized 5′ then 3′ distance; up to 500 sampled reads. |
| `read_coverage_fraction` | Unique-read transcript coverage distribution with its summary median marked. |
| `dist_to_5p` | Read distance to the annotated 5′ end, with the configured tolerance marked. |
| `dist_to_3p` | Read distance to the annotated 3′ end, with the configured tolerance marked. |
| `full_length_fraction` | 5′ complete, 3′ complete, and full-length-like fractions among terminal-evaluable unique reads. |

The terminal-distance displays use the sampled 95th percentile, normally capped
near 2 kb, and pool longer distances into the final bin. Read the TSV metrics
when exact per-read distances are needed.

The HTML report embeds PNG figures and can be opened without a web server. The
PDF report includes summary metrics, run parameters, and QC plots for sharing.

## Practical notes

- **Reference compatibility:** a complete chromosome-name mismatch, such as
  `chr1` versus `1`, stops analysis. Check both the reference assembly and naming
  convention when assignment rates are unexpectedly low.
- **Annotation dependence:** completeness is relative to the selected annotated
  isoform. Incomplete transcript models, inaccurate ends, and unresolved isoforms
  can change assignment and terminal calls. Full-length-like is a QC label for
  alignment support, rather than proof of an intact original RNA molecule.
- **Direct RNA data:** inspect 5′ and 3′ completeness separately. A single
  full-length-like fraction can conceal differences between the two ends.
- **Comparing libraries:** keep the reference annotation, filtering, assignment,
  and completeness thresholds consistent across samples. Compare assignment
  yield, coverage fraction, and terminal completeness together.
- **Large BAMs:** alignments and read-level output are streamed. Sample-level
  medians are exact up to 100,000 retained values per metric, then use an online
  quantile estimator; transcript-level medians use constant-memory online
  estimation from the first observation. Histogram samples are capped at
  100,000 values and the read heatmap at 500 reads. Transcript-body arrays are
  allocated for transcripts receiving unique directional evidence; the full
  transcript heatmap can still grow with the number of represented transcripts.

## Support and license

Report bugs or request features through [GitHub Issues](https://github.com/dawangran/IsoComp/issues).
Include the IsoComp version, command, relevant log output, and a small
reproducible input when possible.

IsoComp is released under the [MIT License](LICENSE).
