# Software Availability Plan

This document is the publication-facing availability checklist for IsoComp.

## Current Repository Assets

- Python package metadata: `pyproject.toml`
- Command-line entry point: `isocomp`
- Unit tests: `tests/`
- Continuous integration: `.github/workflows/ci.yml`
- Synthetic data generator: `benchmarks/generate_synthetic_data.py`
- Synthetic truth and sensitivity benchmark: `benchmarks/synthetic_truth_benchmark.py`
- One-command reviewer smoke test: `scripts/run_synthetic_example.sh`

## Reviewer Smoke Test

From a clean checkout:

```bash
python -m pip install -e ".[test]"
pytest
bash scripts/run_synthetic_example.sh example_runs/synthetic_example
```

Expected behavior:

- `pytest` passes without external sequencing data.
- `example_runs/synthetic_example/input/` contains a sorted/indexed synthetic BAM,
  BED12 annotation files, and truth labels.
- `example_runs/synthetic_example/isocomp.synthetic.*` contains IsoComp TSV, JSON,
  and plot outputs.

## Release Checklist

Complete these items before manuscript submission:

1. Confirm the version in `pyproject.toml` and `isocomp/__init__.py`.
2. Run `pytest` on Python 3.10, 3.11, and 3.12 or confirm GitHub Actions passed.
3. Run `bash scripts/run_synthetic_example.sh example_runs/synthetic_example`.
4. Run `python benchmarks/synthetic_truth_benchmark.py --out-dir benchmark_runs/synthetic_truth --no-plots`.
5. Build a source and wheel distribution:

   ```bash
   python -m pip install build
   python -m build
   ```

6. Create a signed or annotated git tag, for example:

   ```bash
   git tag -a v0.1.0 -m "IsoComp v0.1.0"
   git push origin v0.1.0
   ```

7. Create a GitHub release from that tag and attach:

   - source archive
   - wheel or source distribution from `dist/`
   - synthetic example outputs or a compact test-data archive

8. Archive the GitHub release or release archive on Zenodo, Figshare, or Software
   Heritage and record the DOI.
9. Replace all manuscript placeholders for contact, funding, release URL, and DOI.

## Manuscript Availability Text Template

Use this as the starting point for the Bioinformatics Availability and
implementation field after creating the release and archive:

```text
IsoComp v0.1.0 is implemented in Python 3.10 or newer and is released under the
MIT License. Source code, installation instructions, unit tests, synthetic test
data, and scripts to reproduce the synthetic benchmark are available at
https://github.com/dawangran/IsoComp. The exact version used in this manuscript
is archived at [Zenodo/Figshare/Software Heritage DOI to be inserted before
submission].
```

Do not submit the manuscript with the DOI placeholder still present.
