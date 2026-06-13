#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:-example_runs/synthetic_example}"
PYTHON_BIN="${PYTHON:-python}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUT_DIR_ABS="$(cd "${REPO_ROOT}" && mkdir -p "${OUT_DIR}" && cd "${OUT_DIR}" && pwd)"

export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${OUT_DIR_ABS}/matplotlib-cache}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${OUT_DIR_ABS}/xdg-cache}"
mkdir -p "${MPLCONFIGDIR}"
mkdir -p "${XDG_CACHE_HOME}"

cd "${REPO_ROOT}"

"${PYTHON_BIN}" benchmarks/generate_synthetic_data.py \
  --out-dir "${OUT_DIR_ABS}/input" \
  --replicates 2

"${PYTHON_BIN}" -m isocomp.cli \
  --bam "${OUT_DIR_ABS}/input/synthetic.sorted.bam" \
  --annotation "${OUT_DIR_ABS}/input/synthetic.complete.bed12" \
  --out "${OUT_DIR_ABS}/isocomp.synthetic" \
  --bin-num 100 \
  --min-mapq 20 \
  --tss-tol 100 \
  --tes-tol 100 \
  --log-level WARNING \
  --force

"${PYTHON_BIN}" - "${OUT_DIR_ABS}" <<'PY'
from pathlib import Path
import json
import sys

out_dir = Path(sys.argv[1])
summary_path = out_dir / "isocomp.synthetic.sample_summary.tsv"
assignment_path = out_dir / "isocomp.synthetic.read_assignment.tsv"
stats_path = out_dir / "isocomp.synthetic.assignment_stats.json"

print(f"IsoComp synthetic example completed: {out_dir}")
print(f"Sample summary: {summary_path}")
print(f"Read assignments: {assignment_path}")
print(f"Run metadata: {stats_path}")

rows = summary_path.read_text(encoding="utf-8").strip().splitlines()
if len(rows) >= 2:
    header = rows[0].split("\t")
    values = rows[1].split("\t")
    summary = dict(zip(header, values))
    keys = [
        "total_reads",
        "unique_assigned_reads",
        "ambiguous_reads",
        "low_confidence_reads",
        "unassigned_reads",
        "full_length_like_fraction",
    ]
    for key in keys:
        print(f"{key}: {summary.get(key, '')}")

payload = json.loads(stats_path.read_text(encoding="utf-8"))
print(f"IsoComp version: {payload.get('version', 'unknown')}")
PY
