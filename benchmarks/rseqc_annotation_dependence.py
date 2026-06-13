#!/usr/bin/env python3
"""Run and summarize real RSeQC geneBodyCoverage on synthetic annotations.

This script is manuscript-oriented. It compares the same synthetic BAM under
three annotation inputs: truth-expressed, complete, and redundant. The output is
intended for the RSeQC annotation-dependence panel and its source table.
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
RSEQC_LABELS = {
    "expressed_only_annotation": "rseqc.expressed_only_annotation",
    "complete_annotation": "rseqc.complete_annotation",
    "redundant_annotation": "rseqc.redundant_annotation",
}
LEGACY_RSEQC_LABELS = {
    "expressed_only_annotation": "rseqc.expressed_only",
    "complete_annotation": "rseqc.complete",
    "redundant_annotation": "rseqc.redundant",
}


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    configure_plot_cache(args.out_dir)

    input_dir = args.input_dir
    if input_dir is None:
        input_dir = args.out_dir / "input"
        generate_synthetic_input(input_dir, args.replicates, args.python)
    input_dir = input_dir.resolve()

    bam_path = input_dir / "synthetic.sorted.bam"
    annotations = {
        "expressed_only_annotation": input_dir / "synthetic.expressed_only.bed12",
        "complete_annotation": input_dir / "synthetic.complete.bed12",
        "redundant_annotation": input_dir / "synthetic.redundant.bed12",
    }
    missing_inputs = [
        str(path)
        for path in [bam_path, *annotations.values()]
        if not path.exists()
    ]
    if missing_inputs:
        raise SystemExit(
            "Missing synthetic input files. Generate them with "
            "`python benchmarks/generate_synthetic_data.py --out-dir ...`. "
            f"Missing: {', '.join(missing_inputs)}"
        )

    rseqc_output_dir = (args.rseqc_output_dir or args.out_dir).resolve()
    rseqc_output_dir.mkdir(parents=True, exist_ok=True)
    rseqc_prefixes = {
        label: rseqc_output_dir / RSEQC_LABELS[label]
        for label in annotations
    }
    if not args.skip_rseqc:
        rseqc_cmd = resolve_command(args.rseqc_cmd)
        for label, annotation_path in annotations.items():
            run_rseqc(
                rseqc_cmd=rseqc_cmd,
                bam_path=bam_path,
                annotation_path=annotation_path,
                output_prefix=rseqc_prefixes[label],
                image_format=args.image_format,
            )
    else:
        rseqc_prefixes = {
            label: find_existing_rseqc_prefix(label, rseqc_output_dir)
            for label in annotations
        }

    expected_profiles = load_expected_profiles(input_dir / "synthetic.expected_profiles.tsv")
    profile_rows = build_profile_rows(
        expected_profiles=expected_profiles,
        rseqc_prefixes=rseqc_prefixes,
        normalization=args.normalization,
    )
    profile_rows.to_csv(
        args.out_dir / "rseqc_annotation_dependence_profiles.tsv",
        sep="\t",
        index=False,
    )

    metric_rows = summarize_profile_rows(profile_rows, normalization=args.normalization)
    metric_rows.to_csv(
        args.out_dir / "rseqc_annotation_dependence_metrics.tsv",
        sep="\t",
        index=False,
    )

    if not args.no_plot:
        plot_profiles(
            profile_rows,
            args.out_dir / "figure1b_real_rseqc_annotation_dependence.png",
            normalization=args.normalization,
        )

    print(f"Wrote RSeQC annotation-dependence results to {args.out_dir}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run real RSeQC on synthetic expressed, complete, and redundant annotations."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        help="Directory created by benchmarks/generate_synthetic_data.py.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("benchmark_runs/article_synthetic_rseqc"),
        help="Output directory for RSeQC outputs and manuscript-ready summaries.",
    )
    parser.add_argument(
        "--replicates",
        type=int,
        default=20,
        help="Replicates per synthetic read scenario when --input-dir is omitted.",
    )
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python executable used to generate synthetic input when needed.",
    )
    parser.add_argument(
        "--rseqc-output-dir",
        type=Path,
        help=(
            "Directory containing existing RSeQC output prefixes when --skip-rseqc "
            "is used. Defaults to --out-dir."
        ),
    )
    parser.add_argument(
        "--rseqc-cmd",
        default="geneBody_coverage.py",
        help="Path or command name for RSeQC geneBody_coverage.py.",
    )
    parser.add_argument(
        "--skip-rseqc",
        action="store_true",
        help="Parse existing RSeQC outputs instead of launching geneBody_coverage.py.",
    )
    parser.add_argument(
        "--image-format",
        default="png",
        help="Image format passed to RSeQC through -f.",
    )
    parser.add_argument(
        "--normalization",
        choices=["max", "mean"],
        default="max",
        help="Normalization used for the manuscript-ready profile and plot.",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip PNG plot generation.",
    )
    return parser.parse_args()


def configure_plot_cache(out_dir: Path) -> None:
    cache_dir = out_dir / ".plot-cache"
    matplotlib_cache = cache_dir / "matplotlib"
    xdg_cache = cache_dir / "xdg"
    matplotlib_cache.mkdir(parents=True, exist_ok=True)
    xdg_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_cache))
    os.environ.setdefault("XDG_CACHE_HOME", str(xdg_cache))


def generate_synthetic_input(input_dir: Path, replicates: int, python: str) -> None:
    command = [
        python,
        str(REPO_ROOT / "benchmarks" / "generate_synthetic_data.py"),
        "--out-dir",
        str(input_dir),
        "--replicates",
        str(replicates),
    ]
    subprocess.run(command, cwd=REPO_ROOT, check=True)


def resolve_command(command: str) -> str:
    if "/" in command:
        path = Path(command).expanduser()
        if not path.exists():
            raise SystemExit(f"RSeQC command does not exist: {path}")
        return str(path)
    resolved = shutil.which(command)
    if resolved is None:
        raise SystemExit(
            f"Cannot find {command!r} on PATH. Re-run with "
            "`--rseqc-cmd /path/to/geneBody_coverage.py`, or use --skip-rseqc "
            "to parse existing outputs."
        )
    return resolved


def find_existing_rseqc_prefix(label: str, output_dir: Path) -> Path:
    candidates = [
        output_dir / RSEQC_LABELS[label],
        output_dir / LEGACY_RSEQC_LABELS[label],
    ]
    for prefix in candidates:
        if Path(f"{prefix}.geneBodyCoverage.txt").exists():
            return prefix
    return candidates[0]


def run_rseqc(
    *,
    rseqc_cmd: str,
    bam_path: Path,
    annotation_path: Path,
    output_prefix: Path,
    image_format: str,
) -> None:
    command = [
        rseqc_cmd,
        "-i",
        str(bam_path),
        "-r",
        str(annotation_path),
        "-o",
        str(output_prefix),
        "-f",
        image_format,
    ]
    subprocess.run(command, check=True)


def load_expected_profiles(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def build_profile_rows(
    *,
    expected_profiles: pd.DataFrame,
    rseqc_prefixes: dict[str, Path],
    normalization: str,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    if not expected_profiles.empty:
        expected_map = {
            "truth_read_centric": "Truth read-centric",
            "rseqc_style_expressed_only_annotation": "Expected RSeQC-style expressed-only",
            "rseqc_style_complete_annotation": "Expected RSeQC-style complete",
            "rseqc_style_redundant_annotation": "Expected RSeQC-style redundant",
        }
        for profile, label in expected_map.items():
            column = f"{profile}_{normalization}_normalized"
            raw_column = f"{profile}_raw"
            if column not in expected_profiles:
                continue
            for row in expected_profiles.itertuples(index=False):
                rows.append(
                    {
                        "bin": int(row.bin),
                        "profile": profile,
                        "label": label,
                        "source": "expected",
                        "raw_coverage": getattr(row, raw_column),
                        f"{normalization}_normalized_coverage": getattr(row, column),
                        "normalization": normalization,
                    }
                )

    for profile, prefix in rseqc_prefixes.items():
        txt_path = Path(f"{prefix}.geneBodyCoverage.txt")
        if not txt_path.exists():
            raise SystemExit(
                f"Missing RSeQC output: {txt_path}. Run without --skip-rseqc "
                "or check the output prefix."
            )
        raw_values = parse_rseqc_gene_body_txt(txt_path)
        normalized_values = normalize(raw_values, normalization)
        for index, (raw_value, normalized_value) in enumerate(
            zip(raw_values, normalized_values),
            start=1,
        ):
            rows.append(
                {
                    "bin": index,
                    "profile": f"real_rseqc_{profile}",
                    "label": f"Real RSeQC {profile.replace('_', ' ')}",
                    "source": "real_RSeQC",
                    "raw_coverage": raw_value,
                    f"{normalization}_normalized_coverage": normalized_value,
                    "normalization": normalization,
                }
            )
    return pd.DataFrame(rows)


def parse_rseqc_gene_body_txt(path: Path) -> list[float]:
    table = pd.read_csv(path, sep="\t")
    if table.empty:
        raise ValueError(f"RSeQC output has no rows: {path}")
    row = table.iloc[0]
    values = []
    for column in table.columns[1:]:
        values.append(float(row[column]))
    return values


def normalize(values: list[float], method: str) -> list[float]:
    if method == "max":
        denominator = max(values) if values else 0.0
    else:
        denominator = sum(values) / len(values) if values else 0.0
    if denominator <= 0:
        return [0.0 for _ in values]
    return [value / denominator for value in values]


def summarize_profile_rows(profile_rows: pd.DataFrame, *, normalization: str) -> pd.DataFrame:
    y_column = f"{normalization}_normalized_coverage"
    truth = profile_rows.loc[profile_rows["profile"] == "truth_read_centric"]
    if truth.empty:
        return pd.DataFrame()
    truth_values = truth.sort_values("bin")[y_column].tolist()
    rows = []
    for profile, subset in profile_rows.groupby("profile", sort=False):
        values = subset.sort_values("bin")[y_column].tolist()
        rows.append(
            {
                "profile": profile,
                "source": subset["source"].iloc[0],
                "normalization": normalization,
                "mae_vs_truth": mean_absolute_error(values, truth_values),
                "rmse_vs_truth": root_mean_squared_error(values, truth_values),
                "pearson_vs_truth": pearson_correlation(values, truth_values),
                "first_decile": decile_mean(values, 0),
                "last_decile": decile_mean(values, 9),
                "first_to_last_decile_ratio": safe_divide(
                    decile_mean(values, 0),
                    decile_mean(values, 9),
                ),
            }
        )
    return pd.DataFrame(rows)


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
    return safe_divide(numerator, math.sqrt(value_var * truth_var))


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


def plot_profiles(profile_rows: pd.DataFrame, path: Path, *, normalization: str) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipping RSeQC annotation plot")
        return
    if profile_rows.empty:
        return

    y_column = f"{normalization}_normalized_coverage"
    keep_profiles = [
        "truth_read_centric",
        "real_rseqc_expressed_only_annotation",
        "real_rseqc_complete_annotation",
        "real_rseqc_redundant_annotation",
    ]
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for profile in keep_profiles:
        subset = profile_rows.loc[profile_rows["profile"] == profile].sort_values("bin")
        if subset.empty:
            continue
        ax.plot(
            subset["bin"],
            subset[y_column],
            linewidth=2,
            label=subset["label"].iloc[0],
        )
    ax.set_xlabel("Transcript body bin, 5' to 3'")
    ax.set_ylabel(f"{normalization.capitalize()}-normalized coverage")
    ax.set_title("Synthetic real RSeQC annotation-dependence test")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
