#!/usr/bin/env python3
"""Plot benchmark results from an existing comparison.csv produced by analyze_results.py."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def plot_results(frame: pd.DataFrame, output: Path) -> None:
    labels = frame["label"]
    figure, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)

    axes[0, 0].bar(labels, frame["total_wall_time_s"])
    axes[0, 0].set_title("Total wall time")
    axes[0, 0].set_ylabel("Seconds")

    axes[0, 1].bar(labels, frame["neml2_solve_s"])
    axes[0, 1].set_title("NEML2::solve")
    axes[0, 1].set_ylabel("Seconds")

    copy_bottom = frame["cuda_kernel_time_s"].fillna(0)
    axes[1, 0].bar(labels, copy_bottom, label="CUDA kernels")
    axes[1, 0].bar(
        labels,
        frame["h2d_time_s"].fillna(0),
        bottom=copy_bottom,
        label="H2D",
    )
    axes[1, 0].bar(
        labels,
        frame["d2h_time_s"].fillna(0),
        bottom=copy_bottom + frame["h2d_time_s"].fillna(0),
        label="D2H",
    )
    axes[1, 0].set_title("GPU activity from one profile run")
    axes[1, 0].set_ylabel("Seconds")
    axes[1, 0].legend()

    axes[1, 1].plot(labels, frame["newton_iterations"], "o-", label="Newton")
    axes[1, 1].plot(labels, frame["ksp_iterations"], "s-", label="KSP")
    axes[1, 1].set_title("Iteration totals")
    axes[1, 1].set_yscale('log')
    axes[1, 1].set_ylabel("Iterations")
    axes[1, 1].legend()

    for axis in axes.flat:
        axis.tick_params(axis="x", rotation=18)
        axis.grid(axis="y", alpha=0.25)
    figure.savefig(output, dpi=180)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot benchmark results from an existing comparison.csv."
    )
    parser.add_argument(
        "csv_path",
        type=Path,
        nargs="?",
        default=Path("results/comparison.csv"),
        help="Path to comparison.csv (default: results/comparison.csv)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output PNG path (default: <csv_path parent>/comparison.png)",
    )
    args = parser.parse_args()

    frame = pd.read_csv(args.csv_path)
    output = args.output if args.output is not None else args.csv_path.with_name("comparison.png")
    plot_results(frame, output)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
