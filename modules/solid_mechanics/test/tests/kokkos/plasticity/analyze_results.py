#!/usr/bin/env python3
"""Summarize and plot the five NEML2/Kokkos/PETSc benchmark steps.

This version is robust to GNU time files that contain messages such as

    Command exited with non-zero status 1
    95.68

Failed repetitions are counted but excluded from timing averages.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


STEPS = [
    ("step1_plasticity_cpu_neml2", "1. CPU NEML2"),
    ("step2_plasticity_gpu_neml2", "2. GPU NEML2"),
    (
        "step3a_plasticity_cpu_neml2_kokkos_cpu_petsc",
        "3. CPU NEML2 + GPU assembly",
    ),
    (
        "step3b_plasticity_gpu_neml2_kokkos_cpu_petsc",
        "4. GPU NEML2 + GPU assembly",
    ),
    ("step4_plasticity_full_gpu", "5. Full GPU"),
]


FLOAT_LINE = re.compile(
    r"^\s*[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?\s*$"
)


def read_wall_time(path: Path) -> tuple[float, bool]:
    """Read wall time from a GNU time output file.

    Returns
    -------
    wall_time
        The last standalone floating-point number in the file.
    failed
        True if GNU time reported that the command exited nonzero.
    """
    if not path.exists():
        return math.nan, True

    text = path.read_text(errors="replace")

    failed = "Command exited with non-zero status" in text

    values = [
        float(line.strip())
        for line in text.splitlines()
        if FLOAT_LINE.match(line)
    ]

    if not values:
        return math.nan, True

    return values[-1], failed


def perf_section_time(path: Path, wanted: str) -> float:
    """Sum self time for all nodes whose PerfGraph section name matches wanted."""
    if not path.exists():
        return math.nan

    try:
        data = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        return math.nan

    try:
        graph = data["time_steps"][-1]["perf_graph_json"]["graph"]
    except (KeyError, IndexError, TypeError):
        return math.nan

    def visit(value: object) -> float:
        if not isinstance(value, dict):
            return 0.0

        total = 0.0

        for name, node in value.items():
            if not isinstance(node, dict):
                continue

            if name == wanted:
                try:
                    total += float(node.get("time", 0.0))
                except (TypeError, ValueError):
                    pass

            total += visit(node.get("children", {}))

        return total

    return visit(graph)


def iteration_total(path: Path, kind: str) -> float:
    """Sum nonlinear or linear iteration counts from a MOOSE log."""
    if not path.exists():
        return math.nan

    text = path.read_text(errors="replace")

    pattern = re.compile(
        rf"{kind} solve converged.*?iterations?\s+(\d+)",
        re.IGNORECASE,
    )

    values = [int(value) for value in pattern.findall(text)]

    return float(sum(values)) if values else math.nan


def nsys_csv(report_path: Path, report_name: str) -> list[dict[str, str]]:
    """Read one Nsight Systems CSV report from an .nsys-rep file."""
    if not report_path.exists() or shutil.which("nsys") is None:
        return []

    result = subprocess.run(
        [
            "nsys",
            "stats",
            "--report",
            report_name,
            "--format",
            "csv",
            str(report_path),
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    if result.returncode:
        return []

    lines = result.stdout.splitlines()

    header = next(
        (index for index, line in enumerate(lines) if "Total Time" in line),
        None,
    )

    if header is None:
        return []

    return list(csv.DictReader(io.StringIO("\n".join(lines[header:]))))


def seconds(value: str, header: str) -> float:
    """Convert an Nsight Systems time column to seconds."""
    units = {
        "ns": 1e-9,
        "us": 1e-6,
        "ms": 1e-3,
        "s": 1.0,
    }

    match = re.search(r"\((ns|us|ms|s)\)", header)
    scale = units[match.group(1)] if match else 1e-9

    return float(value.replace(",", "")) * scale


def report_time(
    rows: list[dict[str, str]],
    name_pattern: str | None = None,
) -> float:
    """Sum one Nsight Systems report time column, optionally filtered by name."""
    total = 0.0
    found = False

    for row in rows:
        time_key = next((key for key in row if "Total Time" in key), None)

        if time_key is None:
            continue

        if name_pattern is not None:
            name = " ".join(str(value) for value in row.values())

            if re.search(name_pattern, name, re.IGNORECASE) is None:
                continue

        try:
            total += seconds(row[time_key], time_key)
            found = True
        except (TypeError, ValueError):
            pass

    return total if found else math.nan


def mean(values: list[float]) -> float:
    """Mean of all non-NaN values."""
    valid = [value for value in values if not math.isnan(value)]

    return sum(valid) / len(valid) if valid else math.nan


def collect(results_dir: Path) -> pd.DataFrame:
    """Collect benchmark data from one results directory."""
    records = []

    for step, label in STEPS:
        wall_times = []
        neml2_times = []
        newton_iterations = []
        ksp_iterations = []

        total_reps = 0
        failed_reps = 0
        used_reps = 0

        for time_path in sorted(results_dir.glob(f"{step}_rep*.time")):
            total_reps += 1

            stem = time_path.with_suffix("")

            wall_time, failed = read_wall_time(time_path)

            if failed:
                failed_reps += 1
                continue

            if math.isnan(wall_time):
                failed_reps += 1
                continue

            used_reps += 1

            wall_times.append(wall_time)

            neml2_times.append(
                perf_section_time(
                    Path(f"{stem}.perf.json"),
                    "NEML2::solve",
                )
            )

            log_path = Path(f"{stem}.log")

            newton_iterations.append(
                iteration_total(log_path, "Nonlinear")
            )

            ksp_iterations.append(
                iteration_total(log_path, "Linear")
            )

        profile = results_dir / f"{step}_profile.nsys-rep"

        kernels = nsys_csv(profile, "cuda_gpu_kern_sum")
        copies = nsys_csv(profile, "cuda_gpu_mem_time_sum")
        mpi = nsys_csv(profile, "mpi_event_sum")

        records.append(
            {
                "step": step,
                "label": label,
                "total_reps": total_reps,
                "used_reps": used_reps,
                "failed_reps": failed_reps,
                "total_wall_time_s": mean(wall_times),
                "neml2_solve_s": mean(neml2_times),
                "cuda_kernel_time_s": report_time(kernels),
                "mpi_time_s": report_time(mpi),
                "h2d_time_s": report_time(copies, r"HtoD|Host to Device"),
                "d2h_time_s": report_time(copies, r"DtoH|Device to Host"),
                "newton_iterations": mean(newton_iterations),
                "ksp_iterations": mean(ksp_iterations),
            }
        )

    return pd.DataFrame.from_records(records)


def plot_results(frame: pd.DataFrame, output: Path) -> None:
    """Generate summary plots."""
    labels = frame["label"]

    figure, axes = plt.subplots(
        2,
        2,
        figsize=(13, 8),
        constrained_layout=True,
    )

    axes[0, 0].bar(labels, frame["total_wall_time_s"])
    axes[0, 0].set_title("Total wall time")
    axes[0, 0].set_ylabel("Seconds")

    axes[0, 1].bar(labels, frame["neml2_solve_s"])
    axes[0, 1].set_title("NEML2::solve")
    axes[0, 1].set_ylabel("Seconds")

    copy_bottom = frame["cuda_kernel_time_s"].fillna(0.0)

    axes[1, 0].bar(
        labels,
        copy_bottom,
        label="CUDA kernels",
    )

    axes[1, 0].bar(
        labels,
        frame["h2d_time_s"].fillna(0.0),
        bottom=copy_bottom,
        label="H2D",
    )

    axes[1, 0].bar(
        labels,
        frame["d2h_time_s"].fillna(0.0),
        bottom=copy_bottom + frame["h2d_time_s"].fillna(0.0),
        label="D2H",
    )

    axes[1, 0].set_title("GPU activity from one profile run")
    axes[1, 0].set_ylabel("Seconds")
    axes[1, 0].legend()

    axes[1, 1].plot(
        labels,
        frame["newton_iterations"],
        "o-",
        label="Newton",
    )

    axes[1, 1].plot(
        labels,
        frame["ksp_iterations"],
        "s-",
        label="KSP",
    )

    axes[1, 1].set_title("Iteration totals")
    axes[1, 1].set_ylabel("Iterations")
    axes[1, 1].set_yscale('log')
    axes[1, 1].legend()

    for axis in axes.flat:
        axis.tick_params(axis="x", rotation=18)
        axis.grid(axis="y", alpha=0.25)

    figure.savefig(output, dpi=180)
    plt.close(figure)


def print_warnings(frame: pd.DataFrame) -> None:
    """Print useful warnings for failed or missing benchmark data."""
    failed = frame[frame["failed_reps"] > 0]

    if not failed.empty:
        print("\nWarning: some repetitions failed and were excluded from averages.")

        for _, row in failed.iterrows():
            print(
                f"  {row['step']}: "
                f"{int(row['failed_reps'])} failed, "
                f"{int(row['used_reps'])} used, "
                f"{int(row['total_reps'])} total"
            )

    missing = frame[frame["used_reps"] == 0]

    if not missing.empty:
        print("\nWarning: some steps have no successful repetitions.")

        for _, row in missing.iterrows():
            print(f"  {row['step']}")


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "results_dir",
        type=Path,
        nargs="?",
        default=Path("results"),
        help="Directory containing benchmark result files.",
    )

    args = parser.parse_args()

    args.results_dir.mkdir(parents=True, exist_ok=True)

    frame = collect(args.results_dir)

    csv_path = args.results_dir / "comparison.csv"
    plot_path = args.results_dir / "comparison.png"

    frame.to_csv(csv_path, index=False)
    plot_results(frame, plot_path)

    print(frame.to_string(index=False))
    print_warnings(frame)

    print(f"\nWrote {csv_path}")
    print(f"Wrote {plot_path}")
    print("MPI time is not meaningful for these single-rank runs.")


if __name__ == "__main__":
    main()
