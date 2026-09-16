#!/usr/bin/env python3
"""Select the CUDA device with the most free memory for plasticity runners."""

from __future__ import annotations

import argparse
import sys


def _pick_least_loaded_cuda_device(*, reason: str) -> int:
    """Return the physical CUDA index with the most free device memory.

    The selected index is intended for ``CUDA_VISIBLE_DEVICES``.  After the
    launcher exports that variable, CUDA libraries inside the MOOSE process see
    the selected physical GPU as logical ``cuda:0``.  Selecting once before the
    staged run keeps every benchmark repetition on the same device.
    """
    try:
        import torch
    except ImportError as exc:
        raise RuntimeError(
            "CUDA auto-selection requires the NEML2 Python environment with torch installed"
        ) from exc

    if not torch.cuda.is_available():
        raise RuntimeError("CUDA auto-selection requested but CUDA is unavailable")

    n_devices = torch.cuda.device_count()
    if n_devices == 0:
        raise RuntimeError("CUDA is available but no CUDA devices are visible")

    free_by_index: list[tuple[int, int]] = []
    for index in range(n_devices):
        free_bytes, total_bytes = torch.cuda.mem_get_info(index)
        free_bytes = int(free_bytes)
        total_bytes = int(total_bytes)
        free_by_index.append((free_bytes, index))
        print(
            f"[cuda] physical cuda:{index}: "
            f"{free_bytes / 1024**3:.1f} / {total_bytes / 1024**3:.1f} GiB free",
            file=sys.stderr,
        )

    # Most free memory wins; ties go to the lower physical device index.
    free_bytes, chosen = max(free_by_index, key=lambda pair: (pair[0], -pair[1]))
    print(
        f"[cuda] {reason}: auto-selected physical cuda:{chosen} "
        f"({free_bytes / 1024**3:.1f} GiB free); the MOOSE process will see it "
        f"as logical cuda:0 via CUDA_VISIBLE_DEVICES={chosen}",
        file=sys.stderr,
    )
    return chosen


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Print the physical CUDA index with the most free memory."
    )
    parser.add_argument(
        "--reason",
        default="MOOSE plasticity run",
        help="label included in the selection diagnostic",
    )
    args = parser.parse_args()

    # Keep stdout machine-readable for shell command substitution.  All human
    # diagnostics from the selector are written to stderr.
    print(_pick_least_loaded_cuda_device(reason=args.reason))


if __name__ == "__main__":
    main()
