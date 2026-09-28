#!/bin/bash
# One-run N=32 check before the formal LD benchmark.
# Mirror of plasticity/run_n32_smoke.sh. run_benchmarks.sh calls
# analyze_results.py at the end, which writes comparison.csv and
# comparison.png (per-step wall time, NEML2::solve, GPU activity, and
# Newton/KSP iteration totals) into $RESULTS_DIR.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

RESULTS_DIR="$SCRIPT_DIR/results_n32_smoke" \
MESH_N=32 \
REPEATS=1 \
PROFILE=0 \
  "$SCRIPT_DIR/run_benchmarks.sh"
