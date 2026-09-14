#!/bin/bash
# One-run N=64 memory, runtime, and convergence check.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

RESULTS_DIR="$SCRIPT_DIR/results_n64_smoke" \
MESH_N=64 \
REPEATS=1 \
PROFILE=0 \
  "$SCRIPT_DIR/run_benchmarks.sh"
