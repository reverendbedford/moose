#!/bin/bash
# One-run N=32 check before the formal benchmark.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

RESULTS_DIR="$SCRIPT_DIR/results_n32_smoke" \
MESH_N=32 \
REPEATS=1 \
PROFILE=0 \
  "$SCRIPT_DIR/run_benchmarks.sh"
