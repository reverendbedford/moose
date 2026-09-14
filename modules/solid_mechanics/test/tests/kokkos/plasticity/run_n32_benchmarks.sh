#!/bin/bash
# Recommended N=32 benchmark: five clean timings plus one profile per step.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

RESULTS_DIR="$SCRIPT_DIR/results_n32" \
MESH_N=32 \
REPEATS=5 \
PROFILE=1 \
  "$SCRIPT_DIR/run_benchmarks.sh"
