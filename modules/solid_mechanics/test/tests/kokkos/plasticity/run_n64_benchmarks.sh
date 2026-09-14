#!/bin/bash
# N=64 benchmark after the smoke run succeeds.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

RESULTS_DIR="$SCRIPT_DIR/results_n64" \
MESH_N=64 \
REPEATS=3 \
PROFILE=0 \
  "$SCRIPT_DIR/run_benchmarks.sh"
