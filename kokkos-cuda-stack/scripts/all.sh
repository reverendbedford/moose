#!/bin/bash
set -e
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
"$SCRIPT_DIR/build_openmpi.sh"
"$SCRIPT_DIR/build_petsc.sh"
"$SCRIPT_DIR/build_libmesh.sh"
"$SCRIPT_DIR/build_wasp.sh"
"$SCRIPT_DIR/build_moose.sh"
"$SCRIPT_DIR/run_benchmark.sh"
