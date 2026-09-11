#!/bin/bash
# Step 4: verify MOOSE actually dispatches CUDA kernels, then run
# benchmark_cpu_and_gpu.py without --skip-gpu.
#
# nvidia-smi is too coarse for a MOOSE smoke run — a small mesh finishes in
# under a second, well below nvidia-smi's polling interval, so it sees no
# process. `nsys profile --stats=true` reports per-kernel counts even for
# millisecond-lived runs, and any `Moose::Kokkos::` symbol in the output is
# conclusive proof that MOOSE launched CUDA kernels.
set -e

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
. "$SCRIPT_DIR/env.sh"

export PETSC_DIR="$PREFIX"
export PETSC_ARCH=""
export LIBMESH_DIR="$PREFIX"

BENCH="$MOOSE_DIR/modules/solid_mechanics/test/tests/kokkos/linear_elasticity"
LOG="$LOGS/benchmark-$(date +%Y%m%d-%H%M%S).log"
NSYS_OUT="$LOGS/moose-smoke-$(date +%Y%m%d-%H%M%S)"

echo "[run_benchmark] nsys-profiling a smoke run to verify MOOSE launches CUDA kernels..."
/usr/local/cuda/bin/nsys profile --stats=true -o "$NSYS_OUT" \
  mpirun -np 1 "$MOOSE_DIR/modules/solid_mechanics/solid_mechanics-opt" \
    -i "$BENCH/kokkos_material_linear_elasticity.i" \
    --compute-device=cuda --n-threads=1 \
    Mesh/generated/nx=64 Mesh/generated/ny=64 Outputs/exodus=false \
    > "$LOG.smoke" 2>&1

MOOSE_KERNELS=$(grep -c "Moose::Kokkos::" "$LOG.smoke" || true)
echo "[run_benchmark] MOOSE Kokkos kernels launched on GPU in smoke run: $MOOSE_KERNELS"
if [ "$MOOSE_KERNELS" -lt 1 ]; then
  echo "[run_benchmark] ERROR: nsys did not detect any Moose::Kokkos:: CUDA kernels."
  echo "                Smoke log at $LOG.smoke; nsys report at $NSYS_OUT.nsys-rep"
  exit 2
fi

echo "[run_benchmark] Kicking off full benchmark..."
python3 -u "$BENCH/benchmark_cpu_and_gpu.py" 2>&1 | tee "$LOG"
