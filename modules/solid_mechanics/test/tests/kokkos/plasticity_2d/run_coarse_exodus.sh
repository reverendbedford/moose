#!/bin/bash
# Compare the five 2D execution paths while retaining non-preset Kokkos loading.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MOOSE_DIR=$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)

# shellcheck disable=SC1091
source "$MOOSE_DIR/kokkos-cuda-stack/scripts/activate.sh"

EXE=${EXE:-$MOOSE_DIR/modules/solid_mechanics/solid_mechanics-opt}
MESH_N=${MESH_N:-8}
OUTPUT_DIR=${OUTPUT_DIR:-$SCRIPT_DIR/coarse_exodus_n$MESH_N}

if [ ! -x "$EXE" ]; then
  echo "ERROR: solid_mechanics executable not found at $EXE" >&2
  echo "Set EXE=/path/to/solid_mechanics-opt and rerun." >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

CPU_PETSC_ARGS=(
  -vec_type standard
  -nl0_mat_type aij
  -use_gpu_aware_mpi 0
  -snes_converged_reason
  -ksp_converged_reason
)
GPU_PETSC_ARGS=(
  -vec_type kokkos
  -nl0_mat_type aijkokkos
  -use_gpu_aware_mpi 0
  -snes_converged_reason
  -ksp_converged_reason
)

steps=(
  step1_plasticity_cpu_neml2
  step2_plasticity_gpu_neml2
  step3a_plasticity_cpu_neml2_kokkos_cpu_petsc
  step3b_plasticity_gpu_neml2_kokkos_cpu_petsc
  step4_plasticity_full_gpu
)

run_step()
{
  local step=$1
  local input="$SCRIPT_DIR/$step.i"
  local prefix="$OUTPUT_DIR/$step"
  local petsc_args=("${CPU_PETSC_ARGS[@]}")
  local command_args=(-i "$input")

  if [[ "$step" == step3a_* || "$step" == step3b_* || "$step" == step4_* ]]; then
    command_args+=(--compute-device=cuda)
  fi
  if [[ "$step" == step4_* ]]; then
    petsc_args=("${GPU_PETSC_ARGS[@]}")
  fi
  command_args+=(
    "N=$MESH_N"
    "Outputs/exodus=true"
    "Outputs/csv=true"
    "Outputs/file_base=$prefix"
    "${petsc_args[@]}"
  )

  echo
  echo "Running $step (N=$MESH_N)"
  if ! "$EXE" "${command_args[@]}" > "$prefix.log" 2>&1; then
    echo "ERROR: $step failed; last 40 log lines:" >&2
    tail -n 40 "$prefix.log" >&2
    return 1
  fi

  echo "  log: $prefix.log"
  if [ -s "$prefix.csv" ]; then
    echo "  final CSV row:"
    tail -n 1 "$prefix.csv" | sed 's/^/    /'
  fi
}

for step in "${steps[@]}"; do
  run_step "$step"
done

echo
echo "2D coarse solutions written to: $OUTPUT_DIR"
echo "At t=0.005, expect ux_right ~= 0.005 and ux_center ~= 0.0025."
