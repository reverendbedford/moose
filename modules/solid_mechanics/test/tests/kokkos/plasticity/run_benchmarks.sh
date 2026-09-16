#!/bin/bash
# Run timing repetitions separately from one Nsight Systems profile per step.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MOOSE_DIR=$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)

# Source the CUDA-Kokkos stack activation so OpenMPI can find its runtime
# data files. The stack was relocated after build; without OPAL_PREFIX the
# opal_wrapper falls back to a CWD-relative lookup and MPI_Init fails.
# shellcheck disable=SC1091
source "$MOOSE_DIR/kokkos-cuda-stack/scripts/activate.sh"

# Pin the entire staged benchmark to one physical GPU.  If the caller already
# selected CUDA devices (for example through a scheduler or an explicit shell
# export), preserve that choice.  Otherwise choose the GPU with the most free
# memory once, then let CUDA remap it to logical cuda:0 inside every child.
if [ -z "${CUDA_VISIBLE_DEVICES:-}" ]; then
  CUDA_VISIBLE_DEVICES=$(python3 "$SCRIPT_DIR/select_cuda_device.py" \
    --reason "plasticity benchmark")
  export CUDA_VISIBLE_DEVICES
else
  echo "[cuda] honoring existing CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES" >&2
fi

EXE=${EXE:-$MOOSE_DIR/modules/solid_mechanics/solid_mechanics-opt}
RESULTS_DIR=${RESULTS_DIR:-$SCRIPT_DIR/results}
MESH_N=${MESH_N:-16}
REPEATS=${REPEATS:-3}
PROFILE=${PROFILE:-1}

if [ ! -x "$EXE" ]; then
  echo "ERROR: solid_mechanics executable not found at $EXE" >&2
  echo "Set EXE=/path/to/solid_mechanics-opt and rerun." >&2
  exit 1
fi

mkdir -p "$RESULTS_DIR"

CPU_PETSC_ARGS=(
  -vec_type standard
  -nl0_mat_type aij
  -use_gpu_aware_mpi 0
  -snes_converged_reason
  -ksp_converged_reason
  -ksp_view
)
GPU_PETSC_ARGS=(
  -vec_type kokkos
  -nl0_mat_type aijkokkos
  -use_gpu_aware_mpi 0
  -snes_converged_reason
  -ksp_converged_reason
  -ksp_view
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
  local petsc_args=("${CPU_PETSC_ARGS[@]}")
  local device_args=()

  if [[ "$step" == step3a_* || "$step" == step3b_* || "$step" == step4_* ]]; then
    device_args=(--compute-device=cuda)
  fi
  if [[ "$step" == step4_* ]]; then
    petsc_args=("${GPU_PETSC_ARGS[@]}")
  fi

  for rep in $(seq 1 "$REPEATS"); do
    local prefix="$RESULTS_DIR/${step}_rep${rep}"
    echo "Running $step repetition $rep/$REPEATS"
    /usr/bin/time -f '%e' -o "$prefix.time" \
      "$EXE" -i "$input" "${device_args[@]}" \
        "N=$MESH_N" \
        "Outputs/perf_graph_json_file=$prefix.perf.json" \
        "${petsc_args[@]}" \
        --timing \
        > "$prefix.log" 2>&1
  done

  if [ "$PROFILE" = 1 ]; then
    if ! command -v nsys >/dev/null 2>&1; then
      echo "WARNING: nsys not found; skipping the $step GPU/MPI profile." >&2
      return
    fi

    local profile_prefix="$RESULTS_DIR/${step}_profile"
    echo "Profiling $step once with Nsight Systems"
    nsys profile --force-overwrite=true --trace=cuda,nvtx,mpi \
      -o "$profile_prefix" \
      "$EXE" -i "$input" "${device_args[@]}" \
        "N=$MESH_N" \
        "Outputs/perf_graph_json_file=$profile_prefix.perf.json" \
        "${petsc_args[@]}" \
        --timing \
        > "$profile_prefix.log" 2>&1
  fi
}

for step in "${steps[@]}"; do
  run_step "$step"
done

python3 "$SCRIPT_DIR/analyze_results.py" "$RESULTS_DIR"
