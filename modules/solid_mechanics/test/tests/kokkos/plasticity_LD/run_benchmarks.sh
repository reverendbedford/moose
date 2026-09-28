#!/bin/bash
# Large-deformation counterpart of plasticity/run_benchmarks.sh.
# Same shape as the SD script: pick a GPU once, then run the five staged
# inputs REPEATS times each with GNU time and one Nsight Systems profile.
#
# The step names, PETSc argument sets, --compute-device selection, and the
# GNU time / perf-graph / nsys instrumentation are copied verbatim from
# plasticity/run_benchmarks.sh so both suites can be compared with the same
# tooling.  Only the SCRIPT_DIR-derived RESULTS_DIR (plasticity_LD/results by
# default) and the input paths (this directory) change.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MOOSE_DIR=$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)

# shellcheck disable=SC1091
source "$MOOSE_DIR/kokkos-cuda-stack/scripts/activate.sh"

# Reuse the SD select_cuda_device.py helper.
if [ -z "${CUDA_VISIBLE_DEVICES:-}" ]; then
  CUDA_VISIBLE_DEVICES=$(python3 "$MOOSE_DIR/modules/solid_mechanics/test/tests/kokkos/plasticity/select_cuda_device.py" \
    --reason "plasticity_LD benchmark")
  export CUDA_VISIBLE_DEVICES
else
  echo "[cuda] honoring existing CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES" >&2
fi

EXE=${EXE:-$MOOSE_DIR/modules/solid_mechanics/solid_mechanics-opt}
RESULTS_DIR=${RESULTS_DIR:-$SCRIPT_DIR/results}
MESH_N=${MESH_N:-16}
REPEATS=${REPEATS:-3}
PROFILE=${PROFILE:-1}
# Extra Executioner overrides let callers exercise finite kinematics without
# forking the script (for example NUM_STEPS=10 DT=5e-3 -> 5% engineering strain).
NUM_STEPS=${NUM_STEPS:-5}
DT=${DT:-1e-3}

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
  step3_plasticity_cpu_neml2_kokkos_cpu_petsc
  step4_plasticity_gpu_neml2_kokkos_cpu_petsc
  step5_plasticity_full_gpu
)

run_step()
{
  local step=$1
  local input="$SCRIPT_DIR/$step.i"
  local petsc_args=("${CPU_PETSC_ARGS[@]}")
  local device_args=()

  if [[ "$step" == step3_* || "$step" == step4_* || "$step" == step5_* ]]; then
    device_args=(--compute-device=cuda)
  fi
  if [[ "$step" == step5_* ]]; then
    petsc_args=("${GPU_PETSC_ARGS[@]}")
  fi

  for rep in $(seq 1 "$REPEATS"); do
    local prefix="$RESULTS_DIR/${step}_rep${rep}"
    echo "Running $step repetition $rep/$REPEATS"
    /usr/bin/time -f '%e' -o "$prefix.time" \
      "$EXE" -i "$input" "${device_args[@]}" \
        "N=$MESH_N" \
        "Executioner/dt=$DT" \
        "Executioner/num_steps=$NUM_STEPS" \
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
        "Executioner/dt=$DT" \
        "Executioner/num_steps=$NUM_STEPS" \
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
