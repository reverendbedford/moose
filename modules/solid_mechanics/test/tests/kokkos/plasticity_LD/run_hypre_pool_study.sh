#!/bin/bash
# ==============================================================================
# run_hypre_pool_study.sh
#
# Benchmark study for HYPRE GPU Umpire device-pool memory usage and OOM analysis
# on the Step-5 full-GPU N=32 plasticity benchmark.
#
# Evaluates:
#   - default (HYPRE default 4-GiB pool, control case)
#   - 256 MiB (-hypre_umpire_device_pool_size 256)
#   - 512 MiB (-hypre_umpire_device_pool_size 512)
#   - 1024 MiB (-hypre_umpire_device_pool_size 1024)
#
# Keeps all application code, F-bar settings, NEML2 settings, Kokkos settings,
# and PETSc solver options identical to Step 5.
# ==============================================================================
set -u

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MOOSE_DIR=$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)

# Reuse existing benchmark environment activation
# shellcheck disable=SC1091
source "$MOOSE_DIR/kokkos-cuda-stack/scripts/activate.sh"

# Honor or auto-select CUDA device
if [ -z "${CUDA_VISIBLE_DEVICES:-}" ]; then
  CUDA_VISIBLE_DEVICES=$(python3 "$MOOSE_DIR/modules/solid_mechanics/test/tests/kokkos/plasticity/select_cuda_device.py" \
    --reason "HYPRE device-pool study")
  export CUDA_VISIBLE_DEVICES
else
  echo "[cuda] honoring existing CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES" >&2
fi

EXE=${EXE:-$MOOSE_DIR/modules/solid_mechanics/solid_mechanics-opt}
RESULTS_DIR=${RESULTS_DIR:-$SCRIPT_DIR/results_hypre_pool_study}
MESH_N=${MESH_N:-32}
NUM_STEPS=${NUM_STEPS:-5}
DT=${DT:-1e-3}
INPUT_FILE="$SCRIPT_DIR/step5_plasticity_full_gpu.i"

if [ ! -x "$EXE" ]; then
  echo "ERROR: solid_mechanics executable not found at $EXE" >&2
  exit 1
fi

if [ ! -f "$INPUT_FILE" ]; then
  echo "ERROR: Step-5 input file not found at $INPUT_FILE" >&2
  exit 1
fi

mkdir -p "$RESULTS_DIR"

# Exact base PETSc arguments from Step-5 benchmark
BASE_PETSC_ARGS=(
  -vec_type kokkos
  -nl0_mat_type aijkokkos
  -use_gpu_aware_mpi 0
  -snes_converged_reason
  -ksp_converged_reason
  -ksp_view
)

CASES=("default" "256" "512" "1024")

echo "=============================================================================="
echo "Starting HYPRE GPU Umpire Device-Pool Study"
echo "Results directory: $RESULTS_DIR"
echo "Executable:        $EXE"
echo "Input file:        $INPUT_FILE"
echo "Mesh:              N=$MESH_N (32k elements, 107k DOFs)"
echo "Steps:             $NUM_STEPS (dt=$DT)"
echo "Visible GPU:       CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-all}"
echo "Cases:             ${CASES[*]}"
echo "=============================================================================="

for case in "${CASES[@]}"; do
  echo ""
  echo ">>> [Case: $case] Starting at $(date '+%Y-%m-%d %H:%M:%S') <<<"

  # Ensure no orphan benchmark processes remain active
  pkill -9 -f "solid_mechanics-opt.*step5_plasticity_full_gpu" 2>/dev/null || true
  sleep 1

  prefix="$RESULTS_DIR/hypre_pool_${case}"
  log_file="${prefix}.log"
  time_file="${prefix}.time"
  perf_file="${prefix}.perf.json"
  gpu_mem_csv="${prefix}_gpu_memory.csv"
  gpu_before="${prefix}_before.txt"
  gpu_after="${prefix}_after.txt"

  # Record GPU state before run
  nvidia-smi > "$gpu_before" 2>&1 || true

  # Determine PETSc options for this case
  cmd_petsc_args=("${BASE_PETSC_ARGS[@]}")
  if [ "$case" != "default" ]; then
    cmd_petsc_args+=("-hypre_umpire_device_pool_size" "$case")
  fi

  # Start low-overhead background GPU memory sampling (1 Hz)
  (
    echo "timestamp,gpu_index,gpu_name,memory_used_mib,memory_total_mib,gpu_utilization_percent"
    while true; do
      nvidia-smi --query-gpu=timestamp,index,name,memory.used,memory.total,utilization.gpu --format=csv,noheader,nounits 2>/dev/null || true
      sleep 1
    done
  ) > "$gpu_mem_csv" 2>&1 &
  MONITOR_PID=$!

  # Run Step-5 case sequentially
  START_TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
  RUN_EXIT_CODE=0
  /usr/bin/time -f '%e' -o "$time_file" \
    "$EXE" -i "$INPUT_FILE" \
      --compute-device=cuda \
      "N=$MESH_N" \
      "Executioner/dt=$DT" \
      "Executioner/num_steps=$NUM_STEPS" \
      "Outputs/perf_graph_json_file=$perf_file" \
      "${cmd_petsc_args[@]}" \
      --timing \
      > "$log_file" 2>&1 || RUN_EXIT_CODE=$?
  END_TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

  # Stop GPU memory monitor
  kill "$MONITOR_PID" 2>/dev/null || true
  wait "$MONITOR_PID" 2>/dev/null || true

  # Record GPU state after run
  nvidia-smi > "$gpu_after" 2>&1 || true

  echo ">>> [Case: $case] Completed at $END_TIMESTAMP (exit code: $RUN_EXIT_CODE) <<<"
done

# Run parser and summary generator
python3 - << 'EOF' "$RESULTS_DIR" "${CASES[@]}"
import sys
import re
import csv
from pathlib import Path

results_dir = Path(sys.argv[1])
cases = sys.argv[2:]

FLOAT_LINE = re.compile(r"^\s*[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?\s*$")

results = []

for case in cases:
    prefix = results_dir / f"hypre_pool_{case}"
    log_file = prefix.with_suffix(".log")
    time_file = prefix.with_suffix(".time")
    gpu_csv = results_dir / f"hypre_pool_{case}_gpu_memory.csv"

    # Wall time & exit status from GNU time
    wall_seconds = "NA"
    time_failed = False
    if time_file.exists():
        time_text = time_file.read_text(errors="replace")
        if "Command exited with non-zero status" in time_text:
            time_failed = True
        vals = [float(l.strip()) for l in time_text.splitlines() if FLOAT_LINE.match(l)]
        if vals:
            wall_seconds = f"{vals[-1]:.2f}"

    # Log file inspection
    completed = False
    oom = False
    newton_iters = "NA"
    ksp_iters = "NA"
    exit_code = "0" if (not time_failed and wall_seconds != "NA") else "1"
    failure_msg = "NA"

    if log_file.exists():
        log_text = log_file.read_text(errors="replace")
        
        # Check completion
        if "Finished Executing" in log_text and not time_failed:
            completed = True
            exit_code = "0"

        # Check OOM / Umpire runtime_error / cudaMalloc
        oom_match = re.search(r"(!\s*Umpire runtime_error.*?out of memory|cudaMalloc.*?failed with error:\s*out of memory|out of memory)", log_text, re.IGNORECASE)
        if oom_match:
            oom = True
            failure_msg = oom_match.group(0).strip().replace("\n", " ")
        elif not completed:
            # Capture terminating or abort line
            for line in log_text.splitlines():
                if "terminating:" in line or "MPI_ABORT" in line or "Error" in line:
                    failure_msg = line.strip()
                    break

        # Iteration counts
        newton_vals = [int(v) for v in re.findall(r"Nonlinear solve converged.*?iterations?\s+(\d+)", log_text, re.IGNORECASE)]
        if newton_vals:
            newton_iters = str(sum(newton_vals))

        ksp_vals = [int(v) for v in re.findall(r"Linear solve converged.*?iterations?\s+(\d+)", log_text, re.IGNORECASE)]
        if ksp_vals:
            ksp_iters = str(sum(ksp_vals))

    # Peak GPU memory from monitored CSV
    peak_gpu_mem = "NA"
    if gpu_csv.exists():
        try:
            mem_vals = []
            with open(gpu_csv, "r", encoding="utf-8", errors="replace") as f:
                reader = csv.reader(f)
                header = next(reader, None)
                for row in reader:
                    if len(row) >= 4:
                        try:
                            mem_vals.append(float(row[3].strip()))
                        except ValueError:
                            pass
            if mem_vals:
                peak_gpu_mem = f"{max(mem_vals):.0f}"
        except Exception:
            pass

    pool_mib = "4096 (default)" if case == "default" else case
    status = "PASS" if completed else ("OOM" if oom else "FAIL")

    results.append({
        "case": case,
        "pool_mib": pool_mib,
        "status": status,
        "exit_code": exit_code,
        "completed": "true" if completed else "false",
        "oom": "true" if oom else "false",
        "newton_iterations": newton_iters,
        "ksp_iterations": ksp_iters,
        "wall_seconds": wall_seconds,
        "peak_gpu_memory_mib": peak_gpu_mem,
        "failure_msg": failure_msg
    })

# Write CSV summary
csv_file = results_dir / "hypre_pool_study.csv"
with open(csv_file, "w", newline="", encoding="utf-8") as f:
    fieldnames = ["case", "pool_mib", "status", "exit_code", "completed", "oom", "newton_iterations", "ksp_iterations", "wall_seconds", "peak_gpu_memory_mib", "failure_msg"]
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for row in results:
        writer.writerow(row)

# Write human-readable summary
txt_file = results_dir / "hypre_pool_study_summary.txt"
with open(txt_file, "w", encoding="utf-8") as f:
    f.write(f"{'POOL':<12} {'STATUS':<10} {'NEWTON':<10} {'KSP':<10} {'WALL(s)':<12} {'PEAK_GPU(MiB)':<15} {'EXIT':<6}\n")
    f.write("-" * 75 + "\n")
    for r in results:
        f.write(f"{r['case']:<12} {r['status']:<10} {r['newton_iterations']:<10} {r['ksp_iterations']:<10} {r['wall_seconds']:<12} {r['peak_gpu_memory_mib']:<15} {r['exit_code']:<6}\n")

# Print to stdout
print("\n" + "=" * 25 + " HYPRE POOL STUDY REPORT " + "=" * 25)
with open(txt_file, "r") as f:
    print(f.read())
print("Detailed CSV written to:", csv_file)
print("=" * 75)

EOF

echo "HYPRE pool study completed."
