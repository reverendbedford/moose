#!/bin/bash
# Diagnostic runner for the staged step1..step4 NEML2/Kokkos/PETSc benchmarks.
#
# Runs each step ONCE on a small mesh with verbose PETSc + MOOSE timing, and
# prints per-step diagnostics chosen to answer:
#   1. Did the intended PETSc backends actually take effect?
#      In particular for Step 4: is the matrix really 'aijkokkos' and the
#      vector 'kokkos', or did libMesh silently install MATAIJ? '-options_left'
#      lists any option PETSc never consumed -- the smoking gun.
#   2. Did the solver converge, and with the same iteration counts across steps?
#   3. Where was the time spent (MOOSE PerfGraph + PETSc -log_view)?
#   4. Were there any errors / fatals we should surface?
#
# Every step's full stdout+stderr is kept under $RESULTS_DIR/<step>.log so we
# can dig further if any of the summary lines look off. This script does NOT
# use Nsight Systems; it is meant for a fast smoke test / correctness check.
#
# Usage:
#   EXE=/path/to/solid_mechanics-opt MESH_N=8 ./diagnose.sh
#
# Assumes the CUDA-enabled MOOSE stack is already active in the current shell
# (e.g. `source kokkos-cuda-stack/scripts/activate.sh`).

set -uo pipefail  # NOT -e: keep going even if one step fails

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MOOSE_DIR=$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)
EXE=${EXE:-$MOOSE_DIR/modules/solid_mechanics/solid_mechanics-opt}
RESULTS_DIR=${RESULTS_DIR:-$SCRIPT_DIR/diag_results}
MESH_N=${MESH_N:-8}

# Auto-source the CUDA MOOSE stack activation. Even if the user did it in the
# outer shell, that source is a no-op under zsh because activate.sh relies on
# ${BASH_SOURCE[0]} to find itself. This script IS bash (shebang), so BASH_SOURCE
# works here and the export propagates to the executable we launch below.
ACTIVATE="$MOOSE_DIR/kokkos-cuda-stack/scripts/activate.sh"
if [ -r "$ACTIVATE" ]; then
  # shellcheck disable=SC1090
  source "$ACTIVATE"
fi

if [ ! -x "$EXE" ]; then
  echo "ERROR: solid_mechanics executable not found at $EXE" >&2
  echo "Set EXE=/path/to/solid_mechanics-opt and rerun." >&2
  exit 1
fi

# Preflight: NEML2 [eager=true] embeds a CPython interpreter that must be able
# to `import torch`. Catch a broken environment here instead of after four
# identical crashes.
if ! python3 -c "import torch" >/dev/null 2>&1; then
  echo "ERROR: python3 cannot 'import torch' even after sourcing:" >&2
  echo "         $ACTIVATE" >&2
  echo "       VIRTUAL_ENV=${VIRTUAL_ENV:-<unset>}" >&2
  echo "       which python3: $(command -v python3 || echo '(none)')" >&2
  echo "       Check that $MOOSE_DIR/kokkos-cuda-stack/neml2-venv exists" >&2
  echo "       and has torch installed (build_neml2.sh installs it)." >&2
  exit 1
fi

echo "-- Activated: VIRTUAL_ENV=${VIRTUAL_ENV:-<unset>}"
echo "-- python3 : $(command -v python3)"
echo "-- torch   : $(python3 -c 'import torch; print(torch.__version__, "cuda_available="+str(torch.cuda.is_available()))' 2>&1)"

mkdir -p "$RESULTS_DIR"

# PETSc diagnostic options common to every step.
# -ksp_view          : dump KSP/PC/Mat/Vec config (shows resolved types)
# -snes_view         : dump SNES config
# -*_converged_reason: one-line convergence status per solve
# -*_monitor         : per-iteration residual norms (nonlinear + linear)
# -log_view          : PETSc event timing at PetscFinalize
# -options_left      : list any options set but not consumed (critical!)
#
# These are passed on the MOOSE command line, NOT via the PETSC_OPTIONS env
# var. MOOSE's petscSetOptions() (framework/src/utils/PetscSupport.C) calls
# PetscOptionsClear(NULL) before installing its own options, which wipes
# anything PETSC_OPTIONS put in the option database; it then re-adds options
# from the MOOSE command line via addPetscOptionsFromCommandline(). So the
# command line is the only channel that survives.
#
# The SNES/KSP diagnostic flags stay UNPREFIXED: previous runs confirm the
# unprefixed variants (-ksp_view etc.) are consumed and produce output, so
# adding 'nl0_' would just make them show up in 'options you set that were
# not used'.
COMMON_DIAG="-snes_view -snes_converged_reason -snes_monitor \
             -ksp_view  -ksp_converged_reason  -ksp_monitor \
             -log_view -options_left"

# MOOSE guards -mat_type behind a solver-system prefix ('nl0' here); the
# option DB call for matrices is prefixed so unprefixed -mat_type triggers:
#   "Setting option '-mat_type' is not supported without a solver-system
#    prefix. Use '-nl0_mat_type' for assembled libMesh matrices."
# -vec_type, by contrast, is applied by PETSc's VecSetFromOptions() with NO
# system prefix, so it MUST be passed unprefixed; -nl0_vec_type would end
# up ignored (see framework/src/utils/PetscSupport.C).
CPU_OPT="-vec_type standard -nl0_mat_type aij       -use_gpu_aware_mpi 0 $COMMON_DIAG"
GPU_OPT="-vec_type kokkos   -nl0_mat_type aijkokkos -use_gpu_aware_mpi 0 $COMMON_DIAG"

steps=(step1_plasticity_cpu_neml2
       step2_plasticity_gpu_neml2
       step3a_plasticity_cpu_neml2_kokkos_cpu_petsc
       step3b_plasticity_gpu_neml2_kokkos_cpu_petsc
       step4_plasticity_full_gpu)

banner() {
  echo
  printf '=%.0s' {1..72}; echo
  echo "== $*"
  printf '=%.0s' {1..72}; echo
}

env_snapshot() {
  echo "-- Environment"
  echo "  EXE          : $EXE"
  echo "  RESULTS_DIR  : $RESULTS_DIR"
  echo "  MESH_N       : $MESH_N"
  echo "  mpirun       : $(command -v mpirun || echo '(none)')"
  echo "  nsys         : $(command -v nsys  || echo '(none)')"
  if command -v nvidia-smi >/dev/null 2>&1; then
    echo "  GPUs (before run):"
    nvidia-smi --query-gpu=index,name,memory.used,utilization.gpu \
               --format=csv,noheader | sed 's/^/    /'
  else
    echo "  nvidia-smi   : (none)"
  fi
}

diag_from_log() {
  local log=$1

  echo "  -- PETSc types resolved (from -ksp_view header + type lines):"
  # Print object header + type: line in the order they appear. Loose whitespace
  # so nested objects (indented) are captured too.
  local ksp_types
  ksp_types=$(grep -E "^[[:space:]]*(KSP|PC|Mat|Vec) Object|^[[:space:]]+type:" "$log" | head -30)
  if [ -n "$ksp_types" ]; then
    echo "$ksp_types" | sed 's/^/     /'
  else
    echo "     (no ksp_view output found -- did the solver run?)"
  fi

  echo "  -- PETSc options NOT consumed (from -options_left):"
  # PETSc prints "WARNING! There are options you set that were not used:"
  # followed by an indented list, then a blank line.
  local unused
  unused=$(awk '/There are options you set that were not used/{f=1; next}
                f && NF==0{exit}
                f {print}' "$log")
  if [ -n "$unused" ]; then
    echo "$unused" | sed 's/^/     /'
    echo "     >>> These options were IGNORED. For Step 4, seeing -mat_type"
    echo "     >>> or -vec_type here means libMesh hard-set the type and the"
    echo "     >>> run is effectively still CPU-PETSc."
  else
    echo "     (all PETSc options consumed)"
  fi

  echo "  -- Convergence tallies:"
  local nl_ok nl_bad ksp_ok ksp_bad
  # `grep -c` prints "0" AND exits 1 on no match; use `|| true` (not `|| echo 0`)
  # so we don't append a second "0" line.
  nl_ok=$( grep -c "Nonlinear solve converged"        "$log" 2>/dev/null || true)
  nl_bad=$(grep -c "Nonlinear solve did not converge" "$log" 2>/dev/null || true)
  ksp_ok=$( grep -c "Linear solve converged"          "$log" 2>/dev/null || true)
  ksp_bad=$(grep -c "Linear solve did not converge"   "$log" 2>/dev/null || true)
  echo "     nonlinear: converged=${nl_ok:-0}  did_not_converge=${nl_bad:-0}"
  echo "     linear   : converged=${ksp_ok:-0}  did_not_converge=${ksp_bad:-0}"

  echo "  -- MOOSE error block (from '*** ERROR ***' to next separator, capped at 25 lines):"
  local err_block
  err_block=$(awk '/\*\*\* ERROR/{f=1}
                   f{print; n++}
                   f && /^-----/ && n>1 {exit}
                   n>=25 {exit}' "$log")
  if [ -n "$err_block" ]; then
    echo "$err_block" | sed 's/^/     /'
  else
    echo "     (no MOOSE *** ERROR *** block found)"
  fi

  echo "  -- Other fatal-ish lines (first 5, dedup):"
  local other
  other=$(grep -Ei "FATAL|Segmentation|libMesh terminated|traceback|libc\+\+abi|terminate called" "$log" | sort -u | head -5)
  if [ -n "$other" ]; then
    echo "$other" | sed 's/^/     /'
  else
    echo "     (none)"
  fi

  echo "  -- MOOSE time step summary (last 6 status lines):"
  grep -E "Time Step|Solve Converged|Solve Did NOT Converge" "$log" \
    | tail -6 | sed 's/^/     /'

  echo "  -- MOOSE PerfGraph (first 20 lines of --timing block):"
  local pg
  pg=$(awk '/Performance Graph/{f=1} f{print; if(++n>=20) exit}' "$log")
  if [ -n "$pg" ]; then
    echo "$pg" | sed 's/^/     /'
  else
    echo "     (no PerfGraph found in log)"
  fi

  echo "  -- PETSc -log_view Summary of Stages (first 25 lines):"
  local lv
  lv=$(awk '/Summary of Stages/{f=1} f{print; if(++n>=25) exit}' "$log")
  if [ -n "$lv" ]; then
    echo "$lv" | sed 's/^/     /'
  else
    echo "     (no -log_view output found)"
  fi
}

run_step() {
  local step=$1
  local input="$SCRIPT_DIR/$step.i"
  local petsc=$CPU_OPT
  local device_args=()

  if [[ "$step" == step3a_* || "$step" == step3b_* || "$step" == step4_* ]]; then
    device_args=(--compute-device=cuda)
  fi
  if [[ "$step" == step4_* ]]; then
    petsc=$GPU_OPT
  fi

  local prefix="$RESULTS_DIR/$step"
  banner "$step"
  echo "  input        : $input"
  echo "  device flag  : ${device_args[*]:-<none>}"
  echo "  PETSc args   : $petsc"

  # $petsc is intentionally unquoted so its flags become separate argv entries.
  # See the comment on COMMON_DIAG for why these must be on the command line
  # rather than in the PETSC_OPTIONS env var.
  /usr/bin/time -f 'wall=%es user=%Us sys=%Ss maxRSS=%MkB' -o "$prefix.time" \
    "$EXE" -i "$input" "${device_args[@]}" \
      "N=$MESH_N" \
      "Outputs/perf_graph_json_file=$prefix.perf.json" \
      --timing \
      $petsc \
      > "$prefix.log" 2>&1
  local rc=$?

  echo "  exit code    : $rc"
  echo "  timing       : $(cat "$prefix.time" 2>/dev/null || echo '?')"
  diag_from_log "$prefix.log"
}

banner "step1..step4 diagnostic runner"
env_snapshot
for s in "${steps[@]}"; do
  run_step "$s"
done

banner "Summary"
printf "  %-55s %s\n" "step" "wall / user / sys / maxRSS"
for s in "${steps[@]}"; do
  t=$(cat "$RESULTS_DIR/$s.time" 2>/dev/null || echo "?")
  printf "  %-55s %s\n" "$s" "$t"
done

echo
echo "  Raw logs   : $RESULTS_DIR/*.log"
echo
echo "  Step-4 sanity: this should list aijkokkos + kokkos (NOT aij + standard):"
echo "    grep -E '^[[:space:]]*(Mat|Vec) Object|^[[:space:]]+type:' \\"
echo "      $RESULTS_DIR/step4_plasticity_full_gpu.log | head -20"
