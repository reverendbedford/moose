#!/bin/bash
# End-to-end from-scratch CUDA-Kokkos stack build.
#
# Steps and typical wall times (32-core box, MOOSE_JOBS=8):
#   openmpi   ~10 min  (skipped if $PREFIX/bin/ompi_info reports CUDA at the
#                       CURRENT $PREFIX; auto-repairs a stale prefix if the
#                       install was moved -- see build_openmpi.sh)
#   petsc     ~40 min  (always rebuilds; wipes $PETSC_ARCH first)
#   libmesh   ~20 min  (always rebuilds)
#   wasp      ~5 min
#   neml2     ~10-30 min (installs CUDA torch + NEML2 into conda moose env)
#   moose     ~10 min
#   benchmark ~15 min
#
# Flags:
#   --no-neml2   skip the neml2 step and configure MOOSE without --with-neml2.
#                Default: NEML2 is installed and enabled.
#
# Resume policy: if a step fails, do NOT re-run all.sh from scratch (that
# throws away all successful earlier steps -- notably a 40 min PETSc build).
# Run the individual build_*.sh script for the failed step and every step
# after it. The failure message below prints the exact command to paste.

set -e
set -o pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
. "$SCRIPT_DIR/env.sh" >/dev/null

# --- Flag parsing --------------------------------------------------------
WITH_NEML2=1
for arg in "$@"; do
  case "$arg" in
    --no-neml2) WITH_NEML2=0 ;;
    -h|--help)
      cat <<HELP
Usage: $0 [--no-neml2]

  --no-neml2   Skip NEML2 install; configure MOOSE without --with-neml2.
               Default: NEML2 is installed and enabled.
HELP
      exit 0 ;;
    *) echo "$0: unknown flag: $arg (use --help)" >&2; exit 1 ;;
  esac
done

# Passed to build_moose.sh so its ./configure line matches what we built.
export NEML2_SUPPORT=$WITH_NEML2

if [ "$WITH_NEML2" = 1 ]; then
  STEPS=(openmpi petsc libmesh wasp neml2 moose benchmark)
else
  STEPS=(openmpi petsc libmesh wasp moose benchmark)
fi

# --- Submodule init ------------------------------------------------------
# petsc / libmesh / wasp are regular submodules (update=checkout). neml2 is
# marked `update = none` in .gitmodules so a plain `git submodule update
# --init` silently skips it ("Skipping submodule 'framework/contrib/neml2'").
# Override the update strategy for neml2 only, for this invocation only.
SUBS="petsc libmesh framework/contrib/wasp"
GIT_ARGS=()
if [ "$WITH_NEML2" = 1 ]; then
  SUBS="$SUBS framework/contrib/neml2"
  GIT_ARGS+=(-c submodule."framework/contrib/neml2".update=checkout)
fi
echo "[all.sh] git submodule update --init --recursive $SUBS"
git -C "$MOOSE_DIR" "${GIT_ARGS[@]}" submodule update --init --recursive $SUBS

script_for() {
  case "$1" in
    benchmark) echo "$SCRIPT_DIR/run_benchmark.sh" ;;
    *)         echo "$SCRIPT_DIR/build_$1.sh" ;;
  esac
}

fmt_time() {
  local s=$1
  if [ "$s" -ge 60 ]; then
    printf '%dm%02ds' $((s/60)) $((s%60))
  else
    printf '%ds' "$s"
  fi
}

# --- Startup banner -------------------------------------------------------
cat <<EOF
=============================================================================
  MOOSE-Kokkos-CUDA stack build
=============================================================================
  STACK_DIR  : $STACK_DIR
  PREFIX     : $PREFIX
  MOOSE_DIR  : $MOOSE_DIR
  CUDA_DIR   : $CUDA_DIR
  MOOSE_JOBS : $MOOSE_JOBS
  LOGS       : $LOGS
  NEML2      : $([ "$WITH_NEML2" = 1 ] && echo "on (conda moose env)" || echo "off (--no-neml2)")
  Steps      : ${STEPS[*]}
=============================================================================
EOF

t0=$SECONDS

# --- Per-step execution with timing + resume hint on failure -------------
for step in "${STEPS[@]}"; do
  script=$(script_for "$step")
  echo
  echo "[all.sh] === $step (total elapsed: $(fmt_time $((SECONDS - t0)))) ==="
  ts=$SECONDS
  if ! "$script"; then
    dt=$(fmt_time $((SECONDS - ts)))
    latest_log=$(ls -1t "$LOGS/${step}-"*.log 2>/dev/null | head -1)
    remaining=()
    for s in "${STEPS[@]}"; do
      if [ "$s" = "$step" ] || [ "$reached" = 1 ]; then
        remaining+=("$(script_for "$s")")
        reached=1
      fi
    done
    cat >&2 <<EOF

[all.sh] STEP '$step' FAILED after $dt.
[all.sh] latest log: ${latest_log:-'(none written)'}
[all.sh]
[all.sh] To resume WITHOUT redoing successful earlier steps, run:
[all.sh]
[all.sh]   $(printf '%s && \\\n  ' "${remaining[@]}" | sed '$s/ && \\$//')
[all.sh]
[all.sh] (Do NOT re-run all.sh -- build_petsc.sh unconditionally wipes and
[all.sh]  rebuilds \$PETSC_ARCH, so a re-run costs ~40 min of PETSc work.)
EOF
    exit 1
  fi
  echo "[all.sh] $step done in $(fmt_time $((SECONDS - ts)))"
done

# --- Final capability check ---------------------------------------------
total=$(fmt_time $((SECONDS - t0)))
echo
echo "=============================================================================
  ALL STEPS DONE in $total
============================================================================="
EXE="$MOOSE_DIR/modules/solid_mechanics/solid_mechanics-opt"
if [ -x "$EXE" ]; then
  "$EXE" --show-capabilities 2>/dev/null | tail -n +2 | head -n -1 > /tmp/all-cap.json
  WITH_NEML2=$WITH_NEML2 python3 - <<'PY'
import json, os, sys
try:
    d = json.load(open('/tmp/all-cap.json'))
    def v(name): return d.get(name, {}).get('value')
    print(f"  kokkos.value = {v('kokkos')}")
    print(f"  cuda.value   = {v('cuda')}")
    want_neml2 = os.environ.get('WITH_NEML2', '1') == '1'
    if want_neml2:
        print(f"  neml2.value  = {v('neml2')}")
    ok = v('kokkos') and v('kokkos') != 'false' \
         and v('cuda')   and v('cuda')   != 'false'
    if want_neml2:
        ok = ok and v('neml2') and v('neml2') != 'false'
    caps = 'Kokkos+CUDA' + ('+NEML2' if want_neml2 else '')
    print(f"  status       = {'OK -- stack is ' + caps + ' capable' if ok else 'FAILED -- check --show-capabilities'}")
except Exception as e:
    print(f"  (capability parse failed: {e})", file=sys.stderr)
PY
else
  echo "  (solid_mechanics-opt not found at $EXE; skipping capability check)"
fi
