#!/bin/bash
# Step 6: reconfigure MOOSE with --with-kokkos=cuda pointing at the installed
# PETSc+libmesh+WASP, optionally --with-neml2, clean framework, rebuild
# framework + solid_mechanics, capability-check.
#
# NEML2 support: controlled by NEML2_SUPPORT (default 1). When 1, the venv
# $STACK_DIR/neml2-venv/bin is prepended AFTER $PREFIX/bin so:
#   - mpicxx / mpicc / mpif90 stay ours (CUDA-aware, first on PATH)
#   - python3 resolves to the venv's python (3.12, cu124-compatible), so
#     ./configure's `python3 -c "import neml2"` auto-detect succeeds
# The corresponding NEML2 install is done by build_neml2.sh. If NEML2_SUPPORT=1
# but the venv's python3 cannot `import neml2`, this script aborts with a
# clear hint rather than silently building without NEML2.
set -e
set -o pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
. "$SCRIPT_DIR/env.sh"

LOG="$LOGS/moose-$(date +%Y%m%d-%H%M%S).log"
echo "[build_moose] logging to $LOG"

CONFIGURE_ARGS=(--with-kokkos=cuda)

if [ "${NEML2_SUPPORT:-1}" = "1" ]; then
  NEML2_VENV_BIN="$STACK_DIR/neml2-venv/bin"
  if [ ! -x "$NEML2_VENV_BIN/python3" ]; then
    echo "[build_moose] ERROR: NEML2_SUPPORT=1 but venv python3 missing at $NEML2_VENV_BIN" >&2
    echo "[build_moose]        run $SCRIPT_DIR/build_neml2.sh (creates the venv and installs torch+neml2)" >&2
    echo "[build_moose]        or export NEML2_SUPPORT=0 to build without NEML2." >&2
    exit 1
  fi
  export VIRTUAL_ENV="$STACK_DIR/neml2-venv"
  export PATH="$PREFIX/bin:$NEML2_VENV_BIN:$PATH"
  if ! "$NEML2_VENV_BIN/python3" -c 'import neml2' 2>/dev/null; then
    echo "[build_moose] ERROR: NEML2_SUPPORT=1 but 'import neml2' fails in $NEML2_VENV_BIN/python3." >&2
    echo "[build_moose]        run: $SCRIPT_DIR/build_neml2.sh   (or export NEML2_SUPPORT=0 to build without NEML2)" >&2
    exit 1
  fi
  # NEML2 depends on libtorch. MOOSE requires the flag explicitly -- even
  # with --with-neml2 present, ./configure errors out with "NEML2 depends
  # on libtorch. Please enable libtorch support with --with-libtorch"
  # unless --with-libtorch is also passed. Bare --with-libtorch (=yes)
  # triggers auto-detection: LIBTORCH_DIR env var -> python3 sibling of
  # neml2's site-packages -> framework/contrib/pytorch/installed/. All
  # three land on a libtorch install with matching torch/csrc/stable/
  # headers on this box.
  CONFIGURE_ARGS+=(--with-libtorch --with-neml2)
  echo "[build_moose] NEML2 support enabled (venv python3 has neml2 installed)"
else
  echo "[build_moose] NEML2 support disabled (NEML2_SUPPORT=$NEML2_SUPPORT)"
fi

cd "$MOOSE_DIR"
./configure "${CONFIGURE_ARGS[@]}" 2>&1 | tee "$LOG"

# stale headers/lib metadata will confuse the framework build after configure change.
cd "$MOOSE_DIR/framework"
make clean 2>&1 | tee -a "$LOG" || true
make -j "$MOOSE_JOBS" 2>&1 | tee -a "$LOG"

cd "$MOOSE_DIR/modules/solid_mechanics"
make -j "$MOOSE_JOBS" 2>&1 | tee -a "$LOG"

echo
echo "[build_moose] Solid mechanics capability check:"
./solid_mechanics-opt --show-capabilities 2>/dev/null | tail -n +2 | head -n -1 > /tmp/cap.json
python3 -c "
import json
d = json.load(open('/tmp/cap.json'))
print(f\"kokkos.value = {d.get('kokkos',{}).get('value')}\")
print(f\"cuda.value   = {d.get('cuda',{}).get('value')}\")
print(f\"neml2.value  = {d.get('neml2',{}).get('value')}\")
"
