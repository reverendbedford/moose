#!/bin/bash
# Step 6: reconfigure MOOSE with --with-kokkos=cuda pointing at the installed
# PETSc+libmesh+WASP, optionally --with-neml2, clean framework, rebuild
# framework + solid_mechanics, capability-check.
#
# NEML2 support: controlled by NEML2_SUPPORT (default 1). When 1, the conda
# `moose-neml2` env bin is prepended AFTER $PREFIX/bin so:
#   - mpicxx / mpicc / mpif90 stay ours (CUDA-aware, first on PATH)
#   - python3 resolves to moose-neml2 env's Python 3.13, so ./configure's
#     `python3 -c "import neml2"` auto-detect succeeds
# The corresponding NEML2 install is done by build_neml2.sh. If NEML2_SUPPORT=1
# but moose-neml2's python3 cannot `import neml2`, this script aborts with a
# clear hint rather than silently building without NEML2.
set -e
set -o pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
. "$SCRIPT_DIR/env.sh"

LOG="$LOGS/moose-$(date +%Y%m%d-%H%M%S).log"
echo "[build_moose] logging to $LOG"

CONFIGURE_ARGS=(--with-kokkos=cuda)

if [ "${NEML2_SUPPORT:-1}" = "1" ]; then
  CONDA_NEML2_BIN="/home/chenghau.yang/miniforge/envs/moose-neml2/bin"
  if [ ! -x "$CONDA_NEML2_BIN/python3" ]; then
    echo "[build_moose] ERROR: NEML2_SUPPORT=1 but conda moose-neml2 env python3 missing at $CONDA_NEML2_BIN" >&2
    echo "[build_moose]        create it: conda create -y -n moose-neml2 python=3.13 pip" >&2
    exit 1
  fi
  export PATH="$PREFIX/bin:$CONDA_NEML2_BIN:$PATH"
  if ! "$CONDA_NEML2_BIN/python3" -c 'import neml2' 2>/dev/null; then
    echo "[build_moose] ERROR: NEML2_SUPPORT=1 but 'import neml2' fails in conda moose-neml2 env." >&2
    echo "[build_moose]        run: $SCRIPT_DIR/build_neml2.sh   (or export NEML2_SUPPORT=0 to build without NEML2)" >&2
    exit 1
  fi
  CONFIGURE_ARGS+=(--with-neml2)
  echo "[build_moose] NEML2 support enabled (conda moose-neml2 python3 has neml2 installed)"
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
