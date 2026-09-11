#!/bin/bash
# Step 3: reconfigure MOOSE with --with-kokkos=cuda pointing at the new PETSc+libmesh,
# clean framework, rebuild framework + solid_mechanics.
set -e
set -o pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
. "$SCRIPT_DIR/env.sh"

LOG="$LOGS/moose-$(date +%Y%m%d-%H%M%S).log"
echo "[build_moose] logging to $LOG"

export PETSC_DIR="$PREFIX"
export PETSC_ARCH=""
export LIBMESH_DIR="$PREFIX"

cd "$MOOSE_DIR"
./configure --with-kokkos=cuda 2>&1 | tee "$LOG"

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
"
