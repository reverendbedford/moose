#!/bin/bash
# Step 6: reconfigure MOOSE with --with-kokkos=cuda pointing at the installed
# PETSc+libmesh+WASP, optionally --with-neml2, rebuild framework +
# solid_mechanics + test, capability-check both executables.
#
# NEML2 support: controlled by NEML2_SUPPORT (default 1). When 1, the venv
# $STACK_DIR/neml2-venv/bin is prepended AFTER $PREFIX/bin, and miniforge/bin
# AFTER that, so:
#   - mpicxx / mpicc / mpif90 stay ours (CUDA-aware, first on PATH)
#   - python3 resolves to the venv's python (3.12), so ./configure's
#     `python3 -c "import neml2"` auto-detect succeeds
#   - python3-config resolves to miniforge (Py3.12) so MOOSE's -lpython3.12
#     link (for libneml2_eager.so's Python API refs) picks the right version
# The corresponding NEML2 install is done by build_neml2.sh. If NEML2_SUPPORT=1
# but the venv's python3 cannot `import neml2`, this script aborts with a
# clear hint rather than silently building without NEML2.
#
# Stale-triplet cleanup: MOOSE builds embed the compiler's platform triplet
# (e.g. x86_64-pc-linux-gnu vs x86_64-conda-linux-gnu) in every .lo file
# name AND in libtool metadata. MOOSE's `make clean` only removes files
# matching the CURRENT triplet, so an earlier build under a different env
# (e.g. conda `moose`) leaves stale artifacts whose embedded conda paths
# libtool then re-uses when regenerating .la files -- resulting in the
# next MOOSE link line getting -L/miniforge/envs/moose/lib -lpython3.14
# etc., and undefined-reference failures at link time. Auto-detect any
# non-current-triplet artifacts and wipe them before build.
set -e
set -o pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
. "$SCRIPT_DIR/env.sh"

LOG="$LOGS/moose-$(date +%Y%m%d-%H%M%S).log"

EXE_SM="$MOOSE_DIR/modules/solid_mechanics/solid_mechanics-opt"
EXE_TEST="$MOOSE_DIR/test/moose_test-opt"

# --- capability check helper (used by skip-check AND at end) --------------
check_capabilities() {
  # $1 = executable path; prints kokkos/cuda/neml2 values; exits 0 if all ok
  local exe="$1"
  [ -x "$exe" ] || return 1
  "$exe" --show-capabilities 2>/dev/null | tail -n +2 | head -n -1 > /tmp/cap-check.json || return 1
  WITH_NEML2=${NEML2_SUPPORT:-1} python3 - <<'PY'
import json, os, sys
try:    d = json.load(open('/tmp/cap-check.json'))
except: sys.exit(1)
def ok(v): return v and v != 'false'
if not (ok(d.get('kokkos', {}).get('value')) and ok(d.get('cuda', {}).get('value'))):
    sys.exit(1)
if os.environ.get('WITH_NEML2', '1') == '1' and not ok(d.get('neml2', {}).get('value')):
    sys.exit(1)
sys.exit(0)
PY
}

# --- Skip check: both executables work and no stale artifacts -------------
if [ "${FORCE_REBUILD:-0}" != "1" ]; then
  STALE=$(find "$MOOSE_DIR/framework" "$MOOSE_DIR/test" "$MOOSE_DIR/modules/solid_mechanics" \
               -name '*conda-linux-gnu*' 2>/dev/null | head -1)
  if [ -z "$STALE" ] && check_capabilities "$EXE_SM" && check_capabilities "$EXE_TEST"; then
    echo "[build_moose] both solid_mechanics-opt and moose_test-opt already built with matching capabilities; skipping."
    python3 -c "import json; d=json.load(open('/tmp/cap-check.json')); \
      print(f'  kokkos.value = {d.get(\"kokkos\",{}).get(\"value\")}'); \
      print(f'  cuda.value   = {d.get(\"cuda\",{}).get(\"value\")}'); \
      print(f'  neml2.value  = {d.get(\"neml2\",{}).get(\"value\")}')"
    echo "[build_moose]   to force rebuild: FORCE_REBUILD=1 $0   (or rm $EXE_SM)"
    exit 0
  fi
fi

echo "[build_moose] logging to $LOG"

# --- NEML2 support wiring ------------------------------------------------
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
  MINIFORGE_BIN="/home/chenghau.yang/miniforge/bin"
  export PATH="$PREFIX/bin:$NEML2_VENV_BIN:$MINIFORGE_BIN:$PATH"
  if ! "$NEML2_VENV_BIN/python3" -c 'import neml2' 2>/dev/null; then
    echo "[build_moose] ERROR: NEML2_SUPPORT=1 but 'import neml2' fails in $NEML2_VENV_BIN/python3." >&2
    echo "[build_moose]        run: $SCRIPT_DIR/build_neml2.sh   (or export NEML2_SUPPORT=0 to build without NEML2)" >&2
    exit 1
  fi
  CONFIGURE_ARGS+=(--with-libtorch --with-neml2)
  echo "[build_moose] NEML2 support enabled (venv python3 has neml2 installed)"
else
  echo "[build_moose] NEML2 support disabled (NEML2_SUPPORT=$NEML2_SUPPORT)"
fi

# --- Stale-triplet cleanup ----------------------------------------------
# If a previous build used a different platform triplet (e.g. conda's
# x86_64-conda-linux-gnu), MOOSE's `make clean` will not remove its .lo
# files, and libtool will re-embed the conda-era dependency_libs (like
# -L/miniforge/envs/moose/lib -lpython3.14) into the newly-linked .la
# files. That leaks into every downstream link and breaks with undefined
# references. Detect and nuke.
STALE_FILES=$(find "$MOOSE_DIR/framework" "$MOOSE_DIR/test" "$MOOSE_DIR/modules/solid_mechanics" \
                   -name '*conda-linux-gnu*' 2>/dev/null | wc -l)
if [ "$STALE_FILES" -gt 0 ]; then
  echo "[build_moose] detected $STALE_FILES stale conda-triplet build artifacts -- wiping them + all libtool artifacts"
  for d in framework framework/contrib/hit framework/contrib/pcre test modules/solid_mechanics; do
    find "$MOOSE_DIR/$d" \( -name '*.la' -o -name '*.lai' -o -name '*.lo' -o -name '*.o' \) -delete 2>/dev/null || true
    find "$MOOSE_DIR/$d" -name '*conda-linux-gnu*' -delete 2>/dev/null || true
  done
  rm -rf "$MOOSE_DIR/framework/build" "$MOOSE_DIR/test/build" "$MOOSE_DIR/modules/solid_mechanics/build"
  rm -rf "$MOOSE_DIR/framework/.libs" "$MOOSE_DIR/test/lib/.libs" "$MOOSE_DIR/modules/solid_mechanics/lib/.libs"
  rm -f "$MOOSE_DIR/framework/libmoose"*.la "$MOOSE_DIR/framework/libmoose"*.so*
  rm -f "$MOOSE_DIR/framework/contrib/hit/libhit"*.la "$MOOSE_DIR/framework/contrib/hit/libhit"*.so*
  rm -f "$MOOSE_DIR/framework/contrib/pcre/libpcre"*.la "$MOOSE_DIR/framework/contrib/pcre/libpcre"*.so*
  rm -f "$EXE_TEST" "$MOOSE_DIR/test/lib/libmoose_test"*.la "$MOOSE_DIR/test/lib/libmoose_test"*.so* "$MOOSE_DIR/test/lib/dlink.o"
  rm -f "$EXE_SM" "$MOOSE_DIR/modules/solid_mechanics/lib/libsolid_mechanics"*.la "$MOOSE_DIR/modules/solid_mechanics/lib/libsolid_mechanics"*.so*
fi

# --- Configure + build ---------------------------------------------------
cd "$MOOSE_DIR"
./configure "${CONFIGURE_ARGS[@]}" 2>&1 | tee "$LOG"

# make clean framework/ picks up any leftover stale-triplet-agnostic bits.
cd "$MOOSE_DIR/framework"
make clean 2>&1 | tee -a "$LOG" || true
make -j "$MOOSE_JOBS" 2>&1 | tee -a "$LOG"

cd "$MOOSE_DIR/modules/solid_mechanics"
make -j "$MOOSE_JOBS" 2>&1 | tee -a "$LOG"

# Also rebuild test/moose_test-opt so `test/run_tests` works out of the box.
# Same env, same flags, same libtorch/neml2 wiring -- no reason not to.
cd "$MOOSE_DIR/test"
make -j "$MOOSE_JOBS" 2>&1 | tee -a "$LOG"

# --- Capability checks ---------------------------------------------------
echo
echo "[build_moose] --- solid_mechanics-opt capabilities ---"
"$EXE_SM" --show-capabilities 2>/dev/null | tail -n +2 | head -n -1 > /tmp/cap.json
python3 -c "
import json
d = json.load(open('/tmp/cap.json'))
print(f\"  kokkos.value = {d.get('kokkos',{}).get('value')}\")
print(f\"  cuda.value   = {d.get('cuda',{}).get('value')}\")
print(f\"  neml2.value  = {d.get('neml2',{}).get('value')}\")
"

echo
echo "[build_moose] --- moose_test-opt capabilities ---"
"$EXE_TEST" --show-capabilities 2>/dev/null | tail -n +2 | head -n -1 > /tmp/cap-test.json
python3 -c "
import json
d = json.load(open('/tmp/cap-test.json'))
print(f\"  kokkos.value = {d.get('kokkos',{}).get('value')}\")
print(f\"  cuda.value   = {d.get('cuda',{}).get('value')}\")
print(f\"  neml2.value  = {d.get('neml2',{}).get('value')}\")
"
