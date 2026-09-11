#!/bin/bash
# Step 2: build libmesh against the CUDA/Kokkos PETSc installed in $PREFIX.
set -e
set -o pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
. "$SCRIPT_DIR/env.sh"

LOG="$LOGS/libmesh-$(date +%Y%m%d-%H%M%S).log"

# Skip check: if libmesh is already installed at $PREFIX, skip the ~20 min
# rebuild. FORCE_REBUILD=1 (set by all.sh --force) bypasses.
if [ "${FORCE_REBUILD:-0}" != "1" ] \
   && [ -x "$PREFIX/bin/libmesh-config" ] \
   && [ -f "$PREFIX/lib/libmesh_opt.so" ]; then
  echo "[build_libmesh] libmesh already installed at \$PREFIX; skipping."
  echo "[build_libmesh]   libmesh-config : $PREFIX/bin/libmesh-config"
  echo "[build_libmesh]   libmesh_opt.so : $($PREFIX/bin/libmesh-config --version 2>/dev/null | head -1)"
  echo "[build_libmesh]   to force rebuild: FORCE_REBUILD=1 $0   (or rm $PREFIX/bin/libmesh-config)"
  exit 0
fi

echo "[build_libmesh] logging to $LOG"

export PETSC_DIR="$PREFIX"
export PETSC_ARCH=""            # PREFIX-installed PETSc has no arch layer
export LIBMESH_DIR="$PREFIX"    # install libmesh into the same prefix

# NOTE: do NOT set LIBMESH_BUILD_DIR. MOOSE's `scripts/configure_libmesh.sh`
# hardcodes `cd "${SRC_DIR}/build"` regardless of the env var, so overriding
# LIBMESH_BUILD_DIR just makes `update_and_rebuild_libmesh.sh` (which does
# respect it) cd into one dir and configure_libmesh.sh cd into another,
# leaving `make` with no makefile. Let the default `$LIBMESH_SRC_DIR/build`
# be used — the script rm+mkdirs it for us.
unset LIBMESH_BUILD_DIR

# MOOSE ships update_and_rebuild_libmesh.sh; it handles configure + make + install.
# Args:
#   --skip-submodule-update  : don't touch the pinned libmesh submodule
"$MOOSE_DIR/scripts/update_and_rebuild_libmesh.sh" --skip-submodule-update 2>&1 | tee "$LOG"

echo
echo "[build_libmesh] Installed:"
ls -1 "$PREFIX/lib/libmesh"*.so* 2>/dev/null | head
