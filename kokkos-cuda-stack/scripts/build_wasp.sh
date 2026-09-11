#!/bin/bash
# Step 2.5: build WASP into $PREFIX. Needed by MOOSE framework builds.
set -e
set -o pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
. "$SCRIPT_DIR/env.sh"

LOG="$LOGS/wasp-$(date +%Y%m%d-%H%M%S).log"

# Skip check: if WASP is already installed at $PREFIX, skip the ~5 min
# rebuild. FORCE_REBUILD=1 (set by all.sh --force) bypasses.
if [ "${FORCE_REBUILD:-0}" != "1" ] \
   && ls "$PREFIX/lib/libwaspcore.so" >/dev/null 2>&1; then
  echo "[build_wasp] WASP already installed at \$PREFIX; skipping."
  echo "[build_wasp]   libs: $(ls $PREFIX/lib/libwasp*.so 2>/dev/null | wc -l) files"
  echo "[build_wasp]   to force rebuild: FORCE_REBUILD=1 $0   (or rm $PREFIX/lib/libwaspcore.so)"
  exit 0
fi

echo "[build_wasp] logging to $LOG"

export WASP_PREFIX="$PREFIX"

"$MOOSE_DIR/scripts/update_and_rebuild_wasp.sh" --skip-submodule-update 2>&1 | tee "$LOG"

echo
echo "[build_wasp] Installed:"
ls -1 "$PREFIX/lib/libwasp"*.so* 2>/dev/null | head
