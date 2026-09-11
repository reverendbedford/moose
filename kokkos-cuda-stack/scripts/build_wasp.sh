#!/bin/bash
# Step 2.5: build WASP into $PREFIX. Needed by MOOSE framework builds.
set -e
set -o pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
. "$SCRIPT_DIR/env.sh"

LOG="$LOGS/wasp-$(date +%Y%m%d-%H%M%S).log"
echo "[build_wasp] logging to $LOG"

export WASP_PREFIX="$PREFIX"

"$MOOSE_DIR/scripts/update_and_rebuild_wasp.sh" --skip-submodule-update 2>&1 | tee "$LOG"

echo
echo "[build_wasp] Installed:"
ls -1 "$PREFIX/lib/libwasp"*.so* 2>/dev/null | head
