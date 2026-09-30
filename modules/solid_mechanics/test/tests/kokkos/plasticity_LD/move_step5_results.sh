#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SOURCE_DIR=${SOURCE_DIR:-$SCRIPT_DIR/results_step5}
DEST_DIR=${DEST_DIR:-$SCRIPT_DIR/results_n32_smoke}

for extension in log perf.json time; do
  source_file="$SOURCE_DIR/step5_plasticity_full_gpu.$extension"
  destination_file="$DEST_DIR/step5_plasticity_full_gpu_rep1.$extension"

  if [ ! -f "$source_file" ]; then
    echo "ERROR: missing $source_file" >&2
    exit 1
  fi

  mv -f "$source_file" "$destination_file"
done
