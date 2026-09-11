#!/bin/bash
# Step 5: install PyTorch (CUDA) and NEML2 into a dedicated venv.
#
# NEML2 installs as a Python package (pip). Rather than building libtorch
# from source (~1-2 hours) or reusing a conda env, we manage a dedicated
# venv at $STACK_DIR/neml2-venv built off miniforge's BASE python (3.12).
# The main `moose` conda env is on Python 3.14, but CUDA-12.4 torch wheels
# stop at Python 3.13, and driver 550.x on this box is too old for the
# CUDA-12.6 wheels that support 3.14. Miniforge base's 3.12 is the
# newest Python for which cu124 wheels exist AND is already installed
# on this box, so no `conda create` or `apt install` is needed. NEML2's
# C++ artifacts
# (libneml2*.so) are built with cmake using whatever CC/CXX/FC we hand it,
# so we pin them to $PREFIX/bin/mpi* (our from-scratch CUDA-aware OpenMPI)
# to guarantee ABI/MPI compatibility with MOOSE's solid_mechanics-opt.
#
# Idempotency:
#  - torch install: skip if `import torch` works AND torch.cuda.is_available().
#  - NEML2 install: update_and_rebuild_neml2.sh always rebuilds (pip
#    --force-reinstall). This is by design: NEML2 source changes frequently
#    on this branch and the pip wheel cache would otherwise serve stale bits.
#  - submodule init: guarded by presence of framework/contrib/neml2/CMakeLists.txt.
set -e
set -o pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
. "$SCRIPT_DIR/env.sh"

NEML2_VENV="$STACK_DIR/neml2-venv"
NEML2_VENV_BIN="$NEML2_VENV/bin"

# Auto-create the venv on first run. Uses miniforge base's python3 (3.12) as
# the interpreter -- newest cu124-compatible Python already installed on this
# box (main `moose` conda env is on 3.14, unusable with cu124 wheels).
if [ ! -x "$NEML2_VENV_BIN/python3" ]; then
  MINIFORGE_PYTHON="/home/chenghau.yang/miniforge/bin/python3"
  if [ ! -x "$MINIFORGE_PYTHON" ]; then
    echo "[build_neml2] ERROR: miniforge base python not found at $MINIFORGE_PYTHON" >&2
    echo "[build_neml2]        expected miniforge with a base python; adjust MINIFORGE_PYTHON here" >&2
    echo "[build_neml2]        or point the venv at any other Python >= 3.9, <= 3.13 interpreter." >&2
    exit 1
  fi
  echo "[build_neml2] creating venv at $NEML2_VENV (from $MINIFORGE_PYTHON, $($MINIFORGE_PYTHON --version 2>&1))"
  "$MINIFORGE_PYTHON" -m venv "$NEML2_VENV"
fi

# PATH: $PREFIX/bin first so mpicxx/mpicc/mpif90 are ours (CUDA-aware OpenMPI).
# Venv bin next so python3 and pip come from the venv. env.sh already put
# $PREFIX/bin and /usr/local/cuda/bin on PATH; we insert the venv bin between.
# VIRTUAL_ENV mirrors what standard `activate` would export (informational).
export VIRTUAL_ENV="$NEML2_VENV"
export PATH="$PREFIX/bin:$NEML2_VENV_BIN:$PATH"

# Pin compilers for cmake (scikit-build-core reads CC/CXX/FC). Without this,
# cmake would fall back to /usr/bin/cc (system gcc), producing a libneml2.so
# that links against a different libstdc++ than MOOSE's stack.
export CC="$PREFIX/bin/mpicc"
export CXX="$PREFIX/bin/mpicxx"
export FC="$PREFIX/bin/mpif90"

LOG="$LOGS/neml2-$(date +%Y%m%d-%H%M%S).log"
echo "[build_neml2] logging to $LOG"
echo "[build_neml2] python3 -> $(command -v python3) ($(python3 --version 2>&1))"
echo "[build_neml2] pip     -> $(command -v pip)"
echo "[build_neml2] cmake   -> $(command -v cmake) ($(cmake --version 2>&1 | head -1))"
echo "[build_neml2] mpicxx  -> $(command -v mpicxx)"
echo "[build_neml2] CC/CXX/FC pinned to \$PREFIX/bin/mpi*"

# --- 1) NEML2 submodule ------------------------------------------------
# MOOSE .gitmodules marks framework/contrib/neml2 with `update = none`, which
# even a plain `git submodule update --init` respects (silent skip -- would
# show "Skipping submodule 'framework/contrib/neml2'"). Override the update
# strategy with -c submodule.<path>.update=checkout for this invocation only.
NEML2_SRC="$MOOSE_DIR/framework/contrib/neml2"
if [ ! -f "$NEML2_SRC/CMakeLists.txt" ]; then
  echo "[build_neml2] initializing framework/contrib/neml2 submodule (overriding update=none)"
  git -C "$MOOSE_DIR" \
      -c submodule."framework/contrib/neml2".update=checkout \
      submodule update --init --recursive framework/contrib/neml2 2>&1 | tee -a "$LOG"
fi

# --- 2) PyTorch (CUDA 12.4 wheels) -------------------------------------
if python3 -c 'import torch, sys; sys.exit(0 if torch.cuda.is_available() else 1)' 2>/dev/null; then
  echo "[build_neml2] PyTorch with CUDA already installed in $NEML2_VENV:"
  python3 -c 'import torch; print(f"    torch={torch.__version__}  cuda={torch.version.cuda}  cuda_available={torch.cuda.is_available()}")'
else
  echo "[build_neml2] installing PyTorch (CUDA 12.4 wheels) into $NEML2_VENV"
  python3 -m pip install --upgrade pip 2>&1 | tee -a "$LOG"
  python3 -m pip install torch --index-url https://download.pytorch.org/whl/cu124 2>&1 | tee -a "$LOG"
  # Verify.
  if ! python3 -c 'import torch; assert torch.cuda.is_available(), "torch.cuda.is_available() is False"' 2>&1 | tee -a "$LOG"; then
    echo "[build_neml2] ERROR: torch installed but CUDA is not available." >&2
    echo "[build_neml2]        common cause: no CUDA-12.4-compatible wheel for python$(python3 -c 'import sys; print(f\"{sys.version_info.major}.{sys.version_info.minor}\")')." >&2
    exit 1
  fi
  python3 -c 'import torch; print(f"    torch={torch.__version__}  cuda={torch.version.cuda}  cuda_available={torch.cuda.is_available()}")' | tee -a "$LOG"
fi

# --- 3) NEML2 Python build backends ------------------------------------
# update_and_rebuild_neml2.sh checks for these but never installs them
# (--no-deps is used to protect the pinned torch). Install upfront.
# ninja: scikit-build-core defaults to the Ninja generator; without it,
# the wheel build aborts with NinjaNotFoundError. The `ninja` PyPI package
# ships a bundled ninja binary in the venv -- avoids apt install ninja-build.
echo "[build_neml2] ensuring NEML2 Python build deps (scikit-build-core, pybind11, ninja)"
python3 -m pip install --upgrade scikit-build-core pybind11 ninja 2>&1 | tee -a "$LOG"

# --- 4) NEML2 itself ---------------------------------------------------
echo "[build_neml2] pip-installing NEML2 into $NEML2_VENV (via MOOSE's update_and_rebuild_neml2.sh)"
cd "$MOOSE_DIR"
scripts/update_and_rebuild_neml2.sh --skip-submodule-update 2>&1 | tee -a "$LOG"

# --- 5) Report ---------------------------------------------------------
NEML2_PKG=$(python3 -c 'import neml2, os; print(os.path.dirname(neml2.__file__))' 2>/dev/null || echo "(not found)")
echo
echo "[build_neml2] NEML2 Python package: $NEML2_PKG"
echo "[build_neml2] NEML2 C++ headers/libs go into the same venv's include/ and lib/."
echo "[build_neml2] build_moose.sh will pick this up automatically when NEML2_SUPPORT=1."
