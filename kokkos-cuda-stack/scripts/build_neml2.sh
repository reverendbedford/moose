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

# --- 2) libtorch + Python torch built from MOOSE's pinned pytorch commit
# NEML2's C++ code uses very recent torch APIs (torch/csrc/stable/library.h
# and 5-arg torch::inductor::AOTIModelPackageLoader) that don't exist in
# any pip wheel we can install on this box: cu124 wheels cap at torch
# 2.7.0.dev (nightly), and torch versions new enough for NEML2 (2.13-class)
# only ship as cu126+ wheels needing driver 560+ (this box is 550.x). The
# MOOSE-blessed workaround is to build from the pinned pytorch commit at
# framework/contrib/pytorch (currently on release/2.13, July 2026), which
# has those APIs and builds against our CUDA 12.4. See MOOSE's
# scripts/update_and_rebuild_libtorch.sh + scripts/configure_libtorch.sh.
#
# Idempotency: skip the ~60-90 min rebuild when the venv already has a
# torch that (a) is CUDA-capable, (b) has the C++11 ABI, and (c) ships
# the torch/csrc/stable/library.h header NEML2 will #include.
have_correct_torch() {
  python3 - <<'PY' 2>/dev/null
import os, sys
try:
    import torch
except Exception:
    sys.exit(1)
if not torch.cuda.is_available():                                sys.exit(1)
if not getattr(torch._C, "_GLIBCXX_USE_CXX11_ABI", False):       sys.exit(1)
if not os.path.isfile(os.path.join(os.path.dirname(torch.__file__),
                      "include/torch/csrc/stable/library.h")):   sys.exit(1)
sys.exit(0)
PY
}

if have_correct_torch; then
  echo "[build_neml2] torch (CUDA + C++11 ABI + stable API) already installed in $NEML2_VENV:"
  python3 -c 'import torch; print(f"    torch={torch.__version__}  cuda={torch.version.cuda}  cxx11_abi={torch._C._GLIBCXX_USE_CXX11_ABI}")'
else
  echo "[build_neml2] building libtorch + Python torch from MOOSE's pinned pytorch commit"
  echo "[build_neml2] (framework/contrib/pytorch, release/2.13) -- expected ~60-90 min on 32 cores"

  # Init pytorch submodule (.gitmodules marks it `update = none`, same as neml2).
  if [ ! -f "$MOOSE_DIR/framework/contrib/pytorch/setup.py" ]; then
    echo "[build_neml2] initializing framework/contrib/pytorch submodule (overriding update=none)"
    git -C "$MOOSE_DIR" \
        -c submodule."framework/contrib/pytorch".update=checkout \
        submodule update --init --recursive framework/contrib/pytorch 2>&1 | tee -a "$LOG"
  fi

  # Remove any pip-installed torch first so update_and_rebuild_libtorch.sh's
  # pip install of the freshly-built package proceeds cleanly.
  python3 -m pip uninstall -y torch 2>&1 | tee -a "$LOG" || true

  # CUDA_HOME: pytorch's cmake auto-detects CUDA via CUDA_HOME; env.sh sets
  # CUDA_DIR only. Export CUDA_HOME so the build turns on USE_CUDA.
  export CUDA_HOME="$CUDA_DIR"

  # TORCH_CUDA_ARCH_LIST=8.6 restricts CUDA compilation to Ampere sm_86
  # (RTX A5000 on this box). Without this, pytorch compiles kernels for
  # many archs, tripling build time.
  export TORCH_CUDA_ARCH_LIST="8.6"

  cd "$MOOSE_DIR"
  scripts/update_and_rebuild_libtorch.sh --install-python-package 2>&1 | tee -a "$LOG"

  # Verify all three properties, same as the skip-check above.
  if ! have_correct_torch; then
    echo "[build_neml2] ERROR: libtorch build finished but verification failed:" >&2
    echo "[build_neml2]   torch import, cuda.is_available, C++11 ABI, or torch/csrc/stable/library.h" >&2
    exit 1
  fi
  python3 -c 'import torch; print(f"    torch={torch.__version__}  cuda={torch.version.cuda}  cxx11_abi={torch._C._GLIBCXX_USE_CXX11_ABI}")' | tee -a "$LOG"
fi

# --- 3) NEML2 Python build backends AND runtime deps ------------------
# update_and_rebuild_neml2.sh passes --no-deps to protect the pinned torch,
# so it also skips every OTHER dep NEML2 needs. Install the full set here:
#   scikit-build-core  build backend  (chosen by NEML2's pyproject.toml)
#   ninja              build tool     (scikit-build-core default generator)
#   pybind11           C++ bindings   (find_package(pybind11) at cmake time)
#   nmhit>=0.3.6       NEML2-flavored C++ HIT parser (find_package(nmhit)
#                      REQUIRED at cmake time; ships static lib + headers
#                      under <site-packages>/nmhit/{lib,include})
#   pybind11-stubgen   .pyi generator, listed in NEML2's cibuildwheel
#                      before-build hook
# Order matters only in the sense that build_neml2 assumes all these
# succeed before calling update_and_rebuild_neml2.sh.
echo "[build_neml2] ensuring NEML2 build+runtime deps (scikit-build-core, pybind11, ninja, nmhit, pybind11-stubgen)"
python3 -m pip install --upgrade \
  scikit-build-core \
  pybind11 \
  ninja \
  'nmhit>=0.3.6' \
  'pybind11-stubgen>=2.5.5' \
  2>&1 | tee -a "$LOG"

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
