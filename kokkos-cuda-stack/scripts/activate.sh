#!/bin/bash
# Lightweight activation for the moose-kokkos CUDA-Kokkos stack.
# Intended for DAILY USE (running solid_mechanics-opt, `make -j` in
# framework/ or modules/*/) and for sourcing from direnv (.envrc).
#
# Compared to env.sh:
#   - Exports only the USE-side vars that MOOSE Makefiles read at build
#     time and MOOSE reads at runtime: PATH, PETSC_DIR, PETSC_ARCH="",
#     LIBMESH_DIR, WASP_DIR, OPAL_PREFIX.
#   - Does NOT purge existing environment. env.sh purges because autoconf
#     misreads conda leftovers (build_alias, CMAKE_ARGS, CXX_FOR_BUILD,
#     GCC_RANLIB, etc.); that only matters for full-stack builds via
#     scripts/all.sh + scripts/build_*.sh, all of which source env.sh
#     themselves. Daily `make -j` on already-installed deps does not
#     need the purge, and losing conda tools in the shell is annoying.
#   - Does NOT set OMPI_CC/CXX/FC. The default compilers were baked into
#     OpenMPI at build time (gcc-11/g++-11/gfortran-9), so unsetting the
#     wrapper overrides is fine for USE.
#   - Idempotent: safe to source multiple times; PATH is deduplicated.
#
# For full stack rebuilds, use scripts/all.sh or scripts/build_*.sh.
# Do NOT source this file INSTEAD of env.sh for a build.

_activate_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export MOOSE_DIR="${_activate_dir%/kokkos-cuda-stack/scripts}"
export STACK_DIR="$MOOSE_DIR/kokkos-cuda-stack"
export PREFIX="$STACK_DIR/prefix"

# USE-side vars (parallel to env.sh's USE block).
export PETSC_DIR="$PREFIX"
export PETSC_ARCH=""
export LIBMESH_DIR="$PREFIX"
export WASP_DIR="$PREFIX"

# opal_wrapper (mpicc/mpicxx/mpif90) reads its config from
# ${prefix}/share/openmpi/*-wrapper-data.txt where ${prefix} is baked in
# at OpenMPI build time. Currently baked-in path is the pre-relocation
# sibling folder; OPAL_PREFIX makes the wrapper look at the current
# $PREFIX. Redundant after a fresh build_openmpi.sh; harmless.
export OPAL_PREFIX="$PREFIX"

# PATH: prepend $PREFIX/bin (mpicc/mpirun/libmesh-config), add CUDA bin
# (nvcc/nsys). All entries guarded so re-sourcing does not duplicate.
case ":$PATH:" in
  *":$PREFIX/bin:"*) ;;
  *) PATH="$PREFIX/bin:$PATH" ;;
esac
case ":$PATH:" in
  *":/usr/local/cuda/bin:"*) ;;
  *) PATH="$PATH:/usr/local/cuda/bin" ;;
esac

# NEML2 support (only when the venv is present): matches the PATH order
# build_moose.sh uses so `cd anywhere && make -j` produces a binary with
# the same RUNPATH as `scripts/build_moose.sh` would.
#   $NEML2_VENV_BIN  -> python3 (venv, has `import neml2`)
#   miniforge/bin    -> python3-config (miniforge base, Python 3.12,
#                       needed for MOOSE's `python3-config --embed`
#                       -lpython3.12 rather than /usr/bin's Py3.10)
# Both inserted AFTER $PREFIX/bin so mpicxx stays ours.
NEML2_VENV_BIN="$STACK_DIR/neml2-venv/bin"
if [ -x "$NEML2_VENV_BIN/python3" ]; then
  case ":$PATH:" in
    *":$NEML2_VENV_BIN:"*) ;;
    *)
      # Insert right after $PREFIX/bin
      PATH="${PATH/$PREFIX\/bin:/$PREFIX/bin:$NEML2_VENV_BIN:}"
      ;;
  esac
  export VIRTUAL_ENV="$STACK_DIR/neml2-venv"
fi

MINIFORGE_BIN="/home/chenghau.yang/miniforge/bin"
if [ -x "$MINIFORGE_BIN/python3-config" ]; then
  case ":$PATH:" in
    *":$MINIFORGE_BIN:"*) ;;
    *)
      # Insert right after venv bin (or $PREFIX/bin if venv absent)
      if [ -n "${VIRTUAL_ENV:-}" ]; then
        PATH="${PATH/$NEML2_VENV_BIN:/$NEML2_VENV_BIN:$MINIFORGE_BIN:}"
      else
        PATH="${PATH/$PREFIX\/bin:/$PREFIX/bin:$MINIFORGE_BIN:}"
      fi
      ;;
  esac
fi

export PATH

: "${MOOSE_JOBS:=8}"
export MOOSE_JOBS

if [ -n "${PS1:-}" ] && [ -z "${MOOSE_KOKKOS_ACTIVATED:-}" ]; then
  echo "[activate.sh] moose-kokkos CUDA stack activated (PREFIX=$PREFIX)"
fi
export MOOSE_KOKKOS_ACTIVATED=1

unset _activate_dir
