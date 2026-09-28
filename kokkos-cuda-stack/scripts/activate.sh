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

if [ -n "${ZSH_VERSION:-}" ]; then
  _src="${(%):-%x}"
else
  _src="${BASH_SOURCE[0]}"
fi
_activate_dir=$(cd "$(dirname "$_src")" && pwd)
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

# PATH: force ordering $PREFIX/bin : $NEML2_VENV_BIN : $MINIFORGE_BIN : /usr/local/cuda/bin
# in front of whatever the caller inherited. The previous scheme inserted the
# venv path AFTER $PREFIX/bin via a bash substring replace on ${PATH}; when
# the pre-existing PATH did not have $PREFIX/bin at the very front (e.g. a
# conda profile put miniforge/bin first) the substitution silently missed and
# MOOSE's embedded Python booted against miniforge base without torch.
#
# The re-source-safe idempotent rebuild here always ends up with the correct
# ordering, whatever the caller inherited, and duplicates every entry only
# once so repeated sourcing is a no-op.
NEML2_VENV_BIN="$STACK_DIR/neml2-venv/bin"
MINIFORGE_BIN="/home/chenghau.yang/miniforge/bin"

_head="$PREFIX/bin"
if [ -x "$NEML2_VENV_BIN/python3" ]; then
  _head="$_head:$NEML2_VENV_BIN"
  export VIRTUAL_ENV="$STACK_DIR/neml2-venv"
fi
if [ -x "$MINIFORGE_BIN/python3-config" ]; then
  _head="$_head:$MINIFORGE_BIN"
fi
_head="$_head:/usr/local/cuda/bin"

_new_path="$_head"
_orig_ifs="$IFS"
IFS=':'
# Loop over the caller's PATH entries. Word-splitting on `:` gives us the entries
# in both bash and zsh. Guard against unset PATH.
set -f
for _p in ${PATH:-}; do
  case ":$_new_path:" in
    *":$_p:"*) ;;
    *) _new_path="$_new_path:$_p" ;;
  esac
done
set +f
IFS="$_orig_ifs"
export PATH="$_new_path"
unset _head _new_path _orig_ifs _p

: "${MOOSE_JOBS:=8}"
export MOOSE_JOBS

if [ -n "${PS1:-}" ] && [ -z "${MOOSE_KOKKOS_ACTIVATED:-}" ]; then
  echo "[activate.sh] moose-kokkos CUDA stack activated (PREFIX=$PREFIX)"
fi
export MOOSE_KOKKOS_ACTIVATED=1

unset _activate_dir
