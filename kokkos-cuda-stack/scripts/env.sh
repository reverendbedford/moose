#!/bin/bash
# Clean toolchain environment for the from-scratch MOOSE-CUDA build.
# Source this: `. scripts/env.sh` — do NOT execute.
#
# Design: aggressively erase every env var that (a) conda's activation
# stamps in (build_alias, CMAKE_ARGS, CXX_FOR_BUILD, CPP, GCC_RANLIB, all
# CONDA_*, SSL/CURL/REQUESTS_CA_BUNDLE, GSETTINGS_*, XML_CATALOG_FILES),
# (b) autotools/petsc/libmesh read to override compilers, or (c) leak
# conda-provided library paths into linker/preprocessor search. Then
# re-export only the small set of variables we actually want to be set.
# Anything not on the keep-list is unset.
#
# The `configure_libmesh.sh` step showed why this matters: with even one
# conda leftover (`build_alias=x86_64-conda-linux-gnu`), autoconf detected
# the build triplet as conda-linux-gnu and then went looking for
# `x86_64-conda-linux-gnu-mpicc`, which resolved back into the conda env.

KEEP_VARS="HOME USER LOGNAME PWD OLDPWD SHELL TERM DISPLAY XAUTHORITY MAIL LANG LC_ALL TMPDIR"
KEEP_VARS="$KEEP_VARS SSH_AUTH_SOCK SSH_CONNECTION SSH_CLIENT SSH_TTY XDG_RUNTIME_DIR"
KEEP_VARS="$KEEP_VARS MOOSE_CUDA_STACK_KEEP"  # escape hatch

# Extra vars the user may have set intentionally before sourcing env.sh
for extra in ${MOOSE_CUDA_STACK_KEEP:-}; do
  KEEP_VARS="$KEEP_VARS $extra"
done

# Build a Python one-liner that unsets everything not in KEEP_VARS.
# `env -i` inside a running shell is not an option — we need to keep the
# current shell alive for the caller — so we `unset` by name instead.
KEEP_SET=" $(echo $KEEP_VARS | tr -s ' ' | sed 's/ / /g') "
for var in $(compgen -e); do
  case " $KEEP_SET " in
    *" $var "*) : keep ;;
    *) unset "$var" ;;
  esac
done

# Anchor points
export MOOSE_DIR=/home/chenghau.yang/packages/moose-kokkos
export STACK_DIR=${MOOSE_DIR}/kokkos-cuda-stack
export PREFIX=${STACK_DIR}/prefix
export LOGS=${STACK_DIR}/logs

# Deterministic PATH. $PREFIX/bin first so downstream `command -v mpicc` in
# build_libmesh.sh, build_wasp.sh, and update_and_rebuild_libmesh.sh resolves
# to the CUDA-aware OpenMPI wrappers built by build_openmpi.sh. Nonexistent
# entries are harmless (the very first build_openmpi run has no $PREFIX/bin).
# Set BEFORE any external command (e.g. mkdir) is invoked below - the purge
# loop above cleared PATH, so a bare `mkdir` at this point would not resolve.
export PATH=$PREFIX/bin:/usr/local/cuda/bin:/usr/bin:/bin:/sbin:/usr/sbin

# OPAL_PREFIX: OpenMPI's own override for a relocated install. The
# mpicc/mpicxx/mpif90 wrappers are symlinks to opal_wrapper, which reads its
# compile flags from ${prefix}/share/openmpi/<lang>-wrapper-data.txt where
# ${prefix} is baked into the wrapper binary at OpenMPI build time. If the
# stack was moved after OpenMPI was installed, that baked prefix no longer
# exists and every wrapper invocation fails with "Cannot open configuration
# file ...wrapper-data.txt". OPAL_PREFIX makes opal_wrapper look at the
# current $PREFIX instead. On a fresh install it is redundant but harmless.
# build_openmpi.sh detects this stale-prefix condition and forces a rebuild
# so a future stack cleanup can drop this line entirely.
export OPAL_PREFIX="$PREFIX"

mkdir -p "$PREFIX" "$LOGS"

# PETSc / libmesh / MOOSE will use these values
export PETSC_SRC_DIR=${MOOSE_DIR}/petsc
export PETSC_ARCH=arch-scratch-cuda
export LIBMESH_SRC_DIR=${MOOSE_DIR}/libmesh

# CUDA
export CUDA_DIR=/usr/local/cuda

# WASP (built into $PREFIX by scripts/build_wasp.sh)
export WASP_DIR="$PREFIX"

# OpenMPI needs OMPI_FC because system /usr/bin/gfortran symlink is absent.
export OMPI_CC=/usr/bin/gcc-11
export OMPI_CXX=/usr/bin/g++-11
export OMPI_FC=/usr/bin/gfortran-9

# Parallel build. Machine has 32 cores; use 8 to leave headroom.
export MOOSE_JOBS=8

# Announce, if run interactively.
if [ -n "${PS1:-}" ]; then
  echo "[env.sh] purged non-safelist env; PATH=$PATH"
fi
