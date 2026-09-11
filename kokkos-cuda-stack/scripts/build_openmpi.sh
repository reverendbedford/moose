#!/bin/bash
# Step 0: build CUDA-aware OpenMPI 4.1.6 into $PREFIX.
#
# Rationale: system /usr/bin/mpicc wraps an OpenMPI that reports
# `opal_built_with_cuda_support false`. PETSc with -mat_type aijkokkos +
# -vec_type kokkos then aborts (MPI error 76) unless -use_gpu_aware_mpi 0 is
# passed. Building our own CUDA-aware OpenMPI at $PREFIX lets PETSc, libmesh
# and MOOSE all link against it and lets GPU-aware MPI stay on.
#
# Idempotency: if $PREFIX/bin/ompi_info exists and reports CUDA support,
# skip. Delete $PREFIX/bin/ompi_info (or the whole $PREFIX/bin/mpi* set) to
# force a rebuild.
set -e
set -o pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
. "$SCRIPT_DIR/env.sh"

LOG="$LOGS/openmpi-$(date +%Y%m%d-%H%M%S).log"

if [ -x "$PREFIX/bin/ompi_info" ] \
   && "$PREFIX/bin/ompi_info" | grep -q "opal_built_with_cuda_support:true\|MPI extensions:.*cuda"; then
  echo "[build_openmpi] $PREFIX/bin/ompi_info already reports CUDA support; skipping."
  echo "[build_openmpi] delete $PREFIX/bin/ompi_info to force a rebuild."
  exit 0
fi

echo "[build_openmpi] logging to $LOG"

OMPI_VER=4.1.6
OMPI_SRC="$STACK_DIR/src/openmpi-${OMPI_VER}"
OMPI_BUILD="$STACK_DIR/src/openmpi-${OMPI_VER}-build"
OMPI_TARBALL="$STACK_DIR/src/openmpi-${OMPI_VER}.tar.bz2"

mkdir -p "$STACK_DIR/src"

# Fetch tarball only if neither the source tree nor the tarball is present.
if [ ! -d "$OMPI_SRC" ] && [ ! -f "$OMPI_TARBALL" ]; then
  echo "[build_openmpi] fetching openmpi-${OMPI_VER}.tar.bz2" | tee -a "$LOG"
  curl -fL -o "$OMPI_TARBALL" \
    "https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-${OMPI_VER}.tar.bz2" \
    2>&1 | tee -a "$LOG"
fi

# Extract if the source tree is absent.
if [ ! -d "$OMPI_SRC" ]; then
  echo "[build_openmpi] extracting openmpi-${OMPI_VER}.tar.bz2" | tee -a "$LOG"
  tar -C "$STACK_DIR/src" -xjf "$OMPI_TARBALL" 2>&1 | tee -a "$LOG"
fi

# Fresh build tree every time we get here (i.e. every time we're actually
# (re)building). Guards against half-configured state from a killed run.
rm -rf "$OMPI_BUILD"
mkdir -p "$OMPI_BUILD"
cd "$OMPI_BUILD"

# Configure flags match the ones already baked into the installed
# $PREFIX/bin/ompi_info output. Do not change without also updating the README.
"$OMPI_SRC/configure" \
  --prefix="$PREFIX" \
  --with-cuda=/usr/local/cuda \
  --enable-mpi-fortran=usempi \
  --disable-mpi-cxx \
  --disable-oshmem \
  --without-verbs \
  --without-cma \
  CC=/usr/bin/gcc-11 \
  CXX=/usr/bin/g++-11 \
  FC=/usr/bin/gfortran-9 \
  2>&1 | tee -a "$LOG"

make -j "$MOOSE_JOBS" 2>&1 | tee -a "$LOG"
make install 2>&1 | tee -a "$LOG"

echo
echo "[build_openmpi] verifying CUDA support:"
"$PREFIX/bin/ompi_info" | grep -iE 'MPI extensions|Prefix' | head
