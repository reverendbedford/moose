#!/bin/bash
# Step 0: build CUDA-aware OpenMPI 4.1.6 into $PREFIX.
#
# Rationale: system /usr/bin/mpicc wraps an OpenMPI that reports
# `opal_built_with_cuda_support false`. PETSc with -mat_type aijkokkos +
# -vec_type kokkos then aborts (MPI error 76) unless -use_gpu_aware_mpi 0 is
# passed. Building our own CUDA-aware OpenMPI at $PREFIX lets PETSc, libmesh
# and MOOSE all link against it and lets GPU-aware MPI stay on.
#
# Idempotency: if $PREFIX/bin/ompi_info exists AND reports CUDA support AND
# the install was built at (or has been repaired to) the current $PREFIX,
# skip. Delete $PREFIX/bin/ompi_info to force a rebuild.
#
# Stale-prefix self-repair: if $PREFIX was physically moved after OpenMPI was
# installed, opal_wrapper's RUNPATH (an ELF header, not touchable by sed)
# still points at the original prefix; libmpi.la and pkg-config .pc files
# reference it too. Detect via readelf on opal_wrapper's RUNPATH, sed-fix all
# text files under $PREFIX, and force a rebuild so the wrapper binary bakes
# the new prefix. Without this, downstream libtool builds (libmesh contrib)
# fail with "libopen-rte.la is not a valid libtool archive" pointing at the
# vanished old path.
set -e
set -o pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
. "$SCRIPT_DIR/env.sh"

LOG="$LOGS/openmpi-$(date +%Y%m%d-%H%M%S).log"

detect_stale_prefix() {
  # Read opal_wrapper's DT_RUNPATH and extract the first ".../lib" entry as
  # the baked prefix. Empty output = not stale (or wrapper missing).
  [ -x "$PREFIX/bin/opal_wrapper" ] || return 0
  command -v readelf >/dev/null || return 0
  local baked
  baked=$(readelf -d "$PREFIX/bin/opal_wrapper" 2>/dev/null \
          | awk '/R.NPATH/ { gsub(/.*\[|\].*/, ""); print }' \
          | tr ':' '\n' | grep '/lib$' | head -1 | sed 's|/lib$||')
  if [ -n "$baked" ] && [ "$baked" != "$PREFIX" ]; then
    printf '%s\n' "$baked"
  fi
}

if [ -x "$PREFIX/bin/ompi_info" ] \
   && "$PREFIX/bin/ompi_info" | grep -q "opal_built_with_cuda_support:true\|MPI extensions:.*cuda"; then
  stale=$(detect_stale_prefix)
  if [ -z "$stale" ]; then
    echo "[build_openmpi] $PREFIX/bin/ompi_info already reports CUDA support; skipping."
    echo "[build_openmpi] delete $PREFIX/bin/ompi_info to force a rebuild."
    exit 0
  fi
  echo "[build_openmpi] stale prefix detected: opal_wrapper baked at $stale, current \$PREFIX is $PREFIX."
  echo "[build_openmpi] repairing installed text files (.la, .pc, wrapper-data, ...) via sed:"
  n=$(grep -rlI "$stale" "$PREFIX" 2>/dev/null | tee /dev/stderr | wc -l)
  grep -rlI "$stale" "$PREFIX" 2>/dev/null | xargs -r sed -i "s|$stale|$PREFIX|g"
  echo "[build_openmpi] $n text files repaired; forcing OpenMPI rebuild to bake the new prefix into opal_wrapper."
  # Fall through to the rebuild below.
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
