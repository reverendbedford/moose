#!/bin/bash
# Step 1: build PETSc source in the in-tree `petsc/` under PETSC_ARCH=arch-scratch-cuda,
# using the CUDA-aware OpenMPI installed by build_openmpi.sh at $PREFIX/bin
# (wraps gcc-11 + CUDA 12.4). Leaves `arch-moose/` untouched.
set -e
set -o pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
. "$SCRIPT_DIR/env.sh"

LOG="$LOGS/petsc-$(date +%Y%m%d-%H%M%S).log"
echo "[build_petsc] logging to $LOG"

# env.sh sets PETSC_ARCH="" for USING the installed PETSc. During the build
# we need a non-empty arch so PETSc's Makefile has somewhere to place object
# files under $PETSC_SRC_DIR/<arch>. Use a local BUILD_ARCH; do not export.
BUILD_ARCH=arch-scratch-cuda

# Fresh arch dir every run.
rm -rf "$PETSC_SRC_DIR/$BUILD_ARCH"

cd "$PETSC_SRC_DIR"

# We invoke PETSc's configure directly (not update_and_rebuild_petsc.sh) so we
# control every flag. The list mirrors scripts/configure_petsc.sh with:
#   - explicit compilers (mpicc/mpicxx/mpif90)
#   - --prefix pointing to STACK_DIR/prefix so nothing lands in the source tree
#   - CUDA options (arch=86, cuda-dir under /usr/local/cuda)
#   - --download-libceed=0 (was the killer link failure in the conda attempt)
#   - --download-hdf5=1 (no MPI HDF5 in system apt packages)
# Downstream Fortran packages (MUMPS, SCALAPACK, STRUMPACK, HYPRE-fortran)
# need libgfortran.so which lives in /usr/lib/gcc/x86_64-linux-gnu/9/ on this box
# (only libgfortran.so.5 is on the default search path). Bake it into LDFLAGS
# so package makefiles that link with the C driver still resolve gfortran.
export LDFLAGS="-L/usr/lib/gcc/x86_64-linux-gnu/9 ${LDFLAGS:-}"

# Sanity: refuse to build against system MPI. We want the CUDA-aware wrappers
# from $PREFIX/bin so PETSc's Kokkos backend can use GPU-aware MPI without
# -use_gpu_aware_mpi 0. Run scripts/build_openmpi.sh first if this trips.
if [ ! -x "$PREFIX/bin/mpicc" ]; then
  echo "[build_petsc] ERROR: $PREFIX/bin/mpicc missing. Run scripts/build_openmpi.sh first." >&2
  exit 1
fi

python3 ./configure \
  --prefix="$PREFIX" \
  LDFLAGS="$LDFLAGS" \
  --with-cc="$PREFIX/bin/mpicc" \
  --with-cxx="$PREFIX/bin/mpicxx" \
  --with-fc="$PREFIX/bin/mpif90" \
  --with-64-bit-indices \
  --with-cxx-dialect=C++17 \
  --ignoreCxxBoundCheck=1 \
  --with-debugging=no \
  --with-fortran-bindings=0 \
  --with-mpi=1 \
  --with-openmp=1 \
  --with-strict-petscerrorcode=1 \
  --with-shared-libraries=1 \
  --with-sowing=0 \
  --with-x=0 \
  --with-ssl=0 \
  --with-cuda=1 \
  --with-cuda-arch=86 \
  --with-cudac=/usr/local/cuda/bin/nvcc \
  --with-cuda-dir=/usr/local/cuda \
  --with-blas-lib=/usr/lib/x86_64-linux-gnu/libblas.so \
  --with-lapack-lib=/usr/lib/x86_64-linux-gnu/liblapack.so \
  --download-hpddm=1 \
  --download-hypre=1 \
  --download-metis=1 \
  --download-mumps=1 \
  --download-ptscotch=1 \
  --download-parmetis=1 \
  --download-scalapack=1 \
  --download-slepc=1 \
  --download-strumpack=1 \
  --download-superlu_dist=1 \
  --download-kokkos=1 \
  --download-kokkos-commit=4.7.04 \
  --download-kokkos-kernels=1 \
  --download-kokkos-kernels-commit=4.7.04 \
  --download-umpire \
  --download-hdf5=1 \
  --with-hdf5-fortran-bindings=0 \
  --download-zlib=1 \
  --with-libceed=0 \
  --with-make-np="$MOOSE_JOBS" \
  2>&1 | tee "$LOG"

# PETSc after successful configure suggests `make PETSC_DIR=... PETSC_ARCH=... all`.
make PETSC_DIR="$PETSC_SRC_DIR" PETSC_ARCH="$BUILD_ARCH" all 2>&1 | tee -a "$LOG"
make PETSC_DIR="$PETSC_SRC_DIR" PETSC_ARCH="$BUILD_ARCH" install 2>&1 | tee -a "$LOG"

# Quick sanity: what did we get?
echo
echo "[build_petsc] Installed PETSc capabilities:"
grep -E "PETSC_HAVE_CUDA |PETSC_HAVE_KOKKOS|PETSC_PKG_CUDA_MIN_ARCH" \
    "$PREFIX/include/petscconf.h" | head -6
echo
ls -1 "$PREFIX/lib/libpetsc"*.so* "$PREFIX/lib/libkokkos"*.so* 2>/dev/null | head
