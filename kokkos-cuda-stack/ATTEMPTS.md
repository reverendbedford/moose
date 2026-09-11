# From-scratch build log

Chronological journal of what was tried, in what environment, with which result. Each entry: what was run → what happened → why → next move. Prior conda attempts are summarised at the top; new "from scratch" attempts append below.

## Prior conda-based attempts (2026-08-25) — all failed

Full write-up lives in `~/.claude/projects/-home-chenghau-yang-packages-moose-kokkos/memory/project_petsc_kokkos_cpu_only.md`. In one paragraph:

1. **Source PETSc rebuild inside the `moose` conda env** (`PETSC_ARCH=arch-moose-cuda`, keeping `arch-moose` intact). Fought conda's `mpicxx` header-shadowing: fixed Kokkos LAMBDA test with `NVCC_APPEND_FLAGS='-ccbin /usr/bin/g++'`, then Umpire's Allocator.hpp was shadowed by conda's older umpire, then HDF5, etc. Every fix uncovered the next package. Abandoned.
2. **Fresh `moose-cuda` conda env with pre-built `petsc=3.25.4=cuda12_real_hd4db405_0`**: that PETSc has `PETSC_HAVE_CUDA 1` but no `PETSC_HAVE_KOKKOS`. MOOSE's `framework/kokkos.mk` gates on `PETSC_HAVE_KOKKOS`, so the CUDA backend would never turn on. conda-forge only ships `kokkos-cuda12*` up to 4.6.00; MOOSE pins 4.7.04.
3. **`moose-cuda` conda env, petsc/slepc removed, PETSc source rebuild**: got past Kokkos+Umpire+CUDA-LAMBDA; then libceed failed to link (downstream libstdc++ propagation), then HDF5 failed (libz.so.1 not in rpath). Each package's link line uses conda `mpicc`/`mpicxx` wrappers which add `-Wl,-rpath,<env>/lib` selectively and miss dependencies. Abandoned.

## From-scratch attempt (2026-08-26) — in progress

### Environment survey

```
gcc/g++ 11.4.0 (Ubuntu 22.04, /usr/bin/gcc-11, /usr/bin/g++-11)
gfortran 9.5.0 (/usr/bin/gfortran-9 — /usr/bin/gfortran is NOT installed, need OMPI_FC=/usr/bin/gfortran-9)
mpicc / mpicxx / mpif90 — /usr/bin (OpenMPI, wraps system gcc/g++/gfortran)
nvcc 12.4 (/usr/local/cuda/bin/nvcc)
CUDA runtime libs: /usr/lib/x86_64-linux-gnu/libcuda.so.1, libnvidia-ml.so.1
CUDA stubs: /usr/local/cuda/lib64/stubs/libcuda.so, libnvidia-ml.so
Autotools/cmake/perl/python3: all in /usr/bin
libtool binary NOT in PATH (only /usr/bin/libtoolize) — Debian's `libtool` package generates `libtool` per-project via autoreconf, which is fine for PETSc downloads.
zlib1g-dev + libhwloc-dev + libhdf5-dev(non-MPI) present in system.
Disk: 320 TB free on NFS home. RAM: 500 GB. 32 CPU cores.
```

### Design decisions

- **Prefix**: `/home/chenghau.yang/packages/moose-cuda-stack/prefix` (PETSc + libmesh install here; MOOSE stays in-tree).
- **PETSc source**: reuse in-tree `/home/chenghau.yang/packages/moose-kokkos/petsc/` with `PETSC_ARCH=arch-scratch-cuda`. Leaves the existing `arch-moose/` untouched. Uses `--prefix=<prefix>` so the installed layout is a clean sibling of the conda one.
- **libmesh source**: reuse in-tree `/home/chenghau.yang/packages/moose-kokkos/libmesh/`, install to `<prefix>` via `LIBMESH_DIR=<prefix>` + custom `LIBMESH_BUILD_DIR`.
- **MOOSE**: run `./configure --with-kokkos=cuda` with `PETSC_DIR=<prefix>` and `LIBMESH_DIR=<prefix>` set. `make clean` in framework then rebuild.
- **Fortran**: `OMPI_FC=/usr/bin/gfortran-9` — needed because system `/usr/bin/gfortran` symlink is absent.
- **HDF5**: let PETSc `--download-hdf5=1`. The system's serial HDF5 won't work for parallel PETSc.
- **libceed**: DISABLED (`--download-libceed=0 --with-libceed=0`). It's optional for MOOSE and was the specific package that link-failed in the last conda-env attempt. Add back later if needed.

### Attempt 1 — 2026-08-26 00:23

Ran `scripts/build_petsc.sh` with the full download list (openblas, hypre, mumps, kokkos, kokkos-kernels, umpire, hdf5, zlib, hpddm, metis, mumps, ptscotch, parmetis, scalapack, slepc, strumpack, superlu_dist, libceed=off).

**Result:** configure passed; PETSc started downloading and building externalpackages; failed at OpenBLAS.

**Cause:** OpenBLAS shared-library link step:
```
/usr/bin/ld: cannot find -lgfortran: No such file or directory
```
System has `libgfortran-9-dev` which puts `libgfortran.so` at `/usr/lib/gcc/x86_64-linux-gnu/9/libgfortran.so`, but that path is NOT on the default library search path. gcc-11 (which mpicc wraps here) searches `/usr/lib/gcc/x86_64-linux-gnu/11/`, which has no gfortran runtime (system has no gfortran-11). Only `libgfortran.so.5` is at the standard `/usr/lib/x86_64-linux-gnu/`.

**Fix:** replace `--download-openblas` with system BLAS/LAPACK (`libblas.so` + `liblapack.so` are already installed via `libblas-dev` + `liblapack-dev`). This sidesteps the gfortran linkage problem for BLAS entirely; downstream Fortran packages (MUMPS, SCALAPACK, STRUMPACK) still need gfortran but as first-party linkers they add `-L$(gfortran -print-file-name=libgfortran.so)` themselves. Also exported `LDFLAGS=-L/usr/lib/gcc/x86_64-linux-gnu/9` and passed it through to PETSc's configure so downstream packages see the gfortran path too.

### Attempt 2 — 2026-08-26 00:45

Same script with the two fixes above. **Succeeded.**

- Full configure + download + build + install took ~43 minutes with `MOOSE_JOBS=8` on 32-core box.
- `prefix/include/petscconf.h` shows all four flags MOOSE needs:
  ```
  #define PETSC_HAVE_CUDA 1
  #define PETSC_HAVE_KOKKOS 1
  #define PETSC_HAVE_KOKKOS_KERNELS 1
  #define PETSC_PKG_CUDA_MIN_ARCH 86
  ```
  No `PETSC_HAVE_KOKKOS_WITHOUT_GPU` — this build IS Kokkos-with-GPU.
- `prefix/lib/` has `libpetsc.so.3.25.4`, `libkokkoscore.so.4.7.4`, `libkokkoskernels.so.4.7.4`, `libkokkoscontainers`, `libkokkosalgorithms`, and every other downloaded package.

### Attempt 3 — libmesh, first try (2026-08-26 01:29)

Ran `scripts/build_libmesh.sh` with `PETSC_DIR=$PREFIX PETSC_ARCH="" LIBMESH_DIR=$PREFIX LIBMESH_BUILD_DIR=$STACK_DIR/libmesh-build`.

**Result:** immediate failure: `configure_libmesh.sh: line 45: cd: /home/chenghau.yang/packages/moose-kokkos/libmesh/build: No such file or directory`, followed by `make: *** No targets specified`.

**Cause:** MOOSE ships two libmesh scripts:
- `scripts/update_and_rebuild_libmesh.sh` — respects `LIBMESH_BUILD_DIR` (cd's into it after `mkdir -p`)
- `scripts/configure_libmesh.sh` — hardcodes `cd "${SRC_DIR}/build"` at line 45, ignoring `LIBMESH_BUILD_DIR`

Their contract is inconsistent. If you override `LIBMESH_BUILD_DIR`, the rebuild script cd's to your dir, then configure_libmesh cd's *away* to `$SRC/build`, and `make` runs in a dir that has no configure output.

**Fix:** don't override `LIBMESH_BUILD_DIR`. Let the default `$LIBMESH_SRC_DIR/build` be used. The install still goes to `$PREFIX` via `LIBMESH_DIR`.

### Attempt 4 — libmesh, second try (2026-08-26 01:39)

Same script with `unset LIBMESH_BUILD_DIR` fix.

**Result:** libmesh configure died with `checking whether the C compiler works... no`, `configure: error: C compiler cannot create executables`. config.log showed:
- `build system type: x86_64-conda-linux-gnu`
- `found /home/chenghau.yang/miniforge/envs/moose/bin/mpicxx` (a conda binary — despite PATH not containing conda)
- `x86_64-conda-linux-gnu-c++: command not found` (conda's mpicxx trying to invoke its wrapped compiler which is no longer in PATH)

**Cause:** conda's activation stamps several env vars beyond `CONDA_*` that survive an ordinary `unset` list:
- `build_alias=x86_64-conda-linux-gnu` — autoconf reads this and treats the system as `x86_64-conda-linux-gnu`, then looks for `x86_64-conda-linux-gnu-mpicc`
- `CMAKE_ARGS` — points at conda's cross-compile tools
- `CPP=x86_64-conda-linux-gnu-cpp`, `CXX_FOR_BUILD`, `CPP_FOR_BUILD`, `GCC_RANLIB`, …
- `CONDA_TOOLCHAIN_BUILD`, `CONDA_TOOLCHAIN_HOST`
- `GSETTINGS_SCHEMA_DIR_CONDA_BACKUP` (harmless, but symptomatic)

An unset-by-name list is fragile; new conda releases add more vars. Rewrote `env.sh` to purge everything except a small `KEEP_VARS` safelist (`HOME`, `USER`, `PATH` — set explicitly after — `SSH_*`, `LANG`, `TERM`, …). This is a whitelist, so any new conda variable is dropped automatically.

### Attempt 5 — libmesh, third try (2026-08-26 01:39). **Succeeded.**

Same script, whitelisted env. libmesh built all 4 methods (opt, oprof, devel, dbg) into `prefix/lib/libmesh_*.so.0.0.0`. Took ~90 min single-configure across 4 methods, so we get one full-featured library.

### Attempt 6 — MOOSE reconfigure + rebuild (2026-08-26 03:15)

Ran `scripts/build_moose.sh`. Configure passed and framework began compiling; halted at:
```
***ERROR***
WASP does not seem to be available.
Make sure to either run scripts/update_and_rebuild_wasp.sh in your MOOSE directory,
or set WASP_DIR to a valid WASP install
make: *** [/home/chenghau.yang/packages/moose-kokkos/framework/moose.mk:481: wasp_submodule_status] Error 1
```

**Cause:** conda's `moose-wasp` package supplied WASP libs in the `moose` env, and MOOSE's `moose.mk` picks them up via `WASP_DIR=$CONDA_PREFIX` at configure time. Our from-scratch env has no such package. MOOSE ships `framework/contrib/wasp/` as a submodule + `scripts/update_and_rebuild_wasp.sh` for exactly this case.

**Fix:** added `scripts/build_wasp.sh` (calls `update_and_rebuild_wasp.sh` with `WASP_PREFIX=$PREFIX`) and exported `WASP_DIR=$PREFIX` in `env.sh`. Wired into `all.sh` between libmesh and MOOSE.

### Attempt 7 — WASP build (2026-08-26 03:16). **Succeeded.**

`libwaspcore.so.4.4.1`, `libwaspddi`, `libwaspexpr`, `libwasphalite`, etc. installed to `prefix/lib/`.

### Attempt 8 — MOOSE reconfigure + rebuild, retry (2026-08-26 03:16). **Succeeded.**

Same script, WASP now available. Framework and solid_mechanics module both built. `nvlink` printed some `Skipping incompatible '/usr/lib/x86_64-linux-gnu/libpthread.a'` warnings — harmless; those are static archives in a lib dir that also has the .so, and the linker picks the .so.

Final capability check on `solid_mechanics-opt`:
```
kokkos.value = 4.7.4
cuda.value   = 12.4.0
```

Both non-false → MOOSE Kokkos-CUDA backend is enabled.

### Attempt 9 — smoke test, first try (2026-08-26 03:20)

Ran `scripts/run_benchmark.sh`. The nvidia-smi-based smoke gate failed: `nvidia-smi saw MOOSE process on GPU during smoke run: 0`. Script aborted before the real benchmark.

**Cause:** nvidia-smi's default polling interval (~1 s) is coarser than MOOSE's smoke-run wall time (~0.4 s including MPI init). Independent probes proved Kokkos + CUDA are working end-to-end:
- A standalone Kokkos test program (`/tmp/kokkos_big.cpp`, 50 iterations of a 100M-element `RangePolicy<Cuda>` parallel_for) reports `Default Device: Cuda`, `Kokkos::Cuda[0] NVIDIA RTX A5000 : Selected`, and executes in 0.3 s wall time on GPU — but never appears in `nvidia-smi` either (below polling resolution).
- Explicit `Kokkos::View<double*, Kokkos::CudaSpace>("v", 200000000)` with a 15 s sleep DID show up in nvidia-smi: `1732 MiB used` on the RTX A5000.
- `nsys profile --stats=true` on the MOOSE smoke run shows `Moose::Kokkos::` CUDA kernels launched (e.g. `Moose::Kokkos::Dispatcher<...KokkosStressDivergence...ResidualLoop>`, `Moose::Kokkos::FE...`), each in the tens-of-microseconds range. Total GPU kernel time is a few hundred microseconds — much less than nvidia-smi's minimum-detectable duration.

**Fix:** replaced the nvidia-smi check with `nsys profile --stats=true` and a grep for `Moose::Kokkos::` in its output. That is conclusive regardless of kernel duration.

### Attempt 10 — smoke test + full benchmark (2026-08-26 03:26). **Succeeded.**

- nsys detected many `Moose::Kokkos::` CUDA kernels in the smoke run — proof MOOSE dispatched to GPU.
- Full benchmark ran: 5 reps × 5 CPU thread counts (`2 MPI ranks × {1,2,4,8,16} threads/rank`) + 5 reps × 1 GPU config (`1 MPI rank, 1 thread`), mesh 256×256.

Timings (mean wall time, from `cpu_and_gpu.csv`):
| Config | Wall time (s) |
|---|---|
| CPU 2×1 (2 threads) | 6.16 ± 0.29 |
| CPU 2×2 (4 threads) | 6.32 ± 0.13 |
| CPU 2×4 (8 threads) | 6.24 ± 0.30 |
| CPU 2×8 (16 threads) | 6.39 ± 0.12 |
| CPU 2×16 (32 threads) | 6.27 ± 0.34 |
| GPU 1×1 | 7.85 ± 0.04 |

GPU is ~1.25× SLOWER than the fastest CPU config at this problem size.

### Attempt 11 — bigger mesh (2026-08-26 04:15). Ran mesh 512×512, 3 reps + 1 warmup, CSV/PNG in `logs/cpu_and_gpu_512.csv`.

- CPU best (2 MPI × 4 threads): ~29-42 s (very noisy — the machine is shared)
- CPU worst (2 MPI × 16 threads): ~43-45 s
- GPU: 51.3 ± 0.1 s

GPU is still slower. Rationale (worth documenting so you don't chase a build bug):

The `kokkos_material_linear_elasticity.i` input assembles residual + Jacobian in Kokkos kernels (which our nsys smoke test confirmed run on GPU), but the KSP solve is dispatched through PETSc. This PETSc build has CUDA + Kokkos backends available (see `-lkokkoskernels`, `-lcublas`, etc. in the link line), yet the input doesn't set `-mat_type aijkokkos -vec_type kokkos -pc_type ...` or any GPU-preferring solver options. So the linear solve stays on host CPU — only assembly runs on GPU. Assembly is memory-bound; kernel launch + host-device sync for each Newton assembly loop swamps the small on-device compute at this problem size.

**To see GPU pull ahead you need one or both:**
1. Move the KSP solve to GPU: append `-mat_type aijkokkos -vec_type kokkos -pc_type jacobi -ksp_type cg` (or similar) to `--extra-args` in `benchmark_cpu_and_gpu.py`.
2. Grow the mesh to O(10⁶) elements per rank so on-device compute per assembly loop dwarfs host↔device sync overhead.

The stack is correct. This is a MOOSE-input configuration question, not a build one.

### Attempt 12 — try -mat_type aijkokkos (2026-08-26 04:32)

Tested option 1 directly. Result:
```
[0]PETSC ERROR: For MVAPICH2-GDR, you need to set MV2_USE_CUDA=1 ...
[0]PETSC ERROR: For Cray-MPICH, export MPICH_GPU_SUPPORT_ENABLED=1 (see its 'man mpi'); for MPICH, export MPIR_CVAR_ENABLE_GPU=1
MPI_ABORT was invoked ... with errorcode 76.
```

`ompi_info` on the system OpenMPI: `opal_built_with_cuda_support:value:false`. **System OpenMPI is NOT CUDA-aware.** PETSc-with-Kokkos-CUDA requires CUDA-aware MPI when moving vectors/matrices to GPU (so it can pass GPU pointers between ranks without staging through host).

To unlock this path you'd need one of:
- Rebuild OpenMPI with `--with-cuda=/usr/local/cuda`, then rebuild PETSc + libmesh + MOOSE against the new OpenMPI. Big rebuild.
- Switch to `--download-mpich` (or a system CUDA-aware MPICH) in the PETSc configure. Also a big rebuild.
- Stay with GPU-assembly + CPU-KSP as we have now. Fine for cases where assembly dominates (nonlinear materials, many quadrature points per element, complex constitutive laws); poor for linear elasticity where the KSP solve is the bottleneck.

Not attempted in this session — it would double the build time and the user gets a full working GPU-Kokkos MOOSE either way.

### Attempt 13 — run_tests smoke check (2026-08-26 04:35)

Tried `./run_tests --re=kokkos.linear_elasticity --compute-device=cuda -j2` under the from-scratch env. Failed at import: `ModuleNotFoundError: No module named 'pandas'`. System python3 (3.10) has no pandas; conda's `moose` env has it.

**Impact:** the `solid_mechanics-opt` executable itself is fine (the benchmark uses it directly via `mpirun`). Only MOOSE's TestHarness needs pandas. If you want to run the test suite:
- Simplest: `pip install --user pandas packaging pyyaml jinja2 numpy` and rerun. Or `pip install -r <MOOSE_DIR>/python/requirements.txt` if it exists.
- Or: run the test suite under the `moose` conda env (which has all Python deps), but pointed at the from-scratch build by exporting `PATH=$MOOSE_DIR/modules/solid_mechanics:$PATH` and `LD_LIBRARY_PATH=$STACK_DIR/prefix/lib`. Tricky because the conda env's PETSc/libmesh will fight ours.
- Or: don't run TestHarness — `mpirun -np N solid_mechanics-opt -i <input> --compute-device=cuda` works fine.

Not a build regression, just a Python-env note.

## Stack summary (what you can trust from here on)

- `moose` conda env: **unchanged**. Still CPU-only Kokkos, exactly as before.
- `moose-cuda-stack/prefix/`: PETSc 3.25.4 + libmesh + WASP built from source with **system OpenMPI + gcc-11 + CUDA 12.4** (no conda in the toolchain).
- `moose-kokkos/framework/` and `moose-kokkos/modules/solid_mechanics/`: MOOSE reconfigured with `--with-kokkos=cuda`, framework `make clean`-then-rebuilt, module rebuilt.
- `solid_mechanics-opt --show-capabilities`: `kokkos.value = 4.7.4`, `cuda.value = 12.4.0`.
- Kokkos runtime uses `Kokkos::Cuda` as default execution space, kernel launches confirmed via nsys.

## To rerun

```bash
cd /home/chenghau.yang/packages/moose-cuda-stack
# Everything at once (petsc → libmesh → wasp → moose → benchmark)
./scripts/all.sh

# Or individually
./scripts/build_petsc.sh          # ~45 min
./scripts/build_libmesh.sh        # ~90 min (4 methods)
./scripts/build_wasp.sh           # ~2 min
./scripts/build_moose.sh          # ~15 min
./scripts/run_benchmark.sh        # ~15 min (mesh 256)
```

To use the CUDA MOOSE from a shell later:
```bash
. /home/chenghau.yang/packages/moose-cuda-stack/scripts/env.sh
export PETSC_DIR=$PREFIX PETSC_ARCH="" LIBMESH_DIR=$PREFIX WASP_DIR=$PREFIX
$MOOSE_DIR/modules/solid_mechanics/solid_mechanics-opt \
  -i <input.i> --compute-device=cuda
```

To go back to the conda `moose` env's CPU build, just start a fresh shell (don't source `env.sh`) and `conda activate moose`. `./configure && make -j` in framework will pick the conda PETSc/libmesh again and rebuild.

## Files touched outside the stack directory

- `moose-kokkos/framework/libmoose-opt.so*`, `libmoose_kokkos-opt.so*` — rebuilt against the from-scratch stack. Delete `moose-kokkos/framework/build/` and `make clean && make` under conda's `moose` env to go back.
- `moose-kokkos/modules/solid_mechanics/solid_mechanics-opt`, `lib/lib*-opt.so*` — same.
- `moose-kokkos/petsc/arch-scratch-cuda/` — new arch dir; `arch-moose/` is untouched.
- `moose-kokkos/libmesh/build/` — libmesh build tree (used by `configure_libmesh.sh` per attempt #4). Safe to `rm -rf`.
- `moose-kokkos/framework/contrib/wasp/build/` — WASP build tree (safe to `rm -rf`).
- `moose-kokkos/modules/solid_mechanics/test/tests/kokkos/linear_elasticity/cpu_and_gpu.{csv,png}` — updated with real GPU numbers.





