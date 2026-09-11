# MOOSE with Kokkos: build and run

This directory holds a from-scratch **CUDA-capable** MOOSE stack (PETSc +
libmesh + WASP + MOOSE) built against system compilers, plus this README, which
also documents the simpler **CPU-only** Kokkos build that uses the `moose`
conda environment.

Two supported Kokkos backends on this machine:

| Backend | Compilers   | PETSc + libmesh + WASP                         | Where the recipe lives |
|---------|-------------|------------------------------------------------|------------------------|
| CPU     | conda toolchain | `moose` conda env (`/home/chenghau.yang/miniforge/envs/moose`) | [CPU Kokkos build](#cpu-kokkos-build-conda-env) |
| CUDA    | `/usr/bin/gcc-11` + `/usr/local/cuda`      | `kokkos-cuda-stack/prefix/` (this directory)    | [CUDA Kokkos build](#cuda-kokkos-build-from-scratch-stack) |

The two paths are independent. Building one does not disturb the other. What
they share is a single MOOSE source tree; `./configure --with-kokkos=...` at
the MOOSE repo root selects which backend the next `make` compiles for, and
the recipe below shows how to rebuild MOOSE against either stack.

## Layout

The whole stack lives under `moose-kokkos/kokkos-cuda-stack/`. Scripts and docs
are git-tracked; the build artifacts (`prefix/`, `src/`, `logs/`) are
gitignored so a fresh clone of `moose-kokkos` gets everything needed to run
`scripts/all.sh` but not the ~1.4 GB of binaries it produces.

```
moose-kokkos/
  kokkos-cuda-stack/
    README.md          <-- this file                          [tracked]
    ATTEMPTS.md        <-- attempt-by-attempt log             [tracked]
    SUMMARY.md         <-- one-page summary                   [tracked]
    .gitignore         <-- ignores prefix/, src/, logs/       [tracked]
    scripts/                                                  [tracked]
      env.sh              <-- clean toolchain env; source, do not execute
      build_openmpi.sh    <-- CUDA-aware OpenMPI 4.1.6 into $PREFIX/bin (idempotent)
      build_petsc.sh      <-- PETSc 3.25.4 with Kokkos + CUDA
      build_libmesh.sh    <-- libmesh (all four methods: opt, oprof, devel, dbg)
      build_wasp.sh       <-- WASP (MOOSE's HIT parser dependency)
      build_moose.sh      <-- configure MOOSE, clean framework, rebuild framework + solid_mechanics
      run_benchmark.sh    <-- optional: run mesh-scaling benchmark
      all.sh              <-- run build_openmpi -> petsc -> libmesh -> wasp -> moose -> benchmark
    prefix/            <-- install prefix ($PREFIX)           [gitignored]
    src/               <-- OpenMPI source + build tree        [gitignored]
    logs/              <-- build logs                         [gitignored]
```

`env.sh` derives `STACK_DIR=$MOOSE_DIR/kokkos-cuda-stack` and
`PREFIX=$STACK_DIR/prefix`, so moving the MOOSE checkout requires changing
only `MOOSE_DIR` in `env.sh`.

## Prerequisites

Common to both backends:

- MOOSE source tree at `/home/chenghau.yang/packages/moose-kokkos` (this
  branch: `solid_mechanics_kokkos`).
- Enough disk for two PETSc arches: conda's `arch-moose` and the from-scratch
  `arch-scratch-cuda`. Each is ~2-4 GB after third-party downloads.

CPU-only backend also needs:

- Miniforge with the `moose` conda env populated by MOOSE's installer (PETSc
  with Kokkos-CPU already in the env).

CUDA backend also needs:

- System OpenMPI 4.x wrapping `gcc-11` / `g++-11` / `gfortran-9` at
  `/usr/bin/mpicc,mpicxx,mpif90`. (`OMPI_FC` must be set to `/usr/bin/gfortran-9`
  because this box has no bare `/usr/bin/gfortran` symlink.)
- CUDA 12.4 at `/usr/local/cuda` and an Ampere-class GPU. This box has two
  NVIDIA RTX A5000 (sm_86); Kokkos target is `Kokkos_ARCH_AMPERE86`.

---

## CPU Kokkos build (conda env)

Kokkos objects are not built into the default MOOSE configuration; you must
pass `--with-kokkos=cpu` to `./configure` and then clean the framework so the
generated headers pick up the new flag.

```bash
conda activate moose

cd /home/chenghau.yang/packages/moose-kokkos
./configure --with-kokkos=cpu

# stale generated headers/library metadata will confuse the framework build after
# a configure change; skipping `make clean` produces an executable that links
# but reports kokkos.value=false in --show-capabilities.
cd framework && make clean && make -j "$MOOSE_JOBS"

cd ../modules/solid_mechanics && make -j "$MOOSE_JOBS"
```

Verify the resulting `solid_mechanics-opt`:

```bash
./solid_mechanics-opt --show-capabilities \
  | python3 -c 'import json,sys; d=json.load(sys.stdin); print("kokkos.value =", d["kokkos"]["value"])'
```

`kokkos.value` should be a version string (e.g. `4.7.4`), not `false`. If it
is `false`, `./configure --with-kokkos=cpu` was skipped or the framework was
not cleaned after configure. A ~57 KB `solid_mechanics-opt` is normal - the
code lives in `lib/libsolid_mechanics-opt.so*`.

To run a Kokkos input on CPU:

```bash
./solid_mechanics-opt -i <input.i> --compute-device=cpu
```

---

## CUDA Kokkos build (from-scratch stack)

The CUDA path deliberately avoids conda. Conda's `mpicxx` wrapper injects
`-I<env>/include -L<env>/lib` at the *front* of every compilation, which
shadows PETSc's downloaded/built headers (Kokkos, Umpire, HDF5, ...) with
conda-installed copies of the same headers but with different feature flags.
Every earlier CUDA build attempt (see `ATTEMPTS.md`) failed for this reason.
System `/usr/bin/mpicxx` adds only `-I/usr/lib/x86_64-linux-gnu/openmpi/include`;
nothing shadows PETSc's `arch-*/externalpackages`.

### One-shot build

```bash
/home/chenghau.yang/packages/moose-kokkos/kokkos-cuda-stack/scripts/all.sh
```

`all.sh` sources `env.sh` and runs, in order: `build_openmpi.sh` (idempotent -
skips if `$PREFIX/bin/ompi_info` already reports CUDA support), `build_petsc.sh`,
`build_libmesh.sh`, `build_wasp.sh`, `build_moose.sh`, `run_benchmark.sh`.

`env.sh` erases every environment variable not on a small whitelist because
conda's activation stamps many build-hint vars (`build_alias`,
`CXX_FOR_BUILD`, `CMAKE_ARGS`, `CPP`, `GCC_RANLIB`, ...) that autoconf reads
even after `unset` of the obvious ones. The escape hatch
`MOOSE_CUDA_STACK_KEEP="VAR1 VAR2"` lets you preserve extra vars deliberately.

Approximate wall-clock on this box (32-core, `MOOSE_JOBS=8`): OpenMPI ~10 min
(first time only), PETSc ~40 min, libmesh ~20 min, WASP a few min, MOOSE
framework + solid_mechanics ~10 min.

### CUDA-aware OpenMPI

The system OpenMPI is not CUDA-aware:

```bash
/usr/bin/ompi_info | grep opal_built_with_cuda_support   # -> false
```

PETSc runs with `-mat_type aijkokkos -vec_type kokkos` will abort (MPI error
code 76) against a non-CUDA-aware MPI unless GPU-aware MPI is disabled via
`-use_gpu_aware_mpi 0`.

`build_openmpi.sh` (step 0 of `all.sh`) builds OpenMPI 4.1.6 with
`--with-cuda=/usr/local/cuda` into `$PREFIX/bin`. It is idempotent: on
re-runs it checks `$PREFIX/bin/ompi_info` and skips if CUDA support is
already reported. To force a rebuild, delete `$PREFIX/bin/ompi_info`.

`env.sh` prepends `$PREFIX/bin` to `PATH`, so `libmesh`, `WASP`, and
downstream `command -v mpicc` all pick up the CUDA-aware wrappers.
`build_petsc.sh` calls `$PREFIX/bin/{mpicc,mpicxx,mpif90}` by absolute path
and refuses to run if `$PREFIX/bin/mpicc` is missing.

Verify after `build_openmpi.sh`:

```bash
$PREFIX/bin/ompi_info | grep -iE 'MPI extensions|opal_built_with_cuda_support'
# expect: "cuda" listed in MPI extensions
```

### Build MOOSE against the stack (only)

If PETSc/libmesh/WASP are already installed at `$PREFIX` and you only want to
rebuild MOOSE (e.g. after pulling source changes):

```bash
/home/chenghau.yang/packages/moose-kokkos/kokkos-cuda-stack/scripts/build_moose.sh
```

`build_moose.sh` sources `env.sh`, sets `PETSC_DIR=$PREFIX`, `PETSC_ARCH=""`,
`LIBMESH_DIR=$PREFIX`, runs `./configure --with-kokkos=cuda` at the MOOSE
repo root, `make clean` in `framework/`, then `make -j "$MOOSE_JOBS"` in
`framework/` and `modules/solid_mechanics/`. It ends with a capability
check that prints `kokkos.value` and `cuda.value`.

Expected values after a successful CUDA build:

```
kokkos.value = 4.7.4      # or newer Kokkos version
cuda.value   = 12.4.0
```

### Switching back to the conda CPU stack

Open a **fresh shell** (do not source `env.sh`), then follow the
[CPU Kokkos build](#cpu-kokkos-build-conda-env) steps. `make clean` in
`framework/` after re-running `./configure --with-kokkos=cpu` is what makes
this work; the framework's generated headers depend on the configure flag.

---

## Running Kokkos inputs on GPU

### GPU runtime environment

The stack's PETSc lives outside `PETSC_ARCH`, so both must be set on every
GPU-run shell:

```bash
. /home/chenghau.yang/packages/moose-kokkos/kokkos-cuda-stack/scripts/env.sh
export PETSC_DIR=$PREFIX
export PETSC_ARCH=""
export LIBMESH_DIR=$PREFIX
export WASP_DIR=$PREFIX
```

Add on the MOOSE command line:

```bash
./solid_mechanics-opt -i <input.i> --compute-device=cuda
```

### Enabling GPU KSP (two-part recipe)

Assembly on the GPU is not enough by itself; the linear solve also has to run
in Kokkos memory to see the full speedup. That requires two things at once:

**Part 1: input file** - pick a PC/KSP that has a Kokkos path. The default
`ILU` from `SMP` will segfault on `aijkokkos` at the first triangular solve
(backtrace: `PCApply_ILU -> MatSolve_SeqAIJKokkos_LU -> VecGetKokkosView<Kokkos::CudaSpace>`).
Use `jacobi + cg` (or another combination whose components exist for
Kokkos matrices):

```
[Executioner]
  type = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -ksp_type'
  petsc_options_value = 'jacobi   cg'
[]
```

**Part 2: environment variable** - `-mat_type aijkokkos` and `-vec_type kokkos`
cannot go into `petsc_options_iname`. MOOSE rejects the non-prefixed
`-mat_type` with a solver-system-prefix check, and libmesh's primary
`PetscVector::init` calls `VecSetFromOptions` only for the non-prefixed
`-vec_type` at library level. Pass them via the PETSc-wide env var:

```bash
# CUDA-aware MPI available:
export PETSC_OPTIONS="-vec_type kokkos -mat_type aijkokkos"

# CUDA-aware MPI NOT available (system OpenMPI):
export PETSC_OPTIONS="-vec_type kokkos -mat_type aijkokkos -use_gpu_aware_mpi 0"

./solid_mechanics-opt -i <input.i> --compute-device=cuda
```

### Verifying the run actually used the GPU

Use `nsys profile`, not `nvidia-smi`:

```bash
nsys profile --stats=true ./solid_mechanics-opt -i <input.i> --compute-device=cuda 2>&1 \
  | grep -i 'Moose::Kokkos::'
```

`nvidia-smi` polls at ~1 s, but typical MOOSE Kokkos kernels are tens to
hundreds of microseconds each, so `nvidia-smi --query-compute-apps` will
usually see nothing even when the run is entirely on GPU. `nsys` records
every kernel launch and the `--stats=true` summary groups them by name.

---

## Where else to look

- `ATTEMPTS.md` in this directory - full chronology of what failed and why,
  including OpenBLAS/gfortran-path fights, `LIBMESH_BUILD_DIR` divergence
  between `configure_libmesh.sh` and `update_and_rebuild_libmesh.sh`, and the
  `nvidia-smi` false-negative that cost me an hour.
- `modules/solid_mechanics/README.md` - CPU Kokkos build steps in module
  form (this README is the superset).
