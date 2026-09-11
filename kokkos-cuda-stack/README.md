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
      activate.sh         <-- daily use / direnv; light-weight USE vars, no purge
      build_openmpi.sh    <-- CUDA-aware OpenMPI 4.1.6 into $PREFIX/bin (idempotent, self-heals moved installs)
      build_petsc.sh      <-- PETSc 3.25.4 with Kokkos + CUDA
      build_libmesh.sh    <-- libmesh (all four methods: opt, oprof, devel, dbg)
      build_wasp.sh       <-- WASP (MOOSE's HIT parser dependency)
      build_neml2.sh      <-- PyTorch (CUDA 12.4) + NEML2 into conda `moose` env (branch-specific)
      build_moose.sh      <-- configure MOOSE (--with-kokkos=cuda [--with-neml2]), rebuild framework + solid_mechanics
      run_benchmark.sh    <-- optional: run mesh-scaling benchmark
      all.sh              <-- init submodules -> openmpi -> petsc -> libmesh -> wasp -> [neml2] -> moose -> benchmark
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

`all.sh` sources `env.sh`, then runs (in order): `git submodule update --init`
for petsc/libmesh/wasp[+neml2]; `build_openmpi.sh` (idempotent - skips if
`$PREFIX/bin/ompi_info` already reports CUDA support at the current prefix);
`build_petsc.sh`; `build_libmesh.sh`; `build_wasp.sh`; `build_neml2.sh` (only
when NEML2 is on); `build_moose.sh`; `run_benchmark.sh`. NEML2 is on by
default; pass `--no-neml2` to skip it.

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

### Compile workflows

There are four ways to compile, from incremental to full-stack:

| # | Command | Rebuilds | Time | Use when |
|---|---|---|---|---|
| 1 | `. env.sh` + `make -j` in `framework/` or `modules/*/` | Only MOOSE source you edited | seconds-minutes | Everyday work: you edited a `.C`/`.h` |
| 1b | `. env.sh` + `make clean && make -j` in the same dir | Everything in that dir (framework or module) | ~10 min | After a `git pull` that touches many MOOSE files, or after re-running `./configure` |
| 2 | `scripts/build_moose.sh` | All of MOOSE + reruns `./configure --with-kokkos=cuda [--with-neml2]` + `make clean framework/` | ~10 min | You switched compute-device (`cpu` <-> `cuda`), or the Kokkos/NEML2 configure got out of sync |
| 2n | `scripts/build_neml2.sh` | Reinstalls PyTorch (if needed) + NEML2 into conda moose env | ~10-30 min | NEML2 source changed, or PyTorch missing/wrong-CUDA |
| 3 | `scripts/all.sh` (or `scripts/all.sh --no-neml2`) | Everything: submodules + OpenMPI + PETSc + libmesh + WASP + [NEML2] + MOOSE | ~1.5-2 h | Dependency version bump, dep install got corrupted, or first-time setup on a new machine |

Rule of thumb: default to (1). If MOOSE build fails weirdly, try (2). Reach
for (2n) if NEML2 source changed but the stack is otherwise fine. Reach for
(3) only when a dependency (PETSc / libmesh / WASP / OpenMPI) itself needs
to change.

Workflows (1) and (1b) work because `env.sh` exports `PETSC_DIR`,
`PETSC_ARCH=""`, `LIBMESH_DIR`, and `WASP_DIR` (all = `$PREFIX`); MOOSE's
Makefiles read those to locate the installed stack. Do not skip
`. env.sh` -- unset `LIBMESH_DIR` makes MOOSE fall back to
`$MOOSE_DIR/libmesh/installed/`, which is empty on this checkout, and the
build fails with `libmesh-config: not found`.

### NEML2 support

NEML2 provides the GPU-side material update path on this branch. It links
against PyTorch and installs as a Python package. Layout choice on this
box: PyTorch + NEML2 live in the existing `moose` **conda env**
(`/home/chenghau.yang/miniforge/envs/moose`); PETSc / libmesh / WASP /
OpenMPI / MOOSE stay in the from-scratch stack. `build_moose.sh` and
`build_neml2.sh` prepend the conda env `bin/` to `PATH` *after*
`$PREFIX/bin`, so `mpicxx` remains ours (CUDA-aware OpenMPI) while
`python3` comes from the conda env (has NEML2 installed).

Compilers for cmake: `build_neml2.sh` pins `CC=$PREFIX/bin/mpicc`,
`CXX=$PREFIX/bin/mpicxx`, `FC=$PREFIX/bin/mpif90` before invoking pip so
`libneml2.so` links against the same C++ ABI and OpenMPI as MOOSE. See
`build_neml2.sh` for the details.

Toggle: `all.sh` defaults to NEML2 on. Pass `--no-neml2` to skip it and
configure MOOSE without `--with-neml2`. `build_moose.sh` respects the
env var `NEML2_SUPPORT` (default `1`); `NEML2_SUPPORT=0 build_moose.sh`
builds without NEML2.

After (2) or (3), `build_moose.sh` prints a capability check. Expected:

```
kokkos.value = 4.7.4      # or newer Kokkos version
cuda.value   = 12.4.0
neml2.value  = <version>  # only if built --with-neml2
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
