# Solid Mechanics

The MOOSE Solid Mechanics module provides finite-element models for deformation, stress,
inelastic material behavior, fracture, contact, and related mechanics applications. See the
[Solid Mechanics documentation](https://mooseframework.inl.gov/modules/solid_mechanics/index.html)
for usage details and examples.

## Prerequisites

Install MOOSE's compiler and dependency stack by following the
[MOOSE installation instructions](https://mooseframework.inl.gov/getting_started/installation/index.html),
then clone the MOOSE repository. The commands below assume the repository is located at
`<path_to_moose>`.

If MOOSE was installed with Conda, activate its environment before configuring, building, or
testing:

```bash
conda activate moose
```

## Build

For the standard CPU build:

```bash
cd <path_to_moose>/modules/solid_mechanics
make -j4
```

Replace `4` with the number of build jobs appropriate for your system. The resulting optimized
executable is `solid_mechanics-opt`.

### Build With Kokkos

Kokkos objects are not included in the standard configuration. Configure MOOSE from the
repository root before building the module:

```bash
cd <path_to_moose>
./configure --with-kokkos=cpu
cd modules/solid_mechanics
make -j4
```

Use a MOOSE dependency stack whose PETSc installation includes Kokkos support. After changing
the MOOSE configuration or updating dependencies, clean the framework before rebuilding to avoid
stale generated headers or library metadata:

```bash
cd <path_to_moose>/framework
make clean
cd ../modules/solid_mechanics
make -j4
```

Verify the resulting executable reports Kokkos:

```bash
./solid_mechanics-opt --show-capabilities \
  | python3 -c 'import json,sys; print("kokkos.value =", json.load(sys.stdin)["kokkos"]["value"])'
```

`kokkos.value` should be a Kokkos version string (e.g. `4.7.4`), not `false`. If it is `false`,
`./configure --with-kokkos=cpu` did not take effect (usually a missing `framework/make clean`).

### Build With Kokkos-CUDA

The conda `moose` env's PETSc does not include CUDA. A separately maintained from-scratch
stack (PETSc + libmesh + WASP + CUDA-aware OpenMPI) at `<path_to_moose>/kokkos-cuda-stack/`
provides a CUDA-enabled dependency set. Its `README.md` documents:

- how to build the stack (`kokkos-cuda-stack/scripts/all.sh`) and rebuild only MOOSE against it
  (`kokkos-cuda-stack/scripts/build_moose.sh`);
- the runtime environment (`PETSC_DIR`, `PETSC_ARCH=""`, `LIBMESH_DIR`, `WASP_DIR`) needed
  to run `solid_mechanics-opt` against that stack;
- the two-part GPU-KSP recipe (`jacobi + cg` in the input plus
  `PETSC_OPTIONS="-vec_type kokkos -mat_type aijkokkos"` in the environment);
- how to verify GPU dispatch with `nsys profile` (not `nvidia-smi`).

After building against that stack, `solid_mechanics-opt --show-capabilities` should report
both `kokkos.value` and `cuda.value` as version strings. Run with `--compute-device=cuda`.

## Test

Run the complete Solid Mechanics test suite from the module directory:

```bash
cd <path_to_moose>/modules/solid_mechanics
./run_tests -j4
```

For a CPU Kokkos build, select the CPU compute device:

```bash
./run_tests --compute-device=cpu -j4
```

To run only the Kokkos isotropic linear elasticity tests:

```bash
./run_tests --re=kokkos.linear_elasticity --compute-device=cpu -j2
```

## Run

Run an input file with the optimized executable:

```bash
./solid_mechanics-opt -i <input_file.i>
```

For a Kokkos input, also select the configured compute device:

```bash
./solid_mechanics-opt -i <input_file.i> --compute-device=cpu
```
