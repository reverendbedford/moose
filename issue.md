# Plasticity benchmark regression with non-preset Kokkos loading

Branch: `neml2-v3-kokkos-solid-mechanics`

Observed: 2026-09-14 on tree at `34bb5b9855`, built through
`kokkos-cuda-stack/scripts/activate.sh`.

## Scope correction

The problem has only been observed in the staged benchmark under
`modules/solid_mechanics/test/tests/kokkos/plasticity/`.

The existing tests under
`modules/solid_mechanics/test/tests/kokkos/linear_elasticity/` are not affected.
In particular, that directory already covers:

- three-dimensional `KokkosIsotropicElasticity`;
- three-dimensional `KokkosStressDivergence` with Kokkos material properties;
- three-dimensional NEML2 stress and tangent transfer; and
- transient NEML2 perfect plasticity with managed state.

Those passing cases rule out a general failure of 3D Kokkos stress-divergence
assembly, NEML2-to-Kokkos material transfer, or managed plasticity state. The
previous conclusion that the common `KernelGrad` residual path was broken was
therefore too broad.

## Benchmark symptom

The staged benchmark compares the same three-dimensional perfect-plasticity
problem while moving NEML2, assembly, and PETSc between CPU and GPU paths:

| Step | NEML2 | Assembly | PETSc | Observed `ux_center` at `t=0.005` |
| --- | --- | --- | --- | --- |
| 1 | CPU | CPU | CPU | approximately `2.5e-3` |
| 2 | GPU | CPU | CPU | approximately `2.5e-3` |
| 3a | CPU | Kokkos | CPU | `-8.208e-05` |
| 3b | GPU | Kokkos | CPU | `-8.208e-05` |
| 4 | GPU | Kokkos | AIJKokkos | `+4.308e-05` |

The expected center displacement is approximately `2.5e-3`; the controlled
right boundary reaches approximately `5e-3`.

## Most relevant difference

The benchmark's Kokkos configurations had this setting on the controlled
right-face boundary:

```text
preset = false
```

That setting was introduced by commit `34bb5b9855` in Step 3a, Step 3b, and
Step 4. The passing 3D Kokkos inputs in `linear_elasticity/`, including
`neml2_kokkos_perfect_plasticity.i`, use the default `preset = true` with the
same `RealFunctionControl` loading pattern.

This is a substantially narrower difference than the Kokkos stress-divergence
assembly itself. With `preset = false`, the controlled displacement enters only
through the nodal BC residual. With the default preset behavior, it is also
written into the solution before element strain, stress, and residual assembly.
That distinction is especially relevant because `TorchFEInterpolation` builds
the NEML2 strain from the solution field.

Other benchmark-only differences that should remain in the A/B matrix are:

- mesh size (`N=16` by default versus `N=2` in the regression tests);
- explicit GAMG/GMRES options; and
- explicit CPU or Kokkos PETSc vector and matrix types in the runners.

## Current change

Remove `preset = false` from Step 3a, Step 3b, and Step 4 so these benchmarks use
the same controlled Kokkos Dirichlet loading path as the passing 3D Kokkos
plasticity test.

This is not yet proof of the root cause. It is the smallest configuration
change supported by the comparison with the passing test suite.

## Verification

Run the existing regression tests first:

```bash
cd modules/solid_mechanics
PYTHONPATH=$MOOSE_DIR/python:$PYTHONPATH ./run_tests -i tests \
  --re 'kokkos/linear_elasticity\.(material_3d|neml2_material_3d|neml2_perfect_plasticity)$'
```

Then run the coarse benchmark comparison:

```bash
cd modules/solid_mechanics/test/tests/kokkos/plasticity
MESH_N=2 ./run_coarse_exodus.sh
```

At `t=0.005`, verify for all five steps:

- `ux_right` is approximately `0.005`;
- `ux_center` is approximately `0.0025`;
- nonlinear and linear solves do not report divergence; and
- the Exodus displacement fields agree across the five configurations.

If Step 3/4 still fail, test one difference at a time while keeping the passing
`neml2_kokkos_perfect_plasticity.i` input as the baseline:

1. Add the benchmark postprocessors and output settings.
2. Add the benchmark PETSc solver options.
3. Increase `N` from 2 to 4, 8, and 16.
4. Select CPU NEML2 on CUDA Kokkos assembly for the Step 3a transfer path.
5. Select Kokkos PETSc vectors and AIJKokkos for Step 4.

## CUDA follow-up

The local conda environment uses a CPU-only Torch build, so the CUDA paths still
need to be run on a host with CUDA-enabled Torch and NEML2. Run both the 3D
benchmark and the new `plasticity_2d/run_coarse_exodus.sh` matrix, checking:

- Step 2 isolates GPU NEML2 with CPU assembly and PETSc;
- Step 3a isolates CUDA Kokkos assembly with CPU NEML2 and PETSc;
- Step 3b combines GPU NEML2 and CUDA Kokkos assembly with CPU PETSc;
- Step 4 additionally uses Kokkos PETSc vectors and AIJKokkos;
- all steps reach `ux_right` approximately `0.005` and `ux_center`
  approximately `0.0025` at `t=0.005`;
- the nonlinear residual histories show real Newton convergence rather than an
  iteration-zero false convergence; and
- `preset = false` versus the default preset behavior is tested as an explicit
  A/B comparison, including whether Kokkos automatic scaling completes.

Do not interpret Step 3/4 timing results until this correctness comparison
passes.
