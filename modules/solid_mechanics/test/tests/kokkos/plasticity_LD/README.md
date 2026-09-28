# plasticity_LD

Large-deformation counterpart of `../plasticity`. Same script interface,
step layout, and diagnostics; the constitutive path is replaced by the
Kokkos total-Lagrangian NEML2 bridge introduced by commit `3a7455d0`.

## What is shared with ../plasticity

- `analyze_results.py`, `plot_comparison.py`, and `select_cuda_device.py`
  are symlinked so the SD and LD suites can be summarized and plotted the
  same way.
- `run_benchmarks.sh` mirrors the SD driver: same PETSc argument sets, same
  `--compute-device=cuda` selection, same GNU time / perf-graph JSON / Nsight
  Systems capture, same step names (`step1_plasticity_cpu_neml2`, ...).
- Mesh, boundary conditions, loading (`u_x = t` on the right face), and the
  outer solver settings (`gamg`/`gmres`, `dt=1e-3`, `num_steps=5`,
  `nl_abs_tol=1e-10`, `residual_and_jacobian_together=true`) are preserved
  where they remain valid.

## What is different, and why

- **Formulation.** `strain = FINITE`, `formulation = TOTAL`, and
  `volumetric_locking_correction = true` replace the small-strain `SMALL`
  formulation in step1/step2. Steps 3, 4, and 5 wire the low-level Kokkos
  total-Lagrangian objects directly.
- **Solver.** `pc_type = hypre -pc_hypre_type boomeramg -ksp_type gmres`
  with `line_search = none`. The SD suite uses `gamg`/`gmres` for the same
  scalability reason a benchmark wants an iterative solver: LU dominates
  wall time on non-trivial meshes and hides the GPU speedup we are trying
  to measure. `gamg` itself does not work well on the non-symmetric CP LD
  tangent (converges only with backtracking, ~21 Newton per two-step run
  vs 8), so `hypre boomeramg` is the substitute — same Newton count as LU
  with a scalable preconditioner. Alternatives run at N=16 / num_steps=2
  with the step3 config, all producing the same displacement to machine
  precision:

    | pc                                                          | wall (s) | Newton |
    | ----------------------------------------------------------- | -------: | -----: |
    | `hypre boomeramg` + `gmres` (default)                       |       63 |      8 |
    | `ilu -pc_factor_levels 2` + `gmres`                         |       64 |      8 |
    | `lu -pc_factor_mat_solver_type mumps`                       |       80 |      8 |
    | `asm -sub_pc_type lu`                                       |       88 |      8 |
    | `lu`                                                        |       89 |      8 |
    | `gamg` + `gmres` + `snes_linesearch_type bt`                |      193 |     21 |

  LU is the LD reference in this module (every existing Lagrangian test
  under `modules/solid_mechanics/test/tests/lagrangian/` uses it and so
  does `../../neml2/crystal_plasticity/{exact,approx}_kinematics.i`). Use
  it via CLI override when you need bit-identical agreement with those
  references:

      -pc_type lu -pc_hypre_type "" -ksp_type preonly \
      Executioner/petsc_options_iname=-pc_type Executioner/petsc_options_value=lu
- **NEML2 model.** `../../neml2/crystal_plasticity/exact_kinematics_neml2.i`
  replaces `../../neml2/plasticity/perfect_neml2.i`. The SD suite exercises a
  small-strain J2 return map that takes a symmetric strain and returns a
  symmetric Cauchy stress; the LD suite exercises the crystal-plasticity
  exact-kinematics model that takes the deformation gradient and returns a
  non-symmetric PK2 stress with `plastic_deformation_gradient` state.
  Material parameters (elastic constants, yield, hardening) do not map
  between the two constitutive theories, so they are not preserved; the
  loading and mesh are.
- **Input gathering.** `TorchDeformationGradient` replaces
  `TorchSmallStrain`, so NEML2 receives the full non-symmetric
  `F = I + grad(u)`.
- **Output transfer.** `NEML2ToKokkosFullRankTwoMaterialProperty` and
  `NEML2ToKokkosFullRankFourMaterialProperty` replace their symmetric
  counterparts, moving the full PK2 stress and full `dS/dF` from NEML2 to
  Kokkos without Mandel packing or index permutation.
- **PK1 conversion.** `KokkosComputeLagrangianStressCustomPK2` converts
  `(PK2, dS/dF)` to `(PK1, dP/d(grad u))` on device before assembly.
- **Kernel.** `KokkosTotalLagrangianStressDivergence` replaces
  `KokkosStressDivergence` in steps 3a, 3b, 4.
- **Initial state.** `initialize_outputs = 'plastic_deformation_gradient'` +
  `initialize_output_values = 'initial_plastic_defgrad'` and the auxiliary
  `GenericConstantRankTwoTensor` seed the stateful CP model at t = 0.
- **Executioner override knobs.** `run_benchmarks.sh` exposes `NUM_STEPS` and
  `DT` so a caller can crank up loading (e.g. `NUM_STEPS=10 DT=5e-3` for
  5% engineering strain) without editing inputs.

## Reference

The host CPU reference is
`modules/solid_mechanics/test/tests/neml2/crystal_plasticity/exact_kinematics.i`.
Step1 here reproduces it at `N=16` with the SD benchmark's solver stack
(gamg/gmres, dt = 1e-3, 5 steps) and preserves the SD script interface.

## Steps

| Step | Kokkos assembly | NEML2 | PETSc |
| ---- | ---- | ---- | ---- |
| step1 | CPU MOOSE (`ComputeLagrangianStressCustomPK2`) | CPU | CPU |
| step2 | CPU MOOSE (`ComputeLagrangianStressCustomPK2`) | CUDA | CPU |
| step3 | CUDA (`KokkosTotalLagrangianStressDivergence` + Kokkos PK2->PK1) | CPU | CPU |
| step4 | CUDA (`KokkosTotalLagrangianStressDivergence` + Kokkos PK2->PK1) | CUDA | CPU |
| step5 | CUDA (`KokkosTotalLagrangianStressDivergence` + Kokkos PK2->PK1) | CUDA | CUDA (`aijkokkos` + Kokkos vec) |

## Running

```
# smoke, N=8 unit cube, single repetition, no nsys
MESH_N=8 REPEATS=1 PROFILE=0 ./run_benchmarks.sh

# default benchmark (matches SD): N=16, 3 repetitions, one nsys per step
./run_benchmarks.sh

# larger deformation exercising finite kinematics
DT=5e-3 NUM_STEPS=10 ./run_benchmarks.sh
```
