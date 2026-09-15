# Kokkos assembly regression: 3D solid-mechanics residual comes out identically zero

Branch: `neml2-v3-kokkos-solid-mechanics`
Reproduced: 2026-09-14 on tree at `34bb5b9855`, built through
`kokkos-cuda-stack/scripts/activate.sh`, run with `--compute-device=cpu` or
`--compute-device=cuda` (both fail identically).

## Symptom

The plasticity benchmark `modules/solid_mechanics/test/tests/kokkos/plasticity/`
runs to completion and prints `Solve Converged!`, but produces wrong physics:

| Step | `ux_center` at `t=0.005` | Expected |
| --- | --- | --- |
| 1 (CPU MOOSE + CPU NEML2) | `~2.5e-3` | `~2.5e-3` |
| 2 (CPU MOOSE + GPU NEML2) | `~2.5e-3` | `~2.5e-3` |
| 3a (Kokkos + CPU NEML2 + CPU PETSc) | `-8.208e-05` | `~2.5e-3` |
| 3b (Kokkos + GPU NEML2 + CPU PETSc) | `-8.208e-05` | `~2.5e-3` |
| 4  (Kokkos + GPU NEML2 + GPU PETSc)  | `+4.308e-05` | `~2.5e-3` |

Every Kokkos step reports `0 Nonlinear |R| = 0.000000e+00` at every time step
and "converges" without doing any real Newton work. `ux_right` tracks the BC
correctly at every step; only the interior is bogus.

## The regression is **not** in the Torch/NEML2 layer

The bug reproduces without NEML2 or Torch anywhere in the input. Two
pre-existing solid-mechanics tests on this branch already fail EXODIFF against
their own gold files on the current tree:

```
test:kokkos/linear_elasticity.material_3d ..................... FAILED (EXODIFF)
test:kokkos/linear_elasticity.neml2_material_3d ............... FAILED (EXODIFF)
```

Both use pure `KokkosStressDivergence` / `KokkosIsotropicElasticity` with a
constant `KokkosComputeIsotropicElasticity` tangent material. The gold files
have the correct nonzero `disp_x` field; the current code produces
`disp_x = 0` at every node except the boundary.

The corresponding `PetscJacobianTester` variant
(`neml2_material_3d_jacobian`) is a **false positive**: when the residual
identically doesn't depend on the solution, analytical and FD Jacobians are
both zero and agree with each other. That's why CI hasn't caught this.

## What I verified with device-side `printf` instrumentation

Running `modules/solid_mechanics/test/tests/kokkos/linear_elasticity/kokkos_linear_elasticity_3d.i`
(pure Kokkos, no NEML2, no Torch), on a 1×1×1 mesh with `preset=true` on the
right-face DirichletBC:

1. **FE gradients are correct.** `FESystem::operator()` writes
   `_qp_solutions_grad[tag=2](sid, 0)[qp=0]` for the right-boundary element
   with `sum|dof|=0.004`, `grad_disp_x = (0.002, ~0, ~0)`. All three
   displacement variables have their arrays allocated and populated.
2. **`precomputeQpResidual` is called and returns the right stress.** A print
   at the exit of `KokkosIsotropicElasticity::precomputeQpResidual` shows
   `resid = (0.024, ~0, ~0)` for the boundary element under `stress_x`.
3. **`Kernel::computeResidual` dispatches all three kernels.**
   `numBlockElems=8`, `_thread.size=16`, `Dispatcher::parallelFor` runs.
4. **`accumulateTaggedElementalResidual` is never invoked with a nonzero
   `local_re`.** A print at the top of that function (gated only on
   `i == 0`, so it should fire once per dispatched element regardless of
   value) never fires.
5. **The residual vector on the host is all zero after the device→host
   sync.** A print in `Vector::close()` reads
   `size=24 max|dof|=0 nan_count=0` on every `close`.
6. **Under real SNES** (`-snes_monitor -snes_converged_reason`) the sequence
   is: initial `|R|=9.327e-3` → linear solve fails
   (`DIVERGED_PC_FAILED`) → outer solve fails
   (`DIVERGED_FUNCTION_NANORINF`) → MOOSE cuts the time step and retries
   with half the load → same failure → … → eventually `|R|` gets small
   enough that SNES exits at iteration 0 with `|R| ≈ 0` and MOOSE prints
   `Solve Converged!`. That is how the deceptively clean
   `0 Nonlinear |R| = 0.000000e+00` line in the user's log is produced.

So the residual value is being **computed correctly at the QP level** but is
**lost between `precomputeQpResidual` and the tagged residual vector**. The
gap contains only the KernelGrad body lambda's accumulation
`local_re[i] += value * _grad_test.reference(datum, i, qp)` and the
post-body loop in `ResidualObject::computeResidualInternal` that calls
`accumulateTaggedElementalResidual(local_re[i - ib], …)`.

Also observed (unclear if cause or symptom): the KernelGrad inner qp loop's
`printf` inside the loop body fires only for `qp=0` per thread — the loop
never reaches `qp=1..7`. Some invocations at `preR:exit` also return
`resid = (-nan, -nan, -nan)`, which propagates to the Jacobian and is what
trips PETSc's NANORINF.

## What is *not* the bug (already ruled out)

- **`residual_and_jacobian_together = true`.** Flipping it to `false` in
  Step 3a produces bit-identical wrong values.
- **`preset=false` on `disp_x_right`.** The pre-existing tests use the
  default `preset=true` and fail the same way.
- **`TorchFEInterpolation` caching a stale `_petsc_solution` pointer.** My
  first attempt at a fix (re-fetching the pointer on every call) compiles
  and links cleanly but produces bit-identical wrong values. I've reverted
  that fix — the tree is clean.
- **NEML2 tangent conversion.** Verified the Mandel↔tensor factors in
  `NEML2ToKokkosMaterialProperty::computeQpProperties` are correct and
  reproduce the isotropic elastic tensor exactly.
- **`KernelGrad`'s `J^T` factorization math.** Walked the algebra: with
  `datum.J(qp)` = `J^{-1}` and `datum.J(qp).transpose() * (JxW * pre) ·
  grad_test.reference` it does reduce to the correct
  `JxW * pre · grad_test_physical`.
- **`_active_variables` missing a displacement component.** `FESystem::op`
  fires for var=0, 1, 2 with correct DOF sums, so all three disp arrays
  are allocated.

## Where the bug most plausibly lives

The Kokkos assembly framework, in the code path
`Kernel::computeResidual` → `Dispatcher::parallelFor` →
`Kernel::operator()(ResidualLoop, tid, kernel)` →
`kernel.computeResidualInternal(kernel, datum)` →
`KernelGrad::computeResidualInternal` →
`ResidualObject::computeResidualInternal` →
`accumulateTaggedElementalResidual`.

Recent commits on the branch that touch this exact path and are suspect:

- `d199fc241c` — "Support computing residual and Jacobian together" (added
  `computeKokkosResidualAndJacobian`).
- `bee6d18fb4` — "Factor out jacobian application for optimization"
  (moved the `J^T` factor from the test-function side to the trial-function
  side; changed `VariableShapeGradient::operator()` to
  `J * reference(...)`; changed `KernelGrad::computeResidualInternal` to
  use `grad_test.reference` and `J.transpose() * (JxW * pre)`).
- `2e31626298` — "Make hooks consistent with MOOSE for optimized Kokkos
  kernels/BCs" (rewrote `hasUserJacobianHook` and
  `hasUserOffDiagJacobianHook` for `use_precompute_hooks`).

I could not narrow it further in one session — the `KokkosStressDivergence`
and `KokkosIsotropicElasticity` kernels themselves look consistent with the
new hooks, and my device-side prints in header templates behave weirdly
(some fire, some do not, even after wiping the object cache and touching
all `.K` files and rebuilding). Someone with more familiarity with the
Kokkos framework template-dispatch layout should look at the actual
generated code for one of these instantiations.

## Reproducers (both no NEML2, no Torch)

Fastest confirmation:

```
cd modules/solid_mechanics
PYTHONPATH=$MOOSE_DIR/python:$PYTHONPATH ./run_tests -i tests --re "material_3d$"
```

Expected: `2 passed, 2 FAILED (EXODIFF)`.

Direct run of the smaller reproducer:

```
cd modules/solid_mechanics/test/tests/kokkos/linear_elasticity
$MOOSE_DIR/modules/solid_mechanics/solid_mechanics-opt \
  -i kokkos_linear_elasticity_3d.i \
  --compute-device=cpu \
  Executioner/num_steps=1 -use_gpu_aware_mpi 0 \
  -snes_monitor -snes_converged_reason -ksp_converged_reason
```

Look for `DIVERGED_FUNCTION_NANORINF` and `DIVERGED_PC_FAILED` inside the
"converged" step 1.

## Recommendations

1. **Do not run the Step 3/4 benchmarks** in
   `modules/solid_mechanics/test/tests/kokkos/plasticity/` and do not
   compare Step 1/2 timing against them. With this broken assembly the Kokkos
   configurations aren't solving the physics.
2. Report to the Kokkos-MOOSE authors (Namjae Choi is the author of the
   suspect commits above) with the two reproducers.
3. Fix `PetscJacobianTester`'s false positive on zero-residual states so
   this class of regression can't slip through CI again — either sanity-check
   that the analytical Jacobian is not identically zero, or add an Exodiff
   or PointValue guard alongside every `_jacobian` variant.
4. Once the framework fix lands, rerun `run_coarse_exodus.sh` (which is in
   place and correct) and confirm all five configurations produce matching
   physical solutions and nontrivial Newton/KSP iteration counts **before**
   interpreting any timing results.
