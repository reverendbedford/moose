# Rigid body contact in MOOSE

Your job is to implement rigid body contact in MOOSE, up to and including example input files for solving the resulting numerical equations efficiently with PETSc.

## Math

Let $\Gamma_c$ be the contact surface.  In MOOSE each object will provide the contact residual/jacobian/weak form for a pair of surfaces: one on the deformaable finite element domain and the second with a rigid contactor.

Let $g_n(u)$ be the normal gap with $u$ the displacements such that $g_n>0$ is open.  The interpenetrability constraint we are trying to impose is then

$$ g_n(u) \ge 0 \, \mathrm{on} \, \Gamma_c$$

Define the traction imposed on the deformable body due to contact as

$$t_c = -p_n n$$

with $n$ the contact normal and $p_n \ge 0$ the contact pressure.  The constrained problem we're trying to solve is then defined by the conditions

$$g_n(u) \ge 0$$

$$p_n \ge 0$$

$$g_n(u)p_n = 0$$

For rigid frictionless contact we then only need to impose $t_c$ as a Neumann boundary condition on $\Gamma_c$, i.e.

$$ \sigma n = - p_n n $$

Of course the problem is the contact set will vary and we cannot prescribe it a priori.

I would like you to implement this in the usual way, using Lagrange multipliers to represent the contact pressure.  We'll do the active set update with the semismooth Newton method, using Fischer–Burmeister as the complementary function.

I do not want to use a penalty formulation and I would like as robust a method as possible.

## Approach

In general reuse as much of the existing MOOSE contact module objects as possible.  Our eventual goal is to allow the rigid contactor to be describe by a mesh (which would only be used to calculate the contact, no deformable finite elements).

Make all our code changes to the contact module and the solid_mechnics module.  Note: the actual module dependency runs `contact → solid_mechanics`, so `contact-opt` links both modules while `solid_mechanics-opt` does not.  Cross-module tests therefore live under `modules/contact/test/tests/` and run against `contact-opt`.

We need to handle the large deformation case.  I only care about supporting the "new" Lagrangian kernels.

I want exact algorithmic tangents and please do not use AD.

I suggest we first start with a simpler case: the contactor will be describe by some mathematical level set function.  Probably just start with a sphere for Hertzian contact.  After that works move towards the more general geometric formulation.

Structure the problem into a reasonable number of discrete steps.  Each commit should include new objects working towards our goal with a common theme.  Include MOOSE regression tests with each commit.

Keep working until you have tests that:
1. Reproduce the classical solution for spherical rigid contact with an elastic half space
2. Use an inelastic material model for the half-space and demonstrate good convergence
3. Same as #2 but for large deformations.
then take a break and wait for my input prior to moving to phase 2.

Phase 2 will then implement what's needed for general contact surfaces described by meshes.

Please use this file as a log of what you're doing.  You can check it in so we sync if we move machines, but I will be deleting it before we make a pull request with this new feature so don't rely on it for documentation.

---

# Phase 1 plan

## Context

Phase 1 delivers a Lagrange-multiplier rigid-body frictionless contact capability against an analytic level-set contactor (sphere first), integrated with the new-Lagrangian solid_mechanics kernels, with no penalty and no AD.  Phase 1 is done when three regression tests pass in `contact-opt`:

1. Rigid sphere on elastic half-space vs. analytical Hertz (small strain).
2. Same geometry with an inelastic half-space (radial-return plasticity), demonstrating good Newton convergence.
3. Same with the large-deformation UpdatedLagrangian pipeline.

## Decisions locked in

| Decision | Choice |
|---|---|
| Build environment | conda env `moose` — always `conda activate moose` before build/tests |
| Semismooth-Newton driver | **PETSc `SNESVINEWTONSSLS`** with one-sided bounds (`λ ∈ [0, +∞)`); PETSc's SSLS *is* the FB semismooth reformulation |
| LM discretization | **Dual (biorthogonal) P1 LM** on a lower-d block built from the deformable contact sideset (inf-sup stable, diagonal LM mass, matches existing `use_dual=true` mortar pattern) |
| Contactor input | New `LevelSetContactor` UO hierarchy with `SphereContactor` (center, radius) in Phase 1; `PlaneContactor` / `CylinderContactor` / `MeshContactor` are trivial follow-ons |
| Tests binary | `modules/contact/test/tests/rigid_body_contact/` (runs against `contact-opt`, which already links solid_mechanics) |

## Design

### New objects (all under `modules/contact/`)

1. `LevelSetContactor` — `GeneralUserObject`
   - `Real signedDistance(const Point & x) const` (g_LS(x))
   - `RealVectorValue normal(const Point & x) const` (∇g_LS(x))
   - optional `RealTensorValue hessian(const Point & x) const` for exact tangent bits
2. `SphereContactor : LevelSetContactor` — params `center` (Point), `radius` (Real).  Analytical SDF/normal/hessian; unit-tested.
3. `RigidBodyNormalLMContact : LowerDIntegratedBC`
   - Assembles the LM row on the secondary lower-d block:
     R_λ(node i) = ∫_{Γ_c^L} φ_i^dual(x) · g_LS(x + u(x)) dx
   - Diagonal (thanks to dual basis) → row-wise NCP driven by PETSc SSLS via bounds `λ ≥ 0`.
   - Uses a `LevelSetContactor` UO to evaluate g_LS and ∇g_LS at qp (displaced coords when `use_displaced_mesh=true`).
   - Exact `dR_λ/du_j = ∫ φ_i^dual · ∇g_LS · ∂(x+u)/∂u_j dx` (analytic; no AD).
4. `RigidBodyNormalMechanicalContact : LowerDIntegratedBC` (one instance per disp component)
   - Applies the contact traction on the coupled displacement equations:
     R_uk(i) = ∫_{Γ_c^L} λ n_k φ_i dx  (with n = ∇g_LS at qp).
   - Off-diagonal Jacobian couples `u_k` to `λ` and to other `u_l` via ∂n/∂x (when large-def is on).
5. `RigidBodyContactAction` (optional in Phase 1, add if input files get repetitive) — sugar to declare LM variable + bounds + constraint + traction BCs from a single `[Contact/RigidBody/…]` block.

### Reuse without modification

- `LowerDBlockFromSidesetGenerator` — builds the lower-d block Γ_c^L from the deformable contact sideset.
- `MooseVariableBase::_use_dual` — set `use_dual = true` on the LM variable in `[Variables]` to switch its FE to biorthogonal Lagrange.
- `BoundsBase` / `ConstantBounds` — one aux with `bound_type = lower`, `bound_value = 0` on the LM variable (upper defaults to +∞).
- `Moose::PetscSupport::isSNESVI` — already detects `-snes_type vinewtonssls`.
- Solid mechanics: `TotalLagrangianStressDivergence` / `UpdatedLagrangianStressDivergence` (`modules/solid_mechanics/include/kernels/lagrangian/`), `ComputeLagrangianStrain`, `ComputeLagrangianWrappedStress` + objective rate for wrapping small-strain inelastic models (`ComputeMultipleInelasticStress` + `IsotropicPlasticityStressUpdate`).

### Risks / open items

- **Dual basis on a stand-alone lower-d block** (outside `AutomaticMortarGeneration`): 99% confident it works via `use_dual = true` on the LM variable itself, but if libMesh's dual-basis routine assumes the mortar segment mesh, fall back to **P0 LM per lower-d element** for Phase 1 and revisit for Phase 2 (mortar path).  Verify in Commit 2.
- **Consistent tangent contribution from ∂n/∂x when using UL** (large-def): requires the hessian of the level set.  `SphereContactor::hessian` is trivial; the pattern generalizes.
- **`use_displaced_mesh` semantics** for our LM constraint: the level-set gap must be evaluated at deformed positions.  Integrate on the *reference* lower-d block with `x_current = x_ref + u`, so we do not rebuild the lower-d mesh every Newton step.

## Commit sequence (each commit ships code + regression tests)

All work on branch `rigid_body_contact`.

**Commit 1 — Level-set contactor infrastructure.**
- `LevelSetContactor` (base) + `SphereContactor`.
- Test: `modules/contact/test/tests/rigid_body_contact/level_set_contactor/` — probes SDF and normal at a handful of points via aux + CSV postprocessor; gold vs. closed form.

**Commit 2 — Frictionless rigid-body LM constraint + traction (small strain).**
- `RigidBodyNormalLMContact` (LM row on lower-d block) + `RigidBodyNormalMechanicalContact` (traction on disp components).
- Wire dual basis: try `use_dual = true`; if libMesh rejects it on non-mortar lower-d blocks, fall back to `MONOMIAL/CONSTANT` (P0).
- Wire SNESVI via `ConstantBounds` (lower=0) + `-snes_type vinewtonssls`.
- **Test 1 (Phase 1 goal #1):** 2D axisymmetric Hertz sphere on elastic half-space, `TotalLagrangianStressDivergenceAxisymmetricCylindrical` + `ComputeLagrangianLinearElasticStress`.  Gold on contact radius `a` and peak pressure `p₀` from analytical Hertz plus a pressure-profile CSV comparison at a few radii.

**Commit 3 — Inelastic half-space (small strain).**
- No new C++ classes if Commit 2 was well factored; possibly `RigidBodyContactAction` for input ergonomics.
- **Test 2 (Phase 1 goal #2):** same axisymmetric geometry with `ComputeMultipleInelasticStress` + `IsotropicPlasticityStressUpdate` (radial return J2).  Asserts on: solution reproducibility (gold CSV), plastic zone size, and *Newton iteration count vs. load step* PP (semismooth Newton should stay under a small bound each step).

**Commit 4 — Large deformation (UpdatedLagrangian) + inelastic.**
- Any adjustments to LM/traction classes needed for `use_displaced_mesh = true` in UL.  Add off-diagonal Jacobian contributions from ∂n/∂x if not already there.
- **Test 3 (Phase 1 goal #3):** same geometry, `UpdatedLagrangianStressDivergence` + `ComputeLagrangianWrappedStress` (objective rate) + `IsotropicPlasticityStressUpdate`.  Gold plus Newton-iteration PP.

**Commit 5 (optional polish, only if warranted).**
- 3D Hertz sphere-on-half-space verification test using `TotalLagrangianStressDivergence` (Cartesian).
- Documentation pages under `modules/contact/doc/content/source/{userobjects,constraints,bcs}/RigidBody*.md`.

## Critical files

New files (approx paths):
- `modules/contact/include/userobjects/LevelSetContactor.h`
- `modules/contact/include/userobjects/SphereContactor.h`
- `modules/contact/include/bcs/RigidBodyNormalLMContact.h`
- `modules/contact/include/bcs/RigidBodyNormalMechanicalContact.h`
- matching `.C` in `modules/contact/src/…`
- inputs + gold under `modules/contact/test/tests/rigid_body_contact/`

Read for patterns:
- `modules/contact/src/constraints/ComputeWeightedGapLMMechanicalContact.C` (min-NCP for LM row, dual-basis usage)
- `modules/contact/src/constraints/NormalMortarMechanicalContact.C` (traction assembly)
- `framework/src/bcs/LowerDIntegratedBC.C` (LL/LP/PL Jacobian pattern)
- `framework/src/bounds/BoundsBase.C` and `framework/src/utils/PetscSupport.C:1134` (SNESVI plumbing)
- `modules/solid_mechanics/test/tests/lagrangian/cartesian/total/convergence-auto/3D/neumann.i` and `.../lagrangian/materials/convergence/neohookean.i` (TL/UL test patterns)
- `modules/contact/test/tests/hertz_spherical/hertz_contact.i` (legacy penalty Hertz — reference solution to reproduce)
- `modules/contact/test/tests/bouncing-block-contact/frictionless-weighted-gap.i` (lower-d block + dual LM input pattern)

## Verification

Always `conda activate moose` first (per CLAUDE.md §6).  Ask before running a build if we have not already done so in the current turn.

Per commit:

    cd modules/contact
    make -j$(nproc)                                  # rebuilds contact-opt
    ./run_tests --re rigid_body_contact -j$(nproc)   # runs only the new tests

Golden-path checks that must pass to close Phase 1:
- Test 1 (elastic Hertz): analytical `p₀` and `a` within a stated tolerance; pressure-profile CSV diff.
- Test 2 (inelastic small-def): CSV gold + Newton-iteration PP below threshold (e.g. avg ≤ 6/load step).
- Test 3 (inelastic large-def): CSV gold + Newton-iteration PP below threshold.

If dual basis on the stand-alone lower-d block does not work, fall back to P0 LM (see Risks) and note it in the log below; do not switch to a penalty formulation.

## /goal — non-interactive entry point

Paste this to start a fresh Claude Code session working through Phase 1 without further planning:

    /goal Implement Phase 1 of the rigid-body contact plan documented in
    instructions.md. Ground rules:

    - Build environment is conda env `moose` (activate it before any build or
      test invocation).
    - Work on branch `rigid_body_contact`. One commit per plan step (Commits
      1-4, Commit 5 optional). Each commit must include working regression
      tests under `modules/contact/test/tests/rigid_body_contact/` and pass
      `contact-opt` `./run_tests` for the new tests. Do not squash.
    - Follow the design in instructions.md: dual-basis P1 LM on a lower-d
      block, `SNESVINEWTONSSLS` with lower bound 0 on the LM variable, no
      penalty, no AD, exact algorithmic tangents. Reuse `LowerDIntegratedBC`,
      `BoundsBase`, `LowerDBlockFromSidesetGenerator`, and the
      `TotalLagrangianStressDivergence` / `UpdatedLagrangianStressDivergence`
      kernel family. Fall back to P0 LM only if dual basis on the lower-d
      block genuinely does not work — document any fallback in the log below.
    - Follow the AGENTS.md guidelines in CLAUDE.md (Simplicity First, Surgical
      Changes, Goal-Driven Execution). Update the "Log" section of
      instructions.md after each commit (what shipped, any deviations, any
      newly-discovered risks).
    - Stop and wait for user input when: Phase 1 is complete (all three tests
      pass), or a test unexpectedly fails after best-effort debugging, or you
      hit a design ambiguity not covered by instructions.md.

## Log

### Commit 1 — Level-set contactor infrastructure

Shipped:
- `LevelSetContactor` abstract UO in `modules/contact/{include,src}/userobjects/`. Interface: `signedDistance`, `normal`, `hessian` (default zero). Convention: `g_LS > 0` outside the rigid body, `n = grad g_LS`.
- `SphereContactor` — analytic sphere with center + radius; hessian implemented so the class is UL-ready in Commit 4.
- `LevelSetContactorAux` — samples a contactor's signed distance or a normal component into an aux variable, for tests and visualization.
- Regression test `rigid_body_contact/level_set_contactor.sphere` at `modules/contact/test/tests/rigid_body_contact/level_set_contactor/`: 2 × 2 quad grid, `SphereContactor` at (1, 1) r = ½, gold-file compares every node's SDF and normal against the closed form (edges 0.5, corners √2 − ½, center −½).

Deviations from plan: none.

Newly-discovered risks: none.
