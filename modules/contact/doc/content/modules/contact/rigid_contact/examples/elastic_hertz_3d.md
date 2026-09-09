# Elastic Hertz Contact, 3D Quarter Symmetry

## Description

Rigid sphere pressed into an elastic quarter-sphere.  This is the
introductory rigid-body contact example --- displacement-controlled,
small-strain, linear elasticity, and an analytic sphere contactor ---
and is also the workhorse regression case for the formulation.

Symmetry planes are placed at $x = 0$ and $z = 0$, so the model is one
quarter of a sphere-on-sphere Hertzian setup:

- Subdomain `1` is the deformable quarter-sphere, $E = 1.40625\times 10^7$,
  $\nu = 0.25$, radius 2.  Its curved bottom carries sideset `100`; the top
  face carries sideset `2` and is driven by a `FunctionDirichletBC`.
- The mesh file also contains a rigid-indenter block (`1000`) which is
  stripped by a `BlockDeletionGenerator` before the run.  The analytic
  [SphereContactor.md] replaces it.
- The sidesets `1` and `3` are the symmetry planes.

The Hertz reference (with $E^\ast = 1.5\times 10^7$, $R_{\text{eff}} = 1$)
predicts, at $\text{depth}\ d = 0.01$, contact radius
$a = \sqrt{R d} = 0.1$ and peak pressure
$p_0 = 2 E^\ast a / (\pi R) = 9.55\times 10^5$.  The regression output
reproduces both to within mesh discretization error.

## Contact setup

The [`[RigidContact]`](syntax/RigidContact/index.md) sub-block expands to
the full analytic-level-set stack:
[LowerDBlockFromSidesetGenerator](/meshgenerators/LowerDBlockFromSidesetGenerator.md),
[RigidBodyContactSparsity.md], the `normal_lm` variable with
`ConstantBounds`, [RigidBodyNodalNCPKernel.md],
[RigidBodyNormalMechanicalContact.md] on each displacement component, the
problem-coverage relaxation, and an SMP full preconditioner.  The
executioner uses PETSc `SNESVINEWTONSSLS` to enforce $\lambda \ge 0$.

!listing modules/contact/examples/rigid/elastic/hertz_elastic_3d.i block=RigidContact

## Full input

!listing modules/contact/examples/rigid/elastic/hertz_elastic_3d.i

## Suggested visualization

- Contour of `normal_lm` on the contact face (sideset `100`) at the last
  time step.  This shows the disk-shaped contact patch predicted by
  Hertz theory and the pressure falloff toward its edge.
- Contour of the vertical stress `stress_yy` on a diametral cut.

Rendered images will be added in a follow-up commit; the input as it ships
runs in a couple of seconds on a workstation and produces both fields in
the ExodusII output.
