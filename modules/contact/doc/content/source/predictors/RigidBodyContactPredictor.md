# RigidBodyContactPredictor

!syntax description /Executioner/Predictor/RigidBodyContactPredictor

## Description

MOOSE `Predictor` that warm-starts the monolithic Newton solve by
resolving the contact subproblem --- LM DoFs on the contact sideset,
displacement DoFs within
[!param](/Executioner/Predictor/RigidBodyContactPredictor/k_hops) element
hops into the bulk, and (optionally) the load-control scalar --- with a
small local Newton sub-solve.  The initial guess handed to SNES then
already satisfies the local $(u, \lambda, s)$ coupling; the outer solve
only has to correct the far-field.

On sub-solve failure (non-convergence, KSP breakdown, bound-projection
thrashing) the initial guess is reverted to its entry state and only the
outer SNES runs.  On the shipped force-controlled tests the predictor
typically cuts cumulative nonlinear iterations by three to four and
halves wall time.

The predictor inherits the standard `enable` parameter from `MooseObject`
and declares it *controllable*, so the MOOSE Controls system can flip
the predictor on or off at runtime without touching the input file.
See the [theory page](modules/contact/rigid_contact/theory.md#warm-start-predictor)
for the overall role in the formulation.

## Parameters

!syntax parameters /Executioner/Predictor/RigidBodyContactPredictor
