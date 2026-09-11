# KokkosIsotropicElasticity

!if! function=hasCapability('kokkos')

`KokkosIsotropicElasticity` computes the Cartesian stress-divergence term for isotropic,
small-strain linear elasticity directly from the displacement gradients:

\[
\boldsymbol{\sigma} = \lambda\,\mathrm{tr}(\boldsymbol{\epsilon})\boldsymbol{I}
                     + 2\mu\boldsymbol{\epsilon}, \qquad
\boldsymbol{\epsilon} = \tfrac{1}{2}(\nabla\boldsymbol{u} + \nabla\boldsymbol{u}^T).
\]

Add one kernel for each displacement component. Two-dimensional problems use the plane-strain
assumption. The kernel does not support plane stress, curvilinear coordinates, finite strain,
eigenstrain, or spatially varying elastic constants.

## Example Input Syntax

!listing solid_mechanics/test/tests/kokkos/linear_elasticity/kokkos_linear_elasticity.i start=[Kernels] end=[] include-end=true

!syntax parameters /Kernels/KokkosIsotropicElasticity

!syntax inputs /Kernels/KokkosIsotropicElasticity

!syntax children /Kernels/KokkosIsotropicElasticity

!if-end!

!else
!include kokkos/kokkos_warning.md
