# KokkosComputeIsotropicElasticity

!if! function=hasCapability('kokkos')

`KokkosComputeIsotropicElasticity` computes stress and constitutive tangent Kokkos material
properties for Cartesian isotropic small-strain linear elasticity. In two dimensions, the
constitutive response uses the plane-strain assumption.

## Example Input Syntax

!listing solid_mechanics/test/tests/kokkos/linear_elasticity/kokkos_material_linear_elasticity.i start=[Materials] end=[] include-end=true

!syntax parameters /Materials/KokkosComputeIsotropicElasticity

!syntax inputs /Materials/KokkosComputeIsotropicElasticity

!syntax children /Materials/KokkosComputeIsotropicElasticity

!if-end!

!else
!include kokkos/kokkos_warning.md
