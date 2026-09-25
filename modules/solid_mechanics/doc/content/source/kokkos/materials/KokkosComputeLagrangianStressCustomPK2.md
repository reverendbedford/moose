# KokkosComputeLagrangianStressCustomPK2

!if! function=hasCapability('kokkos')

`KokkosComputeLagrangianStressCustomPK2` converts a full second
Piola-Kirchhoff stress $S$ and its derivative with respect to the deformation gradient to first
Piola-Kirchhoff stress:

\[
P_{iJ} = F_{iK} S_{KJ},
\qquad
\frac{\partial P_{iJ}}{\partial F_{kL}} =
\delta_{ik} S_{LJ} + F_{iK}\frac{\partial S_{KJ}}{\partial F_{kL}}.
\]

The input properties are full tensors rather than symmetric Mandel tensors. The material requires
three scalar Cartesian displacement variables on a three-dimensional reference mesh.

## Example Input Syntax

!listing solid_mechanics/test/tests/kokkos/neml2_input_gatherer/exact_kinematics.i block=Materials/pk1

!syntax parameters /Materials/KokkosComputeLagrangianStressCustomPK2

!syntax inputs /Materials/KokkosComputeLagrangianStressCustomPK2

!syntax children /Materials/KokkosComputeLagrangianStressCustomPK2

!if-end!

!else
!include kokkos/kokkos_warning.md
