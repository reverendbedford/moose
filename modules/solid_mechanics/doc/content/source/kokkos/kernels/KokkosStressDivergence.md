# KokkosStressDivergence

!if! function=hasCapability('kokkos')

`KokkosStressDivergence` computes one Cartesian component of the stress-divergence term from
Kokkos material properties containing the stress and constitutive tangent. The stress property
has dimensions $n \times n$, and the tangent property has dimensions
$n \times n \times n \times n$, where $n$ is the number of displacement components.

The tangent is the derivative of stress with respect to symmetric small strain. Add one kernel
for each displacement component. This kernel does not support curvilinear coordinates, finite
strain, or volumetric locking correction.

## Example Input Syntax

!listing solid_mechanics/test/tests/kokkos/linear_elasticity/kokkos_material_linear_elasticity.i start=[Kernels] end=[] include-end=true

!syntax parameters /Kernels/KokkosStressDivergence

!syntax inputs /Kernels/KokkosStressDivergence

!syntax children /Kernels/KokkosStressDivergence

!if-end!

!else
!include kokkos/kokkos_warning.md
