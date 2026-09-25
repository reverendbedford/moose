# TorchDeformationGradient

!if! function=hasCapability('neml2')

`TorchDeformationGradient` computes the full deformation gradient

\[
F_{iJ} = \delta_{iJ} + u_{i,J}
\]

from one to three scalar Cartesian displacement variables and supplies it as a NEML2 model input.
Missing displacement components have zero gradient, while the three-dimensional identity term is
retained.

## Limitations

- The formulation uses gradients on the reference mesh.
- Curvilinear coordinate terms are not included.

## Syntax

!syntax parameters /UserObjects/TorchDeformationGradient

## Example Input Files

!syntax inputs /UserObjects/TorchDeformationGradient

!if-end!

!else

!include neml2/neml2_warning.md
