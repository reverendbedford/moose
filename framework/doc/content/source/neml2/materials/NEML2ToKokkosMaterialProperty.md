# NEML2ToKokkosMaterialProperty

!if! function=hasCapability('neml2')

!alert note
Users are not expected to directly use these objects in an input file. Instead, use the [NEML2 action](syntax/NEML2/index.md).

## Description

This family of objects transfers a NEML2 symmetric tensor output or derivative to a Kokkos material property. The NEML2 executor moves the batched tensor from its compute device to the device used by Kokkos, if necessary, and the Kokkos material converts it from Mandel notation when evaluated.

The following objects are available:

| Class                                           | Kokkos material property rank |
| :---------------------------------------------- | :---------------------------- |
| `NEML2ToKokkosRankTwoMaterialProperty`          | Two                           |
| `NEML2ToKokkosRankFourMaterialProperty`         | Four                          |

The rank-two object expects a trailing NEML2 base shape of `[6]`. The rank-four object expects `[6, 6]`. The supported Mandel component order is `xx`, `yy`, `zz`, `yz`, `xz`, `xy`.

!syntax parameters /Materials/NEML2ToKokkosRankTwoMaterialProperty

!if-end!

!else

!include neml2/neml2_warning.md
