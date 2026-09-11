//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "KokkosMaterial.h"
#include "NEML2OutputInterface.h"
#include "MathUtils.h"

/**
 * Copies a symmetric rank-two or rank-four NEML2 output into a Kokkos material property.
 */
template <unsigned int rank>
class NEML2ToKokkosMaterialProperty : public Moose::Kokkos::Material,
                                      public NEML2OutputInterface
{
  static_assert(rank == 2 || rank == 4, "Only rank-two and rank-four properties are supported.");

public:
  static InputParameters validParams();

  NEML2ToKokkosMaterialProperty(const InputParameters & params);
  NEML2ToKokkosMaterialProperty(const NEML2ToKokkosMaterialProperty & object);

  void meshChanged() override;
  void computeProperties() override;

  template <typename Derived>
  KOKKOS_FUNCTION void computeQpProperties(const unsigned int qp, Datum & datum) const;

private:
  KOKKOS_FUNCTION static unsigned int mandelIndex(unsigned int i, unsigned int j);
  KOKKOS_FUNCTION static Real mandelFactor(unsigned int index);

  Moose::Kokkos::MaterialProperty<Real, rank> _prop;
  Moose::Kokkos::Array<Real> _staged_output;
  Moose::Kokkos::Array<dof_id_type> _element_batch_offset;
  bool _batch_mapping_outdated = true;
  bool _broadcast_output = false;
  int64_t _batch_qp_stride = 0;
};

using NEML2ToKokkosRankTwoMaterialProperty = NEML2ToKokkosMaterialProperty<2>;
using NEML2ToKokkosRankFourMaterialProperty = NEML2ToKokkosMaterialProperty<4>;

template <unsigned int rank>
KOKKOS_FUNCTION unsigned int
NEML2ToKokkosMaterialProperty<rank>::mandelIndex(const unsigned int i, const unsigned int j)
{
  if (i == j)
    return i;
  if ((i == 1 && j == 2) || (i == 2 && j == 1))
    return 3;
  if ((i == 0 && j == 2) || (i == 2 && j == 0))
    return 4;
  return 5;
}

template <unsigned int rank>
KOKKOS_FUNCTION Real
NEML2ToKokkosMaterialProperty<rank>::mandelFactor(const unsigned int index)
{
  return index < 3 ? 1.0 : MathUtils::sqrt2;
}

template <unsigned int rank>
template <typename Derived>
KOKKOS_FUNCTION void
NEML2ToKokkosMaterialProperty<rank>::computeQpProperties(const unsigned int qp,
                                                          Datum & datum) const
{
  auto prop = _prop(datum, qp);
  const auto batch =
      _broadcast_output ? 0 : _element_batch_offset[datum.elemID()] + static_cast<dof_id_type>(qp);

  for (const auto i : make_range(kokkosAssembly().getDimension()))
    for (const auto j : make_range(kokkosAssembly().getDimension()))
    {
      const auto a = mandelIndex(i, j);
      if constexpr (rank == 2)
        prop(i, j) = _staged_output[batch * 6 + a] / mandelFactor(a);
      else
        for (const auto k : make_range(kokkosAssembly().getDimension()))
          for (const auto l : make_range(kokkosAssembly().getDimension()))
          {
            const auto b = mandelIndex(k, l);
            prop(i, j, k, l) = _staged_output[batch * 36 + 6 * a + b] /
                               (mandelFactor(a) * mandelFactor(b));
          }
    }
}
