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
 * Copies a symmetric or full rank-two or rank-four NEML2 output into a Kokkos material property.
 */
template <unsigned int rank, bool full = false>
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
using NEML2ToKokkosFullRankTwoMaterialProperty = NEML2ToKokkosMaterialProperty<2, true>;
using NEML2ToKokkosFullRankFourMaterialProperty = NEML2ToKokkosMaterialProperty<4, true>;

template <unsigned int rank, bool full>
KOKKOS_FUNCTION unsigned int
NEML2ToKokkosMaterialProperty<rank, full>::mandelIndex(const unsigned int i, const unsigned int j)
{
  if (i == j)
    return i;
  if ((i == 1 && j == 2) || (i == 2 && j == 1))
    return 3;
  if ((i == 0 && j == 2) || (i == 2 && j == 0))
    return 4;
  return 5;
}

template <unsigned int rank, bool full>
KOKKOS_FUNCTION Real
NEML2ToKokkosMaterialProperty<rank, full>::mandelFactor(const unsigned int index)
{
  return index < 3 ? 1.0 : MathUtils::sqrt2;
}

template <unsigned int rank, bool full>
template <typename Derived>
KOKKOS_FUNCTION void
NEML2ToKokkosMaterialProperty<rank, full>::computeQpProperties(const unsigned int qp,
                                                                 Datum & datum) const
{
  auto prop = _prop(datum, qp);
  const auto batch =
      _broadcast_output ? 0 : _element_batch_offset[datum.elemID()] + static_cast<dof_id_type>(qp);

  // Cache the mesh dimension in a local before the write loop. Using
  // libMesh::make_range() with either the device-side kokkosAssembly() accessor or a cached
  // local produces an empty range on CUDA in this compilation unit, so the loop body never
  // executes and _prop stays at the value the device View was created with (zero). Raw
  // counter loops against the cached dimension iterate correctly.
  const unsigned int dim = kokkosAssembly().getDimension();
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
    {
      if constexpr (full)
      {
        if constexpr (rank == 2)
          prop(i, j) = _staged_output[batch * 9 + 3 * i + j];
        else
          for (unsigned int k = 0; k < dim; ++k)
            for (unsigned int l = 0; l < dim; ++l)
              prop(i, j, k, l) =
                  _staged_output[batch * 81 + 27 * i + 9 * j + 3 * k + l];
      }
      else
      {
        const auto a = mandelIndex(i, j);
        if constexpr (rank == 2)
          prop(i, j) = _staged_output[batch * 6 + a] / mandelFactor(a);
        else
          for (unsigned int k = 0; k < dim; ++k)
            for (unsigned int l = 0; l < dim; ++l)
            {
              const auto b = mandelIndex(k, l);
              prop(i, j, k, l) = _staged_output[batch * 36 + 6 * a + b] /
                                 (mandelFactor(a) * mandelFactor(b));
            }
      }
    }
}
