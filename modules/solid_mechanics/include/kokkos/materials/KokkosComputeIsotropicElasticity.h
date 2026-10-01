//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "KokkosMaterial.h"

/**
 * Computes isotropic small-strain stress and tangent Kokkos material properties.
 */
class KokkosComputeIsotropicElasticity : public Moose::Kokkos::Material
{
public:
  static InputParameters validParams();

  KokkosComputeIsotropicElasticity(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION void computeQpProperties(const unsigned int qp, Datum & datum) const;

private:
  const unsigned int _ndisp;
  const Real _lambda;
  const Real _mu;
  const Moose::Kokkos::VariableGradient _grad_displacements;
  Moose::Kokkos::MaterialProperty<Real, 2> _stress;
  Moose::Kokkos::MaterialProperty<Real, 4> _tangent;
};

template <typename Derived>
KOKKOS_FUNCTION void
KokkosComputeIsotropicElasticity::computeQpProperties(const unsigned int qp, Datum & datum) const
{
  auto stress = _stress(datum, qp);
  auto tangent = _tangent(datum, qp);
  Real divergence = 0;
  // libMesh::make_range() is host-only and cannot be called from device code; use raw
  // counter loops.
  for (unsigned int i = 0; i < _ndisp; ++i)
    divergence += _grad_displacements(datum, qp, i)(i);

  for (unsigned int i = 0; i < _ndisp; ++i)
    for (unsigned int j = 0; j < _ndisp; ++j)
    {
      stress(i, j) = _lambda * divergence * (i == j) +
                     _mu * (_grad_displacements(datum, qp, i)(j) +
                            _grad_displacements(datum, qp, j)(i));

      for (unsigned int k = 0; k < _ndisp; ++k)
        for (unsigned int l = 0; l < _ndisp; ++l)
          tangent(i, j, k, l) = _lambda * (i == j) * (k == l) +
                                _mu * ((i == k) * (j == l) + (i == l) * (j == k));
    }
}
