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
  for (const auto i : make_range(_ndisp))
    divergence += _grad_displacements(datum, qp, i)(i);

  for (const auto i : make_range(_ndisp))
    for (const auto j : make_range(_ndisp))
    {
      stress(i, j) = _lambda * divergence * (i == j) +
                     _mu * (_grad_displacements(datum, qp, i)(j) +
                            _grad_displacements(datum, qp, j)(i));

      for (const auto k : make_range(_ndisp))
        for (const auto l : make_range(_ndisp))
          tangent(i, j, k, l) = _lambda * (i == j) * (k == l) +
                                _mu * ((i == k) * (j == l) + (i == l) * (j == k));
    }
}
