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
 * Checks full rank-two and rank-four Kokkos properties produced by S = A F B.
 */
class KokkosFullTensorTestMaterial : public Moose::Kokkos::Material
{
public:
  static InputParameters validParams();

  KokkosFullTensorTestMaterial(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION void computeQpProperties(const unsigned int qp, Datum & datum) const;

private:
  Moose::Kokkos::MaterialProperty<Real, 2> _rank_two;
  Moose::Kokkos::MaterialProperty<Real, 4> _rank_four;
  Moose::Kokkos::MaterialProperty<Real, 2> _pk1;
  Moose::Kokkos::MaterialProperty<Real> _rank_two_error;
  Moose::Kokkos::MaterialProperty<Real> _rank_four_error;
  Moose::Kokkos::MaterialProperty<Real> _pk1_error;
  Moose::Kokkos::Array<Real> _A;
  Moose::Kokkos::Array<Real> _B;
  Moose::Kokkos::Array<Real> _F;
};

template <typename Derived>
KOKKOS_FUNCTION void
KokkosFullTensorTestMaterial::computeQpProperties(const unsigned int qp, Datum & datum) const
{
  const auto rank_two = _rank_two(datum, qp);
  const auto rank_four = _rank_four(datum, qp);
  const auto pk1 = _pk1(datum, qp);
  Real rank_two_error = 0;
  Real rank_four_error = 0;
  Real pk1_error = 0;

  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
    {
      Real expected = 0;
      for (unsigned int k = 0; k < 3; ++k)
        for (unsigned int l = 0; l < 3; ++l)
          expected += _A[3 * i + k] * _F[3 * k + l] * _B[3 * l + j];
      const auto error = ::Kokkos::abs(rank_two(i, j) - expected);
      rank_two_error = error > rank_two_error ? error : rank_two_error;

      Real expected_pk1 = 0;
      for (unsigned int K = 0; K < 3; ++K)
        expected_pk1 += _F[3 * i + K] * rank_two(K, j);
      const auto pk1_component_error = ::Kokkos::abs(pk1(i, j) - expected_pk1);
      pk1_error = pk1_component_error > pk1_error ? pk1_component_error : pk1_error;

      for (unsigned int k = 0; k < 3; ++k)
        for (unsigned int l = 0; l < 3; ++l)
        {
          const auto derivative_error =
              ::Kokkos::abs(rank_four(i, j, k, l) - _A[3 * i + k] * _B[3 * l + j]);
          rank_four_error =
              derivative_error > rank_four_error ? derivative_error : rank_four_error;
        }
    }

  _rank_two_error(datum, qp) = rank_two_error;
  _rank_four_error(datum, qp) = rank_four_error;
  _pk1_error(datum, qp) = pk1_error;
}
