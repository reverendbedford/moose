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
 * Converts full PK2 stress and its deformation-gradient derivative to PK1 quantities.
 */
class KokkosComputeLagrangianStressCustomPK2 : public Moose::Kokkos::Material
{
public:
  static InputParameters validParams();

  KokkosComputeLagrangianStressCustomPK2(const InputParameters & parameters);

  void computeProperties() override;

  template <typename Derived>
  KOKKOS_FUNCTION void computeQpProperties(const unsigned int qp, Datum & datum) const;

private:
  const unsigned int _ndisp;
  const Moose::Kokkos::VariableGradient _grad_displacements;
  Moose::Kokkos::MaterialProperty<Real, 2> _pk2_stress;
  Moose::Kokkos::MaterialProperty<Real, 4> _dpk2_dF;
  Moose::Kokkos::MaterialProperty<Real, 2> _pk1_stress;
  Moose::Kokkos::MaterialProperty<Real, 4> _dpk1_d_grad_u;
  bool _need_jacobian;
  Moose::Kokkos::Scalar<bool> _need_jacobian_device;
};

template <typename Derived>
KOKKOS_FUNCTION void
KokkosComputeLagrangianStressCustomPK2::computeQpProperties(const unsigned int qp,
                                                            Datum & datum) const
{
  const auto pk2 = _pk2_stress(datum, qp);
  const auto dpk2_dF = _dpk2_dF(datum, qp);
  auto pk1 = _pk1_stress(datum, qp);
  auto dpk1_d_grad_u = _dpk1_d_grad_u(datum, qp);
  Real F[3][3];

  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      F[i][j] = (i == j) + _grad_displacements(datum, qp, i)(j);

  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int J = 0; J < 3; ++J)
    {
      pk1(i, J) = 0;
      for (unsigned int K = 0; K < 3; ++K)
        pk1(i, J) += F[i][K] * pk2(K, J);

      if (_need_jacobian_device)
        for (unsigned int k = 0; k < 3; ++k)
          for (unsigned int L = 0; L < 3; ++L)
          {
            dpk1_d_grad_u(i, J, k, L) = (i == k) * pk2(L, J);
            for (unsigned int K = 0; K < 3; ++K)
              dpk1_d_grad_u(i, J, k, L) += F[i][K] * dpk2_dF(K, J, k, L);
          }
    }
}
