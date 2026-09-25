//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "KokkosKernelGrad.h"

/**
 * Kokkos total-Lagrangian divergence of a first Piola-Kirchhoff stress.
 */
class KokkosTotalLagrangianStressDivergence : public Moose::Kokkos::KernelGrad
{
  using Real3 = Moose::Kokkos::Real3;

public:
  static InputParameters validParams();

  KokkosTotalLagrangianStressDivergence(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION Real3 precomputeQpResidual(const unsigned int qp, AssemblyDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION Real3 precomputeQpJacobian(const unsigned int j,
                                              const unsigned int qp,
                                              AssemblyDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION Real3 precomputeQpOffDiagJacobian(const unsigned int j,
                                                     const unsigned int jvar,
                                                     const unsigned int qp,
                                                     AssemblyDatum & datum) const;

private:
  KOKKOS_FUNCTION Real3 jacobian(const unsigned int displacement_component,
                                  const Real3 & grad_phi,
                                  const unsigned int qp,
                                  AssemblyDatum & datum) const;

  const unsigned int _component;
  const unsigned int _ndisp;
  Moose::Kokkos::MaterialProperty<Real, 2> _pk1_stress;
  Moose::Kokkos::MaterialProperty<Real, 4> _dpk1_d_grad_u;
  Moose::Kokkos::Array<unsigned int> _displacement_var_ids;
};

template <typename Derived>
KOKKOS_FUNCTION Moose::Kokkos::Real3
KokkosTotalLagrangianStressDivergence::precomputeQpResidual(const unsigned int qp,
                                                            AssemblyDatum & datum) const
{
  Real3 residual(0);
  const auto pk1 = _pk1_stress(datum, qp);
  for (unsigned int J = 0; J < 3; ++J)
    residual(J) = pk1(_component, J);
  return residual;
}

template <typename Derived>
KOKKOS_FUNCTION Moose::Kokkos::Real3
KokkosTotalLagrangianStressDivergence::precomputeQpJacobian(const unsigned int j,
                                                            const unsigned int qp,
                                                            AssemblyDatum & datum) const
{
  return jacobian(_component, _grad_phi(datum, j, qp), qp, datum);
}

template <typename Derived>
KOKKOS_FUNCTION Moose::Kokkos::Real3
KokkosTotalLagrangianStressDivergence::precomputeQpOffDiagJacobian(
    const unsigned int j,
    const unsigned int jvar,
    const unsigned int qp,
    AssemblyDatum & datum) const
{
  for (unsigned int component = 0; component < _ndisp; ++component)
    if (_displacement_var_ids[component] == jvar)
      return jacobian(component, _grad_phi(datum, j, qp), qp, datum);
  return Real3(0);
}

KOKKOS_FUNCTION Moose::Kokkos::Real3
KokkosTotalLagrangianStressDivergence::jacobian(const unsigned int displacement_component,
                                                const Real3 & grad_phi,
                                                const unsigned int qp,
                                                AssemblyDatum & datum) const
{
  Real3 result(0);
  const auto tangent = _dpk1_d_grad_u(datum, qp);
  for (unsigned int J = 0; J < 3; ++J)
    for (unsigned int L = 0; L < 3; ++L)
      result(J) += tangent(_component, J, displacement_component, L) * grad_phi(L);
  return result;
}
