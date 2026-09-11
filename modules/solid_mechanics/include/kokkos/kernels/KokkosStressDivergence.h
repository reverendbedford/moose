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
 * Kokkos kernel for the Cartesian divergence of a stress material property.
 */
class KokkosStressDivergence : public Moose::Kokkos::KernelGrad
{
  using Real3 = Moose::Kokkos::Real3;

public:
  static InputParameters validParams();

  KokkosStressDivergence(const InputParameters & parameters);

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
  Moose::Kokkos::MaterialProperty<Real, 2> _stress;
  Moose::Kokkos::MaterialProperty<Real, 4> _tangent;
  Moose::Kokkos::Array<unsigned int> _displacement_var_ids;
};

template <typename Derived>
KOKKOS_FUNCTION Moose::Kokkos::Real3
KokkosStressDivergence::precomputeQpResidual(const unsigned int qp, AssemblyDatum & datum) const
{
  Real3 residual(0);
  const auto stress = _stress(datum, qp);
  for (const auto j : make_range(_ndisp))
    residual(j) = stress(_component, j);
  return residual;
}

template <typename Derived>
KOKKOS_FUNCTION Moose::Kokkos::Real3
KokkosStressDivergence::precomputeQpJacobian(const unsigned int j,
                                            const unsigned int qp,
                                            AssemblyDatum & datum) const
{
  return jacobian(_component, _grad_phi(datum, j, qp), qp, datum);
}

template <typename Derived>
KOKKOS_FUNCTION Moose::Kokkos::Real3
KokkosStressDivergence::precomputeQpOffDiagJacobian(const unsigned int j,
                                                   const unsigned int jvar,
                                                   const unsigned int qp,
                                                   AssemblyDatum & datum) const
{
  for (const auto component : make_range(_ndisp))
    if (_displacement_var_ids[component] == jvar)
      return jacobian(component, _grad_phi(datum, j, qp), qp, datum);

  return Real3(0);
}

KOKKOS_FUNCTION Moose::Kokkos::Real3
KokkosStressDivergence::jacobian(const unsigned int displacement_component,
                                const Real3 & grad_phi,
                                const unsigned int qp,
                                AssemblyDatum & datum) const
{
  Real3 result(0);
  const auto tangent = _tangent(datum, qp);
  for (const auto j : make_range(_ndisp))
    for (const auto l : make_range(_ndisp))
      result(j) += 0.5 * (tangent(_component, j, displacement_component, l) +
                          tangent(_component, j, l, displacement_component)) *
                   grad_phi(l);
  return result;
}
