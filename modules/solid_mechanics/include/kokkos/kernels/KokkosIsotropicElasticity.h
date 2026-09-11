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
 * Kokkos kernel for Cartesian isotropic small-strain elasticity.
 */
class KokkosIsotropicElasticity : public Moose::Kokkos::KernelGrad
{
  using Real3 = Moose::Kokkos::Real3;

public:
  static InputParameters validParams();

  KokkosIsotropicElasticity(const InputParameters & parameters);

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
                                 const Real3 & grad_phi) const;

  const unsigned int _component;
  const unsigned int _ndisp;
  const Real _lambda;
  const Real _mu;
  const Moose::Kokkos::VariableGradient _grad_displacements;
  Moose::Kokkos::Array<unsigned int> _displacement_var_ids;
};

template <typename Derived>
KOKKOS_FUNCTION Moose::Kokkos::Real3
KokkosIsotropicElasticity::precomputeQpResidual(const unsigned int qp,
                                                AssemblyDatum & datum) const
{
  Real3 residual(0);
  Real divergence = 0;

  for (const auto i : make_range(_ndisp))
  {
    const auto grad_disp = _grad_displacements(datum, qp, i);
    divergence += grad_disp(i);
    residual(i) = _mu * (grad_disp(_component) + _grad_displacements(datum, qp, _component)(i));
  }

  residual(_component) += _lambda * divergence;
  return residual;
}

template <typename Derived>
KOKKOS_FUNCTION Moose::Kokkos::Real3
KokkosIsotropicElasticity::precomputeQpJacobian(const unsigned int j,
                                                const unsigned int qp,
                                                AssemblyDatum & datum) const
{
  return jacobian(_component, _grad_phi(datum, j, qp));
}

template <typename Derived>
KOKKOS_FUNCTION Moose::Kokkos::Real3
KokkosIsotropicElasticity::precomputeQpOffDiagJacobian(const unsigned int j,
                                                       const unsigned int jvar,
                                                       const unsigned int qp,
                                                       AssemblyDatum & datum) const
{
  for (const auto i : make_range(_ndisp))
    if (_displacement_var_ids[i] == jvar)
      return jacobian(i, _grad_phi(datum, j, qp));

  return Real3(0);
}

KOKKOS_FUNCTION Moose::Kokkos::Real3
KokkosIsotropicElasticity::jacobian(const unsigned int displacement_component,
                                    const Real3 & grad_phi) const
{
  Real3 result(0);
  for (const auto i : make_range(_ndisp))
    result(i) = _mu * ((displacement_component == _component ? grad_phi(i) : 0) +
                       (i == displacement_component ? grad_phi(_component) : 0));

  result(_component) += _lambda * grad_phi(displacement_component);
  return result;
}
