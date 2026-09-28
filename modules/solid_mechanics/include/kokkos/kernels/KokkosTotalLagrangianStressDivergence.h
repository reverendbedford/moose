//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "KokkosKernelGrad.h"
#include "KokkosHeader.h"

/**
 * Kokkos total-Lagrangian divergence of a first Piola-Kirchhoff stress.
 *
 * When `stabilize_strain = false` (default), the kernel is a plain KernelGrad and the assembled
 * residual and Jacobian use only per-qp material properties. This preserves the original
 * unstabilized behavior bit-for-bit.
 *
 * When `stabilize_strain = true`, the kernel expects the upstream `KokkosComputeFbarAverage`
 * material to publish element-constant F_avg_invT and the stress material
 * (`KokkosComputeLagrangianStressCustomPK2` with matching `stabilize_strain = true`) to publish
 * per-qp `A(qp)` and `gamma_over_3(qp)`. Under total-mode F-bar the non-local Jacobian factors
 * to a rank-one perturbation in F_avg:
 *   delta_PK1(qp)[i, J] = A(qp)[i, J] * (gamma/3)(qp) * (F_avg_invT(beta, :) . avg_grad_phi_j),
 * where `beta` is the trial-variable displacement component (the diagonal case uses
 * `beta = _component`; off-diagonal uses the component whose variable equals `datum.jvar()`).
 * The kernel overrides `computeJacobianInternal` and `computeOffDiagJacobianInternal` to add
 * this term to the local per-qp Jacobian that is already computed by the material's F-bar
 * chain. The residual is unchanged from the unstabilized version: total-mode F-bar's residual
 * add-on lives only in `F_bar_mode = incremental`, which the Kokkos path does not implement.
 *
 * Restrictions: requires exactly 3 displacement components on a 3D mesh, and elements with
 * <= MAX_CACHED_DOF trial DOFs per variable (checked at first launch).
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

  /// Override the base's parallel Jacobian body so the F-bar branch can run an element-level
  /// pre-pass (avg_grad_phi[j]) before the per-(j, qp) inner loop. When `_stabilize_strain` is
  /// false, this delegates to the same qp-local formulation as the base KernelGrad.
  template <typename Derived>
  KOKKOS_FUNCTION void computeJacobianInternal(const Derived & kernel, AssemblyDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION void computeOffDiagJacobianInternal(const Derived & kernel,
                                                       AssemblyDatum & datum) const;

private:
  KOKKOS_FUNCTION Real3 jacobian(const unsigned int displacement_component,
                                  const Real3 & grad_phi,
                                  const unsigned int qp,
                                  AssemblyDatum & datum) const;

  /// Common F-bar-aware Jacobian body used by both computeJacobianInternal and its off-diag
  /// counterpart. `trial_component` selects the trial displacement component (equal to
  /// `_component` for the diagonal case). Runs the base's qp-local Jacobian body when
  /// `_stabilize_strain_device` is false so the unstabilized path is untouched.
  template <typename Derived>
  KOKKOS_FUNCTION void computeJacobianBody(const Derived & kernel,
                                            AssemblyDatum & datum,
                                            const unsigned int trial_component) const;

  const unsigned int _component;
  const unsigned int _ndisp;
  Moose::Kokkos::MaterialProperty<Real, 2> _pk1_stress;
  Moose::Kokkos::MaterialProperty<Real, 4> _dpk1_d_grad_u;
  Moose::Kokkos::Array<unsigned int> _displacement_var_ids;

  /// If true, add the cross-QP non-local F-bar Jacobian contribution.
  bool _stabilize_strain;
  Moose::Kokkos::Scalar<bool> _stabilize_strain_device;
  /// F-bar auxiliary material properties consumed only when `_stabilize_strain` is true.
  Moose::Kokkos::MaterialProperty<Real, 2> _F_avg_invT;
  Moose::Kokkos::MaterialProperty<Real, 2> _A;
  Moose::Kokkos::MaterialProperty<Real, 0> _gamma_over_3;
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

template <typename Derived>
KOKKOS_FUNCTION void
KokkosTotalLagrangianStressDivergence::computeJacobianInternal(const Derived & kernel,
                                                               AssemblyDatum & datum) const
{
  computeJacobianBody<Derived>(kernel, datum, _component);
}

template <typename Derived>
KOKKOS_FUNCTION void
KokkosTotalLagrangianStressDivergence::computeOffDiagJacobianInternal(const Derived & kernel,
                                                                      AssemblyDatum & datum) const
{
  // Look up the trial displacement component from the coupled variable number.
  unsigned int trial_component = _ndisp; // invalid sentinel
  for (unsigned int c = 0; c < _ndisp; ++c)
    if (_displacement_var_ids[c] == datum.jvar())
    {
      trial_component = c;
      break;
    }
  // If jvar is not one of the displacements (e.g. temperature coupling), fall through to the
  // base's qp-local off-diag body -- the F-bar term does not apply.
  if (trial_component == _ndisp)
  {
    Moose::Kokkos::KernelGrad::computeOffDiagJacobianInternal<Derived>(kernel, datum);
    return;
  }
  computeJacobianBody<Derived>(kernel, datum, trial_component);
}

template <typename Derived>
KOKKOS_FUNCTION void
KokkosTotalLagrangianStressDivergence::computeJacobianBody(const Derived & kernel,
                                                           AssemblyDatum & datum,
                                                           const unsigned int trial_component)
    const
{
  // Fast path: no F-bar, defer to the base's implementation. The base's static dispatch on
  // `Derived` still routes to our precompute* methods above.
  if (!_stabilize_strain_device)
  {
    if (trial_component == _component)
      Moose::Kokkos::KernelGrad::computeJacobianInternal<Derived>(kernel, datum);
    else
      Moose::Kokkos::KernelGrad::computeOffDiagJacobianInternal<Derived>(kernel, datum);
    return;
  }

  // ---- F-bar element pre-pass ----
  // Guard the local avg_grad_phi[] buffer size. HEX27 (27 DOFs / variable) fits in
  // MAX_CACHED_DOF = 30. Any higher-order 3D FE space would overflow this local storage.
  const unsigned int n_jdofs = datum.n_jdofs();
  if (n_jdofs > Moose::Kokkos::MAX_CACHED_DOF)
    ::Kokkos::abort(
        "KokkosTotalLagrangianStressDivergence stabilize_strain path: element has more trial "
        "DOFs per variable than MAX_CACHED_DOF; F-bar element pre-pass buffer would overflow.");

  Real3 avg_grad_phi[Moose::Kokkos::MAX_CACHED_DOF];
  for (unsigned int j = 0; j < n_jdofs; ++j)
    avg_grad_phi[j] = Real3(0);
  Real vol = 0.0;
  for (unsigned int qp = 0; qp < datum.n_qps(); ++qp)
  {
    const Real w = datum.JxW(qp);
    vol += w;
    for (unsigned int j = 0; j < n_jdofs; ++j)
    {
      const Real3 g = _grad_phi(datum, j, qp);
      avg_grad_phi[j](0) += w * g(0);
      avg_grad_phi[j](1) += w * g(1);
      avg_grad_phi[j](2) += w * g(2);
    }
  }
  const Real inv_vol = 1.0 / vol;
  for (unsigned int j = 0; j < n_jdofs; ++j)
  {
    avg_grad_phi[j](0) *= inv_vol;
    avg_grad_phi[j](1) *= inv_vol;
    avg_grad_phi[j](2) *= inv_vol;
  }

  // F_avg_invT is element-constant; read at qp = 0. Cache the trial component's row.
  const auto F_avg_invT_prop = _F_avg_invT(datum, 0);
  const Real3 F_avg_invT_row(F_avg_invT_prop(trial_component, 0),
                             F_avg_invT_prop(trial_component, 1),
                             F_avg_invT_prop(trial_component, 2));

  Moose::Kokkos::ResidualObject::computeJacobianInternal(
      datum,
      [&](Real * local_ke, const unsigned int ib, const unsigned int ie, const unsigned int j)
      {
        // scalar_contr = F_avg_invT[trial_component, :] . avg_grad_phi[j], per (element, j).
        // Column-invariant across qp, so hoisted out of the qp loop.
        const Real scalar_contr = F_avg_invT_row(0) * avg_grad_phi[j](0) +
                                  F_avg_invT_row(1) * avg_grad_phi[j](1) +
                                  F_avg_invT_row(2) * avg_grad_phi[j](2);

        for (unsigned int qp = 0; qp < datum.n_qps(); ++qp)
        {
          // Local qp-wise Jacobian: existing precomputeQpJacobian (which already includes the
          // F-bar dF_stab/dF_ust chain, computed by the stress material in Phase 2).
          Real3 val =
              (trial_component == _component)
                  ? kernel.template precomputeQpJacobian<Derived>(j, qp, datum)
                  : jacobian(trial_component, _grad_phi(datum, j, qp), qp, datum);

          // Non-local F-bar contribution: delta_PK1(alpha, :) = A(alpha, :) * (gamma/3) * scalar_contr
          const auto A_prop = _A(datum, qp);
          const Real g3 = _gamma_over_3(datum, qp);
          const Real coeff = g3 * scalar_contr;
          val(0) += coeff * A_prop(_component, 0);
          val(1) += coeff * A_prop(_component, 1);
          val(2) += coeff * A_prop(_component, 2);

          const Real3 v = datum.J(qp).transpose() * (datum.JxW(qp) * val);
          for (unsigned int i = ib; i < ie; ++i)
            local_ke[i] += v * _grad_test.reference(datum, i, qp);
        }
      });
}
