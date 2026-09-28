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
 *
 * When `stabilize_strain = false` (the default), the material operates in the original
 * unstabilized mode: `_dpk2_dF` is treated as `dPK2/dF_ust`, and PK1 = F_ust * PK2, so
 *   dPK1_local(i,J,k,L) = del_{ik} PK2(L,J) + F_ust(i,K) * dPK2/dF_ust(K,J,k,L).
 *
 * When `stabilize_strain = true`, the material assumes NEML2 saw the total-mode F-bar
 * stabilized deformation gradient F_stab (see `TorchDeformationGradient` with matching
 * `stabilize_strain = true`), so `_dpk2_dF` is now `dPK2/dF_stab`. The PK1 wrap still uses
 * F_ust (matching the host `ComputeLagrangianStressCustomPK2.C:49`
 * `_pk1_stress[_qp] = _F_ust[_qp] * _pk2[_qp]`), so:
 *   dF_stab/dF_ust = gamma * I^(4) - (gamma/3) * F_ust (x) F_ust^{-T},
 *   dPK1_local(i,J,k,L) = del_{ik} PK2(L,J)
 *                       + gamma * F_ust(i,K) * dPK2/dF_stab(K,J,k,L)
 *                       - (gamma/3) * A(i,J) * F_ust^{-T}(k,L),
 * where A(i,J) = F_ust(i,K) * dPK2/dF_stab(K,J,m,n) * F_ust(m,n) is a per-qp R2 tensor and
 * gamma = cbrt(det(F_avg)/det(F_ust)).
 *
 * A and gamma are also declared as material properties so the total-Lagrangian kernel can
 * add the non-local F-bar coupling term
 *   delta_pk1_nl(i,J) = (gamma/3) * A(i,J) * ( F_avg^{-T} : delta_F_avg ),
 * which factors as rank-1 in F_avg for total-mode F-bar. Off by default so the unstabilized
 * path incurs no extra work.
 *
 * Requires 3D Cartesian meshes with 3 displacement variables.
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

  /// If true, treat _dpk2_dF as dPK2/dF_stab and apply the F-bar chain in the Jacobian. Not
  /// declared `const` because `Moose::Kokkos::Scalar<T>` binds to a mutable reference (see the
  /// existing `_need_jacobian` pattern); the value is written only once at construction.
  bool _stabilize_strain;
  Moose::Kokkos::Scalar<bool> _stabilize_strain_device;
  /// F-bar (total-mode) element averages, only consumed when _stabilize_strain is true.
  /// Non-owned; provided by an upstream `KokkosComputeFbarAverage` with constant_on = 'ELEMENT'.
  Moose::Kokkos::MaterialProperty<Real, 2> _F_avg;
  Moose::Kokkos::MaterialProperty<Real, 2> _F_avg_invT;
  /// Per-qp F-bar auxiliaries published for the kernel's non-local Jacobian term.
  /// `_A(i,J) = F_ust(i,K) dPK2/dF_stab(K,J,m,n) F_ust(m,n)`; `_gamma_over_3 = gamma / 3`.
  /// Only meaningful when _stabilize_strain is true.
  Moose::Kokkos::MaterialProperty<Real, 2> _A;
  Moose::Kokkos::MaterialProperty<Real, 0> _gamma_over_3;
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

  // Unstabilized F_ust = I + grad(u). Same value whether or not F-bar is on (F_stab enters the
  // stress chain only through NEML2's dpk2_dF).
  Real F[3][3];
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      F[i][j] = (i == j) + _grad_displacements(datum, qp, i)(j);

  // PK1 = F_ust * PK2 in all cases; matches host `ComputeLagrangianStressCustomPK2.C:49`.
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int J = 0; J < 3; ++J)
    {
      pk1(i, J) = 0;
      for (unsigned int K = 0; K < 3; ++K)
        pk1(i, J) += F[i][K] * pk2(K, J);
    }

  if (!_need_jacobian_device)
    return;

  if (!_stabilize_strain_device)
  {
    // Original unstabilized Jacobian: dPK1(i,J)/d(grad u)(k,L) = del_{ik} PK2(L,J)
    //                                                          + F_ust(i,K) * dPK2/dF(K,J,k,L).
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int J = 0; J < 3; ++J)
        for (unsigned int k = 0; k < 3; ++k)
          for (unsigned int L = 0; L < 3; ++L)
          {
            dpk1_d_grad_u(i, J, k, L) = (i == k) * pk2(L, J);
            for (unsigned int K = 0; K < 3; ++K)
              dpk1_d_grad_u(i, J, k, L) += F[i][K] * dpk2_dF(K, J, k, L);
          }
    return;
  }

  // ---- F-bar total-mode Jacobian chain ----
  // gamma = cbrt(det(F_avg)/det(F_ust)).
  const Real detF = F[0][0] * (F[1][1] * F[2][2] - F[1][2] * F[2][1]) -
                    F[0][1] * (F[1][0] * F[2][2] - F[1][2] * F[2][0]) +
                    F[0][2] * (F[1][0] * F[2][1] - F[1][1] * F[2][0]);
  const auto F_avg = _F_avg(datum, qp); // element-constant; qp arg ignored by storage
  const Real detFavg =
      F_avg(0, 0) * (F_avg(1, 1) * F_avg(2, 2) - F_avg(1, 2) * F_avg(2, 1)) -
      F_avg(0, 1) * (F_avg(1, 0) * F_avg(2, 2) - F_avg(1, 2) * F_avg(2, 0)) +
      F_avg(0, 2) * (F_avg(1, 0) * F_avg(2, 1) - F_avg(1, 1) * F_avg(2, 0));
  // Sign-preserving real cube root so intermediate Newton iterates with det F < 0 at some qp
  // produce a finite value SNES can line-search away from, rather than propagating NaN. Matches
  // the tolerance policy of the torch-side TorchDeformationGradient F-bar.
  const Real ratio = detFavg / detF;
  const Real gamma =
      (ratio < 0.0 ? -1.0 : 1.0) * ::Kokkos::pow(::Kokkos::fabs(ratio), 1.0 / 3.0);

  // F_ust^{-T} = cofactor(F_ust) / det(F_ust). Reuse the cofactor arithmetic explicitly.
  Real Cof[3][3];
  Cof[0][0] = F[1][1] * F[2][2] - F[1][2] * F[2][1];
  Cof[0][1] = -(F[1][0] * F[2][2] - F[1][2] * F[2][0]);
  Cof[0][2] = F[1][0] * F[2][1] - F[1][1] * F[2][0];
  Cof[1][0] = -(F[0][1] * F[2][2] - F[0][2] * F[2][1]);
  Cof[1][1] = F[0][0] * F[2][2] - F[0][2] * F[2][0];
  Cof[1][2] = -(F[0][0] * F[2][1] - F[0][1] * F[2][0]);
  Cof[2][0] = F[0][1] * F[1][2] - F[0][2] * F[1][1];
  Cof[2][1] = -(F[0][0] * F[1][2] - F[0][2] * F[1][0]);
  Cof[2][2] = F[0][0] * F[1][1] - F[0][1] * F[1][0];
  const Real inv_detF = 1.0 / detF;
  Real Fust_invT[3][3];
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      Fust_invT[i][j] = Cof[i][j] * inv_detF;

  // A(i, J) = F_ust(i, K) * dPK2/dF_stab(K, J, m, n) * F_ust(m, n).
  auto A_prop = _A(datum, qp);
  Real Aij[3][3];
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int J = 0; J < 3; ++J)
    {
      Real s = 0.0;
      for (unsigned int K = 0; K < 3; ++K)
        for (unsigned int m = 0; m < 3; ++m)
          for (unsigned int n = 0; n < 3; ++n)
            s += F[i][K] * dpk2_dF(K, J, m, n) * F[m][n];
      Aij[i][J] = s;
      A_prop(i, J) = s;
    }

  // gamma / 3, published as scalar for the kernel's non-local term.
  const Real g3 = gamma / 3.0;
  _gamma_over_3(datum, qp) = g3;

  // Local Jacobian: dPK1(i,J)/d(grad u)(k,L) =
  //     del_{ik} PK2(L,J)
  //   + gamma * F_ust(i,K) * dPK2/dF_stab(K,J,k,L)
  //   - (gamma/3) * A(i,J) * F_ust^{-T}(k,L)
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int J = 0; J < 3; ++J)
      for (unsigned int k = 0; k < 3; ++k)
        for (unsigned int L = 0; L < 3; ++L)
        {
          Real v = (i == k) * pk2(L, J);
          for (unsigned int K = 0; K < 3; ++K)
            v += gamma * F[i][K] * dpk2_dF(K, J, k, L);
          v -= g3 * Aij[i][J] * Fust_invT[k][L];
          dpk1_d_grad_u(i, J, k, L) = v;
        }
}
