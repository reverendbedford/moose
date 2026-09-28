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
 * Element-constant Kokkos material publishing the F-bar (total-mode) element averages needed by
 * KokkosComputeLagrangianStressCustomPK2 and KokkosTotalLagrangianStressDivergence to build a
 * consistent stabilized Jacobian.
 *
 * Publishes (as `constant_on = element` properties):
 *   * F_avg      : the reference-configuration element average of F_ust = I + grad(u), weighted
 *                  by JxW * coord (Kokkos Assembly's `getJxW`, which already includes the
 *                  coordinate-transform factor).
 *   * F_avg_invT : (F_avg)^{-T}, needed by the kernel's non-local Jacobian scalar contraction.
 *
 * Matches the F_bar_mode = "total" branch of the host ComputeLagrangianStrainBase (see
 * `ComputeLagrangianStrainBase.C:830-887`), which is what
 * `Physics/SolidMechanics/QuasiStatic volumetric_locking_correction = true` selects by default.
 *
 * Restrictions:
 *   * 3D Cartesian only (three displacement variables, mesh dimension 3). The cube-root ratio in
 *     the total-mode F-bar formula assumes a 3D volume ratio.
 */
class KokkosComputeFbarAverage : public Moose::Kokkos::Material
{
public:
  static InputParameters validParams();

  KokkosComputeFbarAverage(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION void computeQpProperties(const unsigned int qp, Datum & datum) const;

private:
  const unsigned int _ndisp;
  const Moose::Kokkos::VariableGradient _grad_displacements;
  /// Element-averaged deformation gradient F_avg.
  Moose::Kokkos::MaterialProperty<Real, 2> _F_avg;
  /// (F_avg)^{-T}. Consumed by the TL kernel's non-local Jacobian.
  Moose::Kokkos::MaterialProperty<Real, 2> _F_avg_invT;
};

template <typename Derived>
KOKKOS_FUNCTION void
KokkosComputeFbarAverage::computeQpProperties(const unsigned int /*qp*/, Datum & datum) const
{
  // Material framework calls this once per element when `constant_on = element` is set (see
  // `Material::operator()(ElementCompute, ...)` in framework/include/kokkos/materials/KokkosMaterial.h).
  // We must therefore loop over every qp of the current element internally to build F_avg.
  auto F_avg = _F_avg(datum, 0);
  auto F_avg_invT = _F_avg_invT(datum, 0);

  // Accumulator for the volume-weighted sum of F_ust and the total volume; kept as scalars until
  // the final divide, so we allocate no rank-two work tensor.
  Real F_sum[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  Real vol = 0.0;

  const unsigned int nqp = datum.n_qps();
  for (unsigned int p = 0; p < nqp; ++p)
  {
    const Real w = datum.JxW(p);
    vol += w;
    // F_ust[i, j] = delta_{ij} + (grad u_i)_j
    for (unsigned int i = 0; i < 3; ++i)
    {
      const auto grad_i = _grad_displacements(datum, p, i);
      for (unsigned int j = 0; j < 3; ++j)
        F_sum[i][j] += w * ((i == j ? 1.0 : 0.0) + grad_i(j));
    }
  }

  const Real inv_vol = 1.0 / vol;
  Real Fa[3][3];
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
    {
      Fa[i][j] = F_sum[i][j] * inv_vol;
      F_avg(i, j) = Fa[i][j];
    }

  // Compute F_avg^{-T} via the 3x3 cofactor / determinant formula. On-device inversion; avoids
  // needing to call into any linear-algebra helper. Cofactor matrix C: C_{ij} is the (i, j)
  // cofactor. inv(F)_{ij} = C_{ji} / det. inv(F)^T_{ij} = C_{ij} / det.
  Real C[3][3];
  C[0][0] = Fa[1][1] * Fa[2][2] - Fa[1][2] * Fa[2][1];
  C[0][1] = -(Fa[1][0] * Fa[2][2] - Fa[1][2] * Fa[2][0]);
  C[0][2] = Fa[1][0] * Fa[2][1] - Fa[1][1] * Fa[2][0];
  C[1][0] = -(Fa[0][1] * Fa[2][2] - Fa[0][2] * Fa[2][1]);
  C[1][1] = Fa[0][0] * Fa[2][2] - Fa[0][2] * Fa[2][0];
  C[1][2] = -(Fa[0][0] * Fa[2][1] - Fa[0][1] * Fa[2][0]);
  C[2][0] = Fa[0][1] * Fa[1][2] - Fa[0][2] * Fa[1][1];
  C[2][1] = -(Fa[0][0] * Fa[1][2] - Fa[0][2] * Fa[1][0]);
  C[2][2] = Fa[0][0] * Fa[1][1] - Fa[0][1] * Fa[1][0];
  const Real detF = Fa[0][0] * C[0][0] + Fa[0][1] * C[0][1] + Fa[0][2] * C[0][2];
  const Real inv_det = 1.0 / detF;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      F_avg_invT(i, j) = C[i][j] * inv_det;
}
