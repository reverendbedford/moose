//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#ifdef NEML2_ENABLED

#include "TorchPreKernel.h"

/**
 * Computes the full deformation gradient from displacement gradients for use as a NEML2 input.
 *
 * When `stabilize_strain = true`, applies the total-mode multiplicative F-bar volumetric locking
 * correction before sending the deformation gradient to NEML2. Specifically, at each element the
 * unstabilized F_ust = I + grad(u) is averaged over the element's quadrature points using
 * TorchAssembly's reference-configuration JxWxT weights,
 *
 *   F_avg = (sum_qp F_ust[qp] * JxWxT[qp]) / (sum_qp JxWxT[qp]),
 *
 * a per-qp scaling gamma = cbrt(det(F_avg) / det(F_ust[qp])) is formed, and NEML2 receives
 *
 *   F_stab[qp] = gamma[qp] * F_ust[qp].
 *
 * This matches the F_bar_mode = "total" branch of the host ComputeLagrangianStrain calculator,
 * which is what `Physics/SolidMechanics/QuasiStatic volumetric_locking_correction = true` selects
 * by default (see QuasiStaticSolidMechanicsPhysicsBase's default
 * `volumetric_locking_correction_mode = total`). The PK1 wrap on the Kokkos material side
 * continues to use F_ust, matching the host `PK1 = F_ust * PK2(F_stab)` convention in
 * ComputeLagrangianStressCustomPK2.
 *
 * The unstabilized default (`stabilize_strain = false`) preserves the original behavior exactly:
 * NEML2 sees F_ust directly.
 */
class TorchDeformationGradient : public TorchPreKernel
{
public:
  static InputParameters validParams();

  TorchDeformationGradient(const InputParameters & parameters);

protected:
  void forward() override;

  /// Displacement gradients (each is a batched (nelem*nqp, 3) tensor)
  ///@{
  const at::Tensor * _grad_disp_x;
  const at::Tensor * _grad_disp_y;
  const at::Tensor * _grad_disp_z;
  ///@}

  /// If true, apply the F-bar (total-mode) volumetric locking correction to the deformation
  /// gradient sent to NEML2. Off by default; preserves the original behavior.
  const bool _stabilize_strain;
};

#endif // NEML2_ENABLED
