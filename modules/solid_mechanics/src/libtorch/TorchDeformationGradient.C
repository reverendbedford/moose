//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#ifdef NEML2_ENABLED

#include "TorchDeformationGradient.h"

registerMooseObject("SolidMechanicsApp", TorchDeformationGradient);

InputParameters
TorchDeformationGradient::validParams()
{
  InputParameters params = TorchPreKernel::validParams();
  params.addClassDescription(
      "Calculates the full deformation gradient from one to three displacement variables. "
      "Optionally applies the total-mode F-bar volumetric locking correction so the F sent to "
      "NEML2 matches the host `Physics/SolidMechanics/QuasiStatic volumetric_locking_correction = "
      "true` formulation.");
  params.addRequiredParam<std::vector<NonlinearVariableName>>(
      "displacements", "The displacements to use to calculate the deformation gradient.");
  params.addParam<bool>(
      "stabilize_strain",
      false,
      "If true, apply the total-mode multiplicative F-bar volumetric locking correction to F "
      "before it is passed to NEML2. F_stab[qp] = cbrt(det(F_avg)/det(F_ust[qp])) * F_ust[qp] "
      "where F_avg is the element average of F_ust weighted by TorchAssembly's JxWxT. Off by "
      "default preserves the original unstabilized behavior. Pair with the matching "
      "`stabilize_strain = true` option on `KokkosComputeLagrangianStressCustomPK2` and "
      "`KokkosTotalLagrangianStressDivergence` for a consistent Jacobian.");
  return params;
}

TorchDeformationGradient::TorchDeformationGradient(const InputParameters & parameters)
  : TorchPreKernel(parameters), _stabilize_strain(getParam<bool>("stabilize_strain"))
{
  const auto & disp_vars = getParam<std::vector<NonlinearVariableName>>("displacements");
  if (disp_vars.size() < 1 || disp_vars.size() > 3)
    paramError("displacements",
               "TorchDeformationGradient requires 1 to 3 displacement variables, got ",
               disp_vars.size(),
               ".");
  // The F-bar total-mode formula assumes a 3D deformation gradient (cbrt of a 3D volume ratio).
  // Rather than silently doing the wrong thing for 1D/2D, require the full 3D triple.
  if (_stabilize_strain && disp_vars.size() != 3)
    paramError("stabilize_strain",
               "The F-bar total-mode volumetric correction is only implemented for the full 3D "
               "deformation gradient (3 displacement variables). Got ",
               disp_vars.size(),
               ". Set stabilize_strain=false or provide three displacement variables.");

  _grad_disp_x = &_fe.getGradient(disp_vars[0]);
  _grad_disp_y = disp_vars.size() >= 2 ? &_fe.getGradient(disp_vars[1]) : nullptr;
  _grad_disp_z = disp_vars.size() >= 3 ? &_fe.getGradient(disp_vars[2]) : nullptr;
}

void
TorchDeformationGradient::forward()
{
  const auto & dux = *_grad_disp_x;
  const auto duy = _grad_disp_y ? *_grad_disp_y : at::zeros_like(dux);
  const auto duz = _grad_disp_z ? *_grad_disp_z : at::zeros_like(dux);

  // F_ust[batch, 3, 3] = I + [[dux], [duy], [duz]], where each row is one component's gradient.
  const auto F_ust = at::stack({dux, duy, duz}, -2) + at::eye(3, dux.options());

  if (!_stabilize_strain)
  {
    _output = F_ust;
    return;
  }

  // ---- Total-mode F-bar volumetric locking correction ----
  // Match the host formulation in ComputeLagrangianStrainBase::computeDeformationGradient:
  //   F_avg = sum_qp(F_ust * JxW * coord) / sum_qp(JxW * coord)   (per element)
  //   gamma[qp] = cbrt(det(F_avg) / det(F_ust[qp]))                (per qp)
  //   F_stab[qp] = gamma[qp] * F_ust[qp]                            (per qp)
  //
  // TorchFEInterpolation gathers gradients into the (nelem, nqp, 3) layout consumed by
  // TorchFEM::interpolate, so `at::stack(..., -2)` above already produces the (nelem, nqp, 3, 3)
  // tensor the F-bar reduction wants -- no reshape needed.
  const auto & JxWxT = _assembly.JxWxT(); // (nelem, nqp)
  if (F_ust.size(0) != JxWxT.size(0) || F_ust.size(1) != JxWxT.size(1))
    mooseError("TorchDeformationGradient: F_ust shape (",
               F_ust.size(0),
               ", ",
               F_ust.size(1),
               ", ...) does not match TorchAssembly JxWxT shape (",
               JxWxT.size(0),
               ", ",
               JxWxT.size(1),
               "). F-bar requires uniform (element, quadrature) topology across the mesh.");

  // Weighted element sum of F_ust, then divide by total weight.
  //   w has shape (nelem, nqp, 1, 1); F_avg has shape (nelem, 3, 3).
  const auto w = JxWxT.unsqueeze(-1).unsqueeze(-1);
  const auto sum_w = JxWxT.sum(/*dim=*/1); // (nelem,)
  const auto F_avg = (F_ust * w).sum(/*dim=*/1) / sum_w.unsqueeze(-1).unsqueeze(-1);

  // Determinants; det_F_ust broadcast per qp against the element-scalar det_F_avg.
  const auto det_F_ust = at::linalg_det(F_ust); // (nelem, nqp)
  const auto det_F_avg = at::linalg_det(F_avg); // (nelem,)

  // gamma = cbrt(det_F_avg / det_F_ust). torch has no direct cbrt; use the sign-preserving
  // real cube-root sign(x) * |x|^(1/3) so intermediate Newton iterates that transiently make F
  // non-physical (det <= 0 at some qp -- e.g., the aggressive linear-disp initial guess used
  // by the jacobian test) produce a finite value SNES can line-search away from, matching the
  // host's tolerance for non-physical iterates rather than aborting with a NaN. `at::pow` on a
  // fractional exponent applied to a negative base returns NaN per IEEE 754.
  const auto ratio = det_F_avg.unsqueeze(-1) / det_F_ust;                        // (nelem, nqp)
  const auto gamma = at::sign(ratio) * at::pow(at::abs(ratio), 1.0 / 3.0);       // (nelem, nqp)

  _output = gamma.unsqueeze(-1).unsqueeze(-1) * F_ust; // (nelem, nqp, 3, 3)
}

#endif // NEML2_ENABLED
