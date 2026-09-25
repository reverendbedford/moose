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
      "Calculates the full deformation gradient from one to three displacement variables.");
  params.addRequiredParam<std::vector<NonlinearVariableName>>(
      "displacements", "The displacements to use to calculate the deformation gradient.");
  return params;
}

TorchDeformationGradient::TorchDeformationGradient(const InputParameters & parameters)
  : TorchPreKernel(parameters)
{
  const auto & disp_vars = getParam<std::vector<NonlinearVariableName>>("displacements");
  if (disp_vars.size() < 1 || disp_vars.size() > 3)
    paramError("displacements",
               "TorchDeformationGradient requires 1 to 3 displacement variables, got ",
               disp_vars.size(),
               ".");

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

  _output = at::stack({dux, duy, duz}, -2) + at::eye(3, dux.options());
}

#endif // NEML2_ENABLED
