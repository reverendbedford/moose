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
 */
class TorchDeformationGradient : public TorchPreKernel
{
public:
  static InputParameters validParams();

  TorchDeformationGradient(const InputParameters & parameters);

protected:
  void forward() override;

  /// Displacement gradients
  ///@{
  const at::Tensor * _grad_disp_x;
  const at::Tensor * _grad_disp_y;
  const at::Tensor * _grad_disp_z;
  ///@}
};

#endif // NEML2_ENABLED
