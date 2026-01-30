//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

// Forward Declarations
class Function;

/**
 * Just calculate properties for output...
 */
class FluxOutputMaterial : public Material
{
public:
  static InputParameters validParams();

  FluxOutputMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  const VariableGradient & _grad_temp;

  std::vector<const Function *> _targets;
  const ADMaterialProperty<Real> & _thermal_conductivity;

  ADMaterialProperty<RealVectorValue> & _total_temperature_gradient;
  ADMaterialProperty<RealVectorValue> & _flux;
};