//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADHeatConduction.h"

class ADHeatConductionHomogenized : public ADHeatConduction
{
public:
  static InputParameters validParams();

  ADHeatConductionHomogenized(const InputParameters & parameters);

protected:
  virtual ADRealVectorValue precomputeQpResidual() override;

  std::vector<const Function*> _targets;
};
