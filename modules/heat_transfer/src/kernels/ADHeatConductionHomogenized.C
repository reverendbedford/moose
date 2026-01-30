//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADHeatConductionHomogenized.h"

#include "Function.h"

registerMooseObject("HeatTransferApp", ADHeatConductionHomogenized);

InputParameters
ADHeatConductionHomogenized::validParams()
{
  InputParameters params = ADHeatConduction::validParams();
  params.addRequiredParam<std::vector<FunctionName>>(
      "targets", "Imposed macroscale thermal gradient in each direction.");
  return params;
}

ADHeatConductionHomogenized::ADHeatConductionHomogenized(const InputParameters & parameters)
  : ADHeatConduction(parameters)
{
    // Targets to hit
  const std::vector<FunctionName> & fnames = parameters.get<std::vector<FunctionName>>("targets");
  if (fnames.size() != _mesh.dimension())
    paramError("targets", "Number of target functions must equal the problem dimension.");
  for (const auto & fname : fnames)
    _targets.push_back(&this->getFunctionByName(fname));
}

ADRealVectorValue
ADHeatConductionHomogenized::precomputeQpResidual()
{
  RealVectorValue imposed_gradient;
  for (unsigned int i = 0; i < _mesh.dimension(); ++i)
    imposed_gradient(i) = _targets[i]->value(_t, _q_point[_qp]);
  return _thermal_conductivity[_qp] * imposed_gradient +  ADHeatConduction::precomputeQpResidual();
}
