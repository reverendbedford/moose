//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FluxOutputMaterial.h"
#include "Function.h"


registerMooseObject("HeatTransferApp", FluxOutputMaterial);
InputParameters
FluxOutputMaterial::validParams()
{
    InputParameters params = Material::validParams();

    params.addCoupledVar("temp", "coupled temperature");
    params.addRequiredParam<std::vector<FunctionName>>(
        "targets", "functions giving the targets to hit for constraint types that are not none.");

    params.addParam<MaterialPropertyName>("thermal_conductivity",
                                            "thermal_conductivity",
                                            "the name of the thermal conductivity material property");

    params.addClassDescription("Calculate thermal material values for output");
    return params;
}

FluxOutputMaterial::FluxOutputMaterial(const InputParameters & parameters)
  : Material(parameters),
  _grad_temp(coupledGradient("temp")),
  _thermal_conductivity(getADMaterialProperty<Real>("thermal_conductivity")),
  _total_temperature_gradient(declareADProperty<RealVectorValue>("total_temperature_gradient")),
  _flux(declareADProperty<RealVectorValue>("flux"))
{
    // Targets to hit
    const std::vector<FunctionName> & fnames = parameters.get<std::vector<FunctionName>>("targets");
    if (fnames.size() != _mesh.dimension())
        paramError("targets", "Number of target functions must equal the problem dimension.");
    for (const auto & fname : fnames)
        _targets.push_back(&this->getFunctionByName(fname));
}

void
FluxOutputMaterial::computeQpProperties()
{
    RealVectorValue imposed_gradient;
    for (unsigned int i = 0; i < _mesh.dimension(); ++i)
        imposed_gradient(i) = _targets[i]->value(_t, _q_point[_qp]);
    _total_temperature_gradient[_qp] = imposed_gradient + _grad_temp[_qp];
    _flux[_qp] = _thermal_conductivity[_qp] * _total_temperature_gradient[_qp];
}