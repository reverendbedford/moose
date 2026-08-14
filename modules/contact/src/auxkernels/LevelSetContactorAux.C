//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LevelSetContactorAux.h"
#include "LevelSetContactor.h"

registerMooseObject("ContactApp", LevelSetContactorAux);

InputParameters
LevelSetContactorAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Sample a LevelSetContactor's signed distance or normal component.");
  params.addRequiredParam<UserObjectName>("contactor", "The LevelSetContactor to sample.");
  MooseEnum quantity("signed_distance normal_x normal_y normal_z", "signed_distance");
  params.addParam<MooseEnum>("quantity", quantity, "Which scalar quantity to output.");
  return params;
}

LevelSetContactorAux::LevelSetContactorAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _contactor(getUserObject<LevelSetContactor>("contactor")),
    _quantity(getParam<MooseEnum>("quantity"))
{
}

Real
LevelSetContactorAux::computeValue()
{
  const Point & x = isNodal() ? static_cast<const Point &>(*_current_node) : _q_point[_qp];

  if (_quantity == "signed_distance")
    return _contactor.signedDistance(x);

  const RealVectorValue n = _contactor.normal(x);
  if (_quantity == "normal_x")
    return n(0);
  if (_quantity == "normal_y")
    return n(1);
  return n(2);
}
