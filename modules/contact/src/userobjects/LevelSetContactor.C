//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LevelSetContactor.h"

InputParameters
LevelSetContactor::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription("Base class for a rigid contactor described implicitly by a "
                             "signed-distance (level-set) function.");
  params.addCoupledVar("offset_variable",
                       "Optional Scalar variable representing the rigid body's translation along "
                       "`load_direction`.  When present, every query point x is evaluated as "
                       "x - s * load_direction, so a companion RigidBodyLoadControl can drive s "
                       "under a prescribed external force.");
  params.addParam<Point>(
      "load_direction",
      "Unit vector giving the direction the rigid body moves when its offset scalar "
      "grows positive.  Required if `offset_variable` is set.");
  return params;
}

LevelSetContactor::LevelSetContactor(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _offset_scalar_value(nullptr),
    _offset_scalar_number(libMesh::invalid_uint),
    _load_direction()
{
  if (isCoupledScalar("offset_variable"))
  {
    _offset_scalar_value = &coupledScalarValue("offset_variable");
    _offset_scalar_number = coupledScalar("offset_variable");

    if (!isParamValid("load_direction"))
      paramError("load_direction", "Must be set when `offset_variable` is set.");
    _load_direction = getParam<Point>("load_direction");
    const Real n = _load_direction.norm();
    if (n < TOLERANCE)
      paramError("load_direction", "Must be a nonzero vector.");
    _load_direction /= n; // normalize so callers can treat as a unit vector
  }
  else if (isParamValid("load_direction"))
    paramError("load_direction", "Only meaningful when `offset_variable` is set.");
}

Real
LevelSetContactor::offset() const
{
  return _offset_scalar_number == libMesh::invalid_uint ? 0.0 : (*_offset_scalar_value)[0];
}
