//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SphereContactor.h"

registerMooseObject("ContactApp", SphereContactor);

InputParameters
SphereContactor::validParams()
{
  InputParameters params = LevelSetContactor::validParams();
  params.addClassDescription(
      "Rigid spherical contactor with analytic signed distance and gradient.");
  params.addRequiredParam<Point>("center", "Center of the rigid sphere.");
  params.addRequiredRangeCheckedParam<Real>(
      "radius", "radius > 0", "Radius of the rigid sphere.");
  return params;
}

SphereContactor::SphereContactor(const InputParameters & parameters)
  : LevelSetContactor(parameters),
    _center(getParam<Point>("center")),
    _radius(getParam<Real>("radius"))
{
}

Real
SphereContactor::signedDistance(const Point & x) const
{
  return (x - _center).norm() - _radius;
}

RealVectorValue
SphereContactor::normal(const Point & x) const
{
  const RealVectorValue r = x - _center;
  const Real rn = r.norm();
  // Fall back to +z at the (measure-zero) center to keep the field well-defined
  // rather than dividing by zero. Contact quadrature points never sit exactly
  // at the sphere center in any physical setup.
  if (rn == 0.0)
    return RealVectorValue(0.0, 0.0, 1.0);
  return r / rn;
}

RealTensorValue
SphereContactor::hessian(const Point & x) const
{
  const RealVectorValue r = x - _center;
  const Real rn = r.norm();
  if (rn == 0.0)
    return RealTensorValue();

  const RealVectorValue n = r / rn;
  // H = (I - n n^T) / |r|
  RealTensorValue H;
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      H(i, j) = ((i == j ? 1.0 : 0.0) - n(i) * n(j)) / rn;
  return H;
}
