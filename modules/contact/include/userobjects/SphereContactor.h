//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
#include "LevelSetContactor.h"

class SphereContactor : public LevelSetContactor
{
public:
  static InputParameters validParams();
  SphereContactor(const InputParameters &);
  virtual Real signedDistance(const Point &) const override;
  virtual RealVectorValue normal(const Point &) const override;
  virtual RealTensorValue hessian(const Point &) const override;

protected:
  const Point _center;
  const Real _radius;
};
