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

/**
 * Rigid spherical contactor.
 *
 * Signed distance:  \f$ g_{LS}(x) = \| x - c \| - R \f$.
 * Outward normal:   \f$ n(x) = (x - c) / \| x - c \| \f$.
 * Hessian:          \f$ H(x) = (I - n \otimes n) / \| x - c \| \f$.
 *
 * With this sign convention \f$ g_{LS} > 0 \f$ outside the sphere (open gap).
 */
class SphereContactor : public LevelSetContactor
{
public:
  static InputParameters validParams();

  SphereContactor(const InputParameters & parameters);

  virtual Real signedDistance(const Point & x) const override;
  virtual RealVectorValue normal(const Point & x) const override;
  virtual RealTensorValue hessian(const Point & x) const override;

protected:
  /// Sphere center.
  const Point _center;
  /// Sphere radius.
  const Real _radius;
};
