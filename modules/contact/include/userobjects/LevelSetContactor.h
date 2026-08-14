//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralUserObject.h"

#include "libmesh/point.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"

/**
 * Base class for a rigid contactor described implicitly by a signed-distance
 * (level-set) function \f$ g_{LS}(x) \f$.
 *
 * The convention is \f$ g_{LS}(x) > 0 \f$ outside the rigid body (open gap) and
 * \f$ g_{LS}(x) < 0 \f$ inside it (penetration). The outward normal from the
 * rigid body at \f$ x \f$ is \f$ \nabla g_{LS}(x) \f$ (normalized when needed).
 *
 * A derived class implements the analytic level set, its gradient, and (for
 * exact algorithmic tangents in the large-deformation case) its Hessian.
 */
class LevelSetContactor : public GeneralUserObject
{
public:
  static InputParameters validParams();

  LevelSetContactor(const InputParameters & parameters);

  /// Analytic signed distance evaluated at a spatial point.
  virtual Real signedDistance(const Point & x) const = 0;

  /// Analytic gradient of the signed distance (unit-normalized outward normal
  /// for a true SDF).
  virtual RealVectorValue normal(const Point & x) const = 0;

  /// Analytic Hessian of the signed distance. Default is zero (correct for
  /// planes; overridden by curved contactors such as spheres).
  virtual RealTensorValue hessian(const Point & /*x*/) const { return RealTensorValue(); }

  /// GeneralUserObject hooks are no-ops: the contactor is a pure function of
  /// position.
  virtual void initialize() override final {}
  virtual void execute() override final {}
  virtual void finalize() override final {}
};
