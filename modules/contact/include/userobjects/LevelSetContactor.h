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
 * (level-set) function g_LS(x). Sign convention: g_LS > 0 outside the rigid
 * body (open gap), g_LS < 0 inside it. Outward normal n = grad(g_LS).
 */
class LevelSetContactor : public GeneralUserObject
{
public:
  static InputParameters validParams();
  LevelSetContactor(const InputParameters & parameters);

  /// Bundled result of all three quantities at a query point.  Callers that
  /// need more than one component at the same point should prefer queryAt()
  /// over the per-quantity accessors, because concrete contactors (notably
  /// SurfaceMeshContactor) can share expensive per-point work (KDTree lookup,
  /// point-in-solid classification) across the three quantities.
  struct Query
  {
    Real gap;
    RealVectorValue normal;
    RealTensorValue hessian;
  };

  virtual Real signedDistance(const Point & x) const = 0;
  virtual RealVectorValue normal(const Point & x) const = 0;
  virtual RealTensorValue hessian(const Point &) const { return RealTensorValue(); }

  /// Default: forward to the three virtual accessors.  Concrete contactors
  /// that can share per-point work should override this and, where practical,
  /// implement signedDistance()/normal()/hessian() as one-line wrappers around
  /// queryAt() so callers cannot accidentally pick the slow path.
  virtual Query queryAt(const Point & x) const
  {
    return {signedDistance(x), normal(x), hessian(x)};
  }

  virtual void initialize() override final {}
  virtual void execute() override final {}
  virtual void finalize() override final {}
};
