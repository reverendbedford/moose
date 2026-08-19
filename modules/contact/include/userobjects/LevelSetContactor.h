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
 *
 * Optional load-control: users can attach a Scalar variable `offset_variable`
 * and a unit `load_direction`; every public query point x is transformed to
 * `x - s * load_direction` before being handed to the concrete geometry.
 * This lets a companion RigidBodyLoadControl scalar kernel drive the rigid
 * body's translation along the load direction under a prescribed force.
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

  // Public non-virtual API.  Applies the offset transform (if any) and
  // delegates to the concrete *Raw method below.
  Real signedDistance(const Point & x) const { return signedDistanceRaw(transformed(x)); }
  RealVectorValue normal(const Point & x) const { return normalRaw(transformed(x)); }
  RealTensorValue hessian(const Point & x) const { return hessianRaw(transformed(x)); }
  Query queryAt(const Point & x) const { return queryAtRaw(transformed(x)); }

  /// Whether an offset scalar has been attached to this contactor.
  bool hasOffset() const { return _offset_scalar_number != libMesh::invalid_uint; }
  /// The (unit) direction the rigid body translates in when its offset scalar
  /// grows positive.  Only meaningful when hasOffset() is true.
  const Point & loadDirection() const { return _load_direction; }
  /// The Scalar variable's number in the nonlinear system, for coupling.
  unsigned int offsetVariableNumber() const { return _offset_scalar_number; }
  /// Current value of the offset scalar (i.e. the rigid body's translation
  /// in the load direction).  Returns 0 when no offset scalar is attached.
  Real offset() const;

  virtual void initialize() override final {}
  virtual void execute() override final {}
  virtual void finalize() override final {}

protected:
  /// Raw geometry hooks — concrete contactors implement these; they see the
  /// query point in the contactor's own (untransformed) frame.
  virtual Real signedDistanceRaw(const Point & x) const = 0;
  virtual RealVectorValue normalRaw(const Point & x) const = 0;
  virtual RealTensorValue hessianRaw(const Point &) const { return RealTensorValue(); }
  /// Default: forward to the three raw accessors.  Concrete contactors that
  /// can share per-point work should override this.
  virtual Query queryAtRaw(const Point & x) const
  {
    return {signedDistanceRaw(x), normalRaw(x), hessianRaw(x)};
  }

private:
  /// x' = x - s * load_direction (identity if no offset scalar is attached).
  Point transformed(const Point & x) const
  {
    return _offset_scalar_number == libMesh::invalid_uint
               ? x
               : Point(x - (*_offset_scalar_value)[0] * _load_direction);
  }

  const VariableValue * _offset_scalar_value;
  unsigned int _offset_scalar_number;
  Point _load_direction;
};
