//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RigidBodyNodalNCPKernel.h"
#include "LevelSetContactor.h"
#include "MooseMesh.h"

registerMooseObject("ContactApp", RigidBodyNodalNCPKernel);

InputParameters
RigidBodyNodalNCPKernel::validParams()
{
  InputParameters params = NodalKernel::validParams();
  params.addClassDescription("Node-wise min-NCP for rigid-body frictionless contact "
                             "against an analytic level-set contactor. Companion to "
                             "RigidBodyNormalMechanicalContact for the traction.");
  params.addRequiredParam<UserObjectName>("contactor",
                                          "LevelSetContactor supplying signedDistance / normal.");
  params.addRequiredCoupledVar("displacements",
                               "Displacement variables in order (x, y[, z]).");
  params.addRangeCheckedParam<Real>(
      "c",
      1.0,
      "c > 0",
      "NCP scaling on the gap (balances units of lambda and gap for the FB reformulation).");
  return params;
}

RigidBodyNodalNCPKernel::RigidBodyNodalNCPKernel(const InputParameters & parameters)
  : NodalKernel(parameters),
    _contactor(getUserObject<LevelSetContactor>("contactor")),
    _c(getParam<Real>("c")),
    _ndisp(coupledComponents("displacements")),
    _disp(_ndisp),
    _disp_num(_ndisp)
{
  if (_ndisp != _mesh.dimension())
    paramError("displacements",
               "Number of displacement components must match mesh dimension (",
               _mesh.dimension(),
               ").");
  for (const auto k : make_range(_ndisp))
  {
    _disp[k] = &coupledValue("displacements", k);
    _disp_num[k] = coupled("displacements", k);
  }
}

Point
RigidBodyNodalNCPKernel::deformedNode() const
{
  Point x = *_current_node;
  for (const auto k : make_range(_ndisp))
    x(k) += (*_disp[k])[_qp];
  return x;
}

Real
RigidBodyNodalNCPKernel::physicalGap() const
{
  return _contactor.signedDistance(deformedNode());
}

Real
RigidBodyNodalNCPKernel::computeQpResidual()
{
  return std::min(_u[_qp], _c * physicalGap());
}

Real
RigidBodyNodalNCPKernel::computeQpJacobian()
{
  // lambda-branch active -> R = lambda, dR/dlambda = 1.  gap-branch: dR/dlambda = 0.
  if (_u[_qp] <= _c * physicalGap())
    return 1.0;
  return 0.0;
}

Real
RigidBodyNodalNCPKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Only the gap-branch contributes off-diagonal (to the disp components at
  // the current node).  On the lambda-branch, R does not depend on disp.
  if (_c * physicalGap() >= _u[_qp])
    return 0.0;
  for (const auto k : make_range(_ndisp))
    if (jvar == _disp_num[k])
      return _c * _contactor.normal(deformedNode())(k);
  return 0.0;
}
