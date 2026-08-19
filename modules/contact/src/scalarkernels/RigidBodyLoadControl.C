//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RigidBodyLoadControl.h"

#include "Assembly.h"
#include "Function.h"
#include "LevelSetContactor.h"
#include "MooseMesh.h"
#include "MooseVariableScalar.h"
#include "NodalArea.h"
#include "SystemBase.h"

registerMooseObject("ContactApp", RigidBodyLoadControl);

InputParameters
RigidBodyLoadControl::validParams()
{
  InputParameters params = NodalScalarKernel::validParams();
  params.addClassDescription("Load-control equation for a rigid-body contactor.  Enforces "
                             "F(t) = Sum_i w_i * lambda_i * (n_i . load_direction) at each "
                             "converged Newton state, where the scalar variable is the rigid "
                             "body's translation along `load_direction`.");
  params.addRequiredParam<FunctionName>("force",
                                        "Function returning the applied force magnitude F(t) "
                                        "along the contactor's load_direction.");
  params.addRequiredParam<UserObjectName>(
      "contactor",
      "LevelSetContactor.  Must have offset_variable = this ScalarKernel's `variable`.");
  params.addRequiredParam<UserObjectName>(
      "nodal_area", "NodalArea UserObject providing tributary weights w_i on the contact sideset.");
  params.addRequiredCoupledVar("lm_variable",
                               "The Lagrange multiplier field variable (nodal, on the sideset's "
                               "lower-d block).");
  params.addRequiredCoupledVar("displacements", "Displacement variables in order (x, y[, z]).");
  params.addRangeCheckedParam<Real>(
      "c",
      1.0,
      "c > 0",
      "NCP scaling on the gap.  Must match the `c` used by the companion "
      "RigidBodyNodalNCPKernel so the assembled d R_lambda / d s block is correct.");
  params.addRangeCheckedParam<Real>(
      "spring_stiffness",
      0.0,
      "spring_stiffness >= 0",
      "Optional linear spring reacting the rigid body's offset (adds K_s*s to R_s and K_s to "
      "the (s, s) Jacobian diagonal).  Regularizes the otherwise singular scalar diagonal "
      "when the load-control equation is used off the equilibrium manifold.  Set to a small "
      "fraction of the material's tangent stiffness times a characteristic contact area.");
  return params;
}

RigidBodyLoadControl::RigidBodyLoadControl(const InputParameters & parameters)
  : NodalScalarKernel(parameters),
    _force(getFunction("force")),
    _contactor(getUserObject<LevelSetContactor>("contactor")),
    _nodal_area(getUserObject<NodalArea>("nodal_area")),
    _c(getParam<Real>("c")),
    _spring_stiffness(getParam<Real>("spring_stiffness")),
    _lm_var_num(coupled("lm_variable")),
    _lambda(coupledValue("lm_variable")),
    _ndisp(coupledComponents("displacements")),
    _disp_var_num(_ndisp),
    _disp(_ndisp)
{
  if (!_contactor.hasOffset())
    paramError("contactor",
               "Load control requires the contactor to have `offset_variable` set — "
               "otherwise there is no rigid-body DoF to solve for.");
  if (_contactor.offsetVariableNumber() != _var.number())
    paramError("contactor",
               "The contactor's `offset_variable` (variable number ",
               _contactor.offsetVariableNumber(),
               ") must match this kernel's `variable` (number ",
               _var.number(),
               ").");
  for (const auto k : make_range(_ndisp))
  {
    _disp[k] = &coupledValue("displacements", k);
    _disp_var_num[k] = coupled("displacements", k);
  }
}

Point
RigidBodyLoadControl::deformedNode(std::size_t k) const
{
  const Node & node = _mesh.getMesh().node_ref(_node_ids[k]);
  Point x = node;
  for (const auto d : make_range(_ndisp))
    x(d) += (*_disp[d])[k];
  return x;
}

void
RigidBodyLoadControl::computeResidual()
{
  const Real F = _force.value(_t, Point());
  Real reaction = 0.0;
  for (const auto k : index_range(_node_ids))
  {
    const Node * node = _mesh.getMesh().node_ptr(_node_ids[k]);
    const Real w = _nodal_area.nodalArea(node);
    const auto q = _contactor.queryAt(deformedNode(k));
    reaction += w * _lambda[k] * (q.normal * _contactor.loadDirection());
  }

  prepareVectorTag(_assembly, _var.number());
  // Physical spring: reacts the rigid body's motion, F_spring = -K*s along
  // load_direction.  Signed as -K*s in R_s = F(t) + F_spring - reaction.
  _local_re(0) = F - reaction - _spring_stiffness * _contactor.offset();
  assignTaggedLocalResidual();
}

void
RigidBodyLoadControl::computeJacobian()
{
  const auto & ld = _contactor.loadDirection();

  // (scalar_row, scalar_col): dR_s/ds = _spring_stiffness (0 by default =
  // pure load control; a nonzero spring provides a well-conditioned
  // regularization when the transpose (LM_row, s_col) block does not by
  // itself give the scalar column enough rank).
  prepareMatrixTag(_assembly, _var.number(), _var.number());
  for (const auto i : make_range(_local_ke.m()))
    for (const auto j : make_range(_local_ke.n()))
      _local_ke(i, j) = 0.0;
  _local_ke(0, 0) = -_spring_stiffness; // d(-K*s)/ds = -K
  assignTaggedLocalMatrix();

  // (scalar_row, lm_col): dR_s / dlambda_j = -w_j * (n_j . load_dir).
  // (scalar_row, disp_k_col): dR_s / d disp_k(j) is a hessian-order term
  //   (d n_j / d disp_k = H_{k,l}); left at 0 for MVP.
  // (lm_row, scalar_col): dR_lambda_j / ds = -c * (n_j . load_dir)  on gap
  //   branch, 0 on lambda branch.  Filled here because NodalKernel has no
  //   scalar off-diagonal hook.
  const auto N = _node_ids.size();

  // Precompute per-node reaction contributions and branch states so both
  // block fills below share the same query results.
  std::vector<Real> n_dot_ld(N);
  std::vector<bool> gap_branch_active(N);
  for (const auto k : index_range(_node_ids))
  {
    const auto q = _contactor.queryAt(deformedNode(k));
    n_dot_ld[k] = q.normal * ld;
    gap_branch_active[k] = _c * q.gap < _lambda[k];
  }

  // (scalar_row, lm_col)
  prepareMatrixTag(_assembly, _var.number(), _lm_var_num);
  for (const auto i : make_range(_local_ke.m()))
    for (const auto j : make_range(_local_ke.n()))
      _local_ke(i, j) = 0.0;
  for (const auto k : index_range(_node_ids))
  {
    const Node * node = _mesh.getMesh().node_ptr(_node_ids[k]);
    const Real w = _nodal_area.nodalArea(node);
    _local_ke(0, k) = -w * n_dot_ld[k];
  }
  assignTaggedLocalMatrix();

  // (lm_row, scalar_col).  Fills the entries the NCP kernel can't reach
  // (NodalKernel has no scalar off-diagonal hook).  Without these entries
  // the scalar column of the Jacobian is empty and δs is under-determined,
  // making the linear system rank-deficient.
  prepareMatrixTag(_assembly, _lm_var_num, _var.number());
  for (const auto i : make_range(_local_ke.m()))
    for (const auto j : make_range(_local_ke.n()))
      _local_ke(i, j) = 0.0;
  for (const auto k : index_range(_node_ids))
    if (gap_branch_active[k])
      _local_ke(k, 0) = -_c * n_dot_ld[k];
  assignTaggedLocalMatrix();
}
