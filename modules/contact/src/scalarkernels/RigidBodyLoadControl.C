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
  params.addClassDescription(
      "Load-control constraint for a rigid-body contactor: enforces "
      "F(t) = Sum_i w_i * lambda_i * (n_i . direction) so that the scalar "
      "variable driving the contactor's translation in `direction` is "
      "determined by the target contact reaction.");
  params.addRequiredParam<FunctionName>(
      "force", "Function returning the target integrated contact reaction F(t) along direction.");
  params.addRequiredParam<UserObjectName>(
      "contactor",
      "LevelSetContactor whose translation in `direction` is driven by this kernel's `variable`.");
  params.addRequiredParam<UserObjectName>(
      "nodal_area", "NodalArea UO providing tributary weights w_i on the contact sideset.");
  params.addRequiredCoupledVar(
      "lm_variable",
      "The Lagrange multiplier field variable (nodal, on the contact sideset's lower-d block).");
  params.addRequiredCoupledVar("displacements", "Displacement variables in order (x, y[, z]).");
  params.addRequiredParam<Point>(
      "direction",
      "Unit vector along which the integrated normal contact reaction is measured.  Must be "
      "aligned with a Cartesian axis and match the contactor axis driven by this kernel's "
      "`variable`.");
  params.addRangeCheckedParam<Real>(
      "c",
      1.0,
      "c > 0",
      "NCP scaling on the gap.  MUST match the `c` used by the companion "
      "RigidBodyNodalNCPKernel so the (LM_row, scalar_col) transpose Jacobian block is correct.");
  return params;
}

RigidBodyLoadControl::RigidBodyLoadControl(const InputParameters & parameters)
  : NodalScalarKernel(parameters),
    _force(getFunction("force")),
    _contactor(getUserObject<LevelSetContactor>("contactor")),
    _nodal_area(getUserObject<NodalArea>("nodal_area")),
    _direction(getParam<Point>("direction")),
    _c(getParam<Real>("c")),
    _lm_var_num(coupled("lm_variable")),
    _lambda(coupledValue("lm_variable")),
    _ndisp(coupledComponents("displacements")),
    _disp_var_num(_ndisp),
    _disp(_ndisp)
{
  const Real n = _direction.norm();
  if (n < TOLERANCE)
    paramError("direction", "Must be a nonzero vector.");
  _direction /= n;

  // Verify the direction is aligned with a Cartesian axis, and that the
  // contactor's translation input on that axis is the scalar we are the
  // kernel of.  Anything else means the user's plumbing is inconsistent.
  unsigned int axis = libMesh::invalid_uint;
  for (const auto k : {0u, 1u, 2u})
    if (std::abs(std::abs(_direction(k)) - 1.0) < TOLERANCE)
      axis = k;
  if (axis == libMesh::invalid_uint)
    paramError("direction",
               "Must be aligned with a Cartesian axis (a signed unit vector on x, y, or z).");
  const unsigned int contactor_scalar = _contactor.translationScalarNumber(axis);
  if (contactor_scalar == libMesh::invalid_uint)
    paramError("contactor",
               "The contactor's translation on axis ",
               axis,
               " must be driven by a Scalar variable (via disp_",
               std::array<const char *, 3>{"x", "y", "z"}[axis],
               "_scalar) for this kernel to close the load-control loop.");
  if (contactor_scalar != _var.number())
    paramError("variable",
               "This kernel's `variable` (number ",
               _var.number(),
               ") must match the Scalar variable driving the contactor's axis-",
               axis,
               " translation (number ",
               contactor_scalar,
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
    reaction += w * _lambda[k] * (q.normal * _direction);
  }

  prepareVectorTag(_assembly, _var.number());
  _local_re(0) = reaction - F;
  assignTaggedLocalResidual();
}

void
RigidBodyLoadControl::computeJacobian()
{
  // Precompute per-node normal projections (also decides which branch of
  // the NCP each node is on, for the transpose block).
  const auto N = _node_ids.size();
  std::vector<Real> n_dot_dir(N);
  std::vector<bool> gap_branch(N);
  for (const auto k : index_range(_node_ids))
  {
    const auto q = _contactor.queryAt(deformedNode(k));
    n_dot_dir[k] = q.normal * _direction;
    gap_branch[k] = _c * q.gap < _lambda[k];
  }

  // (scalar_row, scalar_col): dR_s/ds = 0 in this formulation (F(t) does
  // not depend on s and the reaction depends on s only indirectly through
  // the geometric term dn/ds, which is a hessian order term we drop).
  // PETSc's sparsity still needs the entry, so assemble an explicit zero.
  prepareMatrixTag(_assembly, _var.number(), _var.number());
  for (const auto i : make_range(_local_ke.m()))
    for (const auto j : make_range(_local_ke.n()))
      _local_ke(i, j) = 0.0;
  assignTaggedLocalMatrix();

  // (scalar_row, lambda_col): dR_s / dlambda_j = -w_j * (n_j . direction).
  prepareMatrixTag(_assembly, _var.number(), _lm_var_num);
  for (const auto i : make_range(_local_ke.m()))
    for (const auto j : make_range(_local_ke.n()))
      _local_ke(i, j) = 0.0;
  for (const auto k : index_range(_node_ids))
  {
    const Node * node = _mesh.getMesh().node_ptr(_node_ids[k]);
    const Real w = _nodal_area.nodalArea(node);
    _local_ke(0, k) = w * n_dot_dir[k];
  }
  assignTaggedLocalMatrix();

  // (lambda_row, scalar_col): the transpose block.  For each LM node on the
  // gap branch of min(lambda, c*g), R_lambda = c * g_LS(x - s*direction), so
  // dR_lambda / ds = c * grad(g_LS) . (-direction) = -c * (n . direction).
  // Lambda-branch nodes have R_lambda = lambda (independent of s).
  //
  // MOOSE's ScalarKernel dispatch (`addJacobianOffDiagScalar`) only fills
  // (scalar_row × field_col) blocks — it never fills (field_row × scalar_col).
  // NodalKernel likewise has no scalar off-diagonal hook.  The idiomatic
  // workaround (MortarScalarBase pattern) is direct assembly via
  // TaggingInterface::addJacobian with explicit row/column DoF indices,
  // bypassing the tagged-block dispatch.
  const auto & lm_var = _sys.getVariable(_tid, _lm_var_num);
  const auto & scalar_dofs = _var.dofIndices();
  std::vector<dof_id_type> lm_dofs(N);
  for (const auto k : index_range(_node_ids))
    lm_dofs[k] = _mesh.getMesh().node_ref(_node_ids[k]).dof_number(
        _sys.number(), _lm_var_num, /*comp=*/0);

  DenseMatrix<Real> ke_transpose(N, 1);
  for (const auto k : index_range(_node_ids))
    if (gap_branch[k])
      ke_transpose(k, 0) = -_c * n_dot_dir[k];
  addJacobian(_assembly, ke_transpose, lm_dofs, scalar_dofs, lm_var.scalingFactor());
}
