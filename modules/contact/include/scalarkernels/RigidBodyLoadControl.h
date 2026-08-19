//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "NodalScalarKernel.h"

class Function;
class LevelSetContactor;
class NodalArea;

/**
 * Load-control equation for a rigid-body contactor.
 *
 * Solves the Scalar residual
 *
 *   R_s = F(t) - Sum_i w_i * lambda_i * (n_i . load_direction) = 0
 *
 * where the sum is over LM nodes on the contact sideset's lower-d block,
 * w_i is the tributary nodal area (from a companion NodalArea UO), lambda_i
 * is the nodal Lagrange multiplier (contact pressure), n_i is the outward
 * contactor normal at the deformed node position, and `load_direction` is
 * the unit vector of the applied force.  The scalar unknown `s` is the
 * rigid body's translation along `load_direction`, and the contactor
 * transforms all query points by `x - s * load_direction`.
 *
 * The Jacobian assembles both the scalar-row-times-field-column blocks
 * (d R_s / d lambda_i, d R_s / d disp_k(i)) and the transpose block
 * d R_lambda_i / d s (which the NCP kernel cannot fill because MOOSE's
 * NodalKernel base class has no scalar off-diagonal hook).
 */
class RigidBodyLoadControl : public NodalScalarKernel
{
public:
  static InputParameters validParams();
  RigidBodyLoadControl(const InputParameters & parameters);

  virtual void computeResidual() override;
  virtual void computeJacobian() override;

private:
  /// Deformed position of the k-th LM node (undeformed node + displacement).
  Point deformedNode(std::size_t k) const;

  const Function & _force;
  const LevelSetContactor & _contactor;
  const NodalArea & _nodal_area;
  const Real _c;
  const Real _spring_stiffness;

  const unsigned int _lm_var_num;
  const VariableValue & _lambda; ///< Per-node LM values after reinitNodes.

  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var_num;
  std::vector<const VariableValue *> _disp;
};
