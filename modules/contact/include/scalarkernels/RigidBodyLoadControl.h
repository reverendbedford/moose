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
 * Load-control constraint for a rigid-body contactor.
 *
 * The scalar variable `s` is the contactor's translation along `direction`
 * (via the contactor's `disp_x/y/z_scalar` inputs).  Newton adjusts `s` so
 * the integrated normal contact reaction equals a prescribed target:
 *
 *   R_s = Sum_i w_i * lambda_i * (n_i . direction) - F(t) = 0
 *
 * where the sum is over LM nodes on the contact sideset, w_i is the
 * tributary nodal area (from a companion NodalArea UO), lambda_i is the
 * nodal contact pressure, and n_i is the contactor's outward normal at
 * the deformed node position.
 *
 * The Jacobian assembles:
 *  * (scalar_row, lambda_col) : -w_j * (n_j . direction)   — direct.
 *  * (lambda_row, scalar_col) : -c * (n_j . direction) on gap-branch nodes,
 *    0 elsewhere — transpose block that MOOSE's NodalKernel framework can't
 *    fill (NodalKernel has no scalar off-diagonal hook), so this kernel
 *    fills it directly.  The `c` parameter MUST match the companion
 *    RigidBodyNodalNCPKernel's `c` for the Jacobian to be consistent.
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
  Point _direction;
  const Real _c;

  const unsigned int _lm_var_num;
  const VariableValue & _lambda; ///< Per-node LM values after reinitNodes.

  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var_num;
  std::vector<const VariableValue *> _disp;
};
