//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "NodalKernel.h"

class LevelSetContactor;

/**
 * Rigid-body frictionless contact - node-wise NCP enforced against an
 * analytic level-set contactor. Applied at each Lagrange-multiplier DoF
 * living on the deformable contact sideset's lower-d block.
 *
 *   R_i = min( λ_i,  c · g_LS(x_i + u_i) )
 *
 * Complementarity is completed by PETSc `SNESVINEWTONSSLS` + `ConstantBounds`
 * enforcing λ ≥ 0. The nodal min-NCP puts a nonzero {0, 1} on the MOOSE
 * assembled Jacobian's LM diagonal (same structural benefit that mortar's
 * `enforceConstraintOnDof` provides). No mortar segment mesh, no dual-basis
 * integration, no AD.
 *
 * Jacobian:
 *   λ-branch  (λ ≤ c · g):  ∂R/∂λ = 1,  ∂R/∂disp_k = 0
 *   g-branch  (c · g < λ):  ∂R/∂λ = 0,  ∂R/∂disp_k = c · n_k(x + u)
 *
 * Companion class: RigidBodyNormalMechanicalContact applies the
 * -λ · n · φ_test traction on the coupled displacement equations.
 */
class RigidBodyNodalNCPKernel : public NodalKernel
{
public:
  static InputParameters validParams();
  RigidBodyNodalNCPKernel(const InputParameters &);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  Point deformedNode() const;
  Real physicalGap() const;

  const LevelSetContactor & _contactor;
  const Real _c;
  const unsigned int _ndisp;
  std::vector<const VariableValue *> _disp;
  std::vector<unsigned int> _disp_num;
};
