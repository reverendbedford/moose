//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"

class HomogenizationConstraint : public Kernel
{
 public:
  static InputParameters validParams();

  HomogenizationConstraint(const InputParameters & parameters);
 
  virtual void initialSetup();

  virtual void computeResidual();
  virtual void computeOffDiagJacobianScalar(unsigned int jvar);

 protected:
  virtual Real computeQpResidual() {return 0;};
  virtual Real computeQpResidualBase();
  virtual Real computeQpResidualScalar();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  virtual Real computeConstraint();

 protected:
  bool _ld;
  unsigned int _component;

  unsigned int _ndisp;

  std::vector<unsigned int> _disp_nums;
  std::vector<MooseVariable *> _disp_vars;

  unsigned int _lambda_num;
  VariableValue & _lambda;

  unsigned int _index_i;
  unsigned int _index_j;

  const Function & _target;

  enum class HomogenizationType {Cauchy} _type;

  const MaterialProperty<RankTwoTensor> &_stress;
  const MaterialProperty<RankFourTensor> &_material_jacobian;
};
