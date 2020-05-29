//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HomogenizationConstraint.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseVariableFE.h"
#include "MooseVariableScalar.h"
#include "SystemBase.h"
#include "Function.h"

#include "libmesh/quadrature.h"

registerMooseObject("TensorMechanicsApp", HomogenizationConstraint);

InputParameters
HomogenizationConstraint::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredCoupledVar("lambda", "Lagrange multiplier");
  params.addRequiredParam<unsigned int>("component",
                                        "Which direction this kernel acts in");
  params.addRequiredCoupledVar("displacements", "The displacement components");
  
  MooseEnum homogenizationType("cauchy");
  params.addRequiredParam<MooseEnum>("quantity", homogenizationType, 
                                     "Which quantity to homogenize");

  params.addRequiredParam<unsigned int>("index_i",
                                        "First index of controlled quantity");
  params.addRequiredParam<unsigned int>("index_j",
                                        "Second index of controlled quantity");
  params.addRequiredParam<FunctionName>("target",
                                        "Function giving the target to hit");

  return params;
}

HomogenizationConstraint::HomogenizationConstraint(const InputParameters & parameters)
  : Kernel(parameters),
    _ld(getParam<bool>("use_displaced_mesh")),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _disp_nums(_ndisp),
    _disp_vars(_ndisp),
    _lambda_num(coupledScalar("lambda")), 
    _lambda(coupledScalarValue("lambda")),
    _index_i(getParam<unsigned int>("index_i")),
    _index_j(getParam<unsigned int>("index_j")),
    _target(getFunction("target")),
    _type(getParam<MooseEnum>("quantity").getEnum<HomogenizationType>()),
    _stress(getMaterialPropertyByName<RankTwoTensor>("stress")),
    _material_jacobian(
        getMaterialPropertyByName<RankFourTensor>("Jacobian_mult"))
{
}

void
HomogenizationConstraint::initialSetup() {
  for (unsigned int i = 0; i < _ndisp; i++) {
    _disp_nums[i] = coupled("displacements", i);
    _disp_vars[i] = getVar("displacements", i);
  }
}

void
HomogenizationConstraint::computeResidual()
{
  // The contribution for the base kernel
  DenseVector<Number> & re_base = _assembly.residualBlock(_var.number());

  // The contribution for the scalar kernel (size 1)
  DenseVector<Number> & re_scalar = _assembly.residualBlock(_lambda_num);

  for (_qp = 0; _qp < _qrule->n_points(); _qp++) {
    re_scalar(0) += _JxW[_qp] * _coord[_qp] * computeQpResidualScalar();
    for (_i = 0; _i < _test.size(); _i++) {
      re_base(_i) += _JxW[_qp] * _coord[_qp] * computeQpResidualBase();
    }
  }
}

Real
HomogenizationConstraint::computeQpResidualBase()
{
  return _lambda[0] * computeConstraint();
}

Real HomogenizationConstraint::computeConstraint()
{
  Real value = 0.0;
  for (unsigned int l = 0; l < _ndisp; l++) {
    value -= _material_jacobian[_qp](_index_i,_index_j, _component, l) *
        _grad_test[_i][_qp](l);
  }
  return value;
}

Real
HomogenizationConstraint::computeQpResidualScalar()
{
  return -(_stress[_qp](_index_i, _index_j) - _target.value(_t, _q_point[_qp]))
      / _ndisp;
}

void
HomogenizationConstraint::computeOffDiagJacobianScalar(unsigned int jvar)
{
  DenseMatrix<Number> & ken = _assembly.jacobianBlock(_var.number(), jvar);
  DenseMatrix<Number> & kne = _assembly.jacobianBlock(jvar, _var.number());
  MooseVariableScalar & jv = _sys.getScalarVariable(_tid, jvar);

  for (_i = 0; _i < _test.size(); _i++)
    for (_j = 0; _j < jv.order(); _j++)
      for (_qp = 0; _qp < _qrule->n_points(); _qp++)
      {
        Real value = _JxW[_qp] * _coord[_qp] * computeQpOffDiagJacobian(jvar);
        ken(_i, _j) += value;
        kne(_j, _i) += value;
      }
}

Real
HomogenizationConstraint::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _lambda_num)
    return computeConstraint();
  else
    return 0.0;
}
