//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DummyLagrange.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseVariableScalar.h"
#include "Function.h"

registerMooseObject("TensorMechanicsApp", DummyLagrange);

InputParameters
DummyLagrange::validParams()
{
  InputParameters params = ScalarKernel::validParams();

  return params;
}

DummyLagrange::DummyLagrange(const InputParameters & parameters)
  : ScalarKernel(parameters)
{
}

void
DummyLagrange::reinit()
{
}

void
DummyLagrange::computeResidual()
{
  /*
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  for (_i = 0; _i < re.size(); _i++)
    re(_i) += computeQpResidual();
  */
  return; // yes, do nothing, I take care of this in the other class...
}

void
DummyLagrange::computeJacobian()
{
  // Amusingly you need an explicit zero or PETSC whines
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
  for (_i = 0; _i < ke.m(); _i++)
    ke(_i, _i) += 0.0;
}

// Off-diagonal Jacobian dealt with in the other class...
