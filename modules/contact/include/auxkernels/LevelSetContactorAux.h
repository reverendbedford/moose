//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"
#include "MooseEnum.h"

class LevelSetContactor;

/**
 * Auxiliary kernel that samples a LevelSetContactor UserObject at each
 * quadrature or nodal point. Selected quantity is either the signed distance
 * or one Cartesian component of the outward normal.
 *
 * Primary use: regression testing and visualization of an analytic contactor.
 */
class LevelSetContactorAux : public AuxKernel
{
public:
  static InputParameters validParams();

  LevelSetContactorAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// The contactor being sampled.
  const LevelSetContactor & _contactor;

  /// Which scalar quantity to output.
  const MooseEnum _quantity;
};
