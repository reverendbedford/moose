//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "NEML2ModelExecutor.h"

/**
 * Checks executor evaluation counts and verifies state preservation across batch remapping.
 */
class TestNEML2ModelExecutor : public NEML2ModelExecutor
{
public:
  static InputParameters validParams();

  TestNEML2ModelExecutor(const InputParameters & params);

#ifdef NEML2_ENABLED
  void execute() override;

protected:
  bool solve(bool compute_derivative) override;
  void remapState() override;

private:
  unsigned int _solve_calls;
#endif
};
