//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "IterativeMultiAppSolve.h"
#include "NonlinearSystem.h"

class SteffensenSolve : public IterativeMultiAppSolve
{
public:
  SteffensenSolve(Executioner * ex);

  static InputParameters validParams();

  /// Allocate storage for secondary transformed objects
  virtual void allocateStorageForSecondaryTransformed() override final
  {
    // Store a copy of the previous solution here
    _problem.getNonlinearSystemBase().addVector("secondary_xn_m1", false, PARALLEL);
    _problem.getNonlinearSystemBase().addVector("secondary_fxn_m1", false, PARALLEL);

    // Allocate storage for the previous postprocessor values
    _secondary_transformed_pps_values.resize(_secondary_transformed_pps.size());
    for (size_t i = 0; i < _secondary_transformed_pps.size(); i++)
      _secondary_transformed_pps_values[i].resize(2);
  }

private:
  /// Save the variable values as a SubApp
  virtual void savePreviousVariableValuesAsSubApp() override final;

  /// Save the postprocessor values as a SubApp
  virtual void savePreviousPostprocessorValuesAsSubApp() override final;

  /// Whether to use the coupling algorithm (relaxed Picard, Secant, ...) instead of Picard
  virtual bool useCouplingAlgorithmUpdate(bool as_main_app) override final;

  /// Save the previous variables and postprocessors as the main application
  virtual void savePreviousValuesAsMainApp() override final;

  /// Compute the new value of the coupling postprocessors based on the coupling algorithm selected
  virtual void transformPostprocessorsAsMainApp() override final;

  /// Compute the new value of the coupling postprocessors based on the coupling algorithm selected as a SubApp
  virtual void transformPostprocessorsAsSubApp() override final;

  /// Compute the new variable values based on the coupling algorithm selected
  virtual void
  transformVariablesAsMainApp(const std::set<dof_id_type> & transformed_dofs) override final;

  /// Compute the new variable values based on the coupling algorithm selected as a SubApp
  virtual void transformVariablesAsSubApp(
      const std::set<dof_id_type> & secondary_transformed_dofs) override final;

  /// Print the convergence history of the coupling, at every coupling iteration
  virtual void printCouplingConvergenceHistory() override final;
};
