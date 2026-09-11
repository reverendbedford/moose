//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "InputParameters.h"
#include "MooseTypes.h"

#ifdef NEML2_ENABLED
#include <c10/core/Device.h>
#endif

class NEML2ModelExecutor;
class UserObjectInterface;

#ifdef NEML2_ENABLED
namespace at
{
class Tensor;
}
#endif

/**
 * Common access to an output, input derivative, or parameter derivative of a NEML2 model.
 */
class NEML2OutputInterface
{
public:
  static InputParameters validParams();

  NEML2OutputInterface(const InputParameters & params, const UserObjectInterface & interface);

#ifdef NEML2_ENABLED
  /// Whether the executor has produced output for the current evaluation.
  bool outputReady() const;

  /// The selected output or derivative tensor.
  const at::Tensor & neml2Output() const { return _value; }

  /// Starting batch index for a libMesh element ID.
  std::size_t neml2BatchIndex(dof_id_type elem_id) const;

protected:
  /// NEML2 model compute device used by this output consumer.
  const c10::Device & neml2Device() const;

  /// Require this output to be produced on the consumer's device.
  void requestNEML2OutputDevice(const c10::Device & device) const;

  /// Select the requested output after validating the derivative parameters.
  static const at::Tensor & getNEML2Output(const InputParameters & params,
                                           const NEML2ModelExecutor & executor);

  /// User object managing execution of the NEML2 model.
  const NEML2ModelExecutor & _execute_neml2_model;

  /// Selected output or derivative tensor.
  const at::Tensor & _value;
#endif
};
