//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "NEML2OutputInterface.h"

#include "NEML2ModelExecutor.h"
#include "NEML2Utils.h"
#include "UserObjectInterface.h"

InputParameters
NEML2OutputInterface::validParams()
{
  auto params = emptyInputParameters();
  params.addRequiredParam<UserObjectName>("neml2_executor",
                                           "User object managing the execution of the NEML2 model.");
  params.addRequiredParam<std::string>("from_neml2", "NEML2 output variable to read from");
  params.addParam<std::string>(
      "neml2_input_derivative",
      "If supplied return the derivative of the NEML2 output variable with respect to this input");
  params.addParam<std::string>(
      "neml2_parameter_derivative",
      "If supplied return the derivative of the NEML2 output variable with respect to this NEML2 "
      "model parameter");
  return params;
}

NEML2OutputInterface::NEML2OutputInterface(const InputParameters & params,
                                           const UserObjectInterface & interface)
#ifdef NEML2_ENABLED
  : _execute_neml2_model(interface.getUserObject<NEML2ModelExecutor>("neml2_executor")),
    _value(getNEML2Output(params, _execute_neml2_model))
#endif
{
#ifndef NEML2_ENABLED
  (void)params;
  (void)interface;
#endif
  NEML2Utils::assertNEML2Enabled();
}

#ifdef NEML2_ENABLED
const at::Tensor &
NEML2OutputInterface::getNEML2Output(const InputParameters & params,
                                     const NEML2ModelExecutor & executor)
{
  if (params.isParamValid("neml2_input_derivative") &&
      params.isParamValid("neml2_parameter_derivative"))
    params.paramError("neml2_parameter_derivative",
                      "An input derivative and parameter derivative cannot both be requested.");

  const auto & output = params.get<std::string>("from_neml2");
  if (params.isParamValid("neml2_input_derivative"))
    return executor.getOutputDerivative(output,
                                        params.get<std::string>("neml2_input_derivative"));
  if (params.isParamValid("neml2_parameter_derivative"))
    return executor.getOutputParameterDerivative(
        output, params.get<std::string>("neml2_parameter_derivative"));
  return executor.getOutput(output);
}

bool
NEML2OutputInterface::outputReady() const
{
  return _execute_neml2_model.outputReady();
}

const c10::Device &
NEML2OutputInterface::neml2Device() const
{
  return _execute_neml2_model.device();
}

void
NEML2OutputInterface::requestNEML2OutputDevice(const c10::Device & device) const
{
  _execute_neml2_model.setOutputDevice(device);
}

std::size_t
NEML2OutputInterface::neml2BatchIndex(const dof_id_type elem_id) const
{
  return _execute_neml2_model.getBatchIndex(elem_id);
}
#endif
