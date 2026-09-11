//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "TestNEML2ModelExecutor.h"

registerMooseObject("SolidMechanicsTestApp", TestNEML2ModelExecutor);

InputParameters
TestNEML2ModelExecutor::validParams()
{
  return NEML2ModelExecutor::validParams();
}

TestNEML2ModelExecutor::TestNEML2ModelExecutor(const InputParameters & params)
  : NEML2ModelExecutor(params)
#ifdef NEML2_ENABLED
    ,
    _solve_calls(0)
#endif
{
}

#ifdef NEML2_ENABLED
void
TestNEML2ModelExecutor::execute()
{
  _solve_calls = 0;
  NEML2ModelExecutor::execute();

  if (_solve_calls > 1)
    mooseError("NEML2 model was evaluated ", _solve_calls, " times in one executor dispatch");
  if (_solve_calls == 1)
    mooseInfo("NEML2 model was evaluated once in this executor dispatch");
}

bool
TestNEML2ModelExecutor::solve(const bool compute_derivative)
{
  ++_solve_calls;
  return NEML2ModelExecutor::solve(compute_derivative);
}

void
TestNEML2ModelExecutor::remapState()
{
  if (!_state_committed || _state_vars.empty())
  {
    NEML2ModelExecutor::remapState();
    return;
  }

  const auto old_states = _state_vars;
  const auto old_indices = _state_batch_indices;
  const auto new_indices = _batch_index_generator.getBatchIndexMap();
  NEML2ModelExecutor::remapState();

  std::map<std::size_t, dof_id_type> old_offsets;
  std::map<std::size_t, dof_id_type> new_offsets;
  for (const auto & [elem_id, offset] : old_indices)
    old_offsets[offset] = elem_id;
  for (const auto & [elem_id, offset] : new_indices)
    new_offsets[offset] = elem_id;

  const auto batchSize =
      [](const auto & offsets, const std::size_t offset, const std::size_t total_batch_size)
  {
    const auto next = offsets.upper_bound(offset);
    return next == offsets.end() ? total_batch_size - offset : next->first - offset;
  };

  const auto & input_names = model().input_names();
  const auto & input_shapes = model().input_base_shapes();
  unsigned int retained = 0;
  unsigned int initialized = 0;
  for (const auto & [name, old_state] : old_states)
  {
    const auto & state = libmesh_map_find(_state_vars, name);
    const auto name_it = std::find(input_names.begin(), input_names.end(), name);
    mooseAssert(name_it != input_names.end(), "State variable is not a model input");
    const auto base_ndim =
        static_cast<int64_t>(input_shapes[std::distance(input_names.begin(), name_it)].size());
    if (old_state.dim() == base_ndim)
    {
      if (!torch::equal(state, old_state))
        mooseError("NEML2 remapping changed unbatched state variable '", name, "'");
      continue;
    }

    for (const auto & [elem_id, new_offset] : new_indices)
    {
      const auto old_it = old_indices.find(elem_id);
      const auto new_qp = batchSize(new_offsets, new_offset, state.size(0));
      const auto new_slice = state.narrow(0, new_offset, new_qp);
      if (old_it == old_indices.end())
      {
        if (torch::count_nonzero(new_slice).item<int64_t>() != 0)
          mooseError("NEML2 remapping did not zero-initialize state variable '",
                     name,
                     "' for new element ",
                     elem_id);
        ++initialized;
        continue;
      }

      const auto old_qp = batchSize(old_offsets, old_it->second, old_state.size(0));
      if (!torch::equal(new_slice, old_state.narrow(0, old_it->second, old_qp)))
        mooseError("NEML2 remapping changed state variable '",
                   name,
                   "' for retained element ",
                   elem_id);
      ++retained;
    }
  }

  mooseInfo("Verified NEML2 state remap with ",
            retained,
            " retained and ",
            initialized,
            " initialized state slices");
}
#endif
