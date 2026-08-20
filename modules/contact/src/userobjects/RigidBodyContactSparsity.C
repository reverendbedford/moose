//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RigidBodyContactSparsity.h"

#include "MooseMesh.h"
#include "MooseVariableFE.h"
#include "NonlinearSystemBase.h"
#include "SystemBase.h"

#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"

registerMooseObject("ContactApp", RigidBodyContactSparsity);

InputParameters
RigidBodyContactSparsity::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription(
      "Preallocates the cross-node (LM, disp) Jacobian coupling that "
      "RigidBodyNormalMechanicalContact writes on the contact sideset's "
      "lower-d elements, but MOOSE's default sparsity computation misses.");
  params.addRequiredParam<VariableName>("lm_variable",
                                        "The Lagrange multiplier field variable on the contact "
                                        "lower-d block (same as the NCP kernel's `variable`).");
  params.addRequiredParam<std::vector<VariableName>>(
      "displacements", "Displacement variables in order (x, y[, z]).");
  // Attach in the constructor, which runs before es().init() where the
  // sparsity pattern is computed.  execute_on defaults to NONE (nothing to
  // do at runtime).
  params.set<ExecFlagEnum>("execute_on") = EXEC_NONE;
  return params;
}

RigidBodyContactSparsity::RigidBodyContactSparsity(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _lm_var_num(_fe_problem.getVariable(0, getParam<VariableName>("lm_variable")).number()),
    _ndisp(getParam<std::vector<VariableName>>("displacements").size()),
    _disp_var_num(_ndisp)
{
  const auto & disp_names = getParam<std::vector<VariableName>>("displacements");
  for (const auto k : make_range(_ndisp))
    _disp_var_num[k] = _fe_problem.getVariable(0, disp_names[k]).number();

  // Attach ourselves as an AugmentSparsityPattern on the nonlinear system's
  // DofMap.  The (function, object) slots are separate on libMesh's DofMap,
  // so we coexist with MOOSE's existing extraSparsity function callback.
  auto & nl = _fe_problem.getNonlinearSystemBase(/*sys_num=*/0);
  nl.dofMap().attach_extra_sparsity_object(*this);
}

void
RigidBodyContactSparsity::augment_sparsity_pattern(
    libMesh::SparsityPattern::Graph & sparsity,
    std::vector<libMesh::dof_id_type> & n_nz,
    std::vector<libMesh::dof_id_type> & n_oz)
{
  auto & nl = _fe_problem.getNonlinearSystemBase(/*sys_num=*/0);
  const auto & dof_map = nl.dofMap();
  auto & mesh = _fe_problem.mesh();
  const auto proc = processor_id();
  const auto first_dof_on_proc = dof_map.first_dof(proc);
  const auto end_dof_on_proc = dof_map.end_dof(proc);
  const auto n_dofs_on_proc = dof_map.n_local_dofs();
  const auto n_dofs_not_on_proc = dof_map.n_dofs() - dof_map.n_local_dofs();

  const auto & lm_var = nl.getVariable(0, _lm_var_num);
  const std::set<SubdomainID> & lm_blocks = lm_var.blockIDs();

  // For each lower-d block element, gather ALL DoFs from LM + every disp
  // component (from BOTH the lower-d element itself AND its higher-d
  // parent).  Mark every pair as coupled.  This closes the preallocation
  // gap where MOOSE's LowerDIntegratedBC assembles cross-node
  // (LM,disp), (disp,LM), and (disp_i_on_primary, disp_j_on_lower)
  // blocks that libMesh's default element-based sparsity misses.
  std::vector<libMesh::dof_id_type> all_dofs;
  std::vector<libMesh::dof_id_type> di;

  auto gather = [&](const libMesh::Elem * elem)
  {
    for (const auto k : make_range(_ndisp))
    {
      di.clear();
      dof_map.dof_indices(elem, di, _disp_var_num[k]);
      all_dofs.insert(all_dofs.end(), di.begin(), di.end());
    }
    di.clear();
    dof_map.dof_indices(elem, di, _lm_var_num);
    all_dofs.insert(all_dofs.end(), di.begin(), di.end());
  };

  for (const auto * elem : mesh.getMesh().active_local_element_ptr_range())
  {
    if (!lm_blocks.count(elem->subdomain_id()))
      continue;

    all_dofs.clear();
    gather(elem);
    // Also gather from the interior (higher-d) parent so we capture
    // (higher-d test) × (lower-d phi) blocks the LowerDIntegratedBC
    // creates via computeLowerDOffDiagJacobian(PrimaryLower, ...).
    if (const auto * parent = elem->interior_parent())
      gather(parent);

    std::sort(all_dofs.begin(), all_dofs.end());
    all_dofs.erase(std::unique(all_dofs.begin(), all_dofs.end()), all_dofs.end());
    if (all_dofs.empty())
      continue;

    for (const auto r : all_dofs)
    {
      if (r < first_dof_on_proc || r >= end_dof_on_proc)
        continue;
      const auto local = r - first_dof_on_proc;
      auto & row = sparsity[local];
      const auto old_size = row.size();
      for (const auto c : all_dofs)
      {
        if (!std::binary_search(row.begin(), row.begin() + old_size, c))
        {
          row.push_back(c);
          if (c < first_dof_on_proc || c >= end_dof_on_proc)
          {
            if (n_oz[local] < n_dofs_not_on_proc)
              n_oz[local]++;
          }
          else
          {
            if (n_nz[local] < n_dofs_on_proc)
              n_nz[local]++;
          }
        }
      }
      std::sort(row.begin() + old_size, row.end());
      std::inplace_merge(row.begin(), row.begin() + old_size, row.end());
    }
  }
}
