//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MergeCertainNodes.h"
#include "MooseMesh.h"
#include "Conversion.h"
#include "MooseMeshUtils.h"
#include "CastUniquePointer.h"

#include "libmesh/elem.h"

registerMooseObject("MooseApp", MergeCertainNodes);

InputParameters
MergeCertainNodes::validParams()
{
  InputParameters params = MeshGenerator::validParams();

  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  params.addParam<std::vector<std::vector<Real>>>(
      "coord",
      {},
      "The nodes with coordinates you want to be in the "
      "nodeset. Separate multple coords with ';'.");
  params.addParam<Real>(
      "tolerance", TOLERANCE, "The tolerance in which two nodes are considered identical");
  params.addClassDescription(
      "Creates a new node set and a new boundary made with the nodes the user provides.");

  return params;
}

MergeCertainNodes::MergeCertainNodes(const InputParameters & parameters)
  : MeshGenerator(parameters),
    _input(getMesh("input")),
    _coords(getParam<std::vector<std::vector<Real>>>("coord")),
    _tol(getParam<Real>("tolerance"))
{
}

std::unique_ptr<MeshBase>
MergeCertainNodes::generate()
{
  std::unique_ptr<MeshBase> mesh = std::move(_input);

  auto pl = mesh->sub_point_locator();
  pl->set_close_to_point_tol(_tol);
  pl->enable_out_of_mesh_mode();

  // Setup
  unsigned int num_fixed_nodes = 0;
  std::unordered_set<dof_id_type> nodes_removed;

  auto dim = mesh->mesh_dimension();

  // loop on points
  for (auto & pt : _coords)
  {
    // Find candidate elements
    Point p;
    if (pt.size() < dim)
      paramError("coord",
                 "Coordinate ",
                 Moose::stringify(pt),
                 " does not have enough components for a ",
                 dim,
                 "D mesh.");

    if (pt.size() > 3)
      paramError("coord",
                 "Coordinate ",
                 Moose::stringify(pt),
                 " has too many components. Did you maybe forget to separate multiple coordinates "
                 "with a ';'?");

    for (unsigned int j = 0; j < pt.size(); ++j)
      p(j) = pt[j];

    // will contain all elements that are close to the point
    std::set<const Elem *> elements;
    (*pl)(p, elements);
    if (elements.size() == 0)
      continue;

    // choose a reference node, this one will be kept while the others will be removed
    auto ref_node = pl->locate_node(p);
    if (!ref_node)
      continue;

    // loop on elements
    for (auto & elem : elements)
    {
      // loop on nodes in this element
      for (auto & elem_node : elem->node_ref_range())
      {
        // this node has already been removed
        if (nodes_removed.count(elem_node.id()))
          continue;

        // this is the reference node, we don't want to merge it
        if (elem_node.id() == ref_node->id())
          continue;

        // check if the node is close to the point
        if (elem_node.absolute_fuzzy_equals(p, _tol))
        {
          // merge the nodes
          const_cast<Elem *>(elem)->set_node(elem->get_node_index(&elem_node)) =
              const_cast<Node *>(ref_node);
          nodes_removed.insert(elem_node.id());

          num_fixed_nodes++;
          if (num_fixed_nodes < 10)
            _console << "Merging nodes " << *ref_node << " and " << elem_node << std::endl;
          else if (num_fixed_nodes == 10)
            _console << "Node merging will now proceed silently." << std::endl;
        }
      }
    }
  }

  if (mesh->allow_renumbering())
    mesh->renumber_nodes_and_elements();
  else
  {
    mesh->remove_orphaned_nodes();
    mesh->update_parallel_id_counts();
  }

  mesh->set_isnt_prepared();
  return dynamic_pointer_cast<MeshBase>(mesh);
}