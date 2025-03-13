//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MeshGenerator.h"

/*
 * Mesh generator to merge certain duplicate nodes defined by coordinates
 */
class MergeCertainNodes : public MeshGenerator
{
public:
  static InputParameters validParams();

  MergeCertainNodes(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  /// mesh to modify
  std::unique_ptr<MeshBase> & _input;

  /// the coordinates of the nodes to merge
  const std::vector<std::vector<Real>> _coords;
  /// the tolerance for merging nodes
  const Real _tol;
};