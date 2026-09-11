//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "KokkosIsotropicElasticity.h"

registerKokkosResidualObject("SolidMechanicsApp", KokkosIsotropicElasticity);

InputParameters
KokkosIsotropicElasticity::validParams()
{
  InputParameters params = KernelGrad::validParams();
  params.addClassDescription(
      "Computes the stress divergence for Cartesian isotropic small-strain elasticity.");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "component", "component < 3", "Cartesian displacement component for this kernel.");
  params.addRequiredCoupledVar("displacements", "Displacement variables in Cartesian order.");
  params.addRequiredRangeCheckedParam<Real>(
      "youngs_modulus", "youngs_modulus > 0", "Young's modulus.");
  params.addRequiredRangeCheckedParam<Real>(
      "poissons_ratio",
      "poissons_ratio > -1 & poissons_ratio < 0.5",
      "Poisson's ratio.");
  return params;
}

KokkosIsotropicElasticity::KokkosIsotropicElasticity(const InputParameters & parameters)
  : KernelGrad(parameters),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _lambda(getParam<Real>("youngs_modulus") * getParam<Real>("poissons_ratio") /
            ((1 + getParam<Real>("poissons_ratio")) *
             (1 - 2 * getParam<Real>("poissons_ratio")))),
    _mu(getParam<Real>("youngs_modulus") / (2 * (1 + getParam<Real>("poissons_ratio")))),
    _grad_displacements(kokkosCoupledGradients("displacements")),
    _displacement_var_ids(_ndisp)
{
  if (_ndisp != _mesh.dimension())
    paramError("displacements", "The number of displacement variables must match the mesh dimension.");
  if (_component >= _ndisp)
    paramError("component", "The component must identify one of the supplied displacement variables.");
  if (_var.number() != coupled("displacements", _component))
    paramError("variable", "The kernel variable must match the displacement selected by component.");

  for (const auto i : make_range(_ndisp))
    _displacement_var_ids[i] = coupled("displacements", i);
  _displacement_var_ids.copyToDevice();
}
