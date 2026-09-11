//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "KokkosComputeIsotropicElasticity.h"

registerKokkosMaterial("SolidMechanicsApp", KokkosComputeIsotropicElasticity);

InputParameters
KokkosComputeIsotropicElasticity::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Computes isotropic small-strain stress and tangent Kokkos material properties.");
  params.addRequiredCoupledVar("displacements", "Displacement variables in Cartesian order.");
  params.addRequiredRangeCheckedParam<Real>(
      "youngs_modulus", "youngs_modulus > 0", "Young's modulus.");
  params.addRequiredRangeCheckedParam<Real>(
      "poissons_ratio",
      "poissons_ratio > -1 & poissons_ratio < 0.5",
      "Poisson's ratio.");
  params.addParam<MaterialPropertyName>("stress", "stress", "Name of the stress material property.");
  params.addParam<MaterialPropertyName>(
      "tangent", "Jacobian_mult", "Name of the constitutive tangent material property.");
  return params;
}

KokkosComputeIsotropicElasticity::KokkosComputeIsotropicElasticity(
    const InputParameters & parameters)
  : Material(parameters),
    _ndisp(coupledComponents("displacements")),
    _lambda(getParam<Real>("youngs_modulus") * getParam<Real>("poissons_ratio") /
            ((1 + getParam<Real>("poissons_ratio")) *
             (1 - 2 * getParam<Real>("poissons_ratio")))),
    _mu(getParam<Real>("youngs_modulus") / (2 * (1 + getParam<Real>("poissons_ratio")))),
    _grad_displacements(kokkosCoupledGradients("displacements")),
    _stress(declareKokkosProperty<Real, 2>("stress", {_ndisp, _ndisp})),
    _tangent(declareKokkosProperty<Real, 4>("tangent", {_ndisp, _ndisp, _ndisp, _ndisp}))
{
  if (_ndisp != _mesh.dimension())
    paramError("displacements", "The number of displacement variables must match the mesh dimension.");
}
