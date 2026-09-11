//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "KokkosStressDivergence.h"

registerKokkosResidualObject("SolidMechanicsApp", KokkosStressDivergence);

InputParameters
KokkosStressDivergence::validParams()
{
  InputParameters params = KernelGrad::validParams();
  params.addClassDescription(
      "Computes Cartesian stress divergence from Kokkos stress and tangent material properties.");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "component", "component < 3", "Cartesian displacement component for this kernel.");
  params.addRequiredCoupledVar("displacements", "Displacement variables in Cartesian order.");
  params.addParam<MaterialPropertyName>("stress", "stress", "Name of the stress material property.");
  params.addParam<MaterialPropertyName>(
      "tangent", "Jacobian_mult", "Name of the constitutive tangent material property.");
  return params;
}

KokkosStressDivergence::KokkosStressDivergence(const InputParameters & parameters)
  : KernelGrad(parameters),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _stress(getKokkosMaterialProperty<Real, 2>("stress")),
    _tangent(getKokkosMaterialProperty<Real, 4>("tangent")),
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
