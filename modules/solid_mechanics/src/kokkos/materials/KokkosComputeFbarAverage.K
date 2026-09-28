//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "KokkosComputeFbarAverage.h"

registerKokkosMaterial("SolidMechanicsApp", KokkosComputeFbarAverage);

InputParameters
KokkosComputeFbarAverage::validParams()
{
  auto params = Material::validParams();
  params.addClassDescription(
      "Element-constant Kokkos material publishing the F-bar (total-mode) element averages "
      "F_avg and F_avg^{-T} needed by the stabilized Kokkos total-Lagrangian NEML2 path. "
      "Matches the F_bar_mode = 'total' branch of the host ComputeLagrangianStrain calculator.");
  params.addRequiredCoupledVar("displacements", "Displacement variables in Cartesian order.");
  params.addParam<MaterialPropertyName>(
      "F_avg", "F_avg", "Name of the element-averaged deformation gradient material property.");
  params.addParam<MaterialPropertyName>(
      "F_avg_invT",
      "F_avg_invT",
      "Name of the element-averaged inverse-transpose deformation gradient material property.");

  // Force the property to be element-constant. The framework's `constant_on` parameter is added
  // by KokkosMaterialBase; setting the default here means downstream materials automatically
  // access F_avg via the per-element storage index, and `computeQpProperties` is only invoked
  // once per element rather than once per qp.
  params.set<MooseEnum>("constant_on") = "ELEMENT";
  return params;
}

KokkosComputeFbarAverage::KokkosComputeFbarAverage(const InputParameters & parameters)
  : Material(parameters),
    _ndisp(coupledComponents("displacements")),
    _grad_displacements(kokkosCoupledGradients("displacements")),
    _F_avg(declareKokkosProperty<Real, 2>("F_avg", {3, 3})),
    _F_avg_invT(declareKokkosProperty<Real, 2>("F_avg_invT", {3, 3}))
{
  if (_mesh.dimension() != 3 || _ndisp != 3)
    paramError("displacements",
               "KokkosComputeFbarAverage requires exactly three displacement variables on a 3D "
               "mesh (the total-mode F-bar formula uses a 3D volume ratio).");

  if (parameters.get<MooseEnum>("constant_on") != "ELEMENT")
    paramError("constant_on",
               "KokkosComputeFbarAverage requires constant_on = 'ELEMENT'. The F_avg tensor is a "
               "single per-element quantity; per-qp storage would compute the same value at every "
               "qp and use n_qps times as much memory.");
}
