//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "KokkosComputeLagrangianStressCustomPK2.h"
#include "FEProblemBase.h"

registerKokkosMaterial("SolidMechanicsApp", KokkosComputeLagrangianStressCustomPK2);

InputParameters
KokkosComputeLagrangianStressCustomPK2::validParams()
{
  auto params = Material::validParams();
  params.addClassDescription(
      "Converts a full PK2 stress and its deformation-gradient derivative to PK1 stress and its "
      "displacement-gradient derivative.");
  params.addRequiredCoupledVar("displacements", "Displacement variables in Cartesian order.");
  params.addRequiredParam<MaterialPropertyName>("custom_pk2_stress", "Name of the full PK2 stress.");
  params.addRequiredParam<MaterialPropertyName>(
      "custom_pk2_jacobian", "Name of the full PK2 derivative with respect to F.");
  params.addParam<MaterialPropertyName>("pk1_stress", "pk1_stress", "Name of the PK1 stress.");
  params.addParam<MaterialPropertyName>("dpk1_d_grad_u",
                                        "dpk1_d_grad_u",
                                        "Name of the PK1 displacement-gradient derivative.");
  return params;
}

KokkosComputeLagrangianStressCustomPK2::KokkosComputeLagrangianStressCustomPK2(
    const InputParameters & parameters)
  : Material(parameters),
    _ndisp(coupledComponents("displacements")),
    _grad_displacements(kokkosCoupledGradients("displacements")),
    _pk2_stress(getKokkosMaterialProperty<Real, 2>("custom_pk2_stress")),
    _dpk2_dF(getKokkosMaterialProperty<Real, 4>("custom_pk2_jacobian")),
    _pk1_stress(declareKokkosProperty<Real, 2>("pk1_stress", {3, 3})),
    _dpk1_d_grad_u(declareKokkosProperty<Real, 4>("dpk1_d_grad_u", {3, 3, 3, 3})),
    _need_jacobian(false),
    _need_jacobian_device(_need_jacobian)
{
  if (_mesh.dimension() != 3 || _ndisp != 3)
    paramError("displacements", "Exactly three displacement variables on a 3D mesh are required.");
}

void
KokkosComputeLagrangianStressCustomPK2::computeProperties()
{
  _need_jacobian = _fe_problem.currentlyComputingJacobian() ||
                   _fe_problem.currentlyComputingResidualAndJacobian();
  Material::computeProperties();
}
