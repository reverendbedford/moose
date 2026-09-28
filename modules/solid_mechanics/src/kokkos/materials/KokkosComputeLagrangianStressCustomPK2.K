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
      "displacement-gradient derivative. Optionally applies the total-mode F-bar Jacobian chain "
      "to match the host `ComputeLagrangianStressCustomPK2` when the strain calculator has "
      "`stabilize_strain = true`.");
  params.addRequiredCoupledVar("displacements", "Displacement variables in Cartesian order.");
  params.addRequiredParam<MaterialPropertyName>("custom_pk2_stress", "Name of the full PK2 stress.");
  params.addRequiredParam<MaterialPropertyName>(
      "custom_pk2_jacobian", "Name of the full PK2 derivative with respect to F.");
  params.addParam<MaterialPropertyName>("pk1_stress", "pk1_stress", "Name of the PK1 stress.");
  params.addParam<MaterialPropertyName>("dpk1_d_grad_u",
                                        "dpk1_d_grad_u",
                                        "Name of the PK1 displacement-gradient derivative.");
  params.addParam<bool>(
      "stabilize_strain",
      false,
      "If true, chain the total-mode F-bar dF_stab/dF_ust into the local PK1 Jacobian and "
      "publish A(qp) and gamma/3 for the non-local Jacobian assembled by the "
      "KokkosTotalLagrangianStressDivergence kernel. Requires KokkosComputeFbarAverage to publish "
      "F_avg (and F_avg_invT for the kernel). Pair with matching `stabilize_strain = true` on the "
      "kernel and on `TorchDeformationGradient` (which sends F_stab to NEML2). Off by default "
      "preserves the original unstabilized behavior bit-for-bit.");
  params.addParam<MaterialPropertyName>(
      "F_avg", "F_avg", "Name of the element-averaged deformation gradient (F-bar input).");
  params.addParam<MaterialPropertyName>(
      "F_avg_invT",
      "F_avg_invT",
      "Name of the inverse-transpose element-averaged deformation gradient (declared for pairing "
      "with the kernel, not read here).");
  params.addParam<MaterialPropertyName>(
      "fbar_A",
      "fbar_A",
      "Name of the F-bar A(qp) rank-two tensor published for the kernel's non-local Jacobian.");
  params.addParam<MaterialPropertyName>(
      "fbar_gamma_over_3",
      "fbar_gamma_over_3",
      "Name of the F-bar gamma/3 scalar published for the kernel's non-local Jacobian.");
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
    _need_jacobian_device(_need_jacobian),
    _stabilize_strain(getParam<bool>("stabilize_strain")),
    _stabilize_strain_device(_stabilize_strain),
    // Element-constant F-bar inputs (F_avg, F_avg_invT) and per-qp outputs (A, gamma/3) are
    // declared / fetched only when stabilization is on; otherwise leave them defaulted so no
    // upstream KokkosComputeFbarAverage is required for the unstabilized path.
    _F_avg(_stabilize_strain ? getKokkosMaterialProperty<Real, 2>("F_avg")
                             : Moose::Kokkos::MaterialProperty<Real, 2>()),
    _F_avg_invT(_stabilize_strain ? getKokkosMaterialProperty<Real, 2>("F_avg_invT")
                                  : Moose::Kokkos::MaterialProperty<Real, 2>()),
    _A(_stabilize_strain ? declareKokkosProperty<Real, 2>("fbar_A", {3, 3})
                         : Moose::Kokkos::MaterialProperty<Real, 2>()),
    _gamma_over_3(_stabilize_strain ? declareKokkosProperty<Real, 0>("fbar_gamma_over_3")
                                    : Moose::Kokkos::MaterialProperty<Real, 0>())
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
