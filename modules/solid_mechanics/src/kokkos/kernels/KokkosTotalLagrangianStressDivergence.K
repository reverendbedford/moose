//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "KokkosTotalLagrangianStressDivergence.h"

registerKokkosResidualObject("SolidMechanicsApp", KokkosTotalLagrangianStressDivergence);

InputParameters
KokkosTotalLagrangianStressDivergence::validParams()
{
  auto params = KernelGrad::validParams();
  params.addClassDescription(
      "Computes total-Lagrangian stress divergence from full Kokkos PK1 stress and tangent "
      "material properties. Optionally adds the total-mode F-bar non-local Jacobian coupling "
      "when `stabilize_strain = true`.");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "component", "component < 3", "Cartesian displacement component for this kernel.");
  params.addRequiredCoupledVar("displacements", "Displacement variables in Cartesian order.");
  params.addParam<MaterialPropertyName>(
      "pk1_stress", "pk1_stress", "Name of the first Piola-Kirchhoff stress.");
  params.addParam<MaterialPropertyName>(
      "dpk1_d_grad_u", "dpk1_d_grad_u", "Name of the PK1 displacement-gradient derivative.");
  params.addParam<bool>(
      "stabilize_strain",
      false,
      "If true, add the total-mode F-bar cross-quadrature Jacobian coupling. Requires "
      "KokkosComputeFbarAverage (publishing F_avg_invT) and KokkosComputeLagrangianStressCustomPK2 "
      "with matching `stabilize_strain = true` (publishing A(qp) and gamma/3(qp)). Off by "
      "default preserves the original per-qp assembly path bit-for-bit.");
  params.addParam<MaterialPropertyName>(
      "F_avg_invT",
      "F_avg_invT",
      "Name of the inverse-transpose element-averaged deformation gradient material property.");
  params.addParam<MaterialPropertyName>(
      "fbar_A", "fbar_A", "Name of the F-bar A(qp) rank-two tensor.");
  params.addParam<MaterialPropertyName>(
      "fbar_gamma_over_3", "fbar_gamma_over_3", "Name of the F-bar gamma/3 scalar.");
  return params;
}

KokkosTotalLagrangianStressDivergence::KokkosTotalLagrangianStressDivergence(
    const InputParameters & parameters)
  : KernelGrad(parameters),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _pk1_stress(getKokkosMaterialProperty<Real, 2>("pk1_stress")),
    _dpk1_d_grad_u(getKokkosMaterialProperty<Real, 4>("dpk1_d_grad_u")),
    _displacement_var_ids(_ndisp),
    _stabilize_strain(getParam<bool>("stabilize_strain")),
    _stabilize_strain_device(_stabilize_strain),
    _F_avg_invT(_stabilize_strain ? getKokkosMaterialProperty<Real, 2>("F_avg_invT")
                                  : Moose::Kokkos::MaterialProperty<Real, 2>()),
    _A(_stabilize_strain ? getKokkosMaterialProperty<Real, 2>("fbar_A")
                         : Moose::Kokkos::MaterialProperty<Real, 2>()),
    _gamma_over_3(_stabilize_strain ? getKokkosMaterialProperty<Real, 0>("fbar_gamma_over_3")
                                    : Moose::Kokkos::MaterialProperty<Real, 0>())
{
  if (_mesh.dimension() != 3 || _ndisp != 3)
    paramError("displacements", "Exactly three displacement variables on a 3D mesh are required.");
  if (_component >= _ndisp)
    paramError("component", "The component must identify one of the supplied displacement variables.");
  if (_var.number() != coupled("displacements", _component))
    paramError("variable", "The kernel variable must match the displacement selected by component.");

  for (const auto i : make_range(_ndisp))
    _displacement_var_ids[i] = coupled("displacements", i);
  _displacement_var_ids.copyToDevice();
}
