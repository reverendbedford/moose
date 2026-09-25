//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "KokkosFullTensorTestMaterial.h"

registerKokkosMaterial("SolidMechanicsTestApp", KokkosFullTensorTestMaterial);

InputParameters
KokkosFullTensorTestMaterial::validParams()
{
  auto params = Material::validParams();
  params.addRequiredParam<MaterialPropertyName>("rank_two", "Full rank-two property to check.");
  params.addRequiredParam<MaterialPropertyName>("rank_four", "Full rank-four property to check.");
  params.addRequiredParam<MaterialPropertyName>("pk1", "Full PK1 property to check.");
  params.addRequiredParam<std::vector<Real>>("A", "Left multiplication tensor in row-major order.");
  params.addRequiredParam<std::vector<Real>>("B", "Right multiplication tensor in row-major order.");
  params.addRequiredParam<std::vector<Real>>("F", "Input tensor in row-major order.");
  params.addClassDescription("Checks all components of full rank-two and rank-four Kokkos material "
                             "properties against an analytical tensor product.");
  return params;
}

KokkosFullTensorTestMaterial::KokkosFullTensorTestMaterial(const InputParameters & parameters)
  : Material(parameters),
    _rank_two(getKokkosMaterialProperty<Real, 2>("rank_two")),
    _rank_four(getKokkosMaterialProperty<Real, 4>("rank_four")),
    _pk1(getKokkosMaterialProperty<Real, 2>("pk1")),
    _rank_two_error(declareKokkosProperty<Real>("rank_two_error")),
    _rank_four_error(declareKokkosProperty<Real>("rank_four_error")),
    _pk1_error(declareKokkosProperty<Real>("pk1_error")),
    _A(9),
    _B(9),
    _F(9)
{
  for (const auto * name : {"A", "B", "F"})
    if (getParam<std::vector<Real>>(name).size() != 9)
      paramError(name, "Exactly nine row-major tensor components are required.");

  for (const auto i : make_range(9))
  {
    _A[i] = getParam<std::vector<Real>>("A")[i];
    _B[i] = getParam<std::vector<Real>>("B")[i];
    _F[i] = getParam<std::vector<Real>>("F")[i];
  }
  _A.copyToDevice();
  _B.copyToDevice();
  _F.copyToDevice();
}
