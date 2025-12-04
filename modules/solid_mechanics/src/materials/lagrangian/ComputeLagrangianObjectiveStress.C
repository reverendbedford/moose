//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeLagrangianObjectiveStress.h"

#include "FactorizedRankTwoTensor.h"

InputParameters
ComputeLagrangianObjectiveStress::validParams()
{
  InputParameters params = ComputeLagrangianStressCauchy::validParams();

  params.addClassDescription("Stress update based on the small (engineering) stress");

  MooseEnum objectiveRate("truesdell jaumann green_naghdi rashid", "truesdell");
  params.addParam<MooseEnum>(
      "objective_rate", objectiveRate, "Which type of objective integration to use");

  return params;
}

ComputeLagrangianObjectiveStress::ComputeLagrangianObjectiveStress(
    const InputParameters & parameters)
  : ComputeLagrangianStressCauchy(parameters),
    _small_stress(declareProperty<RankTwoTensor>(_base_name + "small_stress")),
    _small_stress_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "small_stress")),
    _small_jacobian(declareProperty<RankFourTensor>(_base_name + "small_jacobian")),
    _cauchy_stress_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "cauchy_stress")),
    _mechanical_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _strain_increment(getMaterialPropertyByName<RankTwoTensor>(_base_name + "strain_increment")),
    _spatial_velocity_gradient_increment(
        getMaterialPropertyByName<RankTwoTensor>(_base_name + "spatial_velocity_increment")),
    _vorticity_increment(
        getMaterialPropertyByName<RankTwoTensor>(_base_name + "vorticity_increment")),
    _def_grad(getMaterialPropertyByName<RankTwoTensor>(_base_name + "deformation_gradient")),
    _def_grad_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "deformation_gradient")),
    _rate(getParam<MooseEnum>("objective_rate").getEnum<ObjectiveRate>()),
    _polar_decomp(_rate == ObjectiveRate::GreenNaghdi || _rate == ObjectiveRate::Rashid),
    _rotation(_polar_decomp ? &declareProperty<RankTwoTensor>(_base_name + "rotation") : nullptr),
    _rotation_old(_polar_decomp ? &getMaterialPropertyOld<RankTwoTensor>(_base_name + "rotation")
                                : nullptr),
    _d_rotation_d_def_grad(
        _polar_decomp ? &declareProperty<RankFourTensor>(derivativePropertyName(
                            _base_name + "rotation", {_base_name + "deformation_gradient"}))
                      : nullptr),
    _stretch(_polar_decomp ? &declareProperty<RankTwoTensor>(_base_name + "stretch") : nullptr)
{
}

void
ComputeLagrangianObjectiveStress::initQpStatefulProperties()
{
  ComputeLagrangianStressBase::initQpStatefulProperties();

  _small_stress[_qp].zero();
  _cauchy_stress[_qp].zero();

  if (_polar_decomp)
    (*_rotation)[_qp] = RankTwoTensor::Identity();
}

void
ComputeLagrangianObjectiveStress::computeQpCauchyStress()
{
  computeQpSmallStress();

  if (!_large_kinematics)
  {
    _cauchy_stress[_qp] = _small_stress[_qp];
    _cauchy_jacobian[_qp] = _small_jacobian[_qp];
  }
  else
  {
    // If large_kinematics = true, do the objective integration
    RankTwoTensor dS = _small_stress[_qp] - _small_stress_old[_qp];

    if (_rate == ObjectiveRate::Truesdell)
      std::tie(_cauchy_stress[_qp], _cauchy_jacobian[_qp]) = objectiveUpdateTruesdell(dS);
    else if (_rate == ObjectiveRate::Jaumann)
      std::tie(_cauchy_stress[_qp], _cauchy_jacobian[_qp]) = objectiveUpdateJaumann(dS);
    else if (_rate == ObjectiveRate::GreenNaghdi)
      std::tie(_cauchy_stress[_qp], _cauchy_jacobian[_qp]) = objectiveUpdateGreenNaghdi(dS);
    else if (_rate == ObjectiveRate::Rashid)
      std::tie(_cauchy_stress[_qp], _cauchy_jacobian[_qp]) = objectiveUpdateRashid(dS);
    else
      mooseError("Internal error: unsupported objective rate.");
  }
}

std::tuple<RankTwoTensor, RankFourTensor>
ComputeLagrangianObjectiveStress::objectiveUpdateTruesdell(const RankTwoTensor & dS)
{
  // Update the Cauchy stress
  auto [S, Jinv] =
      advectStress(_cauchy_stress_old[_qp] + dS, _spatial_velocity_gradient_increment[_qp]);

  // Get the appropriate tangent tensor
  RankFourTensor U = stressAdvectionDerivative(S);

  return {S, cauchyJacobian(Jinv, U)};
}

std::tuple<RankTwoTensor, RankFourTensor>
ComputeLagrangianObjectiveStress::objectiveUpdateJaumann(const RankTwoTensor & dS)
{
  usingTensorIndices(i, j, k, l);

  // Update the Cauchy stress
  auto [S, Jinv] = advectStress(_cauchy_stress_old[_qp] + dS, _vorticity_increment[_qp]);

  // Get the appropriate tangent tensor
  RankTwoTensor I = RankTwoTensor::Identity();
  RankFourTensor ddW_ddL = 0.5 * (I.times<i, k, j, l>(I) - I.times<i, l, j, k>(I));
  RankFourTensor U = stressAdvectionDerivative(S) * ddW_ddL;

  return {S, cauchyJacobian(Jinv, U)};
}

std::tuple<RankTwoTensor, RankFourTensor>
ComputeLagrangianObjectiveStress::objectiveUpdateGreenNaghdi(const RankTwoTensor & dS)
{
  usingTensorIndices(i, j, k, l, m);

  // The kinematic tensor for the Green-Naghdi rate is
  // Omega = dot(R) R^T
  polarDecomposition();
  RankTwoTensor I = RankTwoTensor::Identity();
  RankTwoTensor dR = (*_rotation)[_qp] * (*_rotation_old)[_qp].transpose() - I;
  RankTwoTensor dO = dR * _inv_df[_qp];

  // Update the Cauchy stress
  auto [S, Jinv] = advectStress(_cauchy_stress_old[_qp] + dS, dO);

  // Get the appropriate tangent tensor
  RankFourTensor d_R_d_F = (*_d_rotation_d_def_grad)[_qp];
  RankFourTensor d_F_d_dL = _inv_df[_qp].inverse().times<i, k, l, j>(_def_grad[_qp]);
  RankTwoTensor T = (*_rotation_old)[_qp].transpose() * _inv_df[_qp];
  RankFourTensor d_dO_d_dL =
      T.times<m, j, i, m, k, l>(d_R_d_F * d_F_d_dL) - dR.times<i, k, j, l>(I);
  RankFourTensor U = stressAdvectionDerivative(S) * d_dO_d_dL;

  return {S, cauchyJacobian(Jinv, U)};
}

std::tuple<RankTwoTensor, RankFourTensor>
ComputeLagrangianObjectiveStress::objectiveUpdateRashid(const RankTwoTensor & dS)
{
  usingTensorIndices(i, j, k, l, m);

  // Rashid does a nonlinear update of the form sigma_new = r (sigma_old + dS) r^T
  // with r = R_new R_old^T
  polarDecomposition(true);
  RankTwoTensor dR = (*_rotation)[_qp];

  auto S = dR * (_cauchy_stress_old[_qp] + dS) * dR.transpose();

  // RankFourTensor Jinv = dR.times<i, k, j, l>(dR);
  //  WIP
  auto jac = RankFourTensor::Identity();

  return {S, jac};
}

std::tuple<RankTwoTensor, RankFourTensor>
ComputeLagrangianObjectiveStress::advectStress(const RankTwoTensor & S0,
                                               const RankTwoTensor & dQ) const
{
  RankFourTensor J = updateTensor(dQ);
  RankFourTensor Jinv = J.inverse();
  RankTwoTensor S = Jinv * S0;
  return {S, Jinv};
}

RankFourTensor
ComputeLagrangianObjectiveStress::updateTensor(const RankTwoTensor & dQ) const
{
  auto I = RankTwoTensor::Identity();
  usingTensorIndices(i, j, k, l);
  return (1.0 + dQ.trace()) * I.times<i, k, j, l>(I) - dQ.times<i, k, j, l>(I) -
         I.times<i, k, j, l>(dQ);
}

RankFourTensor
ComputeLagrangianObjectiveStress::stressAdvectionDerivative(const RankTwoTensor & S) const
{
  auto I = RankTwoTensor::Identity();
  usingTensorIndices(i, j, k, l);
  return S.times<i, j, k, l>(I) - I.times<i, k, l, j>(S) - S.times<i, l, j, k>(I);
}

RankFourTensor
ComputeLagrangianObjectiveStress::cauchyJacobian(const RankFourTensor & Jinv,
                                                 const RankFourTensor & U) const
{
  return Jinv * (_small_jacobian[_qp] - U);
}

void
ComputeLagrangianObjectiveStress::polarDecomposition(bool incremental)
{
  RankTwoTensor use_F;
  if (incremental)
    use_F = _def_grad[_qp] * _def_grad_old[_qp].inverse();
  else
    use_F = _def_grad[_qp];

  FactorizedRankTwoTensor C = use_F.transpose() * use_F;
  (*_stretch)[_qp] = MathUtils::sqrt(C).get();
  RankTwoTensor Uinv = MathUtils::sqrt(C).inverse().get();
  (*_rotation)[_qp] = use_F * Uinv;

  // Derivative of rotation w.r.t. the deformation gradient
  RankTwoTensor I = RankTwoTensor::Identity();
  RankTwoTensor Y = (*_stretch)[_qp].trace() * I - (*_stretch)[_qp];
  RankTwoTensor Z = (*_rotation)[_qp] * Y;
  RankTwoTensor O = Z * (*_rotation)[_qp].transpose();
  usingTensorIndices(i, j, k, l);
  (*_d_rotation_d_def_grad)[_qp] = (O.times<i, k, l, j>(Y) - Z.times<i, l, k, j>(Z)) / Y.det();
}
