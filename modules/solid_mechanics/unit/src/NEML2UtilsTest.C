//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "gtest/gtest.h"

#include "NEML2Utils.h"

#ifdef NEML2_ENABLED
TEST(NEML2Utils, from_blob_Real)
{
  MooseArray<Real> data(3);
  data[0] = 1.0;
  data[1] = 2.0;
  data[2] = 3.0;

  auto tensor = NEML2Utils::from_blob(data);
  ASSERT_TRUE(tensor.defined());
  ASSERT_TRUE(tensor.dim() == 1);
  ASSERT_TRUE(tensor.batch_dim() == 1);
  ASSERT_TRUE(tensor.base_dim() == 0);
  ASSERT_TRUE(tensor.size(0) == 3);

  for (neml2::Size n : index_range(data))
    EXPECT_NEAR(tensor.index({n}).item<Real>(), data[n], 1e-12);
}

TEST(NEML2Utils, from_blob_RealVectorValue)
{
  MooseArray<RealVectorValue> data(3);
  data[0] = RealVectorValue(1.0, 2.0, 3.0);
  data[1] = RealVectorValue(4.0, 5.0, 6.0);
  data[2] = RealVectorValue(7.0, 8.0, 9.0);

  auto tensor = NEML2Utils::from_blob(data);
  ASSERT_TRUE(tensor.defined());
  ASSERT_TRUE(tensor.dim() == 2);
  ASSERT_TRUE(tensor.batch_dim() == 1);
  ASSERT_TRUE(tensor.base_dim() == 1);
  ASSERT_TRUE(tensor.size(0) == 3);
  ASSERT_TRUE(tensor.size(1) == 3);

  for (neml2::Size n : index_range(data))
    for (neml2::Size i : make_range(3))
      EXPECT_NEAR(tensor.index({n, i}).item<Real>(), data[n](i), 1e-12);
}

TEST(NEML2Utils, from_blob_RankTwoTensor)
{
  MooseArray<RankTwoTensor> data(2);
  data[0] = RankTwoTensor(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
  data[1] = RankTwoTensor(-1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0, -9.0);

  auto tensor = NEML2Utils::from_blob(data);
  ASSERT_TRUE(tensor.defined());
  ASSERT_TRUE(tensor.dim() == 3);
  ASSERT_TRUE(tensor.batch_dim() == 1);
  ASSERT_TRUE(tensor.base_dim() == 2);
  ASSERT_TRUE(tensor.size(0) == 2);
  ASSERT_TRUE(tensor.size(1) == 3);
  ASSERT_TRUE(tensor.size(2) == 3);

  for (neml2::Size n : index_range(data))
    for (neml2::Size i : make_range(3))
      for (neml2::Size j : make_range(3))
        EXPECT_NEAR(tensor.index({n, i, j}).item<Real>(), data[n](i, j), 1e-12);
}

TEST(NEML2Utils, from_blob_SymmetricRankTwoTensor)
{
  MooseArray<SymmetricRankTwoTensor> data(2);
  data[0] = SymmetricRankTwoTensor(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  data[1] = SymmetricRankTwoTensor(-1.0, -2.0, -3.0, -4.0, -5.0, -6.0);

  auto tensor = NEML2Utils::from_blob(data);
  ASSERT_TRUE(tensor.defined());
  ASSERT_TRUE(tensor.dim() == 2);
  ASSERT_TRUE(tensor.batch_dim() == 1);
  ASSERT_TRUE(tensor.base_dim() == 1);
  ASSERT_TRUE(tensor.size(0) == 2);
  ASSERT_TRUE(tensor.size(1) == 6);

  for (neml2::Size n : index_range(data))
    for (neml2::Size i : make_range(6))
      EXPECT_NEAR(tensor.index({n, i}).item<Real>(), data[n](i), 1e-12);
}

TEST(NEML2Utils, copyTensorToMooseArray_Real)
{
  MooseArray<Real> data(9);
  const auto tensor = torch::tensor({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}},
                                    torch::TensorOptions().dtype(torch::kFloat64));
  NEML2Utils::copyTensorToMooseArray(tensor, data);

  const auto tensor_flat = tensor.reshape({9});
  for (neml2::Size n : index_range(data))
    EXPECT_NEAR(tensor_flat.index({n}).item<Real>(), data[n], 1e-12);
}

TEST(NEML2Utils, copyTensorToMooseArray_RealVectorValue)
{
  MooseArray<RealVectorValue> data(3);
  const auto tensor = torch::tensor({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}},
                                    torch::TensorOptions().dtype(torch::kFloat64));
  NEML2Utils::copyTensorToMooseArray(tensor, data);

  for (neml2::Size n : index_range(data))
    for (neml2::Size i : make_range(3))
      EXPECT_NEAR(tensor.index({n, i}).item<Real>(), data[n](i), 1e-12);
}

TEST(NEML2Utils, copyTensorToMooseArray_RankTwoTensor)
{
  MooseArray<RankTwoTensor> data(2);
  const auto tensor = torch::tensor({{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}},
                                     {{-1.0, -2.0, -3.0}, {-4.0, -5.0, -6.0}, {-7.0, -8.0, -9.0}}},
                                    torch::TensorOptions().dtype(torch::kFloat64));
  NEML2Utils::copyTensorToMooseArray(tensor, data);

  for (neml2::Size n : index_range(data))
    for (neml2::Size i : make_range(3))
      for (neml2::Size j : make_range(3))
        EXPECT_NEAR(tensor.index({n, i, j}).item<Real>(), data[n](i, j), 1e-12);
}

TEST(NEML2Utils, copyTensorToMooseArray_SymmetricRankTwoTensor)
{
  MooseArray<SymmetricRankTwoTensor> data(2);
  const auto tensor =
      torch::tensor({{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}, {-1.0, -2.0, -3.0, -4.0, -5.0, -6.0}},
                    torch::TensorOptions().dtype(torch::kFloat64));
  NEML2Utils::copyTensorToMooseArray(tensor, data);

  for (neml2::Size n : index_range(data))
    for (neml2::Size i : make_range(6))
      EXPECT_NEAR(tensor.index({n, i}).item<Real>(), data[n](i), 1e-12);
}

TEST(NEML2Utils, copyTensorToMooseArray_SymmetricRankFourTensor)
{
  MooseArray<SymmetricRankFourTensor> data(2);
  const auto tensor = torch::tensor({{{1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
                                      {2.0, 3.0, 4.0, 5.0, 6.0, 7.0},
                                      {3.0, 4.0, 5.0, 6.0, 7.0, 8.0},
                                      {4.0, 5.0, 6.0, 7.0, 8.0, 9.0},
                                      {5.0, 6.0, 7.0, 8.0, 9.0, 10.0},
                                      {6.0, 7.0, 8.0, 9.0, 10.0, 11.0}},
                                     {{-1.0, -2.0, -3.0, -4.0, -5.0, -6.0},
                                      {-2.0, -3.0, -4.0, -5.0, -6.0, -7.0},
                                      {-3.0, -4.0, -5.0, -6.0, -7.0, -8.0},
                                      {-4.0, -5.0, -6.0, -7.0, -8.0, -9.0},
                                      {-5.0, -6.0, -7.0, -8.0, -9.0, -10.0},
                                      {-6.0, -7.0, -8.0, -9.0, -10.0, -11.0}}},
                                    torch::TensorOptions().dtype(torch::kFloat64));
  NEML2Utils::copyTensorToMooseArray(tensor, data);

  for (neml2::Size n : index_range(data))
    for (neml2::Size i : make_range(6))
      for (neml2::Size j : make_range(6))
        EXPECT_NEAR(tensor.index({n, i, j}).item<Real>(), data[n](i, j), 1e-12);
}
#endif
