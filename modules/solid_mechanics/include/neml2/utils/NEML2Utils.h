//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#ifdef NEML2_ENABLED

#include "neml2/misc/parser_utils.h"
#include "neml2/tensors/tensors.h"
#include "neml2/models/LabeledAxisAccessor.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "SymmetricRankTwoTensor.h"
#include "SymmetricRankFourTensor.h"
#include "MaterialProperty.h"

#endif

#include "InputParameters.h"
#include "MooseArray.h"

class MooseObject;
class Action;
class SubProblem;

namespace NEML2Utils
{

#ifdef NEML2_ENABLED
/// Assert that the NEML2 variable name sits on either the forces or the state subaxis
void assertVariable(const neml2::VariableName &);

/// Assert that the NEML2 variable name sits on either the old_forces or the old_state subaxis
void assertOldVariable(const neml2::VariableName &);

/// Parse a raw string into NEML2 variable name
neml2::VariableName parseVariableName(const std::string &);

template <typename T>
struct Layout
{
};
template <>
struct Layout<Real>
{
  static constexpr std::array<neml2::Size, 0> shape{};
  static constexpr std::array<neml2::Size, 1> strides{1};
};
template <>
struct Layout<RealVectorValue>
{
  static constexpr std::array<neml2::Size, 1> shape{3};
  static constexpr std::array<neml2::Size, 2> strides{3, 1};
};
template <>
struct Layout<RankTwoTensor>
{
  static constexpr std::array<neml2::Size, 2> shape{3, 3};
  static constexpr std::array<neml2::Size, 3> strides{9, 3, 1};
};
template <>
struct Layout<SymmetricRankTwoTensor>
{
  static constexpr std::array<neml2::Size, 1> shape{6};
  static constexpr std::array<neml2::Size, 2> strides{6, 1};
};

/**
 * @brief Mapping from MooseArray to neml2::Tensor without copying the data
 *
 * This method is used in gatherers which gather data from MOOSE as input variables to the NEML2
 * material model. So in theory, we only need to overload MOOSE types that can potentially be used
 * as input variables.
 */
template <typename T>
neml2::Tensor
from_blob(const MooseArray<T> & data)
{
  // The const_cast is fine because torch works with non-const ptr so that it can optionally handle
  // deallocation. But we are not going to let torch do that.
  const auto torch_tensor =
      torch::from_blob(const_cast<T *>(data.data()),
                       neml2::utils::add_shapes(data.size(), Layout<T>::shape),
                       Layout<T>::strides,
                       torch::TensorOptions().dtype(torch::kFloat64));
  return neml2::Tensor(torch_tensor, 1);
}

/**
 * @brief Mapping from MooseArray to neml2::Tensor without copying the data
 *
 * Similar to the other from_blob method, but this one is used for a vector of MooseArray.
 */
template <typename T>
neml2::Tensor
from_blob(const std::vector<MooseArray<T>> & data)
{
  std::vector<torch::Tensor> tensors(data.size());
  std::transform(data.begin(),
                 data.end(),
                 tensors.begin(),
                 [](const MooseArray<T> & array) { return from_blob(array); });
  return neml2::Tensor(torch::stack(tensors), 2);
}

template <typename T>
void
copyTensorToMooseArray(const torch::Tensor & src, MooseArray<T> & dest)
{
  mooseAssert(src.numel() == dest.size() * Layout<T>::strides[0],
              "Cannot copy neml2::Tensor into a MooseArray<T> with different number of elements.");

  // memcpy reinterpret the data as unsigned char
  const std::size_t n_unsigned_char = src.numel() * sizeof(Real) / sizeof(unsigned char);

  // This assumes the neml2::Tensor and MooseArray<T> has same layout, for example both row-major (T
  // = RankTwoTensor). If the layouts are different, we may need to reshape the neml2::Tensor before
  // memcpy.
  std::memcpy(dest.data(), src.contiguous().data_ptr(), n_unsigned_char);
}

/// Convert a MOOSE data structure to its NEML2 counterpart
template <typename T>
neml2::Tensor toNEML2(const T &);

/**
 * Convert a wrapped (batched) MOOSE data structure to its NEML2 counterpart
 * The wrapper should implement size()
 */
template <typename T>
neml2::Tensor toNEML2Batched(const T & data);

// @{ Template specializations
template <>
neml2::Tensor toNEML2(const Real & v);
template <>
neml2::Tensor toNEML2(const RankTwoTensor & r2t);
template <>
neml2::Tensor toNEML2(const std::vector<Real> & v);
// @}

/// Convert a NEML2 data structure to its MOOSE counterpart
template <typename T>
T toMOOSE(const neml2::Tensor &);

// @{ Template specializations
template <>
Real toMOOSE(const neml2::Tensor & t);
template <>
SymmetricRankTwoTensor toMOOSE(const neml2::Tensor & t);
template <>
std::vector<Real> toMOOSE(const neml2::Tensor & t);
template <>
SymmetricRankFourTensor toMOOSE(const neml2::Tensor & t);
// @}

static std::string NEML2_help_message = R""""(
==============================================================================
To debug NEML2 related issues:
1. Build and run MOOSE in dbg mode.
2. Re-run the simulation using the dbg executable, and often times
   NEML2 will provide a more helpful error message.
3. If the error message is not helpful, or if there is still no error message,
   run the simulation through a debugger: See
   https://mooseframework.inl.gov/application_development/debugging.html
4. If the issue is due to a NEML2 bug, feel free to report it at
   https://github.com/applied-material-modeling/neml2/issues
==============================================================================
)"""";

////////////////////////////////////////////////////////////////////////////////
// Implementations
////////////////////////////////////////////////////////////////////////////////
template <typename T>
neml2::Tensor
toNEML2Batched(const T & data)
{
  std::vector<torch::Tensor> res(data.size());
  for (const auto i : index_range(data))
    res[i] = toNEML2<typename T::value_type>(data[i]);
  return neml2::Tensor(torch::stack(res, 0), 1);
}

#endif // NEML2_ENABLED

/// Determine whether the NEML2 material model should be evaluated
bool shouldCompute(const SubProblem &);

/**
 * Augment docstring if NEML2 is not enabled
 */
std::string docstring(const std::string & desc);

/**
 * Assert that NEML2 is enabled. A MooseError is raised if NEML2 is not enabled.
 */
void assertNEML2Enabled();

} // namespace NEML2Utils
