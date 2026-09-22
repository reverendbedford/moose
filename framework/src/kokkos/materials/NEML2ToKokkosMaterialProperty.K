//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "NEML2ToKokkosMaterialProperty.h"
#include "NEML2Utils.h"

#include <algorithm>

#if defined(NEML2_ENABLED) && defined(MOOSE_ENABLE_KOKKOS_GPU) && defined(KOKKOS_ENABLE_SYCL)
// c10::xpu device queries, used to match the Kokkos SYCL device to a torch XPU device index.
#include <c10/xpu/XPUFunctions.h>
#endif

#ifdef NEML2_ENABLED
namespace
{
/**
 * Torch device backing the Kokkos execution space.
 *
 * NEML2 outputs must be staged on this device before Kokkos can consume them: the unmanaged Views
 * built over the output tensors in computeProperties() are only addressable from the Kokkos
 * execution space if the tensor memory lives on the same device.
 */
at::Device
kokkosTorchDevice()
{
#ifndef MOOSE_ENABLE_KOKKOS_GPU
  return at::Device(at::kCPU);
#elif defined(KOKKOS_ENABLE_CUDA)
  return at::Device(at::kCUDA, ExecSpace().cuda_device());
#elif defined(KOKKOS_ENABLE_SYCL)
  // A USM pointer is only valid within the SYCL context it was allocated in, so an unmanaged View
  // over a torch allocation is addressable from the Kokkos queue only if the two share a context.
  // Both torch and Kokkos take the platform default context today, so this holds; check it rather
  // than let deep_copy() fault if that ever stops being true.
  const auto & queue = ExecSpace().sycl_queue();
  if (queue.get_context() != c10::xpu::get_device_context())
    ::mooseError("Kokkos and torch are using different SYCL contexts, so NEML2 output tensors are "
                 "not addressable from the Kokkos execution space.");

  // torch indexes its XPU devices independently of Kokkos, so map through the SYCL device itself
  // rather than assuming the two enumerations agree.
  const auto & kokkos_device = queue.get_device();
  for (c10::DeviceIndex i = 0; i < c10::xpu::device_count(); ++i)
    if (c10::xpu::get_raw_device(i) == kokkos_device)
      return at::Device(at::kXPU, i);

  ::mooseError("The SYCL device Kokkos is running on is not among the ",
               static_cast<int>(c10::xpu::device_count()),
               " XPU devices visible to torch.");
#else
  ::mooseError("NEML2ToKokkosMaterialProperty supports GPU execution only with the CUDA and SYCL "
               "Kokkos backends.");
#endif
}
}
#endif

registerKokkosMaterial("MooseApp", NEML2ToKokkosRankTwoMaterialProperty);
registerKokkosMaterial("MooseApp", NEML2ToKokkosRankFourMaterialProperty);

template <unsigned int rank>
InputParameters
NEML2ToKokkosMaterialProperty<rank>::validParams()
{
  auto params = Material::validParams();
  params += NEML2OutputInterface::validParams();
  params.addClassDescription("Provides a symmetric rank-" + std::to_string(rank) +
                             " NEML2 output or derivative as a Kokkos material property.");
  params.addRequiredParam<MaterialPropertyName>(
      "to_moose", "Kokkos material property used to store the NEML2 output or derivative.");
  return params;
}

template <unsigned int rank>
NEML2ToKokkosMaterialProperty<rank>::NEML2ToKokkosMaterialProperty(const InputParameters & params)
  : Material(params),
    NEML2OutputInterface(params, *this),
    _prop(declareKokkosProperty<Real, rank>("to_moose",
                                            std::vector<unsigned int>(rank, _mesh.dimension())))
{
  NEML2Utils::assertNEML2Enabled();
#ifdef NEML2_ENABLED
  // NEML2 may compute on another device, but its batched outputs must be transferred to the
  // device where Kokkos consumes them.
  requestNEML2OutputDevice(kokkosTorchDevice());
#endif
}

template <unsigned int rank>
NEML2ToKokkosMaterialProperty<rank>::NEML2ToKokkosMaterialProperty(
    const NEML2ToKokkosMaterialProperty & object)
  : Material(object),
    NEML2OutputInterface(object),
    _prop(object._prop),
    _staged_output(object._staged_output),
    _element_batch_offset(object._element_batch_offset),
    _batch_mapping_outdated(object._batch_mapping_outdated),
    _broadcast_output(object._broadcast_output),
    _batch_qp_stride(object._batch_qp_stride)
{
}

template <unsigned int rank>
void
NEML2ToKokkosMaterialProperty<rank>::meshChanged()
{
  _batch_mapping_outdated = true;
}

template <unsigned int rank>
void
NEML2ToKokkosMaterialProperty<rank>::computeProperties()
{
#ifdef NEML2_ENABLED
  if (!outputReady())
    return;

  constexpr unsigned int num_components = rank == 2 ? 6 : 36;
  const auto num_local_elements = kokkosMesh().getNumLocalElements();

  if (!neml2Output().defined())
  {
    if (_staged_output.size() != num_components || !_broadcast_output)
    {
      _compute_dispatcher.reset();
      _staged_output.create(num_components);
      _broadcast_output = true;
    }
#ifdef MOOSE_ENABLE_KOKKOS_GPU
    ::Kokkos::deep_copy(_staged_output.deviceView(), 0.0);
#else
    std::fill_n(_staged_output.hostData(), num_components, 0);
#endif
  }
  else
  {
    const auto base_ndim = static_cast<int64_t>(rank / 2);
    const auto batch_ndim = neml2Output().dim() - base_ndim;
    const auto expected_shape = rank == 2 ? std::vector<int64_t>{6} : std::vector<int64_t>{6, 6};
    const auto actual_shape = neml2Output().sizes().slice(batch_ndim);
    if (!actual_shape.equals(expected_shape))
      mooseError("Cannot copy NEML2 tensor with shape ",
                 neml2Output().sizes(),
                 " into a symmetric rank-",
                 rank,
                 " Kokkos material property; expected base shape ",
                 expected_shape,
                 ".");
    if (neml2Output().dtype() != at::kDouble)
      mooseError("Cannot copy NEML2 tensor with dtype ",
                 neml2Output().dtype(),
                 " into a Real Kokkos material property.");

    if (batch_ndim < 0 || batch_ndim > 2)
      mooseError("Cannot copy NEML2 tensor with shape ",
                 neml2Output().sizes(),
                 " into a Kokkos material property; expected zero, one, or two batch dimensions.");

    const auto kokkos_device = kokkosTorchDevice();
    if (neml2Output().device() != kokkos_device)
      mooseError("Cannot copy NEML2 output on ",
                 neml2Output().device(),
                 " into a Kokkos material property on ",
                 kokkos_device,
                 ". Set the NEML2 executor output_device to ",
                 kokkos_device,
                 ".");

    const auto output = neml2Output().contiguous();
    const auto batch_size = static_cast<std::size_t>(output.numel() / num_components);
    const bool broadcast_output = batch_ndim == 0;
    if (_staged_output.size() != static_cast<std::size_t>(output.numel()) ||
        _broadcast_output != broadcast_output)
    {
      _compute_dispatcher.reset();
      _staged_output.create(output.numel());
      _broadcast_output = broadcast_output;
    }
#ifdef MOOSE_ENABLE_KOKKOS_GPU
    const auto output_view =
        ::Kokkos::View<const Real *, MemSpace, ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>(
            output.template data_ptr<Real>(), output.numel());
    ::Kokkos::deep_copy(_staged_output.deviceView(), output_view);
#else
    std::copy_n(output.template data_ptr<Real>(), output.numel(), _staged_output.hostData());
#endif

    const auto batch_qp_stride = batch_ndim == 2 ? output.size(1) : 0;
    if (_batch_mapping_outdated || _element_batch_offset.size() != num_local_elements ||
        _batch_qp_stride != batch_qp_stride)
    {
      if (_element_batch_offset.size() != num_local_elements)
      {
        _compute_dispatcher.reset();
        _element_batch_offset.create(num_local_elements);
      }

      for (const auto elem : _mesh.getMesh().active_local_element_ptr_range())
      {
        const auto elem_id = kokkosMesh().getContiguousElementID(elem);
        _element_batch_offset.hostData()[elem_id] =
            batch_ndim == 1 ? neml2BatchIndex(elem->id()) : elem_id * batch_qp_stride;

        const auto num_qps = kokkosAssembly().getNumQpsHost(elem_id);
        if (!_broadcast_output && _element_batch_offset.hostData()[elem_id] + num_qps > batch_size)
          mooseError("Element ",
                     elem->id(),
                     " maps outside NEML2 output tensor with shape ",
                     output.sizes(),
                     ".");
      }

      _element_batch_offset.copyToDevice();
      _batch_mapping_outdated = false;
      _batch_qp_stride = batch_qp_stride;
    }
  }

#ifndef MOOSE_ENABLE_KOKKOS_GPU
  _staged_output.copyToDevice();
#endif
  Material::computeProperties();
#else
  NEML2Utils::assertNEML2Enabled();
#endif
}

template class NEML2ToKokkosMaterialProperty<2>;
template class NEML2ToKokkosMaterialProperty<4>;
