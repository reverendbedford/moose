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
#ifdef MOOSE_ENABLE_KOKKOS_GPU
#ifdef KOKKOS_ENABLE_CUDA
  // Keep NEML2 computation and output retrieval on the device where Kokkos consumes the output
  const at::Device kokkos_device(at::kCUDA, ExecSpace().cuda_device());
  const auto & neml2_device = neml2Device();
  if (neml2_device.type() != kokkos_device.type() ||
      (neml2_device.has_index() && neml2_device.index() != kokkos_device.index()))
    paramError("neml2_executor",
               "The NEML2 model compute device ",
               neml2_device,
               " does not match the Kokkos device ",
               kokkos_device,
               ".");
  requestNEML2OutputDevice(kokkos_device);
#else
  mooseError("NEML2ToKokkosMaterialProperty supports GPU execution only with CUDA.");
#endif
#else
  const at::Device kokkos_device(at::kCPU);
  if (neml2Device() != kokkos_device)
    paramError("neml2_executor",
               "The NEML2 model compute device ",
               neml2Device(),
               " does not match the Kokkos device ",
               kokkos_device,
               ".");
  requestNEML2OutputDevice(kokkos_device);
#endif
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
#ifdef KOKKOS_ENABLE_CUDA
    ::Kokkos::deep_copy(_staged_output.deviceView(), 0.0);
#else
    mooseError("NEML2ToKokkosMaterialProperty supports GPU execution only with CUDA.");
#endif
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

#ifdef MOOSE_ENABLE_KOKKOS_GPU
#ifdef KOKKOS_ENABLE_CUDA
    if (!neml2Output().is_cuda())
      mooseError("Cannot copy NEML2 output on ",
                 neml2Output().device(),
                 " into a CUDA Kokkos material property. Set the NEML2 executor output_device to "
                 "the CUDA device used by Kokkos.");
    const auto kokkos_device = ExecSpace().cuda_device();
    if (neml2Output().device().index() != kokkos_device)
      mooseError("Cannot copy NEML2 output on ",
                 neml2Output().device(),
                 " into a Kokkos material property on cuda:",
                 kokkos_device,
                 ". The NEML2 and Kokkos CUDA device indices must match.");
#else
    mooseError("NEML2ToKokkosMaterialProperty supports GPU execution only with CUDA.");
#endif
#else
    if (!neml2Output().is_cpu())
      mooseError("Cannot copy NEML2 output on ",
                 neml2Output().device(),
                 " into a CPU Kokkos material property. Set the NEML2 executor output_device to "
                 "cpu.");
#endif

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
#ifdef KOKKOS_ENABLE_CUDA
    const auto output_view =
        ::Kokkos::View<const Real *, MemSpace, ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>(
            output.template data_ptr<Real>(), output.numel());
    ::Kokkos::deep_copy(_staged_output.deviceView(), output_view);
#endif
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
