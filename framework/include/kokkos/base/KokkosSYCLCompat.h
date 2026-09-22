//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// Force-included by kokkos.mk ahead of every other header when building the
// Kokkos sources with the SYCL backend.
//
// DPC++ only defines SYCL_EXTERNAL during the host pass; in the device pass it
// deliberately leaves the macro to the user and spells its own uses
// __DPCPP_SYCL_EXTERNAL. Third-party headers (libTorch's
// torch/headeronly/macros/Macros.h among them) declare device-side externs with
// a bare SYCL_EXTERNAL guarded on __SYCL_DEVICE_ONLY__, which then fails to
// parse. Defining it up front to the attribute DPC++ would have used keeps
// __DPCPP_SYCL_EXTERNAL identical to its default while making the documented
// spelling available in both passes.
#if defined(__SYCL_DEVICE_ONLY__) && !defined(SYCL_EXTERNAL)
#define SYCL_EXTERNAL __attribute__((sycl_device))
#endif
