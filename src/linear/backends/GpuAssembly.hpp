// P10 — compile-time selector for the Kokkos -> PETSc assembly bridge option.
//
// `FREHG2_GPU_ASSEMBLY` is 1 iff the build uses a Kokkos GPU backend (CUDA / HIP / SYCL),
// which is the only situation in which PETSc can consume device-resident COO buffers via
// `MatSetValuesCOO` / `VecGetArrayAndMemType` without host staging (Option B, Task 10.3.3/4).
//
// On the macOS OpenMP/Serial build the macro is 0, so every device path below is compiled
// out and the backend uses the host-staged Option A (`MatSetValues` on HostSpace), which is
// bit-identical to the P9 reference. The macro is derived purely from the Kokkos
// configuration macros (pulled in via <Kokkos_Core.hpp>), so it tracks the actual Kokkos
// device the binary was built for — there is no separate runtime switch.
//
// This header lives under src/linear/backends/ (the only directory where PETSc/solver types
// are permitted by the solver seam); physics code never includes it.
#ifndef FREHG2_LINEAR_BACKENDS_GPU_ASSEMBLY_HPP
#define FREHG2_LINEAR_BACKENDS_GPU_ASSEMBLY_HPP

#include <Kokkos_Core.hpp>

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_SYCL)
#define FREHG2_GPU_ASSEMBLY 1
#else
#define FREHG2_GPU_ASSEMBLY 0
#endif

#endif  // FREHG2_LINEAR_BACKENDS_GPU_ASSEMBLY_HPP
