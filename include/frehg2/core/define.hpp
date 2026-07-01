// Frehg2 core type definitions (P1.3).
//
// Central place for the floating-point typedef, Kokkos View aliases (device + host),
// index type, NODATA sentinel, and the small scoped enums shared across modules.
#ifndef FREHG2_CORE_DEFINE_HPP
#define FREHG2_CORE_DEFINE_HPP

#include <cstdint>

#include <Kokkos_Core.hpp>

namespace frehg2 {

// Floating-point precision for the whole model. Double precision is required for the
// linear-solver tolerances used by the SWE/RE gates.
using real = double;

// Signed index type for global/local cell numbering.
using index_t = std::int32_t;

// Sentinel for missing raster / NODATA values.
inline constexpr real NO_DATA_REAL = -9999.0;

// Execution / memory spaces. On macOS the default execution space is Serial or OpenMP;
// on a CUDA build it becomes Cuda. Device views live in the default space's memory;
// host mirrors live in HostSpace.
using DeviceExec = Kokkos::DefaultExecutionSpace;
using DeviceSpace = Kokkos::DefaultExecutionSpace::memory_space;
using HostSpace = Kokkos::HostSpace;

// Device View aliases (LayoutRight = row-major, C order).
using RealArr1D = Kokkos::View<real*, Kokkos::LayoutRight, DeviceSpace>;
using RealArr2D = Kokkos::View<real**, Kokkos::LayoutRight, DeviceSpace>;
using RealArr3D = Kokkos::View<real***, Kokkos::LayoutRight, DeviceSpace>;
using IntArr1D = Kokkos::View<int*, Kokkos::LayoutRight, DeviceSpace>;
using IntArr2D = Kokkos::View<int**, Kokkos::LayoutRight, DeviceSpace>;

// Host-space View aliases (explicit HostSpace; used for I/O staging and the
// Kokkos->PETSc CPU bridge via MatSetValues on HostSpace data).
using RealArr1DHost = Kokkos::View<real*, Kokkos::LayoutRight, HostSpace>;
using RealArr2DHost = Kokkos::View<real**, Kokkos::LayoutRight, HostSpace>;
using RealArr3DHost = Kokkos::View<real***, Kokkos::LayoutRight, HostSpace>;
using IntArr1DHost = Kokkos::View<int*, Kokkos::LayoutRight, HostSpace>;
using IntArr2DHost = Kokkos::View<int**, Kokkos::LayoutRight, HostSpace>;

// High-level simulation mode (which physics modules are active).
enum class SimMode { SW_ONLY, GW_ONLY, COUPLED, SOLUTE };

// Boundary-condition kinds (unified surface + subsurface vocabulary). Maps onto the
// legacy GW integer codes (0->NEUMANN/zero-flux, 1->FIXED_HEAD, 2->FIXED_FLUX) and the
// position-based SW behaviours documented in docs/legacy_audit/bc_code_reference.md.
enum class BCType { DIRICHLET, NEUMANN, FREEFLOW, ZEROGRADIENT, FIXED_HEAD, FIXED_FLUX };

}  // namespace frehg2

#endif  // FREHG2_CORE_DEFINE_HPP
