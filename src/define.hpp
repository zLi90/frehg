#ifndef FREHG_DEFINE_HPP
#define FREHG_DEFINE_HPP

#include <Kokkos_Core.hpp>
#include <cmath>
#include <cstdint>

// ============================================================================
//                               PRECISION SETTINGS
// ============================================================================

// Define the floating-point precision. 
// Using an alias allows you to switch between double and float easily.
// Standard hydrology models use double, but mixed-precision is useful for GPUs.
using Scalar = double;

// Backward compatibility alias (used in existing code)
using Real = Scalar;

// Define the integer type for mesh indexing and loop counters.
// 'int' is usually sufficient, but 'int64_t' may be needed for massive meshes.
// Can be changed to 'long long' or 'std::size_t' for very large meshes.
using Ordinal = int;

// Define the unsigned ordinal type for array sizes and offsets
using SizeType = std::size_t;

// ============================================================================
//                               KOKKOS VIEW DEFINITIONS
// ============================================================================

// Kokkos Views are reference-counted, multi-dimensional arrays.
// They automatically handle memory layout (Row-Major vs Column-Major) 
// based on the compilation target (CPU vs GPU) for optimal performance.

// 1. Device Views (For computation on GPU or multicore CPU)
// ----------------------------------------------------------------------------
// T corresponds to data type (Scalar or Ordinal)
template <typename T>
using View1D = Kokkos::View<T*>;

template <typename T>
using View2D = Kokkos::View<T**>;

template <typename T>
using View3D = Kokkos::View<T***>;

template <typename T>
using View4D = Kokkos::View<T****>;

// 2. Host Mirror Views (For I/O, MPI communication, and initialization)
// ----------------------------------------------------------------------------
// These views always reside in Host memory (RAM).
// Use Kokkos::deep_copy(host_view, device_view) to transfer data.

template <typename T>
using View1DHost = typename View1D<T>::HostMirror;

template <typename T>
using View2DHost = typename View2D<T>::HostMirror;

template <typename T>
using View3DHost = typename View3D<T>::HostMirror;

template <typename T>
using View4DHost = typename View4D<T>::HostMirror;

// ============================================================================
//                           EXECUTION & MEMORY SPACES
// ============================================================================

// Execution Spaces - Define where kernels run
// ----------------------------------------------------------------------------
// DefaultExecutionSpace adapts to available backends (CUDA, OpenMP, Serial, etc.)
using ExecSpace = Kokkos::DefaultExecutionSpace;

// Host execution space (always available, runs on CPU)
// Uses Serial as fallback, but can be OpenMP if enabled
#ifdef KOKKOS_ENABLE_OPENMP
using HostExecSpace = Kokkos::OpenMP;
#else
using HostExecSpace = Kokkos::Serial;
#endif

// Memory Spaces - Define where data resides
// ----------------------------------------------------------------------------
// Default memory space (matches DefaultExecutionSpace)
using MemSpace = Kokkos::DefaultExecutionSpace::memory_space;

// Host memory space (always available, CPU RAM)
using HostMemSpace = Kokkos::HostSpace;

// Device memory space (GPU memory if available, otherwise same as Host)
#ifdef KOKKOS_ENABLE_CUDA
using DeviceMemSpace = Kokkos::CudaSpace;
#elif defined(KOKKOS_ENABLE_HIP)
using DeviceMemSpace = Kokkos::HIPSpace;
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
using DeviceMemSpace = Kokkos::Experimental::OpenMPTargetSpace;
#else
using DeviceMemSpace = Kokkos::HostSpace;  // Fallback to host
#endif

// Unified Memory Space (for systems with unified memory, e.g., CUDA UVM)
#ifdef KOKKOS_ENABLE_CUDA
using UnifiedMemSpace = Kokkos::CudaUVMSpace;
#else
using UnifiedMemSpace = Kokkos::HostSpace;  // Fallback to host
#endif

// ============================================================================
//                           PARALLEL EXECUTION POLICIES
// ============================================================================

// Range policy for 1D parallel loops
// Usage: Kokkos::parallel_for(RangePolicy(0, N), KOKKOS_LAMBDA(int i) { ... });
using RangePolicy = Kokkos::RangePolicy<ExecSpace>;

// MDRange policy for multi-dimensional loops (2D, 3D)
// Usage: Kokkos::parallel_for(MDRangePolicy<Rank<2>>({0,0}, {nx,ny}), 
//                              KOKKOS_LAMBDA(int i, int j) { ... });
template <int Rank>
using MDRangePolicy = Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<Rank>>;

// Team policy for hierarchical parallelism (thread teams)
// Usage: Kokkos::parallel_for(TeamPolicy(N, Kokkos::AUTO), 
//                              KOKKOS_LAMBDA(const TeamMember& team) { ... });
using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;

// ============================================================================
//                           CONSTANTS & GLOBALS
// ============================================================================

namespace Constants {

    // --- Physical Constants ---
    // Using constexpr allows the compiler to optimize these fully.
    
    // Gravitational acceleration (m/s^2)
    static constexpr Scalar g = 9.80665;
    
    // Density of water (kg/m^3) - Standard reference
    static constexpr Scalar rho_w = 1000.0;
    
    // Density of air (kg/m^3)
    static constexpr Scalar rho_a = 1.225;
    
    // Von Karman constant
    static constexpr Scalar kappa = 0.41;

    // --- Numerical Constants ---
    
    // Pi
    static constexpr Scalar pi = 3.14159265358979323846;
    
    // Small number to prevent division by zero
    static constexpr Scalar epsilon = 1.0e-12;
    
    // Large number for initialization
    static constexpr Scalar big_num = 1.0e20;
    
    // Tolerance for nonlinear solvers
    static constexpr Scalar tol_newton = 1.0e-6;
    
    // --- Conversion Factors ---
    static constexpr Scalar rad_to_deg = 180.0 / pi;
    static constexpr Scalar deg_to_rad = pi / 180.0;
}

// ============================================================================
//                           KOKKOS UTILITIES
// ============================================================================

// Note: KOKKOS_LAMBDA is already defined by Kokkos itself
// For CPU backends (Serial, OpenMP), it expands to [=]
// For GPU backends (CUDA, HIP), it includes device annotations
// Usage: Kokkos::parallel_for(N, KOKKOS_LAMBDA(int i) { ... });

#endif // FREHG_DEFINE_HPP