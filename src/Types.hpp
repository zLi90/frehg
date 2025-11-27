#ifndef TYPES_HPP
#define TYPES_HPP

#include <Kokkos_Core.hpp>

// -----------------------------------------------------------------------------
// 1. Precision Settings
// -----------------------------------------------------------------------------
// Change this to 'float' if you want to run mixed-precision tests later
using Real = double;
using Ordinal = int; // For indices

// -----------------------------------------------------------------------------
// 2. Templated View Aliases (The "Versatile" Types)
// -----------------------------------------------------------------------------
// We use a template <typename T> so these can hold int, double, or bool.

// 1D View: For packed state variables (Pressure[cell_id], Salinity[cell_id])
template <typename T>
using View1D = Kokkos::View<T*>;

// 2D View: For connectivity tables (Neighbors[cell_id][direction])
template <typename T>
using View2D = Kokkos::View<T**>;

// 3D View: Only used for legacy I/O or temporary structured data
template <typename T>
using View3D = Kokkos::View<T***>;

// -----------------------------------------------------------------------------
// 3. Host Mirrors (For CPU <-> GPU Communication)
// -----------------------------------------------------------------------------
template <typename T>
using View1D_Host = typename View1D<T>::HostMirror;

template <typename T>
using View2D_Host = typename View2D<T>::HostMirror;

template <typename T>
using View3D_Host = typename View3D<T>::HostMirror;

#endif // TYPES_HPP