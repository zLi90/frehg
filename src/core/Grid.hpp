// Frehg2 structured grid with halo-padded flat indexing (P2.1).
//
// Frehg2 stores fields in a halo-padded flat layout with one ghost cell on every
// horizontal side (x, y). The vertical direction (z, subsurface) has NO horizontal-style
// halo padding; z neighbors are on-rank (see MpiComm::haloExchange3D). This differs from
// legacy Frehg's "interior cells first, ghosts appended" map-based layout; the two are
// bridged for regression by LegacyIndexAdapter (include/frehg2/core/LegacyIndexAdapter.hpp).
//
// Index math is KOKKOS_INLINE_FUNCTION so it is callable from device kernels (Grid is a
// small value type captured by copy).
#ifndef FREHG2_CORE_GRID_HPP
#define FREHG2_CORE_GRID_HPP

#include <tuple>

#include "frehg2/core/define.hpp"
#include "frehg2/core/types.hpp"

namespace frehg2 {

class Grid {
 public:
  Grid() = default;
  Grid(int nx, int ny, int nz, real dx, real dy, real dz, real dz_incre = 1.0);
  explicit Grid(const DomainParams& p);

  KOKKOS_INLINE_FUNCTION int nx() const { return nx_; }
  KOKKOS_INLINE_FUNCTION int ny() const { return ny_; }
  KOKKOS_INLINE_FUNCTION int nz() const { return nz_; }
  KOKKOS_INLINE_FUNCTION int nzActual() const { return nz_; }
  KOKKOS_INLINE_FUNCTION real dx() const { return dx_; }
  KOKKOS_INLINE_FUNCTION real dy() const { return dy_; }
  KOKKOS_INLINE_FUNCTION real dz() const { return dz_; }
  KOKKOS_INLINE_FUNCTION real dzIncre() const { return dz_incre_; }

  // Storage extents.
  KOKKOS_INLINE_FUNCTION int nxHalo() const { return nx_ + 2; }
  KOKKOS_INLINE_FUNCTION int nyHalo() const { return ny_ + 2; }

  // Total storage for a 3D subsurface field: horizontal halo plus one ghost layer on
  // each vertical face (k in [-1, nz] addressable for GW BC storage).
  KOKKOS_INLINE_FUNCTION int nCell() const {
    return (nx_ + 2) * (ny_ + 2) * (nz_ + 2);
  }
  // Physical surface cells, no halo.
  KOKKOS_INLINE_FUNCTION int nSurfaceCell() const { return nx_ * ny_; }
  // 2D surface field storage with halo.
  KOKKOS_INLINE_FUNCTION int nSurfaceStorageCell() const { return (nx_ + 2) * (ny_ + 2); }
  // Physical 3D cells only.
  KOKKOS_INLINE_FUNCTION int nActiveCell() const { return nx_ * ny_ * nz_; }

  // Flat 3D index with horizontal and vertical halo offset:
  // (i+1) + (j+1)*(nx+2) + (k+1)*(nx+2)*(ny+2), k in [-1, nz].
  KOKKOS_INLINE_FUNCTION int getIndex(int i, int j, int k) const {
    return (i + 1) + (j + 1) * (nx_ + 2) + (k + 1) * (nx_ + 2) * (ny_ + 2);
  }
  // Flat 2D surface index with halo offset: (i+1) + (j+1)*(nx+2).
  KOKKOS_INLINE_FUNCTION int getSurfaceIndex(int i, int j) const {
    return (i + 1) + (j + 1) * (nx_ + 2);
  }

  // Inverse of getIndex: returns physical (i,j,k) (halo cells yield i==-1 or i==nx, etc.).
  KOKKOS_INLINE_FUNCTION void getIJK(int flat, int& i, int& j, int& k) const {
    const int plane = (nx_ + 2) * (ny_ + 2);
    k = flat / plane - 1;
    const int rem = flat - (k + 1) * plane;
    j = rem / (nx_ + 2) - 1;
    i = rem % (nx_ + 2) - 1;
  }
  std::tuple<int, int, int> getIJK(int flat) const {
    int i = 0;
    int j = 0;
    int k = 0;
    getIJK(flat, i, j, k);
    return {i, j, k};
  }

  // Inverse of getSurfaceIndex.
  KOKKOS_INLINE_FUNCTION void getIJ(int flat, int& i, int& j) const {
    j = flat / (nx_ + 2) - 1;
    i = flat % (nx_ + 2) - 1;
  }

  KOKKOS_INLINE_FUNCTION bool isActive(int i, int j, int k) const {
    return i >= 0 && i < nx_ && j >= 0 && j < ny_ && k >= 0 && k < nz_;
  }
  KOKKOS_INLINE_FUNCTION bool isSurfaceActive(int i, int j) const {
    return i >= 0 && i < nx_ && j >= 0 && j < ny_;
  }

 private:
  int nx_ = 1;
  int ny_ = 1;
  int nz_ = 1;
  real dx_ = 1.0;
  real dy_ = 1.0;
  real dz_ = 1.0;
  real dz_incre_ = 1.0;
};

}  // namespace frehg2

#endif  // FREHG2_CORE_GRID_HPP
