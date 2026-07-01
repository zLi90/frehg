// Bridges Frehg2 halo-padded flat storage to legacy Frehg interior-only ordering (P2.1).
//
// Legacy Frehg stores interior cells first in sequential order and appends ghost cells at
// the end, with neighbor access via precomputed maps. Frehg2 stores fields halo-padded
// (one ghost on each horizontal side). To compare Frehg2 fields against legacy output you
// MUST reorder Frehg2 storage into legacy interior order with these helpers; a direct
// element-by-element index comparison would be wrong.
//
// Conventions (from docs/legacy_audit/index_conventions.md):
//   legacy surface interior index   : i + j*nx
//   legacy subsurface interior index: (i + j*nx)*nz + k   (k=0 is top layer)
// Frehg2 storage indices use the halo-padded formulas (see Grid):
//   surface storage index   : (i+1) + (j+1)*(nx+2)
//   subsurface storage index: (i+1) + (j+1)*(nx+2) + k*(nx+2)*(ny+2)
#ifndef FREHG2_CORE_LEGACY_INDEX_ADAPTER_HPP
#define FREHG2_CORE_LEGACY_INDEX_ADAPTER_HPP

#include "frehg2/core/define.hpp"

namespace frehg2 {

class LegacyIndexAdapter {
 public:
  // Legacy interior orderings.
  static int legacySurfaceInteriorIndex(int i, int j, int nx) { return i + j * nx; }
  static int legacySubsurfaceInteriorIndex(int i, int j, int k, int nx, int nz) {
    return (i + j * nx) * nz + k;
  }

  // Frehg2 halo-padded storage indices (kept here so the adapter is self-contained).
  static int frehg2SurfaceStorageIndex(int i, int j, int nx) {
    return (i + 1) + (j + 1) * (nx + 2);
  }
  static int frehg2SubsurfaceStorageIndex(int i, int j, int k, int nx, int ny) {
    return (i + 1) + (j + 1) * (nx + 2) + (k + 1) * (nx + 2) * (ny + 2);
  }

  // Halo-padded surface field -> legacy interior order (length nx*ny).
  static RealArr1DHost toLegacySurfaceOrder(const RealArr1DHost& halo, int nx, int ny) {
    RealArr1DHost out("legacy_surface", static_cast<size_t>(nx) * ny);
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        out(legacySurfaceInteriorIndex(i, j, nx)) =
            halo(frehg2SurfaceStorageIndex(i, j, nx));
      }
    }
    return out;
  }

  // Legacy interior order -> halo-padded surface field (length (nx+2)*(ny+2); ghosts 0).
  static RealArr1DHost fromLegacySurfaceOrder(const RealArr1DHost& legacy, int nx, int ny) {
    RealArr1DHost out("halo_surface", static_cast<size_t>(nx + 2) * (ny + 2));
    Kokkos::deep_copy(out, 0.0);
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        out(frehg2SurfaceStorageIndex(i, j, nx)) =
            legacy(legacySurfaceInteriorIndex(i, j, nx));
      }
    }
    return out;
  }

  // Halo-padded subsurface field -> legacy interior order (length nx*ny*nz).
  static RealArr1DHost toLegacySubsurfaceOrder(const RealArr1DHost& halo, int nx, int ny,
                                               int nz) {
    RealArr1DHost out("legacy_subsurface", static_cast<size_t>(nx) * ny * nz);
    for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          out(legacySubsurfaceInteriorIndex(i, j, k, nx, nz)) =
              halo(frehg2SubsurfaceStorageIndex(i, j, k, nx, ny));
        }
      }
    }
    return out;
  }

  // Legacy interior order -> halo-padded subsurface field (length (nx+2)*(ny+2)*nz).
  static RealArr1DHost fromLegacySubsurfaceOrder(const RealArr1DHost& legacy, int nx,
                                                 int ny, int nz) {
    RealArr1DHost out("halo_subsurface", static_cast<size_t>(nx + 2) * (ny + 2) * (nz + 2));
    Kokkos::deep_copy(out, 0.0);
    for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          out(frehg2SubsurfaceStorageIndex(i, j, k, nx, ny)) =
              legacy(legacySubsurfaceInteriorIndex(i, j, k, nx, nz));
        }
      }
    }
    return out;
  }
};

}  // namespace frehg2

#endif  // FREHG2_CORE_LEGACY_INDEX_ADAPTER_HPP
