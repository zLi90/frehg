// Spatially variable soil parameters (P13.3.1, extended to fully heterogeneous 3D in P23).
//
// A SoilMap holds (a) a list of soil CLASSES — each a full `SoilParams` tuple (Ksat, porosity
// theta_s, residual theta_r, van Genuchten alpha/n, specific storage, ...) — and (b) a CLASS
// INDEX raster over this rank's owned cells. Two layouts are supported:
//   * per-COLUMN (2D, P13): index keyed by `(i, j)` in `[0,nx) x [0,ny)`; replicated across all
//     `k` (the column class). `nz_ == 0`.
//   * per-CELL (3D, P23, like SERGHEI): index keyed by `(i, j, k)` in
//     `[0,nx) x [0,ny) x [0,nz)`, row-major `i + j*nx + k*nx*ny`. `nz_ > 0`.
// The RE solver looks up the per-cell class via `SoilMap::classAt(i, j, k)`.
//
// Legacy `frehg` has uniform soil only; this is a Frehg2 feature. The uniform case (one class,
// every cell index 0) reproduces the P5 single-`SoilParams` path bit-for-bit, and a 2D map
// reproduces the P13 per-column path bit-for-bit.
//
// Storage is per-rank and Kokkos/PETSc-free: a plain `std::vector<int>` index + a small vector of
// `SoilParams`. The map is static (built once at startup), so there is no per-step assembly.
#ifndef FREHG2_SOIL_SOIL_MAP_HPP
#define FREHG2_SOIL_SOIL_MAP_HPP

#include <string>
#include <vector>

#include "frehg2/re/SoilParams.hpp"

namespace frehg2 {

class SoilMap {
 public:
  SoilMap() = default;

  // Replace the class table (index 0..nClasses-1). Throws if empty.
  void setClasses(std::vector<SoilParams> classes);

  // Set the per-column class index over a local nx x ny grid (row-major i + j*nx). Every entry
  // must be a valid class id in [0, nClasses). Throws on size mismatch or out-of-range id.
  void setClassIndex(int nx, int ny, std::vector<int> class_idx);

  // Set the fully heterogeneous per-cell (3D) class index over a local nx x ny x nz grid
  // (row-major i + j*nx + k*nx*ny). Every entry must be a valid class id in [0, nClasses).
  // Throws on size mismatch or out-of-range id. (P23, SERGHEI-style per-cell soil.)
  void setClassIndex3D(int nx, int ny, int nz, std::vector<int> class_idx);

  // Uniform single-class map over nx x ny (all cells -> the one class). Convenience for the
  // legacy/uniform path and tests.
  void setUniform(int nx, int ny, const SoilParams& cls);

  // Load a per-cell class-index raster from a whitespace/comma-separated text file containing
  // exactly nx*ny non-negative integers in row-major (i + j*nx) order. Classes must already be
  // set. Handles 1k–1M rows. Throws with the path + counts on any mismatch.
  void loadIndexFromCSV(const std::string& path, int nx, int ny);

  int nClasses() const { return static_cast<int>(classes_.size()); }
  int nx() const { return nx_; }
  int ny() const { return ny_; }
  int nz() const { return nz_; }
  // True iff this is a fully heterogeneous per-cell (3D) map; false for the per-column (2D) map.
  bool is3D() const { return nz_ > 0; }
  bool empty() const { return classes_.empty() || class_idx_.empty(); }

  // Class id for column (i, j); callers should pass owned-interior indices.
  int classIdAt(int i, int j) const {
    return class_idx_[static_cast<size_t>(i) + static_cast<size_t>(j) * static_cast<size_t>(nx_)];
  }
  // Class id for cell (i, j, k). For a 2D map k is ignored (column class is replicated).
  int classIdAt(int i, int j, int k) const {
    if (nz_ == 0) return classIdAt(i, j);
    return class_idx_[static_cast<size_t>(i) + static_cast<size_t>(j) * static_cast<size_t>(nx_) +
                      static_cast<size_t>(k) * static_cast<size_t>(nx_) *
                          static_cast<size_t>(ny_)];
  }
  // Soil parameters for column (i, j).
  const SoilParams& classAt(int i, int j) const {
    return classes_[static_cast<size_t>(classIdAt(i, j))];
  }
  // Soil parameters for cell (i, j, k) (per-cell for a 3D map; column class for a 2D map).
  const SoilParams& classAt(int i, int j, int k) const {
    return classes_[static_cast<size_t>(classIdAt(i, j, k))];
  }
  const SoilParams& classById(int id) const { return classes_[static_cast<size_t>(id)]; }

  // True iff every cell maps to the same class with identical parameters (the bit-identical
  // uniform path). Used by the RE solver only as a diagnostic.
  bool isUniform() const;

 private:
  std::vector<SoilParams> classes_;
  std::vector<int> class_idx_;  // 2D: length nx_*ny_ (i+j*nx_); 3D: nx_*ny_*nz_ (i+j*nx_+k*nx_*ny_)
  int nx_ = 0;
  int ny_ = 0;
  int nz_ = 0;  // 0 => per-column (2D) map; >0 => fully heterogeneous per-cell (3D) map
};

}  // namespace frehg2

#endif  // FREHG2_SOIL_SOIL_MAP_HPP
