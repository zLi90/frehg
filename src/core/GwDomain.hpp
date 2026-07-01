// Subsurface (3D) static domain data (P2.2.2).
#ifndef FREHG2_CORE_GW_DOMAIN_HPP
#define FREHG2_CORE_GW_DOMAIN_HPP

#include "core/Grid.hpp"
#include "frehg2/core/define.hpp"

namespace frehg2 {

// Subsurface static fields. Vertical layering is geometric: dz3d[k] = dz * dz_incre^k,
// with k=0 the TOP layer (docs/legacy_audit/index_conventions.md). soilID and z3d are
// per active 3D cell (length Grid::nActiveCell() = nx*ny*nz).
class GwDomain {
 public:
  // bot_z is the elevation of the bottom of the deepest layer (k = nz-1).
  explicit GwDomain(const Grid& grid, real bot_z = 0.0);

  const Grid& grid() const { return grid_; }
  real botZ() const { return bot_z_; }

  RealArr1D dz3d;    // length nz: geometric layer thicknesses (k=0 top)
  IntArr1D soilID;   // length nActiveCell: soil type per cell
  RealArr1D z3d;     // length nActiveCell: cell-center elevation (top to bottom)

 private:
  Grid grid_;
  real bot_z_ = 0.0;
};

}  // namespace frehg2

#endif  // FREHG2_CORE_GW_DOMAIN_HPP
