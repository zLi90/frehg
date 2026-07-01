// Surface (2D) static domain data (P2.2.1).
#ifndef FREHG2_CORE_DOMAIN_HPP
#define FREHG2_CORE_DOMAIN_HPP

#include "core/Grid.hpp"
#include "frehg2/core/define.hpp"

namespace frehg2 {

// Per-surface-cell static fields (length Grid::nSurfaceCell() = nx*ny, no halo, per the
// P2.2 acceptance). z is bed elevation, area is cell plan area, actMask is 1 for active
// cells (0 otherwise), roughness is Manning's n.
class Domain {
 public:
  explicit Domain(const Grid& grid);

  const Grid& grid() const { return grid_; }

  RealArr1D z;
  RealArr1D area;
  RealArr1D actMask;
  RealArr1D roughness;

 private:
  Grid grid_;
};

}  // namespace frehg2

#endif  // FREHG2_CORE_DOMAIN_HPP
