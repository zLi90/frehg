// Surface-water state (P2.2.3).
#ifndef FREHG2_CORE_STATE_HPP
#define FREHG2_CORE_STATE_HPP

#include "core/Grid.hpp"
#include "frehg2/core/define.hpp"

namespace frehg2 {

// Surface-water state fields (length Grid::nSurfaceCell()). eta (water surface elevation)
// is the PRIMARY unknown solved by the linear system; h = max(0, eta - z) is derived.
class State {
 public:
  explicit State(const Grid& grid);

  const Grid& grid() const { return grid_; }

  RealArr1D eta;  // water surface elevation (primary unknown)
  RealArr1D u;    // depth-averaged x velocity
  RealArr1D v;    // depth-averaged y velocity
  RealArr1D h;    // water depth = max(0, eta - z)
  RealArr1D qss;  // surface-subsurface exchange flux (filled by Coupling)
  // Surface solute concentration (P8). Halo-padded (Grid::getSurfaceIndex) so the transport
  // kernels can read neighbor cells; default 0. Driven by SoluteStepper (wired in P16).
  RealArr1D conc;

 private:
  Grid grid_;
};

}  // namespace frehg2

#endif  // FREHG2_CORE_STATE_HPP
