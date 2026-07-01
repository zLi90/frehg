// Subsurface (groundwater) state (P2.2.4).
#ifndef FREHG2_CORE_GW_STATE_HPP
#define FREHG2_CORE_GW_STATE_HPP

#include "core/Grid.hpp"
#include "frehg2/core/define.hpp"

namespace frehg2 {

// Subsurface state fields (length Grid::nActiveCell()). The GW algorithm is a predictor-
// corrector (PCA): one implicit head solve (predictor) then a flux-based water-content
// corrector (docs/legacy_audit/state_variables.md). h/hn = current/previous head,
// wc/wcn = current/previous water content, h_pred/wc_pred = predictor intermediates.
class GwState {
 public:
  explicit GwState(const Grid& grid);

  const Grid& grid() const { return grid_; }

  RealArr1D h;    // hydraulic head (current)
  RealArr1D hn;   // hydraulic head (previous)
  RealArr1D wc;   // water content (current)
  RealArr1D wcn;  // water content (previous)

  RealArr1D Kx;  // face conductivity x
  RealArr1D Ky;  // face conductivity y
  RealArr1D Kz;  // face conductivity z

  RealArr1D qx;  // Darcy flux at x face
  RealArr1D qy;  // Darcy flux at y face
  RealArr1D qz;  // Darcy flux at z face

  RealArr1D h_pred;   // predictor-stage head
  RealArr1D wc_pred;  // predictor-stage water content

  // Subsurface solute concentration (P8). Halo-padded (Grid::getIndex) for neighbor access;
  // default 0. Driven by SoluteStepper (wired in P16).
  RealArr1D conc;

 private:
  Grid grid_;
};

}  // namespace frehg2

#endif  // FREHG2_CORE_GW_STATE_HPP
