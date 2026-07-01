// Apply flexible initial conditions to SWE/RE solver fields (P14.3.1 / 14.3.4).
#ifndef FREHG2_IC_IC_APPLY_HPP
#define FREHG2_IC_IC_APPLY_HPP

#include "frehg2/core/define.hpp"
#include "ic/ICSpec.hpp"
#include "io/Config.hpp"

namespace frehg2 {

class MpiComm;
class SweSolver;
class ReSolver;
class PolygonIndex;

struct ICApplyContext {
  int gnx = 0;
  int gny = 0;
  int gnz = 1;
  double dx = 1.0;
  double dy = 1.0;
  double dz = 1.0;
  double x0 = 0.0;
  double y0 = 0.0;
  double botz = 0.0;
  const MpiComm* mc = nullptr;
  const Config* config = nullptr;
  PolygonIndex* region_index = nullptr;
};

// Apply all parsed IC specs. When ic.use_restart is true, only restart metadata is recorded;
// the Orchestrator loads the checkpoint after solvers are wired.
//
// `surface_conc` / `subsurface_conc` are the canonical halo-padded solute concentration views
// (State::conc indexed by Grid::getSurfaceIndex, GwState::conc by Grid::getIndex). When non-empty
// (solute enabled and the matching module on), the parsed solute ICs (constant/raster/formula/
// polygon) are written into them; empty views skip solute IC. Defaulted so flow-only callers and
// the existing tests are unchanged.
void applyInitialConditions(const InitialConditionsConfig& ic, ICApplyContext ctx, SweSolver* swe,
                            ReSolver* re, RealArr1DHost surface_conc = RealArr1DHost(),
                            RealArr1DHost subsurface_conc = RealArr1DHost());

// Build polygon index for IC regions (no-op when regions empty).
void buildICRegionIndex(const InitialConditionsConfig& ic, ICApplyContext ctx);

}  // namespace frehg2

#endif  // FREHG2_IC_IC_APPLY_HPP
