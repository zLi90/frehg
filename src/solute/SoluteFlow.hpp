// Flow fields read by the solute transport solver (P8).
//
// The solute solver is a passive scalar carried by an externally-computed flow field. In a
// production run (P16) these views alias the SWE/RE solver fields; in the P8 unit tests they
// are synthesized directly. All arrays are HOST, halo-padded in the same layout the rest of
// Frehg2 uses (Grid::getSurfaceIndex / Grid::getIndex), so the transport kernels can read
// neighbor / face values without copies. Empty views mean "that domain is not driven".
//
// Face conventions (staggered, matching the SWE/RE face layout):
//   surface  u(c) : x-velocity on the EAST face of cell c  (between c and (i+1,j)) [m/s]
//   surface  v(c) : y-velocity on the NORTH face of cell c (between c and (i,j+1)) [m/s]
//   subsurf qx(c) : Darcy velocity on the EAST face  [m/s]
//   subsurf qy(c) : Darcy velocity on the NORTH face [m/s]
//   subsurf qz(c) : Darcy velocity on the DOWN face  (between k and k+1) [m/s]
#ifndef FREHG2_SOLUTE_SOLUTE_FLOW_HPP
#define FREHG2_SOLUTE_SOLUTE_FLOW_HPP

#include "frehg2/core/define.hpp"

namespace frehg2 {

struct SoluteFlow {
  // Surface (length Grid::nSurfaceStorageCell()).
  RealArr1DHost u;
  RealArr1DHost v;
  RealArr1DHost depth;

  // Subsurface (length Grid::nCell()).
  RealArr1DHost qx;
  RealArr1DHost qy;
  RealArr1DHost qz;
  RealArr1DHost wc;

  bool hasSurface() const { return depth.extent(0) > 0; }
  bool hasSubsurface() const { return wc.extent(0) > 0; }
};

}  // namespace frehg2

#endif  // FREHG2_SOLUTE_SOLUTE_FLOW_HPP
