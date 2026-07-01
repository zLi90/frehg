// Explicit, operator-split advection for the passive scalar (P8.3.2).
//
// Conservative finite-volume transport of a cell-centered concentration by a staggered face
// velocity field (SoluteFlow). First-order upwind (default) or minmod-limited MUSCL (2nd
// order), selected by SoluteParams::advection_scheme. Closed (zero-flux) outer boundaries,
// so the discrete scheme conserves the cell-sum of concentration exactly.
//
// Each function computes the maximum advective Courant number over the domain FIRST. If it
// exceeds SoluteParams::cfl_max the concentration is left UNCHANGED and the (out-of-range)
// Courant number is returned so the caller can refuse the step. Otherwise the update is
// applied and the (in-range) Courant number is returned.
#ifndef FREHG2_SOLUTE_ADVECTION_HPP
#define FREHG2_SOLUTE_ADVECTION_HPP

#include "core/Grid.hpp"
#include "frehg2/core/define.hpp"
#include "solute/SoluteFlow.hpp"
#include "solute/SoluteParams.hpp"

namespace frehg2 {

// Surface (2D) advection of a halo-padded concentration field. Returns max Courant number.
real advectSurface(RealArr1DHost& conc, const SoluteFlow& flow, const Grid& grid,
                   const SoluteParams& p, real dt);

// Subsurface (3D) advection of a halo-padded concentration field. Returns max Courant number.
real advectSubsurface(RealArr1DHost& conc, const SoluteFlow& flow, const Grid& grid,
                      const SoluteParams& p, real dt);

}  // namespace frehg2

#endif  // FREHG2_SOLUTE_ADVECTION_HPP
