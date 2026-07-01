// Solute source / sink terms (P8.3.4): rainfall, infiltration mixing, first-order decay.
//
// All operate in place on halo-padded host concentration fields. The rainfall and
// infiltration operators are written in mass-conservative MIXING form (the bounded
// generalization of the legacy first-order add): a wet cell that already holds the inflow
// concentration is left unchanged, so a constant rain of concentration c_rain drives the
// surface concentration to the steady value c_rain.
#ifndef FREHG2_SOLUTE_SOURCE_SINK_HPP
#define FREHG2_SOLUTE_SOURCE_SINK_HPP

#include "core/Grid.hpp"
#include "frehg2/core/define.hpp"
#include "solute/SoluteParams.hpp"

namespace frehg2 {

// Rainfall source on the surface. For each wet cell (depth > min_depth) with rain rate
// `rain` [m/s], mixes incoming water of concentration c_rain into the standing water:
//   C <- (depth*C + dt*rain*c_rain) / (depth + dt*rain).
void applyRainfall(RealArr1DHost& conc_surf, const RealArr1DHost& depth, const Grid& grid,
                   real rain, const SoluteParams& p, real dt);

// Mass-conservative SW<->GW interface solute exchange (P16-completion). Moves dissolved mass
// with the water volume that the coupling just exchanged across the surface / top-GW-cell
// interface, co-located with the water exchange so the post-exchange water heights are exact.
//
// Per column, `vex` is the signed exchanged water height [m] (per unit area):
//   vex > 0  infiltration (SW -> GW): the infiltrating water carries the SURFACE concentration
//            into the top soil cell. Surface conc is unchanged (the water already left at the
//            surface concentration, so removing it does not change it); the top soil cell mixes:
//              C_sub <- (C_sub*(h_sub - vex) + C_surf*vex) / h_sub.
//   vex < 0  seepage (GW -> SW): the seeping water carries the TOP-SOIL concentration into the
//            surface water. Subsurface conc is unchanged; the surface mixes (s = -vex):
//              C_surf <- (C_surf*(depth - s) + C_sub*s) / depth.
//
// `depth_surf` and `h_sub_water` are the CURRENT (post-water-exchange) surface ponded depth and
// top-cell water height (wc*dz), both per unit area. With these the pair (water move + this
// exchange) conserves Sum(C*water) to machine precision (proved in test_source_sink). Columns
// with a non-positive donor height are skipped (nothing to carry).
void applyInterfaceExchange(RealArr1DHost& conc_surf, RealArr1DHost& conc_sub,
                            const RealArr1DHost& depth_surf, const RealArr1DHost& h_sub_water,
                            const RealArr1DHost& vex, const Grid& grid);

// First-order decay over the active surface cells: C <- C * exp(-k_decay*dt).
void applyDecaySurface(RealArr1DHost& conc_surf, const Grid& grid, real k_decay, real dt);

// First-order decay over the active subsurface cells: C <- C * exp(-k_decay*dt).
void applyDecaySubsurface(RealArr1DHost& conc_sub, const Grid& grid, real k_decay, real dt);

}  // namespace frehg2

#endif  // FREHG2_SOLUTE_SOURCE_SINK_HPP
