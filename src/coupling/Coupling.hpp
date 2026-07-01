// Synchronous surface-water / groundwater coupling (P6) — Frehg's original sync coupling
// from legacy/frehg/src/solve.c + the z-back top-face exchange of darcy_flux()
// (legacy/frehg/src/subroutines.c:70-112).
//
// The asynchronous SERGHEI-style coupling is P11. This module is the simplest, fully
// mass-conservative synchronous exchange: at each coupled step the SWE solver advances,
// the SW<->GW exchange flux is computed per column from the updated surface state and the
// current top GW head, and the SAME flux volume is moved out of one domain and into the
// other before the GW solver advances. Conservation is structural (what leaves SW enters
// GW), independent of either solver's internal tolerances.
//
// No PETSc types appear here: the coupling drives the SweSolver / ReSolver through their
// public field/state accessors and lets each solver own its own LinearSolver/SparseSystem
// seam. Per DP3, the coupling loops stay plain serial/host; on-node Kokkos/GPU is the
// global P9/P10 pass.
#ifndef FREHG2_COUPLING_COUPLING_HPP
#define FREHG2_COUPLING_COUPLING_HPP

#include "frehg2/core/define.hpp"

namespace frehg2 {

class SweSolver;
class ReSolver;

struct CouplingParams {
  real dx = 1.0;
  real dy = 1.0;
  real visc = 1.0;        // relative viscosity at the z+ face (legacy r_visczp; 1 if non-baroclinic)
  real min_depth = 1.0e-8;  // ponded-water threshold (legacy param->min_dept)
};

class Coupling {
 public:
  explicit Coupling(const CouplingParams& p) : params_(p) {}

  const CouplingParams& params() const { return params_; }

  // Per-column SW<->GW exchange flux [m^3/s] across the surface/top-GW interface.
  //
  //   q = Az * kface * visc * (h_gw_top - h_surface) / delta,
  //       Az = dx*dy,  kface = Ksz (saturated, surface water present),  delta = 0.5*dz_top.
  //
  // Sign convention (legacy darcy_flux z-back; plan P6.1):
  //   q > 0  -> seepage      (GW -> SW): top GW head exceeds the surface head.
  //   q < 0  -> infiltration (SW -> GW): surface head exceeds the top GW head.
  //   q == 0 at h_gw_top == h_surface (total-head form; no separate gravity term, so the
  //          exchange vanishes exactly at hydrostatic balance and the sign tests are clean).
  //
  // When there is no ponded surface water (dept <= min_depth) an infiltration gradient has
  // nothing to draw from, so q is clamped to >= 0 (only seepage out of the GW is possible).
  // The magnitude is set by Ksz (saturated conductivity) — infiltration is K_sat-limited.
  real columnExchangeRate(real h_gw_top, real h_surface, real dept, real Ksz,
                          real dz_top) const;

  // Full grid: one independent flux per surface column (no cross-column coupling). Writes
  // `q_owned`, length nx*ny, row-major owned index (i + j*nx). The interface heads are
  // datum-free (legacy groundwater.c:804): the surface ponded depth (eta-bottom) is the head
  // imposed above the top GW cell, compared against that cell's pressure head.
  void computeExchangeRates(const SweSolver& swe, const ReSolver& re,
                            RealArr1DHost& q_owned) const;

  // Apply the exchange to the surface: eta += q*dt/Az per column (seepage raises the surface,
  // infiltration lowers it). Updates dept for immediate consistency; the next SWE step
  // re-derives depth/geometry from eta.
  void applyExchangeToSurface(const RealArr1DHost& q_owned, SweSolver& swe, real dt) const;

  // Apply the exchange to the groundwater: top-cell wc -= q*dt/V per column (seepage removes
  // GW water, infiltration adds it). Keeps the top-cell head consistent with the new water
  // content and clamps to [theta_r, theta_s].
  void applyExchangeToGroundwater(const RealArr1DHost& q_owned, ReSolver& re, real dt) const;

  // Compute + conservatively limit + apply the SW<->GW exchange (steps 2-4 of
  // stepCoupled) for the coupling interval `dt`, WITHOUT advancing either solver. Used by
  // the Orchestrator to subcycle the groundwater solver: the exchange volume is moved once
  // per surface step, then the GW solver takes one or more substeps to catch up. Returns
  // the total exchanged volume [m^3] this interval (signed: + = net GW->SW seepage).
  //
  // When `q_out` is non-null it is overwritten with the per-column LIMITED flux [m^3/s]
  // (length nx*ny, owned index i + j*nx, sign as columnExchangeRate) that was actually
  // applied to the water. The solute coupling (Orchestrator) uses this to move dissolved
  // mass with exactly the water volume that crossed the interface, so the SW<->GW solute
  // exchange is mass-conservative by construction.
  real exchange(SweSolver& swe, ReSolver& re, real dt, RealArr1DHost* q_out = nullptr) const;

  // One synchronous coupled step (legacy solve.c sync_coupling==1 ordering):
  //   1. swe.advanceStep(rain, evap)
  //   2. q = computeExchangeRates(swe, re), conservatively clamped by both domains' available
  //      water so the exchanged volume never drives a negative depth or wc outside
  //      [theta_r, theta_s]
  //   3. applyExchangeToSurface(q)   (h_sw += q*dt/Az)
  //   4. applyExchangeToGroundwater(q) (wc_top -= q*dt/V)
  //   5. re.advanceStep()
  // Returns the total exchanged volume [m^3] this step (signed: + = net GW->SW seepage).
  real stepCoupled(SweSolver& swe, ReSolver& re, real rain_rate, real evap_rate,
                   real dt) const;

 private:
  // Clamp `q_owned` in place so |q*dt| never exceeds the water available on the donor side:
  // infiltration (q<0) by surface ponded volume, seepage (q>0) by GW drainable volume.
  void limitExchange(RealArr1DHost& q_owned, const SweSolver& swe, const ReSolver& re,
                     real dt) const;

  CouplingParams params_;
};

}  // namespace frehg2

#endif  // FREHG2_COUPLING_COUPLING_HPP
