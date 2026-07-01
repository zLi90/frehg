#include "coupling/Coupling.hpp"

#include <algorithm>
#include <cmath>

#include "frehg2/core/ParallelFor.hpp"
#include "re/ReSolver.hpp"
#include "re/VanGenuchten.hpp"
#include "swe/SweSolver.hpp"

namespace frehg2 {

real Coupling::columnExchangeRate(real h_gw_top, real h_surface, real dept, real Ksz,
                                  real dz_top) const {
  const real Az = params_.dx * params_.dy;
  const real delta = 0.5 * dz_top;
  const real kface = Ksz;  // surface water present -> saturated top-face conductivity
  real q = Az * kface * params_.visc * (h_gw_top - h_surface) / delta;
  // No ponded surface water: an infiltration gradient (q<0) has nothing to draw from, so
  // only seepage out of the groundwater (q>=0) is physical.
  if (dept <= params_.min_depth && q < 0.0) q = 0.0;
  return q;
}

void Coupling::computeExchangeRates(const SweSolver& swe, const ReSolver& re,
                                    RealArr1DHost& q_owned) const {
  const Grid sg = swe.grid();
  const Grid gg = re.grid();
  const int nx = sg.nx();
  const int ny = sg.ny();
  const SweFields& sf = swe.fields();
  const GwFields& gf = re.fields();
  const SoilParams& soil = re.params().soil;
  const real Ksz = soil.Ks_z;
  const real Az = params_.dx * params_.dy;
  const real visc = params_.visc;
  const real min_depth = params_.min_depth;
  RealArr1DHost eta = sf.eta;
  RealArr1DHost bottom = sf.bottom;
  RealArr1DHost h = gf.h;
  RealArr1DHost dz3d = gf.dz3d;
  RealArr1DHost q = q_owned;

  // Datum-free interface heads (legacy groundwater.c:804): the surface ponded depth `dept`
  // is the head imposed in the ghost above the top GW cell, compared against that cell's
  // pressure head `h`. Both are pressure-like heads at the interface, so no elevation datum
  // is needed and the flux vanishes when they balance. This is the inlined z-back top-face
  // form of columnExchangeRate (per-cell, disjoint writes -> Kokkos parallel-for, P9).
  parallelForSurface("coupling_exchange_rate", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const int si = sg.getSurfaceIndex(i, j);
    const int ti = gg.getIndex(i, j, 0);
    const real dept = Kokkos::max(eta(si) - bottom(si), 0.0);
    const real dz_top = dz3d(ti);
    real qv = Az * Ksz * visc * (h(ti) - dept) / (0.5 * dz_top);
    if (dept <= min_depth && qv < 0.0) qv = 0.0;
    q(static_cast<size_t>(i + j * nx)) = qv;
  });
}

void Coupling::limitExchange(RealArr1DHost& q_owned, const SweSolver& swe, const ReSolver& re,
                             real dt) const {
  const Grid sg = swe.grid();
  const Grid gg = re.grid();
  const int nx = sg.nx();
  const int ny = sg.ny();
  const SweFields& sf = swe.fields();
  const GwFields& gf = re.fields();
  const SoilParams& soil = re.params().soil;
  const real Az = params_.dx * params_.dy;
  const real theta_s = soil.theta_s;
  const real theta_r = soil.theta_r;
  RealArr1DHost eta = sf.eta;
  RealArr1DHost bottom = sf.bottom;
  RealArr1DHost wc = gf.wc;
  RealArr1DHost dz3d = gf.dz3d;
  RealArr1DHost q = q_owned;

  parallelForSurface("coupling_limit_exchange", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const size_t oi = static_cast<size_t>(i + j * nx);
    real qv = q(oi);
    if (qv == 0.0) return;
    const int si = sg.getSurfaceIndex(i, j);
    const int ti = gg.getIndex(i, j, 0);
    const real V = Az * dz3d(ti);
    const real wc_top = wc(ti);
    if (qv < 0.0) {
      // Infiltration SW -> GW: limited by surface ponded volume and GW headroom.
      const real dept = Kokkos::max(eta(si) - bottom(si), 0.0);
      const real q_min_surf = -dept * Az / dt;
      const real q_min_gw = -(theta_s - wc_top) * V / dt;
      qv = Kokkos::max(qv, Kokkos::max(q_min_surf, q_min_gw));
    } else {
      // Seepage GW -> SW: limited by GW drainable volume.
      const real q_max_gw = (wc_top - theta_r) * V / dt;
      qv = Kokkos::min(qv, q_max_gw);
    }
    q(oi) = qv;
  });
}

void Coupling::applyExchangeToSurface(const RealArr1DHost& q_owned, SweSolver& swe,
                                      real dt) const {
  const Grid sg = swe.grid();
  const int nx = sg.nx();
  const int ny = sg.ny();
  SweFields& sf = swe.fields();
  const real Az = params_.dx * params_.dy;
  RealArr1DHost eta = sf.eta;
  RealArr1DHost bottom = sf.bottom;
  RealArr1DHost dept = sf.dept;
  RealArr1DHost q = q_owned;

  parallelForSurface("coupling_apply_surface", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const int si = sg.getSurfaceIndex(i, j);
    const real dEta = q(static_cast<size_t>(i + j * nx)) * dt / Az;
    eta(si) += dEta;
    dept(si) = Kokkos::max(eta(si) - bottom(si), 0.0);
  });
}

void Coupling::applyExchangeToGroundwater(const RealArr1DHost& q_owned, ReSolver& re,
                                          real dt) const {
  const Grid gg = re.grid();
  const int nx = gg.nx();
  const int ny = gg.ny();
  GwFields& gf = re.fields();
  const SoilParams soil = re.params().soil;
  const real Az = params_.dx * params_.dy;
  RealArr1DHost wc = gf.wc;
  RealArr1DHost h = gf.h;
  RealArr1DHost dz3d = gf.dz3d;
  RealArr1DHost q = q_owned;

  parallelForSurface("coupling_apply_gw", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const int ti = gg.getIndex(i, j, 0);
    const real V = Az * dz3d(ti);
    // q>0 seepage removes GW water; q<0 infiltration adds it.
    real wc_new = wc(ti) - q(static_cast<size_t>(i + j * nx)) * dt / V;
    wc_new = Kokkos::min(Kokkos::max(wc_new, soil.theta_r), soil.theta_s);
    wc(ti) = wc_new;
    h(ti) = VanGenuchten::headFromWaterContent(soil, wc_new);
  });
}

real Coupling::exchange(SweSolver& swe, ReSolver& re, real dt, RealArr1DHost* q_out) const {
  const int nx = swe.grid().nx();
  const int ny = swe.grid().ny();
  RealArr1DHost q("q_exchange", static_cast<size_t>(nx * ny));
  computeExchangeRates(swe, re, q);
  limitExchange(q, swe, re, dt);
  // Export the limited flux BEFORE applying it (apply only mutates solver fields, not q) so the
  // caller's solute coupling can move dissolved mass with exactly the volume that crosses.
  if (q_out != nullptr) {
    if (q_out->extent(0) != q.extent(0))
      *q_out = RealArr1DHost("coupling_flux", q.extent(0));
    Kokkos::deep_copy(*q_out, q);
  }
  applyExchangeToSurface(q, swe, dt);
  applyExchangeToGroundwater(q, re, dt);

  real total = 0.0;
  for (int idx = 0; idx < nx * ny; ++idx) total += q(static_cast<size_t>(idx)) * dt;
  return total;
}

real Coupling::stepCoupled(SweSolver& swe, ReSolver& re, real rain_rate, real evap_rate,
                           real dt) const {
  swe.advanceStep(rain_rate, evap_rate);
  const real total = exchange(swe, re, dt);
  re.advanceStep();
  return total;
}

}  // namespace frehg2
