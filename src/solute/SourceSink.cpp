#include "solute/SourceSink.hpp"

#include <cmath>

#include "frehg2/core/ParallelFor.hpp"

namespace frehg2 {

void applyRainfall(RealArr1DHost& conc_surf, const RealArr1DHost& depth, const Grid& grid,
                   real rain, const SoluteParams& p, real dt) {
  if (rain <= 0.0 || dt <= 0.0) return;
  const int nx = grid.nx();
  const int ny = grid.ny();
  const real added = dt * rain;  // rainfall water height entering this step
  const Grid g = grid;
  const real c_rain = p.c_rain;
  const real min_depth = p.min_depth;
  RealArr1DHost cs = conc_surf;
  RealArr1DHost dep = depth;
  parallelForSurface("solute_rainfall", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const size_t c = static_cast<size_t>(g.getSurfaceIndex(i, j));
    const real h = dep(c);
    if (h <= min_depth) return;
    cs(c) = (h * cs(c) + added * c_rain) / (h + added);
  });
}

void applyInterfaceExchange(RealArr1DHost& conc_surf, RealArr1DHost& conc_sub,
                            const RealArr1DHost& depth_surf, const RealArr1DHost& h_sub_water,
                            const RealArr1DHost& vex, const Grid& grid) {
  const int nx = grid.nx();
  const int ny = grid.ny();
  const Grid g = grid;
  RealArr1DHost cs = conc_surf;
  RealArr1DHost cg = conc_sub;
  RealArr1DHost ds = depth_surf;
  RealArr1DHost hg = h_sub_water;
  RealArr1DHost ve = vex;
  parallelForSurface("solute_interface_exchange", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const size_t csurf = static_cast<size_t>(g.getSurfaceIndex(i, j));
    const size_t ctop = static_cast<size_t>(g.getIndex(i, j, 0));
    const real v = ve(csurf);
    if (v > 0.0) {
      // Infiltration SW -> GW: carry the surface concentration into the top soil cell.
      const real h = hg(csurf);
      if (h > 0.0) {
        real h0 = h - v;  // pre-exchange soil water height (>= 0 by the coupling limiter)
        if (h0 < 0.0) h0 = 0.0;
        cg(ctop) = (cg(ctop) * h0 + cs(csurf) * v) / h;
      }
    } else if (v < 0.0) {
      // Seepage GW -> SW: carry the top-soil concentration into the surface water.
      const real s = -v;
      const real d = ds(csurf);
      if (d > 0.0) {
        real d0 = d - s;  // pre-exchange surface depth (>= 0 by the coupling limiter)
        if (d0 < 0.0) d0 = 0.0;
        cs(csurf) = (cs(csurf) * d0 + cg(ctop) * s) / d;
      }
    }
  });
}

void applyDecaySurface(RealArr1DHost& conc_surf, const Grid& grid, real k_decay, real dt) {
  if (k_decay <= 0.0 || dt <= 0.0) return;
  const real fac = std::exp(-k_decay * dt);
  const int nx = grid.nx();
  const int ny = grid.ny();
  const Grid g = grid;
  RealArr1DHost cs = conc_surf;
  parallelForSurface("solute_decay_surf", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    cs(static_cast<size_t>(g.getSurfaceIndex(i, j))) *= fac;
  });
}

void applyDecaySubsurface(RealArr1DHost& conc_sub, const Grid& grid, real k_decay, real dt) {
  if (k_decay <= 0.0 || dt <= 0.0) return;
  const real fac = std::exp(-k_decay * dt);
  const int nx = grid.nx();
  const int ny = grid.ny();
  const int nz = grid.nz();
  const Grid g = grid;
  RealArr1DHost cs = conc_sub;
  parallelForVolume("solute_decay_subs", nx, ny, nz, KOKKOS_LAMBDA(int i, int j, int k) {
    cs(static_cast<size_t>(g.getIndex(i, j, k))) *= fac;
  });
}

}  // namespace frehg2
