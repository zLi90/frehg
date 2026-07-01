#include "bc/PolygonBC.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

#include "swe/SweSolver.hpp"

namespace frehg2 {

PolygonBC::PolygonBC(std::vector<BcRegion> regions, PolygonIndex index, const MpiComm* mc)
    : regions_(std::move(regions)), index_(std::move(index)) {
  global_count_ = index_.globalColumnCounts(mc);
}

real PolygonBC::applySurface(SweSolver& swe, real dt) const {
  if (regions_.empty()) return 0.0;
  SweFields& f = swe.fields();
  const Grid& g = swe.grid();
  const int nx = g.nx();
  const int ny = g.ny();
  const real dx = g.dx();
  const real dy = g.dy();
  const real area = dx * dy;
  const real gravity = swe.params().gravity;

  real net_outflow = 0.0;
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const int c = g.getSurfaceIndex(i, j);
      const int p = index_.lookup(c);
      if (p < 0) continue;
      const BcRegion& r = regions_[static_cast<size_t>(p)];
      const real depth = std::max(f.eta(c) - f.bottom(c), static_cast<real>(0.0));

      switch (r.kind) {
        case BcKind::Depth: {
          // Prescribe the ponded depth (override the solved eta). Net volume removed = drop*area.
          const real h_bc = std::max(r.rate, static_cast<real>(0.0));
          const real old_eta = f.eta(c);
          f.eta(c) = f.bottom(c) + h_bc;
          net_outflow += (old_eta - f.eta(c)) * area;
          break;
        }
        case BcKind::Discharge: {
          // Prescribed volumetric flux Q [m^3/s] (>0 outflow), spread evenly over the region's
          // global cell count. Outflow is limited by the available ponded volume; inflow (Q<0)
          // is unlimited (it adds water).
          const int gc = global_count_[static_cast<size_t>(p)];
          if (gc <= 0) break;
          real vol = r.rate * dt / static_cast<real>(gc);
          const real avail = depth * area;
          if (vol > avail) vol = avail;
          f.eta(c) -= vol / area;
          net_outflow += vol;
          break;
        }
        case BcKind::Critical: {
          // Critical-depth weir outflow: q = sqrt(g) h^{3/2} per unit width, taken over a
          // critical section of width min(dx, dy). Limited by the available ponded volume.
          const real q = std::sqrt(gravity) * std::pow(depth, static_cast<real>(1.5)) *
                         std::min(dx, dy);
          real vol = q * dt;
          const real avail = depth * area;
          if (vol > avail) vol = avail;
          f.eta(c) -= vol / area;
          net_outflow += vol;
          break;
        }
      }
    }

  // Refresh depth + subgrid geometry so the next step starts from a consistent state (mirrors
  // the end-of-advanceStep refresh after the explicit rain/evap update).
  swe.updateDepth();
  swe.updateGeometry();
  return net_outflow;
}

}  // namespace frehg2
