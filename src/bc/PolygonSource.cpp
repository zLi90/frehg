#include "bc/PolygonSource.hpp"

#include <algorithm>
#include <utility>

#include "re/ReSolver.hpp"
#include "swe/SweSolver.hpp"

namespace frehg2 {

PolygonSource::PolygonSource(std::vector<SourceRegion> regions, PolygonIndex index,
                             const MpiComm* mc)
    : regions_(std::move(regions)), index_(std::move(index)) {
  global_count_ = index_.globalColumnCounts(mc);
  for (const SourceRegion& r : regions_) {
    if (r.kind == SourceKind::InflowRate || r.kind == SourceKind::RainfallRate)
      has_surface_ = true;
    if (r.kind == SourceKind::ExtractionWell) has_subsurface_ = true;
  }
}

real PolygonSource::applySurface(SweSolver& swe, real dt) const {
  if (!has_surface_) return 0.0;
  SweFields& f = swe.fields();
  const Grid& g = swe.grid();
  const int nx = g.nx();
  const int ny = g.ny();
  const real area = g.dx() * g.dy();

  real net_added = 0.0;
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const int c = g.getSurfaceIndex(i, j);
      const int p = index_.lookup(c);
      if (p < 0) continue;
      const SourceRegion& r = regions_[static_cast<size_t>(p)];
      if (r.kind == SourceKind::InflowRate) {
        const int gc = global_count_[static_cast<size_t>(p)];
        if (gc <= 0) continue;
        const real vol = r.rate * dt / static_cast<real>(gc);  // m^3 to this cell
        f.eta(c) += vol / area;
        net_added += vol;
      } else if (r.kind == SourceKind::RainfallRate) {
        const real dvol = r.rate * dt * area;
        f.eta(c) += r.rate * dt;
        net_added += dvol;
      }
    }

  swe.updateDepth();
  swe.updateGeometry();
  return net_added;
}

real PolygonSource::applySubsurface(ReSolver& re, real dt) const {
  if (!has_subsurface_) return 0.0;
  GwFields& f = re.fields();
  const Grid& g = re.grid();
  const int nx = g.nx();
  const int ny = g.ny();
  const int kdeep = g.nz() - 1;
  const real dx = g.dx();
  const real dy = g.dy();

  real net_removed = 0.0;
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const int surf = g.getSurfaceIndex(i, j);
      const int p = index_.lookup(surf);
      if (p < 0) continue;
      const SourceRegion& r = regions_[static_cast<size_t>(p)];
      if (r.kind != SourceKind::ExtractionWell) continue;
      const int gc = global_count_[static_cast<size_t>(p)];
      if (gc <= 0) continue;
      const int c = g.getIndex(i, j, kdeep);
      const real cell_vol = dx * dy * f.dz3d(c);
      if (cell_vol <= 0.0) continue;
      // Per-cell residual moisture (respects a P13/P23 heterogeneous soil map; falls back to the
      // uniform soil when no map is attached). Using params().soil here would clamp against the
      // wrong theta_r under a soil map.
      const real theta_r = re.soilParamsAt(i, j, kdeep).theta_r;
      const real extractable = std::max(f.wc(c) - theta_r, static_cast<real>(0.0)) * cell_vol;
      real vol = r.rate * dt / static_cast<real>(gc);  // m^3 to remove from this column (>0)
      if (vol > extractable) vol = extractable;
      if (vol < 0.0) vol = 0.0;  // a "well" only extracts
      f.wc(c) -= vol / cell_vol;
      net_removed += vol;
    }
  return net_removed;
}

}  // namespace frehg2
