#include "monitoring/LineFlux.hpp"

#include <algorithm>
#include <cmath>

namespace frehg2 {

namespace {

double sampleSurfaceComponent(const LineFluxSpec& line, double px, double py, int gnx, int gny,
                             double dx, double dy, double x0, double y0, const SweSolver& swe,
                             const Grid& grid, const MpiComm* mc) {
  const int gi = std::max(0, std::min(gnx - 1, static_cast<int>(std::floor((px - x0) / dx))));
  const int gj = std::max(0, std::min(gny - 1, static_cast<int>(std::floor((py - y0) / dy))));
  int li = gi;
  int lj = gj;
  if (mc != nullptr) {
    const auto [lx, ly, owner] = mc->globalToLocal(gi, gj);
    if (owner != mc->rank()) return 0.0;
    li = lx;
    lj = ly;
  }
  const int c = grid.getSurfaceIndex(li, lj);
  const SweFields& f = swe.fields();
  if (line.field == "u" || line.field == "uu") return f.uu(c);
  if (line.field == "v" || line.field == "vv") return f.vv(c);
  if (line.field == "flux_u") return f.Fu(c);
  if (line.field == "flux_v") return f.Fv(c);
  // Default: normal component q·n.
  const double tx = line.p1[0] - line.p0[0];
  const double ty = line.p1[1] - line.p0[1];
  const double len = std::hypot(tx, ty);
  if (len <= 0.0) return 0.0;
  const double nx = -ty / len;
  const double ny = tx / len;
  return f.uu(c) * nx + f.vv(c) * ny;
}

double sampleGwComponent(const LineFluxSpec& line, double px, double py, int gnx, int gny,
                         int gk, const ReSolver& re, const Grid& grid, const MpiComm* mc) {
  const int gi = std::max(0, std::min(gnx - 1, static_cast<int>(std::floor(px))));
  const int gj = std::max(0, std::min(gny - 1, static_cast<int>(std::floor(py))));
  int li = gi;
  int lj = gj;
  if (mc != nullptr) {
    const auto [lx, ly, owner] = mc->globalToLocal(gi, gj);
    if (owner != mc->rank()) return 0.0;
    li = lx;
    lj = ly;
  }
  const int c = grid.getIndex(li, lj, gk);
  const GwFields& f = re.fields();
  if (line.field == "qx") return f.qx(c);
  if (line.field == "qy") return f.qy(c);
  if (line.field == "qz") return f.qz(c);
  return f.qx(c);
}

}  // namespace

double integrateLineFlux(const LineFluxSpec& line, const SweSolver* swe, const ReSolver* re,
                         const Grid& swe_grid, const Grid& gw_grid, int gnx, int gny, double dx,
                         double dy, double x0, double y0, const MpiComm* mc) {
  const double x_a = line.p0[0];
  const double y_a = line.p0[1];
  const double x_b = line.p1[0];
  const double y_b = line.p1[1];
  const double seg_len = std::hypot(x_b - x_a, y_b - y_a);
  if (seg_len <= 0.0) return 0.0;

  const int nseg =
      std::max(2, static_cast<int>(std::ceil(seg_len / std::min(dx, dy))) * 2);
  double sum = 0.0;
  double prev = 0.0;
  for (int s = 0; s <= nseg; ++s) {
    const double t = static_cast<double>(s) / static_cast<double>(nseg);
    const double px = x_a + t * (x_b - x_a);
    const double py = y_a + t * (y_b - y_a);
    double val = 0.0;
    if (line.subsurface) {
      if (re == nullptr) return 0.0;
      val = sampleGwComponent(line, px, py, gnx, gny, 0, *re, gw_grid, mc);
    } else {
      if (swe == nullptr) return 0.0;
      val = sampleSurfaceComponent(line, px, py, gnx, gny, dx, dy, x0, y0, *swe, swe_grid, mc);
    }
    if (s == 0) {
      prev = val;
      continue;
    }
    const double dl = seg_len / static_cast<double>(nseg);
    sum += 0.5 * (prev + val) * dl;
    prev = val;
  }
  return sum;
}

}  // namespace frehg2
