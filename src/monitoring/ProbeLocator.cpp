#include "monitoring/ProbeLocator.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace frehg2 {

namespace {

int clampIndex(int v, int n) { return std::max(0, std::min(n - 1, v)); }

void xyzToGlobal(double x, double y, double z, int gnx, int gny, int gnz, double dx, double dy,
                 double dz, double x0, double y0, double botz, int& gi, int& gj, int& gk) {
  gi = clampIndex(static_cast<int>(std::floor((x - x0) / dx)), gnx);
  gj = clampIndex(static_cast<int>(std::floor((y - y0) / dy)), gny);
  const double z_top = botz;
  const double z_bot = botz - static_cast<double>(gnz) * dz;
  const double z_clamped = std::max(z_bot + 0.5 * dz, std::min(z_top - 0.5 * dz, z));
  gk = clampIndex(static_cast<int>(std::floor((z_top - z_clamped) / dz)), gnz);
}

}  // namespace

void buildProbeLocations(MonitorBundle& bundle, int gnx, int gny, int gnz, double dx, double dy,
                         double dz, double x0, double y0, double botz, const MpiComm* mc) {
  if (gnx <= 0 || gny <= 0 || gnz <= 0)
    throw std::runtime_error("ProbeLocator: domain dimensions must be positive");

  for (ProbeSpec& p : bundle.probes) {
    if (p.coords_from_xyz) {
      xyzToGlobal(p.x, p.y, p.z, gnx, gny, p.subsurface ? gnz : 1, dx, dy, dz, x0, y0, botz,
                  p.gi, p.gj, p.gk);
    }
    p.gi = clampIndex(p.gi, gnx);
    p.gj = clampIndex(p.gj, gny);
    if (p.subsurface) {
      p.gk = clampIndex(p.gk, gnz);
    } else {
      p.gk = 0;
    }

    if (mc != nullptr) {
      const auto [li, lj, owner] = mc->globalToLocal(p.gi, p.gj);
      p.owner_rank = owner;
      p.local_i = li;
      p.local_j = lj;
      p.local_k = p.gk;
    } else {
      p.owner_rank = 0;
      p.local_i = p.gi;
      p.local_j = p.gj;
      p.local_k = p.gk;
    }
  }
}

}  // namespace frehg2
