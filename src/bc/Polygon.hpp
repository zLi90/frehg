// Polygon geometry for polygon-based BC / source-sink regions (P12.3.1).
//
// A Polygon is a closed 2D ring of (x, y) vertices in *global* domain coordinates. Membership
// is the even-odd ray-casting test (one O(n) pass per query), with a bounding-box fast reject.
// This is a brand-new Frehg2 feature — legacy `frehg` has rectangular BCs only, so there is no
// legacy code to port here.
//
// Header-only and free of Kokkos/PETSc/MPI: it is plain geometry consumed by PolygonIndex.
#ifndef FREHG2_BC_POLYGON_HPP
#define FREHG2_BC_POLYGON_HPP

#include <array>
#include <limits>
#include <string>
#include <vector>

namespace frehg2 {

struct Polygon {
  std::string name;
  std::vector<std::array<double, 2>> vertices;  // ring, in order; first != last (auto-closed)

  // Even-odd ray-casting point-in-polygon test. Casts a +x ray from (x, y) and counts edge
  // crossings; odd => inside. The ring is treated as closed (last vertex connects to first), so
  // concave rings (L-shapes) and explicit "donut" rings expressed as a single self-touching ring
  // are handled by the parity rule. A bounding-box check rejects far points in O(1).
  bool contains(double x, double y) const {
    const size_t n = vertices.size();
    if (n < 3) return false;
    double xmin = std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    double ymin = xmin;
    double ymax = xmax;
    for (const auto& v : vertices) {
      xmin = v[0] < xmin ? v[0] : xmin;
      xmax = v[0] > xmax ? v[0] : xmax;
      ymin = v[1] < ymin ? v[1] : ymin;
      ymax = v[1] > ymax ? v[1] : ymax;
    }
    if (x < xmin || x > xmax || y < ymin || y > ymax) return false;

    bool inside = false;
    for (size_t i = 0, j = n - 1; i < n; j = i++) {
      const double xi = vertices[i][0], yi = vertices[i][1];
      const double xj = vertices[j][0], yj = vertices[j][1];
      // Does the horizontal ray at height y cross edge (j -> i)?
      const bool straddles = (yi > y) != (yj > y);
      if (straddles) {
        const double x_cross = xj + (y - yj) * (xi - xj) / (yi - yj);
        if (x < x_cross) inside = !inside;
      }
    }
    return inside;
  }
};

}  // namespace frehg2

#endif  // FREHG2_BC_POLYGON_HPP
