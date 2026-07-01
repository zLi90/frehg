#include "bc/PolygonIndex.hpp"

#include <tuple>

#include "core/MpiComm.hpp"

namespace frehg2 {

void PolygonIndex::build(const std::vector<Polygon>& polys, const Grid& grid, const MpiComm* mc,
                         double x0, double y0) {
  grid_ = grid;
  nx_ = grid.nx();
  ny_ = grid.ny();
  n_polygons_ = static_cast<int>(polys.size());
  cell_to_polygon_.assign(static_cast<size_t>(grid.nSurfaceStorageCell()), -1);
  if (n_polygons_ == 0) return;

  const double dx = grid.dx();
  const double dy = grid.dy();
  for (int j = 0; j < ny_; ++j)
    for (int i = 0; i < nx_; ++i) {
      int gi = i;
      int gj = j;
      if (mc != nullptr) std::tie(gi, gj) = mc->localToGlobal(i, j);
      const double x = x0 + (static_cast<double>(gi) + 0.5) * dx;
      const double y = y0 + (static_cast<double>(gj) + 0.5) * dy;
      for (int p = 0; p < n_polygons_; ++p) {
        if (polys[static_cast<size_t>(p)].contains(x, y)) {
          cell_to_polygon_[static_cast<size_t>(grid.getSurfaceIndex(i, j))] = p;
          break;  // first match wins
        }
      }
    }
}

int PolygonIndex::lookup(int surf_idx) const {
  if (surf_idx < 0 || surf_idx >= static_cast<int>(cell_to_polygon_.size())) return -1;
  return cell_to_polygon_[static_cast<size_t>(surf_idx)];
}

int PolygonIndex::countInPolygon(int p) const {
  int n = 0;
  for (int j = 0; j < ny_; ++j)
    for (int i = 0; i < nx_; ++i)
      if (cell_to_polygon_[static_cast<size_t>(grid_.getSurfaceIndex(i, j))] == p) ++n;
  return n;
}

std::vector<int> PolygonIndex::globalColumnCounts(const MpiComm* mc) const {
  std::vector<int> counts(static_cast<size_t>(n_polygons_), 0);
  for (int p = 0; p < n_polygons_; ++p) counts[static_cast<size_t>(p)] = countInPolygon(p);
  if (mc != nullptr && mc->size() > 1 && n_polygons_ > 0) {
    MPI_Allreduce(MPI_IN_PLACE, counts.data(), n_polygons_, MPI_INT, MPI_SUM, mc->comm());
  }
  return counts;
}

}  // namespace frehg2
