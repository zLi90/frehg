// Offline cell -> polygon map for polygon BC / source regions (P12.3.2).
//
// Built once at startup, before the first step. For each of this rank's OWNED surface columns
// (i, j) in [0, localNx) x [0, localNy), the column centroid is mapped to global coordinates
//   x = x0 + (gi + 0.5) dx,  y = y0 + (gj + 0.5) dy
// (gi, gj from MpiComm::localToGlobal; serial = identity), and tested against each polygon in
// order. The FIRST polygon whose ring contains the centroid wins (documented precedence for
// overlaps). Complexity is O(n_owned_cells * n_polygons), acceptable for n_polygons < 100
// (P12.4); a future P22 note may add an R-tree.
//
// The map is keyed by the halo-padded surface storage index (Grid::getSurfaceIndex), so callers
// look up exactly the index they iterate with. Subsurface (well) sources reuse this column map
// and address the deepest cell of the matched column.
#ifndef FREHG2_BC_POLYGON_INDEX_HPP
#define FREHG2_BC_POLYGON_INDEX_HPP

#include <vector>

#include "bc/Polygon.hpp"
#include "core/Grid.hpp"

namespace frehg2 {

class MpiComm;

class PolygonIndex {
 public:
  // Build the column -> polygon map. `grid` is this rank's LOCAL grid (owned interior + halo);
  // `mc` may be null for the serial/full-domain case. `x0`/`y0` are the global domain origin
  // (centroid of cell (0,0) is at (x0 + 0.5 dx, y0 + 0.5 dy)).
  void build(const std::vector<Polygon>& polys, const Grid& grid, const MpiComm* mc, double x0,
             double y0);

  // Polygon index in [0, nPolygons) for the column at halo-padded surface index `surf_idx`,
  // or -1 if the column lies in no polygon (or `surf_idx` is a ghost/out-of-range cell).
  int lookup(int surf_idx) const;

  int nPolygons() const { return n_polygons_; }
  int nLocalNx() const { return nx_; }
  int nLocalNy() const { return ny_; }
  // Number of owned columns mapped to polygon `p` (this rank only).
  int countInPolygon(int p) const;
  // Global owned-column count per polygon (MPI_Allreduce over `mc`; local when mc is null/serial).
  // This is the denominator for evenly distributing a region's inflow/extraction across all its
  // cells, correct even when a polygon straddles multiple ranks.
  std::vector<int> globalColumnCounts(const MpiComm* mc) const;
  const Grid& grid() const { return grid_; }
  bool empty() const { return n_polygons_ == 0; }

 private:
  std::vector<int> cell_to_polygon_;  // length grid.nSurfaceStorageCell(); -1 if none
  int nx_ = 0;
  int ny_ = 0;
  int n_polygons_ = 0;
  Grid grid_;
};

}  // namespace frehg2

#endif  // FREHG2_BC_POLYGON_INDEX_HPP
