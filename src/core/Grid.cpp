#include "core/Grid.hpp"

#include <stdexcept>
#include <string>

namespace frehg2 {

namespace {
void validate(int nx, int ny, int nz) {
  if (nx < 1 || ny < 1 || nz < 1) {
    throw std::runtime_error(
        "Grid: nx, ny, nz must all be >= 1 (got nx=" + std::to_string(nx) +
        ", ny=" + std::to_string(ny) + ", nz=" + std::to_string(nz) + ")");
  }
}
}  // namespace

Grid::Grid(int nx, int ny, int nz, real dx, real dy, real dz, real dz_incre)
    : nx_(nx), ny_(ny), nz_(nz), dx_(dx), dy_(dy), dz_(dz), dz_incre_(dz_incre) {
  validate(nx, ny, nz);
  if (dx <= 0.0 || dy <= 0.0 || dz <= 0.0) {
    throw std::runtime_error("Grid: dx, dy, dz must all be > 0");
  }
}

Grid::Grid(const DomainParams& p)
    : Grid(p.nx, p.ny, p.nz, p.dx, p.dy, p.dz, p.dz_incre) {}

}  // namespace frehg2
