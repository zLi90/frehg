// P8.3.2: explicit operator-split advection.
//   - At unit Courant number the conservative upwind scheme is an EXACT one-cell shift, so a
//     advected Gaussian pulse matches the analytically shifted pulse to machine precision.
//   - The CFL safety check refuses to advect (and leaves the field untouched) when the
//     Courant number exceeds cfl_max.
//   - The scheme conserves the cell-sum of concentration on a closed domain.
#define FREHG2_TEST_IMPL
#define FREHG2_TEST_USE_KOKKOS
#include "frehg2_test.hpp"

#include <cmath>

#include <Kokkos_Core.hpp>

#include "core/Grid.hpp"
#include "solute/Advection.hpp"
#include "solute/SoluteFlow.hpp"
#include "solute/SoluteParams.hpp"

using namespace frehg2;

namespace {

real gaussian(real x, real x0, real sd) {
  const real z = (x - x0) / sd;
  return std::exp(-0.5 * z * z);
}

real sumInterior(const RealArr1DHost& c, const Grid& g) {
  real s = 0.0;
  for (int j = 0; j < g.ny(); ++j)
    for (int i = 0; i < g.nx(); ++i)
      s += c(static_cast<size_t>(g.getSurfaceIndex(i, j)));
  return s;
}

}  // namespace

TEST_CASE("advection: unit-Courant x-advection is an exact one-cell shift") {
  const int nx = 24, ny = 1;
  Grid grid(nx, ny, 1, 1.0, 1.0, 1.0);
  SoluteParams p;  // upwind, cfl_max = 1
  SoluteFlow flow;
  flow.u = RealArr1DHost("u", grid.nSurfaceStorageCell());
  flow.v = RealArr1DHost("v", grid.nSurfaceStorageCell());
  flow.depth = RealArr1DHost("d", grid.nSurfaceStorageCell());
  RealArr1DHost c("c", grid.nSurfaceStorageCell());

  const real u0 = 1.0, dt = 1.0;  // Courant = u0*dt/dx = 1
  for (size_t i = 0; i < flow.u.extent(0); ++i) {
    flow.u(i) = u0;
    flow.v(i) = 0.0;
    flow.depth(i) = 1.0;
  }
  const real x0 = 6.0, sd = 1.6;
  for (int i = 0; i < nx; ++i) c(static_cast<size_t>(grid.getSurfaceIndex(i, 0))) = gaussian(i, x0, sd);

  const int nsteps = 6;
  for (int s = 0; s < nsteps; ++s) {
    const real cfl = advectSurface(c, flow, grid, p, dt);
    REQUIRE(cfl == Approx(1.0).margin(1e-14));
  }
  // After nsteps unit-Courant shifts, cell i holds the original value of cell i-nsteps
  // (inflow from the left boundary is zero).
  for (int i = 0; i < nx; ++i) {
    const real expected = (i - nsteps >= 0) ? gaussian(i - nsteps, x0, sd) : 0.0;
    REQUIRE(c(static_cast<size_t>(grid.getSurfaceIndex(i, 0))) == Approx(expected).margin(1e-12));
  }
}

TEST_CASE("advection: unit-Courant y-advection is an exact one-cell shift") {
  const int nx = 1, ny = 24;
  Grid grid(nx, ny, 1, 1.0, 1.0, 1.0);
  SoluteParams p;
  SoluteFlow flow;
  flow.u = RealArr1DHost("u", grid.nSurfaceStorageCell());
  flow.v = RealArr1DHost("v", grid.nSurfaceStorageCell());
  flow.depth = RealArr1DHost("d", grid.nSurfaceStorageCell());
  RealArr1DHost c("c", grid.nSurfaceStorageCell());
  const real v0 = 1.0, dt = 1.0;
  for (size_t i = 0; i < flow.v.extent(0); ++i) {
    flow.u(i) = 0.0;
    flow.v(i) = v0;
    flow.depth(i) = 1.0;
  }
  const real y0 = 7.0, sd = 1.5;
  for (int j = 0; j < ny; ++j) c(static_cast<size_t>(grid.getSurfaceIndex(0, j))) = gaussian(j, y0, sd);
  const int nsteps = 5;
  for (int s = 0; s < nsteps; ++s) advectSurface(c, flow, grid, p, dt);
  for (int j = 0; j < ny; ++j) {
    const real expected = (j - nsteps >= 0) ? gaussian(j - nsteps, y0, sd) : 0.0;
    REQUIRE(c(static_cast<size_t>(grid.getSurfaceIndex(0, j))) == Approx(expected).margin(1e-12));
  }
}

TEST_CASE("advection: CFL safety check refuses and leaves the field unchanged") {
  const int nx = 16, ny = 1;
  Grid grid(nx, ny, 1, 1.0, 1.0, 1.0);
  SoluteParams p;  // cfl_max = 1
  SoluteFlow flow;
  flow.u = RealArr1DHost("u", grid.nSurfaceStorageCell());
  flow.v = RealArr1DHost("v", grid.nSurfaceStorageCell());
  flow.depth = RealArr1DHost("d", grid.nSurfaceStorageCell());
  RealArr1DHost c("c", grid.nSurfaceStorageCell());
  for (size_t i = 0; i < flow.u.extent(0); ++i) {
    flow.u(i) = 1.0;
    flow.depth(i) = 1.0;
  }
  RealArr1DHost c0("c0", grid.nSurfaceStorageCell());
  for (int i = 0; i < nx; ++i) {
    const real val = gaussian(i, 8.0, 2.0);
    c(static_cast<size_t>(grid.getSurfaceIndex(i, 0))) = val;
    c0(static_cast<size_t>(grid.getSurfaceIndex(i, 0))) = val;
  }
  const real dt_big = 1.5;  // Courant = 1.5 > 1
  const real cfl = advectSurface(c, flow, grid, p, dt_big);
  REQUIRE(cfl > p.cfl_max);
  for (int i = 0; i < nx; ++i)
    REQUIRE(c(static_cast<size_t>(grid.getSurfaceIndex(i, 0))) ==
            Approx(c0(static_cast<size_t>(grid.getSurfaceIndex(i, 0)))).margin(0.0));
}

TEST_CASE("advection: closed-domain mass conservation (sub-unit Courant)") {
  const int nx = 30, ny = 1;
  Grid grid(nx, ny, 1, 1.0, 1.0, 1.0);
  SoluteParams p;
  SoluteFlow flow;
  flow.u = RealArr1DHost("u", grid.nSurfaceStorageCell());
  flow.v = RealArr1DHost("v", grid.nSurfaceStorageCell());
  flow.depth = RealArr1DHost("d", grid.nSurfaceStorageCell());
  RealArr1DHost c("c", grid.nSurfaceStorageCell());
  for (size_t i = 0; i < flow.u.extent(0); ++i) {
    flow.u(i) = 0.7;  // Courant = 0.7 < 1
    flow.depth(i) = 1.0;
  }
  for (int i = 0; i < nx; ++i)
    c(static_cast<size_t>(grid.getSurfaceIndex(i, 0))) = gaussian(i, 8.0, 2.0);
  const real m0 = sumInterior(c, grid);
  for (int s = 0; s < 10; ++s) advectSurface(c, flow, grid, p, 1.0);
  const real m1 = sumInterior(c, grid);
  REQUIRE(std::fabs(m1 - m0) < 1e-12 * std::fabs(m0));
}
