// P12.3.3 / P12.3.4 acceptance: PolygonBC and PolygonSource applied to real solver state.
// Discharge removes exactly Q*dt; depth overrides eta; critical removes the weir volume; inflow
// adds exactly Q*dt; rainfall adds r*dt over masked cells; extraction removes Q*dt from the
// deepest GW cell (clamped at theta_r). Mass bookkeeping is asserted to machine precision.
#define FREHG2_TEST_IMPL
#define FREHG2_TEST_USE_KOKKOS
#include <cmath>
#include <vector>

#include "bc/PolygonBC.hpp"
#include "bc/PolygonIndex.hpp"
#include "bc/PolygonSource.hpp"
#include "core/Grid.hpp"
#include "frehg2_test.hpp"
#include "re/ReSolver.hpp"
#include "swe/SweSolver.hpp"

using namespace frehg2;

namespace {
// Column polygon over global x-cell i=icol (dx=10, origin 0): covers x in [icol*10, (icol+1)*10].
Polygon columnPolygon(int icol, double dx) {
  const double x0 = icol * dx, x1 = (icol + 1) * dx;
  Polygon p;
  p.vertices = {{x0, -50.0}, {x1, -50.0}, {x1, 1.0e6}, {x0, 1.0e6}};
  return p;
}

void initSwe(SweSolver& swe, double init_eta) {
  SweParams p;
  p.gravity = 9.81;
  p.min_depth = 1.0e-8;
  p.dt = 1.0;
  swe.setParams(p);
  swe.setBathymetryConstant(0.0);
  swe.initializeState(init_eta);
}

double etaAt(const SweSolver& swe, int i, int j) {
  return swe.fields().eta(swe.grid().getSurfaceIndex(i, j));
}
}  // namespace

TEST_CASE("polygon BC: discharge removes exactly Q*dt (outflow == prescribed Q)") {
  const double dx = 10.0, dy = 10.0, dt = 2.0, Q = 1.0;
  Grid g(4, 4, 1, dx, dy, 0.1);
  SweSolver swe(g, nullptr);
  initSwe(swe, 2.0);

  BcRegion r;
  r.kind = BcKind::Discharge;
  r.rate = Q;
  r.polygon = columnPolygon(3, dx);  // last column, 4 cells
  PolygonIndex idx;
  idx.build({r.polygon}, g, nullptr, 0.0, 0.0);
  PolygonBC bc({r}, idx, nullptr);

  const double out = bc.applySurface(swe, dt);
  REQUIRE(out == Approx(Q * dt).margin(1e-9));  // outflow matches prescribed Q to < 1e-6
  // Each of the 4 tagged cells lost Q*dt/4 volume => depth drop (Q*dt/4)/(dx*dy).
  const double drop = (Q * dt / 4.0) / (dx * dy);
  REQUIRE(etaAt(swe, 3, 0) == Approx(2.0 - drop).margin(1e-12));
  REQUIRE(etaAt(swe, 3, 2) == Approx(2.0 - drop).margin(1e-12));
  REQUIRE(etaAt(swe, 0, 0) == Approx(2.0).margin(1e-12));  // untagged unchanged
}

TEST_CASE("polygon BC: depth overrides eta; outflow clamped to available water") {
  const double dx = 10.0, dt = 1.0;
  Grid g(4, 4, 1, dx, dx, 0.1);
  SweSolver swe(g, nullptr);
  initSwe(swe, 2.0);

  BcRegion r;
  r.kind = BcKind::Depth;
  r.rate = 0.5;  // prescribe 0.5 m depth
  r.polygon = columnPolygon(0, dx);
  PolygonIndex idx;
  idx.build({r.polygon}, g, nullptr, 0.0, 0.0);
  PolygonBC bc({r}, idx, nullptr);
  const double net = bc.applySurface(swe, dt);

  REQUIRE(etaAt(swe, 0, 1) == Approx(0.5).margin(1e-12));  // bottom 0 + h 0.5
  REQUIRE(net == Approx((2.0 - 0.5) * dx * dx * 4.0).margin(1e-9));  // volume drained

  // Discharge larger than available is clamped (cannot drive depth negative).
  SweSolver swe2(g, nullptr);
  initSwe(swe2, 0.01);  // only 0.01 m available
  BcRegion big;
  big.kind = BcKind::Discharge;
  big.rate = 1.0e6;
  big.polygon = columnPolygon(0, dx);
  PolygonIndex idx2;
  idx2.build({big.polygon}, g, nullptr, 0.0, 0.0);
  PolygonBC bc2({big}, idx2, nullptr);
  const double out = bc2.applySurface(swe2, dt);
  REQUIRE(out == Approx(0.01 * dx * dx * 4.0).margin(1e-9));  // only available water leaves
  REQUIRE(etaAt(swe2, 0, 0) == Approx(0.0).margin(1e-12));     // drained to bottom
}

TEST_CASE("polygon BC: critical-depth weir removes sqrt(g) h^{3/2} min(dx,dy) dt") {
  const double dx = 10.0, dt = 0.5, h = 1.0;
  Grid g(4, 4, 1, dx, dx, 0.1);
  SweSolver swe(g, nullptr);
  initSwe(swe, h);  // depth = 1 m (bottom 0)
  BcRegion r;
  r.kind = BcKind::Critical;
  r.polygon = columnPolygon(3, dx);
  PolygonIndex idx;
  idx.build({r.polygon}, g, nullptr, 0.0, 0.0);
  PolygonBC bc({r}, idx, nullptr);
  const double out = bc.applySurface(swe, dt);
  const double per_cell = std::sqrt(9.81) * std::pow(h, 1.5) * dx * dt;  // q*dt, width=min(dx,dy)
  REQUIRE(out == Approx(per_cell * 4.0).margin(1e-9));
}

TEST_CASE("polygon source: inflow adds Q*dt; rainfall adds r*dt over masked cells") {
  const double dx = 10.0, dt = 2.0, Q = 3.0;
  Grid g(4, 4, 1, dx, dx, 0.1);
  SweSolver swe(g, nullptr);
  initSwe(swe, 1.0);

  SourceRegion in;
  in.kind = SourceKind::InflowRate;
  in.rate = Q;
  in.polygon = columnPolygon(0, dx);
  PolygonIndex idx;
  idx.build({in.polygon}, g, nullptr, 0.0, 0.0);
  PolygonSource src({in}, idx, nullptr);
  const double added = src.applySurface(swe, dt);
  REQUIRE(added == Approx(Q * dt).margin(1e-9));
  REQUIRE(etaAt(swe, 0, 0) == Approx(1.0 + (Q * dt / 4.0) / (dx * dx)).margin(1e-12));
  REQUIRE(etaAt(swe, 2, 0) == Approx(1.0).margin(1e-12));  // outside the inflow column

  SweSolver swe2(g, nullptr);
  initSwe(swe2, 1.0);
  SourceRegion rain;
  rain.kind = SourceKind::RainfallRate;
  rain.rate = 1.0e-3;  // m/s
  rain.polygon = columnPolygon(1, dx);
  PolygonIndex idx2;
  idx2.build({rain.polygon}, g, nullptr, 0.0, 0.0);
  PolygonSource src2({rain}, idx2, nullptr);
  const double r_added = src2.applySurface(swe2, dt);
  REQUIRE(etaAt(swe2, 1, 0) == Approx(1.0 + 1.0e-3 * dt).margin(1e-12));
  REQUIRE(r_added == Approx(1.0e-3 * dt * dx * dx * 4.0).margin(1e-9));
}

TEST_CASE("polygon source: extraction well removes Q*dt from the deepest GW cell, clamped") {
  const double dx = 10.0, dz = 0.1, dt = 1.0, Q = 0.5;
  Grid g3(4, 4, 5, dx, dx, dz);
  ReSolver re(g3, nullptr);
  ReParams rp;
  rp.dx = dx;
  rp.dy = dx;
  rp.dz = dz;
  rp.soil.theta_s = 0.4;
  rp.soil.theta_r = 0.05;
  re.setParams(rp);
  re.initializeUniformColumn(0.3);  // wc = 0.3 everywhere

  SourceRegion well;
  well.kind = SourceKind::ExtractionWell;
  well.rate = Q;
  well.polygon = columnPolygon(2, dx);  // one column of cells (4 in y)
  PolygonIndex idx;
  idx.build({well.polygon}, g3, nullptr, 0.0, 0.0);
  PolygonSource src({well}, idx, nullptr);

  const double removed = src.applySubsurface(re, dt);
  const double cell_vol = dx * dx * dz;
  REQUIRE(removed == Approx(Q * dt).margin(1e-9));
  // deepest cell (k=4) of a tagged column dropped by (Q*dt/4)/cell_vol.
  const int cdeep = g3.getIndex(2, 0, 4);
  REQUIRE(re.fields().wc(cdeep) == Approx(0.3 - (Q * dt / 4.0) / cell_vol).margin(1e-12));
  // shallower cells untouched; cells outside the column untouched.
  REQUIRE(re.fields().wc(g3.getIndex(2, 0, 0)) == Approx(0.3).margin(1e-12));
  REQUIRE(re.fields().wc(g3.getIndex(0, 0, 4)) == Approx(0.3).margin(1e-12));

  // A huge extraction is clamped at theta_r (cannot pump the cell below residual).
  ReSolver re2(g3, nullptr);
  re2.setParams(rp);
  re2.initializeUniformColumn(0.3);
  SourceRegion big = well;
  big.rate = 1.0e9;
  PolygonIndex idx2;
  idx2.build({big.polygon}, g3, nullptr, 0.0, 0.0);
  PolygonSource src2({big}, idx2, nullptr);
  const double removed_big = src2.applySubsurface(re2, dt);
  REQUIRE(re2.fields().wc(cdeep) == Approx(0.05).margin(1e-12));  // clamped at theta_r
  const double capacity = (0.3 - 0.05) * cell_vol * 4.0;          // 4 columns in y
  REQUIRE(removed_big == Approx(capacity).margin(1e-9));
}
