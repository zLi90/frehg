// P14 manufactured gate: constant / raster / formula / polygon IC dispatch on small grids.
#define FREHG2_TEST_IMPL
#define FREHG2_TEST_USE_KOKKOS
#include "frehg2_test.hpp"

#include <cmath>
#include <fstream>
#include <string>

#include "core/Grid.hpp"
#include "ic/ICApply.hpp"
#include "ic/ICConfig.hpp"
#include "io/Config.hpp"
#include "re/ReSolver.hpp"
#include "swe/SweSolver.hpp"

using namespace frehg2;

namespace {

constexpr const char* kTmp = FREHG2_IO_TMP;

RealArr1DHost surfaceEtaPhysical(const SweSolver& swe) {
  const Grid& g = swe.grid();
  const int nx = g.nx(), ny = g.ny();
  RealArr1DHost eta("eta", static_cast<size_t>(nx * ny));
  const auto& f = swe.fields();
  const real off = swe.params().offset;
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const int c = g.getSurfaceIndex(i, j);
      eta(static_cast<size_t>(i + j * nx)) = f.eta(c) - off;
    }
  return eta;
}

void writeListRaster(const std::string& path, int nx, int ny, double fill) {
  std::ofstream out(path);
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) out << (fill + 0.01 * static_cast<double>(i + j * nx)) << '\n';
}

ICApplyContext baseCtx(const Config& config, int gnx, int gny, int gnz) {
  ICApplyContext ctx;
  ctx.gnx = gnx;
  ctx.gny = gny;
  ctx.gnz = gnz;
  ctx.dx = config.getOr<double>("domain.dx", 1.0);
  ctx.dy = config.getOr<double>("domain.dy", 1.0);
  ctx.dz = config.getOr<double>("domain.dz", 1.0);
  ctx.x0 = config.getOr<double>("domain.x0", 0.0);
  ctx.y0 = config.getOr<double>("domain.y0", 0.0);
  ctx.botz = config.getOr<double>("domain.botz", 0.0);
  ctx.config = &config;
  return ctx;
}

}  // namespace

TEST_CASE("constant IC matches scalar initializeState path") {
  const std::string yaml =
      "schema_version: '2.0'\n"
      "modules: {surface_water: true, groundwater: false, solute: false}\n"
      "time: {dt: 1.0, t_end: 1.0}\n"
      "domain: {nx: 4, ny: 3, nz: 1, dx: 10.0, dy: 10.0, dz: 0.1, botz: -1.0}\n"
      "initial_conditions: {surface_water: {eta: 0.25}}\n";
  Config config = Config::fromString(yaml, kTmp);
  InitialConditionsConfig ic = parseInitialConditions(config, 0.0);
  ICApplyContext ctx = baseCtx(config, 4, 3, 1);

  Grid g(4, 3, 1, 10.0, 10.0, 0.1);
  SweSolver swe(g, nullptr);
  SweParams p;
  p.offset = 1.0;
  swe.setParams(p);
  RealArr1DHost bed("bed", 12);
  for (size_t k = 0; k < 12; ++k) bed(k) = 0.0;
  swe.setBathymetry(bed);

  applyInitialConditions(ic, ctx, &swe, nullptr);
  const RealArr1DHost eta = surfaceEtaPhysical(swe);
  for (size_t k = 0; k < eta.extent(0); ++k)
    REQUIRE(eta(k) == Approx(0.25).margin(1.0e-12));
}

TEST_CASE("raster list IC loads per-cell eta") {
  const std::string raster = std::string(kTmp) + "/ic_eta.list";
  writeListRaster(raster, 3, 2, 0.5);
  const std::string yaml =
      "schema_version: '2.0'\n"
      "modules: {surface_water: true, groundwater: false, solute: false}\n"
      "time: {dt: 1.0, t_end: 1.0}\n"
      "domain: {nx: 3, ny: 2, nz: 1, dx: 5.0, dy: 5.0, dz: 0.1, botz: 0.0}\n"
      "initial_conditions:\n"
      "  surface_water:\n"
      "    eta: {type: raster, file: ic_eta.list, format: list}\n";
  Config config = Config::fromString(yaml, kTmp);
  InitialConditionsConfig ic = parseInitialConditions(config, 0.0);
  ICApplyContext ctx = baseCtx(config, 3, 2, 1);

  Grid g(3, 2, 1, 5.0, 5.0, 0.1);
  SweSolver swe(g, nullptr);
  SweParams p;
  swe.setParams(p);
  RealArr1DHost bed("bed", 6);
  for (size_t k = 0; k < 6; ++k) bed(k) = 0.0;
  swe.setBathymetry(bed);

  applyInitialConditions(ic, ctx, &swe, nullptr);
  const RealArr1DHost eta = surfaceEtaPhysical(swe);
  for (int j = 0; j < 2; ++j)
    for (int i = 0; i < 3; ++i) {
      const double expect = 0.5 + 0.01 * static_cast<double>(i + j * 3);
      REQUIRE(eta(static_cast<size_t>(i + j * 3)) == Approx(expect).margin(1.0e-12));
    }
}

TEST_CASE("formula IC evaluates x+y on cell centroids") {
  const std::string yaml =
      "schema_version: '2.0'\n"
      "modules: {surface_water: true, groundwater: false, solute: false}\n"
      "time: {dt: 1.0, t_end: 1.0}\n"
      "domain: {nx: 2, ny: 2, nz: 1, dx: 1.0, dy: 1.0, dz: 0.1, botz: 0.0, x0: 0.0, y0: 0.0}\n"
      "initial_conditions:\n"
      "  surface_water:\n"
      "    eta: {type: formula, formula: 'x+y'}\n";
  Config config = Config::fromString(yaml, kTmp);
  InitialConditionsConfig ic = parseInitialConditions(config, 0.0);
  ICApplyContext ctx = baseCtx(config, 2, 2, 1);

  Grid g(2, 2, 1, 1.0, 1.0, 0.1);
  SweSolver swe(g, nullptr);
  swe.setParams(SweParams{});
  RealArr1DHost bed("bed", 4);
  for (size_t k = 0; k < 4; ++k) bed(k) = 0.0;
  swe.setBathymetry(bed);

  applyInitialConditions(ic, ctx, &swe, nullptr);
  const RealArr1DHost eta = surfaceEtaPhysical(swe);
  // Centroids: (0.5,0.5)->1, (1.5,0.5)->2, (0.5,1.5)->2, (1.5,1.5)->3
  REQUIRE(eta(0) == Approx(1.0).margin(1.0e-12));
  REQUIRE(eta(1) == Approx(2.0).margin(1.0e-12));
  REQUIRE(eta(2) == Approx(2.0).margin(1.0e-12));
  REQUIRE(eta(3) == Approx(3.0).margin(1.0e-12));
}

TEST_CASE("polygon IC assigns region values on columns") {
  const std::string yaml =
      "schema_version: '2.0'\n"
      "modules: {surface_water: true, groundwater: false, solute: false}\n"
      "time: {dt: 1.0, t_end: 1.0}\n"
      "domain: {nx: 4, ny: 1, nz: 1, dx: 1.0, dy: 1.0, dz: 0.1, botz: 0.0, x0: 0.0, y0: 0.0}\n"
      "initial_conditions:\n"
      "  regions:\n"
      "    - {name: left, value: 1.0, vertices: [[0,0],[2,0],[2,1],[0,1]]}\n"
      "    - {name: right, value: 3.0, vertices: [[2,0],[4,0],[4,1],[2,1]]}\n"
      "  surface_water:\n"
      "    eta: {type: polygon, default: 0.0}\n";
  Config config = Config::fromString(yaml, kTmp);
  InitialConditionsConfig ic = parseInitialConditions(config, 0.0);
  ICApplyContext ctx = baseCtx(config, 4, 1, 1);

  Grid g(4, 1, 1, 1.0, 1.0, 0.1);
  SweSolver swe(g, nullptr);
  swe.setParams(SweParams{});
  RealArr1DHost bed("bed", 4);
  for (size_t k = 0; k < 4; ++k) bed(k) = 0.0;
  swe.setBathymetry(bed);

  applyInitialConditions(ic, ctx, &swe, nullptr);
  const RealArr1DHost eta = surfaceEtaPhysical(swe);
  REQUIRE(eta(0) == Approx(1.0).margin(1.0e-12));
  REQUIRE(eta(1) == Approx(1.0).margin(1.0e-12));
  REQUIRE(eta(2) == Approx(3.0).margin(1.0e-12));
  REQUIRE(eta(3) == Approx(3.0).margin(1.0e-12));
}

TEST_CASE("head raster IC sets GW pressure head field") {
  const std::string raster = std::string(kTmp) + "/ic_head.list";
  {
    std::ofstream out(raster);
    for (int k = 0; k < 3; ++k) out << (-1.0 - 0.1 * k) << '\n';
  }
  const std::string yaml =
      "schema_version: '2.0'\n"
      "modules: {surface_water: false, groundwater: true, solute: false}\n"
      "time: {dt: 1.0, t_end: 1.0}\n"
      "domain: {nx: 1, ny: 1, nz: 3, dx: 1.0, dy: 1.0, dz: 0.01, botz: -1.0}\n"
      "soil: {types: [{id: 0, theta_s: 0.4, theta_r: 0.1, vg: {alpha: 1.0, n: 2.0}, "
      "k_sat: {x: 0, y: 0, z: 1e-6}}]}\n"
      "initial_conditions:\n"
      "  groundwater:\n"
      "    h: {type: raster, file: ic_head.list, format: list}\n";
  Config config = Config::fromString(yaml, kTmp);
  InitialConditionsConfig ic = parseInitialConditions(config, 0.1);
  ICApplyContext ctx = baseCtx(config, 1, 1, 3);

  Grid g(1, 1, 3, 1.0, 1.0, 0.01);
  ReSolver re(g, nullptr);
  ReParams rp;
  rp.soil.theta_s = 0.4;
  rp.soil.theta_r = 0.1;
  rp.soil.alpha = 1.0;
  rp.soil.n = 2.0;
  rp.soil.Ks_z = 1.0e-6;
  re.setParams(rp);

  applyInitialConditions(ic, ctx, nullptr, &re);
  const Grid& rg = re.grid();
  for (int k = 0; k < 3; ++k) {
    const double expect = -1.0 - 0.1 * k;
    REQUIRE(re.fields().h(rg.getIndex(0, 0, k)) == Approx(expect).margin(1.0e-10));
  }
}
