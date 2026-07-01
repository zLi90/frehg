// P8.3.4: source/sink terms.
//   - Constant rain at concentration c_rain drives the surface concentration to c_rain.
//   - Infiltration mixing produces the water-height-weighted average concentration.
//   - First-order decay follows C(t) = C0 * exp(-k*t).
#define FREHG2_TEST_IMPL
#define FREHG2_TEST_USE_KOKKOS
#include "frehg2_test.hpp"

#include <algorithm>
#include <cmath>

#include <Kokkos_Core.hpp>

#include "core/Grid.hpp"
#include "solute/SoluteParams.hpp"
#include "solute/SourceSink.hpp"

using namespace frehg2;

TEST_CASE("source/sink: constant rain drives surface conc to c_rain") {
  for (real c_rain : {5.0, 35.0}) {
    Grid grid(1, 1, 1, 1.0, 1.0, 1.0);
    SoluteParams p;
    p.c_rain = c_rain;
    RealArr1DHost conc("conc", grid.nSurfaceStorageCell());
    RealArr1DHost depth("depth", grid.nSurfaceStorageCell());
    const size_t c = static_cast<size_t>(grid.getSurfaceIndex(0, 0));
    conc(c) = 0.0;
    depth(c) = 0.1;
    const real rain = 1.0e-3, dt = 1.0;
    for (int s = 0; s < 4000; ++s) applyRainfall(conc, depth, grid, rain, p, dt);
    REQUIRE(conc(c) == Approx(c_rain).margin(1e-3));
    // Steady state is exactly fixed: applying rain at c_rain to water already at c_rain
    // leaves it unchanged.
    conc(c) = c_rain;
    applyRainfall(conc, depth, grid, rain, p, dt);
    REQUIRE(conc(c) == Approx(c_rain).margin(1e-12));
  }
}

TEST_CASE("source/sink: SW<->GW interface exchange conserves solute mass (both directions)") {
  // The interface exchange is applied right after the water exchange. We mimic that pairing: the
  // water move shifts a height `vex` from one side to the other (concentrations held), then
  // applyInterfaceExchange moves the dissolved mass. The pair must conserve Sum(C*water).
  // vex > 0 = infiltration (SW -> GW); vex < 0 = seepage (GW -> SW).
  struct Case {
    real cs, cg, ds0, hg0, vex;
  };
  const Case cases[] = {
      {10.0, 0.0, 0.60, 0.040, 0.010},   // infiltration, clean surface slug into clean soil
      {2.0, 7.0, 0.50, 0.100, 0.030},    // infiltration, both sides non-zero
      {3.0, 9.0, 0.40, 0.200, -0.050},   // seepage, soil slug into surface
      {5.0, 1.0, 0.30, 0.150, -0.020},   // seepage, both sides non-zero
  };
  for (const auto& tc : cases) {
    Grid grid(1, 1, 1, 1.0, 1.0, 1.0);
    RealArr1DHost cs("cs", grid.nSurfaceStorageCell());
    RealArr1DHost cg("cg", grid.nCell());
    RealArr1DHost depth("depth", grid.nSurfaceStorageCell());
    RealArr1DHost hsub("hsub", grid.nSurfaceStorageCell());
    RealArr1DHost vex("vex", grid.nSurfaceStorageCell());
    const size_t csurf = static_cast<size_t>(grid.getSurfaceIndex(0, 0));
    const size_t ctop = static_cast<size_t>(grid.getIndex(0, 0, 0));
    cs(csurf) = tc.cs;
    cg(ctop) = tc.cg;
    // Water move (concentrations held): surface loses vex, soil gains vex.
    const real ds1 = tc.ds0 - tc.vex;
    const real hg1 = tc.hg0 + tc.vex;
    depth(csurf) = ds1;
    hsub(csurf) = hg1;
    vex(csurf) = tc.vex;
    const double m_before = static_cast<double>(tc.cs) * tc.ds0 + static_cast<double>(tc.cg) * tc.hg0;
    applyInterfaceExchange(cs, cg, depth, hsub, vex, grid);
    const double m_after =
        static_cast<double>(cs(csurf)) * ds1 + static_cast<double>(cg(ctop)) * hg1;
    REQUIRE(m_after == Approx(m_before).margin(1e-12));
    // Resulting concentrations stay bounded by the two endpoints (no spurious over/undershoot).
    const real lo = std::min(tc.cs, tc.cg);
    const real hi = std::max(tc.cs, tc.cg);
    REQUIRE(cs(csurf) >= lo - 1e-12);
    REQUIRE(cs(csurf) <= hi + 1e-12);
    REQUIRE(cg(ctop) >= lo - 1e-12);
    REQUIRE(cg(ctop) <= hi + 1e-12);
  }
}

TEST_CASE("source/sink: pure decay follows exp(-k t)") {
  Grid grid(1, 1, 1, 1.0, 1.0, 1.0);
  RealArr1DHost conc("conc", grid.nSurfaceStorageCell());
  const size_t c = static_cast<size_t>(grid.getSurfaceIndex(0, 0));
  const real C0 = 7.0, k = 0.1, dt = 2.0;
  conc(c) = C0;
  for (int s = 1; s <= 20; ++s) {
    applyDecaySurface(conc, grid, k, dt);
    const real expected = C0 * std::exp(-k * dt * s);
    REQUIRE(conc(c) == Approx(expected).margin(1e-10));
  }
}
