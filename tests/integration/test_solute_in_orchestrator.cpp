// P16.3.1 acceptance: the standalone P8 solute solver is wired into Orchestrator::step() under
// `solute.enabled: true`. A surface-water run with constant rainfall concentration must:
//   - raise the surface concentration where rain fell (conc non-zero),
//   - write a monitor CSV with a "C" probe column at output_interval cadence, and
//   - conserve total solute mass: final == initial + rain input (to 1e-9 relative).
//
// To make the mass balance exact we use a flat, closed, uniformly-rained wet region. The SWE
// solver excludes rainfall on the global y+ row (evaprain), so that row is made a DRY wall via
// the bathymetry; the remaining wet cells all receive identical rain, stay at equal depth, and
// develop no flow -> the advection (which conserves cell-sum) moves nothing and Sum(C*depth) is
// conserved to machine precision.
#include <petscksp.h>

#include <Kokkos_Core.hpp>
#include <cmath>
#include <fstream>
#include <mpi.h>
#include <sstream>
#include <string>

#include "core/Grid.hpp"
#include "driver/Orchestrator.hpp"
#include "frehg2_test.hpp"
#include "io/Config.hpp"
#include "swe/SweFields.hpp"
#include "swe/SweSolver.hpp"

using namespace frehg2;

namespace {
const char* kTmp = FREHG2_IO_TMP;

// nx x ny domain, rows 0..ny-2 flat (botz 0), the global y+ row a high dry wall.
std::string writeWallBathymetry(int nx, int ny) {
  const std::string path = std::string(kTmp) + "/solute_bath.txt";
  std::ofstream f(path);
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) f << (j == ny - 1 ? 100.0 : 0.0) << "\n";
  return path;
}

std::string soluteConfig(const std::string& out_h5, const std::string& bath, int nx, int ny,
                         double dt, double t_end, long long max_steps, double rain,
                         double c_rain, double out_interval) {
  std::ostringstream s;
  s << "schema_version: '2.0'\n"
    << "simulation: {id: solute_test, mode: surface_water}\n"
    << "domain: {nx: " << nx << ", ny: " << ny << ", nz: 1, dx: 10.0, dy: 10.0, dz: 0.1,\n"
    << "  botz: 0.0, bathymetry: {from_file: true, file: " << bath << ", format: list}}\n"
    << "time: {dt: " << dt << ", t_end: " << t_end << ", max_steps: " << max_steps
    << ", output_interval: " << out_interval << "}\n"
    << "modules: {surface_water: true, groundwater: false, solute: true}\n"
    << "surface_water: {gravity: 9.81, manning: 0.02, min_depth: 1.0e-8,\n"
    << "  viscosity: {x: 0.0, y: 0.0}, h_diffusion_ref: 0.1, waterfall_depth: 1.0e-8}\n"
    << "initial_conditions: {surface_water: {eta: 0.6}}\n"
    << "sources: {surface: {rainfall: {from_file: false, rate: " << rain
    << "}, evaporation: {rate: 0.0}}}\n"
    << "solute: {enabled: true, c_rain: " << c_rain << ", k_decay: 0.0, D: 0.0,\n"
    << "  advection_scheme: upwind, diffusion_scheme: none, cfl_max: 1.0}\n"
    << "monitors: {probes: [{name: center, i: 0, j: 1, fields: [C, eta]}]}\n"
    << "output: {format: hdf5, filename: " << out_h5 << ", io_mode: serial_gather}\n";
  return s.str();
}

int countDataRows(const std::string& csv) {
  std::ifstream in(csv);
  std::string line;
  int rows = 0;
  bool header = true;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (header) {
      header = false;
      continue;
    }
    ++rows;
  }
  return rows;
}
}  // namespace

TEST_CASE("solute runs through the Orchestrator: conc, monitor CSV, mass balance") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  // nx=1 closed column (the b1-sw-style closed basin the legacy SWE BCs are tuned for); the
  // global y+ row is a dry wall so the wet cells are uniformly rained and stay quiescent.
  const int nx = 1, ny = 4;
  const std::string bath = writeWallBathymetry(nx, ny);
  const std::string out = std::string(kTmp) + "/solute_orch.h5";
  const double dt = 1.0, t_end = 100.0, rain = 1.0e-4, c_rain = 10.0;

  Orchestrator orch;
  Config cfg = Config::fromString(
      soluteConfig(out, bath, nx, ny, dt, t_end, 100, rain, c_rain, 10.0), "");
  orch.initialize(cfg);
  REQUIRE(orch.soluteEnabled());
  orch.run();

  // The uniformly-rained wet column stays quiescent (velocities ~0), so the upwind advection
  // moves nothing and Sum(C*depth) is conserved to machine precision.
  {
    const SweFields& f = orch.swe()->fields();
    const Grid& g = orch.swe()->grid();
    double vmax = 0.0;
    for (int j = 0; j < ny - 1; ++j)
      for (int i = 0; i < nx; ++i) {
        const int c = g.getSurfaceIndex(i, j);
        vmax = std::max(vmax, std::abs(static_cast<double>(f.vv(c))));
      }
    REQUIRE(vmax < 1.0e-9);
  }

  // (1) Concentration is non-zero where rain fell.
  REQUIRE(orch.maxSurfaceConc() > 1.0e-3);

  // (2) Monitor CSV exists with a C column and data rows at output cadence.
  const std::string csv = std::string(kTmp) + "/monitors/solute_test.csv";
  std::ifstream hin(csv);
  REQUIRE(hin.good());
  std::string header;
  std::getline(hin, header);
  REQUIRE(header.find("center.C") != std::string::npos);
  REQUIRE(countDataRows(csv) >= 10);

  // (3) Mass balance: initial (0) + rain input == final, to 1e-9 relative.
  const double initial = orch.soluteInitialMass();
  const double added = orch.soluteAddedByRain();
  const double final_mass = orch.soluteFinalMass();
  REQUIRE(initial == Approx(0.0).margin(1e-30));
  REQUIRE(added > 0.0);
  const double expected = initial + added;
  const double rel = std::abs(final_mass - expected) / expected;
  std::fprintf(stderr, "  solute mass: initial=%.6g added=%.6g final=%.6g rel_err=%.3e\n",
               initial, added, final_mass, rel);
  REQUIRE(rel < 1.0e-9);
}

TEST_CASE("solute.enabled:false is a full no-op") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const int nx = 4, ny = 4;
  const std::string bath = writeWallBathymetry(nx, ny);
  const std::string out = std::string(kTmp) + "/solute_off.h5";
  std::string yaml = soluteConfig(out, bath, nx, ny, 1.0, 10.0, 10, 1.0e-4, 10.0, 0.0);
  // Flip solute off.
  const std::string needle = "solute: {enabled: true";
  yaml.replace(yaml.find(needle), needle.size(), "solute: {enabled: false");

  Orchestrator orch;
  Config cfg = Config::fromString(yaml, "");
  orch.initialize(cfg);
  REQUIRE(!orch.soluteEnabled());
  orch.run();
  REQUIRE(orch.maxSurfaceConc() == Approx(0.0).margin(1e-30));
  REQUIRE(orch.soluteFinalMass() == Approx(0.0).margin(1e-30));
}

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  MPI_Init(&argc, &argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  int rc = frehg2test::runAll();
  int grc = rc;
  MPI_Allreduce(MPI_IN_PLACE, &grc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  PetscFinalize();
  MPI_Finalize();
  Kokkos::finalize();
  return grc;
}
