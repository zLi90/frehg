// P7.2 acceptance: the grid comes entirely from YAML. Arbitrary surface sizes and vertical
// layer counts work through the unified driver, and an invalid grid (nx=0) produces a clear
// error rather than a crash.
#include <petscksp.h>

#include <Kokkos_Core.hpp>
#include <cmath>
#include <fstream>
#include <mpi.h>
#include <sstream>
#include <string>
#include <vector>

#include "core/Grid.hpp"
#include "driver/Orchestrator.hpp"
#include "frehg2_test.hpp"
#include "integration/orch_test_util.hpp"
#include "io/Config.hpp"
#include "re/ReSolver.hpp"
#include "swe/SweSolver.hpp"

using namespace frehg2;
using namespace frehg2::orch_test;

namespace {
const char* kTmp = FREHG2_IO_TMP;

// A bathymetry-from-file SWE config (nx x ny) reading `bath_file` with the given format.
std::string bathConfig(const std::string& out_h5, int nx, int ny, const std::string& bath_file,
                       const std::string& bath_fmt, double init_eta) {
  std::ostringstream s;
  s << "schema_version: '2.0'\n"
    << "simulation: {id: bath_test, mode: surface_water}\n"
    << "domain: {nx: " << nx << ", ny: " << ny
    << ", nz: 1, dx: 10.0, dy: 10.0, dz: 0.1, botz: 0.0,\n"
    << "  bathymetry: {from_file: true, file: " << bath_file << ", format: " << bath_fmt
    << "}}\n"
    << "time: {dt: 1.0, t_end: 3.0, max_steps: 3, output_interval: 0}\n"
    << "modules: {surface_water: true, groundwater: false, solute: false}\n"
    << "surface_water: {gravity: 9.81, manning: 0.02, min_depth: 1.0e-8,\n"
    << "  viscosity: {x: 1.0e-6, y: 1.0e-6}, h_diffusion_ref: 0.1, waterfall_depth: 1.0e-8}\n"
    << "initial_conditions: {surface_water: {eta: " << init_eta << "}}\n"
    << "sources: {surface: {rainfall: {from_file: false, rate: 0.0}, evaporation: {rate: 0.0}}}\n"
    << "output: {format: hdf5, filename: " << out_h5 << ", io_mode: serial_gather}\n";
  return s.str();
}

RealArr1DHost ownedDepth(const SweSolver& swe) {
  const Grid& g = swe.grid();
  const int nx = g.nx(), ny = g.ny();
  RealArr1DHost d("d", static_cast<size_t>(nx * ny));
  const auto& f = swe.fields();
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const int c = g.getSurfaceIndex(i, j);
      d(static_cast<size_t>(i + j * nx)) = std::max(f.eta(c) - f.bottom(c), 0.0);
    }
  return d;
}
}

TEST_CASE("arbitrary surface grid sizes run through the driver") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  for (int n : {5, 10, 50}) {
    const std::string out = std::string(kTmp) + "/grid_sw_" + std::to_string(n) + ".h5";
    Config cfg = Config::fromString(swConfig(out, n, n, 1.0, 2.0, 2, 0.5, 0.0), "");
    Orchestrator orch;
    orch.initialize(cfg);
    REQUIRE(orch.swe()->grid().nx() == n);
    REQUIRE(orch.swe()->grid().ny() == n);
    orch.run();
    REQUIRE(orch.stepCount() == 2);
  }
}

TEST_CASE("arbitrary vertical layer counts run through the driver") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  for (int nz : {1, 3, 10, 100}) {
    const std::string out = std::string(kTmp) + "/grid_gw_" + std::to_string(nz) + ".h5";
    // Adaptive GW dt -> use a large t_end so max_steps governs the step count.
    Config cfg = Config::fromString(gwConfig(out, nz, 1.0e-4, 1.0e9, 10, 0.1), "");
    Orchestrator orch;
    orch.initialize(cfg);
    REQUIRE(orch.re()->grid().nz() == nz);
    orch.run();
    REQUIRE(orch.stepCount() == 10);
  }
}

TEST_CASE("ESRI raster DEM bathymetry matches an equivalent plain-list run") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const int nx = 4, ny = 4;
  // Elevation in model orientation: E[gj][gi], gj increasing northward.
  std::vector<std::vector<double>> E(ny, std::vector<double>(nx));
  for (int gj = 0; gj < ny; ++gj)
    for (int gi = 0; gi < nx; ++gi) E[gj][gi] = 0.05 * (gi + 1) + 0.1 * gj;

  // Plain list: order gi + gj*nx.
  const std::string plain = std::string(kTmp) + "/bath_plain.txt";
  {
    std::ofstream o(plain);
    for (int gj = 0; gj < ny; ++gj)
      for (int gi = 0; gi < nx; ++gi) o << E[gj][gi] << "\n";
  }
  // ESRI raster: rows north-first (row 0 == gj = ny-1).
  const std::string ras = std::string(kTmp) + "/bath_dem.asc";
  {
    std::ofstream o(ras);
    o << "ncols " << nx << "\nnrows " << ny << "\nxllcorner 0.0\nyllcorner 0.0\n"
      << "cellsize 10.0\nNODATA_value -9999\n";
    for (int r = 0; r < ny; ++r) {
      const int gj = ny - 1 - r;
      for (int gi = 0; gi < nx; ++gi) o << E[gj][gi] << (gi + 1 < nx ? " " : "\n");
    }
  }

  RealArr1DHost from_plain, from_raster;
  {
    Config c = Config::fromString(
        bathConfig(std::string(kTmp) + "/bath_plain.h5", nx, ny, plain, "list", 0.9), "");
    Orchestrator orch;
    orch.initialize(c);
    orch.run();
    from_plain = ownedDepth(*orch.swe());
  }
  {
    Config c = Config::fromString(
        bathConfig(std::string(kTmp) + "/bath_raster.h5", nx, ny, ras, "raster", 0.9), "");
    Orchestrator orch;
    orch.initialize(c);
    orch.run();
    from_raster = ownedDepth(*orch.swe());
  }
  const double l2 = relL2(from_raster, from_plain);
  std::fprintf(stderr, "  raster-vs-plain bathymetry rel-L2 = %.3e\n", l2);
  REQUIRE(l2 < 1.0e-12);
}

TEST_CASE("raster DEM with NODATA notches runs (masking deferred, cells filled)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const int nx = 4, ny = 4;
  const std::string ras = std::string(kTmp) + "/bath_nodata.asc";
  {
    std::ofstream o(ras);
    o << "ncols " << nx << "\nnrows " << ny << "\nxllcorner 0.0\nyllcorner 0.0\n"
      << "cellsize 10.0\nNODATA_value -9999\n";
    for (int r = 0; r < ny; ++r) {
      for (int gi = 0; gi < nx; ++gi) {
        const bool notch = (r == 0 && gi == 0);  // a single NODATA corner
        o << (notch ? -9999.0 : 0.1 * (gi + 1)) << (gi + 1 < nx ? " " : "\n");
      }
    }
  }
  Config c = Config::fromString(
      bathConfig(std::string(kTmp) + "/bath_nodata.h5", nx, ny, ras, "raster", 0.9), "");
  Orchestrator orch;
  orch.initialize(c);
  orch.run();
  // Run completes and state stays finite (NODATA filled with min valid elevation).
  const RealArr1DHost d = ownedDepth(*orch.swe());
  for (size_t i = 0; i < d.extent(0); ++i) REQUIRE(std::isfinite(d(i)));
  REQUIRE(orch.stepCount() == 3);
}

TEST_CASE("invalid grid (nx=0) gives a clear error, not a crash") {
  bool threw = false;
  std::string msg;
  try {
    // nx=0 must be rejected by config validation before any solver is built.
    Config cfg = Config::fromString(swConfig(std::string(kTmp) + "/bad.h5", 0, 5, 1.0, 1.0,
                                             1, 0.5, 0.0),
                                    "");
    Orchestrator orch;
    orch.initialize(cfg);
  } catch (const std::runtime_error& e) {
    threw = true;
    msg = e.what();
  }
  REQUIRE(threw);
  REQUIRE(msg.find("domain.nx") != std::string::npos);
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
