// P7.5 acceptance: every successful run writes simulation_summary.txt containing the required
// keys. This test runs a short simulation through the Orchestrator, parses the summary, and
// verifies each required key is present (and a couple of values are sane).
#include <petscksp.h>

#include <Kokkos_Core.hpp>
#include <cmath>
#include <fstream>
#include <mpi.h>
#include <sstream>
#include <string>

#include "driver/Orchestrator.hpp"
#include "frehg2_test.hpp"
#include "integration/orch_test_util.hpp"
#include "io/Config.hpp"

using namespace frehg2;
using namespace frehg2::orch_test;

namespace {
const char* kTmp = FREHG2_IO_TMP;

std::string readFile(const std::string& p) {
  std::ifstream f(p);
  std::stringstream ss;
  ss << f.rdbuf();
  return ss.str();
}

// Extract the numeric value following `key ` in a "key value\n" summary file.
double summaryValue(const std::string& text, const std::string& key) {
  const std::string needle = key + " ";
  const size_t p = text.find("\n" + needle) != std::string::npos
                       ? text.find("\n" + needle) + 1
                       : text.find(needle);
  if (p == std::string::npos) return std::nan("");
  const size_t b = p + needle.size();
  return std::stod(text.substr(b));
}
}  // namespace

TEST_CASE("simulation_summary.txt is written with the required keys") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const std::string out = std::string(kTmp) + "/summary.h5";
  Config cfg = Config::fromString(swConfig(out, 5, 5, 1.0, 4.0, 4, 0.6, 1.0e-5), "");
  Orchestrator orch;
  orch.initialize(cfg);
  orch.run();

  const std::string text = readFile(orch.summaryPath());
  REQUIRE(!text.empty());
  for (const char* key : {"frehg2_version", "simulation_id", "modules", "mpi_ranks",
                          "kokkos_execution_space", "total_runtime_seconds",
                          "output_intervals_completed", "wall_clock_hours"}) {
    REQUIRE(text.find(key) != std::string::npos);
  }
  // Spot-check content fidelity.
  REQUIRE(text.find("simulation_id \"sw_test\"") != std::string::npos);
  REQUIRE(text.find("mpi_ranks 1") != std::string::npos);
  REQUIRE(text.find("surface_water=true") != std::string::npos);
  REQUIRE(text.find("groundwater=false") != std::string::npos);
  // P21 perf breakdown keys (always written).
  REQUIRE(text.find("perf_sw_ksp_seconds") != std::string::npos);
  REQUIRE(text.find("perf_sw_assembly_seconds") != std::string::npos);
  REQUIRE(text.find("perf_cells_touched") != std::string::npos);
  // Flow-only runs must NOT emit the solute mass-balance block.
  REQUIRE(text.find("solute_initial_mass") == std::string::npos);

  // P19 water mass-balance budget: keys present and the SW-only closed budget closes
  // (final == initial + rain_in - outflow, with no GW/outlet here) to better than 1%.
  for (const char* key : {"water_initial_volume", "water_final_volume", "water_rain_volume_in",
                          "water_polygon_outflow_volume", "coupling_exchange_volume"}) {
    REQUIRE(text.find(key) != std::string::npos);
  }
  REQUIRE(text.find("groundwater_present false") != std::string::npos);
  REQUIRE(text.find("surface_water_present true") != std::string::npos);
  const double v0 = summaryValue(text, "water_initial_volume");
  const double vf = summaryValue(text, "water_final_volume");
  const double rin = summaryValue(text, "water_rain_volume_in");
  const double vout = summaryValue(text, "water_polygon_outflow_volume");
  REQUIRE(v0 > 0.0);          // init eta=0.6 over a flat bed => standing water
  REQUIRE(rin > 0.0);         // rainfall was applied
  const double expected = v0 + rin - vout;
  REQUIRE(std::abs(vf - expected) <= 1.0e-2 * expected);
}

namespace {
std::string soluteSummaryConfig(const std::string& out_h5) {
  std::ostringstream s;
  s << "schema_version: '2.0'\n"
    << "simulation: {id: sol_sum, mode: surface_water}\n"
    << "domain: {nx: 4, ny: 4, nz: 1, dx: 10.0, dy: 10.0, dz: 0.1, botz: 0.0}\n"
    << "time: {dt: 1.0, t_end: 10.0, max_steps: 10, output_interval: 0}\n"
    << "modules: {surface_water: true, groundwater: false, solute: true}\n"
    << "surface_water: {gravity: 9.81, manning: 0.02, min_depth: 1.0e-8,\n"
    << "  viscosity: {x: 1.0e-6, y: 1.0e-6}, h_diffusion_ref: 0.1, waterfall_depth: 1.0e-8}\n"
    << "initial_conditions: {surface_water: {eta: 0.6}}\n"
    << "sources: {surface: {rainfall: {from_file: false, rate: 1.0e-4}, evaporation: {rate: 0.0}}}\n"
    << "solute: {enabled: true, c_rain: 10.0, k_decay: 0.0, D: 0.0,\n"
    << "  advection_scheme: upwind, diffusion_scheme: none, cfl_max: 1.0}\n"
    << "output: {format: hdf5, filename: " << out_h5 << ", io_mode: serial_gather}\n";
  return s.str();
}
}  // namespace

TEST_CASE("simulation_summary.txt includes the solute mass-balance block when solute is on") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const std::string out = std::string(kTmp) + "/sol_summary.h5";
  Config cfg = Config::fromString(soluteSummaryConfig(out), "");
  Orchestrator orch;
  orch.initialize(cfg);
  orch.run();

  const std::string text = readFile(orch.summaryPath());
  for (const char* key : {"solute_initial_mass", "solute_final_mass", "solute_added_by_rain",
                          "solute_decay_loss", "solute_relative_mass_error"}) {
    REQUIRE(text.find(key) != std::string::npos);
  }
  REQUIRE(text.find("solute=true") != std::string::npos);
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
