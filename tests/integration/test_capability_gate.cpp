// Capability gate: the production driver must FAIL FAST on any option that is parsed by the
// schema but not implemented by the realized solvers, instead of silently ignoring it. Each case
// enables one such option and asserts Orchestrator::initialize() throws with an identifying
// message. Positive controls prove the supported paths (incl. bottom fixed-flux and a valid
// output.variables list) still initialize.
#include <petscksp.h>

#include <Kokkos_Core.hpp>
#include <mpi.h>
#include <string>

#include "driver/Orchestrator.hpp"
#include "frehg2_test.hpp"
#include "io/Config.hpp"

using namespace frehg2;

namespace {
const char* kTmp = FREHG2_IO_TMP;

// Minimal SW-only config; `extra` is spliced into the surface_water map / top level as needed.
std::string swYaml(const std::string& sw_extra, const std::string& evap_model = "0",
                   const std::string& top_extra = "", const std::string& out_vars = "") {
  std::string s =
      "schema_version: '2.0'\n"
      "simulation: {id: gate, mode: surface_water}\n"
      "domain: {nx: 3, ny: 3, nz: 1, dx: 10.0, dy: 10.0, dz: 0.1, botz: 0.0}\n"
      "time: {dt: 1.0, t_end: 1.0, max_steps: 1, output_interval: 0}\n"
      "modules: {surface_water: true, groundwater: false, solute: false}\n"
      "surface_water: {gravity: 9.81, manning: 0.02, min_depth: 1.0e-8, "
      "viscosity: {x: 1.0e-6, y: 1.0e-6}, h_diffusion_ref: 0.1, waterfall_depth: 1.0e-8" +
      sw_extra +
      "}\n"
      "initial_conditions: {surface_water: {eta: 0.5}}\n"
      "sources: {surface: {rainfall: {from_file: false, rate: 0.0}, evaporation: {rate: 0.0, "
      "model: " +
      evap_model + "}}}\n" + top_extra +
      "output: {format: hdf5, filename: " + std::string(kTmp) + "/gate_sw.h5, io_mode: "
      "serial_gather" +
      (out_vars.empty() ? std::string("") : (", variables: " + out_vars)) + "}\n";
  return s;
}

// Minimal GW-only config with a configurable bc_type_gw vector and extra groundwater keys.
std::string gwYaml(const std::string& bc_type_gw, const std::string& gw_extra = "",
                   const std::string& solver = "pca") {
  return "schema_version: '2.0'\n"
         "simulation: {id: gate, mode: groundwater}\n"
         "domain: {nx: 1, ny: 1, nz: 6, dx: 1.0, dy: 1.0, dz: 0.01, dz_incre: 1.0, botz: -1.0}\n"
         "time: {dt: 0.1, t_end: 0.1, max_steps: 1, output_interval: 0}\n"
         "modules: {surface_water: false, groundwater: true, solute: false}\n"
         "groundwater: {solver: " + solver +
         ", full_3d: false, adaptive_dt: false, use_corrector: true, "
         "dt_max: 2.0, dt_min: 0.1, co_max: 2.0, specific_storage: 1.0e-5, bc_type_gw: " +
         bc_type_gw + gw_extra +
         "}\n"
         "soil: {types: [{id: 0, theta_s: 0.33, theta_r: 0.0, vg: {alpha: 1.43, n: 1.56}, "
         "k_sat: {x: 0.0, y: 0.0, z: 2.89e-6}}]}\n"
         "initial_conditions: {groundwater: {wc: 0.2}}\n"
         "output: {format: hdf5, filename: " +
         std::string(kTmp) + "/gate_gw.h5, io_mode: serial_gather}\n";
}

bool initThrows(const std::string& yaml, const std::string& expect_substr) {
  bool threw = false;
  try {
    Config cfg = Config::fromString(yaml, kTmp);
    Orchestrator orch;
    orch.initialize(cfg);
  } catch (const std::runtime_error& e) {
    threw = true;
    REQUIRE(std::string(e.what()).find(expect_substr) != std::string::npos);
  }
  return threw;
}

void initOk(const std::string& yaml) {
  Config cfg = Config::fromString(yaml, kTmp);
  Orchestrator orch;
  orch.initialize(cfg);  // must not throw
}
}  // namespace

TEST_CASE("Gate rejects unimplemented surface-water features") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  REQUIRE(initThrows(swYaml(", diffusive_wave: true"), "diffusive_wave"));
  REQUIRE(initThrows(swYaml(", wind: {enabled: true}"), "wind"));
  REQUIRE(initThrows(swYaml(", subgrid: {enabled: true}"), "subgrid"));
  REQUIRE(initThrows(swYaml("", /*evap_model=*/"1"), "evaporation.model"));
  // Baroclinic coupling (parsed under solute) is not implemented.
  REQUIRE(initThrows(swYaml("", "0", "solute: {baroclinic: true}\n"), "baroclinic"));
}

TEST_CASE("Gate rejects unimplemented groundwater BC modes and schemes") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  // Lateral Dirichlet / fixed-flux are not implemented (any non-zero on faces 0..3).
  REQUIRE(initThrows(gwYaml("[1,0,0,0,0,1]"), "bc_type_gw[0]"));
  REQUIRE(initThrows(gwYaml("[0,1,0,0,0,1]"), "bc_type_gw[1]"));
  REQUIRE(initThrows(gwYaml("[0,0,2,0,0,1]"), "bc_type_gw[2]"));
  REQUIRE(initThrows(gwYaml("[0,0,0,2,0,1]"), "bc_type_gw[3]"));
  // Out-of-range vertical mode.
  REQUIRE(initThrows(gwYaml("[0,0,0,0,3,1]"), "bc_type_gw[4]"));
  // Iterative (Picard/Newton) head solves are not implemented.
  REQUIRE(initThrows(gwYaml("[0,0,0,0,0,1]", ", iter_solve: 1"), "iter_solve"));
  REQUIRE(initThrows(gwYaml("[0,0,0,0,0,1]", "", /*solver=*/"newton"), "solver"));
}

TEST_CASE("Gate rejects invalid output.variables and accepts supported ones") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  // Unknown field name.
  REQUIRE(initThrows(swYaml("", "0", "", "[water_depth, bogus_field]"), "unknown field"));
  // Field of a disabled module (qx needs groundwater).
  REQUIRE(initThrows(swYaml("", "0", "", "[water_depth, qx]"), "groundwater"));
  // Supported surface set initializes.
  initOk(swYaml("", "0", "", "[water_depth, eta, u, v]"));
}

TEST_CASE("Gate accepts supported GW BC modes (no-flux/Dirichlet/bottom & top fixed-flux)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  initOk(gwYaml("[0,0,0,0,0,1]"));                       // top Dirichlet (b2-gw style)
  initOk(gwYaml("[0,0,0,0,1,1]"));                       // bottom Dirichlet
  initOk(gwYaml("[0,0,0,0,2,0]", ", qbot: 1.0e-7"));     // bottom fixed-flux (qbot)
  initOk(gwYaml("[0,0,0,0,0,2]", ", qtop: -1.0e-7"));    // top fixed-flux (qtop)
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
