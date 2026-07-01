#include <petscksp.h>

#include <Kokkos_Core.hpp>
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

int countDataRows(const std::string& path) {
  std::ifstream in(path);
  if (!in) return -1;
  std::string line;
  int n = 0;
  while (std::getline(in, line)) ++n;
  return n - 1;  // minus header
}
}  // namespace

TEST_CASE("Orchestrator writes monitor CSV at output_interval") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const std::string out = std::string(kTmp) + "/mon_orch.h5";
  std::ostringstream yaml;
  yaml << "schema_version: '2.0'\n"
       << "simulation: {id: mon_test, mode: surface_water}\n"
       << "domain: {nx: 5, ny: 5, nz: 1, dx: 10.0, dy: 10.0, dz: 0.1, botz: 0.0}\n"
       << "time: {dt: 1.0, t_end: 3.0, max_steps: 3, output_interval: 1.0}\n"
       << "modules: {surface_water: true, groundwater: false, solute: false}\n"
       << "surface_water: {gravity: 9.81, manning: 0.02, min_depth: 1.0e-8,\n"
       << "  viscosity: {x: 1.0e-6, y: 1.0e-6}, h_diffusion_ref: 0.1, waterfall_depth: 1.0e-8}\n"
       << "initial_conditions: {surface_water: {eta: 1.0}}\n"
       << "monitors:\n"
       << "  probes:\n"
       << "    - {name: c, xyz: [25.0, 25.0, 0.0], fields: [eta]}\n"
       << "sources: {surface: {rainfall: {from_file: false, rate: 0.0}, evaporation: {rate: 0.0}}}\n"
       << "output: {format: hdf5, filename: " << out << ", io_mode: serial_gather}\n";

  Config c = Config::fromString(yaml.str(), kTmp);
  Orchestrator orch;
  orch.initialize(c);
  orch.run();

  const std::string csv = std::string(kTmp) + "/monitors/mon_test.csv";
  REQUIRE(countDataRows(csv) == 4);  // t=0,1,2,3
}

TEST_CASE("legacy monitoring.points path produces CSV (b1-style indices)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const std::string out = std::string(kTmp) + "/mon_legacy.h5";
  std::ostringstream yaml;
  yaml << "schema_version: '2.0'\n"
       << "simulation: {id: b1mon, mode: surface_water}\n"
       << "domain: {nx: 1, ny: 10, nz: 1, dx: 80.0, dy: 80.0, dz: 0.1, botz: -3.0}\n"
       << "time: {dt: 1.0, t_end: 2.0, max_steps: 2, output_interval: 1.0}\n"
       << "modules: {surface_water: true, groundwater: false, solute: false}\n"
       << "surface_water: {gravity: 9.81, manning: 0.02, min_depth: 1.0e-8,\n"
       << "  viscosity: {x: 1.0e-6, y: 1.0e-6}, h_diffusion_ref: 0.1, waterfall_depth: 1.0e-8}\n"
       << "initial_conditions: {surface_water: {eta: -2.0}}\n"
       << "monitoring:\n"
       << "  points:\n"
       << "    - {i: 0, j: 4, fields: [depth]}\n"
       << "sources: {surface: {rainfall: {from_file: false, rate: 0.0}, evaporation: {rate: 0.0}}}\n"
       << "output: {format: hdf5, filename: " << out << ", io_mode: serial_gather}\n";

  Config c = Config::fromString(yaml.str(), kTmp);
  Orchestrator orch;
  orch.initialize(c);
  orch.run();

  const std::string csv = std::string(kTmp) + "/monitors/b1mon.csv";
  std::ifstream in(csv);
  REQUIRE(in.good());
  std::string header;
  REQUIRE(std::getline(in, header));
  REQUIRE(header == "time,monitor_0.depth");
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
