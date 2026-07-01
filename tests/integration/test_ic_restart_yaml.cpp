// P14: restart IC type in YAML round-trips through Orchestrator::initialize + run().
#include <petscksp.h>

#include <Kokkos_Core.hpp>
#include <cmath>
#include <iomanip>
#include <mpi.h>
#include <sstream>
#include <string>

#include "driver/Orchestrator.hpp"
#include "frehg2_test.hpp"
#include "integration/orch_test_util.hpp"
#include "io/Config.hpp"
#include "swe/SweSolver.hpp"

using namespace frehg2;
using namespace frehg2::orch_test;

namespace {
const char* kTmp = FREHG2_IO_TMP;

std::string ckptKey(double t) {
  const long long r = std::llround(t);
  if (std::fabs(t - static_cast<double>(r)) < 1.0e-9) return std::to_string(r);
  std::ostringstream os;
  os << std::setprecision(12) << t;
  return os.str();
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
}  // namespace

TEST_CASE("YAML restart IC matches CLI restart path") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const double dt = 5.0;
  const double t_mid = 25.0;
  const double t_end = 50.0;
  const double init_eta = 0.6;

  const std::string outA = std::string(kTmp) + "/ic_yaml_ckpt_a.h5";
  {
    Config c = Config::fromString(swConfig(outA, 5, 5, dt, t_mid, 5, init_eta, 0.0, t_mid), kTmp);
    Orchestrator orch;
    orch.initialize(c);
    orch.run();
  }
  const std::string ckpt = std::string(kTmp) + "/ic_yaml_ckpt_a.ckpt." + ckptKey(t_mid) + ".h5";

  RealArr1DHost cli;
  {
    Config c = Config::fromString(swConfig(outA, 5, 5, dt, t_end, 10, init_eta), kTmp);
    Orchestrator orch;
    orch.initialize(c);
    orch.restart(ckpt, t_mid);
    cli = ownedDepth(*orch.swe());
    REQUIRE(orch.time() == Approx(t_end).margin(1e-9));
  }

  RealArr1DHost yaml_ic;
  {
    std::ostringstream yaml;
    yaml << "schema_version: '2.0'\n"
         << "simulation: {id: ic_restart, mode: surface_water}\n"
         << "domain: {nx: 5, ny: 5, nz: 1, dx: 10.0, dy: 10.0, dz: 0.1, botz: 0.0}\n"
         << "time: {dt: " << dt << ", t_end: " << t_end << ", max_steps: 10, output_interval: 0}\n"
         << "modules: {surface_water: true, groundwater: false, solute: false}\n"
         << "surface_water: {gravity: 9.81, manning: 0.02, min_depth: 1.0e-8,\n"
         << "  viscosity: {x: 1.0e-6, y: 1.0e-6}, h_diffusion_ref: 0.1, waterfall_depth: 1.0e-8}\n"
         << "initial_conditions: {type: restart, file: '" << ckpt << "', time: " << t_mid << "}\n"
         << "sources: {surface: {rainfall: {from_file: false, rate: 0.0}, evaporation: {rate: 0.0}}}\n"
         << "output: {format: hdf5, filename: " << outA << ", io_mode: serial_gather}\n";
    Config c = Config::fromString(yaml.str(), kTmp);
    Orchestrator orch;
    orch.initialize(c);
    orch.run();
    yaml_ic = ownedDepth(*orch.swe());
    REQUIRE(orch.time() == Approx(t_end).margin(1e-9));
  }

  const double l2 = relL2(yaml_ic, cli);
  std::fprintf(stderr, "  YAML restart IC vs CLI restart rel-L2 = %.3e\n", l2);
  REQUIRE(l2 < 1.0e-10);
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
