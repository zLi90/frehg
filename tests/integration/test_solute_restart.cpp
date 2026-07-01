// P16.3.2 acceptance: restart preserves the solute `conc` field. Running a solute-enabled
// surface-water case continuously vs. running halfway, checkpointing, restarting, and finishing
// must reproduce the concentration field to L2 < 1e-12 (same executable / backend).
//
// A flat (no-wall) basin with rainfall is used on purpose: the rain-excluded y+ row creates a
// small depth gradient -> real flow -> the upwind advection actually redistributes solute, so
// the test exercises the full transport path (not a quiescent field).
#include <petscksp.h>

#include <Kokkos_Core.hpp>
#include <cmath>
#include <mpi.h>
#include <sstream>
#include <string>

#include "driver/Orchestrator.hpp"
#include "frehg2_test.hpp"
#include "io/Config.hpp"

using namespace frehg2;

namespace {
const char* kTmp = FREHG2_IO_TMP;

double relL2(const RealArr1DHost& a, const RealArr1DHost& b) {
  double num = 0.0, den = 0.0;
  for (size_t i = 0; i < a.extent(0); ++i) {
    const double d = a(i) - b(i);
    num += d * d;
    den += b(i) * b(i);
  }
  if (den == 0.0) return std::sqrt(num);
  return std::sqrt(num) / std::sqrt(den);
}

std::string soluteConfig(const std::string& out_h5, double dt, double t_end, long long max_steps,
                         double dt_checkpoint) {
  std::ostringstream s;
  s << "schema_version: '2.0'\n"
    << "simulation: {id: solute_restart, mode: surface_water}\n"
    << "domain: {nx: 5, ny: 5, nz: 1, dx: 10.0, dy: 10.0, dz: 0.1, botz: 0.0}\n"
    << "time: {dt: " << dt << ", t_end: " << t_end << ", max_steps: " << max_steps
    << ", output_interval: 0, dt_checkpoint: " << dt_checkpoint << ", max_checkpoints: 4}\n"
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

TEST_CASE("solute restart reproduces the conc field (L2 < 1e-12)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const double dt = 1.0;
  const double t_mid = 50.0;
  const double t_end = 100.0;

  // Continuous reference.
  RealArr1DHost cont;
  {
    Config c = Config::fromString(
        soluteConfig(std::string(kTmp) + "/sol_cont.h5", dt, t_end, 100, 0.0), "");
    Orchestrator orch;
    orch.initialize(c);
    orch.run();
    cont = orch.soluteSurfaceConcOwned();
  }
  REQUIRE(cont.extent(0) > 0);

  // Phase A: run to the midpoint and checkpoint there.
  const std::string outA = std::string(kTmp) + "/sol_phaseA.h5";
  {
    Config c = Config::fromString(soluteConfig(outA, dt, t_mid, 50, t_mid), "");
    Orchestrator orch;
    orch.initialize(c);
    orch.run();
  }
  const std::string ckpt = std::string(kTmp) + "/sol_phaseA.ckpt.50.h5";

  // Phase B: restart from the checkpoint and finish.
  RealArr1DHost restarted;
  {
    Config c = Config::fromString(soluteConfig(outA, dt, t_end, 100, 0.0), "");
    Orchestrator orch;
    orch.initialize(c);
    orch.restart(ckpt, t_mid);
    REQUIRE(orch.time() == Approx(t_end).margin(1e-9));
    restarted = orch.soluteSurfaceConcOwned();
  }

  const double l2 = relL2(restarted, cont);
  std::fprintf(stderr, "  solute restart-vs-continuous rel-L2 = %.3e\n", l2);
  REQUIRE(l2 < 1.0e-12);
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
