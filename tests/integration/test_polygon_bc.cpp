// P12.4 BLOCKING GATE: a synthetic channel with a polygon outflow BC, driven through the
// production Orchestrator (build polygon index at init, apply every step). The realized per-step
// outflow equals the prescribed discharge Q, so the accumulated outflow over the run matches
// Q * t_end to < 1e-6. A second case adds a matching polygon inflow source and asserts the
// injected and drained volumes both equal rate * t_end. These exercise the full main()->run()
// path (polygons are applied inside Orchestrator::step, every step, not just at init).
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

// SW-only flat channel (dx=dy=10, bottom 0, ponded init_eta). `polygon_block` appends the
// boundaries:/sources: YAML. No scalar `sources:` map is emitted (rainfall defaults to 0), so a
// polygon `sources:` sequence does not collide with a duplicate key.
std::string channelConfig(const std::string& out, int nx, int ny, double dt, double t_end,
                          long long max_steps, double init_eta, const std::string& polygon_block) {
  std::string s;
  s += "schema_version: '2.0'\n";
  s += "simulation: {id: channel, mode: surface_water}\n";
  s += "domain: {nx: " + std::to_string(nx) + ", ny: " + std::to_string(ny) +
       ", nz: 1, dx: 10.0, dy: 10.0, dz: 0.1, botz: 0.0, x0: 0.0, y0: 0.0}\n";
  s += "time: {dt: " + std::to_string(dt) + ", t_end: " + std::to_string(t_end) +
       ", max_steps: " + std::to_string(max_steps) + ", output_interval: 0}\n";
  s += "modules: {surface_water: true, groundwater: false, solute: false}\n";
  s += "surface_water: {gravity: 9.81, manning: 0.02, min_depth: 1.0e-8,\n";
  s += "  viscosity: {x: 1.0e-6, y: 1.0e-6}, h_diffusion_ref: 0.1, waterfall_depth: 1.0e-8}\n";
  s += "initial_conditions: {surface_water: {eta: " + std::to_string(init_eta) + "}}\n";
  s += "output: {format: hdf5, filename: " + out + ", io_mode: serial_gather}\n";
  s += polygon_block;
  return s;
}
}  // namespace

TEST_CASE("polygon BC gate: channel outflow matches prescribed discharge Q to < 1e-6") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;  // gate runs on a single rank

  const double dt = 1.0, t_end = 5.0, Q = 0.5;
  const long long nsteps = 5;
  // Outlet polygon over the last x-column (nx=20 => centroid x=195 in [189,201]).
  const std::string block =
      "boundaries:\n"
      "  - {name: outlet, type: bc_discharge,\n"
      "     vertices: [[189,-10],[201,-10],[201,50],[189,50]], rate: 0.5}\n";
  Config cfg = Config::fromString(
      channelConfig(std::string(kTmp) + "/poly_channel.h5", 20, 4, dt, t_end, nsteps, 2.0, block),
      "");

  Orchestrator orch;
  orch.initialize(cfg);
  REQUIRE(orch.hasPolygonBc());
  REQUIRE_FALSE(orch.hasPolygonSource());
  orch.run();

  REQUIRE(orch.stepCount() == nsteps);
  // The gate: accumulated outflow == Q * t_end (each step drains exactly Q*dt; water is ample).
  REQUIRE(orch.polygonOutflowVolume() == Approx(Q * t_end).margin(1e-6));
}

TEST_CASE("polygon source+BC gate: inflow injected and outflow drained both equal rate*t_end") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const double dt = 1.0, t_end = 5.0, rate = 0.4;
  const long long nsteps = 5;
  const std::string block =
      "boundaries:\n"
      "  - {name: outlet, type: bc_discharge,\n"
      "     vertices: [[189,-10],[201,-10],[201,50],[189,50]], rate: 0.4}\n"
      "sources:\n"
      "  - {name: inlet, type: inflow_rate,\n"
      "     vertices: [[-10,-10],[11,-10],[11,50],[-10,50]], rate: 0.4}\n";
  Config cfg = Config::fromString(
      channelConfig(std::string(kTmp) + "/poly_channel2.h5", 20, 4, dt, t_end, nsteps, 2.0, block),
      "");

  Orchestrator orch;
  orch.initialize(cfg);
  REQUIRE(orch.hasPolygonBc());
  REQUIRE(orch.hasPolygonSource());
  orch.run();

  REQUIRE(orch.polygonInflowVolume() == Approx(rate * t_end).margin(1e-6));
  REQUIRE(orch.polygonOutflowVolume() == Approx(rate * t_end).margin(1e-6));
}

TEST_CASE("polygon gate: no polygons => zero accumulators (b1-sw/b2-gw safe)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  Config cfg = Config::fromString(
      channelConfig(std::string(kTmp) + "/poly_none.h5", 8, 4, 1.0, 3.0, 3, 2.0, ""), "");
  Orchestrator orch;
  orch.initialize(cfg);
  REQUIRE_FALSE(orch.hasPolygonBc());
  REQUIRE_FALSE(orch.hasPolygonSource());
  orch.run();
  REQUIRE(orch.polygonOutflowVolume() == Approx(0.0).margin(1e-15));
  REQUIRE(orch.polygonInflowVolume() == Approx(0.0).margin(1e-15));
  REQUIRE(orch.polygonWellVolume() == Approx(0.0).margin(1e-15));
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
