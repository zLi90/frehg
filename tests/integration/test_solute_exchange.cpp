// P16-completion acceptance: the SW<->GW interface solute exchange is wired into the coupled
// production path. When the coupling moves water across the surface / top-GW-cell interface, the
// dissolved solute is carried with that exact volume (Orchestrator::applyCouplingSoluteExchange,
// co-located with the water exchange), so total solute mass is conserved to machine precision.
//
// Setup: a closed, quiescent column (nx=1, ny with a dry y+ wall, nz=1). A surface solute slug
// (initial_conditions.solute.surface) infiltrates into the clean groundwater as the coupling
// drives SW->GW exchange. No rain, no decay, no diffusion, no flow -> the only solute motion is
// the interface exchange, so Sum(C*depth*area) + Sum(C*wc*vol) must be invariant.
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

using namespace frehg2;

namespace {
const char* kTmp = FREHG2_IO_TMP;

// nx x ny domain: rows 0..ny-2 flat (botz 0), the global y+ row a high dry wall (matches the SWE
// rainfall/closed-basin convention used by the other solute tests).
std::string writeWallBathymetry(int nx, int ny) {
  const std::string path = std::string(kTmp) + "/solute_ex_bath.txt";
  std::ofstream f(path);
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) f << (j == ny - 1 ? 100.0 : 0.0) << "\n";
  return path;
}

std::string coupledSoluteConfig(const std::string& out_h5, const std::string& bath, int nx, int ny,
                                double dt, double t_end, long long max_steps, double init_eta,
                                double init_wc, double c_surface) {
  std::ostringstream s;
  s << "schema_version: '2.0'\n"
    << "simulation: {id: solute_exchange, mode: coupled}\n"
    << "domain: {nx: " << nx << ", ny: " << ny << ", nz: 1, dx: 10.0, dy: 10.0, dz: 0.05,\n"
    << "  dz_incre: 1.0, botz: 0.0, follow_terrain: false,\n"
    << "  bathymetry: {from_file: true, file: " << bath << ", format: list}}\n"
    << "time: {dt: " << dt << ", t_end: " << t_end << ", max_steps: " << max_steps
    << ", output_interval: 0}\n"
    << "modules: {surface_water: true, groundwater: true, solute: true}\n"
    << "coupling: {mode: sync, surface_dt: " << dt << ", groundwater_dt: " << dt << "}\n"
    << "surface_water: {gravity: 9.81, manning: 0.02, min_depth: 1.0e-8,\n"
    << "  viscosity: {x: 0.0, y: 0.0}, h_diffusion_ref: 0.1, waterfall_depth: 1.0e-8}\n"
    << "groundwater: {solver: pca, full_3d: false, adaptive_dt: false, use_corrector: true,\n"
    << "  use_vg: true, use_mvg: false, air_entry_value: -0.02, dt_max: " << dt
    << ", dt_min: " << dt << ", co_max: 2.0, specific_storage: 1.0e-5,\n"
    << "  bc_type_gw: [0, 0, 0, 0, 0, 0]}\n"
    << "soil: {types: [{id: 0, theta_s: 0.4, theta_r: 0.08,\n"
    << "  vg: {alpha: 1.0, n: 2.0}, k_sat: {x: 0.0, y: 0.0, z: 1.0e-6}}]}\n"
    << "initial_conditions: {surface_water: {eta: " << init_eta << "},\n"
    << "  groundwater: {wc: " << init_wc << "},\n"
    << "  solute: {surface: " << c_surface << ", subsurface: 0.0}}\n"
    << "sources: {surface: {rainfall: {from_file: false, rate: 0.0}, evaporation: {rate: 0.0}}}\n"
    << "solute: {enabled: true, c_rain: 0.0, k_decay: 0.0, D: 0.0,\n"
    << "  advection_scheme: upwind, diffusion_scheme: none, cfl_max: 1.0}\n"
    << "output: {format: hdf5, filename: " << out_h5 << ", io_mode: serial_gather}\n";
  return s.str();
}
}  // namespace

TEST_CASE("coupled SW<->GW solute exchange conserves total mass to machine precision") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;  // serial conservation gate (restart/IC paths are serial)

  const int nx = 1, ny = 4;
  const std::string bath = writeWallBathymetry(nx, ny);
  const std::string out = std::string(kTmp) + "/solute_exchange.h5";
  const double dt = 1.0, t_end = 200.0;
  const double init_eta = 0.6, init_wc = 0.2, c_surface = 10.0;

  Orchestrator orch;
  Config cfg = Config::fromString(
      coupledSoluteConfig(out, bath, nx, ny, dt, t_end, 200, init_eta, init_wc, c_surface), "");
  orch.initialize(cfg);
  REQUIRE(orch.soluteEnabled());
  REQUIRE(orch.coupled());

  // The solute IC must have been applied (parsed-but-unwired before this change): the surface
  // starts with the prescribed slug, the subsurface clean.
  REQUIRE(orch.maxSurfaceConc() == Approx(c_surface).margin(1e-12));
  REQUIRE(orch.maxSubsurfaceConc() == Approx(0.0).margin(1e-30));

  orch.run();

  // soluteInitialMass()/soluteFinalMass() are captured by the run loop (initial at run start, so
  // it reflects the applied IC); read them after run().
  const double initial = orch.soluteInitialMass();
  const double final_mass = orch.soluteFinalMass();
  REQUIRE(initial > 0.0);

  // Exchange actually happened: net infiltration moved water SW->GW (signed exchange < 0) and the
  // solute followed (subsurface conc now non-zero, surface slug partly drained).
  REQUIRE(orch.exchangeVolume() < 0.0);
  REQUIRE(orch.maxSubsurfaceConc() > 1.0e-3);

  // Conservation: with no rain / decay / diffusion and a quiescent field, the total solute mass
  // is invariant -> final == initial to machine precision.
  const double rel = std::abs(final_mass - initial) / initial;
  std::fprintf(stderr,
               "  coupled solute exchange: initial=%.10g final=%.10g exch_vol=%.6g rel_err=%.3e\n",
               initial, final_mass, static_cast<double>(orch.exchangeVolume()), rel);
  REQUIRE(rel < 1.0e-10);
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
