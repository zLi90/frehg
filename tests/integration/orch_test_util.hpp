// Shared helpers for the P7 Orchestrator integration tests: small in-memory YAML configs
// and L2 utilities. Kept header-only so each test links just the orchestrator library.
#pragma once

#include <cmath>
#include <sstream>
#include <string>

namespace frehg2 {
namespace orch_test {

inline double relL2(const RealArr1DHost& a, const RealArr1DHost& b) {
  double num = 0.0, den = 0.0;
  const size_t n = a.extent(0);
  for (size_t i = 0; i < n; ++i) {
    const double d = a(i) - b(i);
    num += d * d;
    den += b(i) * b(i);
  }
  if (den == 0.0) return std::sqrt(num);
  return std::sqrt(num) / std::sqrt(den);
}

// Surface-water-only config (constant flat bottom + ponded water + optional rainfall).
inline std::string swConfig(const std::string& out_h5, int nx, int ny, double dt,
                            double t_end, long long max_steps, double init_eta,
                            double rain_rate = 0.0, double dt_checkpoint = 0.0) {
  std::ostringstream s;
  s << "schema_version: '2.0'\n"
    << "simulation: {id: sw_test, mode: surface_water}\n"
    << "domain: {nx: " << nx << ", ny: " << ny
    << ", nz: 1, dx: 10.0, dy: 10.0, dz: 0.1, botz: 0.0}\n"
    << "time: {dt: " << dt << ", t_end: " << t_end << ", max_steps: " << max_steps
    << ", output_interval: 0, dt_checkpoint: " << dt_checkpoint << ", max_checkpoints: 4}\n"
    << "modules: {surface_water: true, groundwater: false, solute: false}\n"
    << "surface_water: {gravity: 9.81, manning: 0.02, min_depth: 1.0e-8,\n"
    << "  viscosity: {x: 1.0e-6, y: 1.0e-6}, h_diffusion_ref: 0.1, waterfall_depth: 1.0e-8}\n"
    << "initial_conditions: {surface_water: {eta: " << init_eta << "}}\n"
    << "sources: {surface: {rainfall: {from_file: false, rate: " << rain_rate
    << "}, evaporation: {rate: 0.0}}}\n"
    << "output: {format: hdf5, filename: " << out_h5 << ", io_mode: serial_gather}\n";
  return s.str();
}

// Groundwater-only config (b2-gw-like infiltration column).
inline std::string gwConfig(const std::string& out_h5, int nz, double dt, double t_end,
                            long long max_steps, double init_wc) {
  std::ostringstream s;
  s << "schema_version: '2.0'\n"
    << "simulation: {id: gw_test, mode: groundwater}\n"
    << "domain: {nx: 1, ny: 1, nz: " << nz
    << ", dx: 1.0, dy: 1.0, dz: 0.01, dz_incre: 1.0, botz: -1.0, follow_terrain: false}\n"
    << "time: {dt: " << dt << ", t_end: " << t_end << ", max_steps: " << max_steps
    << ", output_interval: 0}\n"
    << "modules: {surface_water: false, groundwater: true, solute: false}\n"
    << "groundwater: {solver: pca, full_3d: false, adaptive_dt: true, use_corrector: true,\n"
    << "  use_vg: true, use_mvg: false, air_entry_value: -0.02, dt_max: 2.0, dt_min: "
    << dt << ", co_max: 2.0, specific_storage: 1.0e-5, bc_type_gw: [0,0,0,0,0,1]}\n"
    << "soil: {types: [{id: 0, theta_s: 0.33, theta_r: 0.0,\n"
    << "  vg: {alpha: 1.43, n: 1.56}, k_sat: {x: 0.0, y: 0.0, z: 2.89e-6}}]}\n"
    << "initial_conditions: {groundwater: {wc: " << init_wc << "}}\n"
    << "output: {format: hdf5, filename: " << out_h5 << ", io_mode: serial_gather}\n";
  return s.str();
}

// Coupled SW+GW config.
inline std::string coupledConfig(const std::string& out_h5, int nx, int ny, int nz,
                                 double surface_dt, double t_end, long long max_steps,
                                 double init_eta, double init_wc,
                                 const std::string& coupling_mode = "sync") {
  std::ostringstream s;
  s << "schema_version: '2.0'\n"
    << "simulation: {id: coupled_test, mode: coupled}\n"
    << "domain: {nx: " << nx << ", ny: " << ny << ", nz: " << nz
    << ", dx: 10.0, dy: 10.0, dz: 0.05, dz_incre: 1.0, botz: -0.5, follow_terrain: false}\n"
    << "time: {dt: " << surface_dt << ", t_end: " << t_end << ", max_steps: " << max_steps
    << ", output_interval: 0}\n"
    << "modules: {surface_water: true, groundwater: true, solute: false}\n"
    << "coupling: {mode: " << coupling_mode << ", surface_dt: " << surface_dt
    << ", groundwater_dt: " << surface_dt << "}\n"
    << "surface_water: {gravity: 9.81, manning: 0.02, min_depth: 1.0e-8,\n"
    << "  viscosity: {x: 1.0e-6, y: 1.0e-6}, h_diffusion_ref: 0.1, waterfall_depth: 1.0e-8}\n"
    << "groundwater: {solver: pca, full_3d: false, adaptive_dt: true, use_corrector: true,\n"
    << "  use_vg: true, use_mvg: false, air_entry_value: -0.02, dt_max: 2.0, dt_min: 0.001,\n"
    << "  co_max: 2.0, specific_storage: 1.0e-5, bc_type_gw: [0,0,0,0,0,0]}\n"
    << "soil: {types: [{id: 0, theta_s: 0.4, theta_r: 0.08,\n"
    << "  vg: {alpha: 1.0, n: 2.0}, k_sat: {x: 0.0, y: 0.0, z: 1.0e-6}}]}\n"
    << "initial_conditions: {surface_water: {eta: " << init_eta << "},\n"
    << "  groundwater: {wc: " << init_wc << "}}\n"
    << "sources: {surface: {rainfall: {from_file: false, rate: 0.0}, evaporation: {rate: 0.0}}}\n"
    << "output: {format: hdf5, filename: " << out_h5 << ", io_mode: serial_gather}\n";
  return s.str();
}

// b1-sw config that reuses the repository's real bathymetry + rainfall files (config_dir is
// the benchmark directory). Lets the Orchestrator reproduce the P4 direct path exactly.
inline std::string b1Config(const std::string& out_h5, double init_eta, double dt,
                            double t_end, long long max_steps, double dt_checkpoint = 0.0) {
  std::ostringstream s;
  s << "schema_version: '2.0'\n"
    << "simulation: {id: b1-sw, mode: surface_water, legacy_finput: b1-input/}\n"
    << "domain: {nx: 1, ny: 10, nz: 1, dx: 80.0, dy: 80.0, dz: 0.1, botz: -3.0,\n"
    << "  bathymetry: {from_file: true, file: bath}}\n"
    << "time: {dt: " << dt << ", t_end: " << t_end << ", max_steps: " << max_steps
    << ", output_interval: 0, dt_checkpoint: " << dt_checkpoint << ", max_checkpoints: 4}\n"
    << "modules: {surface_water: true, groundwater: false, solute: false}\n"
    << "surface_water: {gravity: 9.81, manning: 0.019, min_depth: 1.0e-8,\n"
    << "  viscosity: {x: 1.0e-6, y: 1.0e-6}, h_diffusion_ref: 0.1, waterfall_depth: 1.0e-8}\n"
    << "initial_conditions: {surface_water: {eta: " << init_eta << "}}\n"
    << "sources: {surface: {rainfall: {from_file: true, file: rain}, evaporation: {rate: 0.0}}}\n"
    << "output: {format: hdf5, filename: " << out_h5 << ", io_mode: serial_gather}\n";
  return s.str();
}

}  // namespace orch_test
}  // namespace frehg2
