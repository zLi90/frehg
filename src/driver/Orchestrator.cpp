#include "driver/Orchestrator.hpp"

#include <Kokkos_Core.hpp>
#include <mpi.h>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "bc/Polygon.hpp"
#include "bc/PolygonBC.hpp"
#include "bc/PolygonConfig.hpp"
#include "bc/PolygonIndex.hpp"
#include "bc/PolygonSource.hpp"
#include "core/GwState.hpp"
#include "core/MpiComm.hpp"
#include "core/State.hpp"
#include "frehg2/core/ParallelFor.hpp"
#include "coupling/Coupling.hpp"
#include "driver/AsyncPipeline.hpp"
#include "ic/ICApply.hpp"
#include "ic/ICConfig.hpp"
#include "io/AsciiRaster.hpp"
#include "frehg2/io/build_info.hpp"
#include "frehg2/linear/LinearSolverFactory.hpp"
#include "frehg2/linear/SolverConfig.hpp"
#include "frehg2/perf/Timer.hpp"
#include "io/Config.hpp"
#include "io/Hdf5Reader.hpp"
#include "io/Sha256.hpp"
#include "linear/DomainDecomposition.hpp"
#include "monitoring/MonitorConfig.hpp"
#include "monitoring/MonitorWriter.hpp"
#include "re/ReSolver.hpp"
#include "re/VanGenuchten.hpp"
#include "soil/SoilConfig.hpp"
#include "soil/SoilMap.hpp"
#include "solute/SoluteFlow.hpp"
#include "solute/SoluteParams.hpp"
#include "solute/SoluteStepper.hpp"
#include "solute/SourceSink.hpp"
#include "swe/SweFields.hpp"
#include "swe/SweSolver.hpp"
#include "re/GwFields.hpp"

namespace frehg2 {

namespace {
// Legacy-exact piecewise-linear time-series interpolation (matches the b1-sw runner so the
// Orchestrator reproduces the P4 direct path bit-for-bit). `t`/`val` are paired samples.
double interpSeries(const std::vector<double>& t, const std::vector<double>& val, double tc) {
  const int n = static_cast<int>(t.size());
  if (n == 0) return 0.0;
  if (tc == 0.0) return val[0];
  int ind = 1;
  if (ind < n && tc > t[ind]) {
    ind += 1;
    while (ind < n && tc > t[ind]) ind += 1;
  }
  if (ind >= n) return val[static_cast<size_t>(n - 1)];
  return val[static_cast<size_t>(ind - 1)] +
         (val[static_cast<size_t>(ind)] - val[static_cast<size_t>(ind - 1)]) *
             (tc - t[static_cast<size_t>(ind - 1)]) /
             (t[static_cast<size_t>(ind)] - t[static_cast<size_t>(ind - 1)]);
}

std::vector<double> readDoublesFile(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("Orchestrator: cannot open data file '" + path + "'");
  std::vector<double> v;
  double x;
  while (in >> x) v.push_back(x);
  return v;
}

// P17 production schema gate: enforce the frozen v2 schema at the single production entry
// point (the Orchestrator). schema_version must be exactly "2.0" and every required top-level
// section must be present. Errors name the offending key and point at the migration tool.
void validateSchemaV2(const Config& config) {
  const std::string src = config.source();
  if (!config.has("schema_version")) {
    throw std::runtime_error(
        "schema error: required top-level 'schema_version' missing (in " + src +
        "). This is the frozen v2 schema; run tools/migrate_yaml_v1_to_v2 to upgrade.");
  }
  const std::string ver = config.getOr<std::string>("schema_version", "");
  if (ver != "2.0") {
    throw std::runtime_error("schema error: schema_version='" + ver +
                             "' is not supported; the production schema is '2.0' (in " + src +
                             "). Run tools/migrate_yaml_v1_to_v2 to upgrade.");
  }
  for (const char* sect : {"simulation", "domain", "time", "modules", "output"}) {
    if (!config.has(sect)) {
      throw std::runtime_error(std::string("schema error: required top-level section '") + sect +
                               "' missing (in " + src + ")");
    }
  }
}

// Capability gate: fail fast on any option that is PARSED by the schema but NOT IMPLEMENTED by
// the realized solvers. Without this, such options are silent no-ops (the config asks for physics
// the model never performs), which is exactly the class of "looks wired but isn't" bug we forbid.
// Every check below is keyed on a feature being ENABLED/selected; the default (disabled) paths are
// untouched, so all existing benchmarks are unaffected. Errors name the key and the reason.
void validateCapabilities(const Config& config, bool sw_enabled, bool gw_enabled,
                          bool solute_enabled) {
  const std::string src = config.source();
  auto reject = [&src](const std::string& key, const std::string& why) {
    throw std::runtime_error("unsupported configuration: '" + key + "' " + why + " (in " + src +
                             "). Disable it or remove the key; see docs/capabilities.md.");
  };

  // ---- Surface-water features that are declared in the schema but not implemented -------------
  if (sw_enabled) {
    if (config.getOr<bool>("surface_water.diffusive_wave", false))
      reject("surface_water.diffusive_wave",
             "is not implemented (only the dynamic-wave SWE, difuwave==0, is realized)");
    if (config.getOr<bool>("surface_water.wind.enabled", false))
      reject("surface_water.wind.enabled", "wind forcing is not implemented");
    if (config.getOr<bool>("surface_water.subgrid.enabled", false))
      reject("surface_water.subgrid.enabled", "subgrid bathymetry is not implemented");
    // Evaporation: only a constant/zero rate is applied; the legacy evaporation MODELS are not.
    const int evap_model = config.getOr<int>("sources.surface.evaporation.model", 0);
    if (evap_model != 0)
      reject("sources.surface.evaporation.model",
             "only constant-rate evaporation (model 0) is implemented; model=" +
                 std::to_string(evap_model) + " is not");
  }

  // ---- Density-driven (baroclinic) coupling: parsed but not implemented -----------------------
  if (config.getOr<bool>("solute.baroclinic", false))
    reject("solute.baroclinic", "baroclinic (density-driven) coupling is not implemented");

  // ---- Groundwater scheme + boundary conditions ----------------------------------------------
  if (gw_enabled) {
    const std::string gw_solver = config.getOr<std::string>("groundwater.solver", "pca");
    if (gw_solver != "pca")
      reject("groundwater.solver",
             "only the PCA predictor-corrector scheme ('pca') is implemented; '" + gw_solver +
                 "' (iterative Picard/Newton) is not");
    if (config.getOr<int>("groundwater.iter_solve", 0) != 0)
      reject("groundwater.iter_solve",
             "iterative (Picard/Newton) head solves are not implemented; use the PCA scheme "
             "(iter_solve: 0)");

    // GW boundary-condition modes. Realized: lateral faces (x-, x+, y-, y+) = no-flux (0) only;
    // bottom (z+) and top (z-) = no-flux (0) / Dirichlet (1) / fixed-flux (2). Lateral Dirichlet
    // and lateral fixed-flux are NOT implemented (the value would be silently ignored), so they
    // are rejected rather than mishandled. (legacy bctype_GW indices 0..5.)
    const YAML::Node bc = config.root()["groundwater"]["bc_type_gw"];
    if (bc && bc.IsSequence()) {
      const char* face[6] = {"x-", "x+", "y-", "y+", "z+ (bottom)", "z- (top)"};
      for (int i = 0; i < 6 && i < static_cast<int>(bc.size()); ++i) {
        const int mode = bc[static_cast<size_t>(i)].as<int>();
        if (i < 4 && mode != 0)
          reject(std::string("groundwater.bc_type_gw[") + std::to_string(i) + "] (" + face[i] + ")",
                 "lateral Dirichlet/fixed-flux GW boundaries are not implemented; only no-flux (0) "
                 "is supported on lateral faces");
        if (i >= 4 && (mode < 0 || mode > 2))
          reject(std::string("groundwater.bc_type_gw[") + std::to_string(i) + "] (" + face[i] + ")",
                 "must be 0 (no-flux), 1 (Dirichlet), or 2 (fixed-flux)");
      }
    }
  }

  (void)solute_enabled;
}
}  // namespace

Orchestrator::Orchestrator() = default;
Orchestrator::~Orchestrator() = default;

real Orchestrator::nextMultipleAfter(real interval, real t) const {
  if (interval <= 0.0) return std::numeric_limits<real>::max();
  const real k = std::floor(t / interval + 1.0e-9) + 1.0;
  return k * interval;
}

void Orchestrator::mkParentDir(const std::string& path) const {
  if (rank_ != 0 || path.empty()) return;
  const std::filesystem::path p(path);
  const std::filesystem::path dir = p.parent_path();
  if (!dir.empty()) std::filesystem::create_directories(dir);
}

void Orchestrator::initialize(const Config& config) {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &size_);

  // ---- P17: production schema v2 gate -------------------------------------------------
  // The production driver runs ONLY the frozen v2 schema. A missing/incorrect schema_version
  // or a missing required top-level section is a hard error (use tools/migrate_yaml_v1_to_v2
  // to upgrade a v1/experimental YAML). Component-level parsers (e.g. SoluteParams::fromConfig)
  // do not enforce this, so unit tests can still exercise partial configs.
  validateSchemaV2(config);

  sim_id_ = config.getOr<std::string>("simulation.id", "frehg2");
  mode_ = config.getOr<std::string>("simulation.mode", "coupled");

  sw_enabled_ = config.getOr<bool>("modules.surface_water", false);
  gw_enabled_ = config.getOr<bool>("modules.groundwater", false);
  solute_enabled_ = config.getOr<bool>("modules.solute", false);

  // Capability gate (run right after the schema gate): reject any enabled option the realized
  // solvers do not implement, so a config can never silently request physics that never runs.
  validateCapabilities(config, sw_enabled_, gw_enabled_, solute_enabled_);

  // simulation.mode selects the coupling pattern; warn on a mode/module mismatch.
  if (mode_ == "coupled" && (!sw_enabled_ || !gw_enabled_)) {
    if (rank_ == 0) {
      std::cerr << "frehg2: warning: simulation.mode='coupled' but "
                << (!sw_enabled_ ? "surface_water" : "groundwater")
                << " module is disabled; running uncoupled.\n";
    }
  }
  if (mode_ == "surface_water" && !sw_enabled_)
    throw std::runtime_error("Orchestrator: mode 'surface_water' requires modules.surface_water");
  if (mode_ == "groundwater" && !gw_enabled_)
    throw std::runtime_error("Orchestrator: mode 'groundwater' requires modules.groundwater");
  const std::string coupling_mode_str = config.getOr<std::string>("coupling.mode", "");
  coupled_ = sw_enabled_ && gw_enabled_ &&
             (mode_ == "coupled" || coupling_mode_str == "sync" ||
              coupling_mode_str == "sequential" || coupling_mode_str == "async");
  coupling_mode_ =
      (coupling_mode_str == "async") ? CouplingMode::Async : CouplingMode::Sequential;
  if (coupling_mode_ == CouplingMode::Async && coupled_ && size_ > 1 && rank_ == 0) {
    std::cerr << "frehg2: warning: coupling.mode='async' is validated on a single rank; running "
                 "the (numerically identical) synchronous coupling on "
              << size_
              << " ranks. Multi-rank async awaits the P10 PetscSubcomm split (Task 10.3.7).\n";
  }

  // ---- Grid (entirely from YAML; Config already validated nx,ny,nz > 0) -----------------
  nx_ = config.get<int>("domain.nx");
  ny_ = config.get<int>("domain.ny");
  nz_ = config.get<int>("domain.nz");
  dx_ = config.getOr<double>("domain.dx", 1.0);
  dy_ = config.getOr<double>("domain.dy", 1.0);
  dz_ = config.getOr<double>("domain.dz", 1.0);
  dz_incre_ = config.getOr<double>("domain.dz_incre", 1.0);
  botz_ = config.getOr<double>("domain.botz", 0.0);
  if (dx_ <= 0.0 || dy_ <= 0.0 || dz_ <= 0.0)
    throw std::runtime_error("Orchestrator: domain.dx/dy/dz must be positive");

  // ---- Time / coupling parameters ------------------------------------------------------
  dt_ = config.get<double>("time.dt");
  t_end_ = config.get<double>("time.t_end");
  max_steps_ = config.getOr<long long>("time.max_steps", 0);
  output_interval_ = config.getOr<double>("time.output_interval", 0.0);
  dt_checkpoint_ = config.getOr<double>("time.dt_checkpoint", 0.0);
  max_checkpoints_ = config.getOr<int>("time.max_checkpoints", 2);
  if (max_checkpoints_ < 1) max_checkpoints_ = 1;
  surface_dt_ = config.getOr<double>("coupling.surface_dt", dt_);
  groundwater_dt_ = config.getOr<double>("coupling.groundwater_dt", dt_);
  evap_rate_ = config.getOr<double>("sources.surface.evaporation.rate", 0.0);

  // ---- MPI process grid ----------------------------------------------------------------
  int mpi_nx = 1, mpi_ny = 1;
  if (size_ > 1) {
    mpi_nx = config.getOr<int>("domain.mpi.nx", 1);
    mpi_ny = config.getOr<int>("domain.mpi.ny", 1);
    if (mpi_nx * mpi_ny != size_) {
      // Default decomposition: split along y (matches the b1-sw NX=1 convention).
      mpi_nx = 1;
      mpi_ny = size_;
    }
  }
  mc_ = std::make_unique<MpiComm>(nx_, ny_, mpi_nx, mpi_ny, MPI_COMM_WORLD);

  buildGrids(config);
  loadRainfall(config);
  if (sw_enabled_) buildSurfaceWater(config);
  if (gw_enabled_) buildGroundwater(config);
  if (coupled_) {
    CouplingParams cp;
    cp.dx = dx_;
    cp.dy = dy_;
    cp.visc = 1.0;
    cp.min_depth = config.getOr<double>("surface_water.min_depth", 1.0e-8);
    coupling_ = std::make_unique<Coupling>(cp);
  }
  buildPolygons(config);
  buildSolute(config);

  config_sha256_ = sha256Hex(config.resolvedString());
  git_sha_ = build_info::kGitSha;
  setupInitialConditions(config);
  buildOutput(config);
  buildMonitors(config);

  // P19 water mass-balance: snapshot the post-IC water volume; rain inflow accumulates from 0.
  water_initial_volume_ = computeWaterVolume();
  water_rain_in_ = 0.0;

  // P21: attach the run-wide perf recorder to the flow solvers.
  if (swe_) swe_->setPerfRecorder(&perf_);
  if (re_) re_->setPerfRecorder(&perf_);
}

void Orchestrator::buildGrids(const Config&) {
  const int lnx = mc_->localNx();
  const int lny = mc_->localNy();
  swe_grid_ = Grid(lnx, lny, 1, dx_, dy_, dz_);
  gw_grid_ = Grid(lnx, lny, nz_, dx_, dy_, dz_, dz_incre_);
}

void Orchestrator::loadRainfall(const Config& config) {
  rain_from_file_ = config.getOr<bool>("sources.surface.rainfall.from_file", false);
  rain_rate_const_ = config.getOr<double>("sources.surface.rainfall.rate", 0.0);
  if (!rain_from_file_) return;
  const std::string file =
      config.resolvePath(config.getOr<std::string>("sources.surface.rainfall.file", "rain"));
  const std::vector<double> flat = readDoublesFile(file);
  if (flat.empty() || flat.size() % 2 != 0)
    throw std::runtime_error("Orchestrator: rainfall file '" + file +
                             "' must contain an even number of (time, rate) values");
  const size_t npairs = flat.size() / 2;
  rain_t_.resize(npairs);
  rain_v_.resize(npairs);
  for (size_t k = 0; k < npairs; ++k) {
    rain_t_[k] = flat[2 * k];
    rain_v_[k] = flat[2 * k + 1];
  }
}

real Orchestrator::rainAt(real t) const {
  if (!sw_enabled_) return 0.0;
  if (rain_from_file_) return interpSeries(rain_t_, rain_v_, t);
  return rain_rate_const_;
}

void Orchestrator::buildSurfaceWater(const Config& config) {
  swe_ = std::make_unique<SweSolver>(swe_grid_, mc_.get());

  SweParams p;
  p.gravity = config.getOr<double>("surface_water.gravity", 9.81);
  p.manning = config.getOr<double>("surface_water.manning", 0.0);
  p.min_depth = config.getOr<double>("surface_water.min_depth", 1.0e-8);
  p.visc_x = config.getOr<double>("surface_water.viscosity.x", 0.0);
  p.visc_y = config.getOr<double>("surface_water.viscosity.y", 0.0);
  p.hD = config.getOr<double>("surface_water.h_diffusion_ref", 0.1);
  p.wtfh = config.getOr<double>("surface_water.waterfall_depth", 1.0e-8);
  p.dt = dt_;

  // Global bathymetry (physical, row-major gi + gj*nx). From file or constant bottom.
  const int gnx = nx_;
  const int gny = ny_;
  std::vector<double> bath_global(static_cast<size_t>(gnx) * static_cast<size_t>(gny), botz_);
  const bool bath_from_file = config.getOr<bool>("domain.bathymetry.from_file", false);
  if (bath_from_file) {
    const std::string file =
        config.resolvePath(config.getOr<std::string>("domain.bathymetry.file", "bath"));
    // Format selection: ESRI ASCII raster DEM vs. a plain list of values. "auto" infers a
    // raster from a .asc extension; "raster"/"esri" force it; anything else is a plain list.
    std::string fmt = config.getOr<std::string>("domain.bathymetry.format", "auto");
    bool as_raster = (fmt == "raster" || fmt == "esri" || fmt == "ascii_raster");
    if (fmt == "auto")
      as_raster = file.size() >= 4 && file.compare(file.size() - 4, 4, ".asc") == 0;

    if (as_raster) {
      const AsciiRaster ras = AsciiRaster::fromFile(file);
      if (ras.ncols() != gnx || ras.nrows() != gny)
        throw std::runtime_error(
            "Orchestrator: raster DEM '" + file + "' is " + std::to_string(ras.ncols()) + "x" +
            std::to_string(ras.nrows()) + " but domain is nx=" + std::to_string(gnx) +
            ", ny=" + std::to_string(gny));
      const RealArr1DHost& d = ras.data();
      // Lowest valid elevation, used to fill NODATA cells (see deferral note below).
      double vmin = std::numeric_limits<double>::max();
      int n_nodata = 0;
      for (size_t i = 0; i < d.extent(0); ++i)
        if (!ras.isNoData(static_cast<int>(i), d(i))) vmin = std::min(vmin, double(d(i)));
      if (vmin == std::numeric_limits<double>::max()) vmin = botz_;
      // ESRI rows run north (row 0) to south; map so model gj increases northward.
      for (int gj = 0; gj < gny; ++gj)
        for (int gi = 0; gi < gnx; ++gi) {
          const size_t src = static_cast<size_t>(gi + (gny - 1 - gj) * gnx);
          double v = d(src);
          if (ras.isNoData(static_cast<int>(src), v)) {
            ++n_nodata;
            v = vmin;  // keep the cell active for now (masking deferred, see below)
          }
          bath_global[static_cast<size_t>(gi + gj * gnx)] = v;
        }
      // NOTE: true inactive-cell masking for NODATA notches (no-flux walls, assembly skip)
      // is deferred to the domain-masking / polygon-BC phase. The legacy-exact SWE/RE
      // solvers have no inactive-cell concept yet, so NODATA cells are filled with the
      // minimum valid elevation and remain active rather than being silently dropped.
      if (n_nodata > 0 && rank_ == 0)
        std::cerr << "frehg2: warning: raster DEM '" << file << "' has " << n_nodata
                  << " NODATA cell(s); inactive-cell masking is deferred to a later phase, "
                     "filling them with the minimum valid elevation (cells stay active).\n";
    } else {
      const std::vector<double> raw = readDoublesFile(file);
      if (raw.size() == static_cast<size_t>(gnx) * static_cast<size_t>(gny)) {
        for (size_t i = 0; i < raw.size(); ++i) bath_global[i] = raw[i];
      } else if (gnx == 1 && raw.size() == static_cast<size_t>(gny)) {
        for (int j = 0; j < gny; ++j)
          bath_global[static_cast<size_t>(j)] = raw[static_cast<size_t>(j)];
      } else {
        throw std::runtime_error("Orchestrator: bathymetry file '" + file + "' has " +
                                 std::to_string(raw.size()) + " values; expected " +
                                 std::to_string(gnx * gny) + " (nx*ny) or " +
                                 std::to_string(gny) + " (column, nx==1)");
      }
    }
  }
  // Legacy vertical datum: offset = -min(bottom) (MANDATORY for dry-cell singularity rows).
  double bmin = bath_global[0];
  for (double b : bath_global) bmin = std::min(bmin, b);
  p.offset = -bmin;
  swe_->setParams(p);

  // This rank's owned bathymetry block (physical, local i + j*localNx).
  const int lnx = mc_->localNx();
  const int lny = mc_->localNy();
  RealArr1DHost bed_local("bed_local", static_cast<size_t>(lnx) * static_cast<size_t>(lny));
  for (int lj = 0; lj < lny; ++lj)
    for (int li = 0; li < lnx; ++li) {
      const int gi = mc_->i0() + li;
      const int gj = mc_->j0() + lj;
      bed_local(static_cast<size_t>(li + lj * lnx)) =
          bath_global[static_cast<size_t>(gi + gj * gnx)] + p.offset;
    }
  swe_->setBathymetry(bed_local);

  dd2_ = std::make_unique<Decomp2D>(*mc_);
  SolverConfig cfg;
  cfg.ksp_type = config.getOr<std::string>("surface_water.solver.ksp_type", "cg");
  cfg.pc_type = config.getOr<std::string>("surface_water.solver.pc_type", "jacobi");
  cfg.rtol = config.getOr<double>("surface_water.solver.rtol", 1.0e-12);
  sw_linear_ = makeLinearSolver(config.getOr<std::string>("solver.backend", "petsc"), cfg);
  swe_->attachSolver(*sw_linear_, *dd2_);
}

void Orchestrator::buildGroundwater(const Config& config) {
  re_ = std::make_unique<ReSolver>(gw_grid_, mc_.get());

  ReParams rp;
  const YAML::Node soil_types = config.root()["soil"]["types"];
  if (!soil_types || !soil_types.IsSequence() || soil_types.size() == 0)
    throw std::runtime_error("Orchestrator: groundwater enabled but soil.types is empty");
  // Class 0 drives the uniform path (rp.soil); the full class list feeds the optional SoilMap.
  rp.soil = parseSoilClasses(config)[0];

  rp.dx = dx_;
  rp.dy = dy_;
  rp.dz = dz_;
  rp.botz = botz_;
  rp.dt = dt_;
  rp.dt_min = config.getOr<double>("groundwater.dt_min", 1.0e-4);
  rp.dt_max = config.getOr<double>("groundwater.dt_max", 2.0);
  rp.co_max = config.getOr<double>("groundwater.co_max", 2.0);
  rp.adaptive_dt = config.getOr<bool>("groundwater.adaptive_dt", true);
  rp.use_corrector = config.getOr<bool>("groundwater.use_corrector", true);
  rp.use_full3d = config.getOr<bool>("groundwater.full_3d", false);
  rp.follow_terrain = config.getOr<bool>("domain.follow_terrain", false);
  rp.htop = config.getOr<double>("groundwater.htop", 0.0);
  rp.hbot = config.getOr<double>("groundwater.hbot", 0.0);
  rp.qtop = config.getOr<double>("groundwater.qtop", 0.0);
  rp.qbot = config.getOr<double>("groundwater.qbot", 0.0);
  const YAML::Node bc = config.root()["groundwater"]["bc_type_gw"];
  if (bc && bc.IsSequence()) {
    for (int i = 0; i < 6 && i < static_cast<int>(bc.size()); ++i)
      rp.bc_type[i] = bc[static_cast<size_t>(i)].as<int>();
  }
  re_->setParams(rp);

  // Non-uniform soil (P13): attach a SoilMap ONLY when soil.map.from_file is set. With no class
  // map the RE solver keeps the uniform rp.soil path, so b2-gw stays bit-identical.
  buildSoilMap(config);

  // Partial-width top fixed-flux recharge (legacy bctype_GW[5]==2 / qtop, e.g. b3-kirkland's
  // polygontop). `groundwater.recharge` carries a `rate` [m/s] (negative = downward infiltration)
  // and a `polygon` ring [[x,y],...] in global coords; owned columns whose centroid lies inside
  // the ring get the rate, the rest get 0 (=> no-flux top). No-op when the section is absent, so
  // every existing GW config (b2-gw etc.) is unchanged.
  buildRecharge(config);

  dd3_ = std::make_unique<Decomp3D>(*mc_, nz_);
  SolverConfig cfg;
  cfg.ksp_type = config.getOr<std::string>("groundwater.solver.ksp_type", "cg");
  cfg.pc_type = config.getOr<std::string>("groundwater.solver.pc_type", "sor");
  cfg.rtol = config.getOr<double>("groundwater.solver.rtol", 1.0e-10);
  gw_linear_ = makeLinearSolver(config.getOr<std::string>("solver.backend", "petsc"), cfg);
  re_->attachSolver(*gw_linear_, *dd3_);
}

void Orchestrator::buildSoilMap(const Config& config) {
  // Only build a SoilMap when a per-cell class map is configured (soil.map.from_file). Otherwise
  // the RE solver keeps the uniform rp.soil path (bit-identical to P5 / b2-gw).
  if (!config.getOr<bool>("soil.map.from_file", false)) return;

  std::vector<SoilParams> classes = parseSoilClasses(config);
  const int n_class = static_cast<int>(classes.size());

  const int gnx = nx_;
  const int gny = ny_;

  // Load one 2D class-index layer (raster or list) into a global gnx*gny int vector
  // (row-major gi + gj*gnx). Shared by the per-column (2D) and per-cell (3D) paths.
  auto loadLayer = [&](const std::string& file) -> std::vector<int> {
    std::string fmt = config.getOr<std::string>("soil.map.format", "auto");
    bool as_raster = (fmt == "raster" || fmt == "esri" || fmt == "ascii_raster");
    if (fmt == "auto")
      as_raster = file.size() >= 4 && file.compare(file.size() - 4, 4, ".asc") == 0;
    auto toClassId = [&](double v) -> int {
      const int id = static_cast<int>(std::lround(v));
      if (id < 0 || id >= n_class)
        throw std::runtime_error("Orchestrator: soil class id " + std::to_string(id) + " in '" +
                                 file + "' is out of range [0," + std::to_string(n_class) +
                                 ") (have " + std::to_string(n_class) + " soil.types)");
      return id;
    };
    std::vector<int> g(static_cast<size_t>(gnx) * static_cast<size_t>(gny), 0);
    if (as_raster) {
      const AsciiRaster ras = AsciiRaster::fromFile(file);
      if (ras.ncols() != gnx || ras.nrows() != gny)
        throw std::runtime_error("Orchestrator: soil class raster '" + file + "' is " +
                                 std::to_string(ras.ncols()) + "x" + std::to_string(ras.nrows()) +
                                 " but domain is nx=" + std::to_string(gnx) +
                                 ", ny=" + std::to_string(gny));
      const RealArr1DHost& d = ras.data();
      // ESRI rows run north (row 0) to south; map so model gj increases northward (as bathymetry).
      for (int gj = 0; gj < gny; ++gj)
        for (int gi = 0; gi < gnx; ++gi) {
          const size_t src = static_cast<size_t>(gi + (gny - 1 - gj) * gnx);
          g[static_cast<size_t>(gi + gj * gnx)] = toClassId(d(src));
        }
    } else {
      const std::vector<double> raw = readDoublesFile(file);
      if (raw.size() != static_cast<size_t>(gnx) * static_cast<size_t>(gny))
        throw std::runtime_error("Orchestrator: soil class file '" + file + "' has " +
                                 std::to_string(raw.size()) + " values; expected nx*ny=" +
                                 std::to_string(gnx * gny));
      for (size_t i = 0; i < raw.size(); ++i) g[i] = toClassId(raw[i]);
    }
    return g;
  };

  const int lnx = mc_->localNx();
  const int lny = mc_->localNy();
  soil_map_ = std::make_unique<SoilMap>();
  soil_map_->setClasses(std::move(classes));

  // Fully heterogeneous per-cell (3D) soil (P23, SERGHEI-style): `soil.map.layers` is a
  // sequence of nz per-layer class files (top k=0 -> bottom k=nz-1). Otherwise the existing
  // per-column (2D) `soil.map.file` path is used (bit-identical to P13).
  const YAML::Node layers = config.root()["soil"]["map"]["layers"];
  if (layers && layers.IsSequence() && layers.size() > 0) {
    if (static_cast<int>(layers.size()) != nz_)
      throw std::runtime_error("Orchestrator: soil.map.layers has " +
                               std::to_string(layers.size()) + " layers but domain nz=" +
                               std::to_string(nz_) + " (need one class file per vertical layer)");
    std::vector<int> cls_local(static_cast<size_t>(lnx) * static_cast<size_t>(lny) *
                                   static_cast<size_t>(nz_),
                               0);
    for (int k = 0; k < nz_; ++k) {
      const std::string lf = config.resolvePath(layers[static_cast<size_t>(k)].as<std::string>());
      const std::vector<int> g = loadLayer(lf);
      for (int lj = 0; lj < lny; ++lj)
        for (int li = 0; li < lnx; ++li) {
          const int gi = mc_->i0() + li;
          const int gj = mc_->j0() + lj;
          cls_local[static_cast<size_t>(li + lj * lnx + k * lnx * lny)] =
              g[static_cast<size_t>(gi + gj * gnx)];
        }
    }
    soil_map_->setClassIndex3D(lnx, lny, nz_, std::move(cls_local));
  } else {
    const std::string file =
        config.resolvePath(config.getOr<std::string>("soil.map.file", "soil_class"));
    const std::vector<int> cls_global = loadLayer(file);
    std::vector<int> cls_local(static_cast<size_t>(lnx) * static_cast<size_t>(lny), 0);
    for (int lj = 0; lj < lny; ++lj)
      for (int li = 0; li < lnx; ++li) {
        const int gi = mc_->i0() + li;
        const int gj = mc_->j0() + lj;
        cls_local[static_cast<size_t>(li + lj * lnx)] =
            cls_global[static_cast<size_t>(gi + gj * gnx)];
      }
    soil_map_->setClassIndex(lnx, lny, std::move(cls_local));
  }
  re_->setSoilMap(soil_map_.get());
}

void Orchestrator::buildRecharge(const Config& config) {
  const YAML::Node rech = config.root()["groundwater"]["recharge"];
  if (!rech || !rech.IsMap()) return;
  const double rate = rech["rate"] ? rech["rate"].as<double>() : 0.0;
  const YAML::Node verts = rech["polygon"];

  Polygon poly;
  if (verts && verts.IsSequence()) {
    for (const auto& v : verts) {
      if (v.IsSequence() && v.size() >= 2)
        poly.vertices.push_back({v[0].as<double>(), v[1].as<double>()});
    }
  }

  const double x0 = config.getOr<double>("domain.x0", 0.0);
  const double y0 = config.getOr<double>("domain.y0", 0.0);
  const int lnx = mc_->localNx();
  const int lny = mc_->localNy();
  std::vector<real> q(static_cast<size_t>(lnx) * static_cast<size_t>(lny), 0.0);
  const bool whole_top = poly.vertices.size() < 3;  // no polygon => uniform recharge everywhere
  for (int lj = 0; lj < lny; ++lj)
    for (int li = 0; li < lnx; ++li) {
      const int gi = mc_->i0() + li;
      const int gj = mc_->j0() + lj;
      const double xc = x0 + (gi + 0.5) * dx_;
      const double yc = y0 + (gj + 0.5) * dy_;
      if (whole_top || poly.contains(xc, yc))
        q[static_cast<size_t>(li + lj * lnx)] = static_cast<real>(rate);
    }
  re_->setTopFluxField(q);

  if (rank_ == 0) {
    std::cerr << "frehg2: groundwater top recharge active (rate=" << rate
              << " m/s, " << (whole_top ? "whole top" : "polygon region") << ").\n";
  }
}

void Orchestrator::buildPolygons(const Config& config) {
  // Parse polygon BC / source regions (empty when the sections are absent => exact b1-sw/b2-gw
  // regression). The column index uses this rank's local grid + global origin (default 0,0).
  PolygonRegions pr = parsePolygonRegions(config);
  const double x0 = config.getOr<double>("domain.x0", 0.0);
  const double y0 = config.getOr<double>("domain.y0", 0.0);
  const Grid& grid = sw_enabled_ ? swe_grid_ : gw_grid_;

  std::vector<Polygon> bpolys;
  bpolys.reserve(pr.boundaries.size());
  for (const BcRegion& r : pr.boundaries) bpolys.push_back(r.polygon);
  PolygonIndex bidx;
  bidx.build(bpolys, grid, mc_.get(), x0, y0);
  poly_bc_ = std::make_unique<PolygonBC>(pr.boundaries, std::move(bidx), mc_.get());

  std::vector<Polygon> spolys;
  spolys.reserve(pr.sources.size());
  for (const SourceRegion& r : pr.sources) spolys.push_back(r.polygon);
  PolygonIndex sidx;
  sidx.build(spolys, grid, mc_.get(), x0, y0);
  poly_source_ = std::make_unique<PolygonSource>(pr.sources, std::move(sidx), mc_.get());

  if (rank_ == 0 && (!poly_bc_->empty() || !poly_source_->empty())) {
    std::cerr << "frehg2: polygon regions active: " << poly_bc_->nRegions()
              << " boundary, " << pr.sources.size() << " source.\n";
  }
}

void Orchestrator::setupInitialConditions(const Config& config) {
  const real default_wc = re_ ? re_->params().soil.theta_r : 0.0;
  const InitialConditionsConfig ic = parseInitialConditions(config, default_wc);

  if (ic.use_restart) {
    const std::string path = config.resolvePath(ic.restart_file);
    loadCheckpointState(path, ic.restart_time);
    return;
  }

  ICApplyContext ctx;
  ctx.gnx = nx_;
  ctx.gny = ny_;
  ctx.gnz = nz_;
  ctx.dx = dx_;
  ctx.dy = dy_;
  ctx.dz = dz_;
  ctx.x0 = config.getOr<double>("domain.x0", 0.0);
  ctx.y0 = config.getOr<double>("domain.y0", 0.0);
  ctx.botz = botz_;
  ctx.mc = mc_.get();
  ctx.config = &config;
  // Pass the canonical conc views when solute is enabled so the parsed solute ICs are applied
  // (P16-completion). Host aliases share storage with State/GwState (host build).
  RealArr1DHost surf_conc, subs_conc;
  if (solute_ && sw_enabled_) surf_conc = sw_state_->conc;
  if (solute_ && gw_enabled_) subs_conc = gw_state_->conc;
  applyInitialConditions(ic, ctx, swe_.get(), re_.get(), surf_conc, subs_conc);
}

void Orchestrator::loadCheckpointState(const std::string& checkpoint_file,
                                       real checkpoint_time) {
  if (size_ > 1)
    throw std::runtime_error(
        "Orchestrator checkpoint load is supported on a single rank (serial)");
  if (!swe_ && !re_)
    throw std::runtime_error("Orchestrator checkpoint load requires initialize() first");

  RestartState rs = Hdf5Reader::readCheckpoint(checkpoint_file, config_sha256_);

  auto restore = [&](std::vector<std::pair<std::string, RealArr1DHost*>> views,
                     const std::string& pfx) {
    for (auto& nv : views) {
      auto it = rs.extra.find(pfx + nv.first);
      if (it == rs.extra.end())
        throw std::runtime_error("Orchestrator checkpoint load: missing field '" + pfx +
                                 nv.first + "'");
      if (it->second.extent(0) != nv.second->extent(0))
        throw std::runtime_error("Orchestrator checkpoint load: field '" + pfx + nv.first +
                                 "' size mismatch (grid differs from checkpoint)");
      Kokkos::deep_copy(*nv.second, it->second);
    }
  };
  if (sw_enabled_) {
    if (!rs.has_sw)
      throw std::runtime_error("Orchestrator checkpoint load: checkpoint has no SW state");
    restore(swe_->fields().namedViews(), "sw_");
  }
  if (gw_enabled_) {
    if (!rs.has_gw)
      throw std::runtime_error("Orchestrator checkpoint load: checkpoint has no GW state");
    restore(re_->fields().namedViews(), "gw_");
    re_->primeAdaptiveDt(rs.dt);
  }
  if (solute_) {
    if (!rs.has_solute)
      throw std::runtime_error(
          "Orchestrator checkpoint load: solute enabled but checkpoint has no solute state");
    if (sw_enabled_) {
      auto it = rs.extra.find("solute_sw_conc");
      if (it == rs.extra.end())
        throw std::runtime_error("Orchestrator checkpoint load: missing 'solute_sw_conc'");
      Kokkos::deep_copy(sw_state_->conc, it->second);
    }
    if (gw_enabled_) {
      auto it = rs.extra.find("solute_gw_conc");
      if (it == rs.extra.end())
        throw std::runtime_error("Orchestrator checkpoint load: missing 'solute_gw_conc'");
      Kokkos::deep_copy(gw_state_->conc, it->second);
    }
  }

  t_ = checkpoint_time;
  gw_time_ = checkpoint_time;
  next_output_t_ = nextMultipleAfter(output_interval_, t_);
  next_checkpoint_t_ = nextMultipleAfter(dt_checkpoint_, t_);
  solute_added_by_rain_ = 0.0;
  solute_decay_loss_ = 0.0;
  resuming_from_checkpoint_ = true;
}

void Orchestrator::buildSolute(const Config& config) {
  // Fully no-op unless solute.enabled: no State/GwState allocation, no solute LinearSolver,
  // no PETSc assembly. solute_ stays null so applySolute() / computeSoluteMass() short-circuit.
  const SoluteParams sp = SoluteParams::fromConfig(config);
  // solute.enabled is authoritative (the plan's P16 key); it overrides modules.solute so a
  // run with solute.enabled:false does NO solute work even if modules.solute was left true.
  solute_enabled_ = sp.enabled;
  if (!sp.enabled) return;
  if (!sw_enabled_ && !gw_enabled_)
    throw std::runtime_error("Orchestrator: solute.enabled requires surface_water or groundwater");

  solute_substeps_ = std::max(1, config.getOr<int>("solute.substeps", 1));
  solute_start_time_ = config.getOr<double>("solute.start_time", 0.0);

  // The conc fields live in the canonical P2 State/GwState (Grid::getSurfaceIndex / getIndex).
  sw_state_ = std::make_unique<State>(swe_grid_);
  gw_state_ = std::make_unique<GwState>(gw_grid_);

  // The stepper drives whichever domains are populated by the per-step SoluteFlow.
  solute_ = std::make_unique<SoluteStepper>(gw_enabled_ ? gw_grid_ : swe_grid_, sp);

  // Implicit diffusion goes through dedicated LinearSolver backends (separate KSP/matrix from
  // the SWE/RE solves). Only built when implicit diffusion is actually requested (D > 0).
  if (sp.diffusion_scheme == "implicit" && sp.D > 0.0) {
    SolverConfig cfg;
    cfg.ksp_type = config.getOr<std::string>("solute.solver.ksp_type", "cg");
    cfg.pc_type = config.getOr<std::string>("solute.solver.pc_type", "jacobi");
    cfg.rtol = config.getOr<double>("solute.solver.rtol", 1.0e-12);
    const std::string backend = config.getOr<std::string>("solver.backend", "petsc");
    if (sw_enabled_ && dd2_) {
      solute_surf_linear_ = makeLinearSolver(backend, cfg);
      solute_->attachSurfaceDiffusion(*solute_surf_linear_, *dd2_);
    }
    if (gw_enabled_ && dd3_) {
      solute_subs_linear_ = makeLinearSolver(backend, cfg);
      solute_->attachSubsurfaceDiffusion(*solute_subs_linear_, *dd3_);
    }
  }
}

double Orchestrator::maxSurfaceConc() const {
  if (!solute_ || !sw_enabled_) return 0.0;
  const int lnx = mc_->localNx();
  const int lny = mc_->localNy();
  const Grid g = swe_grid_;
  RealArr1DHost conc = sw_state_->conc;
  double m = 0.0;
  Kokkos::parallel_reduce(
      "solute_maxconc_surf",
      Kokkos::MDRangePolicy<LoopExec, Kokkos::Rank<2>>({0, 0}, {lnx, lny}),
      KOKKOS_LAMBDA(int li, int lj, double& lmax) {
        const double v = conc(g.getSurfaceIndex(li, lj));
        if (v > lmax) lmax = v;
      },
      Kokkos::Max<double>(m));
  return m;
}

double Orchestrator::maxSubsurfaceConc() const {
  if (!solute_ || !gw_enabled_) return 0.0;
  const int lnx = mc_->localNx();
  const int lny = mc_->localNy();
  const Grid g = gw_grid_;
  RealArr1DHost conc = gw_state_->conc;
  double m = 0.0;
  Kokkos::parallel_reduce(
      "solute_maxconc_subs",
      Kokkos::MDRangePolicy<LoopExec, Kokkos::Rank<3>>({0, 0, 0}, {lnx, lny, nz_}),
      KOKKOS_LAMBDA(int li, int lj, int k, double& lmax) {
        const double v = conc(g.getIndex(li, lj, k));
        if (v > lmax) lmax = v;
      },
      Kokkos::Max<double>(m));
  return m;
}

RealArr1DHost Orchestrator::soluteSurfaceConcOwned() const {
  if (!solute_ || !sw_enabled_) return RealArr1DHost();
  const int lnx = mc_->localNx();
  const int lny = mc_->localNy();
  RealArr1DHost out("conc_surf", static_cast<size_t>(lnx) * static_cast<size_t>(lny));
  const Grid g = swe_grid_;
  RealArr1DHost conc = sw_state_->conc;
  parallelForSurface("solute_pack_surf", lnx, lny, KOKKOS_LAMBDA(int li, int lj) {
    out(static_cast<size_t>(li + lj * lnx)) = conc(g.getSurfaceIndex(li, lj));
  });
  return out;
}

RealArr1DHost Orchestrator::soluteSubsurfaceConcOwned() const {
  if (!solute_ || !gw_enabled_) return RealArr1DHost();
  const int lnx = mc_->localNx();
  const int lny = mc_->localNy();
  RealArr1DHost out("conc_subs", static_cast<size_t>(lnx) * static_cast<size_t>(lny) *
                                     static_cast<size_t>(nz_));
  const Grid g = gw_grid_;
  RealArr1DHost conc = gw_state_->conc;
  parallelForVolume("solute_pack_subs", lnx, lny, nz_, KOKKOS_LAMBDA(int li, int lj, int k) {
    out(static_cast<size_t>(li + lj * lnx + k * lnx * lny)) = conc(g.getIndex(li, lj, k));
  });
  return out;
}

double Orchestrator::computeSoluteMass() const {
  if (!solute_) return 0.0;
  const int lnx = mc_->localNx();
  const int lny = mc_->localNy();
  const double area = dx_ * dy_;
  double local = 0.0;
  if (sw_enabled_ && swe_) {
    const SweFields& f = swe_->fields();
    const Grid g = swe_grid_;
    RealArr1DHost conc = sw_state_->conc;
    RealArr1DHost eta = f.eta;
    RealArr1DHost bottom = f.bottom;
    double sw_sum = 0.0;
    Kokkos::parallel_reduce(
        "solute_mass_surf",
        Kokkos::MDRangePolicy<LoopExec, Kokkos::Rank<2>>({0, 0}, {lnx, lny}),
        KOKKOS_LAMBDA(int li, int lj, double& acc) {
          const int c = g.getSurfaceIndex(li, lj);
          const double depth = Kokkos::max(0.0, eta(c) - bottom(c));
          acc += conc(c) * depth * area;
        },
        sw_sum);
    local += sw_sum;
  }
  if (gw_enabled_ && re_) {
    const GwFields& f = re_->fields();
    const double vol = dx_ * dy_ * dz_;
    const Grid g = gw_grid_;
    RealArr1DHost conc = gw_state_->conc;
    RealArr1DHost wc = f.wc;
    double gw_sum = 0.0;
    Kokkos::parallel_reduce(
        "solute_mass_subs",
        Kokkos::MDRangePolicy<LoopExec, Kokkos::Rank<3>>({0, 0, 0}, {lnx, lny, nz_}),
        KOKKOS_LAMBDA(int li, int lj, int k, double& acc) {
          const int c = g.getIndex(li, lj, k);
          acc += conc(c) * wc(c) * vol;
        },
        gw_sum);
    local += gw_sum;
  }
  double global = local;
  if (size_ > 1) MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, mc_->comm());
  return global;
}

double Orchestrator::computeWaterVolume() const {
  const int lnx = mc_->localNx();
  const int lny = mc_->localNy();
  const double area = dx_ * dy_;
  double local = 0.0;
  if (sw_enabled_ && swe_) {
    const SweFields& f = swe_->fields();
    const Grid g = swe_grid_;
    RealArr1DHost eta = f.eta, bottom = f.bottom;
    double sw_sum = 0.0;
    Kokkos::parallel_reduce(
        "water_vol_surf",
        Kokkos::MDRangePolicy<LoopExec, Kokkos::Rank<2>>({0, 0}, {lnx, lny}),
        KOKKOS_LAMBDA(int li, int lj, double& acc) {
          const int c = g.getSurfaceIndex(li, lj);
          acc += Kokkos::max(0.0, eta(c) - bottom(c)) * area;
        },
        sw_sum);
    local += sw_sum;
  }
  if (gw_enabled_ && re_) {
    const GwFields& f = re_->fields();
    const double vol = dx_ * dy_ * dz_;
    const Grid g = gw_grid_;
    RealArr1DHost wc = f.wc;
    double gw_sum = 0.0;
    Kokkos::parallel_reduce(
        "water_vol_subs",
        Kokkos::MDRangePolicy<LoopExec, Kokkos::Rank<3>>({0, 0, 0}, {lnx, lny, nz_}),
        KOKKOS_LAMBDA(int li, int lj, int k, double& acc) {
          acc += wc(g.getIndex(li, lj, k)) * vol;
        },
        gw_sum);
    local += gw_sum;
  }
  double global = local;
  if (size_ > 1) MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, mc_->comm());
  return global;
}

void Orchestrator::applySolute(real dt, real rain) {
  if (!solute_ || dt <= 0.0) return;
  if (t_ + dt <= solute_start_time_) return;
  perf::ScopedTimer t(&perf_, perf::Region::Solute);

  const SoluteParams& sp = solute_->params();
  const int lnx = mc_->localNx();
  const int lny = mc_->localNy();
  const double area = dx_ * dy_;

  // Build the operator-split flow snapshot. Surface advection/mixing use the START-of-step depth
  // (= post-step depth minus the rain just added by the SWE solver), so the rain solute exactly
  // balances the rain water and Sum(C*depth) is conserved. Velocities/Darcy fluxes are taken at
  // the current solver state (explicit splitting).
  const int gny = ny_;
  SoluteFlow flow;
  if (sw_enabled_ && swe_) {
    const SweFields& f = swe_->fields();
    const size_t ns = static_cast<size_t>(swe_grid_.nSurfaceStorageCell());
    flow.u = RealArr1DHost("sol_u", ns);
    flow.v = RealArr1DHost("sol_v", ns);
    flow.depth = RealArr1DHost("sol_depth", ns);
    Kokkos::deep_copy(flow.u, f.uu);
    Kokkos::deep_copy(flow.v, f.vv);
    // Reconstruct the start-of-step depth (= post-step depth minus the rain the SWE solver just
    // added) so the rain solute exactly balances the rain water. The SWE adds rain to every
    // interior cell EXCEPT the global y+ boundary row (evaprain, globalJ != gny-1), so that row
    // had no rain added and its start-of-step depth equals its current depth.
    for (int lj = 0; lj < lny; ++lj) {
      const int gj = mc_->j0() + lj;
      const bool got_rain = (gj != gny - 1);
      for (int li = 0; li < lnx; ++li) {
        const int c = swe_grid_.getSurfaceIndex(li, lj);
        const double post = std::max<double>(0.0, f.eta(c) - f.bottom(c));
        flow.depth(static_cast<size_t>(c)) =
            got_rain ? std::max<double>(0.0, post - rain * dt) : post;
      }
    }
  }
  if (gw_enabled_ && re_) {
    const GwFields& f = re_->fields();
    const size_t nc = static_cast<size_t>(gw_grid_.nCell());
    flow.qx = RealArr1DHost("sol_qx", nc);
    flow.qy = RealArr1DHost("sol_qy", nc);
    flow.qz = RealArr1DHost("sol_qz", nc);
    flow.wc = RealArr1DHost("sol_wc", nc);
    Kokkos::deep_copy(flow.qx, f.qx);
    Kokkos::deep_copy(flow.qy, f.qy);
    Kokkos::deep_copy(flow.qz, f.qz);
    Kokkos::deep_copy(flow.wc, f.wc);
  }

  // Rainfall solute source (mass-conservative mixing) applied HERE, not inside the stepper, so
  // it can mirror the SWE rain exclusion at the global y+ row exactly (the stepper takes a scalar
  // rain and would otherwise add solute to a row that received no water). Order is preserved:
  // rain is the first solute operator, so the stepper is then called with rain = 0.
  double rain_added_step = 0.0;
  if (sw_enabled_ && rain > 0.0 && sp.c_rain != 0.0) {
    double local_rain = 0.0;
    const double dw = rain * dt;  // rain water height added to each rained cell
    for (int lj = 0; lj < lny; ++lj) {
      const int gj = mc_->j0() + lj;
      if (gj == gny - 1) continue;  // SWE excludes rain on the global y+ row
      for (int li = 0; li < lnx; ++li) {
        const int c = swe_grid_.getSurfaceIndex(li, lj);
        const double d0 = flow.depth(static_cast<size_t>(c));
        if (d0 <= sp.min_depth) continue;  // skip cells dry at the start of the step
        const double cnew = (d0 * sw_state_->conc(c) + dw * sp.c_rain) / (d0 + dw);
        sw_state_->conc(c) = cnew;
        local_rain += dw * sp.c_rain * area;
      }
    }
    rain_added_step = local_rain;
    if (size_ > 1)
      MPI_Allreduce(&local_rain, &rain_added_step, 1, MPI_DOUBLE, MPI_SUM, mc_->comm());
    solute_added_by_rain_ += rain_added_step;
  }

  const double mass_before = (sp.k_decay > 0.0) ? computeSoluteMass() : 0.0;

  // Substep the solute over the flow step (solute dt = flow dt / substeps). Refusal (CFL) is
  // retried with progressively finer substepping; warn once if even the finest is refused. Rain
  // is already applied above, so the stepper runs with rain = 0 (it still does infiltration
  // mixing / decay / advection / diffusion).
  int nsub = solute_substeps_;
  const int nsub_cap = 4096;
  bool done = false;
  while (!done) {
    const real sub_dt = dt / static_cast<real>(nsub);
    bool refused = false;
    for (int s = 0; s < nsub; ++s) {
      const StepResult res = solute_->step(*sw_state_, *gw_state_, flow, sub_dt, 0.0);
      if (!res.ok) {
        refused = true;
        break;
      }
    }
    if (!refused) {
      done = true;
    } else if (nsub < nsub_cap) {
      nsub *= 2;  // finer substeps to satisfy the advective CFL
    } else {
      if (!solute_cfl_warned_ && rank_ == 0) {
        std::cerr << "frehg2: warning: solute advective CFL exceeded solute.cfl_max even at "
                  << nsub << " substeps; concentration left unchanged this step.\n";
        solute_cfl_warned_ = true;
      }
      done = true;
    }
  }

  if (sp.k_decay > 0.0) {
    const double mass_after = computeSoluteMass();
    // mass_before already includes this step's rain; the stepper's only mass sink is first-order
    // decay (advection/diffusion conserve mass), so the drop is the decay loss.
    solute_decay_loss_ += std::max<double>(0.0, mass_before - mass_after);
  }
}

void Orchestrator::buildMonitors(const Config& config) {
  if (output_interval_ <= 0.0) return;
  const MonitorBundle bundle = parseMonitors(config, sw_enabled_, gw_enabled_);
  if (bundle.probes.empty() && bundle.lines.empty()) return;

  monitor_writer_ = std::make_unique<MonitorWriter>();
  const double x0 = config.getOr<double>("domain.x0", 0.0);
  const double y0 = config.getOr<double>("domain.y0", 0.0);
  monitor_writer_->configure(bundle, out_dir_, sim_id_, nx_, ny_, nz_, dx_, dy_, dz_, x0, y0,
                             botz_, mc_.get());
  monitor_writer_->open(resuming_from_checkpoint_);
}

void Orchestrator::writeMonitors(real time) {
  if (!monitor_writer_ || !monitor_writer_->active()) return;
  const RealArr1DHost* sc = (solute_ && sw_state_) ? &sw_state_->conc : nullptr;
  const RealArr1DHost* gc = (solute_ && gw_state_) ? &gw_state_->conc : nullptr;
  monitor_writer_->writeRow(time, swe_.get(), re_.get(), swe_grid_, gw_grid_, sc, gc);
}

void Orchestrator::buildOutput(const Config& config) {
  out_format_ = config.getOr<std::string>("output.format", "hdf5");
  out_filename_ = config.getOr<std::string>("output.filename", "");
  io_mode_ = ioModeFromString(config.getOr<std::string>("output.io_mode", "serial_gather"));

  // Output field selection. When `output.variables` is absent we write the historical default
  // set (water_depth+eta for SW, head+moisture for GW, plus concentration when solute is on) so
  // existing configs are byte-for-byte unchanged. When present, the list is honored and validated:
  // each name must map to a field of an ENABLED module, otherwise it is a hard error (no silent
  // drop). Canonical field names + accepted aliases are resolved here into sw_/gw_out_vars_.
  sw_out_vars_.clear();
  gw_out_vars_.clear();
  const YAML::Node vars = config.root()["output"]["variables"];
  if (!vars || !vars.IsSequence() || vars.size() == 0) {
    if (sw_enabled_) sw_out_vars_ = {"water_depth", "eta"};
    if (gw_enabled_) gw_out_vars_ = {"head", "moisture"};
    if (solute_enabled_) {
      if (sw_enabled_) sw_out_vars_.push_back("concentration");
      if (gw_enabled_) gw_out_vars_.push_back("concentration");
    }
  } else {
    auto lower = [](std::string s) {
      for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
      return s;
    };
    auto has = [](const std::vector<std::string>& v, const std::string& n) {
      return std::find(v.begin(), v.end(), n) != v.end();
    };
    auto addUnique = [&has](std::vector<std::string>& v, const std::string& n) {
      if (!has(v, n)) v.push_back(n);
    };
    for (const auto& node : vars) {
      const std::string raw = node.as<std::string>();
      const std::string name = lower(raw);
      // Surface fields.
      if (name == "water_depth" || name == "depth") {
        if (!sw_enabled_)
          throw std::runtime_error("output.variables: '" + raw +
                                   "' requires modules.surface_water to be enabled");
        addUnique(sw_out_vars_, "water_depth");
      } else if (name == "eta" || name == "wse" || name == "water_surface_elevation") {
        if (!sw_enabled_)
          throw std::runtime_error("output.variables: '" + raw +
                                   "' requires modules.surface_water to be enabled");
        addUnique(sw_out_vars_, "eta");
      } else if (name == "u" || name == "uu" || name == "velocity_x" || name == "velocity_u") {
        if (!sw_enabled_)
          throw std::runtime_error("output.variables: '" + raw +
                                   "' requires modules.surface_water to be enabled");
        addUnique(sw_out_vars_, "velocity_u");
      } else if (name == "v" || name == "vv" || name == "velocity_y" || name == "velocity_v") {
        if (!sw_enabled_)
          throw std::runtime_error("output.variables: '" + raw +
                                   "' requires modules.surface_water to be enabled");
        addUnique(sw_out_vars_, "velocity_v");
        // Subsurface fields.
      } else if (name == "head" || name == "pressure_head") {
        if (!gw_enabled_)
          throw std::runtime_error("output.variables: '" + raw +
                                   "' requires modules.groundwater to be enabled");
        addUnique(gw_out_vars_, "head");
      } else if (name == "moisture" || name == "water_content" || name == "wc" ||
                 name == "theta") {
        if (!gw_enabled_)
          throw std::runtime_error("output.variables: '" + raw +
                                   "' requires modules.groundwater to be enabled");
        addUnique(gw_out_vars_, "moisture");
      } else if (name == "qx" || name == "qy" || name == "qz") {
        if (!gw_enabled_)
          throw std::runtime_error("output.variables: '" + raw +
                                   "' requires modules.groundwater to be enabled");
        addUnique(gw_out_vars_, name);
        // Concentration applies to whichever enabled module carries solute.
      } else if (name == "concentration" || name == "conc") {
        if (!solute_enabled_)
          throw std::runtime_error("output.variables: '" + raw +
                                   "' requires modules.solute (or solute.enabled) to be enabled");
        if (sw_enabled_) addUnique(sw_out_vars_, "concentration");
        if (gw_enabled_) addUnique(gw_out_vars_, "concentration");
      } else {
        throw std::runtime_error(
            "output.variables: unknown field '" + raw +
            "'. Supported: water_depth, eta, u, v (surface); head, moisture, qx, qy, qz "
            "(subsurface); concentration (solute).");
      }
    }
  }

  if (out_filename_.empty()) {
    out_dir_ = ".";
    summary_path_ = "./simulation_summary.txt";
    return;  // no field-output writer; run still writes the summary.
  }
  const std::string resolved = config.resolvePath(out_filename_);
  out_dir_ = std::filesystem::path(resolved).parent_path().string();
  if (out_dir_.empty()) out_dir_ = ".";
  summary_path_ = out_dir_ + "/simulation_summary.txt";
  mkParentDir(resolved);
  MPI_Barrier(MPI_COMM_WORLD);

  // I/O layout: this rank's owned cells mapped to global physical indices.
  const int lnx = mc_->localNx();
  const int lny = mc_->localNy();
  layout_.nx = nx_;
  layout_.ny = ny_;
  layout_.nz = gw_enabled_ ? nz_ : 1;
  layout_.surf_global_idx =
      IntArr1DHost("surf_idx", static_cast<size_t>(lnx) * static_cast<size_t>(lny));
  for (int lj = 0; lj < lny; ++lj)
    for (int li = 0; li < lnx; ++li) {
      const int gi = mc_->i0() + li;
      const int gj = mc_->j0() + lj;
      layout_.surf_global_idx(static_cast<size_t>(li + lj * lnx)) = gi + gj * nx_;
    }
  if (gw_enabled_) {
    layout_.subs_global_idx = IntArr1DHost(
        "subs_idx", static_cast<size_t>(lnx) * static_cast<size_t>(lny) * static_cast<size_t>(nz_));
    size_t e = 0;
    for (int k = 0; k < nz_; ++k)
      for (int lj = 0; lj < lny; ++lj)
        for (int li = 0; li < lnx; ++li) {
          const int gi = mc_->i0() + li;
          const int gj = mc_->j0() + lj;
          layout_.subs_global_idx(e++) = gi + gj * nx_ + k * nx_ * ny_;
        }
  } else {
    layout_.subs_global_idx = IntArr1DHost("subs_idx", 0);
  }

  writer_ = makeOutputWriter(out_format_, MPI_COMM_WORLD, layout_, io_mode_);
  RunMetadata meta;
  meta.title = config.getOr<std::string>("simulation.title", sim_id_);
  meta.frehg2_version = build_info::kFrehg2Version;
  meta.git_sha = git_sha_;
  meta.git_dirty = build_info::kGitDirty;
  meta.config_sha256 = config_sha256_;
  meta.build_type = build_info::kBuildType;
  meta.compiler = build_info::kCompiler;
  meta.kokkos_backend = Kokkos::DefaultExecutionSpace::name();
  meta.solver_backend = "petsc";
  meta.mpi_ranks = size_;
  writer_->openFile(resolved, meta);
  writer_->writeDomain(gw_enabled_ ? gw_grid_ : swe_grid_);
}

void Orchestrator::swAdvanceTo(real t_target) {
  const real rain = rainAt(t_target);
  const real dt = t_target - t_;  // the surface step this call advances
  {
    std::lock_guard<std::mutex> lock(solve_mu_);
    swe_->advanceStep(rain, evap_rate_);
  }
  // P19 water budget: accumulate this rank's rain inflow. The SWE solver (evaprain) adds rain to
  // every owned cell EXCEPT the global y+ boundary row (globalJ == ny_-1), so mirror that mask.
  if (rain != 0.0 && dt > 0.0) {
    int rain_rows = 0;
    for (int lj = 0; lj < mc_->localNy(); ++lj)
      if (mc_->j0() + lj != ny_ - 1) ++rain_rows;
    water_rain_in_ += static_cast<double>(rain) * dt * dx_ * dy_ *
                      static_cast<double>(mc_->localNx()) * rain_rows;
  }
  // Surface polygon source/BC are explicit post-solve updates (like the legacy explicit rain in
  // evaprain). Folding them here keeps sequential step() and the async swAdvance callback in sync.
  applySurfacePolygons(dt);
}

void Orchestrator::gwCatchUpTo(real t_target) {
  const real eps = 1.0e-9 * std::max<real>(1.0, t_end_);
  const real t_start = gw_time_;
  while (gw_time_ < t_target - eps) {
    {
      std::lock_guard<std::mutex> lock(solve_mu_);
      re_->advanceStep();
    }
    gw_time_ += re_->dtUsed();
    ++gw_substeps_;
  }
  // Subsurface extraction wells: explicit sink over the groundwater time actually advanced this
  // window (0 on a drain() no-op). Folding here keeps wells consistent in sequential and async.
  const real advanced = gw_time_ - t_start;
  if (advanced > 0.0) applySubsurfacePolygons(advanced);
}

void Orchestrator::applyCouplingExchange(real dt) {
  perf::ScopedTimer t(&perf_, perf::Region::Coupling);
  // Capture the per-column flux only when we will move solute with it (solute coupled to both
  // domains); otherwise the plain water-only exchange is used.
  const bool move_solute = solute_ && sw_enabled_ && gw_enabled_;
  RealArr1DHost* qcap = move_solute ? &coupling_flux_ : nullptr;
  exchange_volume_ += coupling_->exchange(*swe_, *re_, dt, qcap);
  if (move_solute) applyCouplingSoluteExchange(dt);
}

void Orchestrator::applyCouplingSoluteExchange(real dt) {
  // Build the per-column signed exchanged height vex (+ infiltration / - seepage) from the
  // captured flux, plus the CURRENT post-exchange surface depth and top-cell water height, then
  // move solute conservatively. Heights use the uniform dz_ basis of computeSoluteMass so the
  // (water move + this transfer) conserves Sum(C*water) to machine precision. Plain host loop
  // (per-column, like the rain-mixing accounting in applySolute).
  const int lnx = mc_->localNx();
  const int lny = mc_->localNy();
  const double area = dx_ * dy_;
  const SweFields& sf = swe_->fields();
  const GwFields& gf = re_->fields();
  RealArr1DHost vex("sol_vex", static_cast<size_t>(swe_grid_.nSurfaceStorageCell()));
  RealArr1DHost depth("sol_ex_depth", static_cast<size_t>(swe_grid_.nSurfaceStorageCell()));
  RealArr1DHost hsub("sol_ex_hsub", static_cast<size_t>(swe_grid_.nSurfaceStorageCell()));
  for (int lj = 0; lj < lny; ++lj)
    for (int li = 0; li < lnx; ++li) {
      const int cs = swe_grid_.getSurfaceIndex(li, lj);
      const int ct = gw_grid_.getIndex(li, lj, 0);
      const double q = coupling_flux_(static_cast<size_t>(li + lj * lnx));  // m^3/s, >0 seepage
      vex(static_cast<size_t>(cs)) = -q * dt / area;  // +ve infiltration (SW->GW)
      depth(static_cast<size_t>(cs)) =
          std::max(0.0, static_cast<double>(sf.eta(cs)) - static_cast<double>(sf.bottom(cs)));
      hsub(static_cast<size_t>(cs)) = static_cast<double>(gf.wc(ct)) * dz_;
    }
  // Host aliases of the canonical (device) conc storage (same space on the host build); writes
  // through them update the underlying State/GwState conc.
  RealArr1DHost surf_conc = sw_state_->conc;
  RealArr1DHost subs_conc = gw_state_->conc;
  applyInterfaceExchange(surf_conc, subs_conc, depth, hsub, vex, gw_grid_);
}

void Orchestrator::applySurfacePolygons(real dt) {
  if (!sw_enabled_ || !swe_) return;
  perf::ScopedTimer t(&perf_, perf::Region::Polygon);
  if (poly_source_ && poly_source_->hasSurface())
    polygon_inflow_volume_ += poly_source_->applySurface(*swe_, dt);
  if (poly_bc_ && !poly_bc_->empty())
    polygon_outflow_volume_ += poly_bc_->applySurface(*swe_, dt);
}

void Orchestrator::applySubsurfacePolygons(real dt) {
  if (!gw_enabled_ || !re_) return;
  perf::ScopedTimer t(&perf_, perf::Region::Polygon);
  if (poly_source_ && poly_source_->hasSubsurface())
    polygon_well_volume_ += poly_source_->applySubsurface(*re_, dt);
}

bool Orchestrator::hasPolygonBc() const { return poly_bc_ != nullptr && !poly_bc_->empty(); }
bool Orchestrator::hasPolygonSource() const {
  return poly_source_ != nullptr && !poly_source_->empty();
}

real Orchestrator::step(real dt_requested) {
  if (coupled_) {
    // Synchronous (Gauss-Seidel) coupled step: SW advance -> exchange -> GW catch-up. The async
    // pipeline (runLoopAsyncCoupled) reorders these across steps for overlap but executes the
    // identical operation sequence, so step() stays the synchronous, fully-advanced semantics
    // every external caller (and every direct test) relies on, in both Sequential and Async mode.
    const real dt_sw = dt_requested;
    const real t_new = t_ + dt_sw;
    swAdvanceTo(t_new);
    applyCouplingExchange(dt_sw);
    gwCatchUpTo(t_new);
    applySolute(dt_sw, rainAt(t_new));  // operator-split transport after the flow step
    t_ = t_new;
    return dt_sw;
  }
  if (sw_enabled_) {
    const real t_new = t_ + dt_requested;
    swAdvanceTo(t_new);  // includes surface polygon source/BC
    applySolute(dt_requested, rainAt(t_new));
    t_ = t_new;
    return dt_requested;
  }
  // Groundwater-only: the RE solver self-adapts dt; the requested dt is advisory.
  re_->advanceStep();
  const real used = re_->dtUsed();
  applySubsurfacePolygons(used);  // extraction wells over the step actually taken
  applySolute(used, 0.0);
  t_ += used;
  gw_time_ = t_;
  ++gw_substeps_;
  return used;
}

void Orchestrator::writeFieldOutputs(real time) {
  if (!writer_) return;
  perf::ScopedTimer t(&perf_, perf::Region::Io);
  const int lnx = mc_->localNx();
  const int lny = mc_->localNy();
  const size_t nsurf = static_cast<size_t>(lnx) * static_cast<size_t>(lny);
  auto want = [](const std::vector<std::string>& v, const char* n) {
    return std::find(v.begin(), v.end(), n) != v.end();
  };
  if (sw_enabled_) {
    const SweFields& f = swe_->fields();
    const Grid g = swe_grid_;
    if (want(sw_out_vars_, "water_depth")) {
      RealArr1DHost depth("depth_owned", nsurf);
      RealArr1DHost fe = f.eta, fb = f.bottom;
      parallelForSurface("io_pack_depth", lnx, lny, KOKKOS_LAMBDA(int li, int lj) {
        const int c = g.getSurfaceIndex(li, lj);
        depth(static_cast<size_t>(li + lj * lnx)) = Kokkos::max(fe(c) - fb(c), 0.0);
      });
      writer_->writeSurfaceField("water_depth", time, depth, "m");
    }
    if (want(sw_out_vars_, "eta")) {
      RealArr1DHost eta("eta_owned", nsurf);
      RealArr1DHost fe = f.eta;
      parallelForSurface("io_pack_eta", lnx, lny, KOKKOS_LAMBDA(int li, int lj) {
        eta(static_cast<size_t>(li + lj * lnx)) = fe(g.getSurfaceIndex(li, lj));
      });
      writer_->writeSurfaceField("eta", time, eta, "m");
    }
    if (want(sw_out_vars_, "velocity_u")) {
      RealArr1DHost u("u_owned", nsurf);
      RealArr1DHost fu = f.uu;
      parallelForSurface("io_pack_u", lnx, lny, KOKKOS_LAMBDA(int li, int lj) {
        u(static_cast<size_t>(li + lj * lnx)) = fu(g.getSurfaceIndex(li, lj));
      });
      writer_->writeSurfaceField("velocity_u", time, u, "m/s");
    }
    if (want(sw_out_vars_, "velocity_v")) {
      RealArr1DHost v("v_owned", nsurf);
      RealArr1DHost fv = f.vv;
      parallelForSurface("io_pack_v", lnx, lny, KOKKOS_LAMBDA(int li, int lj) {
        v(static_cast<size_t>(li + lj * lnx)) = fv(g.getSurfaceIndex(li, lj));
      });
      writer_->writeSurfaceField("velocity_v", time, v, "m/s");
    }
  }
  if (gw_enabled_) {
    const size_t n = nsurf * static_cast<size_t>(nz_);
    const GwFields& f = re_->fields();
    const Grid g = gw_grid_;
    if (want(gw_out_vars_, "head")) {
      RealArr1DHost head("head_owned", n);
      RealArr1DHost fh = f.h;
      parallelForVolume("io_pack_head", lnx, lny, nz_, KOKKOS_LAMBDA(int li, int lj, int k) {
        head(static_cast<size_t>(li + lj * lnx + k * lnx * lny)) = fh(g.getIndex(li, lj, k));
      });
      writer_->writeSubsurfaceField("head", time, head, "m");
    }
    if (want(gw_out_vars_, "moisture")) {
      RealArr1DHost moist("moist_owned", n);
      RealArr1DHost fwc = f.wc;
      parallelForVolume("io_pack_moist", lnx, lny, nz_, KOKKOS_LAMBDA(int li, int lj, int k) {
        moist(static_cast<size_t>(li + lj * lnx + k * lnx * lny)) = fwc(g.getIndex(li, lj, k));
      });
      writer_->writeSubsurfaceField("moisture", time, moist, "");
    }
    // Darcy face volumetric fluxes [m^3/s] (qx/qy are 0 outside the full-3D path; qz is always
    // the corrector's vertical flux). Written only when requested.
    const char* qnames[3] = {"qx", "qy", "qz"};
    const RealArr1DHost qsrc[3] = {f.qx, f.qy, f.qz};
    for (int qi = 0; qi < 3; ++qi) {
      if (!want(gw_out_vars_, qnames[qi])) continue;
      RealArr1DHost q("q_owned", n);
      RealArr1DHost src = qsrc[qi];
      parallelForVolume("io_pack_flux", lnx, lny, nz_, KOKKOS_LAMBDA(int li, int lj, int k) {
        q(static_cast<size_t>(li + lj * lnx + k * lnx * lny)) = src(g.getIndex(li, lj, k));
      });
      writer_->writeSubsurfaceField(qnames[qi], time, q, "m^3/s");
    }
  }
  if (solute_) {
    if (sw_enabled_ && want(sw_out_vars_, "concentration")) {
      RealArr1DHost cs("conc_surf_owned", nsurf);
      const Grid g = swe_grid_;
      RealArr1DHost conc = sw_state_->conc;
      parallelForSurface("io_pack_conc_surf", lnx, lny, KOKKOS_LAMBDA(int li, int lj) {
        cs(static_cast<size_t>(li + lj * lnx)) = conc(g.getSurfaceIndex(li, lj));
      });
      writer_->writeSurfaceField("concentration", time, cs, "kg/m^3");
    }
    if (gw_enabled_ && want(gw_out_vars_, "concentration")) {
      const size_t n = nsurf * static_cast<size_t>(nz_);
      RealArr1DHost cg("conc_subs_owned", n);
      const Grid g = gw_grid_;
      RealArr1DHost conc = gw_state_->conc;
      parallelForVolume("io_pack_conc_subs", lnx, lny, nz_, KOKKOS_LAMBDA(int li, int lj, int k) {
        cg(static_cast<size_t>(li + lj * lnx + k * lnx * lny)) = conc(g.getIndex(li, lj, k));
      });
      writer_->writeSubsurfaceField("concentration", time, cg, "kg/m^3");
    }
  }
  last_output_t_ = time;
  ++outputs_written_;
  writeMonitors(time);
}

void Orchestrator::writeCheckpointNow() {
  if (!writer_) return;
  perf::ScopedTimer t(&perf_, perf::Region::Io);
  CheckpointState cp;
  cp.time = t_;
  cp.dt = gw_enabled_ ? re_->dt() : dt_;
  cp.has_sw = sw_enabled_;
  cp.has_gw = gw_enabled_;
  cp.has_solute = solute_ != nullptr;
  cp.config_sha256 = config_sha256_;
  cp.git_sha = git_sha_;
  cp.storage_layout = "internal-with-halo";

  auto deepCopy = [](const RealArr1DHost& src, const std::string& nm) {
    RealArr1DHost dst(nm, src.extent(0));
    Kokkos::deep_copy(dst, src);
    return dst;
  };
  if (sw_enabled_) {
    SweFields& f = swe_->fields();
    cp.eta = deepCopy(f.eta, "eta");
    cp.u = deepCopy(f.uu, "u");
    cp.v = deepCopy(f.vv, "v");
    for (auto& nv : f.namedViews()) cp.extra["sw_" + nv.first] = deepCopy(*nv.second, nv.first);
  }
  if (gw_enabled_) {
    GwFields& f = re_->fields();
    cp.h = deepCopy(f.h, "h");
    cp.hn = deepCopy(f.hn, "hn");
    cp.wc = deepCopy(f.wc, "wc");
    cp.wcn = deepCopy(f.wcn, "wcn");
    for (auto& nv : f.namedViews()) cp.extra["gw_" + nv.first] = deepCopy(*nv.second, nv.first);
  }
  if (solute_) {
    // Solute conc is full halo-padded state; store in extra for bit-exact restart and also set
    // the documented cp.C (surface preferred, else subsurface) for the schema / h5py interop.
    if (sw_enabled_) cp.extra["solute_sw_conc"] = deepCopy(sw_state_->conc, "conc");
    if (gw_enabled_) cp.extra["solute_gw_conc"] = deepCopy(gw_state_->conc, "conc");
    cp.C = deepCopy(sw_enabled_ ? sw_state_->conc : gw_state_->conc, "C");
  }
  writer_->writeCheckpoint(cp);
}

void Orchestrator::runLoop() {
  // Async coupling (P11) is the double-buffered pipeline; everything else (sequential coupling,
  // SW-only, GW-only) uses the straight-line loop. Multi-rank async falls back to sequential
  // (warned at initialize()) until the P10 PetscSubcomm split makes concurrent solves safe.
  // Solute transport (P16) requires a synchronized flow state for operator splitting (it runs
  // inside the synchronous step()), so it forces the sequential loop. The async pipeline (P11)
  // lags the GW window, which has no defined transport semantics yet.
  if (asyncActive() && !solute_) {
    runLoopAsyncCoupled();
  } else {
    runLoopSequential();
  }
}

void Orchestrator::runLoopSequential() {
  const auto wall_start = std::chrono::steady_clock::now();
  const real eps = 1.0e-9 * std::max<real>(1.0, t_end_);
  if (solute_) solute_initial_mass_ = computeSoluteMass();

  while (t_ < t_end_ - eps && (max_steps_ <= 0 || step_count_ < max_steps_)) {
    step(dt_);
    ++step_count_;
    if (t_ >= next_output_t_ - eps) {
      writeFieldOutputs(t_);
      next_output_t_ = nextMultipleAfter(output_interval_, t_);
    }
    if (dt_checkpoint_ > 0.0 && t_ >= next_checkpoint_t_ - eps) {
      writeCheckpointNow();
      next_checkpoint_t_ = nextMultipleAfter(dt_checkpoint_, t_);
    }
  }

  // Always emit a final snapshot at the end state (if not already written there).
  if (writer_ && last_output_t_ < t_ - eps) {
    writeFieldOutputs(t_);
  }

  if (solute_) solute_final_mass_ = computeSoluteMass();
  const auto wall_end = std::chrono::steady_clock::now();
  const double wall_s =
      std::chrono::duration<double>(wall_end - wall_start).count();
  if (writer_) writer_->close();
  if (monitor_writer_) monitor_writer_->close();
  writeSummary(wall_s);
}

void Orchestrator::runLoopAsyncCoupled() {
  const auto wall_start = std::chrono::steady_clock::now();
  const real eps = 1.0e-9 * std::max<real>(1.0, t_end_);

  AsyncPipeline::Callbacks cb;
  cb.swAdvance = [this](real t_target) { swAdvanceTo(t_target); };
  cb.gwCatchUp = [this](real t_target) { gwCatchUpTo(t_target); };
  cb.applyExchange = [this](real dt) { applyCouplingExchange(dt); };
  AsyncPipeline pipeline(std::move(cb));

  while (t_ < t_end_ - eps && (max_steps_ <= 0 || step_count_ < max_steps_)) {
    // SW advances [t_, t_+dt_] now; the groundwater catch-up of the previous window overlaps on
    // a worker thread. After this returns, SW is at t_+dt_ and GW is one window behind.
    pipeline.advanceWindow(t_, dt_);
    t_ += dt_;
    ++step_count_;

    // Output / checkpoint read BOTH domains at the same time, so drain the trailing GW window
    // first to synchronize them. drain() is idempotent if no window is pending.
    if (t_ >= next_output_t_ - eps) {
      pipeline.drain();
      writeFieldOutputs(t_);
      next_output_t_ = nextMultipleAfter(output_interval_, t_);
    }
    if (dt_checkpoint_ > 0.0 && t_ >= next_checkpoint_t_ - eps) {
      pipeline.drain();
      writeCheckpointNow();
      next_checkpoint_t_ = nextMultipleAfter(dt_checkpoint_, t_);
    }
  }

  // Finish the trailing GW window so the end state is fully synchronized before final output.
  pipeline.drain();
  if (writer_ && last_output_t_ < t_ - eps) {
    writeFieldOutputs(t_);
  }

  const auto wall_end = std::chrono::steady_clock::now();
  const double wall_s =
      std::chrono::duration<double>(wall_end - wall_start).count();
  if (writer_) writer_->close();
  if (monitor_writer_) monitor_writer_->close();
  writeSummary(wall_s);
}

void Orchestrator::run() {
  if (!resuming_from_checkpoint_) {
    t_ = 0.0;
    gw_time_ = 0.0;
    step_count_ = 0;
    gw_substeps_ = 0;
    outputs_written_ = 0;
    exchange_volume_ = 0.0;
    polygon_outflow_volume_ = 0.0;
    polygon_inflow_volume_ = 0.0;
    polygon_well_volume_ = 0.0;
    solute_added_by_rain_ = 0.0;
    solute_decay_loss_ = 0.0;
    next_output_t_ = output_interval_ > 0.0 ? output_interval_ : t_end_ + 1.0;
    next_checkpoint_t_ = dt_checkpoint_ > 0.0 ? dt_checkpoint_ : t_end_ + 1.0;
    writeFieldOutputs(0.0);  // initial state
  } else {
    resuming_from_checkpoint_ = false;
  }
  runLoop();
}

void Orchestrator::restart(const std::string& checkpoint_file, real checkpoint_time) {
  loadCheckpointState(checkpoint_file, checkpoint_time);
  runLoop();
}

void Orchestrator::writeSummary(double wall_seconds) {
  // Water mass-balance reductions: ALL ranks must participate before the rank-0 write guard.
  const double water_final = computeWaterVolume();  // already MPI-reduced (global)
  double rain_in = water_rain_in_, outflow = polygon_outflow_volume_,
         inflow = polygon_inflow_volume_, well = polygon_well_volume_, exch = exchange_volume_;
  if (size_ > 1) {
    double loc[5] = {rain_in, outflow, inflow, well, exch};
    double glb[5] = {0, 0, 0, 0, 0};
    MPI_Allreduce(loc, glb, 5, MPI_DOUBLE, MPI_SUM, mc_->comm());
    rain_in = glb[0]; outflow = glb[1]; inflow = glb[2]; well = glb[3]; exch = glb[4];
  }

  perf_.mpiReduce(mc_->comm());

  if (rank_ != 0) return;
  std::ofstream out(summary_path_);
  if (!out) return;
  auto boolWord = [](bool b) { return b ? "true" : "false"; };
  out << "frehg2_version \"" << build_info::kFrehg2Version << " (" << git_sha_ << ")\"\n";
  out << "simulation_id \"" << sim_id_ << "\"\n";
  out << "modules {surface_water=" << boolWord(sw_enabled_)
      << ", groundwater=" << boolWord(gw_enabled_)
      << ", solute=" << boolWord(solute_enabled_) << "}\n";
  out << "mpi_ranks " << size_ << "\n";
  out << "kokkos_execution_space \"" << Kokkos::DefaultExecutionSpace::name() << "\"\n";
  out << "total_runtime_seconds " << wall_seconds << "\n";
  out << "simulated_time_seconds " << static_cast<double>(t_) << "\n";
  out << "time_steps_completed " << step_count_ << "\n";
  out << "gw_substeps_completed " << gw_substeps_ << "\n";
  out << "output_intervals_completed " << outputs_written_ << "\n";

  // P19 water mass-balance budget. For an SW-only run with no groundwater flux BCs the budget is
  // closable: final == initial + rain_in + polygon_inflow - polygon_outflow (the harness checks
  // the residual). When groundwater is present the GW flux-BC inflow is not instrumented, so the
  // residual is informational only (harness flags `groundwater_present`).
  out << "water_initial_volume " << water_initial_volume_ << "\n";
  out << "water_final_volume " << water_final << "\n";
  out << "water_rain_volume_in " << rain_in << "\n";
  out << "water_polygon_inflow_volume " << inflow << "\n";
  out << "water_polygon_outflow_volume " << outflow << "\n";
  out << "water_polygon_well_volume " << well << "\n";
  out << "coupling_exchange_volume " << exch << "\n";
  out << "groundwater_present " << boolWord(gw_enabled_) << "\n";
  out << "surface_water_present " << boolWord(sw_enabled_) << "\n";
  if (solute_enabled_) {
    const double expected = solute_initial_mass_ + solute_added_by_rain_ - solute_decay_loss_;
    const double denom = std::max(std::abs(solute_initial_mass_ + solute_added_by_rain_), 1.0e-300);
    const double rel_err = std::abs(solute_final_mass_ - expected) / denom;
    out << "solute_initial_mass " << solute_initial_mass_ << "\n";
    out << "solute_final_mass " << solute_final_mass_ << "\n";
    out << "solute_added_by_rain " << solute_added_by_rain_ << "\n";
    out << "solute_decay_loss " << solute_decay_loss_ << "\n";
    out << "solute_relative_mass_error " << rel_err << "\n";
  }
  perf_.writeSummaryLines(out, step_count_, gw_substeps_);
  out << "wall_clock_hours " << (wall_seconds / 3600.0) << "\n";
}

}  // namespace frehg2
