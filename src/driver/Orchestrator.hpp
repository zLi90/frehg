// General-purpose production driver (P7).
//
// The Orchestrator is the ONLY production code path from main() to the physics. It builds a
// Grid, model-owned MPI decomposition, the backend-agnostic linear solver(s), the SWE / RE
// solvers, and the synchronous Coupling entirely from a single YAML Config, then runs the
// coupled time loop with output, checkpointing, restart, and a run summary. Every feature
// (BCs, ICs, sources, output, coupling) is reachable from run(); nothing is wired only in a
// test fixture.
//
// Architecture note (deviation from the P7 plan sketch): the plan listed Domain/GwDomain/
// State/GwState members, but in the realized P4/P5 design each solver owns its own
// halo-padded fields (SweFields / GwFields). The Orchestrator therefore owns Grid, MpiComm,
// Decomp2D/Decomp3D, the LinearSolver backends, the solvers, and the Coupling — matching the
// actual solver architecture. No PETSc type appears here: the linear backend is created via
// the backend-agnostic makeLinearSolver() factory.
#ifndef FREHG2_DRIVER_ORCHESTRATOR_HPP
#define FREHG2_DRIVER_ORCHESTRATOR_HPP

#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "core/Grid.hpp"
#include "frehg2/core/define.hpp"
#include "frehg2/io/OutputWriter.hpp"
#include "frehg2/perf/PerfRecorder.hpp"

namespace frehg2 {

class Config;
class MpiComm;
class Decomp2D;
class Decomp3D;
class LinearSolver;
class SweSolver;
class ReSolver;
class Coupling;
class PolygonBC;
class PolygonSource;
class SoilMap;
class MonitorWriter;
class SoluteStepper;
class State;
class GwState;

class Orchestrator {
 public:
  // SW<->GW coupling execution mode (P11). Sequential is the P6 blocking-sync path (the
  // regression-safe default for b1-sw / b2-gw); Async is the P11 double-buffered pipeline,
  // selected by the YAML key `coupling.mode: async`. Async is bit-identical to Sequential.
  enum class CouplingMode { Sequential, Async };

  Orchestrator();
  ~Orchestrator();

  Orchestrator(const Orchestrator&) = delete;
  Orchestrator& operator=(const Orchestrator&) = delete;

  // Build every component from the resolved config. Throws std::runtime_error on an invalid
  // config (missing required field, file not found, inconsistent module/mode selection).
  void initialize(const Config& config);

  // Advance one coupled time step (the requested dt is the surface step; adaptive GW may
  // take several substeps to catch up). Returns the surface dt actually used.
  real step(real dt_requested);

  // Run to completion: step() loop until t >= t_end (or max_steps), writing output at
  // output_interval, checkpoints at dt_checkpoint, and simulation_summary.txt at the end.
  void run();

  // Restart from a checkpoint file: requires initialize() first (to build the components),
  // then loads the full solver state, sets the clock to checkpoint_time, and runs to t_end.
  void restart(const std::string& checkpoint_file, real checkpoint_time);

  // --- Diagnostics / test accessors -------------------------------------------------------
  const SweSolver* swe() const { return swe_.get(); }
  const ReSolver* re() const { return re_.get(); }
  bool swEnabled() const { return swe_ != nullptr; }
  bool gwEnabled() const { return re_ != nullptr; }
  bool coupled() const { return coupling_ != nullptr; }
  CouplingMode couplingMode() const { return coupling_mode_; }
  // Whether the run() loop will actually use the async pipeline (async requested AND coupled
  // AND single rank; multi-rank async falls back to synchronous coupling until P10).
  bool asyncActive() const {
    return coupling_mode_ == CouplingMode::Async && coupling_ != nullptr && size_ == 1;
  }
  real time() const { return t_; }
  real gwTime() const { return gw_time_; }
  long long stepCount() const { return step_count_; }
  long long gwSubstepCount() const { return gw_substeps_; }
  int outputsWritten() const { return outputs_written_; }
  const std::string& summaryPath() const { return summary_path_; }
  // Net signed exchange volume [m^3] accumulated over the run (+ = GW->SW seepage).
  real exchangeVolume() const { return exchange_volume_; }
  // Polygon BC / source diagnostics accumulated over the run (P12). Outflow = volume removed by
  // boundary regions; inflow = surface source volume added; well = subsurface volume extracted.
  real polygonOutflowVolume() const { return polygon_outflow_volume_; }
  real polygonInflowVolume() const { return polygon_inflow_volume_; }
  real polygonWellVolume() const { return polygon_well_volume_; }
  bool hasPolygonBc() const;
  bool hasPolygonSource() const;

  // Solute diagnostics (P16). soluteEnabled() reflects the resolved solute.enabled.
  bool soluteEnabled() const { return solute_ != nullptr; }
  double soluteInitialMass() const { return solute_initial_mass_; }
  double soluteFinalMass() const { return solute_final_mass_; }
  double soluteAddedByRain() const { return solute_added_by_rain_; }
  double soluteDecayLoss() const { return solute_decay_loss_; }
  // Max concentration over OWNED cells (per-rank; for single-rank tests). 0 when solute disabled.
  double maxSurfaceConc() const;
  double maxSubsurfaceConc() const;
  // Owned (no-halo) concentration snapshots for tests (row-major gi+gj*nx / +k*nx*ny). Empty
  // views when solute is disabled or the corresponding module is off.
  RealArr1DHost soluteSurfaceConcOwned() const;
  RealArr1DHost soluteSubsurfaceConcOwned() const;

 private:
  void buildGrids(const Config& config);
  void buildSurfaceWater(const Config& config);
  void buildGroundwater(const Config& config);
  // Build + attach the non-uniform SoilMap (P13). No-op (leaves the RE solver on the uniform
  // path) when no per-cell soil class map is configured.
  void buildSoilMap(const Config& config);
  void buildRecharge(const Config& config);
  void buildPolygons(const Config& config);
  // Apply flexible initial conditions once at init (P14). Restart IC loads checkpoint state.
  void setupInitialConditions(const Config& config);
  void loadCheckpointState(const std::string& checkpoint_file, real checkpoint_time);
  // Build the solute transport stepper (P16). No-op (leaves solute_ null) unless
  // `solute.enabled` is true, so a flow-only run does NO solute work / PETSc assembly.
  void buildSolute(const Config& config);
  // Run one operator-split solute step after the flow step, using the start-of-step flow
  // snapshot. No-op when solute_ is null. Accumulates the rain/decay mass-balance terms.
  void applySolute(real dt, real rain);
  // Total solute mass over OWNED cells: surface Sum(C*depth*area) + subsurface Sum(C*wc*vol),
  // reduced across ranks. Returns 0 when solute disabled.
  double computeSoluteMass() const;
  // Total water volume over OWNED cells [m^3]: surface Sum(depth*area) + subsurface Sum(wc*vol),
  // reduced across ranks. Used for the P19 water mass-balance budget in simulation_summary.txt.
  double computeWaterVolume() const;
  void buildOutput(const Config& config);
  void buildMonitors(const Config& config);
  void writeMonitors(real time);
  void loadRainfall(const Config& config);
  real rainAt(real t) const;
  void mkParentDir(const std::string& path) const;
  void writeFieldOutputs(real time);
  void writeCheckpointNow();
  void writeSummary(double wall_seconds);
  void runLoop();
  void runLoopSequential();
  void runLoopAsyncCoupled();
  // Per-domain coupling operations shared by the sequential step() and the async pipeline so
  // both paths execute the identical Gauss-Seidel sequence (SW advance -> exchange -> GW
  // catch-up). swAdvanceTo / gwCatchUpTo take solve_mu_ around their KSP solves (PETSc on a
  // shared communicator is not thread-safe); applyCouplingExchange runs after the join (no lock).
  void swAdvanceTo(real t_target);
  void gwCatchUpTo(real t_target);
  void applyCouplingExchange(real dt);
  // Move dissolved solute mass with the water the coupling just exchanged (co-located with
  // applyCouplingExchange so the post-exchange water heights are exact -> machine-precision
  // conservation). No-op unless solute is enabled on both coupled domains.
  void applyCouplingSoluteExchange(real dt);
  // Apply polygon source/BC regions (no-op when none). Surface regions run after the surface
  // solve (folded into swAdvanceTo); subsurface wells run after the GW catch-up.
  void applySurfacePolygons(real dt);
  void applySubsurfacePolygons(real dt);
  real nextMultipleAfter(real interval, real t) const;

  // Resolved configuration / run parameters.
  std::string sim_id_;
  std::string mode_;
  bool sw_enabled_ = false;
  bool gw_enabled_ = false;
  bool solute_enabled_ = false;
  bool coupled_ = false;
  CouplingMode coupling_mode_ = CouplingMode::Sequential;

  int nx_ = 0, ny_ = 0, nz_ = 0;
  real dx_ = 1.0, dy_ = 1.0, dz_ = 1.0, dz_incre_ = 1.0, botz_ = 0.0;

  real dt_ = 1.0;
  real t_end_ = 0.0;
  long long max_steps_ = 0;  // 0 => no cap (t_end governs)
  real output_interval_ = 0.0;
  real dt_checkpoint_ = 0.0;
  int max_checkpoints_ = 2;

  real surface_dt_ = 0.0;
  real groundwater_dt_ = 0.0;
  real evap_rate_ = 0.0;

  // Rainfall forcing.
  bool rain_from_file_ = false;
  real rain_rate_const_ = 0.0;
  std::vector<double> rain_t_;
  std::vector<double> rain_v_;

  // Output.
  std::string out_format_ = "hdf5";
  std::string out_filename_;
  std::string out_dir_;
  std::string summary_path_;
  IoMode io_mode_ = IoMode::SerialGather;
  std::string config_sha256_;
  std::string git_sha_;
  std::vector<std::string> sw_out_vars_;
  std::vector<std::string> gw_out_vars_;

  // Clock / counters.
  real t_ = 0.0;
  real gw_time_ = 0.0;
  long long step_count_ = 0;
  long long gw_substeps_ = 0;
  int outputs_written_ = 0;
  real next_output_t_ = 0.0;
  real next_checkpoint_t_ = 0.0;
  real last_output_t_ = -1.0e30;
  bool resuming_from_checkpoint_ = false;
  real exchange_volume_ = 0.0;

  // Components (model-owned; no PETSc DMDA, no solver-library type here).
  std::unique_ptr<MpiComm> mc_;
  Grid swe_grid_;
  Grid gw_grid_;
  std::unique_ptr<Decomp2D> dd2_;
  std::unique_ptr<Decomp3D> dd3_;
  std::unique_ptr<LinearSolver> sw_linear_;
  std::unique_ptr<LinearSolver> gw_linear_;
  std::unique_ptr<SweSolver> swe_;
  std::unique_ptr<ReSolver> re_;
  std::unique_ptr<SoilMap> soil_map_;
  std::unique_ptr<Coupling> coupling_;
  std::unique_ptr<PolygonBC> poly_bc_;
  std::unique_ptr<PolygonSource> poly_source_;
  std::unique_ptr<OutputWriter> writer_;
  std::unique_ptr<MonitorWriter> monitor_writer_;
  IoLayout layout_;

  // Solute transport (P16). The conc fields live in sw_state_/gw_state_ (the canonical P2
  // State/GwState conc storage); the stepper drives them through its own diffusion backends.
  std::unique_ptr<SoluteStepper> solute_;
  std::unique_ptr<State> sw_state_;
  std::unique_ptr<GwState> gw_state_;
  std::unique_ptr<LinearSolver> solute_surf_linear_;
  std::unique_ptr<LinearSolver> solute_subs_linear_;
  // Per-column SW<->GW exchange flux [m^3/s] captured from the coupling (owned index i+j*nx),
  // used to move solute with the exchanged water. Allocated only when solute is coupled.
  RealArr1DHost coupling_flux_;
  int solute_substeps_ = 1;             // solute.substeps (>=1): solute dt = flow dt / substeps
  real solute_start_time_ = 0.0;        // solute.start_time: transport begins once t >= this
  bool solute_cfl_warned_ = false;      // emit the CFL-refusal warning at most once
  // Mass-balance accounting (written to simulation_summary.txt).
  double solute_initial_mass_ = 0.0;
  double solute_final_mass_ = 0.0;
  double solute_added_by_rain_ = 0.0;
  double solute_decay_loss_ = 0.0;

  real polygon_outflow_volume_ = 0.0;
  real polygon_inflow_volume_ = 0.0;
  real polygon_well_volume_ = 0.0;

  // Water mass-balance accounting (P19; written to simulation_summary.txt). Initial volume is
  // snapshotted after IC; rain inflow is accumulated as the SWE solver applies it.
  double water_initial_volume_ = 0.0;
  double water_rain_in_ = 0.0;

  int rank_ = 0;
  int size_ = 1;

  // Serializes the SWE and RE KSP solves when the async pipeline runs them on separate threads
  // (PETSc KSPSolve on a shared communicator is not thread-safe). Unused on the sequential path.
  std::mutex solve_mu_;

  // P21 performance instrumentation (always active; written to simulation_summary.txt).
  PerfRecorder perf_;
};

}  // namespace frehg2

#endif  // FREHG2_DRIVER_ORCHESTRATOR_HPP
