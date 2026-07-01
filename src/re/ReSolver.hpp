// Predictor-corrector Richards-equation solver (P5) — legacy-exact port of the b2-gw /
// iter_solve==0, use_corrector==1 path from legacy/frehg/src/groundwater.c.
#ifndef FREHG2_RE_RE_SOLVER_HPP
#define FREHG2_RE_RE_SOLVER_HPP

#include <memory>
#include <vector>

#include "core/Grid.hpp"
#include "core/MpiComm.hpp"
#include "frehg2/core/define.hpp"
#include "re/GwFields.hpp"
#include "re/VanGenuchten.hpp"

namespace frehg2 {

class Config;
class LinearSolver;
class SparseSystem;
class Decomp3D;
class SoilMap;
class PerfRecorder;

struct ReParams {
  SoilParams soil;
  real dx = 1.0;
  real dy = 1.0;
  real dz = 0.01;
  real botz = -1.0;
  real dt = 1.0e-4;
  real dt_min = 1.0e-4;
  real dt_max = 2.0;
  real co_max = 2.0;
  bool adaptive_dt = true;
  bool use_corrector = true;
  bool use_full3d = false;
  bool follow_terrain = false;
  bool sim_shallowwater = false;
  bool baroclinic = false;
  // legacy bctype_GW[0..5]: x-, x+, y-, y+, z+ (bottom), z- (top)
  int bc_type[6] = {0, 0, 0, 0, 0, 1};
  real htop = 0.0;   // top (z-) Dirichlet head, used when bc_type[5]==1
  real hbot = 0.0;   // bottom (z+) Dirichlet head, used when bc_type[4]==1
  real qtop = 0.0;   // top fixed flux (m/s, <0 = downward infiltration), used when bc_type[5]==2
  real qbot = 0.0;   // bottom fixed flux (m/s, >0 = upward inflow), used when bc_type[4]==2
};

class ReSolver {
 public:
  explicit ReSolver(const Grid& grid, const MpiComm* mc = nullptr);
  ~ReSolver();

  void setParams(const ReParams& p) { params_ = p; }
  const ReParams& params() const { return params_; }

  // Attach a spatially variable soil map (P13). The map's local nx,ny must match this rank's
  // grid. When unset, every cell uses the uniform `params_.soil` (bit-identical to P5). Must be
  // set BEFORE initializeUniformColumn so per-column initial heads use the per-column soil.
  void setSoilMap(const SoilMap* map) { soil_map_ = map; }
  const SoilMap* soilMap() const { return soil_map_; }

  // Public per-cell soil accessor (the heterogeneous-soil-aware value for cell (i,j,k), or the
  // uniform params_.soil when no map is attached). External consumers (e.g. polygon extraction
  // wells) MUST use this rather than params().soil so they respect a P13/P23 soil map.
  const SoilParams& soilParamsAt(int i, int j, int k) const { return soilAt(i, j, k); }

  void attachSolver(LinearSolver& solver, Decomp3D& dd);

  // P21 performance instrumentation (optional; Orchestrator attaches its recorder).
  void setPerfRecorder(PerfRecorder* rec) { perf_ = rec; }

  // Uniform column geometry + init_wc IC (legacy initialize.c init_wc path).
  void initializeUniformColumn(real init_wc);

  // P14 flexible IC: build column geometry only (no wc/h), then set fields and finalize.
  void initializeGeometry();
  void setInitialWaterContent(const std::vector<double>& wc_local);
  void setInitialHead(const std::vector<double>& h_local);
  void finalizeInitialState();

  // Advance one legacy PCA groundwater step; updates params_.dt when adaptive_dt.
  void advanceStep();

  const GwFields& fields() const { return fields_; }
  GwFields& fields() { return fields_; }
  const Grid& grid() const { return grid_; }
  real dt() const { return params_.dt; }
  // dt applied in the last completed advanceStep() (before adaptiveTimeStep).
  real dtUsed() const { return dt_used_; }
  // Force the next step's dt (used for legacy dt-sequence replay parity tests).
  void setDt(real dt) { params_.dt = dt; }

  // Restore the adaptive-dt state on restart (P7.4): at a step boundary dtn_ == params_.dt,
  // so priming both reproduces the continuous run's next-step dt exactly.
  void primeAdaptiveDt(real dt) {
    params_.dt = dt;
    dtn_ = dt;
  }

  // Per-column top fixed-flux (legacy `qtop`, used only when bc_type[5]==2). `q_local_interior`
  // holds one value per OWNED column (size nx*ny, indexed i + j*nx); the legacy sign convention
  // is negative = downward infiltration into the domain (b3-kirkland fixed-flux). When unset, the
  // uniform scalar `params_.qtop` applies to every column. This is the partial-width recharge BC
  // (legacy reads it from a per-top-cell polygon); columns outside the recharge region get 0.
  void setTopFluxField(const std::vector<real>& q_local_interior);

 private:
  int s(int i, int j, int k) const { return grid_.getIndex(i, j, k); }
  // Per-cell soil parameters: the SoilMap class for cell (i,j,k), or the uniform params_.soil
  // when no map is attached. A 2D (per-column) map ignores k (P13); a 3D map is fully
  // heterogeneous per cell (P23, SERGHEI-style). (i,j,k) are clamped into the owned interior so
  // ghost neighbor lookups at domain edges resolve to the adjacent interior cell.
  const SoilParams& soilAt(int i, int j, int k) const;
  bool active(int i, int j, int k) const;
  bool isTop(int k) const { return k == 0; }
  bool isBottom(int k) const { return k == grid_.nz() - 1; }
  real cellVolume(int i, int j, int k) const;
  real dzCell(int k) const;

  void buildColumnGeometry();
  void savePreviousStepFields();
  void updateDiagnostics();
  void computeKFace();
  void baroclinicFace();
  void groundwaterMatCoeff();
  void groundwaterRhs();
  void assembleAndSolve();
  void enforceHeadBc();
  void groundwaterFlux();
  void checkRoom();
  void updateWaterContent();
  void reallocateWaterContent();
  void clampWaterContent();
  void enforceMoistureBc();
  void adaptiveTimeStep();

  Grid grid_;
  ReParams params_;
  GwFields fields_;
  const MpiComm* mc_ = nullptr;
  const SoilMap* soil_map_ = nullptr;

  LinearSolver* solver_ = nullptr;
  Decomp3D* dd_ = nullptr;
  std::unique_ptr<SparseSystem> sys_;
  RealArr1D b_owned_;
  RealArr1D x_owned_;
  real dtn_ = 0.0;
  real dt_used_ = 0.0;
  PerfRecorder* perf_ = nullptr;

  // Per-column top fixed-flux field (legacy qtop), indexed by grid_.getSurfaceIndex(i,j).
  // Empty/unset => use the uniform params_.qtop. Only consulted when bc_type[5]==2.
  RealArr1DHost qtop_field_;
  bool has_qtop_field_ = false;
};

}  // namespace frehg2

#endif  // FREHG2_RE_RE_SOLVER_HPP
