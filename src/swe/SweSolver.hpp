// Semi-implicit shallow-water solver (P4) — legacy-exact port of the b1-sw path from
// legacy/frehg/src/shallowwater.c (+ geometry from initialize.c).
//
// Bring-up follows the mandated P4.S sequence (4a serial scalar core → 4e MPI). This unit
// currently implements stage 4a foundations: initialization, bathymetry face geometry, cell
// /face depth (`update_depth`), and wet-area/volume geometry (`Vs/Asz/Asx/Asy/...`). It
// owns no PETSc types — the linear solve goes through the SparseSystem/LinearSolver
// interface (wired in the matrix-assembly gate). SWE local loops stay plain serial/host;
// on-node Kokkos/GPU is deferred to the global P9/P10 passes (DP3).
#ifndef FREHG2_SWE_SWE_SOLVER_HPP
#define FREHG2_SWE_SWE_SOLVER_HPP

#include <memory>

#include "core/Grid.hpp"
#include "core/MpiComm.hpp"
#include "frehg2/core/define.hpp"
#include "swe/SweFields.hpp"

namespace frehg2 {

class Config;
class LinearSolver;
class SparseSystem;
class Decomp2D;
class PerfRecorder;

struct SweParams {
  real gravity = 9.81;
  real manning = 0.0;
  real min_depth = 1.0e-8;
  real visc_x = 0.0;
  real visc_y = 0.0;
  real hD = 0.1;       // thin-layer depth threshold (drag exponent switch)
  real hE = 0.0;       // evaporation depth threshold (legacy hE)
  real wtfh = 1.0e-8;  // waterfall depth threshold
  real dt = 1.0;
  real offset = 0.0;   // vertical datum offset added to init_eta (legacy offset[0])
};

class SweSolver {
 public:
  // `local_grid` spans this rank's owned interior (localNx x localNy) plus halo. When
  // `mc` is null (unit tests without MPI), the grid is treated as the full global domain.
  explicit SweSolver(const Grid& local_grid, const MpiComm* mc = nullptr);
  ~SweSolver();  // out-of-line: unique_ptr<SparseSystem> (incomplete in this header)

  // Read SW parameters + bathymetry + initial eta from config, then set up the full
  // initial state (eta, depth, face geometry). Bathymetry source:
  //   domain.bathymetry.constant (real)  OR  a prior setBathymetry()/setBathymetryConstant().
  void initialize(const Config& config);

  // Lower-level setup used by tests and by initialize().
  void setParams(const SweParams& p) { params_ = p; }
  void setBathymetryConstant(real z);
  // physical bathymetry, length nSurfaceCell, row-major idx = i + j*nx.
  void setBathymetry(const RealArr1DHost& physical);
  // Set eta = init_eta + offset everywhere (clamped to bottom), u=v=0, then refresh
  // depth and geometry. Bathymetry must already be set.
  void initializeState(real init_eta);

  // P14 flexible IC: write owned-cell physical eta / velocities, then finalize depth+geometry.
  void setInitialEtaPhysical(const RealArr1DHost& physical_eta);
  void setInitialVelocityConstant(real u, real v);
  void setInitialVelocityPhysical(const RealArr1DHost& u, const RealArr1DHost& v);
  void finalizeInitialState();

  // Attach the backend-agnostic linear solver + 2D decomposition (one-rank during 4a-4d,
  // multi-rank from 4e). The decomposition must describe the same nx,ny as `grid`.
  void attachSolver(LinearSolver& solver, Decomp2D& dd);

  // P21 performance instrumentation (optional; Orchestrator attaches its recorder).
  void setPerfRecorder(PerfRecorder* rec) { perf_ = rec; }

  // Advance one semi-implicit SWE step (legacy solve_shallowwater + shallowwater_velocity +
  // update_subgrid_variable). Requires attachSolver().
  void advanceStep(real rain_rate, real evap_rate);

  // Legacy-exact geometry/depth refresh (callable each step).
  void updateBottomFaces();  // boundary_bath: bottom ghosts + bottomXP/bottomYP
  void updateDepth();        // update_depth: dept, deptx, depty
  void updateGeometry();     // update_subgrid_variable (use_subgrid==0): Vs/Asz/Asx/Asy/...

  const Grid& grid() const { return grid_; }
  const SweParams& params() const { return params_; }
  SweFields& fields() { return fields_; }
  const SweFields& fields() const { return fields_; }

  // CFL diagnostic: max( (|u|+sqrt(g*h)) * dt / min(dx,dy) ) over wet cells (P4.3.3).
  real maxCfl() const;

  // Gather owned `dept` onto rank 0 in global row-major order (gi + gj*globalNx).
  // Requires `mc` and a prior attachSolver() call. No-op when not parallel.
  void gatherDeptGlobal(RealArr1DHost& global_dept) const;

 private:
  // Halo-padded surface index (local i in [-1,localNx], j in [-1,localNy]).
  int s(int i, int j) const { return grid_.getSurfaceIndex(i, j); }

  int globalNx() const { return mc_ ? mc_->globalNx() : grid_.nx(); }
  int globalNy() const { return mc_ ? mc_->globalNy() : grid_.ny(); }
  int globalI(int li) const { return mc_ ? mc_->i0() + li : li; }
  int globalJ(int lj) const { return mc_ ? mc_->j0() + lj : lj; }
  bool parallel() const { return mc_ != nullptr && mc_->size() > 1; }

  // Fill ghost cells by zero-gradient mirror at domain outer edges only (MPI partition
  // edges are filled by haloExchange).
  void fillDomainEdgeGhosts(RealArr1DHost& f) const;
  // Serial path: zero-gradient on all four domain boundaries.
  void fillGhostsZeroGradient(RealArr1DHost& f) const;
  void exchangeHalo(RealArr1DHost& f) const;
  void exchangeStepFieldsPreSolve();
  void exchangeMomentumCoeffs();
  void exchangeDepthFields();
  void exchangeSubgridGeometry();
  void exchangeStepFieldsPostVelocity();

  // Legacy-exact substeps (b1-sw path, difuwave==0).
  void enforceSurfBc();      // enforce_surf_bc: clamp eta>=bottom + zero-gradient eta ghosts
  void momentumSource();     // momentum_source: Ex,Ey,Dx,Dy
  void shallowwaterRhs();    // shallowwater_rhs: Srhs
  void shallowwaterMatCoeff();  // shallowwater_mat_coeff: Sxp/Sxm/Syp/Sym/Sct + BC/singularity
  void assembleAndSolve();   // build_shallowwater_system + solve -> eta (via SparseSystem)
  void cflLimiter();         // cfl_limiter: 1-cell wetting + small-depth removal
  void evaprain(real rain_rate, real evap_rate);  // evaprain: rain/evap + dept refresh
  void updateDragCoef();     // update_drag_coef: CDx,CDy
  void waterfallLocation();  // waterfall_location: wtfx,wtfy
  void updateVelocity();     // update_velocity: uu,vv + limiters + Fu,Fv,cflx,cfly
  void enforceVeloBc();      // enforce_velo_bc (zero-gradient ghosts; tidal disabled)
  void interpVelocity();     // interp_velocity: uy,vx

  Grid grid_;
  SweParams params_;
  SweFields fields_;

  const MpiComm* mc_ = nullptr;

  // Linear-solver wiring (no PETSc types here; seam holds).
  LinearSolver* solver_ = nullptr;
  Decomp2D* dd_ = nullptr;
  std::unique_ptr<SparseSystem> sys_;
  RealArr1D b_owned_;
  RealArr1D x_owned_;
  PerfRecorder* perf_ = nullptr;
};

}  // namespace frehg2

#endif  // FREHG2_SWE_SWE_SOLVER_HPP
