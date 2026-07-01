// Standalone solute step driver (P8.3.5).
//
// Drives one operator-split solute transport step on the surface (State) and/or subsurface
// (GwState) concentration fields, in the legacy order: source/sink -> advect -> diffuse ->
// record. This is the SOLVER ONLY; it is invoked from unit-test drivers in P8 and is wired
// into the Orchestrator in P16 (this header intentionally does NOT depend on Orchestrator).
//
// State.conc / GwState.conc are the canonical (device) concentration fields. The stepper
// bridges them to host work buffers, runs the plain-loop kernels (Kokkos-ification is the
// global P9 pass), and writes the result back only when the step is accepted, so a CFL
// refusal leaves the input state untouched.
#ifndef FREHG2_SOLUTE_SOLUTE_STEPPER_HPP
#define FREHG2_SOLUTE_SOLUTE_STEPPER_HPP

#include <memory>

#include "core/GwState.hpp"
#include "core/Grid.hpp"
#include "core/State.hpp"
#include "frehg2/core/define.hpp"
#include "solute/Diffusion.hpp"
#include "solute/SoluteFlow.hpp"
#include "solute/SoluteParams.hpp"

namespace frehg2 {

class LinearSolver;
class Decomp2D;
class Decomp3D;

struct StepResult {
  bool ok = true;             // false if the advective CFL exceeded cfl_max (step refused)
  real max_cfl = 0.0;          // max advective Courant number this step
  real max_diffusion_cfl = 0.0;  // diffusion number D*dt/dx_min^2 (informational; implicit)
};

class SoluteStepper {
 public:
  SoluteStepper(const Grid& grid, const SoluteParams& p);
  ~SoluteStepper();

  // Attach the implicit-diffusion backends (optional; needed only when diffusion_scheme is
  // "implicit" and D > 0). Surface uses Decomp2D, subsurface uses Decomp3D.
  void attachSurfaceDiffusion(LinearSolver& solver, Decomp2D& dd);
  void attachSubsurfaceDiffusion(LinearSolver& solver, Decomp3D& dd);

  const SoluteParams& params() const { return p_; }

  // Advance one step. `rain` is the current surface rainfall rate [m/s] (0 if none). The
  // domains driven are determined by which SoluteFlow fields are populated.
  StepResult step(State& sw, GwState& gw, const SoluteFlow& flow, real dt, real rain = 0.0);

 private:
  Grid grid_;
  SoluteParams p_;
  std::unique_ptr<DiffusionSolver> surf_diff_;
  std::unique_ptr<DiffusionSolver> subs_diff_;
};

}  // namespace frehg2

#endif  // FREHG2_SOLUTE_SOLUTE_STEPPER_HPP
