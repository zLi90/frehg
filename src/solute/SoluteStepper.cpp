#include "solute/SoluteStepper.hpp"

#include <Kokkos_Core.hpp>
#include <algorithm>

#include "frehg2/core/ParallelFor.hpp"
#include "solute/Advection.hpp"
#include "solute/SourceSink.hpp"

namespace frehg2 {

SoluteStepper::SoluteStepper(const Grid& grid, const SoluteParams& p)
    : grid_(grid), p_(p) {}

SoluteStepper::~SoluteStepper() = default;

void SoluteStepper::attachSurfaceDiffusion(LinearSolver& solver, Decomp2D& dd) {
  surf_diff_ = std::make_unique<DiffusionSolver>(solver, dd, grid_);
}

void SoluteStepper::attachSubsurfaceDiffusion(LinearSolver& solver, Decomp3D& dd) {
  subs_diff_ = std::make_unique<DiffusionSolver>(solver, dd, grid_);
}

StepResult SoluteStepper::step(State& sw, GwState& gw, const SoluteFlow& flow, real dt,
                               real rain) {
  StepResult r;
  if (!p_.enabled) return r;  // no-op: no transport, no PETSc assembly

  const bool do_surf = flow.hasSurface();
  const bool do_subs = flow.hasSubsurface();

  // Bridge canonical (device) state -> host work buffers.
  RealArr1DHost cs("solute_cs", sw.conc.extent(0));
  RealArr1DHost cg("solute_cg", gw.conc.extent(0));
  Kokkos::deep_copy(cs, sw.conc);
  Kokkos::deep_copy(cg, gw.conc);

  // (1) Source / sink: rainfall -> decay. The SW<->GW interface solute exchange is NOT done
  // here: it is co-located with the coupling's water exchange in the Orchestrator
  // (applyInterfaceExchange) so the post-exchange water heights are exact and the transfer is
  // mass-conservative. The stepper only does bulk in-domain transport (rain, decay, advection,
  // diffusion).
  if (do_surf) applyRainfall(cs, flow.depth, grid_, rain, p_, dt);
  if (do_surf) applyDecaySurface(cs, grid_, p_.k_decay, dt);
  if (do_subs) applyDecaySubsurface(cg, grid_, p_.k_decay, dt);

  // (2) Advection (explicit). CFL refusal leaves the canonical state untouched.
  real cfl = 0.0;
  if (do_surf) cfl = std::max(cfl, advectSurface(cs, flow, grid_, p_, dt));
  if (do_subs) cfl = std::max(cfl, advectSubsurface(cg, flow, grid_, p_, dt));
  r.max_cfl = cfl;
  if (cfl > p_.cfl_max) {
    r.ok = false;
    return r;  // discard host buffers; State / GwState unchanged
  }

  // (3) Diffusion (implicit) through the LinearSolver seam.
  if (p_.diffusion_scheme == "implicit" && p_.D > 0.0) {
    if (do_surf && surf_diff_) surf_diff_->solve(cs, p_.D, dt);
    if (do_subs && subs_diff_) subs_diff_->solve(cg, p_.D, dt);
  }

  // (4) Record: commit host buffers back to the canonical (device) fields.
  if (do_surf) Kokkos::deep_copy(sw.conc, cs);
  if (do_subs) Kokkos::deep_copy(gw.conc, cg);

  real dmin = grid_.dx();
  dmin = std::min(dmin, grid_.dy());
  if (do_subs) dmin = std::min(dmin, grid_.dz());
  r.max_diffusion_cfl = p_.D * dt / (dmin * dmin);
  return r;
}

}  // namespace frehg2
