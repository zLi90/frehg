// Solute-transport parameters (P8.3.6).
//
// Backend-neutral configuration for the passive-scalar advection-diffusion solver. Parsed
// from the YAML `solute:` block. When `enabled` is false the whole P8 code path is a no-op
// (P16 turns it on through the Orchestrator).
#ifndef FREHG2_SOLUTE_SOLUTE_PARAMS_HPP
#define FREHG2_SOLUTE_SOLUTE_PARAMS_HPP

#include <string>

#include "frehg2/core/define.hpp"

namespace frehg2 {

class Config;

struct SoluteParams {
  bool enabled = false;             // solute.enabled  (false => no-op, no PETSc assembly)
  real c_rain = 0.0;                // solute.c_rain   rainfall concentration
  real k_decay = 0.0;              // solute.k_decay  first-order decay constant [1/T]
  real D = 1.0e-9;                  // solute.D        isotropic diffusion/dispersion [m^2/s]
  std::string advection_scheme = "upwind";    // "upwind" (1st order) | "muscl" (2nd order)
  std::string diffusion_scheme = "implicit";  // "implicit" (PETSc) | "none"
  real cfl_max = 1.0;              // advective CFL ceiling; step refuses above this
  real min_depth = 1.0e-8;         // surface wet threshold for sources / advection

  // Parse from a Config (the `solute:` block); missing keys take the defaults above.
  static SoluteParams fromConfig(const Config& cfg);
};

}  // namespace frehg2

#endif  // FREHG2_SOLUTE_SOLUTE_PARAMS_HPP
