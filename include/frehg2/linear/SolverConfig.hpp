// Backend-neutral linear-solver configuration (P2.5.2).
//
// Populated from YAML / PETSc options; carries NO library-specific calls. The strings map
// onto PETSc KSP/PC type names today and onto an equivalent for any future backend.
#ifndef FREHG2_LINEAR_SOLVER_CONFIG_HPP
#define FREHG2_LINEAR_SOLVER_CONFIG_HPP

#include <string>

#include "frehg2/core/define.hpp"

namespace frehg2 {

struct SolverConfig {
  std::string ksp_type = "cg";    // Krylov method (cg for SPD, gmres for non-symmetric)
  std::string pc_type = "gamg";   // preconditioner (algebraic multigrid by default)
  real rtol = 1e-8;               // relative residual tolerance
  real atol = 1e-50;              // absolute residual tolerance (PETSc default)
  int max_it = 10000;             // maximum Krylov iterations
};

}  // namespace frehg2

#endif  // FREHG2_LINEAR_SOLVER_CONFIG_HPP
