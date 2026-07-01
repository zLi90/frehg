#include "frehg2/linear/LinearSolverFactory.hpp"

#include <stdexcept>

#include "linear/backends/PetscLinearSolver.hpp"

namespace frehg2 {

std::unique_ptr<LinearSolver> makeLinearSolver(const std::string& backend,
                                               const SolverConfig& cfg) {
  if (backend.empty() || backend == "petsc") {
    return std::make_unique<PetscLinearSolver>(cfg);
  }
  throw std::runtime_error("makeLinearSolver: unknown / not-built linear-solver backend '" +
                           backend + "' (only 'petsc' is sanctioned)");
}

}  // namespace frehg2
