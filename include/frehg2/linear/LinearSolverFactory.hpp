// Backend-agnostic LinearSolver factory (P7).
//
// The production driver (Orchestrator) must create the concrete linear-solver backend
// without ever including a solver-library header. This factory returns a
// std::unique_ptr<LinearSolver> for the requested backend name; the only sanctioned
// backend today is "petsc" (PetscLinearSolver, the default). The implementation lives in
// src/linear/backends/ — the single place permitted to include <petscksp.h> — so callers
// (driver, physics) stay free of any PETSc dependency.
#ifndef FREHG2_LINEAR_LINEAR_SOLVER_FACTORY_HPP
#define FREHG2_LINEAR_LINEAR_SOLVER_FACTORY_HPP

#include <memory>
#include <string>

#include "frehg2/linear/LinearSolver.hpp"
#include "frehg2/linear/SolverConfig.hpp"

namespace frehg2 {

// Construct the linear-solver backend selected by `backend` (default "petsc"). Throws
// std::runtime_error for an unknown / not-built backend so misconfiguration fails loud.
std::unique_ptr<LinearSolver> makeLinearSolver(const std::string& backend,
                                               const SolverConfig& cfg);

}  // namespace frehg2

#endif  // FREHG2_LINEAR_LINEAR_SOLVER_FACTORY_HPP
