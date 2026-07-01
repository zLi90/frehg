// Implicit, operator-split diffusion for the passive scalar (P8.3.3).
//
// Solves (I - dt * D * L) C^{n+1} = C^n where L is the standard 5-point (surface) / 7-point
// (subsurface) Neumann Laplacian on the uniform grid. The matrix is SPD and is assembled and
// solved through the backend-agnostic SparseSystem / LinearSolver interface (no PETSc types
// here; the seam check covers src/solute). With zero-flux (Neumann) outer boundaries the
// operator has zero row sums, so the implicit update conserves the cell-sum of concentration.
#ifndef FREHG2_SOLUTE_DIFFUSION_HPP
#define FREHG2_SOLUTE_DIFFUSION_HPP

#include <memory>

#include "core/Grid.hpp"
#include "frehg2/core/define.hpp"

namespace frehg2 {

class LinearSolver;
class SparseSystem;
class Decomp2D;
class Decomp3D;

class DiffusionSolver {
 public:
  // Surface (2D, 5-point) diffusion.
  DiffusionSolver(LinearSolver& solver, Decomp2D& dd, const Grid& grid);
  // Subsurface (3D, 7-point) diffusion.
  DiffusionSolver(LinearSolver& solver, Decomp3D& dd, const Grid& grid);
  ~DiffusionSolver();

  // In-place implicit diffusion of a halo-padded host concentration field over one step.
  // No-op when D <= 0 or dt <= 0.
  void solve(RealArr1DHost& conc, real D, real dt);

 private:
  void solve2D(RealArr1DHost& conc, real D, real dt);
  void solve3D(RealArr1DHost& conc, real D, real dt);

  LinearSolver* solver_ = nullptr;
  Decomp2D* dd2_ = nullptr;
  Decomp3D* dd3_ = nullptr;
  Grid grid_;
  std::unique_ptr<SparseSystem> sys_;
  RealArr1D b_owned_;
  RealArr1D x_owned_;
};

}  // namespace frehg2

#endif  // FREHG2_SOLUTE_DIFFUSION_HPP
