#include "solute/Diffusion.hpp"

#include <Kokkos_Core.hpp>
#include <tuple>

#include "frehg2/core/ParallelFor.hpp"
#include "frehg2/linear/LinearSolver.hpp"
#include "frehg2/linear/SparseSystem.hpp"
#include "linear/DomainDecomposition.hpp"

namespace frehg2 {

DiffusionSolver::DiffusionSolver(LinearSolver& solver, Decomp2D& dd, const Grid& grid)
    : solver_(&solver), dd2_(&dd), grid_(grid) {
  sys_ = solver.createSystem(dd);
  const int nloc = dd.ownedRowCount();
  b_owned_ = RealArr1D("solute_diff_b", static_cast<size_t>(nloc));
  x_owned_ = RealArr1D("solute_diff_x", static_cast<size_t>(nloc));
}

DiffusionSolver::DiffusionSolver(LinearSolver& solver, Decomp3D& dd, const Grid& grid)
    : solver_(&solver), dd3_(&dd), grid_(grid) {
  sys_ = solver.createSystem(dd);
  const int nloc = dd.ownedRowCount();
  b_owned_ = RealArr1D("solute_diff_b", static_cast<size_t>(nloc));
  x_owned_ = RealArr1D("solute_diff_x", static_cast<size_t>(nloc));
}

DiffusionSolver::~DiffusionSolver() = default;

void DiffusionSolver::solve(RealArr1DHost& conc, real D, real dt) {
  if (D <= 0.0 || dt <= 0.0) return;
  if (dd2_ != nullptr) {
    solve2D(conc, D, dt);
  } else {
    solve3D(conc, D, dt);
  }
}

void DiffusionSolver::solve2D(RealArr1DHost& conc, real D, real dt) {
  const int gnx = grid_.nx();
  const int gny = grid_.ny();
  const real cx = D * dt / (grid_.dx() * grid_.dx());
  const real cy = D * dt / (grid_.dy() * grid_.dy());
  const auto corners = dd2_->localCorners();
  const int i0 = std::get<0>(corners), j0 = std::get<1>(corners);
  const int ni = std::get<2>(corners), nj = std::get<3>(corners);
  const int rstart = dd2_->ownershipRange().first;

  sys_->zeroEntries();
  sys_->beginAssembly();
  auto bh = Kokkos::create_mirror_view(b_owned_);
  for (int gj = j0; gj < j0 + nj; ++gj)
    for (int gi = i0; gi < i0 + ni; ++gi) {
      const int grow = dd2_->globalRow(gi, gj);
      int cols[kMaxStencil];
      int ncols = 0;
      dd2_->stencilColumns(gi, gj, cols, ncols);
      const int rowW = (gi - 1 >= 0) ? dd2_->globalRow(gi - 1, gj) : -1;
      const int rowE = (gi + 1 < gnx) ? dd2_->globalRow(gi + 1, gj) : -1;
      const int rowS = (gj - 1 >= 0) ? dd2_->globalRow(gi, gj - 1) : -1;
      const int rowN = (gj + 1 < gny) ? dd2_->globalRow(gi, gj + 1) : -1;
      real diag = 1.0;
      if (rowW >= 0) diag += cx;
      if (rowE >= 0) diag += cx;
      if (rowS >= 0) diag += cy;
      if (rowN >= 0) diag += cy;
      real vals[kMaxStencil];
      for (int cc = 0; cc < ncols; ++cc) {
        const int col = cols[cc];
        if (col == grow)
          vals[cc] = diag;
        else if (col == rowW || col == rowE)
          vals[cc] = -cx;
        else if (col == rowS || col == rowN)
          vals[cc] = -cy;
        else
          vals[cc] = 0.0;
      }
      sys_->addRow(grow, ncols, cols, vals);
      bh(static_cast<size_t>(grow - rstart)) =
          conc(static_cast<size_t>(grid_.getSurfaceIndex(gi, gj)));
    }
  Kokkos::deep_copy(b_owned_, bh);
  sys_->endAssembly();

  solver_->setup(*sys_);
  solver_->solve(*sys_, b_owned_, x_owned_);

  auto xh = Kokkos::create_mirror_view(x_owned_);
  Kokkos::deep_copy(xh, x_owned_);
  // Solution unpack: per-cell, disjoint writes -> Kokkos parallel-for (P9).
  const Grid g = grid_;
  Decomp2D* dd = dd2_;
  RealArr1DHost cc = conc;
  const int rs = rstart;
  parallelForSurface("solute_diff2d_unpack", ni, nj, KOKKOS_LAMBDA(int li, int lj) {
    const int gi = i0 + li;
    const int gj = j0 + lj;
    const int grow = dd->globalRow(gi, gj);
    cc(static_cast<size_t>(g.getSurfaceIndex(gi, gj))) = xh(static_cast<size_t>(grow - rs));
  });
}

void DiffusionSolver::solve3D(RealArr1DHost& conc, real D, real dt) {
  const int gnx = grid_.nx();
  const int gny = grid_.ny();
  const int gnz = grid_.nz();
  const real cx = D * dt / (grid_.dx() * grid_.dx());
  const real cy = D * dt / (grid_.dy() * grid_.dy());
  const real cz = D * dt / (grid_.dz() * grid_.dz());
  const auto corners = dd3_->localCorners();
  const int i0 = std::get<0>(corners), j0 = std::get<1>(corners);
  const int k0 = std::get<2>(corners);
  const int ni = std::get<3>(corners), nj = std::get<4>(corners), nk = std::get<5>(corners);
  const int rstart = dd3_->ownershipRange().first;

  sys_->zeroEntries();
  sys_->beginAssembly();
  auto bh = Kokkos::create_mirror_view(b_owned_);
  for (int gk = k0; gk < k0 + nk; ++gk)
    for (int gj = j0; gj < j0 + nj; ++gj)
      for (int gi = i0; gi < i0 + ni; ++gi) {
        const int grow = dd3_->globalRow(gi, gj, gk);
        int cols[kMaxStencil];
        int ncols = 0;
        dd3_->stencilColumns(gi, gj, gk, cols, ncols);
        const int rowW = (gi - 1 >= 0) ? dd3_->globalRow(gi - 1, gj, gk) : -1;
        const int rowE = (gi + 1 < gnx) ? dd3_->globalRow(gi + 1, gj, gk) : -1;
        const int rowS = (gj - 1 >= 0) ? dd3_->globalRow(gi, gj - 1, gk) : -1;
        const int rowN = (gj + 1 < gny) ? dd3_->globalRow(gi, gj + 1, gk) : -1;
        const int rowB = (gk - 1 >= 0) ? dd3_->globalRow(gi, gj, gk - 1) : -1;
        const int rowT = (gk + 1 < gnz) ? dd3_->globalRow(gi, gj, gk + 1) : -1;
        real diag = 1.0;
        if (rowW >= 0) diag += cx;
        if (rowE >= 0) diag += cx;
        if (rowS >= 0) diag += cy;
        if (rowN >= 0) diag += cy;
        if (rowB >= 0) diag += cz;
        if (rowT >= 0) diag += cz;
        real vals[kMaxStencil];
        for (int cc = 0; cc < ncols; ++cc) {
          const int col = cols[cc];
          if (col == grow)
            vals[cc] = diag;
          else if (col == rowW || col == rowE)
            vals[cc] = -cx;
          else if (col == rowS || col == rowN)
            vals[cc] = -cy;
          else if (col == rowB || col == rowT)
            vals[cc] = -cz;
          else
            vals[cc] = 0.0;
        }
        sys_->addRow(grow, ncols, cols, vals);
        bh(static_cast<size_t>(grow - rstart)) =
            conc(static_cast<size_t>(grid_.getIndex(gi, gj, gk)));
      }
  Kokkos::deep_copy(b_owned_, bh);
  sys_->endAssembly();

  solver_->setup(*sys_);
  solver_->solve(*sys_, b_owned_, x_owned_);

  auto xh = Kokkos::create_mirror_view(x_owned_);
  Kokkos::deep_copy(xh, x_owned_);
  // Solution unpack: per-cell, disjoint writes -> Kokkos parallel-for (P9).
  const Grid g = grid_;
  Decomp3D* dd = dd3_;
  RealArr1DHost cc = conc;
  const int rs = rstart;
  parallelForVolume("solute_diff3d_unpack", ni, nj, nk, KOKKOS_LAMBDA(int li, int lj, int lk) {
    const int gi = i0 + li;
    const int gj = j0 + lj;
    const int gk = k0 + lk;
    const int grow = dd->globalRow(gi, gj, gk);
    cc(static_cast<size_t>(g.getIndex(gi, gj, gk))) = xh(static_cast<size_t>(grow - rs));
  });
}

}  // namespace frehg2
