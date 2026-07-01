#include "swe/SweSolver.hpp"

#include <Kokkos_Core.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <mpi.h>

#include "frehg2/core/ParallelFor.hpp"
#include "frehg2/linear/LinearSolver.hpp"
#include "frehg2/linear/SparseSystem.hpp"
#include "frehg2/perf/PerfRecorder.hpp"
#include "frehg2/perf/Timer.hpp"
#include "io/Config.hpp"
#include "linear/DomainDecomposition.hpp"

namespace frehg2 {

SweSolver::SweSolver(const Grid& local_grid, const MpiComm* mc)
    : grid_(local_grid), fields_(local_grid), mc_(mc) {}

SweSolver::~SweSolver() = default;

void SweSolver::fillGhostsZeroGradient(RealArr1DHost& f) const {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  for (int i = 0; i < nx; ++i) {
    f(s(i, -1)) = f(s(i, 0));
    f(s(i, ny)) = f(s(i, ny - 1));
  }
  for (int j = 0; j < ny; ++j) {
    f(s(-1, j)) = f(s(0, j));
    f(s(nx, j)) = f(s(nx - 1, j));
  }
}

void SweSolver::fillDomainEdgeGhosts(RealArr1DHost& f) const {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  if (!parallel()) {
    fillGhostsZeroGradient(f);
    return;
  }
  if (mc_->leftRank() == MPI_PROC_NULL) {
    for (int j = 0; j < ny; ++j) f(s(-1, j)) = f(s(0, j));
  }
  if (mc_->rightRank() == MPI_PROC_NULL) {
    for (int j = 0; j < ny; ++j) f(s(nx, j)) = f(s(nx - 1, j));
  }
  if (mc_->downRank() == MPI_PROC_NULL) {
    for (int i = 0; i < nx; ++i) f(s(i, -1)) = f(s(i, 0));
  }
  if (mc_->upRank() == MPI_PROC_NULL) {
    for (int i = 0; i < nx; ++i) f(s(i, ny)) = f(s(i, ny - 1));
  }
}

void SweSolver::exchangeHalo(RealArr1DHost& f) const {
  if (!parallel() || dd_ == nullptr) return;
  RealArr1D dev("swe_halo", f.extent(0));
  Kokkos::deep_copy(dev, f);
  dd_->haloExchange(dev);
  Kokkos::deep_copy(f, dev);
}

void SweSolver::exchangeStepFieldsPreSolve() {
  if (!parallel()) return;
  auto& f = fields_;
  exchangeHalo(f.eta);
  exchangeHalo(f.dept);
  exchangeHalo(f.uu);
  exchangeHalo(f.vv);
  exchangeHalo(f.uy);
  exchangeHalo(f.vx);
  exchangeHalo(f.Vs);
  exchangeHalo(f.Asx);
  exchangeHalo(f.Asy);
  exchangeHalo(f.Asz);
  exchangeHalo(f.Vsx);
  exchangeHalo(f.Vsy);
  exchangeHalo(f.Aszx);
  exchangeHalo(f.Aszy);
  exchangeHalo(f.Ex);
  exchangeHalo(f.Ey);
  exchangeHalo(f.Dx);
  exchangeHalo(f.Dy);
  exchangeHalo(f.CDx);
  exchangeHalo(f.CDy);
}

void SweSolver::exchangeMomentumCoeffs() {
  if (!parallel()) return;
  auto& f = fields_;
  exchangeHalo(f.Ex);
  exchangeHalo(f.Ey);
  exchangeHalo(f.Dx);
  exchangeHalo(f.Dy);
}

void SweSolver::exchangeDepthFields() {
  if (!parallel()) return;
  auto& f = fields_;
  exchangeHalo(f.dept);
  exchangeHalo(f.deptx);
  exchangeHalo(f.depty);
}

void SweSolver::exchangeSubgridGeometry() {
  if (!parallel()) return;
  auto& f = fields_;
  exchangeHalo(f.Vs);
  exchangeHalo(f.Vsx);
  exchangeHalo(f.Vsy);
  exchangeHalo(f.Asx);
  exchangeHalo(f.Asy);
  exchangeHalo(f.Asz);
  exchangeHalo(f.Aszx);
  exchangeHalo(f.Aszy);
}

void SweSolver::exchangeStepFieldsPostVelocity() {
  if (!parallel()) return;
  auto& f = fields_;
  exchangeHalo(f.uu);
  exchangeHalo(f.vv);
  exchangeHalo(f.Vs);
  exchangeHalo(f.Asx);
  exchangeHalo(f.Asy);
  exchangeHalo(f.Asz);
  exchangeHalo(f.Vsx);
  exchangeHalo(f.Vsy);
  exchangeHalo(f.Aszx);
  exchangeHalo(f.Aszy);
}

void SweSolver::setBathymetryConstant(real z) {
  const int n = grid_.nSurfaceStorageCell();
  RealArr1DHost bottom = fields_.bottom;
  parallelForRange("swe_set_bathy_const", n, KOKKOS_LAMBDA(int idx) { bottom(idx) = z; });
  updateBottomFaces();
}

void SweSolver::setBathymetry(const RealArr1DHost& physical) {
  const int lnx = grid_.nx();
  const int lny = grid_.ny();
  const int gnx = globalNx();
  for (int j = 0; j < lny; ++j) {
    for (int i = 0; i < lnx; ++i) {
      const int gi = globalI(i);
      const int gj = globalJ(j);
      fields_.bottom(s(i, j)) = physical(static_cast<size_t>(gi + gj * gnx));
    }
  }
  updateBottomFaces();
  if (parallel()) exchangeHalo(fields_.bottom);
}

void SweSolver::updateBottomFaces() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  auto& bottom = fields_.bottom;
  auto& bXP = fields_.bottomXP;
  auto& bYP = fields_.bottomYP;

  fillDomainEdgeGhosts(bottom);
  if (parallel()) exchangeHalo(bottom);
  {
    const Grid g = grid_;
    RealArr1DHost bot = bottom;
    RealArr1DHost bxp = bXP;
    RealArr1DHost byp = bYP;
    parallelForSurface("swe_bottom_faces", nx, ny, KOKKOS_LAMBDA(int i, int j) {
      const int c = g.getSurfaceIndex(i, j);
      bxp(c) = bot(c);
      if (bot(g.getSurfaceIndex(i + 1, j)) > bot(c)) bxp(c) = bot(g.getSurfaceIndex(i + 1, j));
      byp(c) = bot(c);
      if (bot(g.getSurfaceIndex(i, j + 1)) > bot(c)) byp(c) = bot(g.getSurfaceIndex(i, j + 1));
    });
  }
  // Face-bed ghosts at domain outer edges only (partition edges use haloExchange).
  if (!parallel() || mc_->downRank() == MPI_PROC_NULL) {
    for (int i = 0; i < nx; ++i) {
      bYP(s(i, -1)) = bottom(s(i, 0));
      if (bottom(s(i, -1)) > bottom(s(i, 0))) bYP(s(i, -1)) = bottom(s(i, -1));
    }
  }
  if (!parallel() || mc_->upRank() == MPI_PROC_NULL) {
    for (int i = 0; i < nx; ++i) bYP(s(i, ny)) = bYP(s(i, ny - 1));
  }
  if (!parallel() || mc_->leftRank() == MPI_PROC_NULL) {
    for (int j = 0; j < ny; ++j) {
      bXP(s(-1, j)) = bottom(s(0, j));
      if (bottom(s(-1, j)) > bottom(s(0, j))) bXP(s(-1, j)) = bottom(s(-1, j));
    }
  }
  if (!parallel() || mc_->rightRank() == MPI_PROC_NULL) {
    for (int j = 0; j < ny; ++j) bXP(s(nx, j)) = bXP(s(nx - 1, j));
  }
  if (parallel()) {
    exchangeHalo(bottom);
    exchangeHalo(bXP);
    exchangeHalo(bYP);
  }
}

void SweSolver::updateDepth() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const real min_dept = params_.min_depth;
  auto& eta = fields_.eta;
  auto& bottom = fields_.bottom;
  auto& dept = fields_.dept;
  auto& deptx = fields_.deptx;
  auto& depty = fields_.depty;

  const Grid g = grid_;
  // Snap tiny interior depth to the bed (update_depth lines 939-944).
  {
    RealArr1DHost e = eta;
    RealArr1DHost b = bottom;
    parallelForSurface("swe_depth_snap", nx, ny, KOKKOS_LAMBDA(int i, int j) {
      const int c = g.getSurfaceIndex(i, j);
      const real diff = e(c) - b(c);
      if (diff > 0.0 && diff < min_dept) e(c) = b(c);
    });
  }

  // Center depth over ALL cells incl. ghosts (uses eta ghosts as last set by enforceSurfBc).
  const int n = grid_.nSurfaceStorageCell();
  {
    RealArr1DHost e = eta;
    RealArr1DHost b = bottom;
    RealArr1DHost dp = dept;
    RealArr1DHost dpx = deptx;
    RealArr1DHost dpy = depty;
    parallelForRange("swe_depth_center", n, KOKKOS_LAMBDA(int idx) {
      dp(idx) = e(idx) - b(idx);
      if (dp(idx) <= min_dept) dp(idx) = 0.0;
      dpx(idx) = 0.0;
      dpy(idx) = 0.0;
    });
  }

  {
    RealArr1DHost e = eta;
    RealArr1DHost b = bottom;
    RealArr1DHost dpx = deptx;
    RealArr1DHost dpy = depty;
    parallelForSurface("swe_depth_faces", nx, ny, KOKKOS_LAMBDA(int i, int j) {
      const int c = g.getSurfaceIndex(i, j);
      real eta_hi = e(c), bot_hi = b(c);
      if (e(g.getSurfaceIndex(i + 1, j)) > eta_hi) eta_hi = e(g.getSurfaceIndex(i + 1, j));
      if (b(g.getSurfaceIndex(i + 1, j)) > bot_hi) bot_hi = b(g.getSurfaceIndex(i + 1, j));
      dpx(c) = eta_hi - bot_hi;
      eta_hi = e(c);
      bot_hi = b(c);
      if (e(g.getSurfaceIndex(i, j + 1)) > eta_hi) eta_hi = e(g.getSurfaceIndex(i, j + 1));
      if (b(g.getSurfaceIndex(i, j + 1)) > bot_hi) bot_hi = b(g.getSurfaceIndex(i, j + 1));
      dpy(c) = eta_hi - bot_hi;
    });
  }

  if (!parallel() || mc_->downRank() == MPI_PROC_NULL) {
    for (int i = 0; i < nx; ++i) {
      real eta_hi = eta(s(i, 0)), bot_hi = bottom(s(i, 0));
      if (eta(s(i, -1)) > eta_hi) eta_hi = eta(s(i, -1));
      if (bottom(s(i, -1)) > bot_hi) bot_hi = bottom(s(i, -1));
      depty(s(i, -1)) = eta_hi - bot_hi;
    }
  }
  if (!parallel() || mc_->upRank() == MPI_PROC_NULL) {
    for (int i = 0; i < nx; ++i) {
      real eta_hi = eta(s(i, ny - 1)), bot_hi = bottom(s(i, ny - 1));
      if (eta(s(i, ny)) > eta_hi) eta_hi = eta(s(i, ny));
      if (bottom(s(i, ny)) > bot_hi) bot_hi = bottom(s(i, ny));
      depty(s(i, ny)) = eta_hi - bot_hi;
    }
  }
  if (!parallel() || mc_->leftRank() == MPI_PROC_NULL) {
    for (int j = 0; j < ny; ++j) {
      real eta_hi = eta(s(0, j)), bot_hi = bottom(s(0, j));
      if (eta(s(-1, j)) > eta_hi) eta_hi = eta(s(-1, j));
      if (bottom(s(-1, j)) > bot_hi) bot_hi = bottom(s(-1, j));
      deptx(s(-1, j)) = eta_hi - bot_hi;
    }
  }
  if (!parallel() || mc_->rightRank() == MPI_PROC_NULL) {
    for (int j = 0; j < ny; ++j) {
      real eta_hi = eta(s(nx - 1, j)), bot_hi = bottom(s(nx - 1, j));
      if (eta(s(nx, j)) > eta_hi) eta_hi = eta(s(nx, j));
      if (bottom(s(nx, j)) > bot_hi) bot_hi = bottom(s(nx, j));
      deptx(s(nx, j)) = eta_hi - bot_hi;
    }
  }
  {
    RealArr1DHost dpx = deptx;
    RealArr1DHost dpy = depty;
    parallelForRange("swe_depth_clamp", n, KOKKOS_LAMBDA(int idx) {
      if (dpx(idx) < 0.0) dpx(idx) = 0.0;
      if (dpy(idx) < 0.0) dpy(idx) = 0.0;
    });
  }
}

void SweSolver::updateGeometry() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const real dx = grid_.dx();
  const real dy = grid_.dy();
  auto& dept = fields_.dept;
  auto& deptx = fields_.deptx;
  auto& depty = fields_.depty;
  auto& Vs = fields_.Vs;
  auto& Vsn = fields_.Vsn;
  auto& Asx = fields_.Asx;
  auto& Asy = fields_.Asy;
  auto& Asz = fields_.Asz;
  auto& Aszx = fields_.Aszx;
  auto& Aszy = fields_.Aszy;
  auto& Vsx = fields_.Vsx;
  auto& Vsy = fields_.Vsy;

  const int n = grid_.nSurfaceStorageCell();
  const Grid g = grid_;
  const int gnx1 = globalNx();
  const int gny1 = globalNy();
  {
    RealArr1DHost vsn = Vsn, vs = Vs, asz = Asz, asx = Asx, asy = Asy;
    RealArr1DHost dp = dept, dpx = deptx, dpy = depty;
    parallelForRange("swe_geom_cell", n, KOKKOS_LAMBDA(int idx) {
      vsn(idx) = vs(idx);
      vs(idx) = dp(idx) * dx * dy;
      asz(idx) = dp(idx) > 0.0 ? dx * dy : 0.0;
      asx(idx) = dpx(idx) * dy;
      asy(idx) = dpy(idx) * dx;
      if (gnx1 == 1) asx(idx) = 0.0;
      if (gny1 == 1) asy(idx) = 0.0;
    });
  }
  {
    RealArr1DHost vs = Vs, asz = Asz, vsx = Vsx, vsy = Vsy, aszx = Aszx, aszy = Aszy;
    parallelForSurface("swe_geom_face", nx, ny, KOKKOS_LAMBDA(int i, int j) {
      const int c = g.getSurfaceIndex(i, j);
      vsx(c) = 0.5 * (vs(c) + vs(g.getSurfaceIndex(i + 1, j)));
      vsy(c) = 0.5 * (vs(c) + vs(g.getSurfaceIndex(i, j + 1)));
      aszx(c) = 0.5 * (asz(c) + asz(g.getSurfaceIndex(i + 1, j)));
      aszy(c) = 0.5 * (asz(c) + asz(g.getSurfaceIndex(i, j + 1)));
    });
  }
  // Ghost geometry: domain outer edges only; partition edges come from haloExchange
  // (exchangeStepFieldsPostVelocity in advanceStep).
  fillDomainEdgeGhosts(Vs);
  fillDomainEdgeGhosts(Asx);
  fillDomainEdgeGhosts(Asy);
  fillDomainEdgeGhosts(Asz);
  if (!parallel() || mc_->downRank() == MPI_PROC_NULL) {
    for (int i = 0; i < nx; ++i) {
      Vsy(s(i, -1)) = Vs(s(i, 0));
      Aszy(s(i, -1)) = Asz(s(i, 0));
    }
  }
  if (!parallel() || mc_->leftRank() == MPI_PROC_NULL) {
    for (int j = 0; j < ny; ++j) {
      Vsx(s(-1, j)) = Vs(s(0, j));
      Aszx(s(-1, j)) = Asz(s(0, j));
    }
  }
}

void SweSolver::enforceSurfBc() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  auto& eta = fields_.eta;
  auto& bottom = fields_.bottom;
  const Grid g = grid_;
  RealArr1DHost e = eta;
  RealArr1DHost b = bottom;
  parallelForSurface("swe_enforce_surf_bc", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const int c = g.getSurfaceIndex(i, j);
    if (e(c) < b(c)) e(c) = b(c);
  });
  fillDomainEdgeGhosts(eta);
}

void SweSolver::momentumSource() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const real dx = grid_.dx();
  const real dy = grid_.dy();
  const real dt = params_.dt;
  const real viscx = params_.visc_x;
  const real viscy = params_.visc_y;
  const Grid g = grid_;
  RealArr1DHost uu_ = fields_.uu, vv_ = fields_.vv, vx_ = fields_.vx, uy_ = fields_.uy;
  RealArr1DHost cflx = fields_.cflx, cfly = fields_.cfly;
  RealArr1DHost Vsx = fields_.Vsx, Vsy = fields_.Vsy, Asx = fields_.Asx, Asy = fields_.Asy;
  RealArr1DHost Aszx = fields_.Aszx, Aszy = fields_.Aszy, CDx = fields_.CDx, CDy = fields_.CDy;
  RealArr1DHost Dx = fields_.Dx, Dy = fields_.Dy, Ex = fields_.Ex, Ey = fields_.Ey;
  parallelForSurface("swe_momentum_source", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const int c = g.getSurfaceIndex(i, j);
    const int iM = g.getSurfaceIndex(i - 1, j), iP = g.getSurfaceIndex(i + 1, j);
    const int jM = g.getSurfaceIndex(i, j - 1), jP = g.getSurfaceIndex(i, j + 1);
    const real uu = uu_(c), vv = vv_(c), vx = vx_(c), uy = uy_(c);

    real advX = (0.5 / dx) * ((uu + Kokkos::fabs(uu)) * (uu - uu_(iM)) +
                              (uu - Kokkos::fabs(uu)) * (uu_(iP) - uu)) +
                (0.5 / dy) * ((vx + Kokkos::fabs(vx)) * (uu - uu_(jM)) +
                              (vx - Kokkos::fabs(vx)) * (uu_(jP) - uu));
    real advY = (0.5 / dx) * ((uy + Kokkos::fabs(uy)) * (vv - vv_(iM)) +
                              (uy - Kokkos::fabs(uy)) * (vv_(iP) - vv)) +
                (0.5 / dy) * ((vv + Kokkos::fabs(vv)) * (vv - vv_(jM)) +
                              (vv - Kokkos::fabs(vv)) * (vv_(jP) - vv));
    if (uu == 0.0 || cflx(c) > 0.7)
      advX = 0.0;
    else if (cflx(c) > 0.5)
      advX = advX * (0.7 - cflx(c)) / (0.7 - 0.5);
    if (vv == 0.0 || cfly(c) > 0.7)
      advY = 0.0;
    else if (cfly(c) > 0.5)
      advY = advY * (0.7 - cfly(c)) / (0.7 - 0.5);

    real difX = 0.0, difY = 0.0;
    if (Vsx(c) > 0.0) {
      difX = (viscx / Vsx(c) / dx) * (Asx(c) * (uu_(iP) - uu) - Asx(c) * (uu - uu_(iM))) +
             (viscy / Vsx(c) / dy) * (Asy(c) * (uu_(jP) - uu) - Asy(jM) * (uu - uu_(jM)));
    }
    if (Vsy(c) > 0.0) {
      difY = (viscx / Vsy(c) / dx) * (Asx(c) * (vv_(iP) - vv) - Asx(iM) * (vv - vv_(iM))) +
             (viscy / Vsy(c) / dy) * (Asy(c) * (vv_(jP) - vv) - Asy(c) * (vv - vv_(jM)));
    }

    const real velx = Kokkos::sqrt(uu * uu + vx * vx);
    const real vely = Kokkos::sqrt(uy * uy + vv * vv);
    real facdx = 0.0, facdy = 0.0;
    if (Vsx(c) > 0.0) facdx = Aszx(c) / Vsx(c);
    if (Vsy(c) > 0.0) facdy = Aszy(c) / Vsy(c);

    Dx(c) = 1.0 / (0.5 * dt * CDx(c) * velx * facdx + 1.0);
    Dy(c) = 1.0 / (0.5 * dt * CDy(c) * vely * facdy + 1.0);
    Ex(c) = (uu + dt * (difX - advX)) * Dx(c);
    Ey(c) = (vv + dt * (difY - advY)) * Dy(c);
  });
}

void SweSolver::shallowwaterRhs() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const real dt = params_.dt;
  const real hE = params_.hE;
  const Grid g = grid_;
  RealArr1DHost Srhs = fields_.Srhs, eta = fields_.eta, Asz = fields_.Asz, Asx = fields_.Asx;
  RealArr1DHost Ex = fields_.Ex, Asy = fields_.Asy, Ey = fields_.Ey, bottom = fields_.bottom;
  RealArr1DHost evap = fields_.evap;
  parallelForSurface("swe_rhs", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const int c = g.getSurfaceIndex(i, j);
    const int iM = g.getSurfaceIndex(i - 1, j), jM = g.getSurfaceIndex(i, j - 1);
    Srhs(c) = eta(c) * Asz(c) -
              dt * (Asx(c) * Ex(c) - Asx(iM) * Ex(iM) + Asy(c) * Ey(c) - Asy(jM) * Ey(jM));
    if (eta(c) - bottom(c) > hE) Srhs(c) -= evap(c) * Asz(c) * dt;
  });
}

void SweSolver::shallowwaterMatCoeff() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const real dx = grid_.dx();
  const real dy = grid_.dy();
  const real dt = params_.dt;
  const real coef = params_.gravity * dt * dt;
  const Grid g = grid_;
  const int i0g = mc_ ? mc_->i0() : 0;
  const int j0g = mc_ ? mc_->j0() : 0;
  const int gnx = globalNx();
  const int gny = globalNy();
  const bool par = parallel();
  const int px = mc_ ? mc_->px() : 0;
  const int py = mc_ ? mc_->py() : 0;
  const int mnx = mc_ ? mc_->mpiNx() : 1;
  const int mny = mc_ ? mc_->mpiNy() : 1;
  RealArr1DHost Sxp = fields_.Sxp, Sxm = fields_.Sxm, Syp = fields_.Syp, Sym = fields_.Sym;
  RealArr1DHost Sct = fields_.Sct, Srhs = fields_.Srhs, Asz = fields_.Asz;
  RealArr1DHost Asx = fields_.Asx, Asy = fields_.Asy, Dx = fields_.Dx, Dy = fields_.Dy;
  RealArr1DHost Vsx = fields_.Vsx, Vsy = fields_.Vsy, dept = fields_.dept, eta = fields_.eta;
  RealArr1DHost uu = fields_.uu, vv = fields_.vv;
  parallelForSurface("swe_mat_coeff", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const int c = g.getSurfaceIndex(i, j);
    const int iM = g.getSurfaceIndex(i - 1, j), jM = g.getSurfaceIndex(i, j - 1);
    Sxp(c) = Sxm(c) = Syp(c) = Sym(c) = 0.0;
    if (Vsx(c) > 0.0) Sxp(c) = coef * Asx(c) * Asx(c) * Dx(c) / Vsx(c);
    if (Vsx(iM) > 0.0) Sxm(c) = coef * Asx(iM) * Asx(iM) * Dx(iM) / Vsx(iM);
    if (Vsy(c) > 0.0) Syp(c) = coef * Asy(c) * Asy(c) * Dy(c) / Vsy(c);
    if (Vsy(jM) > 0.0) Sym(c) = coef * Asy(jM) * Asy(jM) * Dy(jM) / Vsy(jM);
    Sct(c) = Asz(c) + Sxp(c) + Sxm(c) + Syp(c) + Sym(c);

    if (dept(c) == 0.0) {
      Sct(c) = dx * dy;
      Srhs(c) = eta(c) * dx * dy;
      if (uu(c) == 0.0 && uu(iM) == 0.0 && vv(c) == 0.0 && vv(jM) == 0.0) {
        Sxp(c) = Sxm(c) = Syp(c) = Sym(c) = 0.0;
      }
    } else if (Sct(c) == 0.0) {
      Sct(c) = dx * dy;
      Srhs(c) = eta(c) * dx * dy;
    }
    // x-boundary (legacy shallowwater_mat_coeff with MPI partition folding).
    const int gi = i0g + i;
    const int gj = j0g + j;
    if (gi == 0) {
      if (!par || px == 0) Sct(c) -= Sxm(c);
      else Srhs(c) += Sxm(c) * eta(g.getSurfaceIndex(i - 1, j));
      if (gnx == 1) Sct(c) -= Sxp(c);
    } else if (gi == gnx - 1) {
      if (!par || px == mnx - 1) Sct(c) -= Sxp(c);
      else Srhs(c) += Sxp(c) * eta(g.getSurfaceIndex(i + 1, j));
    }
    if (gj == 0) {
      if (!par || py == 0) Sct(c) -= Sym(c);
      else Srhs(c) += Sym(c) * eta(g.getSurfaceIndex(i, j - 1));
      if (gny == 1) Sct(c) -= Syp(c);
    } else if (gj == gny - 1) {
      if (!par || py == mny - 1) Sct(c) -= Syp(c);
      else Srhs(c) += Syp(c) * eta(g.getSurfaceIndex(i, j + 1));
    }
  });
}

void SweSolver::assembleAndSolve() {
  perf::ScopedTimer t_asm(perf_, perf::Region::SwAssembly);
  const int gnx = globalNx();
  const int gny = globalNy();
  const auto corners = dd_->localCorners();
  const int i0 = std::get<0>(corners), j0 = std::get<1>(corners);
  const int ni = std::get<2>(corners), nj = std::get<3>(corners);
  const int rstart = dd_->ownershipRange().first;
  auto& f = fields_;

  sys_->zeroEntries();
  sys_->beginAssembly();
  RealArr1DHost b_host = b_owned_;
  int64_t bytes = 0;
  for (int gj = j0; gj < j0 + nj; ++gj)
    for (int gi = i0; gi < i0 + ni; ++gi) {
      const int lc = s(gi - i0, gj - j0);
      const int grow = dd_->globalRow(gi, gj);
      int cols[kMaxStencil];
      int ncols = 0;
      dd_->stencilColumns(gi, gj, cols, ncols);
      const int rowW = (gi - 1 >= 0) ? dd_->globalRow(gi - 1, gj) : -1;
      const int rowE = (gi + 1 < gnx) ? dd_->globalRow(gi + 1, gj) : -1;
      const int rowS = (gj - 1 >= 0) ? dd_->globalRow(gi, gj - 1) : -1;
      const int rowN = (gj + 1 < gny) ? dd_->globalRow(gi, gj + 1) : -1;
      real vals[kMaxStencil];
      for (int cc = 0; cc < ncols; ++cc) {
        const int col = cols[cc];
        if (col == grow)
          vals[cc] = f.Sct(lc);
        else if (col == rowW)
          vals[cc] = -f.Sxm(lc);
        else if (col == rowE)
          vals[cc] = -f.Sxp(lc);
        else if (col == rowS)
          vals[cc] = -f.Sym(lc);
        else if (col == rowN)
          vals[cc] = -f.Syp(lc);
        else
          vals[cc] = 0.0;
      }
      sys_->addRow(grow, ncols, cols, vals);
      b_host(static_cast<size_t>(grow - rstart)) = f.Srhs(lc);
      bytes += static_cast<int64_t>(ncols) * (2 * static_cast<int>(sizeof(int)) + static_cast<int>(sizeof(real)));
    }
  if (b_host.data() != b_owned_.data()) Kokkos::deep_copy(b_owned_, b_host);
  if (perf_) perf_->counters().addBytes(bytes + static_cast<int64_t>(ni) * nj * static_cast<int>(sizeof(real)));
  sys_->endAssembly();

  {
    perf::ScopedTimer t_ksp(perf_, perf::Region::SwKsp);
    // PetscLinearSolver::solve() calls setup() internally; avoid a duplicate setup here (P21).
    solver_->solve(*sys_, b_owned_, x_owned_);
    if (perf_) perf_->counters().addKspIters(solver_->getIterationCount());
  }

  RealArr1DHost x_host = x_owned_;
  const Grid g = grid_;
  Decomp2D* dd = dd_;
  RealArr1DHost eta = f.eta;
  const int rs = rstart;
  parallelForSurface("swe_solve_unpack", ni, nj, KOKKOS_LAMBDA(int li, int lj) {
    const int gi = i0 + li;
    const int gj = j0 + lj;
    const int grow = dd->globalRow(gi, gj);
    eta(g.getSurfaceIndex(li, lj)) = x_host(static_cast<size_t>(grow - rs));
  });
}

void SweSolver::cflLimiter() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const real min_dept = params_.min_depth;
  const Grid g = grid_;
  RealArr1DHost eta = fields_.eta, bottom = fields_.bottom, dept = fields_.dept;
  RealArr1DHost cfl_active = fields_.cfl_active;
  parallelForSurface("swe_cfl_wet", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const int c = g.getSurfaceIndex(i, j);
    const real diff = eta(c) - bottom(c);
    if (dept(c) <= 0.0 && diff > 0.0) {
      int wet = 0;
      if (dept(g.getSurfaceIndex(i + 1, j)) > 0.0) wet = 1;
      else if (dept(g.getSurfaceIndex(i - 1, j)) > 0.0) wet = 1;
      else if (dept(g.getSurfaceIndex(i, j + 1)) > 0.0) wet = 1;
      else if (dept(g.getSurfaceIndex(i, j - 1)) > 0.0) wet = 1;
      if (wet == 0) {
        eta(c) = bottom(c);
        cfl_active(c) = 1.0;
      }
    }
  });
  parallelForSurface("swe_cfl_smalldepth", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const int c = g.getSurfaceIndex(i, j);
    const real diff = eta(c) - bottom(c);
    if (diff > 0.0 && diff < min_dept) eta(c) = bottom(c);
  });
}

void SweSolver::evaprain(real rain_rate, real evap_rate) {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const real dt = params_.dt;
  const real min_dept = params_.min_depth;
  const real hE = params_.hE;
  const Grid g = grid_;
  const int j0g = mc_ ? mc_->j0() : 0;
  const int gny = globalNy();
  RealArr1DHost eta = fields_.eta, bottom = fields_.bottom, dept = fields_.dept;
  // Rainfall (sim_groundwater==0 path): add to all interior cells except global y+ row.
  parallelForSurface("swe_rain", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    if (j0g + j != gny - 1) eta(g.getSurfaceIndex(i, j)) += rain_rate * dt;
  });
  // Evaporation over layers thicker than hE.
  parallelForSurface("swe_evap", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const int c = g.getSurfaceIndex(i, j);
    if (eta(c) - bottom(c) > hE) eta(c) -= evap_rate * dt;
  });
  // Remove negative/small depth and refresh dept.
  parallelForSurface("swe_evaprain_depth", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const int c = g.getSurfaceIndex(i, j);
    const real diff = eta(c) - bottom(c);
    if (diff < min_dept) eta(c) = bottom(c);
    dept(c) = eta(c) - bottom(c);
  });
}

void SweSolver::updateDragCoef() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const real dx = grid_.dx();
  const real dy = grid_.dy();
  const real coef = params_.gravity * params_.manning * params_.manning;
  const real hD = params_.hD;
  const Grid g = grid_;
  RealArr1DHost Vs = fields_.Vs, CDx = fields_.CDx, CDy = fields_.CDy;
  parallelForSurface("swe_drag", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const int c = g.getSurfaceIndex(i, j);
    if (Vs(c) > 0.0) {
      const real effh = Vs(c) / (dx * dy);
      const real expo = (effh < hD) ? (2.0 / 3.0) : (1.0 / 3.0);
      CDx(c) = coef / Kokkos::pow(effh, expo);
      CDy(c) = coef / Kokkos::pow(effh, expo);
    }
  });
}

void SweSolver::waterfallLocation() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const real wtfh = params_.wtfh;
  const Grid g = grid_;
  const int n = grid_.nSurfaceStorageCell();
  RealArr1DHost wtfx = fields_.wtfx, wtfy = fields_.wtfy, eta = fields_.eta;
  RealArr1DHost bottomXP = fields_.bottomXP, bottomYP = fields_.bottomYP;
  parallelForRange("swe_wtf_reset", n, KOKKOS_LAMBDA(int idx) {
    wtfx(idx) = 0.0;
    wtfy(idx) = 0.0;
  });
  parallelForSurface("swe_wtf_fwd", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const int c = g.getSurfaceIndex(i, j);
    real facd = eta(g.getSurfaceIndex(i + 1, j)) - bottomXP(c);
    if (eta(c) < bottomXP(c) && facd > wtfh) wtfx(c) = -1.0;
    facd = eta(g.getSurfaceIndex(i, j + 1)) - bottomYP(c);
    if (eta(c) < bottomYP(c) && facd > wtfh) wtfy(c) = -1.0;
  });
  parallelForSurface("swe_wtf_bwd", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const int c = g.getSurfaceIndex(i, j);
    const int iM = g.getSurfaceIndex(i - 1, j), jM = g.getSurfaceIndex(i, j - 1);
    real facd = eta(iM) - bottomXP(iM);
    if (eta(c) < bottomXP(iM) && facd > wtfh) wtfx(c) = 1.0;
    facd = eta(jM) - bottomYP(jM);
    if (eta(c) < bottomYP(jM) && facd > wtfh) wtfy(c) = 1.0;
  });
}

void SweSolver::updateVelocity() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const real dx = grid_.dx();
  const real dy = grid_.dy();
  const real dt = params_.dt;
  const real wtfh = params_.wtfh;
  const real coef = params_.gravity * dt;
  const Grid g = grid_;
  const int n = grid_.nSurfaceStorageCell();
  {
    RealArr1DHost un = fields_.un, vn = fields_.vn, uu = fields_.uu, vv = fields_.vv;
    parallelForRange("swe_vel_save", n, KOKKOS_LAMBDA(int idx) {
      un(idx) = uu(idx);
      vn(idx) = vv(idx);
    });
  }
  {
    RealArr1DHost uu = fields_.uu, vv = fields_.vv, Vsx = fields_.Vsx, Vsy = fields_.Vsy;
    RealArr1DHost Asx = fields_.Asx, Asy = fields_.Asy, Ex = fields_.Ex, Ey = fields_.Ey;
    RealArr1DHost Dx = fields_.Dx, Dy = fields_.Dy, eta = fields_.eta;
    parallelForSurface("swe_vel_compute", nx, ny, KOKKOS_LAMBDA(int i, int j) {
      const int c = g.getSurfaceIndex(i, j);
      uu(c) = 0.0;
      vv(c) = 0.0;
      const real effhx = (Vsx(c) > 0.0) ? Asx(c) / Vsx(c) : 0.0;
      const real effhy = (Vsy(c) > 0.0) ? Asy(c) / Vsy(c) : 0.0;
      uu(c) = (Ex(c) - coef * effhx * (eta(g.getSurfaceIndex(i + 1, j)) - eta(c))) * Dx(c);
      vv(c) = (Ey(c) - coef * effhy * (eta(g.getSurfaceIndex(i, j + 1)) - eta(c))) * Dy(c);
    });
  }
  // Wet/dry + CFL velocity limiters write the WEST/SOUTH neighbour faces (uu(iM), vv(jM)),
  // so this pass has a cross-cell write dependency and stays SEQUENTIAL (legacy order).
  {
    auto& f = fields_;
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        const int c = s(i, j);
        const int iM = s(i - 1, j), jM = s(i, j - 1);
        if (f.Asx(c) < wtfh * dy) f.uu(c) = 0.0;
        if (f.Asy(c) < wtfh * dx) f.vv(c) = 0.0;
        if (f.dept(c) < wtfh) {
          if (f.uu(c) > 0.0) f.uu(c) = 0.0;
          if (f.uu(iM) < 0.0) f.uu(iM) = 0.0;
          if (f.vv(c) > 0.0) f.vv(c) = 0.0;
          if (f.vv(jM) < 0.0) f.vv(jM) = 0.0;
        }
        if (f.cfl_active(c) == 1.0) {
          f.uu(c) = 0.0;
          f.uu(iM) = 0.0;
          f.vv(c) = 0.0;
          f.vv(jM) = 0.0;
          f.cfl_active(c) = 0.0;
        }
      }
  }
  {
    RealArr1DHost Fu = fields_.Fu, Fv = fields_.Fv, uu = fields_.uu, vv = fields_.vv;
    RealArr1DHost Asx = fields_.Asx, Asy = fields_.Asy, cflx = fields_.cflx, cfly = fields_.cfly;
    parallelForSurface("swe_vel_flux", nx, ny, KOKKOS_LAMBDA(int i, int j) {
      const int c = g.getSurfaceIndex(i, j);
      Fu(c) = uu(c) * Asx(c);
      Fv(c) = vv(c) * Asy(c);
      cflx(c) = Kokkos::fabs(uu(c) * dt / dx);
      cfly(c) = Kokkos::fabs(vv(c) * dt / dy);
    });
  }
}

void SweSolver::enforceVeloBc() {
  auto& f = fields_;
  fillDomainEdgeGhosts(f.uu);
  fillDomainEdgeGhosts(f.vv);
  exchangeHalo(f.uu);
  exchangeHalo(f.vv);
}

void SweSolver::interpVelocity() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const Grid g = grid_;
  RealArr1DHost uy = fields_.uy, vx = fields_.vx, uu = fields_.uu, vv = fields_.vv;
  parallelForSurface("swe_interp_vel", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const int c = g.getSurfaceIndex(i, j);
    uy(c) = 0.25 * (uu(c) + uu(g.getSurfaceIndex(i - 1, j)) + uu(g.getSurfaceIndex(i, j + 1)) +
                    uu(g.getSurfaceIndex(i - 1, j + 1)));
    vx(c) = 0.25 * (vv(c) + vv(g.getSurfaceIndex(i, j - 1)) + vv(g.getSurfaceIndex(i + 1, j)) +
                    vv(g.getSurfaceIndex(i + 1, j - 1)));
  });
}

void SweSolver::attachSolver(LinearSolver& solver, Decomp2D& dd) {
  solver_ = &solver;
  dd_ = &dd;
  sys_ = solver.createSystem(dd);
  const int nloc = dd.ownedRowCount();
  b_owned_ = RealArr1D("swe_b", static_cast<size_t>(nloc));
  x_owned_ = RealArr1D("swe_x", static_cast<size_t>(nloc));
}

void SweSolver::advanceStep(real rain_rate, real evap_rate) {
  if (perf_) perf_->counters().addCells(static_cast<int64_t>(grid_.nx()) * grid_.ny());
  auto& f = fields_;
  const int n = grid_.nSurfaceStorageCell();
  {
    perf::ScopedTimer t(perf_, perf::Region::SwUpdate);
    RealArr1DHost rain = fields_.rain, evap = fields_.evap, etan = fields_.etan, eta = fields_.eta;
    parallelForRange("swe_step_save", n, KOKKOS_LAMBDA(int idx) {
      rain(idx) = rain_rate;
      evap(idx) = evap_rate;
      etan(idx) = eta(idx);
    });
    exchangeStepFieldsPreSolve();
    enforceSurfBc();
    if (parallel()) {
      exchangeHalo(f.etan);
      exchangeHalo(f.eta);
    }
    momentumSource();
    exchangeMomentumCoeffs();
    shallowwaterRhs();
    shallowwaterMatCoeff();
  }
  assembleAndSolve();
  {
    perf::ScopedTimer t(perf_, perf::Region::SwUpdate);
    enforceSurfBc();
    cflLimiter();
    evaprain(rain_rate, evap_rate);
    updateDepth();
    if (parallel()) exchangeHalo(f.eta);
    updateDepth();
    exchangeDepthFields();
    exchangeSubgridGeometry();
    updateDragCoef();
    waterfallLocation();
    updateVelocity();
    enforceVeloBc();
    if (parallel()) {
      exchangeHalo(f.uu);
      exchangeHalo(f.vv);
    }
    interpVelocity();
    if (parallel()) {
      exchangeHalo(f.uy);
      exchangeHalo(f.vx);
    }
    updateGeometry();
    exchangeStepFieldsPostVelocity();
  }
}

void SweSolver::setInitialEtaPhysical(const RealArr1DHost& physical_eta) {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  if (physical_eta.extent(0) != static_cast<size_t>(nx * ny))
    throw std::runtime_error("SweSolver::setInitialEtaPhysical: expected nx*ny values");
  const Grid g = grid_;
  const real offset = params_.offset;
  RealArr1DHost eta = fields_.eta;
  RealArr1DHost bottom = fields_.bottom;
  RealArr1DHost phys = physical_eta;
  parallelForSurface("swe_ic_eta", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const int c = g.getSurfaceIndex(i, j);
    eta(c) = phys(static_cast<size_t>(i + j * nx)) + offset;
    if (eta(c) < bottom(c)) eta(c) = bottom(c);
  });
}

void SweSolver::setInitialVelocityConstant(real u, real v) {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const Grid g = grid_;
  RealArr1DHost uu = fields_.uu, vv = fields_.vv, CDx = fields_.CDx, CDy = fields_.CDy;
  parallelForSurface("swe_ic_vel_const", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const int c = g.getSurfaceIndex(i, j);
    uu(c) = u;
    vv(c) = v;
    CDx(c) = 0.0;
    CDy(c) = 0.0;
  });
}

void SweSolver::setInitialVelocityPhysical(const RealArr1DHost& u,
                                           const RealArr1DHost& v) {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  if (u.extent(0) != static_cast<size_t>(nx * ny) || v.extent(0) != static_cast<size_t>(nx * ny))
    throw std::runtime_error("SweSolver::setInitialVelocityPhysical: expected nx*ny values");
  const Grid g = grid_;
  RealArr1DHost uu = fields_.uu, vv = fields_.vv, CDx = fields_.CDx, CDy = fields_.CDy;
  RealArr1DHost uin = u, vin = v;
  parallelForSurface("swe_ic_vel_phys", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    const int c = g.getSurfaceIndex(i, j);
    uu(c) = uin(static_cast<size_t>(i + j * nx));
    vv(c) = vin(static_cast<size_t>(i + j * nx));
    CDx(c) = 0.0;
    CDy(c) = 0.0;
  });
}

void SweSolver::finalizeInitialState() {
  enforceSurfBc();
  fillDomainEdgeGhosts(fields_.uu);
  fillDomainEdgeGhosts(fields_.vv);
  if (parallel()) {
    exchangeHalo(fields_.uu);
    exchangeHalo(fields_.vv);
  }
  updateDepth();
  updateGeometry();
  if (parallel()) exchangeStepFieldsPostVelocity();
  const int n = grid_.nSurfaceStorageCell();
  RealArr1DHost etan = fields_.etan, eta = fields_.eta, un = fields_.un, uu = fields_.uu;
  RealArr1DHost vn = fields_.vn, vv = fields_.vv, Vsn = fields_.Vsn, Vs = fields_.Vs;
  parallelForRange("swe_finalize_save", n, KOKKOS_LAMBDA(int idx) {
    etan(idx) = eta(idx);
    un(idx) = uu(idx);
    vn(idx) = vv(idx);
    Vsn(idx) = Vs(idx);
  });
}

void SweSolver::initializeState(real init_eta) {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  RealArr1DHost physical("ic_eta", static_cast<size_t>(nx * ny));
  for (size_t k = 0; k < physical.extent(0); ++k) physical(k) = init_eta;
  setInitialEtaPhysical(physical);
  setInitialVelocityConstant(0.0, 0.0);
  finalizeInitialState();
}

void SweSolver::initialize(const Config& config) {
  SweParams p;
  p.gravity = config.getOr<double>("surface_water.gravity", 9.81);
  p.manning = config.getOr<double>("surface_water.manning", 0.0);
  p.min_depth = config.getOr<double>("surface_water.min_depth", 1.0e-8);
  p.visc_x = config.getOr<double>("surface_water.viscosity.x", 0.0);
  p.visc_y = config.getOr<double>("surface_water.viscosity.y", 0.0);
  p.hD = config.getOr<double>("surface_water.h_diffusion_ref", 0.1);
  p.wtfh = config.getOr<double>("surface_water.waterfall_depth", 1.0e-8);
  p.dt = config.getOr<double>("time.dt", 1.0);
  setParams(p);

  const real botz = config.getOr<double>("domain.botz", 0.0);
  setBathymetryConstant(botz);
  const real init_eta = config.getOr<double>("initial_conditions.surface_water.eta", botz);
  initializeState(init_eta);
}

real SweSolver::maxCfl() const {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const real grav = params_.gravity;
  const real dt = params_.dt;
  const real dmin = std::min(grid_.dx(), grid_.dy());
  const Grid gr = grid_;
  RealArr1DHost dept = fields_.dept, uu = fields_.uu, vv = fields_.vv;
  real cfl = 0.0;
  Kokkos::parallel_reduce(
      "swe_max_cfl", Kokkos::MDRangePolicy<LoopExec, Kokkos::Rank<2>>({0, 0}, {nx, ny}),
      KOKKOS_LAMBDA(int i, int j, real& lmax) {
        const int c = gr.getSurfaceIndex(i, j);
        const real h = dept(c);
        if (h <= 0.0) return;
        const real speed = Kokkos::sqrt(uu(c) * uu(c) + vv(c) * vv(c));
        const real local = (speed + Kokkos::sqrt(grav * h)) * dt / dmin;
        if (local > lmax) lmax = local;
      },
      Kokkos::Max<real>(cfl));
  return cfl;
}

void SweSolver::gatherDeptGlobal(RealArr1DHost& global_dept) const {
  if (mc_ == nullptr || dd_ == nullptr) {
    const int nx = grid_.nx();
    const int ny = grid_.ny();
    const Grid g = grid_;
    RealArr1DHost dept = fields_.dept;
    RealArr1DHost gd = global_dept;
    parallelForSurface("swe_gather_serial", nx, ny, KOKKOS_LAMBDA(int i, int j) {
      gd(static_cast<size_t>(i + j * nx)) = dept(g.getSurfaceIndex(i, j));
    });
    return;
  }
  RealArr1D dev("gather_dept", fields_.dept.extent(0));
  Kokkos::deep_copy(dev, fields_.dept);
  mc_->gatherToRank0(dev, grid_.nx() + 2, grid_.ny() + 2, global_dept);
}

}  // namespace frehg2
