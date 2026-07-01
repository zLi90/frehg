#include "re/ReSolver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "frehg2/core/ParallelFor.hpp"
#include "frehg2/linear/LinearSolver.hpp"
#include "frehg2/linear/SparseSystem.hpp"
#include "frehg2/perf/PerfRecorder.hpp"
#include "frehg2/perf/Timer.hpp"
#include "linear/DomainDecomposition.hpp"
#include "soil/SoilMap.hpp"

namespace frehg2 {

namespace {

// Vertical Darcy flux through the +z (down) face of cell (i,j,k), or the top back face when
// `back`. `soilC` is the own cell's soil; `soilZp` is the cell below (k+1)'s soil (== soilC for
// a 2D/uniform map). Per-cell Ks_z lets a 3D heterogeneous map vary conductivity vertically.
KOKKOS_INLINE_FUNCTION real darcyFluxZ(const ReParams& p, const GwFields& f,
                                       const SoilParams& soilC, const SoilParams& soilZp, int i,
                                       int j, int k, const Grid& grid, bool back) {
  const int c = grid.getIndex(i, j, k);
  const real faceA = p.dx * p.dy;
  if (back && k == 0) {
    const int km = grid.getIndex(i, j, -1);
    const real delta = 0.5 * f.dz3d(c);
    real kface = soilC.Ks_z;
    if (p.bc_type[5] == 0) kface = 0.0;
    // Legacy fixed-flux top (bctype_GW[5]==2): the top back-face conductivity is the cell's own
    // pressure-dependent K (groundwater.c:288 `Kz[icjckM]=Kp=compute_K(h[ii])`), not Ks_z.
    else if (p.bc_type[5] == 2)
      kface = VanGenuchten::conductivityFromHead(soilC, f.h(c), soilC.Ks_z);
    // Legacy darcy_flux z/back: ic=ghost (km), is=cell (c) → q ∝ h(cell) - h(ghost).
    // The "- rho" term is the gravity gradient dz/dz = r_rhozp (1 for non-baroclinic).
    const real hc = f.h(km);
    const real hs = f.h(c);
    const real rho = f.r_rhozp(km);
    const real visc = f.r_visczp(km);
    return faceA * kface * visc * ((hs - hc) / delta - rho);
  }
  const int kp = grid.getIndex(i, j, k + 1);
  real delta = 0.5 * (f.dz3d(c) + f.dz3d(kp));
  if (k == grid.nz() - 1) {
    delta = 0.5 * f.dz3d(c);
    if (p.bc_type[4] == 0) return 0.0;
  }
  const real kc = VanGenuchten::conductivityFromHead(soilC, f.h(c), soilC.Ks_z);
  const real ks = VanGenuchten::conductivityFromHead(soilZp, f.h(kp), soilZp.Ks_z);
  const real kface = 0.5 * (kc + ks);
  const real rho = f.r_rhozp(c);
  const real visc = f.r_visczp(c);
  return faceA * kface * visc * ((f.h(kp) - f.h(c)) / delta - rho);
}
}  // namespace

ReSolver::ReSolver(const Grid& grid, const MpiComm* mc)
    : grid_(grid), fields_(grid), mc_(mc) {
  qtop_field_ = RealArr1DHost("re_qtop", static_cast<size_t>(grid_.nSurfaceStorageCell()));
}

ReSolver::~ReSolver() = default;

void ReSolver::setTopFluxField(const std::vector<real>& q_local_interior) {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  if (static_cast<int>(q_local_interior.size()) != nx * ny)
    throw std::runtime_error("ReSolver::setTopFluxField: expected nx*ny=" +
                             std::to_string(nx * ny) + " values, got " +
                             std::to_string(q_local_interior.size()));
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      qtop_field_(grid_.getSurfaceIndex(i, j)) =
          q_local_interior[static_cast<size_t>(i + j * nx)];
  has_qtop_field_ = true;
}

bool ReSolver::active(int i, int j, int k) const {
  const int idx = s(i, j, k);
  return fields_.actv(idx) > 0.5;
}

const SoilParams& ReSolver::soilAt(int i, int j, int k) const {
  if (soil_map_ == nullptr) return params_.soil;
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const int nz = grid_.nz();
  const int ci = i < 0 ? 0 : (i >= nx ? nx - 1 : i);
  const int cj = j < 0 ? 0 : (j >= ny ? ny - 1 : j);
  const int ck = k < 0 ? 0 : (k >= nz ? nz - 1 : k);
  return soil_map_->classAt(ci, cj, ck);
}

real ReSolver::dzCell(int k) const { return fields_.dz3d(s(0, 0, k)); }

real ReSolver::cellVolume(int i, int j, int k) const {
  return params_.dx * params_.dy * dzCell(k);
}

void ReSolver::buildColumnGeometry() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const int nz = grid_.nz();
  const real dz0 = params_.dz;
  const real botz = params_.botz;
  (void)botz;

  // Legacy cartesian bot1d: top cell bottom at bath - dz (bath=0 for b2-gw flat top).
  RealArr1DHost bot1d("bot1d", static_cast<size_t>(nz));
  bot1d(0) = -dz0;
  for (int k = 1; k < nz; ++k) bot1d(static_cast<size_t>(k)) = bot1d(static_cast<size_t>(k - 1)) - dz0;

  const Grid g = grid_;
  RealArr1DHost bot3d = fields_.bot3d, dz3d = fields_.dz3d, actv = fields_.actv;
  RealArr1DHost r_rho = fields_.r_rho, r_viscxp = fields_.r_viscxp;
  RealArr1DHost r_viscyp = fields_.r_viscyp, r_visczp = fields_.r_visczp;
  parallelForVolume("re_build_geom", nx, ny, nz, KOKKOS_LAMBDA(int i, int j, int k) {
    const int idx = g.getIndex(i, j, k);
    bot3d(idx) = bot1d(static_cast<size_t>(k));
    dz3d(idx) = dz0;
    actv(idx) = 1.0;
    r_rho(idx) = 1.0;
    r_viscxp(idx) = 1.0;
    r_viscyp(idx) = 1.0;
    r_visczp(idx) = 1.0;
  });
}

void ReSolver::initializeGeometry() { buildColumnGeometry(); }

void ReSolver::setInitialWaterContent(const std::vector<double>& wc_local) {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const int nz = grid_.nz();
  const size_t expect = static_cast<size_t>(nx * ny * nz);
  if (wc_local.size() != expect)
    throw std::runtime_error("ReSolver::setInitialWaterContent: expected nx*ny*nz values");
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        const SoilParams& soil = soilAt(i, j, k);
        const int idx = s(i, j, k);
        const double wc = wc_local[static_cast<size_t>(i + j * nx + k * nx * ny)];
        fields_.wc(idx) = wc;
        fields_.h(idx) = VanGenuchten::headFromWaterContent(soil, wc);
        fields_.wcn(idx) = wc;
        fields_.hn(idx) = fields_.h(idx);
      }
}

void ReSolver::setInitialHead(const std::vector<double>& h_local) {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const int nz = grid_.nz();
  const size_t expect = static_cast<size_t>(nx * ny * nz);
  if (h_local.size() != expect)
    throw std::runtime_error("ReSolver::setInitialHead: expected nx*ny*nz values");
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        const SoilParams& soil = soilAt(i, j, k);
        const int idx = s(i, j, k);
        const double h = h_local[static_cast<size_t>(i + j * nx + k * nx * ny)];
        fields_.h(idx) = h;
        fields_.wc(idx) = VanGenuchten::waterContentFromHead(soil, h);
        fields_.wcn(idx) = fields_.wc(idx);
        fields_.hn(idx) = h;
      }
}

void ReSolver::finalizeInitialState() {
  updateDiagnostics();
  dtn_ = params_.dt;
}

void ReSolver::initializeUniformColumn(real init_wc) {
  buildColumnGeometry();
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const int nz = grid_.nz();
  std::vector<double> wc_local(static_cast<size_t>(nx * ny * nz), init_wc);
  setInitialWaterContent(wc_local);
  finalizeInitialState();
}

void ReSolver::savePreviousStepFields() {
  const int n = grid_.nCell();
  RealArr1DHost hn = fields_.hn, h = fields_.h, wcn = fields_.wcn, wc = fields_.wc;
  parallelForRange("re_save_prev", n, KOKKOS_LAMBDA(int idx) {
    hn(idx) = h(idx);
    wcn(idx) = wc(idx);
  });
}

void ReSolver::updateDiagnostics() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const int nz = grid_.nz();
  const Grid g = grid_;
  const ReSolver* self = this;
  RealArr1DHost hwc = fields_.hwc, wch = fields_.wch, ch = fields_.ch;
  RealArr1DHost wc = fields_.wc, h = fields_.h;
  parallelForVolume("re_diagnostics", nx, ny, nz, KOKKOS_LAMBDA(int i, int j, int k) {
    if (!self->active(i, j, k)) return;
    const SoilParams& soil = self->soilAt(i, j, k);
    const int idx = g.getIndex(i, j, k);
    hwc(idx) = VanGenuchten::headFromWaterContent(soil, wc(idx));
    wch(idx) = VanGenuchten::waterContentFromHead(soil, h(idx));
    ch(idx) = VanGenuchten::specificMoistureCapacity(soil, h(idx));
  });
}

void ReSolver::computeKFace() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const int nz = grid_.nz();
  auto& f = fields_;
  const Grid g = grid_;
  const ReSolver* self = this;
  const int bc2 = params_.bc_type[2];
  const bool full3d = params_.use_full3d;
  // Global-boundary ownership for the no-flux wall closure and the main-loop y+ handling.
  // Serial / no-comm => this rank owns every global boundary.
  const bool ownsYp = (mc_ == nullptr) || (mc_->j0() + ny == mc_->globalNy());
  RealArr1DHost h = fields_.h, Kx = fields_.Kx, Ky = fields_.Ky, Kz = fields_.Kz;

  // 3D MPI: lateral conductivity/flux need neighbor heads. Exchange h ghosts before the face
  // means. Gated behind use_full3d so the vertical-only path is byte-for-byte unchanged.
  if (full3d && mc_ != nullptr && mc_->size() > 1)
    mc_->haloExchange3D(h, nx + 2, ny + 2, nz);

  parallelForVolume("re_kface", nx, ny, nz, KOKKOS_LAMBDA(int i, int j, int k) {
    if (!self->active(i, j, k)) return;
    const int c = g.getIndex(i, j, k);
    // Per-cell soil (P13 per-column / P23 per-cell): own cell plus the +x/+y/+z neighbor cells
    // for the face means (arithmetic mean of own/neighbor conductivity, exactly the P5 formula;
    // uniform soil makes the references the same SoilParams, reproducing P5/P13 bit-for-bit).
    const SoilParams& soilC = self->soilAt(i, j, k);
    // Kx
    {
      const int ip = g.getIndex(i + 1, j, k);
      const SoilParams& soilXp = self->soilAt(i + 1, j, k);
      const real Kp = VanGenuchten::conductivityFromHead(soilXp, h(ip), soilXp.Ks_x);
      const real Km = VanGenuchten::conductivityFromHead(soilC, h(c), soilC.Ks_x);
      Kx(c) = 0.5 * (Kp + Km);
    }
    // Ky
    {
      const int jp = g.getIndex(i, j + 1, k);
      const SoilParams& soilYp = self->soilAt(i, j + 1, k);
      const real Kp = VanGenuchten::conductivityFromHead(soilYp, h(jp), soilYp.Ks_y);
      const real Km = VanGenuchten::conductivityFromHead(soilC, h(c), soilC.Ks_y);
      Ky(c) = 0.5 * (Kp + Km);
      // The y+ wall handling fires at the LOCAL last row for the vertical-only path (P5/P13,
      // unchanged) but only at the GLOBAL y+ boundary for the 3D MPI path (interior rank
      // boundaries keep the neighbor-mean computed above).
      if (j == ny - 1 && (!full3d || ownsYp)) {
        if (bc2 == 0) Ky(c) = 0.0;
        else if (bc2 == 1) Ky(c) = VanGenuchten::conductivityFromHead(soilC, h(c), soilC.Ks_y);
      }
    }
    // Kz (downward k+ face between cell k and k+1). Legacy compute_K_face:
    //   if istop[icjckP] -> Kp  (cell BELOW is a top cell: terrain, never on a flat
    //                            column; do NOT confuse with this cell being the top)
    //   else if actv[ii]==0 -> 0 (skipped: inactive cells are continue'd above)
    //   else if icjckP > n3ci -> Km (cell below is out of the interior domain)
    //   else -> 0.5*(Kp+Km)
    {
      // z faces use the own cell (c) and the cell below (k+1). For a 2D/uniform map both
      // resolve to the same SoilParams (bit-identical to P5/P13); for a 3D map this is the
      // per-cell arithmetic mean across the vertical face.
      const int kp = g.getIndex(i, j, k + 1);
      const SoilParams& soilZp = self->soilAt(i, j, k + 1);
      const real Kp = VanGenuchten::conductivityFromHead(soilZp, h(kp), soilZp.Ks_z);
      const real Km = VanGenuchten::conductivityFromHead(soilC, h(c), soilC.Ks_z);
      if (k + 1 < nz && self->isTop(k + 1)) Kz(c) = Kp;
      else if (k + 1 >= nz) Kz(c) = Km;
      else if (!self->active(i, j, k + 1)) Kz(c) = 0.0;
      else Kz(c) = 0.5 * (Kp + Km);
    }
  });

  if (!full3d) {
    // Vertical-only path (P5/P13): legacy local-boundary handling, unchanged (bit-identical).
    // x- boundary faces
    for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < ny; ++j) {
        const int im = s(-1, j, k);
        const int c0 = s(0, j, k);
        const SoilParams& soil0 = soilAt(0, j, k);
        const real Kp = VanGenuchten::conductivityFromHead(soil0, f.h(im), soil0.Ks_x);
        f.Kx(c0) = 0.5 * (Kp + Kp);
        if (params_.bc_type[0] == 0) f.Kx(s(-1, j, k)) = 0.0;
        if (params_.bc_type[1] == 0) f.Kx(s(nx, j, k)) = 0.0;
      }
    }
    // y- boundary faces
    for (int k = 0; k < nz; ++k) {
      for (int i = 0; i < nx; ++i) {
        const int jm = s(i, -1, k);
        const SoilParams& soil0 = soilAt(i, 0, k);
        const real Kp = VanGenuchten::conductivityFromHead(soil0, f.h(jm), soil0.Ks_y);
        f.Ky(s(i, 0, k)) = 0.5 * (Kp + Kp);
        if (params_.bc_type[3] == 0) f.Ky(s(i, -1, k)) = 0.0;
      }
    }
  } else {
    // 3D path (P23): exchange the +x/+y face conductivities so each rank's -x/-y ghost face
    // (used by Gxm/Gym and the lateral flux divergence) holds the shared face value, then close
    // the global no-flux walls at their correct face. Convention: Kx(c) is the +x face of cell c.
    if (mc_ != nullptr && mc_->size() > 1) {
      mc_->haloExchange3D(f.Kx, nx + 2, ny + 2, nz);
      mc_->haloExchange3D(f.Ky, nx + 2, ny + 2, nz);
    }
    const int gnx = mc_ ? mc_->globalNx() : nx;
    const int gny = mc_ ? mc_->globalNy() : ny;
    const int gi0 = mc_ ? mc_->i0() : 0;
    const int gj0 = mc_ ? mc_->j0() : 0;
    const bool ownsXm = (gi0 == 0);
    const bool ownsXp = (gi0 + nx == gnx);
    const bool ownsYm = (gj0 == 0);
    const bool ownsYpL = (gj0 + ny == gny);
    for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < ny; ++j) {
        if (ownsXm && params_.bc_type[0] == 0) f.Kx(s(-1, j, k)) = 0.0;  // x- wall (-x face of i=0)
        if (ownsXp && params_.bc_type[1] == 0) f.Kx(s(nx - 1, j, k)) = 0.0;  // x+ wall
      }
      for (int i = 0; i < nx; ++i) {
        if (ownsYm && params_.bc_type[3] == 0) f.Ky(s(i, -1, k)) = 0.0;  // y- wall (-y face of j=0)
        if (ownsYpL && params_.bc_type[2] == 0) f.Ky(s(i, ny - 1, k)) = 0.0;  // y+ wall
      }
    }
  }
  // bottom k+ boundary
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      f.Kz(s(i, j, nz)) = 0.0;
      if (params_.bc_type[4] == 0) f.Kz(s(i, j, nz - 1)) = 0.0;
    }
  }
  // top surface k- face
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      if (!isTop(0)) continue;
      const int c = s(i, j, 0);
      const int km = s(i, j, -1);
      const SoilParams& soil = soilAt(i, j, 0);
      const real Kp = VanGenuchten::conductivityFromHead(soil, f.h(c), soil.Ks_z);
      if (params_.sim_shallowwater) {
        (void)km;
      } else if (params_.bc_type[5] == 1) {
        f.Kz(km) = soil.Ks_z;
      } else if (params_.bc_type[5] == 0) {
        f.Kz(km) = 0.0;
      } else {
        f.Kz(km) = Kp;
      }
    }
  }

  if (!params_.use_full3d) {
    for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          const int c = s(i, j, k);
          if (f.wc(c) < soilAt(i, j, k).theta_s) {
            f.Kx(c) = 0.0;
            f.Kx(s(i - 1, j, k)) = 0.0;
            f.Ky(c) = 0.0;
            f.Ky(s(i, j - 1, k)) = 0.0;
          }
        }
      }
    }
  }
}

void ReSolver::baroclinicFace() {
  if (params_.baroclinic) return;
  const int n = grid_.nCell();
  RealArr1DHost r_rhoxp = fields_.r_rhoxp, r_rhoyp = fields_.r_rhoyp, r_rhozp = fields_.r_rhozp;
  RealArr1DHost r_viscxp = fields_.r_viscxp, r_viscyp = fields_.r_viscyp, r_visczp = fields_.r_visczp;
  parallelForRange("re_baroclinic", n, KOKKOS_LAMBDA(int idx) {
    r_rhoxp(idx) = 1.0;
    r_rhoyp(idx) = 1.0;
    r_rhozp(idx) = 1.0;
    r_viscxp(idx) = 1.0;
    r_viscyp(idx) = 1.0;
    r_visczp(idx) = 1.0;
  });
}

void ReSolver::groundwaterMatCoeff() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const int nz = grid_.nz();
  const real dx = params_.dx;
  const real dy = params_.dy;
  const real dtg = params_.dt;
  const Grid g = grid_;
  const ReSolver* self = this;
  const int bc2 = params_.bc_type[2], bc3 = params_.bc_type[3];
  const int bc4 = params_.bc_type[4], bc5 = params_.bc_type[5];
  const bool sim_sw = params_.sim_shallowwater;
  RealArr1DHost Gxp = fields_.Gxp, Gxm = fields_.Gxm, Gyp = fields_.Gyp, Gym = fields_.Gym;
  RealArr1DHost Gzp = fields_.Gzp, Gzm = fields_.Gzm, Gct = fields_.Gct;
  RealArr1DHost Kx = fields_.Kx, Ky = fields_.Ky, Kz = fields_.Kz, ch = fields_.ch, wcn = fields_.wcn;
  RealArr1DHost r_rhoxp = fields_.r_rhoxp, r_viscxp = fields_.r_viscxp;
  RealArr1DHost r_rhoyp = fields_.r_rhoyp, r_viscyp = fields_.r_viscyp;
  RealArr1DHost r_rhozp = fields_.r_rhozp, r_visczp = fields_.r_visczp, r_rho = fields_.r_rho;

  parallelForVolume("re_mat_coeff", nx, ny, nz, KOKKOS_LAMBDA(int i, int j, int k) {
    const int c = g.getIndex(i, j, k);
    if (!self->active(i, j, k)) {
      Gxp(c) = Gxm(c) = Gyp(c) = Gym(c) = Gzp(c) = Gzm(c) = 0.0;
      Gct(c) = 1.0;
      return;
    }
    const SoilParams& soil = self->soilAt(i, j, k);
    const int im = g.getIndex(i - 1, j, k);
    const int jm = g.getIndex(i, j - 1, k);
    const int km = g.getIndex(i, j, k - 1);

    // Lateral coefficients carry the face area (Ax = dy*dz, Ay = dx*dz), matching the
    // legacy `* gmap->Ax[ii] * param->dx / dx^2` form and the vertical Gz* face-area form.
    // (b2-gw and the laterally-decoupled MPI gate have Kx=Ky=0, so this is bit-identical there.)
    const real Ax_k = dy * self->dzCell(k);
    const real Ay_k = dx * self->dzCell(k);
    Gxp(c) = -Kx(c) * dtg * r_rhoxp(c) * r_viscxp(c) * Ax_k / dx;
    Gxm(c) = -Kx(im) * dtg * r_rhoxp(im) * r_viscxp(im) * Ax_k / dx;
    Gyp(c) = -Ky(c) * dtg * r_rhoyp(c) * r_viscyp(c) * Ay_k / dy;
    Gym(c) = -Ky(jm) * dtg * r_rhoyp(jm) * r_viscyp(jm) * Ay_k / dy;
    if (j == ny - 1 && bc2 == 1) Gyp(c) = 2.0 * Gyp(c);
    if (j == 0 && bc3 == 1) Gym(c) = 2.0 * Gym(c);

    const real dz_down = (k + 1 < nz) ? self->dzCell(k + 1) : self->dzCell(k);
    const real dzf_p = 0.5 * (self->dzCell(k) + dz_down);
    Gzp(c) = -Kz(c) * dtg * r_rhozp(c) * r_visczp(c) / (self->dzCell(k) * dzf_p);
    Gzp(c) *= self->cellVolume(i, j, k);

    real dzf_m;
    if (self->isTop(k)) dzf_m = 0.5 * self->dzCell(k);
    else dzf_m = 0.5 * (self->dzCell(k) + self->dzCell(k - 1));
    Gzm(c) = -Kz(km) * dtg * r_rhozp(km) * r_visczp(km) / (self->dzCell(k) * dzf_m);
    Gzm(c) *= self->cellVolume(i, j, k);

    Gct(c) = (ch(c) + soil.Ss * wcn(c) / soil.theta_s) * r_rho(c);
    Gct(c) *= self->cellVolume(i, j, k);
    Gct(c) -= (Gxp(c) + Gxm(c) + Gyp(c) + Gym(c));

    if (self->isBottom(k) && bc4 != 1) {
      // no zp term in ct
    } else {
      Gct(c) -= Gzp(c);
    }
    if (self->isTop(k)) {
      if (sim_sw) {
        // shallow-water-coupled top: no zm term in ct
      } else if (bc5 == 1) {
        Gct(c) -= Gzm(c);
      }
    } else {
      Gct(c) -= Gzm(c);
    }
  });
}

void ReSolver::groundwaterRhs() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const int nz = grid_.nz();
  const real dtg = params_.dt;
  const Grid g = grid_;
  const ReSolver* self = this;
  const int bc4 = params_.bc_type[4], bc5 = params_.bc_type[5];
  const bool sim_sw = params_.sim_shallowwater;
  const real htop = params_.htop;
  const real qtop_uniform = params_.qtop;
  const real qbot_uniform = params_.qbot;
  const bool has_qtop = has_qtop_field_;
  // Lateral Dirichlet folding happens only at the GLOBAL domain boundary; interior rank
  // boundaries carry the +/-x,y coupling as real off-diagonal matrix entries (see
  // assembleAndSolve), so folding there would double-count. Serial: global == local.
  const int gnx = mc_ ? mc_->globalNx() : nx;
  const int gny = mc_ ? mc_->globalNy() : ny;
  const int gi0 = mc_ ? mc_->i0() : 0;
  const int gj0 = mc_ ? mc_->j0() : 0;
  RealArr1DHost Grhs = fields_.Grhs, hn = fields_.hn, ch = fields_.ch, wcn = fields_.wcn;
  RealArr1DHost Kz = fields_.Kz, r_rho = fields_.r_rho, r_rhozp = fields_.r_rhozp;
  RealArr1DHost r_visczp = fields_.r_visczp;
  RealArr1DHost Gxp = fields_.Gxp, Gxm = fields_.Gxm, Gyp = fields_.Gyp, Gym = fields_.Gym;
  RealArr1DHost Gzp = fields_.Gzp, Gzm = fields_.Gzm;
  RealArr1DHost qtop_field = qtop_field_;

  parallelForVolume("re_rhs", nx, ny, nz, KOKKOS_LAMBDA(int i, int j, int k) {
    const int c = g.getIndex(i, j, k);
    if (!self->active(i, j, k)) {
      Grhs(c) = hn(c);
      return;
    }
    const SoilParams& soil = self->soilAt(i, j, k);
    const int im = g.getIndex(i - 1, j, k);
    const int jm = g.getIndex(i, j - 1, k);
    const int km = g.getIndex(i, j, k - 1);
    const int ip = g.getIndex(i + 1, j, k);
    const int jp = g.getIndex(i, j + 1, k);
    const int kp = g.getIndex(i, j, k + 1);
    const real V = self->cellVolume(i, j, k);
    const real dzk = self->dzCell(k);

    Grhs(c) = (ch(c) + soil.Ss * wcn(c) / soil.theta_s) * hn(c) * r_rho(c) * V;
    Grhs(c) -= dtg * V * Kz(c) * r_rhozp(c) * r_rhozp(c) * r_visczp(c) / dzk;
    Grhs(c) += dtg * V * Kz(km) * r_rhozp(km) * r_rhozp(km) * r_visczp(km) / dzk;

    if (gi0 + i == gnx - 1) Grhs(c) -= Gxp(c) * hn(ip);
    if (gi0 + i == 0) Grhs(c) -= Gxm(c) * hn(im);
    if (gj0 + j == gny - 1) Grhs(c) -= Gyp(c) * hn(jp);
    if (gj0 + j == 0) Grhs(c) -= Gym(c) * hn(jm);

    if (self->isBottom(k)) {
      if (bc4 == 1) Grhs(c) += Gzp(c) * hn(kp);
      // Fixed-flux bottom (bctype_GW[4]==2 / qbot): mirror of the top fixed-flux term. Replace
      // the +z (down) bottom-face gravity term (subtracted above via Kz(c)) with the prescribed
      // flux qbot. Sign convention: qbot>0 = upward inflow INTO the domain (adds water), matching
      // qz>0 = upward flux from below. Validated by the bottom-flux mass-balance/stability test.
      if (bc4 == 2) {
        Grhs(c) += dtg * r_rhozp(c) * V *
                   (qbot_uniform + Kz(c) * r_rhozp(c) * r_visczp(c)) / dzk;
      }
    }
    if (self->isTop(k) && !sim_sw && bc5 == 1) {
      Grhs(c) -= Gzm(c) * htop;
    }
    // Fixed-flux top (bctype_GW[5]==2, legacy groundwater.c:666-669): replace the top-face
    // gravity term (already added above via Kz(km)) with the prescribed infiltration flux qtop.
    // qtop<0 = downward infiltration. The per-column field (partial-width recharge) overrides the
    // uniform scalar when set; columns outside the recharge region carry qtop=0 (=> no-flux top).
    if (self->isTop(k) && !sim_sw && bc5 == 2) {
      const int km = g.getIndex(i, j, -1);
      const real qcol = has_qtop ? qtop_field(g.getSurfaceIndex(i, j)) : qtop_uniform;
      Grhs(c) += -dtg * r_rhozp(km) * V *
                 (qcol + Kz(km) * r_rhozp(km) * r_visczp(km)) / dzk;
    }
  });
}

void ReSolver::attachSolver(LinearSolver& solver, Decomp3D& dd) {
  solver_ = &solver;
  dd_ = &dd;
  sys_ = solver.createSystem(dd);
  const int nloc = dd.ownedRowCount();
  b_owned_ = RealArr1D("re_b", static_cast<size_t>(nloc));
  x_owned_ = RealArr1D("re_x", static_cast<size_t>(nloc));
}

void ReSolver::assembleAndSolve() {
  perf::ScopedTimer t_asm(perf_, perf::Region::GwAssembly);
  const int gnx = grid_.nx();
  const int gny = grid_.ny();
  const int gnz = grid_.nz();
  const auto corners = dd_->localCorners();
  const int i0 = std::get<0>(corners), j0 = std::get<1>(corners);
  const int k0 = std::get<2>(corners);
  const int ni = std::get<3>(corners), nj = std::get<4>(corners), nk = std::get<5>(corners);
  const int rstart = dd_->ownershipRange().first;
  auto& f = fields_;

  sys_->zeroEntries();
  sys_->beginAssembly();
  RealArr1DHost b_host = b_owned_;
  int64_t bytes = 0;
  for (int gk = k0; gk < k0 + nk; ++gk)
    for (int gj = j0; gj < j0 + nj; ++gj)
      for (int gi = i0; gi < i0 + ni; ++gi) {
        const int lc = s(gi - i0, gj - j0, gk - k0);
        const int grow = dd_->globalRow(gi, gj, gk);
        int cols[kMaxStencil];
        int ncols = 0;
        dd_->stencilColumns(gi, gj, gk, cols, ncols);
        const int rowW = (gi - 1 >= 0) ? dd_->globalRow(gi - 1, gj, gk) : -1;
        const int rowE = (gi + 1 < gnx) ? dd_->globalRow(gi + 1, gj, gk) : -1;
        const int rowS = (gj - 1 >= 0) ? dd_->globalRow(gi, gj - 1, gk) : -1;
        const int rowN = (gj + 1 < gny) ? dd_->globalRow(gi, gj + 1, gk) : -1;
        const int rowB = (gk - 1 >= 0) ? dd_->globalRow(gi, gj, gk - 1) : -1;
        const int rowT = (gk + 1 < gnz) ? dd_->globalRow(gi, gj, gk + 1) : -1;
        real vals[kMaxStencil];
        for (int cc = 0; cc < ncols; ++cc) {
          const int col = cols[cc];
          if (col == grow)
            vals[cc] = f.Gct(lc);
          else if (col == rowW)
            vals[cc] = f.Gxm(lc);
          else if (col == rowE)
            vals[cc] = f.Gxp(lc);
          else if (col == rowS)
            vals[cc] = f.Gym(lc);
          else if (col == rowN)
            vals[cc] = f.Gyp(lc);
          else if (col == rowB)
            vals[cc] = f.Gzm(lc);
          else if (col == rowT)
            vals[cc] = f.Gzp(lc);
          else
            vals[cc] = 0.0;
        }
        sys_->addRow(grow, ncols, cols, vals);
        b_host(static_cast<size_t>(grow - rstart)) = f.Grhs(lc);
        bytes += static_cast<int64_t>(ncols) * (2 * static_cast<int>(sizeof(int)) + static_cast<int>(sizeof(real)));
      }
  if (b_host.data() != b_owned_.data()) Kokkos::deep_copy(b_owned_, b_host);
  if (perf_) perf_->counters().addBytes(bytes + static_cast<int64_t>(ni) * nj * nk * static_cast<int>(sizeof(real)));
  sys_->endAssembly();

  {
    perf::ScopedTimer t_ksp(perf_, perf::Region::GwKsp);
    solver_->solve(*sys_, b_owned_, x_owned_);
    if (perf_) perf_->counters().addKspIters(solver_->getIterationCount());
  }

  RealArr1DHost x_host = x_owned_;
  const Grid g = grid_;
  Decomp3D* dd = dd_;
  RealArr1DHost h = f.h;
  const int rs = rstart;
  parallelForVolume("re_solve_unpack", ni, nj, nk, KOKKOS_LAMBDA(int li, int lj, int lk) {
    const int gi = i0 + li;
    const int gj = j0 + lj;
    const int gk = k0 + lk;
    const int grow = dd->globalRow(gi, gj, gk);
    h(g.getIndex(li, lj, lk)) = x_host(static_cast<size_t>(grow - rs));
  });
}

void ReSolver::enforceHeadBc() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const int nz = grid_.nz();
  auto& f = fields_;
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      f.h(s(-1, j, k)) = f.h(s(0, j, k));
      f.h(s(nx, j, k)) = f.h(s(nx - 1, j, k));
    }
  }
  for (int k = 0; k < nz; ++k) {
    for (int i = 0; i < nx; ++i) {
      f.h(s(i, -1, k)) = f.h(s(i, 0, k));
      f.h(s(i, ny, k)) = f.h(s(i, ny - 1, k));
    }
  }
  if (params_.bc_type[4] == 1) {
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) f.h(s(i, j, nz)) = params_.hbot;
  } else {
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) f.h(s(i, j, nz)) = f.h(s(i, j, nz - 1));
  }
  if (!params_.sim_shallowwater && params_.bc_type[5] == 1) {
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) f.h(s(i, j, -1)) = params_.htop;
  } else {
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) f.h(s(i, j, -1)) = f.h(s(i, j, 0));
  }
}

void ReSolver::groundwaterFlux() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const int nz = grid_.nz();
  const Grid g = grid_;
  const ReSolver* self = this;
  const ReParams p = params_;
  const GwFields f = fields_;
  const bool full3d = params_.use_full3d;
  const real dx = params_.dx, dy = params_.dy;
  RealArr1DHost qx = fields_.qx, qy = fields_.qy, qz = fields_.qz;
  RealArr1DHost h = fields_.h, Kx = fields_.Kx, Ky = fields_.Ky;
  RealArr1DHost r_viscxp = fields_.r_viscxp, r_viscyp = fields_.r_viscyp;
  parallelForVolume("re_flux", nx, ny, nz, KOKKOS_LAMBDA(int i, int j, int k) {
    if (!self->active(i, j, k)) return;
    const int c = g.getIndex(i, j, k);
    if (full3d) {
      // Lateral Darcy flux through the +x and +y faces of cell c (no gravity term laterally;
      // face areas Ax = dy*dz, Ay = dx*dz). Kx(c)/Ky(c) are the face means from computeKFace;
      // closed global walls have Kx=Ky=0 so wall flux is exactly zero. Legacy darcy_flux
      // "x"/"y","none" with cos=1, sin=0, rho=0 for the flat (follow_terrain=false) case.
      const real dzk = self->dzCell(k);
      const int ip = g.getIndex(i + 1, j, k);
      const int jp = g.getIndex(i, j + 1, k);
      qx(c) = (dy * dzk) * Kx(c) * r_viscxp(c) * (h(ip) - h(c)) / dx;
      qy(c) = (dx * dzk) * Ky(c) * r_viscyp(c) * (h(jp) - h(c)) / dy;
    } else {
      qx(c) = 0.0;
      qy(c) = 0.0;
    }
    qz(c) = darcyFluxZ(p, f, self->soilAt(i, j, k), self->soilAt(i, j, k + 1), i, j, k, g, false);
  });
  const int bc5 = params_.bc_type[5];
  const int bc4 = params_.bc_type[4];
  const bool has_qtop = has_qtop_field_;
  const real qtop_uniform = params_.qtop;
  const real qbot_uniform = params_.qbot;
  const real Az = dx * dy;
  const int nzc = nz;
  RealArr1DHost qtop_field = qtop_field_;
  parallelForSurface("re_flux_top", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    if (!self->isTop(0)) return;
    const int km = g.getIndex(i, j, -1);
    if (bc5 == 2) {
      // Fixed-flux top: the top-face flux IS the prescribed recharge (Neumann BC). Sign: the
      // corrector adds water for qz(km)<0 (cf. the b2-gw bc5==1 path), and legacy qtop<0 is
      // downward infiltration, so qz(km) = qtop*Az gives the correct mass addition. This makes
      // the PCA corrector consistent with the predictor's qtop term (legacy's commented `dqep`).
      const real qcol = has_qtop ? qtop_field(g.getSurfaceIndex(i, j)) : qtop_uniform;
      qz(km) = qcol * Az;
    } else {
      const SoilParams& st = self->soilAt(i, j, 0);
      qz(km) = darcyFluxZ(p, f, st, st, i, j, 0, g, true);
    }
  });
  // Fixed-flux bottom (bctype_GW[4]==2 / qbot): override the deepest cell's +z (down) face flux
  // with the prescribed value so the corrector mass change matches the predictor's qbot term.
  // qbot>0 = upward inflow into the domain (qz>0 = upward), so qz(bottomFace) = qbot*Az adds water.
  if (bc4 == 2) {
    parallelForSurface("re_flux_bot", nx, ny, KOKKOS_LAMBDA(int i, int j) {
      const int c = g.getIndex(i, j, nzc - 1);
      qz(c) = qbot_uniform * Az;
    });
  }
  // 3D MPI: the -x/-y ghost back-faces (used as the inflow face by updateWaterContent's
  // divergence) must hold the neighbor's shared +face flux. Exchanging qx/qy fills
  // qx(-1)/qy(-1) from the neighbor's owned last-cell +face, making the divergence conservative
  // across ranks. z is on-rank. Gated behind use_full3d (vertical path is unchanged).
  if (full3d && mc_ != nullptr && mc_->size() > 1) {
    mc_->haloExchange3D(fields_.qx, nx + 2, ny + 2, nz);
    mc_->haloExchange3D(fields_.qy, nx + 2, ny + 2, nz);
  }
}

void ReSolver::checkRoom() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const int nz = grid_.nz();
  const Grid g = grid_;
  const ReSolver* self = this;
  const real dx = params_.dx, dy = params_.dy;
  RealArr1DHost room = fields_.room, wc = fields_.wc;
  parallelForVolume("re_room", nx, ny, nz, KOKKOS_LAMBDA(int i, int j, int k) {
    const int c = g.getIndex(i, j, k);
    if (!self->active(i, j, k))
      room(c) = 0.0;
    else
      room(c) = (self->soilAt(i, j, k).theta_s - wc(c)) * self->dzCell(k) * dx * dy;
  });
}

void ReSolver::updateWaterContent() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const int nz = grid_.nz();
  const real dtg = params_.dt;
  const Grid g = grid_;
  const ReSolver* self = this;
  RealArr1DHost wc = fields_.wc, wcn = fields_.wcn, h = fields_.h, hn = fields_.hn;
  RealArr1DHost qx = fields_.qx, qy = fields_.qy, qz = fields_.qz, r_rho = fields_.r_rho;
  RealArr1DHost r_rhoxp = fields_.r_rhoxp, r_rhoyp = fields_.r_rhoyp, r_rhozp = fields_.r_rhozp;
  parallelForVolume("re_update_wc", nx, ny, nz, KOKKOS_LAMBDA(int i, int j, int k) {
    if (!self->active(i, j, k)) return;
    const SoilParams& soil = self->soilAt(i, j, k);
    const int c = g.getIndex(i, j, k);
    const real V = self->cellVolume(i, j, k);
    const real coeff = V * (r_rho(c) + r_rho(c) * soil.Ss * (h(c) - hn(c)) / soil.theta_s);
    const int imc = g.getIndex(i - 1, j, k);
    const int jmc = g.getIndex(i, j - 1, k);
    const int kmc = g.getIndex(i, j, k - 1);
    const real dqx = dtg * (qx(c) * r_rhoxp(c) - qx(imc) * r_rhoxp(imc));
    const real dqy = dtg * (qy(c) * r_rhoyp(c) - qy(jmc) * r_rhoyp(jmc));
    const real dqz = dtg * (qz(c) * r_rhozp(c) - qz(kmc) * r_rhozp(kmc));
    wc(c) = (wcn(c) * r_rho(c) * V + dqx + dqy + dqz) / coeff;
  });
}

void ReSolver::reallocateWaterContent() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const int nz = grid_.nz();
  auto& f = fields_;

  // 3D MPI: the lateral saturated-neighbor checks below read the +/-x,y ghost wc, so refresh
  // the wc halo first. Gated behind use_full3d (the vertical path does no lateral check).
  if (params_.use_full3d && mc_ != nullptr && mc_->size() > 1)
    mc_->haloExchange3D(f.wc, nx + 2, ny + 2, nz);

  // "Saturated neighbor" compares each neighbor's water content to that neighbor's own theta_s
  // (P13 per-column / P23 per-cell). For uniform soil every cell shares theta_s, reproducing the
  // P5 test bit-for-bit.
  auto adjacentToSat = [&](int i, int j, int k) {
    const int c = s(i, j, k);
    if (f.wc(s(i, j, k + 1)) >= soilAt(i, j, k + 1).theta_s) return true;
    if (f.Kx(c) != 0.0 && f.wc(s(i + 1, j, k)) >= soilAt(i + 1, j, k).theta_s) return true;
    if (f.Kx(s(i - 1, j, k)) != 0.0 && f.wc(s(i - 1, j, k)) >= soilAt(i - 1, j, k).theta_s)
      return true;
    if (f.Ky(c) != 0.0 && f.wc(s(i, j + 1, k)) >= soilAt(i, j + 1, k).theta_s) return true;
    if (f.Ky(s(i, j - 1, k)) != 0.0 && f.wc(s(i, j - 1, k)) >= soilAt(i, j - 1, k).theta_s)
      return true;
    if (isTop(k)) {
      if (!params_.sim_shallowwater && params_.bc_type[5] == 1 && params_.htop >= 0.0) return true;
    } else if (f.wc(s(i, j, k - 1)) >= soilAt(i, j, k - 1).theta_s) {
      return true;
    }
    return false;
  };

  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        if (!active(i, j, k)) continue;
        const SoilParams& soil = soilAt(i, j, k);
        const int c = s(i, j, k);
        f.wch(c) = VanGenuchten::waterContentFromHead(soil, f.h(c));
        f.hwc(c) = VanGenuchten::headFromWaterContent(soil, f.wc(c));
        if (f.wc(c) >= soil.theta_s) {
          f.wc(c) = soil.theta_s;
        } else if (!adjacentToSat(i, j, k)) {
          if (f.wc(c) < 0.9999 * soil.theta_s) f.h(c) = f.hwc(c);
        } else {
          f.wc(c) = f.wch(c);
        }
      }
    }
  }
}

void ReSolver::clampWaterContent() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const int nz = grid_.nz();
  const Grid g = grid_;
  const ReSolver* self = this;
  RealArr1DHost wc = fields_.wc;
  parallelForVolume("re_clamp_wc", nx, ny, nz, KOKKOS_LAMBDA(int i, int j, int k) {
    if (!self->active(i, j, k)) return;
    const SoilParams& soil = self->soilAt(i, j, k);
    const int c = g.getIndex(i, j, k);
    if (wc(c) > soil.theta_s) wc(c) = soil.theta_s;
    else if (wc(c) < soil.theta_r) wc(c) = soil.theta_r;
  });
}

void ReSolver::enforceMoistureBc() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const int nz = grid_.nz();
  auto& f = fields_;
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      f.wc(s(-1, j, k)) = f.wc(s(0, j, k));
      f.wc(s(nx, j, k)) = f.wc(s(nx - 1, j, k));
    }
  }
  for (int k = 0; k < nz; ++k) {
    for (int i = 0; i < nx; ++i) {
      f.wc(s(i, -1, k)) = f.wc(s(i, 0, k));
      f.wc(s(i, ny, k)) = f.wc(s(i, ny - 1, k));
    }
  }
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      f.wc(s(i, j, nz)) = f.wc(s(i, j, nz - 1));
      f.wc(s(i, j, -1)) = f.wc(s(i, j, 0));
    }
  }
}

void ReSolver::adaptiveTimeStep() {
  const int nx = grid_.nx();
  const int ny = grid_.ny();
  const int nz = grid_.nz();
  auto& f = fields_;
  const real r_red = 0.75;
  const real r_inc = 1.25;
  real dq_max = 0.0;
  real dt_Comin = 1.0e8;

  params_.dt = dtn_;
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        if (!active(i, j, k) || isTop(k)) continue;
        const SoilParams& soil = soilAt(i, j, k);
        const int c = s(i, j, k);
        // Face areas perpendicular to each axis: Ax = dy*dz, Ay = dx*dz, Az = dx*dy.
        const real Ax = params_.dy * dzCell(k);
        const real Ay = params_.dx * dzCell(k);
        const real Az = params_.dx * params_.dy;
        const real qin = f.qx(c) / Ax + f.qy(c) / Ay + f.qz(c) / Az;
        const real qou = f.qx(s(i - 1, j, k)) / Ax + f.qy(s(i, j - 1, k)) / Ay +
                         f.qz(s(i, j, k - 1)) / Az;
        const real dq = std::fabs(qin - qou) * params_.dt / dzCell(k);
        dq_max = std::max(dq_max, dq);
        if (f.wc(c) < soil.theta_s) {
          const real dKdwc = VanGenuchten::dKdwc(soil, f.wc(c), soil.Ks_z);
          const real dt_Co = params_.co_max * dzCell(k) / dKdwc;
          dt_Comin = std::min(dt_Comin, dt_Co);
        }
      }
    }
  }
  // The dt criteria are global: reduce the flux-imbalance max and the Courant min across
  // all ranks so every rank advances with the same dt (rank-count equivalence).
  if (mc_ != nullptr && mc_->size() > 1) {
    real dq_local = dq_max, dtco_local = dt_Comin;
    MPI_Allreduce(&dq_local, &dq_max, 1, MPI_DOUBLE, MPI_MAX, mc_->comm());
    MPI_Allreduce(&dtco_local, &dt_Comin, 1, MPI_DOUBLE, MPI_MIN, mc_->comm());
  }
  if (dq_max > 0.02)
    params_.dt *= r_red;
  else if (dq_max < 0.01)
    params_.dt *= r_inc;
  if (params_.dt > dt_Comin) params_.dt = dt_Comin;
  if (params_.dt > params_.dt_max) params_.dt = params_.dt_max;
  if (params_.dt < params_.dt_min) params_.dt = params_.dt_min;
  dtn_ = params_.dt;
}

void ReSolver::advanceStep() {
  dt_used_ = params_.dt;
  if (perf_) perf_->counters().addCells(static_cast<int64_t>(grid_.nx()) * grid_.ny() * grid_.nz());
  {
    perf::ScopedTimer t(perf_, perf::Region::GwUpdate);
    savePreviousStepFields();
    updateDiagnostics();
    computeKFace();
    baroclinicFace();
    groundwaterMatCoeff();
    groundwaterRhs();
  }
  assembleAndSolve();
  {
    perf::ScopedTimer t(perf_, perf::Region::GwUpdate);
    enforceHeadBc();
    computeKFace();
    groundwaterFlux();
    if (params_.use_corrector) {
      checkRoom();
      updateWaterContent();
    } else {
      const int nx = grid_.nx();
      const int ny = grid_.ny();
      const int nz = grid_.nz();
      for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
          for (int i = 0; i < nx; ++i)
            if (active(i, j, k))
              fields_.wc(s(i, j, k)) =
                  VanGenuchten::waterContentFromHead(soilAt(i, j, k), fields_.h(s(i, j, k)));
    }
    updateDiagnostics();
    reallocateWaterContent();
    clampWaterContent();
    enforceMoistureBc();
    if (params_.adaptive_dt) adaptiveTimeStep();
  }
}

}  // namespace frehg2
