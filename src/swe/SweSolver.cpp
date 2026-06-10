#include "swe/SweSolver.hpp"

#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

#ifdef USE_PETSC
#include <petsc.h>
#endif

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace frehg2 {

SweSolver::SweSolver(Grid grid, SweParameters parameters)
    : grid_(std::move(grid)),
      parameters_(parameters)
{
}

const Grid& SweSolver::grid() const noexcept
{
    return grid_;
}

const SweParameters& SweSolver::parameters() const noexcept
{
    return parameters_;
}

void SweSolver::initialize(const Domain& domain, State& state, real initial_eta) const
{
#ifdef USE_KOKKOS
    const auto bed_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), domain.z);
    auto h_h = Kokkos::create_mirror_view(state.h_old);
    auto h_new_h = Kokkos::create_mirror_view(state.h_new);
    auto z_h = Kokkos::create_mirror_view(state.z);

    for (index_t idx = 0; idx < grid_.nSurfaceCellMem(); ++idx) {
        const real depth = std::max<real>(0.0, initial_eta - bed_h(idx));
        h_h(idx) = depth;
        h_new_h(idx) = depth;
        z_h(idx) = bed_h(idx);
    }

    Kokkos::deep_copy(state.h_old, h_h);
    Kokkos::deep_copy(state.h_new, h_new_h);
    Kokkos::deep_copy(state.z, z_h);
    Kokkos::deep_copy(state.hu_old, 0.0);
    Kokkos::deep_copy(state.hu_new, 0.0);
    Kokkos::deep_copy(state.hv_old, 0.0);
    Kokkos::deep_copy(state.hv_new, 0.0);
#else
    (void)domain;
    (void)state;
    (void)initial_eta;
#endif
}

SweLinearSystem SweSolver::assembleLinearSystem(
    const std::vector<real>& eta,
    const std::vector<real>& bed,
    const std::vector<int>& active_mask,
    const std::vector<real>& ex,
    const std::vector<real>& ey,
    const std::vector<real>& dx_factor,
    const std::vector<real>& dy_factor) const
{
    const auto n2ci = static_cast<std::size_t>(grid_.nSurfaceCell());
    const auto n2ct = static_cast<std::size_t>(grid_.nSurfaceCellMem());
    if (eta.size() < n2ct || bed.size() < n2ct || active_mask.size() < n2ct ||
        ex.size() < n2ct || ey.size() < n2ct || dx_factor.size() < n2ct ||
        dy_factor.size() < n2ct) {
        throw std::invalid_argument("SWE input arrays must cover surface memory cells");
    }

    SweLinearSystem system;
    system.dom2mat.assign(n2ci, -1);
    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        if (active_mask[idx] > 0) {
            system.dom2mat[idx] = system.n;
            system.mat2dom.push_back(static_cast<int>(idx));
            ++system.n;
        }
    }
    system.rhs.assign(static_cast<std::size_t>(system.n), 0.0);

    const real area = grid_.dx() * grid_.dy();
    const real coef = parameters_.gravity * parameters_.dt * parameters_.dt;

    for (int row = 0; row < system.n; ++row) {
        const int idx = system.mat2dom[static_cast<std::size_t>(row)];
        const auto logical = grid_.getSurfaceLogicalIndex(static_cast<index_t>(idx));
        const int i = logical[0];
        const int j = logical[1];

        const int im = xm(i, j);
        const int jm = ym(i, j);
        const int ip = xp(i, j);
        const int jp = yp(i, j);

        real sxp = coef * area * dx_factor[static_cast<std::size_t>(idx)];
        real sxm = coef * area * dx_factor[static_cast<std::size_t>(im)];
        real syp = coef * area * dy_factor[static_cast<std::size_t>(idx)];
        real sym = coef * area * dy_factor[static_cast<std::size_t>(jm)];
        real rhs = eta[static_cast<std::size_t>(idx)] * area -
                   parameters_.dt *
                       (grid_.dy() * ex[static_cast<std::size_t>(idx)] -
                        grid_.dy() * ex[static_cast<std::size_t>(im)] +
                        grid_.dx() * ey[static_cast<std::size_t>(idx)] -
                        grid_.dx() * ey[static_cast<std::size_t>(jm)]);
        real sct = area + sxp + sxm + syp + sym;

        if (eta[static_cast<std::size_t>(idx)] <= bed[static_cast<std::size_t>(idx)]) {
            sct = area;
            rhs = eta[static_cast<std::size_t>(idx)] * area;
            sxp = 0.0;
            sxm = 0.0;
            syp = 0.0;
            sym = 0.0;
        }

        if (i == 0) {
            sct -= sxm;
        }
        if (i == grid_.nx() - 1) {
            sct -= sxp;
        }
        if (j == 0) {
            sct -= sym;
        }
        if (j == grid_.ny() - 1) {
            sct -= syp;
        }

        auto addEntry = [&system, row](int domain_col, real value) {
            if (domain_col < 0) {
                return;
            }
            const auto col_domain = static_cast<std::size_t>(domain_col);
            if (col_domain >= system.dom2mat.size()) {
                return;
            }
            const int col = system.dom2mat[col_domain];
            if (col >= 0) {
                system.entries.push_back(SweMatrixEntry{row, col, value});
            }
        };

        if (j > 0) {
            addEntry(jm, -sym);
        }
        if (i > 0) {
            addEntry(im, -sxm);
        }
        addEntry(idx, sct);
        if (i < grid_.nx() - 1) {
            addEntry(ip, -sxp);
        }
        if (j < grid_.ny() - 1) {
            addEntry(jp, -syp);
        }

        system.rhs[static_cast<std::size_t>(row)] = rhs;
    }

    return system;
}

std::vector<real> SweSolver::solveLinearSystem(const SweLinearSystem& system) const
{
    if (system.n <= 0) {
        return {};
    }

#ifdef USE_PETSC
    Mat matrix = nullptr;
    Vec rhs = nullptr;
    Vec solution = nullptr;
    KSP solver = nullptr;

    MatCreateSeqAIJ(PETSC_COMM_SELF, system.n, system.n, 5, nullptr, &matrix);
    for (const auto& entry : system.entries) {
        MatSetValue(matrix, entry.row, entry.col, entry.value, ADD_VALUES);
    }
    MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);

    VecCreateSeq(PETSC_COMM_SELF, system.n, &rhs);
    VecCreateSeq(PETSC_COMM_SELF, system.n, &solution);
    for (int i = 0; i < system.n; ++i) {
        VecSetValue(rhs, i, system.rhs[static_cast<std::size_t>(i)], INSERT_VALUES);
    }
    VecAssemblyBegin(rhs);
    VecAssemblyEnd(rhs);

    KSPCreate(PETSC_COMM_SELF, &solver);
    KSPSetOperators(solver, matrix, matrix);
    KSPSetType(solver, KSPCG);
    PC pc = nullptr;
    KSPGetPC(solver, &pc);
    PCSetType(pc, PCSOR);
    PCSORSetSymmetric(pc, SOR_SYMMETRIC_SWEEP);
    KSPSetTolerances(solver, 1.0e-8, PETSC_DEFAULT, PETSC_DEFAULT, 10000000);
    KSPSetFromOptions(solver);
    KSPSolve(solver, rhs, solution);

    const PetscScalar* raw_solution = nullptr;
    VecGetArrayRead(solution, &raw_solution);
    std::vector<real> values(static_cast<std::size_t>(system.n), 0.0);
    for (int i = 0; i < system.n; ++i) {
        values[static_cast<std::size_t>(i)] = PetscRealPart(raw_solution[i]);
    }
    VecRestoreArrayRead(solution, &raw_solution);

    KSPDestroy(&solver);
    VecDestroy(&solution);
    VecDestroy(&rhs);
    MatDestroy(&matrix);
    return values;
#else
    throw std::runtime_error("PETSc support is required for SWE linear solves");
#endif
}

void SweSolver::applyRainEvap(std::vector<real>& eta, const std::vector<real>& bed) const
{
    if (eta.size() < grid_.nSurfaceCell() || bed.size() < grid_.nSurfaceCell()) {
        throw std::invalid_argument("rain/evap arrays do not cover interior cells");
    }
    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        eta[idx] += (parameters_.rain_rate - parameters_.evap_rate) * parameters_.dt;
        if (eta[idx] - bed[idx] < parameters_.min_depth) {
            eta[idx] = bed[idx];
        }
    }
}

real SweSolver::computeCflDiagnostic(
    const std::vector<real>& h,
    const std::vector<real>& u,
    const std::vector<real>& v) const
{
    real max_speed = 0.0;
    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        const real wave_speed = std::sqrt(parameters_.gravity * std::max<real>(h[idx], 0.0));
        max_speed = std::max(max_speed, std::abs(u[idx]) + wave_speed);
        max_speed = std::max(max_speed, std::abs(v[idx]) + wave_speed);
    }
    if (max_speed == 0.0) {
        return 0.0;
    }
    const real length = std::min(grid_.dx(), grid_.dy());
    return parameters_.dt * max_speed / length;
}

void SweSolver::updateVelocity(
    std::vector<real>& u,
    std::vector<real>& v,
    const std::vector<real>& eta_old,
    const std::vector<real>& eta_new) const
{
    (void)eta_old;
    for (int j = 0; j < grid_.ny(); ++j) {
        for (int i = 0; i < grid_.nx(); ++i) {
            const auto idx = grid_.getSurfaceIndex(i, j);
            if (i < grid_.nx() - 1) {
                u[idx] -= parameters_.dt * parameters_.gravity *
                          (eta_new[xp(i, j)] - eta_new[idx]) / grid_.dx();
            }
            if (j < grid_.ny() - 1) {
                v[idx] -= parameters_.dt * parameters_.gravity *
                          (eta_new[yp(i, j)] - eta_new[idx]) / grid_.dy();
            }
        }
    }
}

LegacySweState SweSolver::initializeLegacyState(const std::vector<real>& bed, real initial_eta) const
{
    const auto n2ci = static_cast<std::size_t>(grid_.nSurfaceCell());
    const auto n2ct = static_cast<std::size_t>(grid_.nSurfaceCellMem());
    if (bed.size() != n2ci) {
        throw std::invalid_argument("legacy SWE bed size must match the interior surface grid");
    }

    LegacySweState state;
    state.bottom.assign(n2ct, 0.0);
    state.bottom_xp.assign(n2ct, 0.0);
    state.bottom_yp.assign(n2ct, 0.0);
    state.active.assign(n2ct, 0);

    state.eta.assign(n2ct, 0.0);
    state.etan.assign(n2ct, 0.0);
    state.dept.assign(n2ct, 0.0);
    state.deptx.assign(n2ct, 0.0);
    state.depty.assign(n2ct, 0.0);

    state.uu.assign(n2ct, 0.0);
    state.un.assign(n2ct, 0.0);
    state.uy.assign(n2ct, 0.0);
    state.vv.assign(n2ct, 0.0);
    state.vn.assign(n2ct, 0.0);
    state.vx.assign(n2ct, 0.0);

    state.Fu.assign(n2ct, 0.0);
    state.Fv.assign(n2ct, 0.0);
    state.Ex.assign(n2ct, 0.0);
    state.Ey.assign(n2ct, 0.0);
    state.Dx.assign(n2ct, 0.0);
    state.Dy.assign(n2ct, 0.0);
    state.CDx.assign(n2ct, 0.0);
    state.CDy.assign(n2ct, 0.0);

    state.Vs.assign(n2ct, 0.0);
    state.Vsn.assign(n2ct, 0.0);
    state.Vflux.assign(n2ci, 0.0);
    state.Vsx.assign(n2ct, 0.0);
    state.Vsy.assign(n2ct, 0.0);
    state.Asz.assign(n2ct, 0.0);
    state.Aszx.assign(n2ct, 0.0);
    state.Aszy.assign(n2ct, 0.0);
    state.Asx.assign(n2ct, 0.0);
    state.Asy.assign(n2ct, 0.0);

    state.Sxp.assign(n2ci, 0.0);
    state.Sxm.assign(n2ci, 0.0);
    state.Syp.assign(n2ci, 0.0);
    state.Sym.assign(n2ci, 0.0);
    state.Sct.assign(n2ci, 0.0);
    state.Srhs.assign(n2ci, 0.0);

    state.wtfx.assign(n2ct, 0.0);
    state.wtfy.assign(n2ct, 0.0);
    state.cflx.assign(n2ci, 0.0);
    state.cfly.assign(n2ci, 0.0);
    state.cfl_active.assign(n2ci, 0.0);

    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        state.bottom[idx] = bed[idx];
        state.eta[idx] = std::max(initial_eta, bed[idx]);
        state.active[idx] = 1;
    }
    fillLegacyGhosts(state);

    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        const auto ij = grid_.getSurfaceLogicalIndex(idx);
        state.bottom_xp[idx] = std::max(state.bottom[idx], state.bottom[static_cast<std::size_t>(xp(ij[0], ij[1]))]);
        state.bottom_yp[idx] = std::max(state.bottom[idx], state.bottom[static_cast<std::size_t>(yp(ij[0], ij[1]))]);
    }
    for (int i = 0; i < grid_.nx(); ++i) {
        const auto south = grid_.getSurfaceIndex(i, 0);
        const auto north = grid_.getSurfaceIndex(i, grid_.ny() - 1);
        const auto ym_ghost = grid_.surfaceGhostIndex(Direction::YM, i, 0);
        const auto yp_ghost = grid_.surfaceGhostIndex(Direction::YP, i, grid_.ny() - 1);
        state.bottom_yp[yp_ghost] = state.bottom_yp[north];
        state.bottom_yp[ym_ghost] = std::max(state.bottom[ym_ghost], state.bottom[south]);
    }
    for (int j = 0; j < grid_.ny(); ++j) {
        const auto west = grid_.getSurfaceIndex(0, j);
        const auto east = grid_.getSurfaceIndex(grid_.nx() - 1, j);
        const auto xm_ghost = grid_.surfaceGhostIndex(Direction::XM, 0, j);
        const auto xp_ghost = grid_.surfaceGhostIndex(Direction::XP, grid_.nx() - 1, j);
        state.bottom_xp[xp_ghost] = state.bottom_xp[east];
        state.bottom_xp[xm_ghost] = std::max(state.bottom[xm_ghost], state.bottom[west]);
    }

    updateLegacyDepth(state);
    updateLegacySubgridVariables(state, false);
    updateLegacyDragCoefficient(state);
    return state;
}

void SweSolver::advanceLegacyStep(LegacySweState& state, real rain_rate, real evap_rate) const
{
    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        state.etan[idx] = state.eta[idx];
    }
    enforceLegacySurfaceBc(state);
    legacyMomentumSource(state);
    enforceLegacyReflectiveBoundaryFluxes(state);
    auto system = assembleLegacyLinearSystem(state);
    const auto solution = solveLinearSystem(system);
    for (int row = 0; row < system.n; ++row) {
        state.eta[static_cast<std::size_t>(system.mat2dom[static_cast<std::size_t>(row)])] =
            solution[static_cast<std::size_t>(row)];
    }
    enforceLegacySurfaceBc(state);
    legacyCflLimiter(state, rain_rate);
    legacyEvapRain(state, rain_rate, evap_rate);
    updateLegacyDepth(state);
    updateLegacyDepth(state);
    updateLegacySubgridVariables(state, true);
    updateLegacyDragCoefficient(state);
    legacyWaterfallLocation(state);
    updateLegacyVelocity(state);
    enforceLegacyReflectiveBoundaryFluxes(state);
    applyLegacyFreeOutflowBc(state);
    enforceLegacyVelocityBc(state);
    interpolateLegacyVelocity(state);
}

void SweSolver::fillLegacyGhosts(LegacySweState& state) const
{
    for (int i = 0; i < grid_.nx(); ++i) {
        const auto south = grid_.getSurfaceIndex(i, 0);
        const auto north = grid_.getSurfaceIndex(i, grid_.ny() - 1);
        const auto ym_ghost = grid_.surfaceGhostIndex(Direction::YM, i, 0);
        const auto yp_ghost = grid_.surfaceGhostIndex(Direction::YP, i, grid_.ny() - 1);
        state.bottom[ym_ghost] = state.bottom[south];
        state.bottom[yp_ghost] = state.bottom[north];
        state.eta[ym_ghost] = state.eta[south];
        state.eta[yp_ghost] = state.eta[north];
    }
    for (int j = 0; j < grid_.ny(); ++j) {
        const auto west = grid_.getSurfaceIndex(0, j);
        const auto east = grid_.getSurfaceIndex(grid_.nx() - 1, j);
        const auto xm_ghost = grid_.surfaceGhostIndex(Direction::XM, 0, j);
        const auto xp_ghost = grid_.surfaceGhostIndex(Direction::XP, grid_.nx() - 1, j);
        state.bottom[xm_ghost] = state.bottom[west];
        state.bottom[xp_ghost] = state.bottom[east];
        state.eta[xm_ghost] = state.eta[west];
        state.eta[xp_ghost] = state.eta[east];
    }
}

void SweSolver::enforceLegacySurfaceBc(LegacySweState& state) const
{
    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        if (state.eta[idx] < state.bottom[idx]) {
            state.eta[idx] = state.bottom[idx];
        }
    }
    for (int i = 0; i < grid_.nx(); ++i) {
        const auto south = grid_.getSurfaceIndex(i, 0);
        const auto north = grid_.getSurfaceIndex(i, grid_.ny() - 1);
        state.eta[grid_.surfaceGhostIndex(Direction::YM, i, 0)] = state.eta[south];
        state.eta[grid_.surfaceGhostIndex(Direction::YP, i, grid_.ny() - 1)] = state.eta[north];
    }
    for (int j = 0; j < grid_.ny(); ++j) {
        const auto west = grid_.getSurfaceIndex(0, j);
        const auto east = grid_.getSurfaceIndex(grid_.nx() - 1, j);
        state.eta[grid_.surfaceGhostIndex(Direction::XM, 0, j)] = state.eta[west];
        state.eta[grid_.surfaceGhostIndex(Direction::XP, grid_.nx() - 1, j)] = state.eta[east];
    }
}

void SweSolver::updateLegacyDepth(LegacySweState& state) const
{
    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        const real diff = state.eta[idx] - state.bottom[idx];
        if (diff > 0.0 && diff < parameters_.min_depth) {
            state.eta[idx] = state.bottom[idx];
        }
    }
    for (index_t idx = 0; idx < grid_.nSurfaceCellMem(); ++idx) {
        state.dept[idx] = state.eta[idx] - state.bottom[idx];
        if (state.dept[idx] <= parameters_.min_depth) {
            state.dept[idx] = 0.0;
        }
        state.deptx[idx] = 0.0;
        state.depty[idx] = 0.0;
    }
    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        const auto ij = grid_.getSurfaceLogicalIndex(idx);
        const auto ip = static_cast<std::size_t>(xp(ij[0], ij[1]));
        const auto jp = static_cast<std::size_t>(yp(ij[0], ij[1]));
        state.deptx[idx] = std::max(state.eta[idx], state.eta[ip]) -
                           std::max(state.bottom[idx], state.bottom[ip]);
        state.depty[idx] = std::max(state.eta[idx], state.eta[jp]) -
                           std::max(state.bottom[idx], state.bottom[jp]);
    }
    for (int i = 0; i < grid_.nx(); ++i) {
        const auto south = grid_.getSurfaceIndex(i, 0);
        const auto north = grid_.getSurfaceIndex(i, grid_.ny() - 1);
        const auto ym_ghost = grid_.surfaceGhostIndex(Direction::YM, i, 0);
        const auto yp_ghost = grid_.surfaceGhostIndex(Direction::YP, i, grid_.ny() - 1);
        state.depty[ym_ghost] = std::max(state.eta[south], state.eta[ym_ghost]) -
                                std::max(state.bottom[south], state.bottom[ym_ghost]);
        state.depty[yp_ghost] = std::max(state.eta[north], state.eta[yp_ghost]) -
                                std::max(state.bottom[north], state.bottom[yp_ghost]);
    }
    for (int j = 0; j < grid_.ny(); ++j) {
        const auto west = grid_.getSurfaceIndex(0, j);
        const auto east = grid_.getSurfaceIndex(grid_.nx() - 1, j);
        const auto xm_ghost = grid_.surfaceGhostIndex(Direction::XM, 0, j);
        const auto xp_ghost = grid_.surfaceGhostIndex(Direction::XP, grid_.nx() - 1, j);
        state.deptx[xm_ghost] = std::max(state.eta[west], state.eta[xm_ghost]) -
                                std::max(state.bottom[west], state.bottom[xm_ghost]);
        state.deptx[xp_ghost] = std::max(state.eta[east], state.eta[xp_ghost]) -
                                std::max(state.bottom[east], state.bottom[xp_ghost]);
    }
    for (index_t idx = 0; idx < grid_.nSurfaceCellMem(); ++idx) {
        if (state.deptx[idx] < 0.0) {
            state.deptx[idx] = 0.0;
        }
        if (state.depty[idx] < 0.0) {
            state.depty[idx] = 0.0;
        }
    }
}

void SweSolver::updateLegacySubgridVariables(LegacySweState& state, bool zero_degenerate_faces) const
{
    const real area = grid_.dx() * grid_.dy();
    for (index_t idx = 0; idx < grid_.nSurfaceCellMem(); ++idx) {
        state.Vsn[idx] = state.Vs[idx];
        state.Vs[idx] = state.dept[idx] * area;
        state.Asz[idx] = state.dept[idx] > 0.0 ? area : 0.0;
        state.Asx[idx] = state.deptx[idx] * grid_.dy();
        state.Asy[idx] = state.depty[idx] * grid_.dx();
        if (zero_degenerate_faces && grid_.nx() == 1) {
            state.Asx[idx] = 0.0;
        }
        if (zero_degenerate_faces && grid_.ny() == 1) {
            state.Asy[idx] = 0.0;
        }
    }
    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        const auto ij = grid_.getSurfaceLogicalIndex(idx);
        const auto ip = static_cast<std::size_t>(xp(ij[0], ij[1]));
        const auto jp = static_cast<std::size_t>(yp(ij[0], ij[1]));
        state.Vsx[idx] = 0.5 * (state.Vs[idx] + state.Vs[ip]);
        state.Vsy[idx] = 0.5 * (state.Vs[idx] + state.Vs[jp]);
        state.Aszx[idx] = 0.5 * (state.Asz[idx] + state.Asz[ip]);
        state.Aszy[idx] = 0.5 * (state.Asz[idx] + state.Asz[jp]);
    }
    if (zero_degenerate_faces) {
        for (int i = 0; i < grid_.nx(); ++i) {
            const auto south = grid_.getSurfaceIndex(i, 0);
            const auto north = grid_.getSurfaceIndex(i, grid_.ny() - 1);
            const auto ym_ghost = grid_.surfaceGhostIndex(Direction::YM, i, 0);
            const auto yp_ghost = grid_.surfaceGhostIndex(Direction::YP, i, grid_.ny() - 1);
            state.Vs[ym_ghost] = state.Vs[south];
            state.Vs[yp_ghost] = state.Vs[north];
            state.Asx[ym_ghost] = state.Asx[south];
            state.Asx[yp_ghost] = state.Asx[north];
            state.Asy[ym_ghost] = state.Asy[south];
            state.Asy[yp_ghost] = state.Asy[north];
            state.Asz[ym_ghost] = state.Asz[south];
            state.Asz[yp_ghost] = state.Asz[north];
            state.Vsy[ym_ghost] = state.Vs[south];
            state.Aszy[ym_ghost] = state.Asz[south];
        }
        for (int j = 0; j < grid_.ny(); ++j) {
            const auto west = grid_.getSurfaceIndex(0, j);
            const auto east = grid_.getSurfaceIndex(grid_.nx() - 1, j);
            const auto xm_ghost = grid_.surfaceGhostIndex(Direction::XM, 0, j);
            const auto xp_ghost = grid_.surfaceGhostIndex(Direction::XP, grid_.nx() - 1, j);
            state.Vs[xm_ghost] = state.Vs[west];
            state.Vs[xp_ghost] = state.Vs[east];
            state.Asx[xm_ghost] = state.Asx[west];
            state.Asx[xp_ghost] = state.Asx[east];
            state.Asy[xm_ghost] = state.Asy[west];
            state.Asy[xp_ghost] = state.Asy[east];
            state.Asz[xm_ghost] = state.Asz[west];
            state.Asz[xp_ghost] = state.Asz[east];
            state.Vsx[xm_ghost] = state.Vs[west];
            state.Aszx[xm_ghost] = state.Asz[west];
        }
    } else {
        for (int i = 0; i < grid_.nx(); ++i) {
            const auto south = grid_.getSurfaceIndex(i, 0);
            const auto ym_ghost = grid_.surfaceGhostIndex(Direction::YM, i, 0);
            state.Vsy[ym_ghost] = state.Vs[south];
            state.Aszy[ym_ghost] = state.Asz[south];
        }
        for (int j = 0; j < grid_.ny(); ++j) {
            const auto west = grid_.getSurfaceIndex(0, j);
            const auto xm_ghost = grid_.surfaceGhostIndex(Direction::XM, 0, j);
            state.Vsx[xm_ghost] = state.Vs[west];
            state.Aszx[xm_ghost] = state.Asz[west];
        }
    }
}

void SweSolver::updateLegacyDragCoefficient(LegacySweState& state) const
{
    const real coef = parameters_.gravity * parameters_.manning * parameters_.manning;
    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        if (state.Vs[idx] > 0.0) {
            const real effh = state.Vs[idx] / (grid_.dx() * grid_.dy());
            const real expo = effh < parameters_.hD ? 2.0 / 3.0 : 1.0 / 3.0;
            state.CDx[idx] = coef / std::pow(effh, expo);
            state.CDy[idx] = coef / std::pow(effh, expo);
        }
    }
}

void SweSolver::legacyMomentumSource(LegacySweState& state) const
{
    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        const auto ij = grid_.getSurfaceLogicalIndex(idx);
        const auto im = static_cast<std::size_t>(xm(ij[0], ij[1]));
        const auto ip = static_cast<std::size_t>(xp(ij[0], ij[1]));
        const auto jm = static_cast<std::size_t>(ym(ij[0], ij[1]));
        const auto jp = static_cast<std::size_t>(yp(ij[0], ij[1]));

        real adv_x = (0.5 / grid_.dx()) *
                         ((state.uu[idx] + std::abs(state.uu[idx])) * (state.uu[idx] - state.uu[im]) +
                          (state.uu[idx] - std::abs(state.uu[idx])) * (state.uu[ip] - state.uu[idx])) +
                     (0.5 / grid_.dy()) *
                         ((state.vx[idx] + std::abs(state.vx[idx])) * (state.uu[idx] - state.uu[jm]) +
                          (state.vx[idx] - std::abs(state.vx[idx])) * (state.uu[jp] - state.uu[idx]));
        real adv_y = (0.5 / grid_.dx()) *
                         ((state.uy[idx] + std::abs(state.uy[idx])) * (state.vv[idx] - state.vv[im]) +
                          (state.uy[idx] - std::abs(state.uy[idx])) * (state.vv[ip] - state.vv[idx])) +
                     (0.5 / grid_.dy()) *
                         ((state.vv[idx] + std::abs(state.vv[idx])) * (state.vv[idx] - state.vv[jm]) +
                          (state.vv[idx] - std::abs(state.vv[idx])) * (state.vv[jp] - state.vv[idx]));
        if (state.uu[idx] == 0.0 || state.cflx[idx] > 0.7) {
            adv_x = 0.0;
        } else if (state.cflx[idx] > 0.5) {
            adv_x *= (0.7 - state.cflx[idx]) / (0.7 - 0.5);
        }
        if (state.vv[idx] == 0.0 || state.cfly[idx] > 0.7) {
            adv_y = 0.0;
        } else if (state.cfly[idx] > 0.5) {
            adv_y *= (0.7 - state.cfly[idx]) / (0.7 - 0.5);
        }

        real dif_x = 0.0;
        real dif_y = 0.0;
        if (state.Vsx[idx] > 0.0) {
            dif_x = (parameters_.viscx / state.Vsx[idx] / grid_.dx()) *
                        (state.Asx[idx] * (state.uu[ip] - state.uu[idx]) -
                         state.Asx[idx] * (state.uu[idx] - state.uu[im])) +
                    (parameters_.viscy / state.Vsx[idx] / grid_.dy()) *
                        (state.Asy[idx] * (state.uu[jp] - state.uu[idx]) -
                         state.Asy[jm] * (state.uu[idx] - state.uu[jm]));
        }
        if (state.Vsy[idx] > 0.0) {
            dif_y = (parameters_.viscx / state.Vsy[idx] / grid_.dx()) *
                        (state.Asx[idx] * (state.vv[ip] - state.vv[idx]) -
                         state.Asx[im] * (state.vv[idx] - state.vv[im])) +
                    (parameters_.viscy / state.Vsy[idx] / grid_.dy()) *
                        (state.Asy[idx] * (state.vv[jp] - state.vv[idx]) -
                         state.Asy[idx] * (state.vv[idx] - state.vv[jm]));
        }

        real facdx = 0.0;
        real facdy = 0.0;
        const real velx = std::sqrt(state.uu[idx] * state.uu[idx] + state.vx[idx] * state.vx[idx]);
        const real vely = std::sqrt(state.uy[idx] * state.uy[idx] + state.vv[idx] * state.vv[idx]);
        if (state.Vsx[idx] > 0.0) {
            facdx = state.Aszx[idx] / state.Vsx[idx];
        }
        if (state.Vsy[idx] > 0.0) {
            facdy = state.Aszy[idx] / state.Vsy[idx];
        }
        state.Dx[idx] = 1.0 / (0.5 * parameters_.dt * state.CDx[idx] * velx * facdx + 1.0);
        state.Dy[idx] = 1.0 / (0.5 * parameters_.dt * state.CDy[idx] * vely * facdy + 1.0);
        state.Ex[idx] = (state.uu[idx] + parameters_.dt * (dif_x - adv_x)) * state.Dx[idx];
        state.Ey[idx] = (state.vv[idx] + parameters_.dt * (dif_y - adv_y)) * state.Dy[idx];
    }
}

bool SweSolver::isFreeOutflowBoundaryCell(index_t cell, real normal_x, real normal_y) const
{
    for (const auto& boundary : parameters_.free_outflow_boundaries) {
        if (boundary.normal_x * normal_x + boundary.normal_y * normal_y <= 0.0) {
            continue;
        }
        if (std::find(boundary.cells.begin(), boundary.cells.end(), cell) != boundary.cells.end()) {
            return true;
        }
    }
    return false;
}

void SweSolver::enforceLegacyReflectiveBoundaryFluxes(LegacySweState& state) const
{
    for (int i = 0; i < grid_.nx(); ++i) {
        const auto south = grid_.getSurfaceIndex(i, 0);
        const auto north = grid_.getSurfaceIndex(i, grid_.ny() - 1);
        const auto ym_ghost = grid_.surfaceGhostIndex(Direction::YM, i, 0);

        if (!isFreeOutflowBoundaryCell(south, 0.0, -1.0)) {
            state.vv[ym_ghost] = 0.0;
            state.Ey[ym_ghost] = 0.0;
            state.Fv[ym_ghost] = 0.0;
        }
        if (!isFreeOutflowBoundaryCell(north, 0.0, 1.0)) {
            state.vv[north] = 0.0;
            state.Ey[north] = 0.0;
            state.Fv[north] = 0.0;
        }
    }

    for (int j = 0; j < grid_.ny(); ++j) {
        const auto west = grid_.getSurfaceIndex(0, j);
        const auto east = grid_.getSurfaceIndex(grid_.nx() - 1, j);
        const auto xm_ghost = grid_.surfaceGhostIndex(Direction::XM, 0, j);

        if (!isFreeOutflowBoundaryCell(west, -1.0, 0.0)) {
            state.uu[xm_ghost] = 0.0;
            state.Ex[xm_ghost] = 0.0;
            state.Fu[xm_ghost] = 0.0;
        }
        if (!isFreeOutflowBoundaryCell(east, 1.0, 0.0)) {
            state.uu[east] = 0.0;
            state.Ex[east] = 0.0;
            state.Fu[east] = 0.0;
        }
    }
}

SweLinearSystem SweSolver::assembleLegacyLinearSystem(LegacySweState& state) const
{
    const auto n2ci = static_cast<std::size_t>(grid_.nSurfaceCell());
    SweLinearSystem system;
    system.n = static_cast<int>(n2ci);
    system.dom2mat.resize(n2ci);
    system.mat2dom.resize(n2ci);
    system.rhs.assign(n2ci, 0.0);
    for (std::size_t idx = 0; idx < n2ci; ++idx) {
        system.dom2mat[idx] = static_cast<int>(idx);
        system.mat2dom[idx] = static_cast<int>(idx);
    }

    const real coef = parameters_.gravity * parameters_.dt * parameters_.dt;
    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        const auto ij = grid_.getSurfaceLogicalIndex(idx);
        const auto im = static_cast<std::size_t>(xm(ij[0], ij[1]));
        const auto jm = static_cast<std::size_t>(ym(ij[0], ij[1]));
        const auto ip = static_cast<std::size_t>(xp(ij[0], ij[1]));
        const auto jp = static_cast<std::size_t>(yp(ij[0], ij[1]));

        state.Srhs[idx] = state.eta[idx] * state.Asz[idx] -
                          parameters_.dt *
                              (state.Asx[idx] * state.Ex[idx] - state.Asx[im] * state.Ex[im] +
                               state.Asy[idx] * state.Ey[idx] - state.Asy[jm] * state.Ey[jm]);
        state.Sxp[idx] = state.Vsx[idx] > 0.0
                             ? coef * state.Asx[idx] * state.Asx[idx] * state.Dx[idx] / state.Vsx[idx]
                             : 0.0;
        state.Sxm[idx] = state.Vsx[im] > 0.0
                             ? coef * state.Asx[im] * state.Asx[im] * state.Dx[im] / state.Vsx[im]
                             : 0.0;
        state.Syp[idx] = state.Vsy[idx] > 0.0
                             ? coef * state.Asy[idx] * state.Asy[idx] * state.Dy[idx] / state.Vsy[idx]
                             : 0.0;
        state.Sym[idx] = state.Vsy[jm] > 0.0
                             ? coef * state.Asy[jm] * state.Asy[jm] * state.Dy[jm] / state.Vsy[jm]
                             : 0.0;
        state.Sct[idx] = state.Asz[idx] + state.Sxp[idx] + state.Sxm[idx] + state.Syp[idx] + state.Sym[idx];
        if (state.dept[idx] == 0.0) {
            state.Sct[idx] = grid_.dx() * grid_.dy();
            state.Srhs[idx] = state.eta[idx] * grid_.dx() * grid_.dy();
            if (state.uu[idx] == 0.0 && state.uu[im] == 0.0 && state.vv[idx] == 0.0 &&
                state.vv[jm] == 0.0) {
                state.Sxp[idx] = 0.0;
                state.Sxm[idx] = 0.0;
                state.Syp[idx] = 0.0;
                state.Sym[idx] = 0.0;
            }
        } else if (state.Sct[idx] == 0.0) {
            state.Sct[idx] = grid_.dx() * grid_.dy();
            state.Srhs[idx] = state.eta[idx] * grid_.dx() * grid_.dy();
        }

        if (ij[0] == 0) {
            state.Sct[idx] -= state.Sxm[idx];
            if (grid_.nx() == 1) {
                state.Sct[idx] -= state.Sxp[idx];
            }
        } else if (ij[0] == grid_.nx() - 1) {
            state.Sct[idx] -= state.Sxp[idx];
        }
        if (ij[1] == 0) {
            state.Sct[idx] -= state.Sym[idx];
            if (grid_.ny() == 1) {
                state.Sct[idx] -= state.Syp[idx];
            }
        } else if (ij[1] == grid_.ny() - 1) {
            state.Sct[idx] -= state.Syp[idx];
        }

        auto add = [&system, idx](int col_domain, real value) {
            if (col_domain >= 0 &&
                static_cast<std::size_t>(col_domain) < system.dom2mat.size()) {
                system.entries.push_back(
                    SweMatrixEntry{static_cast<int>(idx), system.dom2mat[static_cast<std::size_t>(col_domain)], value});
            }
        };
        if (ij[1] > 0) {
            add(static_cast<int>(jm), -state.Sym[idx]);
        }
        if (ij[0] > 0) {
            add(static_cast<int>(im), -state.Sxm[idx]);
        }
        add(static_cast<int>(idx), state.Sct[idx]);
        if (ij[0] < grid_.nx() - 1) {
            add(static_cast<int>(ip), -state.Sxp[idx]);
        }
        if (ij[1] < grid_.ny() - 1) {
            add(static_cast<int>(jp), -state.Syp[idx]);
        }
        system.rhs[idx] = state.Srhs[idx];
    }
    return system;
}

void SweSolver::legacyCflLimiter(LegacySweState& state, real rain_rate) const
{
    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        const auto ij = grid_.getSurfaceLogicalIndex(idx);
        const real diff = state.eta[idx] - state.bottom[idx];
        if (state.dept[idx] <= 0.0 && diff > 0.0 && rain_rate <= 0.0) {
            bool wet = false;
            wet = wet || state.dept[static_cast<std::size_t>(xp(ij[0], ij[1]))] > 0.0;
            wet = wet || state.dept[static_cast<std::size_t>(xm(ij[0], ij[1]))] > 0.0;
            wet = wet || state.dept[static_cast<std::size_t>(yp(ij[0], ij[1]))] > 0.0;
            wet = wet || state.dept[static_cast<std::size_t>(ym(ij[0], ij[1]))] > 0.0;
            if (!wet) {
                state.eta[idx] = state.bottom[idx];
                state.cfl_active[idx] = 1.0;
            }
        }
    }
    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        const real diff = state.eta[idx] - state.bottom[idx];
        if (diff > 0.0 && diff < parameters_.min_depth) {
            state.eta[idx] = state.bottom[idx];
        }
    }
}

void SweSolver::legacyEvapRain(LegacySweState& state, real rain_rate, real evap_rate) const
{
    state.rain_sum += rain_rate * parameters_.dt;
    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        const auto ij = grid_.getSurfaceLogicalIndex(idx);
        if (ij[1] != grid_.ny() - 1) {
            state.eta[idx] += rain_rate * parameters_.dt;
        }
        state.eta[idx] -= evap_rate * parameters_.dt;
    }
    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        const real diff = state.eta[idx] - state.bottom[idx];
        if (diff < 0.0 || diff < parameters_.min_depth) {
            state.eta[idx] = state.bottom[idx];
        }
    }
}

void SweSolver::legacyWaterfallLocation(LegacySweState& state) const
{
    std::fill(state.wtfx.begin(), state.wtfx.end(), 0.0);
    std::fill(state.wtfy.begin(), state.wtfy.end(), 0.0);
    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        const auto ij = grid_.getSurfaceLogicalIndex(idx);
        const auto ip = static_cast<std::size_t>(xp(ij[0], ij[1]));
        const auto im = static_cast<std::size_t>(xm(ij[0], ij[1]));
        const auto jp = static_cast<std::size_t>(yp(ij[0], ij[1]));
        const auto jm = static_cast<std::size_t>(ym(ij[0], ij[1]));
        real facd = state.eta[ip] - state.bottom_xp[idx];
        if (state.eta[idx] < state.bottom_xp[idx] && facd > parameters_.wtfh) {
            state.wtfx[idx] = -1.0;
        }
        facd = state.eta[jp] - state.bottom_yp[idx];
        if (state.eta[idx] < state.bottom_yp[idx] && facd > parameters_.wtfh) {
            state.wtfy[idx] = -1.0;
        }
        facd = state.eta[im] - state.bottom_xp[im];
        if (state.eta[idx] < state.bottom_xp[im] && facd > parameters_.wtfh) {
            state.wtfx[idx] = 1.0;
        }
        facd = state.eta[jm] - state.bottom_yp[jm];
        if (state.eta[idx] < state.bottom_yp[jm] && facd > parameters_.wtfh) {
            state.wtfy[idx] = 1.0;
        }
    }
}

void SweSolver::updateLegacyVelocity(LegacySweState& state) const
{
    const real coef = parameters_.gravity * parameters_.dt;
    for (index_t idx = 0; idx < grid_.nSurfaceCellMem(); ++idx) {
        state.un[idx] = state.uu[idx];
        state.vn[idx] = state.vv[idx];
    }
    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        const auto ij = grid_.getSurfaceLogicalIndex(idx);
        const auto ip = static_cast<std::size_t>(xp(ij[0], ij[1]));
        const auto jp = static_cast<std::size_t>(yp(ij[0], ij[1]));
        const real effhx = state.Vsx[idx] > 0.0 ? state.Asx[idx] / state.Vsx[idx] : 0.0;
        const real effhy = state.Vsy[idx] > 0.0 ? state.Asy[idx] / state.Vsy[idx] : 0.0;
        state.uu[idx] = (state.Ex[idx] - coef * effhx * (state.eta[ip] - state.eta[idx])) * state.Dx[idx];
        state.vv[idx] = (state.Ey[idx] - coef * effhy * (state.eta[jp] - state.eta[idx])) * state.Dy[idx];
    }
    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        const auto ij = grid_.getSurfaceLogicalIndex(idx);
        const auto im = static_cast<std::size_t>(xm(ij[0], ij[1]));
        const auto jm = static_cast<std::size_t>(ym(ij[0], ij[1]));
        if (state.Asx[idx] < parameters_.wtfh * grid_.dy()) {
            state.uu[idx] = 0.0;
        }
        if (state.Asy[idx] < parameters_.wtfh * grid_.dx()) {
            state.vv[idx] = 0.0;
        }
        if (state.dept[idx] < parameters_.wtfh) {
            if (state.uu[idx] > 0.0) {
                state.uu[idx] = 0.0;
            }
            if (state.uu[im] < 0.0) {
                state.uu[im] = 0.0;
            }
            if (state.vv[idx] > 0.0) {
                state.vv[idx] = 0.0;
            }
            if (state.vv[jm] < 0.0) {
                state.vv[jm] = 0.0;
            }
        }
        if (state.cfl_active[idx] == 1.0) {
            state.uu[idx] = 0.0;
            state.uu[im] = 0.0;
            state.vv[idx] = 0.0;
            state.vv[jm] = 0.0;
            state.cfl_active[idx] = 0.0;
        }
    }
    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        state.Fu[idx] = state.uu[idx] * state.Asx[idx];
        state.Fv[idx] = state.vv[idx] * state.Asy[idx];
        state.cflx[idx] = std::abs(state.uu[idx] * parameters_.dt / grid_.dx());
        state.cfly[idx] = std::abs(state.vv[idx] * parameters_.dt / grid_.dy());
    }
}

void SweSolver::enforceLegacyVelocityBc(LegacySweState& state) const
{
    for (int i = 0; i < grid_.nx(); ++i) {
        const auto south = grid_.getSurfaceIndex(i, 0);
        const auto north = grid_.getSurfaceIndex(i, grid_.ny() - 1);
        const auto ym_ghost = grid_.surfaceGhostIndex(Direction::YM, i, 0);
        const auto yp_ghost = grid_.surfaceGhostIndex(Direction::YP, i, grid_.ny() - 1);
        state.uu[ym_ghost] = state.uu[south];
        state.vv[ym_ghost] = state.vv[south];
        state.uu[yp_ghost] = state.uu[north];
        state.vv[yp_ghost] = state.vv[north];
    }
    for (int j = 0; j < grid_.ny(); ++j) {
        const auto west = grid_.getSurfaceIndex(0, j);
        const auto east = grid_.getSurfaceIndex(grid_.nx() - 1, j);
        const auto xm_ghost = grid_.surfaceGhostIndex(Direction::XM, 0, j);
        const auto xp_ghost = grid_.surfaceGhostIndex(Direction::XP, grid_.nx() - 1, j);
        state.uu[xm_ghost] = state.uu[west];
        state.vv[xm_ghost] = state.vv[west];
        state.uu[xp_ghost] = state.uu[east];
        state.vv[xp_ghost] = state.vv[east];
    }
}

void SweSolver::applyLegacyFreeOutflowBc(LegacySweState& state) const
{
    for (const auto& boundary : parameters_.free_outflow_boundaries) {
        for (const auto cell : boundary.cells) {
            if (cell >= state.uu.size()) {
                continue;
            }
            const real speed =
                std::sqrt(state.uu[cell] * state.uu[cell] + state.vv[cell] * state.vv[cell]);
            state.uu[cell] = boundary.normal_x * speed;
            state.vv[cell] = boundary.normal_y * speed;
        }
    }
}

void SweSolver::interpolateLegacyVelocity(LegacySweState& state) const
{
    for (index_t idx = 0; idx < grid_.nSurfaceCell(); ++idx) {
        const auto ij = grid_.getSurfaceLogicalIndex(idx);
        const auto im = static_cast<std::size_t>(xm(ij[0], ij[1]));
        const auto jm = static_cast<std::size_t>(ym(ij[0], ij[1]));
        const auto ip = static_cast<std::size_t>(xp(ij[0], ij[1]));
        const auto jp = static_cast<std::size_t>(yp(ij[0], ij[1]));
        const auto jp_minus_one = static_cast<std::size_t>(static_cast<int>(jp) - 1);
        const auto ip_minus_nx = static_cast<std::size_t>(static_cast<int>(ip) - grid_.nx());
        state.uy[idx] = 0.25 * (state.uu[idx] + state.uu[im] + state.uu[jp] +
                                state.uu[jp_minus_one]);
        state.vx[idx] = 0.25 * (state.vv[idx] + state.vv[jm] + state.vv[ip] +
                                state.vv[ip_minus_nx]);
    }
}

int SweSolver::xp(int i, int j) const
{
    return i == grid_.nx() - 1 ? static_cast<int>(grid_.surfaceGhostIndex(Direction::XP, i, j))
                               : static_cast<int>(grid_.getSurfaceIndex(i + 1, j));
}

int SweSolver::xm(int i, int j) const
{
    return i == 0 ? static_cast<int>(grid_.surfaceGhostIndex(Direction::XM, i, j))
                  : static_cast<int>(grid_.getSurfaceIndex(i - 1, j));
}

int SweSolver::yp(int i, int j) const
{
    return j == grid_.ny() - 1 ? static_cast<int>(grid_.surfaceGhostIndex(Direction::YP, i, j))
                               : static_cast<int>(grid_.getSurfaceIndex(i, j + 1));
}

int SweSolver::ym(int i, int j) const
{
    return j == 0 ? static_cast<int>(grid_.surfaceGhostIndex(Direction::YM, i, j))
                  : static_cast<int>(grid_.getSurfaceIndex(i, j - 1));
}

}  // namespace frehg2
