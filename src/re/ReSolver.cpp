#include "re/ReSolver.hpp"

#include "core/InitialCondition.hpp"
#include "io/AsciiRaster.hpp"

#ifdef USE_PETSC
#include <petsc.h>
#endif

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace frehg2 {

namespace {

constexpr real EPS = 1.0e-7;
constexpr real WC_HEAD_EPS = 1.0e-5;

real mvgWaterContentMax(const ReParameters& parameters)
{
    if (!parameters.use_mvg) {
        return parameters.wcs;
    }
    const real m = 1.0 - 1.0 / parameters.soil_n;
    return parameters.wcr + (parameters.wcs - parameters.wcr) *
                                std::pow(1.0 + std::pow(std::abs(parameters.aev) *
                                                             parameters.soil_a,
                                                         parameters.soil_n),
                                         m);
}

real clamp(real value, real lower, real upper)
{
    return std::max(lower, std::min(value, upper));
}

}  // namespace

ReSolver::ReSolver(Grid grid, ReParameters parameters)
    : grid_(std::move(grid)),
      parameters_(parameters)
{
    if (!parameters_.qtop_surface.empty() &&
        parameters_.qtop_surface.size() != static_cast<std::size_t>(grid_.nSurfaceCell())) {
        throw std::invalid_argument("groundwater qtop_surface size must match surface grid cells");
    }
}

const Grid& ReSolver::grid() const noexcept
{
    return grid_;
}

const ReParameters& ReSolver::parameters() const noexcept
{
    return parameters_;
}

SoilParameters ReSolver::uniformSoil(const ReParameters& parameters)
{
    SoilParameters soil;
    soil.ksx = parameters.ksx;
    soil.ksy = parameters.ksy;
    soil.ksz = parameters.ksz;
    soil.soil_a = parameters.soil_a;
    soil.soil_n = parameters.soil_n;
    soil.wcs = parameters.wcs;
    soil.wcr = parameters.wcr;
    soil.aev = parameters.aev;
    return soil;
}

ReParameters ReSolver::parametersForSoil(
    const ReParameters& parameters,
    const SoilParameters& soil)
{
    ReParameters result = parameters;
    result.ksx = soil.ksx;
    result.ksy = soil.ksy;
    result.ksz = soil.ksz;
    result.soil_a = soil.soil_a;
    result.soil_n = soil.soil_n;
    result.wcs = soil.wcs;
    result.wcr = soil.wcr;
    result.aev = soil.aev;
    return result;
}

std::vector<int> ReSolver::expandSurfaceSoilMap(
    const Grid& grid,
    const std::vector<int>& surface_soil_id)
{
    if (surface_soil_id.size() != static_cast<std::size_t>(grid.nSurfaceCell())) {
        throw std::invalid_argument("surface soil map size must match surface grid cells");
    }
    std::vector<int> soil_id(static_cast<std::size_t>(grid.nCell()), 0);
    for (int j = 0; j < grid.ny(); ++j) {
        for (int i = 0; i < grid.nx(); ++i) {
            const auto surface = grid.getSurfaceIndex(i, j);
            for (int k = 0; k < grid.nz(); ++k) {
                soil_id[grid.getIndex(i, j, k)] = surface_soil_id[surface];
            }
        }
    }
    return soil_id;
}

std::vector<int> ReSolver::readSoilMap(const Grid& grid, const std::string& filename)
{
    const auto values = InitialCondition::readAsciiGroundwaterField(grid, filename);
    std::vector<int> soil_id(values.size(), 0);
    for (std::size_t cell = 0; cell < values.size(); ++cell) {
        soil_id[cell] = static_cast<int>(std::lround(values[cell]));
    }
    return soil_id;
}

std::vector<int> ReSolver::readSurfaceSoilMap(const Grid& grid, const std::string& filename)
{
    AsciiRaster raster;
    raster.read(filename);
    if (raster.ncols() != grid.nx() || raster.nrows() != grid.ny()) {
        throw std::runtime_error("soil map dimensions must match the surface grid");
    }

    std::vector<int> surface_soil_id(static_cast<std::size_t>(grid.nSurfaceCell()), 0);
    for (int j = 0; j < grid.ny(); ++j) {
        for (int i = 0; i < grid.nx(); ++i) {
            const auto idx = grid.getSurfaceIndex(i, j);
            if (raster.active(idx) == 0) {
                throw std::runtime_error("soil map contains NODATA in an active groundwater cell");
            }
            surface_soil_id[idx] = static_cast<int>(raster.value(idx));
        }
    }
    return surface_soil_id;
}

int ReSolver::legacyNz(real bath, real bot_z, real dz, real dz_multiplier)
{
    if (dz <= 0.0 || dz_multiplier <= 0.0) {
        throw std::invalid_argument("groundwater dz and dz multiplier must be positive");
    }
    if (dz_multiplier == 1.0) {
        return static_cast<int>(std::ceil((bath - bot_z) / dz));
    }

    real bot_new = bath - dz;
    real dz_new = dz;
    int nz = 1;
    while (bot_new > bot_z) {
        dz_new *= dz_multiplier;
        bot_new -= dz_new;
        ++nz;
    }
    return nz;
}

real ReSolver::waterContentFromHead(real h, const ReParameters& parameters)
{
    real wc = parameters.wcr;
    if (parameters.use_vg) {
        const real m = 1.0 - 1.0 / parameters.soil_n;
        const real wcm = mvgWaterContentMax(parameters);
        const real s = std::pow(
            1.0 + std::pow(std::abs(parameters.soil_a * h), parameters.soil_n), -m);
        wc = h > parameters.aev ? parameters.wcs : parameters.wcr + (wcm - parameters.wcr) * s;
    } else {
        wc = parameters.wcr + (parameters.wcs - parameters.wcr) * std::exp(0.1634 * h);
    }
    return clamp(wc, parameters.wcr, parameters.wcs);
}

real ReSolver::headFromWaterContent(real wc, const ReParameters& parameters)
{
    if (parameters.use_vg) {
        const real m = 1.0 - 1.0 / parameters.soil_n;
        const real wcm = mvgWaterContentMax(parameters);
        wc = std::max(wc, parameters.wcr + WC_HEAD_EPS);
        if (wc < parameters.wcs) {
            return -(1.0 / parameters.soil_a) *
                   std::pow(std::pow((wcm - parameters.wcr) / (wc - parameters.wcr), 1.0 / m) -
                                1.0,
                            1.0 / parameters.soil_n);
        }
        return 0.0;
    }
    return std::log((wc - parameters.wcr) / (parameters.wcs - parameters.wcr)) / 0.1634;
}

real ReSolver::specificCapacity(real h, const ReParameters& parameters)
{
    const real m = 1.0 - 1.0 / parameters.soil_n;
    const real wcm = mvgWaterContentMax(parameters);
    const real ah = std::abs(parameters.soil_a * h);
    const real numerator = parameters.soil_a * parameters.soil_n * m *
                           (wcm - parameters.wcr) *
                           std::pow(ah, parameters.soil_n - 1.0);
    const real denominator = std::pow(1.0 + std::pow(ah, parameters.soil_n), m + 1.0);
    real c = numerator / denominator;
    if (parameters.use_mvg) {
        if (h > parameters.aev) {
            c = 0.0;
        }
    } else if (h > 0.0) {
        c = 0.0;
    }
    return c;
}

real ReSolver::hydraulicConductivity(real h, real ks, const ReParameters& parameters)
{
    real effective_k = ks;
    if (parameters.use_vg) {
        const real m = 1.0 - 1.0 / parameters.soil_n;
        const real s = std::pow(
            1.0 + std::pow(std::abs(parameters.soil_a * h), parameters.soil_n), -m);
        if (parameters.use_mvg) {
            const real wcm = mvgWaterContentMax(parameters);
            const real numerator =
                1.0 - std::pow(1.0 - std::pow(s * (parameters.wcs - parameters.wcr) /
                                                   (wcm - parameters.wcr),
                                               1.0 / m),
                                m);
            const real denominator =
                1.0 - std::pow(1.0 - std::pow((parameters.wcs - parameters.wcr) /
                                                   (wcm - parameters.wcr),
                                               1.0 / m),
                                m);
            effective_k = denominator == 0.0
                              ? ks
                              : ks * std::pow(s, 0.5) * std::pow(numerator / denominator, 2.0);
        } else {
            effective_k =
                ks * std::pow(s, 0.5) *
                std::pow(1.0 - std::pow(1.0 - std::pow(s, 1.0 / m), m), 2.0);
        }
    } else {
        effective_k = ks * std::exp(0.1634 * h);
    }

    if (effective_k > ks) {
        effective_k = ks;
    }
    if (parameters.use_mvg) {
        if (h > parameters.aev) {
            effective_k = ks;
        }
    } else if (h > 0.0) {
        effective_k = ks;
    }
    return effective_k;
}

real ReSolver::conductivityFace(real h_left, real h_right, real ks, const ReParameters& parameters)
{
    return 0.5 * (hydraulicConductivity(h_left, ks, parameters) +
                  hydraulicConductivity(h_right, ks, parameters));
}

real ReSolver::conductivityDerivativeWc(real wc, real ks, const ReParameters& parameters)
{
    const real m = 1.0 - 1.0 / parameters.soil_n;
    const real wcm = mvgWaterContentMax(parameters);
    if (wc > 0.9999 * parameters.wcs && wc < parameters.wcs) {
        wc = 0.9999 * parameters.wcs;
    }
    const real s = (wc - parameters.wcr) / (parameters.wcs - parameters.wcr);
    if (s <= 0.0) {
        return 0.0;
    }

    real term1 = 0.0;
    real term2 = 0.0;
    if (!parameters.use_mvg) {
        const real term0 = std::pow(1.0 - std::pow(s, 1.0 / m), m);
        term1 = 0.5 * ks * std::pow(s, -0.5) * (1.0 - term0) * (1.0 - term0);
        term2 = 2.0 * ks * std::pow(s, (2.0 - m) / (2.0 * m)) * (1.0 - term0) *
                std::pow(1.0 - std::pow(s, 1.0 / m), m - 1.0);
    } else {
        const real c2 = (parameters.wcs - parameters.wcr) / (wcm - parameters.wcr);
        const real c1 =
            1.0 / std::pow(1.0 - std::pow(1.0 - std::pow(c2, 1.0 / m), m), 2.0);
        const real term0 = std::pow(1.0 - std::pow(c2 * s, 1.0 / m), m);
        term1 = 0.5 * c1 * ks * std::pow(s, -0.5) * (1.0 - term0) * (1.0 - term0);
        term2 = 2.0 * c1 * ks * c2 * std::pow(s, 0.5) *
                std::pow(c2 * s, 1.0 / m - 1.0) * (1.0 - term0) *
                std::pow(1.0 - std::pow(c2 * s, 1.0 / m), m - 1.0);
    }
    return (term1 + term2) / (parameters.wcs - parameters.wcr);
}

LegacyReState ReSolver::initializeLegacyState(real bath, real bot_z) const
{
    (void)bot_z;
    const auto n = static_cast<std::size_t>(grid_.nCell());
    LegacyReState state;
    state.h.assign(n, 0.0);
    state.hn.assign(n, 0.0);
    state.hnm.assign(n, 0.0);
    state.wc.assign(n, parameters_.init_wc);
    state.wcn.assign(n, parameters_.init_wc);
    state.wch.assign(n, parameters_.init_wc);
    state.hwc.assign(n, 0.0);
    state.ch.assign(n, 0.0);
    state.kx.assign(n, 0.0);
    state.ky.assign(n, 0.0);
    state.kz.assign(n, 0.0);
    state.qx.assign(n, 0.0);
    state.qy.assign(n, 0.0);
    state.qz.assign(n, 0.0);
    state.qz_top = 0.0;
    state.qz_top_surface.assign(static_cast<std::size_t>(grid_.nSurfaceCell()), 0.0);
    state.gxp.assign(n, 0.0);
    state.gxm.assign(n, 0.0);
    state.gyp.assign(n, 0.0);
    state.gym.assign(n, 0.0);
    state.gzp.assign(n, 0.0);
    state.gzm.assign(n, 0.0);
    state.gct.assign(n, 0.0);
    state.grhs.assign(n, 0.0);
    state.bot3d.assign(n, 0.0);
    state.dz3d.assign(n, 0.0);
    state.zcntr.assign(n, 0.0);
    state.active.assign(n, 1);
    state.istop.assign(n, 0);
    state.soil_id.assign(n, 0);
    state.dtg = parameters_.dt;

    if (!parameters_.soil_id.empty()) {
        if (parameters_.soil_id.size() != n) {
            throw std::invalid_argument("groundwater soil_id size must match grid cells");
        }
        state.soil_id = parameters_.soil_id;
    }

    for (int j = 0; j < grid_.ny(); ++j) {
        for (int i = 0; i < grid_.nx(); ++i) {
            for (int k = 0; k < grid_.nz(); ++k) {
                const int cell = index(i, j, k);
                state.bot3d[cell] = bath - static_cast<real>(k + 1) * grid_.dz(k);
                state.dz3d[cell] = grid_.dz(k);
                state.zcntr[cell] = state.bot3d[cell] + 0.5 * state.dz3d[cell];
                state.istop[cell] = k == 0 ? 1 : 0;
            }
        }
    }

    for (std::size_t cell_idx = 0; cell_idx < n; ++cell_idx) {
        const auto cell = static_cast<int>(cell_idx);
        const auto cell_params = parametersForCell(state, cell);
        if (parameters_.init_wc >= cell_params.wcr && parameters_.init_wc <= cell_params.wcs) {
            state.wc[cell] = parameters_.init_wc;
            state.h[cell] = parameters_.init_wc == cell_params.wcs
                                ? bath - state.bot3d[cell] - 0.5 * state.dz3d[cell]
                                : headFromWaterContent(state.wc[cell], cell_params);
        } else if (parameters_.init_h <= 0.0) {
            state.h[cell] = parameters_.init_h;
            state.wc[cell] = waterContentFromHead(parameters_.init_h, cell_params);
        } else {
            const real water_table =
                parameters_.init_wt_rel > 0.0 ? bath - parameters_.init_wt_rel : parameters_.init_wt_abs;
            state.h[cell] = water_table - state.bot3d[cell] - 0.5 * state.dz3d[cell];
            state.wc[cell] = waterContentFromHead(state.h[cell], cell_params);
        }
    }

    for (std::size_t cell = 0; cell < n; ++cell) {
        const auto cell_params = parametersForCell(state, static_cast<int>(cell));
        state.hn[cell] = state.h[cell];
        state.hnm[cell] = state.h[cell];
        state.wcn[cell] = state.wc[cell];
        state.hwc[cell] = headFromWaterContent(state.wc[cell], cell_params);
        state.wch[cell] = waterContentFromHead(state.h[cell], cell_params);
        state.ch[cell] = specificCapacity(state.h[cell], cell_params);
        state.kx[cell] = hydraulicConductivity(state.h[cell], cell_params.ksx, cell_params);
        state.ky[cell] = hydraulicConductivity(state.h[cell], cell_params.ksy, cell_params);
        state.kz[cell] = hydraulicConductivity(state.h[cell], cell_params.ksz, cell_params);
    }
    return state;
}

void ReSolver::advanceLegacyStep(LegacyReState& state) const
{
    for (std::size_t cell = 0; cell < state.h.size(); ++cell) {
        const auto cell_params = parametersForCell(state, static_cast<int>(cell));
        state.hnm[cell] = state.hn[cell];
        state.hn[cell] = state.h[cell];
        state.wcn[cell] = state.wc[cell];
        state.hwc[cell] = headFromWaterContent(state.wc[cell], cell_params);
        state.wch[cell] = waterContentFromHead(state.h[cell], cell_params);
        state.ch[cell] = specificCapacity(state.h[cell], cell_params);
    }

    computeConductivityFaces(state);
    const auto system = assemblePredictorSystem(state);
    state.h = solveLinearSystem(system);
    computeConductivityFaces(state);
    computeDarcyFlux(state);
    if (parameters_.use_corrector) {
        updateWaterContent(state);
        enforceLegacyMoistureConsistency(state);
    } else {
        for (std::size_t cell = 0; cell < state.h.size(); ++cell) {
            state.wc[cell] = waterContentFromHead(
                state.h[cell],
                parametersForCell(state, static_cast<int>(cell)));
        }
    }

    for (std::size_t cell = 0; cell < state.h.size(); ++cell) {
        const auto cell_params = parametersForCell(state, static_cast<int>(cell));
        state.wc[cell] = clamp(state.wc[cell], cell_params.wcr, cell_params.wcs);
        state.hwc[cell] = headFromWaterContent(state.wc[cell], cell_params);
        state.wch[cell] = waterContentFromHead(state.h[cell], cell_params);
    }
    if (parameters_.dt_adjust) {
        state.dtg = adaptTimeStep(state);
    }
}

void ReSolver::computeConductivityFaces(LegacyReState& state) const
{
    for (int j = 0; j < grid_.ny(); ++j) {
        for (int i = 0; i < grid_.nx(); ++i) {
            for (int k = 0; k < grid_.nz(); ++k) {
                const int cell = index(i, j, k);
                const auto cell_params = parametersForCell(state, cell);
                if (i < grid_.nx() - 1) {
                    const int east = index(i + 1, j, k);
                    const auto east_params = parametersForCell(state, east);
                    const real k_here =
                        hydraulicConductivity(state.h[cell], cell_params.ksx, cell_params);
                    const real k_east =
                        hydraulicConductivity(state.h[east], east_params.ksx, east_params);
                    state.kx[cell] = 0.5 * (k_here + k_east);
                } else {
                    state.kx[cell] = 0.0;
                }
                if (j < grid_.ny() - 1) {
                    const int north = index(i, j + 1, k);
                    const auto north_params = parametersForCell(state, north);
                    const real k_here =
                        hydraulicConductivity(state.h[cell], cell_params.ksy, cell_params);
                    const real k_north =
                        hydraulicConductivity(state.h[north], north_params.ksy, north_params);
                    state.ky[cell] = 0.5 * (k_here + k_north);
                } else {
                    state.ky[cell] = 0.0;
                }
                if (k == grid_.nz() - 1) {
                    state.kz[cell] =
                        parameters_.bc_type[4] == 0
                            ? 0.0
                            : hydraulicConductivity(state.h[cell], cell_params.ksz, cell_params);
                } else {
                    const int below = zp(cell);
                    const auto below_params = parametersForCell(state, below);
                    const real k_here =
                        hydraulicConductivity(state.h[cell], cell_params.ksz, cell_params);
                    const real k_below =
                        hydraulicConductivity(state.h[below], below_params.ksz, below_params);
                    state.kz[cell] = 0.5 * (k_here + k_below);
                }
            }
        }
    }
}

ReLinearSystem ReSolver::assemblePredictorSystem(LegacyReState& state) const
{
    ReLinearSystem system;
    system.n = static_cast<int>(state.h.size());
    system.rhs.assign(state.h.size(), 0.0);

    for (int j = 0; j < grid_.ny(); ++j) {
        for (int i = 0; i < grid_.nx(); ++i) {
            const int surface_cell = static_cast<int>(grid_.getSurfaceIndex(i, j));
            for (int k = 0; k < grid_.nz(); ++k) {
                const int cell = index(i, j, k);
                const auto cell_params = parametersForCell(state, cell);
                const real volume = grid_.dx() * grid_.dy() * state.dz3d[cell];
                const real storage =
                    (state.ch[cell] + parameters_.ss * state.wcn[cell] / cell_params.wcs) * volume;
                real gxp = 0.0;
                real gxm = 0.0;
                real gyp = 0.0;
                real gym = 0.0;
                if (i < grid_.nx() - 1) {
                    gxp = -state.kx[cell] * state.dtg * volume / (grid_.dx() * grid_.dx());
                }
                if (i > 0) {
                    gxm = -state.kx[index(i - 1, j, k)] * state.dtg * volume /
                          (grid_.dx() * grid_.dx());
                }
                if (j < grid_.ny() - 1) {
                    gyp = -state.ky[cell] * state.dtg * volume / (grid_.dy() * grid_.dy());
                }
                if (j > 0) {
                    gym = -state.ky[index(i, j - 1, k)] * state.dtg * volume /
                          (grid_.dy() * grid_.dy());
                }
                real gzp = 0.0;
                if (k < grid_.nz() - 1) {
                    const real dzf = 0.5 * (state.dz3d[cell] + state.dz3d[zp(cell)]);
                    gzp = -state.kz[cell] * state.dtg / (state.dz3d[cell] * dzf) * volume;
                }

                real k_top = 0.0;
                if (k == 0) {
                    if (parameters_.bc_type[5] == 1) {
                        k_top = cell_params.ksz;
                    } else if (parameters_.bc_type[5] != 0 && parameters_.bc_type[5] != 2) {
                        k_top = hydraulicConductivity(state.h[cell], cell_params.ksz, cell_params);
                    }
                } else {
                    k_top = state.kz[zm(cell)];
                }
                const real dzm = k == 0 ? 0.5 * state.dz3d[cell]
                                        : 0.5 * (state.dz3d[cell] + state.dz3d[zm(cell)]);
                const real gzm = -k_top * state.dtg / (state.dz3d[cell] * dzm) * volume;

                real gct = storage;
                gct -= gxp + gxm + gyp + gym;
                if (k < grid_.nz() - 1 || parameters_.bc_type[4] == 1) {
                    gct -= gzp;
                }
                if (k > 0 || parameters_.bc_type[5] == 1) {
                    gct -= gzm;
                }

                real rhs = storage * state.hn[cell];
                rhs -= state.dtg * volume * state.kz[cell] / state.dz3d[cell];
                rhs += state.dtg * volume * k_top / state.dz3d[cell];
                if (k == grid_.nz() - 1 && parameters_.bc_type[4] == 2) {
                    rhs += volume *
                           (parameters_.qbot + state.dtg * state.kz[cell] / state.dz3d[cell]);
                } else if (k == grid_.nz() - 1 && parameters_.bc_type[4] == 1) {
                    rhs += gzp * parameters_.hbot;
                }
                if (k == 0) {
                    if (parameters_.bc_type[5] == 2) {
                        rhs += -state.dtg * volume * (qtopForSurfaceCell(surface_cell) + k_top) /
                               state.dz3d[cell];
                    } else if (parameters_.bc_type[5] == 1) {
                        rhs -= gzm * parameters_.htop;
                    }
                }

                state.gct[cell] = gct;
                state.gxp[cell] = gxp;
                state.gxm[cell] = gxm;
                state.gyp[cell] = gyp;
                state.gym[cell] = gym;
                state.gzp[cell] = gzp;
                state.gzm[cell] = gzm;
                state.grhs[cell] = rhs;
                system.rhs[cell] = rhs;

                if (i > 0) {
                    system.entries.push_back({cell, index(i - 1, j, k), gxm});
                }
                if (j > 0) {
                    system.entries.push_back({cell, index(i, j - 1, k), gym});
                }
                if (k > 0) {
                    system.entries.push_back({cell, zm(cell), gzm});
                }
                system.entries.push_back({cell, cell, gct});
                if (k < grid_.nz() - 1) {
                    system.entries.push_back({cell, zp(cell), gzp});
                }
                if (j < grid_.ny() - 1) {
                    system.entries.push_back({cell, index(i, j + 1, k), gyp});
                }
                if (i < grid_.nx() - 1) {
                    system.entries.push_back({cell, index(i + 1, j, k), gxp});
                }
            }
        }
    }
    return system;
}

std::vector<real> ReSolver::solveLinearSystem(const ReLinearSystem& system) const
{
    std::vector<real> solution(static_cast<std::size_t>(system.n), 0.0);
#ifdef USE_PETSC
    Mat matrix = nullptr;
    Vec rhs = nullptr;
    Vec x = nullptr;
    KSP ksp = nullptr;
    PC pc = nullptr;

    MatCreate(PETSC_COMM_SELF, &matrix);
    MatSetSizes(matrix, PETSC_DECIDE, PETSC_DECIDE, system.n, system.n);
    MatSetFromOptions(matrix);
    MatSeqAIJSetPreallocation(matrix, 7, nullptr);
    MatSetUp(matrix);
    for (const auto& entry : system.entries) {
        MatSetValue(matrix, entry.row, entry.col, entry.value, ADD_VALUES);
    }
    MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);

    VecCreateSeq(PETSC_COMM_SELF, system.n, &rhs);
    VecDuplicate(rhs, &x);
    for (int row = 0; row < system.n; ++row) {
        VecSetValue(rhs, row, system.rhs[static_cast<std::size_t>(row)], INSERT_VALUES);
    }
    VecAssemblyBegin(rhs);
    VecAssemblyEnd(rhs);

    KSPCreate(PETSC_COMM_SELF, &ksp);
    KSPSetOperators(ksp, matrix, matrix);
    KSPSetType(ksp, KSPPREONLY);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCLU);
    KSPSetTolerances(ksp, 1.0e-8, PETSC_DEFAULT, PETSC_DEFAULT, 1000);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, rhs, x);

    const PetscScalar* values = nullptr;
    VecGetArrayRead(x, &values);
    for (int row = 0; row < system.n; ++row) {
        solution[static_cast<std::size_t>(row)] = values[row];
    }
    VecRestoreArrayRead(x, &values);

    KSPDestroy(&ksp);
    VecDestroy(&x);
    VecDestroy(&rhs);
    MatDestroy(&matrix);
#else
    std::vector<real> dense(static_cast<std::size_t>(system.n * system.n), 0.0);
    for (const auto& entry : system.entries) {
        dense[static_cast<std::size_t>(entry.row * system.n + entry.col)] += entry.value;
    }
    solution = system.rhs;
    for (int pivot = 0; pivot < system.n; ++pivot) {
        const real diag = dense[static_cast<std::size_t>(pivot * system.n + pivot)];
        for (int col = pivot; col < system.n; ++col) {
            dense[static_cast<std::size_t>(pivot * system.n + col)] /= diag;
        }
        solution[static_cast<std::size_t>(pivot)] /= diag;
        for (int row = pivot + 1; row < system.n; ++row) {
            const real factor = dense[static_cast<std::size_t>(row * system.n + pivot)];
            for (int col = pivot; col < system.n; ++col) {
                dense[static_cast<std::size_t>(row * system.n + col)] -=
                    factor * dense[static_cast<std::size_t>(pivot * system.n + col)];
            }
            solution[static_cast<std::size_t>(row)] -=
                factor * solution[static_cast<std::size_t>(pivot)];
        }
    }
    for (int row = system.n - 1; row >= 0; --row) {
        for (int col = row + 1; col < system.n; ++col) {
            solution[static_cast<std::size_t>(row)] -=
                dense[static_cast<std::size_t>(row * system.n + col)] *
                solution[static_cast<std::size_t>(col)];
        }
    }
#endif
    return solution;
}

void ReSolver::computeDarcyFlux(LegacyReState& state) const
{
    const real horizontal_area = grid_.dx() * grid_.dy();
    std::fill(state.qz_top_surface.begin(), state.qz_top_surface.end(), 0.0);
    for (int j = 0; j < grid_.ny(); ++j) {
        for (int i = 0; i < grid_.nx(); ++i) {
            const int surface_cell = static_cast<int>(grid_.getSurfaceIndex(i, j));
            for (int k = 0; k < grid_.nz(); ++k) {
                const int cell = index(i, j, k);
                const real x_area = grid_.dy() * state.dz3d[cell];
                const real y_area = grid_.dx() * state.dz3d[cell];
                if (i < grid_.nx() - 1) {
                    const int east = index(i + 1, j, k);
                    state.qx[cell] =
                        x_area * state.kx[cell] * (state.h[east] - state.h[cell]) / grid_.dx();
                } else {
                    state.qx[cell] = 0.0;
                }
                if (j < grid_.ny() - 1) {
                    const int north = index(i, j + 1, k);
                    state.qy[cell] =
                        y_area * state.ky[cell] * (state.h[north] - state.h[cell]) / grid_.dy();
                } else {
                    state.qy[cell] = 0.0;
                }
                if (k == grid_.nz() - 1) {
                    state.qz[cell] = 0.0;
                } else {
                    const int below = zp(cell);
                    const real dzf = 0.5 * (state.dz3d[cell] + state.dz3d[below]);
                    state.qz[cell] =
                        horizontal_area * state.kz[cell] *
                        ((state.h[below] - state.h[cell]) / dzf - 1.0);
                }
            }

            const int top_cell = index(i, j, 0);
            if (parameters_.bc_type[5] == 2) {
                state.qz_top_surface[surface_cell] =
                    horizontal_area * qtopForSurfaceCell(surface_cell);
            } else if (parameters_.bc_type[5] == 1) {
                const auto top_params = parametersForCell(state, top_cell);
                state.qz_top_surface[surface_cell] =
                    horizontal_area * top_params.ksz *
                    ((state.h[top_cell] - parameters_.htop) / (0.5 * state.dz3d[top_cell]) - 1.0);
            }
        }
    }
    state.qz_top = state.qz_top_surface.empty() ? 0.0 : state.qz_top_surface.front();
}

void ReSolver::updateWaterContent(LegacyReState& state) const
{
    const real horizontal_area = grid_.dx() * grid_.dy();
    for (int j = 0; j < grid_.ny(); ++j) {
        for (int i = 0; i < grid_.nx(); ++i) {
            const int surface_cell = static_cast<int>(grid_.getSurfaceIndex(i, j));
            for (int k = 0; k < grid_.nz(); ++k) {
                const int cell = index(i, j, k);
                const auto cell_params = parametersForCell(state, cell);
                const real volume = horizontal_area * state.dz3d[cell];
                const real coeff =
                    volume *
                    (1.0 + parameters_.ss * (state.h[cell] - state.hn[cell]) / cell_params.wcs);
                const real qx_plus = i < grid_.nx() - 1 ? state.qx[cell] : 0.0;
                const real qx_minus = i > 0 ? state.qx[index(i - 1, j, k)] : 0.0;
                const real qy_plus = j < grid_.ny() - 1 ? state.qy[cell] : 0.0;
                const real qy_minus = j > 0 ? state.qy[index(i, j - 1, k)] : 0.0;
                const real q_down = state.qz[cell];
                const real q_top =
                    k == 0 ? state.qz_top_surface[surface_cell] : state.qz[zm(cell)];
                const real dq = state.dtg *
                                (qx_plus - qx_minus + qy_plus - qy_minus + q_down - q_top);
                state.wc[cell] = (state.wcn[cell] * volume + dq) / coeff;
            }
        }
    }
}

void ReSolver::enforceLegacyMoistureConsistency(LegacyReState& state) const
{
    for (int cell = 0; cell < static_cast<int>(state.wc.size()); ++cell) {
        const auto cell_params = parametersForCell(state, cell);
        state.wch[cell] = waterContentFromHead(state.h[cell], cell_params);
        state.hwc[cell] = headFromWaterContent(state.wc[cell], cell_params);
        if (state.wc[cell] >= cell_params.wcs) {
            state.wc[cell] = cell_params.wcs;
        } else if (state.wc[cell] <= cell_params.wcr + WC_HEAD_EPS) {
            state.wc[cell] = cell_params.wcr;
        } else if (isAdjacentToSaturation(state, cell)) {
            state.wc[cell] = state.wch[cell];
        } else if (state.wc[cell] < 0.9999 * cell_params.wcs) {
            state.h[cell] = state.hwc[cell];
        }
    }
}

real ReSolver::adaptTimeStep(const LegacyReState& state) const
{
    real dq_max = 0.0;
    real dt_co_min = 1.0e8;
    for (int j = 0; j < grid_.ny(); ++j) {
        for (int i = 0; i < grid_.nx(); ++i) {
            for (int k = 0; k < grid_.nz(); ++k) {
                const int cell = index(i, j, k);
                if (isTop(cell)) {
                    continue;
                }
                const real area = grid_.dx() * grid_.dy();
                const real qin = state.qz[cell] / area;
                const real qout = state.qz[zm(cell)] / area;
                dq_max = std::max(dq_max, std::abs(qin - qout) * state.dtg / state.dz3d[cell]);
                const auto cell_params = parametersForCell(state, cell);
                if (state.wc[cell] < cell_params.wcs) {
                    const real dkdwc =
                        conductivityDerivativeWc(state.wc[cell], cell_params.ksz, cell_params);
                    if (dkdwc > 0.0) {
                        dt_co_min =
                            std::min(dt_co_min, parameters_.co_max * state.dz3d[cell] / dkdwc);
                    }
                }
            }
        }
    }

    real next_dt = state.dtg;
    if (dq_max > 0.02) {
        next_dt *= 0.75;
    } else if (dq_max < 0.01) {
        next_dt *= 1.25;
    }
    next_dt = std::min(next_dt, dt_co_min);
    next_dt = std::min(next_dt, parameters_.dt_max);
    next_dt = std::max(next_dt, parameters_.dt_min);
    return next_dt;
}

int ReSolver::index(int i, int j, int k) const
{
    return static_cast<int>(grid_.getIndex(i, j, k));
}

int ReSolver::zp(int cell) const
{
    return cell + 1;
}

int ReSolver::zm(int cell) const
{
    return cell - 1;
}

real ReSolver::qtopForSurfaceCell(int surface_cell) const
{
    if (parameters_.qtop_surface.empty()) {
        return parameters_.qtop;
    }
    return parameters_.qtop_surface[static_cast<std::size_t>(surface_cell)];
}

const SoilParameters& ReSolver::soilForCell(const LegacyReState& state, int cell) const
{
    static const SoilParameters fallback{};
    if (cell < 0 || cell >= static_cast<int>(state.soil_id.size())) {
        throw std::out_of_range("soil cell index is out of range");
    }
    if (parameters_.soil_table.empty()) {
        return fallback;
    }
    const int soil_id = state.soil_id[static_cast<std::size_t>(cell)];
    if (soil_id < 0 || soil_id >= static_cast<int>(parameters_.soil_table.size())) {
        throw std::out_of_range("soil id is outside the VG table");
    }
    return parameters_.soil_table[static_cast<std::size_t>(soil_id)];
}

ReParameters ReSolver::parametersForCell(const LegacyReState& state, int cell) const
{
    if (parameters_.soil_table.empty()) {
        return parameters_;
    }
    return parametersForSoil(parameters_, soilForCell(state, cell));
}

bool ReSolver::isTop(int cell) const
{
    return cell % grid_.nz() == 0;
}

bool ReSolver::isAdjacentToSaturation(const LegacyReState& state, int cell) const
{
    const int k = cell % grid_.nz();
    if (k < grid_.nz() - 1) {
        const int below = zp(cell);
        if (state.wc[below] >= parametersForCell(state, below).wcs) {
            return true;
        }
    }
    if (k > 0) {
        const int above = zm(cell);
        if (state.wc[above] >= parametersForCell(state, above).wcs) {
            return true;
        }
    }
    if (state.wc[cell] >= parametersForCell(state, cell).wcs) {
        return true;
    }
    if (k == 0) {
        return parameters_.bc_type[5] == 1 && parameters_.htop >= 0.0;
    }
    return false;
}

}  // namespace frehg2
