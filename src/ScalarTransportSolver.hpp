#ifndef FREHG_SCALAR_TRANSPORT_SOLVER_HPP
#define FREHG_SCALAR_TRANSPORT_SOLVER_HPP

#include "define.hpp"
#include "Domain.hpp"
#include "StateVariables.hpp"
#include "ActiveCellMesh.hpp"
#include <memory>
#include <cmath>

// ============================================================================
//                      SCALAR TRANSPORT SOLVER
// ============================================================================
// Solves scalar transport (advection-diffusion/dispersion) equations
// Uses TVD scheme for advection and handles coupling between domains

// Forward declarations
class SwDomain;
class GwDomain;
class SwStateVariables;
class GwStateVariables;

// ============================================================================
//                      TVD SUPERBEE LIMITER
// ============================================================================
// 2nd-order TVD scheme with superbee limiter for advection

KOKKOS_INLINE_FUNCTION
Scalar tvd_superbee(Scalar sp, Scalar sc, Scalar sm, Scalar u, Scalar delta, Scalar dt) {
    Scalar phi = 0.0;
    Scalar r1 = 1.0;
    Scalar r2 = 2.0;
    Scalar coef = std::abs(u) * dt / delta;
    
    if (sp != sc) {
        Scalar r = (sc - sm) / (sp - sc);
        if (2.0 * r < 1.0) {
            r1 = 2.0 * r;
        }
        if (r < 2.0) {
            r2 = r;
        }
        if (r1 > 0.0 || r2 > 0.0) {
            if (r1 > r2) {
                phi = r1;
            } else {
                phi = r2;
            }
        }
    }
    
    return sc + 0.5 * phi * (1.0 - coef) * (sp - sc);
}

// ============================================================================
//                      BASE SCALAR TRANSPORT SOLVER
// ============================================================================
// Common functionality for both surface water and groundwater scalar transport

class BaseScalarTransportSolver {
public:
    // --- Scalar Data ---
    // Scalar concentration (per scalar index)
    // For multiple scalars, we'll use a 2D view: View2D<Scalar> scalar_concentration(num_scalars, num_cells)
    // For now, we'll assume a single scalar and use View1D
    
    // --- Parameters ---
    Scalar dt;                  // Time step
    Scalar s_lim_hi;               // Upper limit for scalar concentration
    Scalar s_lim_lo;               // Lower limit for scalar concentration
    bool use_tvd;                  // Use TVD scheme for advection
    
    // Constructor
    BaseScalarTransportSolver(Scalar _dt, 
                              Scalar _s_lim_hi = 200.0,
                              Scalar _s_lim_lo = 0.0,
                              bool _use_tvd = true)
        : dt(_dt), s_lim_hi(_s_lim_hi), s_lim_lo(_s_lim_lo), use_tvd(_use_tvd) {}
    
    virtual ~BaseScalarTransportSolver() = default;
    
    // Apply scalar limiter to keep values within bounds
    KOKKOS_INLINE_FUNCTION
    Scalar apply_limiter(Scalar s, Scalar s_min, Scalar s_max) const {
        if (s > s_max && s_max < s_lim_hi) {
            return s_max;
        } else if (s < s_min && s_min > s_lim_lo) {
            return s_min;
        }
        if (s > s_lim_hi) {
            return s_lim_hi;
        } else if (s < s_lim_lo) {
            return s_lim_lo;
        }
        return s;
    }
};

// ============================================================================
//                      SURFACE WATER SCALAR TRANSPORT (2D)
// ============================================================================

class SwScalarTransportSolver : public BaseScalarTransportSolver {
public:
    // --- Domain and State ---
    SwDomain& domain;
    ActiveCellMesh& active_mesh;
    SwStateVariables& state;
    
    // --- Scalar Arrays ---
    View1D<Scalar> scalar_concentration;      // Scalar concentration (s_surf)
    View1D<Scalar> scalar_mass;              // Scalar mass (sm_surf)
    View1D<Scalar> scalar_old;               // Old scalar concentration
    
    // --- Diffusion Coefficients ---
    Scalar diff_x;                            // Diffusion coefficient in x-direction
    Scalar diff_y;                            // Diffusion coefficient in y-direction
    
    // Constructor
    SwScalarTransportSolver(SwDomain& _domain,
                           ActiveCellMesh& _active_mesh,
                           SwStateVariables& _state,
                           Scalar _dt,
                           Scalar _diff_x = 0.0,
                           Scalar _diff_y = 0.0,
                           Scalar _s_lim_hi = 200.0,
                           Scalar _s_lim_lo = 0.0,
                           bool _use_tvd = true)
        : BaseScalarTransportSolver(_dt, _s_lim_hi, _s_lim_lo, _use_tvd),
          domain(_domain), active_mesh(_active_mesh), state(_state),
          diff_x(_diff_x), diff_y(_diff_y) {
        
        // Allocate scalar arrays
        scalar_concentration = View1D<Scalar>("scalar_concentration", domain.num_cells_total);
        scalar_mass = View1D<Scalar>("scalar_mass", domain.num_cells_total);
        scalar_old = View1D<Scalar>("scalar_old", domain.num_cells_total);
    }
    
    // Main solver function
    void solve() {
        // Save old values
        Kokkos::deep_copy(scalar_old, scalar_concentration);
        
        // Compute scalar mass from concentration
        compute_scalar_mass();
        
        // Advective transport
        compute_advection();
        
        // Diffusive transport
        compute_diffusion();
        
        // Update concentration from mass
        update_concentration();
        
        // Apply limiters
        apply_scalar_limiters();
    }
    
private:
    // Compute scalar mass from concentration
    void compute_scalar_mass() {
        auto _scalar_concentration = scalar_concentration;
        auto _scalar_mass = scalar_mass;
        auto _volume_old = state.volume_old;
        auto _active_mask = domain.active_mask;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask(domain_idx) > 0) {
                    _scalar_mass(domain_idx) = _scalar_concentration(domain_idx) * _volume_old(domain_idx);
                } else {
                    _scalar_mass(domain_idx) = 0.0;
                }
            });
    }
    
    // Compute advective transport using TVD scheme
    void compute_advection() {
        auto _scalar_concentration = scalar_concentration;
        auto _scalar_mass = scalar_mass;
        auto _flux_x = state.flux_x;
        auto _flux_y = state.flux_y;
        auto _area_x = state.area_x;
        auto _area_y = state.area_y;
        auto _active_mask = domain.active_mask;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask(domain_idx) == 0) {
                    return;
                }
                
                // Get neighbors
                Ordinal active_left = _neighbor_left(i);
                Ordinal active_right = _neighbor_right(i);
                Ordinal active_back = _neighbor_back(i);
                Ordinal active_front = _neighbor_front(i);
                
                Ordinal i_left = (active_left >= 0) ? _active_to_domain(active_left) : -1;
                Ordinal i_right = (active_right >= 0) ? _active_to_domain(active_right) : -1;
                Ordinal i_back = (active_back >= 0) ? _active_to_domain(active_back) : -1;
                Ordinal i_front = (active_front >= 0) ? _active_to_domain(active_front) : -1;
                
                // X-direction advection
                Scalar sip = 0.0;  // Scalar at i+1/2 face
                Scalar sim = 0.0;  // Scalar at i-1/2 face
                
                // i+1/2 face
                Scalar flux_x_here = _flux_x(domain_idx);
                if (flux_x_here > 0.0) {
                    // Flow from left to right
                    sip = _scalar_concentration(domain_idx);
                    if (use_tvd && i_right >= 0 && i_left >= 0) {
                        sip = tvd_superbee(_scalar_concentration(i_right), 
                                          _scalar_concentration(domain_idx),
                                          _scalar_concentration(i_left),
                                          flux_x_here / _area_x(domain_idx),
                                          domain.dx, dt);
                    }
                } else if (flux_x_here < 0.0) {
                    // Flow from right to left
                    if (i_right >= 0) {
                        sip = _scalar_concentration(i_right);
                        if (use_tvd && i_right >= 0 && i_left >= 0) {
                            Ordinal i_right_right = -1;
                            if (active_right >= 0) {
                                Ordinal active_right_right = _neighbor_right(active_right);
                                if (active_right_right >= 0) {
                                    i_right_right = _active_to_domain(active_right_right);
                                }
                            }
                            if (i_right_right < 0) i_right_right = i_right;
                            sip = tvd_superbee(_scalar_concentration(domain_idx),
                                              _scalar_concentration(i_right),
                                              _scalar_concentration(i_right_right),
                                              flux_x_here / _area_x(domain_idx),
                                              domain.dx, dt);
                        }
                    } else {
                        sip = _scalar_concentration(domain_idx);
                    }
                }
                
                // i-1/2 face
                if (i_left >= 0) {
                    Scalar flux_x_left = _flux_x(i_left);
                    if (flux_x_left > 0.0) {
                        sim = _scalar_concentration(i_left);
                        if (use_tvd && i_left >= 0 && i_right >= 0) {
                            sim = tvd_superbee(_scalar_concentration(domain_idx),
                                              _scalar_concentration(i_left),
                                              _scalar_concentration(i_right),
                                              flux_x_left / _area_x(i_left),
                                              domain.dx, dt);
                        }
                    } else if (flux_x_left < 0.0) {
                        sim = _scalar_concentration(domain_idx);
                        if (use_tvd && i_left >= 0) {
                            Ordinal i_left_left = -1;
                            if (active_left >= 0) {
                                Ordinal active_left_left = _neighbor_left(active_left);
                                if (active_left_left >= 0) {
                                    i_left_left = _active_to_domain(active_left_left);
                                }
                            }
                            if (i_left_left < 0) i_left_left = i_left;
                            sim = tvd_superbee(_scalar_concentration(i_left),
                                              _scalar_concentration(domain_idx),
                                              _scalar_concentration(i_left_left),
                                              flux_x_left / _area_x(i_left),
                                              domain.dx, dt);
                        }
                    }
                } else {
                    sim = _scalar_concentration(domain_idx);
                }
                
                // Y-direction advection
                Scalar sjp = 0.0;  // Scalar at j+1/2 face
                Scalar sjm = 0.0;  // Scalar at j-1/2 face
                
                // j+1/2 face
                Scalar flux_y_here = _flux_y(domain_idx);
                if (flux_y_here > 0.0) {
                    sjp = _scalar_concentration(domain_idx);
                    if (use_tvd && i_front >= 0 && i_back >= 0) {
                        sjp = tvd_superbee(_scalar_concentration(i_front),
                                          _scalar_concentration(domain_idx),
                                          _scalar_concentration(i_back),
                                          flux_y_here / _area_y(domain_idx),
                                          domain.dy, dt);
                    }
                } else if (flux_y_here < 0.0) {
                    if (i_front >= 0) {
                        sjp = _scalar_concentration(i_front);
                        if (use_tvd && i_front >= 0 && i_back >= 0) {
                            Ordinal i_front_front = -1;
                            if (active_front >= 0) {
                                Ordinal active_front_front = _neighbor_front(active_front);
                                if (active_front_front >= 0) {
                                    i_front_front = _active_to_domain(active_front_front);
                                }
                            }
                            if (i_front_front < 0) i_front_front = i_front;
                            sjp = tvd_superbee(_scalar_concentration(domain_idx),
                                              _scalar_concentration(i_front),
                                              _scalar_concentration(i_front_front),
                                              flux_y_here / _area_y(domain_idx),
                                              domain.dy, dt);
                        }
                    } else {
                        sjp = _scalar_concentration(domain_idx);
                    }
                }
                
                // j-1/2 face
                if (i_back >= 0) {
                    Scalar flux_y_back = _flux_y(i_back);
                    if (flux_y_back > 0.0) {
                        sjm = _scalar_concentration(i_back);
                        if (use_tvd && i_back >= 0 && i_front >= 0) {
                            sjm = tvd_superbee(_scalar_concentration(domain_idx),
                                              _scalar_concentration(i_back),
                                              _scalar_concentration(i_front),
                                              flux_y_back / _area_y(i_back),
                                              domain.dy, dt);
                        }
                    } else if (flux_y_back < 0.0) {
                        sjm = _scalar_concentration(domain_idx);
                        if (use_tvd && i_back >= 0) {
                            Ordinal i_back_back = -1;
                            if (active_back >= 0) {
                                Ordinal active_back_back = _neighbor_back(active_back);
                                if (active_back_back >= 0) {
                                    i_back_back = _active_to_domain(active_back_back);
                                }
                            }
                            if (i_back_back < 0) i_back_back = i_back;
                            sjm = tvd_superbee(_scalar_concentration(i_back),
                                              _scalar_concentration(domain_idx),
                                              _scalar_concentration(i_back_back),
                                              flux_y_back / _area_y(i_back),
                                              domain.dy, dt);
                        }
                    }
                } else {
                    sjm = _scalar_concentration(domain_idx);
                }
                
                // Update scalar mass from advection
                _scalar_mass(domain_idx) += dt * (-flux_x_here * sip + 
                                                  (i_left >= 0 ? _flux_x(i_left) : 0.0) * sim -
                                                  flux_y_here * sjp + 
                                                  (i_back >= 0 ? _flux_y(i_back) : 0.0) * sjm);
            });
    }
    
    // Compute diffusive transport
    void compute_diffusion() {
        auto _scalar_concentration = scalar_concentration;
        auto _scalar_mass = scalar_mass;
        auto _area_x = state.area_x;
        auto _area_y = state.area_y;
        auto _active_mask = domain.active_mask;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask(domain_idx) == 0) {
                    return;
                }
                
                Ordinal active_left = _neighbor_left(i);
                Ordinal active_right = _neighbor_right(i);
                Ordinal active_back = _neighbor_back(i);
                Ordinal active_front = _neighbor_front(i);
                
                Ordinal i_left = (active_left >= 0) ? _active_to_domain(active_left) : -1;
                Ordinal i_right = (active_right >= 0) ? _active_to_domain(active_right) : -1;
                Ordinal i_back = (active_back >= 0) ? _active_to_domain(active_back) : -1;
                Ordinal i_front = (active_front >= 0) ? _active_to_domain(active_front) : -1;
                
                // X-direction diffusion
                Scalar diff_x_flux = 0.0;
                if (i_right >= 0) {
                    diff_x_flux += (diff_x * _area_x(domain_idx) / domain.dx) * 
                                  (_scalar_concentration(i_right) - _scalar_concentration(domain_idx));
                }
                if (i_left >= 0) {
                    diff_x_flux -= (diff_x * _area_x(i_left) / domain.dx) * 
                                  (_scalar_concentration(domain_idx) - _scalar_concentration(i_left));
                }
                
                // Y-direction diffusion
                Scalar diff_y_flux = 0.0;
                if (i_front >= 0) {
                    diff_y_flux += (diff_y * _area_y(domain_idx) / domain.dy) * 
                                  (_scalar_concentration(i_front) - _scalar_concentration(domain_idx));
                }
                if (i_back >= 0) {
                    diff_y_flux -= (diff_y * _area_y(i_back) / domain.dy) * 
                                  (_scalar_concentration(domain_idx) - _scalar_concentration(i_back));
                }
                
                // Update scalar mass from diffusion
                _scalar_mass(domain_idx) += dt * (diff_x_flux + diff_y_flux);
            });
    }
    
    // Update concentration from mass
    void update_concentration() {
        auto _scalar_concentration = scalar_concentration;
        auto _scalar_mass = scalar_mass;
        auto _volume = state.volume;
        auto _active_mask = domain.active_mask;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask(domain_idx) > 0) {
                    if (_volume(domain_idx) > 0.0) {
                        _scalar_concentration(domain_idx) = _scalar_mass(domain_idx) / _volume(domain_idx);
                    } else {
                        _scalar_concentration(domain_idx) = 0.0;
                    }
                } else {
                    _scalar_concentration(domain_idx) = 0.0;
                }
            });
    }
    
    // Apply scalar limiters
    void apply_scalar_limiters() {
        auto _scalar_concentration = scalar_concentration;
        auto _scalar_old = scalar_old;
        auto _area_x = state.area_x;
        auto _active_mask = domain.active_mask;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _active_to_domain = active_mesh.active_to_domain;
        
        // First pass: find min/max in neighborhood
        View1D<Scalar> s_min("s_min", domain.num_cells_total);
        View1D<Scalar> s_max("s_max", domain.num_cells_total);
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask(domain_idx) == 0) {
                    s_min(domain_idx) = s_lim_hi;
                    s_max(domain_idx) = s_lim_lo;
                    return;
                }
                
                s_min(domain_idx) = _scalar_concentration(domain_idx);
                s_max(domain_idx) = _scalar_concentration(domain_idx);
                
                Ordinal active_right = _neighbor_right(i);
                Ordinal active_left = _neighbor_left(i);
                Ordinal active_front = _neighbor_front(i);
                Ordinal active_back = _neighbor_back(i);
                
                if (active_right >= 0 && _area_x(domain_idx) > 0.0) {
                    Ordinal i_right = _active_to_domain(active_right);
                    if (_scalar_concentration(i_right) > s_max(domain_idx)) {
                        s_max(domain_idx) = _scalar_concentration(i_right);
                    }
                    if (_scalar_concentration(i_right) < s_min(domain_idx)) {
                        s_min(domain_idx) = _scalar_concentration(i_right);
                    }
                }
                if (active_left >= 0 && state.area_x(_active_to_domain(active_left)) > 0.0) {
                    Ordinal i_left = _active_to_domain(active_left);
                    if (_scalar_concentration(i_left) > s_max(domain_idx)) {
                        s_max(domain_idx) = _scalar_concentration(i_left);
                    }
                    if (_scalar_concentration(i_left) < s_min(domain_idx)) {
                        s_min(domain_idx) = _scalar_concentration(i_left);
                    }
                }
                if (active_front >= 0 && state.area_y(domain_idx) > 0.0) {
                    Ordinal i_front = _active_to_domain(active_front);
                    if (_scalar_concentration(i_front) > s_max(domain_idx)) {
                        s_max(domain_idx) = _scalar_concentration(i_front);
                    }
                    if (_scalar_concentration(i_front) < s_min(domain_idx)) {
                        s_min(domain_idx) = _scalar_concentration(i_front);
                    }
                }
                if (active_back >= 0 && state.area_y(_active_to_domain(active_back)) > 0.0) {
                    Ordinal i_back = _active_to_domain(active_back);
                    if (_scalar_concentration(i_back) > s_max(domain_idx)) {
                        s_max(domain_idx) = _scalar_concentration(i_back);
                    }
                    if (_scalar_concentration(i_back) < s_min(domain_idx)) {
                        s_min(domain_idx) = _scalar_concentration(i_back);
                    }
                }
            });
        
        // Second pass: apply limiter
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask(domain_idx) > 0) {
                    _scalar_concentration(domain_idx) = apply_limiter(_scalar_concentration(domain_idx),
                                                                       s_min(domain_idx), s_max(domain_idx));
                }
            });
    }
};

// ============================================================================
//                      GROUNDWATER SCALAR TRANSPORT (3D)
// ============================================================================

class GwScalarTransportSolver : public BaseScalarTransportSolver {
public:
    // --- Domain and State ---
    GwDomain& domain;
    ActiveCellMesh& active_mesh;
    GwStateVariables& state;
    
    // --- Scalar Arrays ---
    View1D<Scalar> scalar_concentration;      // Scalar concentration (s_subs)
    View1D<Scalar> scalar_mass;              // Scalar mass (sm_subs)
    View1D<Scalar> scalar_old;               // Old scalar concentration
    
    // --- Dispersion Tensor Components ---
    View1D<Scalar> Dxx, Dxy, Dxz;  // Dispersion tensor components
    View1D<Scalar> Dyx, Dyy, Dyz;
    View1D<Scalar> Dzx, Dzy, Dzz;
    
    // --- Dispersion Parameters ---
    Scalar diff_x;                  // Base diffusion coefficient in x-direction
    Scalar diff_y;                  // Base diffusion coefficient in y-direction
    Scalar diff_z;                  // Base diffusion coefficient in z-direction
    Scalar disp_longitudinal;      // Longitudinal dispersivity
    Scalar disp_transverse;         // Transverse dispersivity
    
    // Constructor
    GwScalarTransportSolver(GwDomain& _domain,
                           ActiveCellMesh& _active_mesh,
                           GwStateVariables& _state,
                           Scalar _dt,
                           Scalar _diff_x = 0.0,
                           Scalar _diff_y = 0.0,
                           Scalar _diff_z = 0.0,
                           Scalar _disp_longitudinal = 0.0,
                           Scalar _disp_transverse = 0.0,
                           Scalar _s_lim_hi = 200.0,
                           Scalar _s_lim_lo = 0.0,
                           bool _use_tvd = true)
        : BaseScalarTransportSolver(_dt, _s_lim_hi, _s_lim_lo, _use_tvd),
          domain(_domain), active_mesh(_active_mesh), state(_state),
          diff_x(_diff_x), diff_y(_diff_y), diff_z(_diff_z),
          disp_longitudinal(_disp_longitudinal), disp_transverse(_disp_transverse) {
        
        // Allocate scalar arrays
        Ordinal num_cells_3d = domain.num_cells_total * domain.nz;
        scalar_concentration = View1D<Scalar>("scalar_concentration", num_cells_3d);
        scalar_mass = View1D<Scalar>("scalar_mass", num_cells_3d);
        scalar_old = View1D<Scalar>("scalar_old", num_cells_3d);
        
        // Allocate dispersion tensor arrays
        Dxx = View1D<Scalar>("Dxx", num_cells_3d);
        Dxy = View1D<Scalar>("Dxy", num_cells_3d);
        Dxz = View1D<Scalar>("Dxz", num_cells_3d);
        Dyx = View1D<Scalar>("Dyx", num_cells_3d);
        Dyy = View1D<Scalar>("Dyy", num_cells_3d);
        Dyz = View1D<Scalar>("Dyz", num_cells_3d);
        Dzx = View1D<Scalar>("Dzx", num_cells_3d);
        Dzy = View1D<Scalar>("Dzy", num_cells_3d);
        Dzz = View1D<Scalar>("Dzz", num_cells_3d);
    }
    
    // Main solver function
    void solve() {
        // Save old values
        Kokkos::deep_copy(scalar_old, scalar_concentration);
        
        // Compute scalar mass from concentration
        compute_scalar_mass();
        
        // Advective transport
        compute_advection();
        
        // Compute dispersion tensor
        compute_dispersion_tensor();
        
        // Dispersive transport
        compute_dispersion();
        
        // Update concentration from mass
        update_concentration();
        
        // Apply limiters (to be implemented later)
        // apply_scalar_limiters();
    }
    
private:
    // Compute scalar mass from concentration
    void compute_scalar_mass() {
        auto _scalar_concentration = scalar_concentration;
        auto _scalar_mass = scalar_mass;
        auto _volume_old = state.volume_old;
        auto _active_mask_3d = domain.active_mask_3d;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) > 0) {
                    _scalar_mass(domain_idx) = _scalar_concentration(domain_idx) * _volume_old(domain_idx);
                } else {
                    _scalar_mass(domain_idx) = 0.0;
                }
            });
    }
    
    // Compute advective transport using TVD scheme (3D)
    void compute_advection() {
        auto _scalar_concentration = scalar_concentration;
        auto _scalar_mass = scalar_mass;
        auto _flux_x = state.flux_x;
        auto _flux_y = state.flux_y;
        auto _flux_z = state.flux_z;
        auto _active_mask_3d = domain.active_mask_3d;
        auto _dz_layers = domain.dz_layers;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _neighbor_bottom = active_mesh.neighbor_bottom;
        auto _neighbor_top = active_mesh.neighbor_top;
        auto _active_to_domain = active_mesh.active_to_domain;
        auto _coord_k = active_mesh.coord_k;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) {
                    return;
                }
                
                // Get neighbors
                Ordinal active_left = _neighbor_left(i);
                Ordinal active_right = _neighbor_right(i);
                Ordinal active_back = _neighbor_back(i);
                Ordinal active_front = _neighbor_front(i);
                Ordinal active_bottom = _neighbor_bottom(i);
                Ordinal active_top = _neighbor_top(i);
                
                Ordinal i_left = (active_left >= 0) ? _active_to_domain(active_left) : -1;
                Ordinal i_right = (active_right >= 0) ? _active_to_domain(active_right) : -1;
                Ordinal i_back = (active_back >= 0) ? _active_to_domain(active_back) : -1;
                Ordinal i_front = (active_front >= 0) ? _active_to_domain(active_front) : -1;
                Ordinal i_bottom = (active_bottom >= 0) ? _active_to_domain(active_bottom) : -1;
                Ordinal i_top = (active_top >= 0) ? _active_to_domain(active_top) : -1;
                
                // Get layer index and thickness
                Ordinal k = _coord_k(i);
                Scalar dz = _dz_layers(k);
                
                // Face areas
                Scalar Ax = domain.dy * dz;  // X-face area
                Scalar Ay = domain.dx * dz;  // Y-face area
                Scalar Az = domain.dx * domain.dy;  // Z-face area
                
                // X-direction advection
                Scalar sip = 0.0;  // Scalar at i+1/2 face
                Scalar sim = 0.0;  // Scalar at i-1/2 face
                
                // i+1/2 face
                Scalar flux_x_here = _flux_x(domain_idx);
                if (flux_x_here > 0.0) {
                    // Flow from left to right
                    sip = _scalar_concentration(domain_idx);
                    if (use_tvd && i_right >= 0 && i_left >= 0) {
                        sip = tvd_superbee(_scalar_concentration(i_right), 
                                          _scalar_concentration(domain_idx),
                                          _scalar_concentration(i_left),
                                          flux_x_here / Ax,
                                          domain.dx, dt);
                    }
                } else if (flux_x_here < 0.0) {
                    // Flow from right to left
                    if (i_right >= 0) {
                        sip = _scalar_concentration(i_right);
                        if (use_tvd && i_right >= 0 && i_left >= 0) {
                            // Get i+2 neighbor for TVD
                            Ordinal i_right_right = -1;
                            if (active_right >= 0) {
                                Ordinal active_right_right = _neighbor_right(active_right);
                                if (active_right_right >= 0) {
                                    i_right_right = _active_to_domain(active_right_right);
                                }
                            }
                            if (i_right_right >= 0) {
                                sip = tvd_superbee(_scalar_concentration(domain_idx),
                                                  _scalar_concentration(i_right),
                                                  _scalar_concentration(i_right_right),
                                                  flux_x_here / Ax,
                                                  domain.dx, dt);
                            } else {
                                sip = _scalar_concentration(i_right);
                            }
                        }
                    } else {
                        sip = _scalar_concentration(domain_idx);
                    }
                }
                
                // i-1/2 face
                if (i_left >= 0) {
                    Scalar flux_x_left = _flux_x(i_left);
                    if (flux_x_left > 0.0) {
                        sim = _scalar_concentration(i_left);
                        if (use_tvd && i_left >= 0 && i_right >= 0) {
                            sim = tvd_superbee(_scalar_concentration(domain_idx),
                                              _scalar_concentration(i_left),
                                              _scalar_concentration(i_right),
                                              flux_x_left / Ax,
                                              domain.dx, dt);
                        }
                    } else if (flux_x_left < 0.0) {
                        sim = _scalar_concentration(domain_idx);
                        if (use_tvd && i_left >= 0) {
                            // Get i-2 neighbor for TVD
                            Ordinal i_left_left = -1;
                            if (active_left >= 0) {
                                Ordinal active_left_left = _neighbor_left(active_left);
                                if (active_left_left >= 0) {
                                    i_left_left = _active_to_domain(active_left_left);
                                }
                            }
                            if (i_left_left >= 0) {
                                sim = tvd_superbee(_scalar_concentration(i_left),
                                                  _scalar_concentration(domain_idx),
                                                  _scalar_concentration(i_left_left),
                                                  flux_x_left / Ax,
                                                  domain.dx, dt);
                            } else {
                                sim = _scalar_concentration(domain_idx);
                            }
                        }
                    }
                } else {
                    sim = _scalar_concentration(domain_idx);
                }
                
                // Y-direction advection
                Scalar sjp = 0.0;  // Scalar at j+1/2 face
                Scalar sjm = 0.0;  // Scalar at j-1/2 face
                
                // j+1/2 face
                Scalar flux_y_here = _flux_y(domain_idx);
                if (flux_y_here > 0.0) {
                    sjp = _scalar_concentration(domain_idx);
                    if (use_tvd && i_front >= 0 && i_back >= 0) {
                        sjp = tvd_superbee(_scalar_concentration(i_front),
                                          _scalar_concentration(domain_idx),
                                          _scalar_concentration(i_back),
                                          flux_y_here / Ay,
                                          domain.dy, dt);
                    }
                } else if (flux_y_here < 0.0) {
                    if (i_front >= 0) {
                        sjp = _scalar_concentration(i_front);
                        if (use_tvd && i_front >= 0 && i_back >= 0) {
                            // Get j+2 neighbor for TVD
                            Ordinal i_front_front = -1;
                            if (active_front >= 0) {
                                Ordinal active_front_front = _neighbor_front(active_front);
                                if (active_front_front >= 0) {
                                    i_front_front = _active_to_domain(active_front_front);
                                }
                            }
                            if (i_front_front >= 0) {
                                sjp = tvd_superbee(_scalar_concentration(domain_idx),
                                                  _scalar_concentration(i_front),
                                                  _scalar_concentration(i_front_front),
                                                  flux_y_here / Ay,
                                                  domain.dy, dt);
                            } else {
                                sjp = _scalar_concentration(i_front);
                            }
                        }
                    } else {
                        sjp = _scalar_concentration(domain_idx);
                    }
                }
                
                // j-1/2 face
                if (i_back >= 0) {
                    Scalar flux_y_back = _flux_y(i_back);
                    if (flux_y_back > 0.0) {
                        sjm = _scalar_concentration(i_back);
                        if (use_tvd && i_back >= 0 && i_front >= 0) {
                            sjm = tvd_superbee(_scalar_concentration(domain_idx),
                                              _scalar_concentration(i_back),
                                              _scalar_concentration(i_front),
                                              flux_y_back / Ay,
                                              domain.dy, dt);
                        }
                    } else if (flux_y_back < 0.0) {
                        sjm = _scalar_concentration(domain_idx);
                        if (use_tvd && i_back >= 0) {
                            // Get j-2 neighbor for TVD
                            Ordinal i_back_back = -1;
                            if (active_back >= 0) {
                                Ordinal active_back_back = _neighbor_back(active_back);
                                if (active_back_back >= 0) {
                                    i_back_back = _active_to_domain(active_back_back);
                                }
                            }
                            if (i_back_back >= 0) {
                                sjm = tvd_superbee(_scalar_concentration(i_back),
                                                  _scalar_concentration(domain_idx),
                                                  _scalar_concentration(i_back_back),
                                                  flux_y_back / Ay,
                                                  domain.dy, dt);
                            } else {
                                sjm = _scalar_concentration(domain_idx);
                            }
                        }
                    }
                } else {
                    sjm = _scalar_concentration(domain_idx);
                }
                
                // Z-direction advection
                Scalar skp = 0.0;  // Scalar at k+1/2 face
                Scalar skm = 0.0;  // Scalar at k-1/2 face
                
                // k+1/2 face
                Scalar flux_z_here = _flux_z(domain_idx);
                if (flux_z_here > 0.0) {
                    if (i_top >= 0) {
                        skp = _scalar_concentration(i_top);
                        if (use_tvd && i_top >= 0 && i_bottom >= 0) {
                            // Get k+2 neighbor for TVD
                            Ordinal i_top_top = -1;
                            if (active_top >= 0) {
                                Ordinal active_top_top = _neighbor_top(active_top);
                                if (active_top_top >= 0) {
                                    i_top_top = _active_to_domain(active_top_top);
                                }
                            }
                            if (i_top_top >= 0) {
                                skp = tvd_superbee(_scalar_concentration(domain_idx),
                                                  _scalar_concentration(i_top),
                                                  _scalar_concentration(i_top_top),
                                                  flux_z_here / Az,
                                                  dz, dt);
                            } else {
                                skp = _scalar_concentration(i_top);
                            }
                        }
                    } else {
                        skp = _scalar_concentration(domain_idx);
                    }
                } else if (flux_z_here < 0.0) {
                    skp = _scalar_concentration(domain_idx);
                    if (use_tvd && i_top >= 0 && i_bottom >= 0) {
                        skp = tvd_superbee(_scalar_concentration(i_top),
                                          _scalar_concentration(domain_idx),
                                          _scalar_concentration(i_bottom),
                                          flux_z_here / Az,
                                          dz, dt);
                    }
                }
                
                // k-1/2 face
                if (i_bottom >= 0) {
                    Scalar flux_z_bottom = _flux_z(i_bottom);
                    if (flux_z_bottom > 0.0) {
                        skm = _scalar_concentration(i_bottom);
                        if (use_tvd && i_bottom >= 0 && i_top >= 0) {
                            skm = tvd_superbee(_scalar_concentration(domain_idx),
                                              _scalar_concentration(i_bottom),
                                              _scalar_concentration(i_top),
                                              flux_z_bottom / Az,
                                              dz, dt);
                        }
                    } else if (flux_z_bottom < 0.0) {
                        skm = _scalar_concentration(domain_idx);
                        if (use_tvd && i_bottom >= 0) {
                            // Get k-2 neighbor for TVD
                            Ordinal i_bottom_bottom = -1;
                            if (active_bottom >= 0) {
                                Ordinal active_bottom_bottom = _neighbor_bottom(active_bottom);
                                if (active_bottom_bottom >= 0) {
                                    i_bottom_bottom = _active_to_domain(active_bottom_bottom);
                                }
                            }
                            if (i_bottom_bottom >= 0) {
                                skm = tvd_superbee(_scalar_concentration(i_bottom),
                                                  _scalar_concentration(domain_idx),
                                                  _scalar_concentration(i_bottom_bottom),
                                                  flux_z_bottom / Az,
                                                  dz, dt);
                            } else {
                                skm = _scalar_concentration(domain_idx);
                            }
                        }
                    }
                } else {
                    skm = _scalar_concentration(domain_idx);
                }
                
                // Update scalar mass from advection
                Scalar fx = flux_x_here * sip - (i_left >= 0 ? _flux_x(i_left) : 0.0) * sim;
                Scalar fy = flux_y_here * sjp - (i_back >= 0 ? _flux_y(i_back) : 0.0) * sjm;
                Scalar fz = flux_z_here * skp - (i_bottom >= 0 ? _flux_z(i_bottom) : 0.0) * skm;
                
                _scalar_mass(domain_idx) += dt * (fx + fy + fz);
            });
    }
    
    // ========================================================================
    // COMPUTE DISPERSION TENSOR
    // ========================================================================
    void compute_dispersion_tensor() {
        auto _Dxx = Dxx;
        auto _Dxy = Dxy;
        auto _Dxz = Dxz;
        auto _Dyx = Dyx;
        auto _Dyy = Dyy;
        auto _Dyz = Dyz;
        auto _Dzx = Dzx;
        auto _Dzy = Dzy;
        auto _Dzz = Dzz;
        auto _water_content_sat = state.water_content_sat;
        auto _flux_x = state.flux_x;
        auto _flux_y = state.flux_y;
        auto _flux_z = state.flux_z;
        auto _active_mask_3d = domain.active_mask_3d;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) {
                    _Dxx(domain_idx) = 0.0;
                    _Dxy(domain_idx) = 0.0;
                    _Dxz(domain_idx) = 0.0;
                    _Dyx(domain_idx) = 0.0;
                    _Dyy(domain_idx) = 0.0;
                    _Dyz(domain_idx) = 0.0;
                    _Dzx(domain_idx) = 0.0;
                    _Dzy(domain_idx) = 0.0;
                    _Dzz(domain_idx) = 0.0;
                    return;
                }
                
                // Base diffusion coefficients (molecular diffusion * water content)
                Scalar wcs = _water_content_sat(domain_idx);
                _Dxx(domain_idx) = diff_x * wcs;
                _Dxy(domain_idx) = 0.0;
                _Dxz(domain_idx) = 0.0;
                _Dyx(domain_idx) = 0.0;
                _Dyy(domain_idx) = diff_y * wcs;
                _Dyz(domain_idx) = 0.0;
                _Dzx(domain_idx) = 0.0;
                _Dzy(domain_idx) = 0.0;
                _Dzz(domain_idx) = diff_z * wcs;
                
                // Add velocity-dependent mechanical dispersion
                Scalar qx = _flux_x(domain_idx);
                Scalar qy = _flux_y(domain_idx);
                Scalar qz = _flux_z(domain_idx);
                Scalar q_abs = std::sqrt(qx * qx + qy * qy + qz * qz);
                
                if (q_abs > 0.0 && (disp_longitudinal > 0.0 || disp_transverse > 0.0)) {
                    Scalar qx2 = qx * qx;
                    Scalar qy2 = qy * qy;
                    Scalar qz2 = qz * qz;
                    Scalar qxqy = qx * qy;
                    Scalar qxqz = qx * qz;
                    Scalar qyqz = qy * qz;
                    
                    // Longitudinal dispersion: along flow direction
                    // Transverse dispersion: perpendicular to flow direction
                    _Dxx(domain_idx) += disp_longitudinal * qx2 / q_abs +
                                       disp_transverse * (qy2 + qz2) / q_abs;
                    
                    _Dxy(domain_idx) += (disp_longitudinal - disp_transverse) * qxqy / q_abs;
                    _Dxz(domain_idx) += (disp_longitudinal - disp_transverse) * qxqz / q_abs;
                    
                    _Dyy(domain_idx) += disp_longitudinal * qy2 / q_abs +
                                       disp_transverse * (qx2 + qz2) / q_abs;
                    
                    _Dyx(domain_idx) += (disp_longitudinal - disp_transverse) * qxqy / q_abs;
                    _Dyz(domain_idx) += (disp_longitudinal - disp_transverse) * qyqz / q_abs;
                    
                    _Dzz(domain_idx) += disp_longitudinal * qz2 / q_abs +
                                       disp_transverse * (qx2 + qy2) / q_abs;
                    
                    _Dzx(domain_idx) += (disp_longitudinal - disp_transverse) * qxqz / q_abs;
                    _Dzy(domain_idx) += (disp_longitudinal - disp_transverse) * qyqz / q_abs;
                }
                
                // Ensure non-negative off-diagonal terms (symmetric tensor)
                if (_Dxy(domain_idx) < 0.0) _Dxy(domain_idx) = 0.0;
                if (_Dxz(domain_idx) < 0.0) _Dxz(domain_idx) = 0.0;
                if (_Dyx(domain_idx) < 0.0) _Dyx(domain_idx) = 0.0;
                if (_Dyz(domain_idx) < 0.0) _Dyz(domain_idx) = 0.0;
                if (_Dzx(domain_idx) < 0.0) _Dzx(domain_idx) = 0.0;
                if (_Dzy(domain_idx) < 0.0) _Dzy(domain_idx) = 0.0;
            });
    }
    
    // ========================================================================
    // COMPUTE DISPERSIVE FLUX (With Cross Terms)
    // ========================================================================
    void compute_dispersion() {
        auto _scalar_concentration = scalar_concentration;
        auto _scalar_mass = scalar_mass;
        auto _Dxx = Dxx;
        auto _Dxy = Dxy;
        auto _Dxz = Dxz;
        auto _Dyx = Dyx;
        auto _Dyy = Dyy;
        auto _Dyz = Dyz;
        auto _Dzx = Dzx;
        auto _Dzy = Dzy;
        auto _Dzz = Dzz;
        auto _conductivity_x = state.conductivity_x;
        auto _conductivity_y = state.conductivity_y;
        auto _conductivity_z = state.conductivity_z;
        auto _active_mask_3d = domain.active_mask_3d;
        auto _dz_layers = domain.dz_layers;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _neighbor_bottom = active_mesh.neighbor_bottom;
        auto _neighbor_top = active_mesh.neighbor_top;
        auto _active_to_domain = active_mesh.active_to_domain;
        auto _coord_k = active_mesh.coord_k;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) {
                    return;
                }
                
                // Get neighbors
                Ordinal active_left = _neighbor_left(i);
                Ordinal active_right = _neighbor_right(i);
                Ordinal active_back = _neighbor_back(i);
                Ordinal active_front = _neighbor_front(i);
                Ordinal active_bottom = _neighbor_bottom(i);
                Ordinal active_top = _neighbor_top(i);
                
                Ordinal i_left = (active_left >= 0) ? _active_to_domain(active_left) : -1;
                Ordinal i_right = (active_right >= 0) ? _active_to_domain(active_right) : -1;
                Ordinal i_back = (active_back >= 0) ? _active_to_domain(active_back) : -1;
                Ordinal i_front = (active_front >= 0) ? _active_to_domain(active_front) : -1;
                Ordinal i_bottom = (active_bottom >= 0) ? _active_to_domain(active_bottom) : -1;
                Ordinal i_top = (active_top >= 0) ? _active_to_domain(active_top) : -1;
                
                // Get layer index and thickness
                Ordinal k = _coord_k(i);
                Scalar dz = _dz_layers(k);
                
                // Face areas (simplified - for terrain-following, would need cosx, cosy)
                Scalar Ax = domain.dy * dz;  // X-face area
                Scalar Ay = domain.dx * dz;  // Y-face area
                Scalar Az = domain.dx * domain.dy;  // Z-face area
                
                // Distances (for terrain-following, would use dx/cosx, dy/cosy)
                Scalar dist_x = domain.dx;
                Scalar dist_y = domain.dy;
                Scalar dist_z = dz;
                
                // X-direction dispersive flux at i+1/2 face
                // Includes: fx (Dxx), fy (Dxy), fz (Dxz)
                Scalar jip_fx = 0.0, jip_fy = 0.0, jip_fz = 0.0;
                Scalar jim_fx = 0.0, jim_fy = 0.0, jim_fz = 0.0;
                
                // i+1/2 face: fx term (Dxx - longitudinal)
                if (i_right >= 0 && _conductivity_x(domain_idx) > 0.0 && 
                    _active_mask_3d(i_right) > 0) {
                    Scalar Dxx_face = 0.5 * (_Dxx(domain_idx) + _Dxx(i_right));
                    jip_fx = Dxx_face * Ax * 
                             (_scalar_concentration(i_right) - _scalar_concentration(domain_idx)) / dist_x;
                }
                
                // i+1/2 face: fy term (Dxy - cross term in y-direction)
                // Need 4-point average at j+1/2 and j-1/2 positions
                if (i_front >= 0 && i_right >= 0 && _conductivity_y(domain_idx) > 0.0 &&
                    _active_mask_3d(i_front) > 0 && _active_mask_3d(i_right) > 0) {
                    // Get i_right_front (i+1, j+1, k)
                    Ordinal active_right_front = (active_right >= 0 && active_front >= 0) ? 
                        _neighbor_front(active_right) : -1;
                    Ordinal i_right_front = (active_right_front >= 0) ? 
                        _active_to_domain(active_right_front) : -1;
                    
                    // Get i_right_back (i+1, j-1, k)
                    Ordinal active_right_back = (active_right >= 0 && active_back >= 0) ? 
                        _neighbor_back(active_right) : -1;
                    Ordinal i_right_back = (active_right_back >= 0) ? 
                        _active_to_domain(active_right_back) : -1;
                    
                    // 4-point average at j+1/2: (i,j,k) + (i+1,j,k) + (i,j+1,k) + (i+1,j+1,k)
                    Scalar syp = 0.0;
                    int count_yp = 0;
                    if (i_right >= 0) { syp += _scalar_concentration(i_right); count_yp++; }
                    if (domain_idx >= 0) { syp += _scalar_concentration(domain_idx); count_yp++; }
                    if (i_front >= 0) { syp += _scalar_concentration(i_front); count_yp++; }
                    if (i_right_front >= 0) { syp += _scalar_concentration(i_right_front); count_yp++; }
                    if (count_yp > 0) syp /= count_yp;
                    
                    // 4-point average at j-1/2: (i,j,k) + (i+1,j,k) + (i,j-1,k) + (i+1,j-1,k)
                    Scalar sym = 0.0;
                    int count_ym = 0;
                    if (i_right >= 0) { sym += _scalar_concentration(i_right); count_ym++; }
                    if (domain_idx >= 0) { sym += _scalar_concentration(domain_idx); count_ym++; }
                    if (i_back >= 0) { sym += _scalar_concentration(i_back); count_ym++; }
                    if (i_right_back >= 0) { sym += _scalar_concentration(i_right_back); count_ym++; }
                    if (count_ym > 0) sym /= count_ym;
                    
                    if (count_yp > 0 && count_ym > 0) {
                        jip_fy = _Dxy(domain_idx) * Ay * (syp - sym) / dist_y;
                    }
                }
                
                // i+1/2 face: fz term (Dxz - cross term in z-direction)
                // Need 4-point average at k+1/2 and k-1/2 positions
                if (i_top >= 0 && i_right >= 0 && _conductivity_z(domain_idx) > 0.0 &&
                    _active_mask_3d(i_top) > 0 && _active_mask_3d(i_right) > 0) {
                    // Get i_right_top (i+1, j, k+1)
                    Ordinal active_right_top = (active_right >= 0 && active_top >= 0) ? 
                        _neighbor_top(active_right) : -1;
                    Ordinal i_right_top = (active_right_top >= 0) ? 
                        _active_to_domain(active_right_top) : -1;
                    
                    // Get i_right_bottom (i+1, j, k-1)
                    Ordinal active_right_bottom = (active_right >= 0 && active_bottom >= 0) ? 
                        _neighbor_bottom(active_right) : -1;
                    Ordinal i_right_bottom = (active_right_bottom >= 0) ? 
                        _active_to_domain(active_right_bottom) : -1;
                    
                    // 4-point average at k+1/2: (i,j,k) + (i+1,j,k) + (i,j,k+1) + (i+1,j,k+1)
                    Scalar szp = 0.0;
                    int count_zp = 0;
                    if (i_right >= 0) { szp += _scalar_concentration(i_right); count_zp++; }
                    if (domain_idx >= 0) { szp += _scalar_concentration(domain_idx); count_zp++; }
                    if (i_top >= 0) { szp += _scalar_concentration(i_top); count_zp++; }
                    if (i_right_top >= 0) { szp += _scalar_concentration(i_right_top); count_zp++; }
                    if (count_zp > 0) szp /= count_zp;
                    
                    // 4-point average at k-1/2: (i,j,k) + (i+1,j,k) + (i,j,k-1) + (i+1,j,k-1)
                    Scalar szm = 0.0;
                    int count_zm = 0;
                    if (i_right >= 0) { szm += _scalar_concentration(i_right); count_zm++; }
                    if (domain_idx >= 0) { szm += _scalar_concentration(domain_idx); count_zm++; }
                    if (i_bottom >= 0) { szm += _scalar_concentration(i_bottom); count_zm++; }
                    if (i_right_bottom >= 0) { szm += _scalar_concentration(i_right_bottom); count_zm++; }
                    if (count_zm > 0) szm /= count_zm;
                    
                    if (count_zp > 0 && count_zm > 0) {
                        jip_fz = _Dxz(domain_idx) * Az * (szp - szm) / dist_z;
                    }
                }
                
                // i-1/2 face: similar computation but from left neighbor's perspective
                if (i_left >= 0 && _conductivity_x(i_left) > 0.0 && 
                    _active_mask_3d(i_left) > 0) {
                    // fx term (Dxx)
                    jim_fx = _Dxx(i_left) * Ax * 
                             (_scalar_concentration(domain_idx) - _scalar_concentration(i_left)) / dist_x;
                    
                    // fy term (Dxy) - cross term
                    if (i_front >= 0 && _active_mask_3d(i_front) > 0) {
                        Ordinal active_left_front = (active_left >= 0 && active_front >= 0) ? 
                            _neighbor_front(active_left) : -1;
                        Ordinal i_left_front = (active_left_front >= 0) ? 
                            _active_to_domain(active_left_front) : -1;
                        Ordinal active_left_back = (active_left >= 0 && active_back >= 0) ? 
                            _neighbor_back(active_left) : -1;
                        Ordinal i_left_back = (active_left_back >= 0) ? 
                            _active_to_domain(active_left_back) : -1;
                        
                        Scalar syp = 0.0, sym = 0.0;
                        int count_yp = 0, count_ym = 0;
                        if (domain_idx >= 0) { syp += _scalar_concentration(domain_idx); count_yp++; }
                        if (i_left >= 0) { syp += _scalar_concentration(i_left); count_yp++; }
                        if (i_front >= 0) { syp += _scalar_concentration(i_front); count_yp++; }
                        if (i_left_front >= 0) { syp += _scalar_concentration(i_left_front); count_yp++; }
                        if (count_yp > 0) syp /= count_yp;
                        
                        if (domain_idx >= 0) { sym += _scalar_concentration(domain_idx); count_ym++; }
                        if (i_left >= 0) { sym += _scalar_concentration(i_left); count_ym++; }
                        if (i_back >= 0) { sym += _scalar_concentration(i_back); count_ym++; }
                        if (i_left_back >= 0) { sym += _scalar_concentration(i_left_back); count_ym++; }
                        if (count_ym > 0) sym /= count_ym;
                        
                        if (count_yp > 0 && count_ym > 0) {
                            jim_fy = _Dxy(i_left) * Ay * (syp - sym) / dist_y;
                        }
                    }
                    
                    // fz term (Dxz) - cross term
                    if (i_top >= 0 && _active_mask_3d(i_top) > 0) {
                        Ordinal active_left_top = (active_left >= 0 && active_top >= 0) ? 
                            _neighbor_top(active_left) : -1;
                        Ordinal i_left_top = (active_left_top >= 0) ? 
                            _active_to_domain(active_left_top) : -1;
                        Ordinal active_left_bottom = (active_left >= 0 && active_bottom >= 0) ? 
                            _neighbor_bottom(active_left) : -1;
                        Ordinal i_left_bottom = (active_left_bottom >= 0) ? 
                            _active_to_domain(active_left_bottom) : -1;
                        
                        Scalar szp = 0.0, szm = 0.0;
                        int count_zp = 0, count_zm = 0;
                        if (domain_idx >= 0) { szp += _scalar_concentration(domain_idx); count_zp++; }
                        if (i_left >= 0) { szp += _scalar_concentration(i_left); count_zp++; }
                        if (i_top >= 0) { szp += _scalar_concentration(i_top); count_zp++; }
                        if (i_left_top >= 0) { szp += _scalar_concentration(i_left_top); count_zp++; }
                        if (count_zp > 0) szp /= count_zp;
                        
                        if (domain_idx >= 0) { szm += _scalar_concentration(domain_idx); count_zm++; }
                        if (i_left >= 0) { szm += _scalar_concentration(i_left); count_zm++; }
                        if (i_bottom >= 0) { szm += _scalar_concentration(i_bottom); count_zm++; }
                        if (i_left_bottom >= 0) { szm += _scalar_concentration(i_left_bottom); count_zm++; }
                        if (count_zm > 0) szm /= count_zm;
                        
                        if (count_zp > 0 && count_zm > 0) {
                            jim_fz = _Dxz(i_left) * Az * (szp - szm) / dist_z;
                        }
                    }
                }
                
                Scalar jip = jip_fx + jip_fy + jip_fz;
                Scalar jim = jim_fx + jim_fy + jim_fz;
                
                // Y-direction dispersive flux at j+1/2 face
                // Includes: fx (Dyx), fy (Dyy), fz (Dyz)
                Scalar jjp_fx = 0.0, jjp_fy = 0.0, jjp_fz = 0.0;
                Scalar jjm_fx = 0.0, jjm_fy = 0.0, jjm_fz = 0.0;
                
                // j+1/2 face: fy term (Dyy - longitudinal)
                if (i_front >= 0 && _conductivity_y(domain_idx) > 0.0 && 
                    _active_mask_3d(i_front) > 0) {
                    Scalar Dyy_face = 0.5 * (_Dyy(domain_idx) + _Dyy(i_front));
                    jjp_fy = Dyy_face * Ay * 
                             (_scalar_concentration(i_front) - _scalar_concentration(domain_idx)) / dist_y;
                }
                
                // j+1/2 face: fx term (Dyx - cross term in x-direction)
                if (i_right >= 0 && i_front >= 0 && _conductivity_x(domain_idx) > 0.0 &&
                    _active_mask_3d(i_right) > 0 && _active_mask_3d(i_front) > 0) {
                    Ordinal active_right_front = (active_right >= 0 && active_front >= 0) ? 
                        _neighbor_front(active_right) : -1;
                    Ordinal i_right_front = (active_right_front >= 0) ? 
                        _active_to_domain(active_right_front) : -1;
                    Ordinal active_left_front = (active_left >= 0 && active_front >= 0) ? 
                        _neighbor_front(active_left) : -1;
                    Ordinal i_left_front = (active_left_front >= 0) ? 
                        _active_to_domain(active_left_front) : -1;
                    
                    // 4-point average at i+1/2: (i,j,k) + (i+1,j,k) + (i,j+1,k) + (i+1,j+1,k)
                    Scalar sxp = 0.0;
                    int count_xp = 0;
                    if (domain_idx >= 0) { sxp += _scalar_concentration(domain_idx); count_xp++; }
                    if (i_right >= 0) { sxp += _scalar_concentration(i_right); count_xp++; }
                    if (i_front >= 0) { sxp += _scalar_concentration(i_front); count_xp++; }
                    if (i_right_front >= 0) { sxp += _scalar_concentration(i_right_front); count_xp++; }
                    if (count_xp > 0) sxp /= count_xp;
                    
                    // 4-point average at i-1/2: (i,j,k) + (i-1,j,k) + (i,j+1,k) + (i-1,j+1,k)
                    Scalar sxm = 0.0;
                    int count_xm = 0;
                    if (domain_idx >= 0) { sxm += _scalar_concentration(domain_idx); count_xm++; }
                    if (i_left >= 0) { sxm += _scalar_concentration(i_left); count_xm++; }
                    if (i_front >= 0) { sxm += _scalar_concentration(i_front); count_xm++; }
                    if (i_left_front >= 0) { sxm += _scalar_concentration(i_left_front); count_xm++; }
                    if (count_xm > 0) sxm /= count_xm;
                    
                    if (count_xp > 0 && count_xm > 0) {
                        jjp_fx = _Dyx(domain_idx) * Ax * (sxp - sxm) / dist_x;
                    }
                }
                
                // j+1/2 face: fz term (Dyz - cross term in z-direction)
                if (i_top >= 0 && i_front >= 0 && _conductivity_z(domain_idx) > 0.0 &&
                    _active_mask_3d(i_top) > 0 && _active_mask_3d(i_front) > 0) {
                    Ordinal active_front_top = (active_front >= 0 && active_top >= 0) ? 
                        _neighbor_top(active_front) : -1;
                    Ordinal i_front_top = (active_front_top >= 0) ? 
                        _active_to_domain(active_front_top) : -1;
                    Ordinal active_front_bottom = (active_front >= 0 && active_bottom >= 0) ? 
                        _neighbor_bottom(active_front) : -1;
                    Ordinal i_front_bottom = (active_front_bottom >= 0) ? 
                        _active_to_domain(active_front_bottom) : -1;
                    
                    // 4-point average at k+1/2: (i,j,k) + (i,j+1,k) + (i,j,k+1) + (i,j+1,k+1)
                    Scalar szp = 0.0;
                    int count_zp = 0;
                    if (domain_idx >= 0) { szp += _scalar_concentration(domain_idx); count_zp++; }
                    if (i_front >= 0) { szp += _scalar_concentration(i_front); count_zp++; }
                    if (i_top >= 0) { szp += _scalar_concentration(i_top); count_zp++; }
                    if (i_front_top >= 0) { szp += _scalar_concentration(i_front_top); count_zp++; }
                    if (count_zp > 0) szp /= count_zp;
                    
                    // 4-point average at k-1/2: (i,j,k) + (i,j+1,k) + (i,j,k-1) + (i,j+1,k-1)
                    Scalar szm = 0.0;
                    int count_zm = 0;
                    if (domain_idx >= 0) { szm += _scalar_concentration(domain_idx); count_zm++; }
                    if (i_front >= 0) { szm += _scalar_concentration(i_front); count_zm++; }
                    if (i_bottom >= 0) { szm += _scalar_concentration(i_bottom); count_zm++; }
                    if (i_front_bottom >= 0) { szm += _scalar_concentration(i_front_bottom); count_zm++; }
                    if (count_zm > 0) szm /= count_zm;
                    
                    if (count_zp > 0 && count_zm > 0) {
                        jjp_fz = _Dyz(domain_idx) * Az * (szp - szm) / dist_z;
                    }
                }
                
                // j-1/2 face: similar computation but from back neighbor's perspective
                if (i_back >= 0 && _conductivity_y(i_back) > 0.0 && 
                    _active_mask_3d(i_back) > 0) {
                    // fy term (Dyy)
                    jjm_fy = _Dyy(i_back) * Ay * 
                             (_scalar_concentration(domain_idx) - _scalar_concentration(i_back)) / dist_y;
                    
                    // fx term (Dyx) - cross term
                    if (i_right >= 0 && _active_mask_3d(i_right) > 0) {
                        Ordinal active_right_back = (active_right >= 0 && active_back >= 0) ? 
                            _neighbor_back(active_right) : -1;
                        Ordinal i_right_back = (active_right_back >= 0) ? 
                            _active_to_domain(active_right_back) : -1;
                        Ordinal active_left_back = (active_left >= 0 && active_back >= 0) ? 
                            _neighbor_back(active_left) : -1;
                        Ordinal i_left_back = (active_left_back >= 0) ? 
                            _active_to_domain(active_left_back) : -1;
                        
                        Scalar sxp = 0.0, sxm = 0.0;
                        int count_xp = 0, count_xm = 0;
                        if (domain_idx >= 0) { sxp += _scalar_concentration(domain_idx); count_xp++; }
                        if (i_right >= 0) { sxp += _scalar_concentration(i_right); count_xp++; }
                        if (i_back >= 0) { sxp += _scalar_concentration(i_back); count_xp++; }
                        if (i_right_back >= 0) { sxp += _scalar_concentration(i_right_back); count_xp++; }
                        if (count_xp > 0) sxp /= count_xp;
                        
                        if (domain_idx >= 0) { sxm += _scalar_concentration(domain_idx); count_xm++; }
                        if (i_left >= 0) { sxm += _scalar_concentration(i_left); count_xm++; }
                        if (i_back >= 0) { sxm += _scalar_concentration(i_back); count_xm++; }
                        if (i_left_back >= 0) { sxm += _scalar_concentration(i_left_back); count_xm++; }
                        if (count_xm > 0) sxm /= count_xm;
                        
                        if (count_xp > 0 && count_xm > 0) {
                            jjm_fx = _Dyx(i_back) * Ax * (sxp - sxm) / dist_x;
                        }
                    }
                    
                    // fz term (Dyz) - cross term
                    if (i_top >= 0 && _active_mask_3d(i_top) > 0) {
                        Ordinal active_back_top = (active_back >= 0 && active_top >= 0) ? 
                            _neighbor_top(active_back) : -1;
                        Ordinal i_back_top = (active_back_top >= 0) ? 
                            _active_to_domain(active_back_top) : -1;
                        Ordinal active_back_bottom = (active_back >= 0 && active_bottom >= 0) ? 
                            _neighbor_bottom(active_back) : -1;
                        Ordinal i_back_bottom = (active_back_bottom >= 0) ? 
                            _active_to_domain(active_back_bottom) : -1;
                        
                        Scalar szp = 0.0, szm = 0.0;
                        int count_zp = 0, count_zm = 0;
                        if (domain_idx >= 0) { szp += _scalar_concentration(domain_idx); count_zp++; }
                        if (i_back >= 0) { szp += _scalar_concentration(i_back); count_zp++; }
                        if (i_top >= 0) { szp += _scalar_concentration(i_top); count_zp++; }
                        if (i_back_top >= 0) { szp += _scalar_concentration(i_back_top); count_zp++; }
                        if (count_zp > 0) szp /= count_zp;
                        
                        if (domain_idx >= 0) { szm += _scalar_concentration(domain_idx); count_zm++; }
                        if (i_back >= 0) { szm += _scalar_concentration(i_back); count_zm++; }
                        if (i_bottom >= 0) { szm += _scalar_concentration(i_bottom); count_zm++; }
                        if (i_back_bottom >= 0) { szm += _scalar_concentration(i_back_bottom); count_zm++; }
                        if (count_zm > 0) szm /= count_zm;
                        
                        if (count_zp > 0 && count_zm > 0) {
                            jjm_fz = _Dyz(i_back) * Az * (szp - szm) / dist_z;
                        }
                    }
                }
                
                Scalar jjp = jjp_fx + jjp_fy + jjp_fz;
                Scalar jjm = jjm_fx + jjm_fy + jjm_fz;
                
                // Z-direction dispersive flux at k+1/2 face
                // Includes: fx (Dzx), fy (Dzy), fz (Dzz)
                Scalar jkp_fx = 0.0, jkp_fy = 0.0, jkp_fz = 0.0;
                Scalar jkm_fx = 0.0, jkm_fy = 0.0, jkm_fz = 0.0;
                
                // k+1/2 face: fz term (Dzz - longitudinal)
                if (i_top >= 0 && _conductivity_z(domain_idx) > 0.0 && 
                    _active_mask_3d(i_top) > 0) {
                    Ordinal k_top = _coord_k(active_top);
                    Scalar dz_top = _dz_layers(k_top);
                    Scalar dist_z_top = 0.5 * (dz + dz_top);
                    Scalar Dzz_face = 0.5 * (_Dzz(domain_idx) + _Dzz(i_top));
                    jkp_fz = Dzz_face * Az * 
                             (_scalar_concentration(i_top) - _scalar_concentration(domain_idx)) / dist_z_top;
                }
                
                // k+1/2 face: fx term (Dzx - cross term in x-direction)
                if (i_right >= 0 && i_top >= 0 && _conductivity_x(domain_idx) > 0.0 &&
                    _active_mask_3d(i_right) > 0 && _active_mask_3d(i_top) > 0) {
                    Ordinal active_right_top = (active_right >= 0 && active_top >= 0) ? 
                        _neighbor_top(active_right) : -1;
                    Ordinal i_right_top = (active_right_top >= 0) ? 
                        _active_to_domain(active_right_top) : -1;
                    Ordinal active_left_top = (active_left >= 0 && active_top >= 0) ? 
                        _neighbor_top(active_left) : -1;
                    Ordinal i_left_top = (active_left_top >= 0) ? 
                        _active_to_domain(active_left_top) : -1;
                    
                    // 4-point average at i+1/2: (i,j,k) + (i+1,j,k) + (i,j,k+1) + (i+1,j,k+1)
                    Scalar sxp = 0.0;
                    int count_xp = 0;
                    if (domain_idx >= 0) { sxp += _scalar_concentration(domain_idx); count_xp++; }
                    if (i_right >= 0) { sxp += _scalar_concentration(i_right); count_xp++; }
                    if (i_top >= 0) { sxp += _scalar_concentration(i_top); count_xp++; }
                    if (i_right_top >= 0) { sxp += _scalar_concentration(i_right_top); count_xp++; }
                    if (count_xp > 0) sxp /= count_xp;
                    
                    // 4-point average at i-1/2: (i,j,k) + (i-1,j,k) + (i,j,k+1) + (i-1,j,k+1)
                    Scalar sxm = 0.0;
                    int count_xm = 0;
                    if (domain_idx >= 0) { sxm += _scalar_concentration(domain_idx); count_xm++; }
                    if (i_left >= 0) { sxm += _scalar_concentration(i_left); count_xm++; }
                    if (i_top >= 0) { sxm += _scalar_concentration(i_top); count_xm++; }
                    if (i_left_top >= 0) { sxm += _scalar_concentration(i_left_top); count_xm++; }
                    if (count_xm > 0) sxm /= count_xm;
                    
                    if (count_xp > 0 && count_xm > 0) {
                        jkp_fx = _Dzx(domain_idx) * Ax * (sxp - sxm) / dist_x;
                    }
                }
                
                // k+1/2 face: fy term (Dzy - cross term in y-direction)
                if (i_front >= 0 && i_top >= 0 && _conductivity_y(domain_idx) > 0.0 &&
                    _active_mask_3d(i_front) > 0 && _active_mask_3d(i_top) > 0) {
                    Ordinal active_front_top = (active_front >= 0 && active_top >= 0) ? 
                        _neighbor_top(active_front) : -1;
                    Ordinal i_front_top = (active_front_top >= 0) ? 
                        _active_to_domain(active_front_top) : -1;
                    Ordinal active_back_top = (active_back >= 0 && active_top >= 0) ? 
                        _neighbor_top(active_back) : -1;
                    Ordinal i_back_top = (active_back_top >= 0) ? 
                        _active_to_domain(active_back_top) : -1;
                    
                    // 4-point average at j+1/2: (i,j,k) + (i,j+1,k) + (i,j,k+1) + (i,j+1,k+1)
                    Scalar syp = 0.0;
                    int count_yp = 0;
                    if (domain_idx >= 0) { syp += _scalar_concentration(domain_idx); count_yp++; }
                    if (i_front >= 0) { syp += _scalar_concentration(i_front); count_yp++; }
                    if (i_top >= 0) { syp += _scalar_concentration(i_top); count_yp++; }
                    if (i_front_top >= 0) { syp += _scalar_concentration(i_front_top); count_yp++; }
                    if (count_yp > 0) syp /= count_yp;
                    
                    // 4-point average at j-1/2: (i,j,k) + (i,j-1,k) + (i,j,k+1) + (i,j-1,k+1)
                    Scalar sym = 0.0;
                    int count_ym = 0;
                    if (domain_idx >= 0) { sym += _scalar_concentration(domain_idx); count_ym++; }
                    if (i_back >= 0) { sym += _scalar_concentration(i_back); count_ym++; }
                    if (i_top >= 0) { sym += _scalar_concentration(i_top); count_ym++; }
                    if (i_back_top >= 0) { sym += _scalar_concentration(i_back_top); count_ym++; }
                    if (count_ym > 0) sym /= count_ym;
                    
                    if (count_yp > 0 && count_ym > 0) {
                        jkp_fy = _Dzy(domain_idx) * Ay * (syp - sym) / dist_y;
                    }
                }
                
                // k-1/2 face: similar computation but from bottom neighbor's perspective
                if (i_bottom >= 0 && _conductivity_z(i_bottom) > 0.0 && 
                    _active_mask_3d(i_bottom) > 0) {
                    Ordinal k_bottom = _coord_k(active_bottom);
                    Scalar dz_bottom = _dz_layers(k_bottom);
                    Scalar dist_z_bottom = 0.5 * (dz + dz_bottom);
                    
                    // fz term (Dzz)
                    jkm_fz = _Dzz(i_bottom) * Az * 
                             (_scalar_concentration(domain_idx) - _scalar_concentration(i_bottom)) / dist_z_bottom;
                    
                    // fx term (Dzx) - cross term
                    if (i_right >= 0 && _active_mask_3d(i_right) > 0) {
                        Ordinal active_right_bottom = (active_right >= 0 && active_bottom >= 0) ? 
                            _neighbor_bottom(active_right) : -1;
                        Ordinal i_right_bottom = (active_right_bottom >= 0) ? 
                            _active_to_domain(active_right_bottom) : -1;
                        Ordinal active_left_bottom = (active_left >= 0 && active_bottom >= 0) ? 
                            _neighbor_bottom(active_left) : -1;
                        Ordinal i_left_bottom = (active_left_bottom >= 0) ? 
                            _active_to_domain(active_left_bottom) : -1;
                        
                        Scalar sxp = 0.0, sxm = 0.0;
                        int count_xp = 0, count_xm = 0;
                        if (domain_idx >= 0) { sxp += _scalar_concentration(domain_idx); count_xp++; }
                        if (i_right >= 0) { sxp += _scalar_concentration(i_right); count_xp++; }
                        if (i_bottom >= 0) { sxp += _scalar_concentration(i_bottom); count_xp++; }
                        if (i_right_bottom >= 0) { sxp += _scalar_concentration(i_right_bottom); count_xp++; }
                        if (count_xp > 0) sxp /= count_xp;
                        
                        if (domain_idx >= 0) { sxm += _scalar_concentration(domain_idx); count_xm++; }
                        if (i_left >= 0) { sxm += _scalar_concentration(i_left); count_xm++; }
                        if (i_bottom >= 0) { sxm += _scalar_concentration(i_bottom); count_xm++; }
                        if (i_left_bottom >= 0) { sxm += _scalar_concentration(i_left_bottom); count_xm++; }
                        if (count_xm > 0) sxm /= count_xm;
                        
                        if (count_xp > 0 && count_xm > 0) {
                            jkm_fx = _Dzx(i_bottom) * Ax * (sxp - sxm) / dist_x;
                        }
                    }
                    
                    // fy term (Dzy) - cross term
                    if (i_front >= 0 && _active_mask_3d(i_front) > 0) {
                        Ordinal active_front_bottom = (active_front >= 0 && active_bottom >= 0) ? 
                            _neighbor_bottom(active_front) : -1;
                        Ordinal i_front_bottom = (active_front_bottom >= 0) ? 
                            _active_to_domain(active_front_bottom) : -1;
                        Ordinal active_back_bottom = (active_back >= 0 && active_bottom >= 0) ? 
                            _neighbor_bottom(active_back) : -1;
                        Ordinal i_back_bottom = (active_back_bottom >= 0) ? 
                            _active_to_domain(active_back_bottom) : -1;
                        
                        Scalar syp = 0.0, sym = 0.0;
                        int count_yp = 0, count_ym = 0;
                        if (domain_idx >= 0) { syp += _scalar_concentration(domain_idx); count_yp++; }
                        if (i_front >= 0) { syp += _scalar_concentration(i_front); count_yp++; }
                        if (i_bottom >= 0) { syp += _scalar_concentration(i_bottom); count_yp++; }
                        if (i_front_bottom >= 0) { syp += _scalar_concentration(i_front_bottom); count_yp++; }
                        if (count_yp > 0) syp /= count_yp;
                        
                        if (domain_idx >= 0) { sym += _scalar_concentration(domain_idx); count_ym++; }
                        if (i_back >= 0) { sym += _scalar_concentration(i_back); count_ym++; }
                        if (i_bottom >= 0) { sym += _scalar_concentration(i_bottom); count_ym++; }
                        if (i_back_bottom >= 0) { sym += _scalar_concentration(i_back_bottom); count_ym++; }
                        if (count_ym > 0) sym /= count_ym;
                        
                        if (count_yp > 0 && count_ym > 0) {
                            jkm_fy = _Dzy(i_bottom) * Ay * (syp - sym) / dist_y;
                        }
                    }
                }
                
                Scalar jkp = jkp_fx + jkp_fy + jkp_fz;
                Scalar jkm = jkm_fx + jkm_fy + jkm_fz;
                
                // Update scalar mass from dispersive flux divergence
                // Divergence = (jip - jim) + (jjp - jjm) + (jkp - jkm)
                _scalar_mass(domain_idx) += dt * ((jip - jim) + (jjp - jjm) + (jkp - jkm));
            });
    }
    
    // Update concentration from mass
    void update_concentration() {
        auto _scalar_concentration = scalar_concentration;
        auto _scalar_mass = scalar_mass;
        auto _volume = state.volume;
        auto _active_mask_3d = domain.active_mask_3d;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) > 0) {
                    if (_volume(domain_idx) > 0.0) {
                        _scalar_concentration(domain_idx) = _scalar_mass(domain_idx) / _volume(domain_idx);
                    } else {
                        _scalar_concentration(domain_idx) = 0.0;
                    }
                } else {
                    _scalar_concentration(domain_idx) = 0.0;
                }
            });
    }
};

#endif // FREHG_SCALAR_TRANSPORT_SOLVER_HPP

