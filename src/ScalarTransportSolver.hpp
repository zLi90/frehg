#ifndef FREHG_SCALAR_TRANSPORT_SOLVER_HPP
#define FREHG_SCALAR_TRANSPORT_SOLVER_HPP

#include "define.hpp"
#include "Domain.hpp"
#include "StateVariables.hpp"
#include "ActiveCellMesh.hpp"
#include "BoundaryConditions.hpp"
#include "SourceSinkTerms.hpp"
#include <memory>
#include <cmath>
#include <limits>

// Use classes from Frehg namespace
using Frehg::SwScalarBoundaryConditionManager;
using Frehg::GwScalarBoundaryConditionManager;
using Frehg::ScalarBcType;
using Frehg::SwSourceSinkManager;
using Frehg::SourceSinkType;

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
class GwScalarTransportSolver;  // Forward declaration for cross-reference

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
    
    // --- Boundary Conditions ---
    SwScalarBoundaryConditionManager* bc_manager_;  // Pointer to boundary condition manager
    Ordinal scalar_index_;                          // Index of this scalar species
    
    // --- Source/Sink Terms ---
    SwSourceSinkManager* source_sink_manager_;  // Pointer to source/sink manager
    
    // Constructor
    SwScalarTransportSolver(SwDomain& _domain,
                           ActiveCellMesh& _active_mesh,
                           SwStateVariables& _state,
                           Scalar _dt,
                           Ordinal _scalar_index = 0,
                           Scalar _diff_x = 0.0,
                           Scalar _diff_y = 0.0,
                           Scalar _s_lim_hi = 200.0,
                           Scalar _s_lim_lo = 0.0,
                           bool _use_tvd = true,
                           SwScalarBoundaryConditionManager* _bc_manager = nullptr,
                           SwSourceSinkManager* _source_sink_manager = nullptr)
        : BaseScalarTransportSolver(_dt, _s_lim_hi, _s_lim_lo, _use_tvd),
          domain(_domain), active_mesh(_active_mesh), state(_state),
          diff_x(_diff_x), diff_y(_diff_y),
          bc_manager_(_bc_manager), scalar_index_(_scalar_index),
          source_sink_manager_(_source_sink_manager) {
        
        // Allocate scalar arrays
        scalar_concentration = View1D<Scalar>("scalar_concentration", domain.num_cells_total);
        scalar_mass = View1D<Scalar>("scalar_mass", domain.num_cells_total);
        scalar_old = View1D<Scalar>("scalar_old", domain.num_cells_total);
    }
    
    // Set boundary condition manager
    void set_boundary_condition_manager(SwScalarBoundaryConditionManager* bc_manager) {
        bc_manager_ = bc_manager;
    }
    
    // Set source/sink manager
    void set_source_sink_manager(SwSourceSinkManager* source_sink_manager) {
        source_sink_manager_ = source_sink_manager;
    }
    
    // Main solver function
    void solve(Scalar current_time = 0.0) {
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
        
        // Apply source/sink terms (dilution/concentration effects)
        apply_source_sink_terms(current_time);
        
        // Exchange scalar with subsurface (if coupled)
        // This will be called from ModelDriver after groundwater solve
        
        // Apply limiters
        apply_scalar_limiters();
        
        // Enforce boundary conditions
        enforce_boundary_conditions(current_time);
    }
    
    // ========================================================================
    // ENFORCE BOUNDARY CONDITIONS (2D)
    // ========================================================================
    void enforce_boundary_conditions(Scalar current_time) {
        if (!bc_manager_) return;
        
        auto _scalar = scalar_concentration;
        const auto& bcs = bc_manager_->get_boundary_conditions(scalar_index_);
        
        for (const auto& bc : bcs) {
            if (bc.type != ScalarBcType::PRESCRIBED_CONCENTRATION) continue;
            
            Scalar bc_value = bc.get_value(current_time);
            Ordinal n_bc_cells = bc.cell_indices.size();
            
            if (n_bc_cells == 0) continue;
            
            // Create device view for BC cell indices
            View1D<Ordinal> d_bc_indices("bc_indices", n_bc_cells);
            auto h_bc_indices = Kokkos::create_mirror_view(d_bc_indices);
            for (Ordinal i = 0; i < n_bc_cells; ++i) {
                h_bc_indices(i) = bc.cell_indices[i];
            }
            Kokkos::deep_copy(d_bc_indices, h_bc_indices);
            
            // Apply BC values on device
            Kokkos::parallel_for(RangePolicy(0, n_bc_cells),
                KOKKOS_LAMBDA (const Ordinal i) {
                    Ordinal cell_idx = d_bc_indices(i);
                    _scalar(cell_idx) = bc_value;
                });
        }
    }
    
    // ========================================================================
    // EXCHANGE SCALAR WITH SUBSURFACE
    // ========================================================================
    // Exchanges scalar mass between surface and subsurface during seepage
    // gw_scalar: Groundwater scalar concentration view
    // gw_scalar_mass: Groundwater scalar mass view
    // gw_state: Groundwater state variables
    // gw_domain: Groundwater domain
    void exchange_scalar_with_subsurface(GwScalarTransportSolver* gw_scalar_solver,
                                        const GwStateVariables& gw_state,
                                        const GwDomain& gw_domain);
    
    // Template version that takes arrays directly (called after GwScalarTransportSolver is defined)
    void exchange_scalar_with_subsurface_impl(View1D<Scalar>& gw_scalar_concentration,
                                              View1D<Scalar>& gw_scalar_mass,
                                              const GwStateVariables& gw_state,
                                              const GwDomain& gw_domain) {
        auto _scalar_concentration = scalar_concentration;
        auto _scalar_mass = scalar_mass;
        auto _seepage_rate = state.seepage_rate;
        auto _depth = state.depth;
        auto _area_top = state.area_top;
        
        // Get subsurface scalar concentration
        auto _gw_scalar_concentration = gw_scalar_concentration;
        auto _gw_scalar_mass = gw_scalar_mass;
        
        // Mesh coupling maps
        auto _gw_to_sw = gw_domain.mesh->gw_to_sw_idx;
        auto _sw_to_gw = domain.mesh->sw_to_gw_idx;
        
        // Create host mirrors
        auto h_scalar = Kokkos::create_mirror_view(_scalar_concentration);
        auto h_scalar_mass = Kokkos::create_mirror_view(_scalar_mass);
        auto h_seepage_rate = Kokkos::create_mirror_view(_seepage_rate);
        auto h_depth = Kokkos::create_mirror_view(_depth);
        auto h_area_top = Kokkos::create_mirror_view(_area_top);
        auto h_gw_scalar = Kokkos::create_mirror_view(_gw_scalar_concentration);
        auto h_gw_scalar_mass = Kokkos::create_mirror_view(_gw_scalar_mass);
        auto h_gw_to_sw = Kokkos::create_mirror_view(_gw_to_sw);
        auto h_sw_to_gw = Kokkos::create_mirror_view(_sw_to_gw);
        
        Kokkos::deep_copy(h_scalar, _scalar_concentration);
        Kokkos::deep_copy(h_scalar_mass, _scalar_mass);
        Kokkos::deep_copy(h_seepage_rate, _seepage_rate);
        Kokkos::deep_copy(h_depth, _depth);
        Kokkos::deep_copy(h_area_top, _area_top);
        Kokkos::deep_copy(h_gw_scalar, _gw_scalar_concentration);
        Kokkos::deep_copy(h_gw_scalar_mass, _gw_scalar_mass);
        Kokkos::deep_copy(h_gw_to_sw, _gw_to_sw);
        Kokkos::deep_copy(h_sw_to_gw, _sw_to_gw);
        
        // Get diffusion coefficient (Dzz) - need to get from domain or state
        // For now, use a default value or get from gw_state if available
        Scalar Dzz = 0.0; // Will need to be passed as parameter or stored in state
        
        // Process each surface cell
        for (Ordinal sw_idx = 0; sw_idx < domain.num_cells_total; ++sw_idx) {
            if (domain.active_mask(sw_idx) == 0) continue;
            
            Scalar qss = h_seepage_rate(sw_idx);
            Scalar dept = h_depth(sw_idx);
            Scalar Aini = domain.dx * domain.dy;
            if (h_area_top(sw_idx) > 0.0) {
                Aini = h_area_top(sw_idx);
            }
            
            // Get corresponding subsurface top cell
            Ordinal gw_top_idx = h_sw_to_gw(sw_idx);
            if (gw_top_idx < 0 || gw_top_idx >= gw_domain.num_cells_3d_total) continue;
            
            // Get subsurface scalar concentration at top cell
            Scalar s_surfkP = h_gw_scalar(gw_top_idx); // Subsurface scalar at top
            Scalar s_surf = h_scalar(sw_idx); // Surface scalar
            
            Scalar sseepage = 0.0;
            
            // Seepage (qss > 0): water comes from subsurface to surface
            if (qss > 0.0) {
                // Scalar doesn't leave subsurface if surface is dry
                if (dept > 0.0) {
                    // Advective + diffusive flux
                    // sseepage = qss * s_surfkP + (2*Dzz/dz) * (s_surfkP - s_surf)
                    // For now, simplified: sseepage = qss * s_surfkP
                    // Full implementation would include diffusive term
                    sseepage = qss * s_surfkP;
                    // TODO: Add diffusive term when Dzz is available
                    // Scalar dz = ...; // Get layer thickness
                    // sseepage = qss * s_surfkP + (2.0 * Dzz / dz) * (s_surfkP - s_surf);
                } else {
                    sseepage = 0.0;
                }
            }
            // Infiltration (qss < 0): water goes from surface to subsurface
            else if (qss < 0.0) {
                // Advective + diffusive flux
                sseepage = qss * s_surf;
                // TODO: Add diffusive term when Dzz is available
                // sseepage = qss * s_surf + (2.0 * Dzz / dz) * (s_surfkP - s_surf);
            }
            // No exchange (qss == 0)
            
            // Update surface scalar mass
            h_scalar_mass(sw_idx) += sseepage * Aini * dt;
        }
        
        // Copy back to device
        Kokkos::deep_copy(_scalar_mass, h_scalar_mass);
    }
    
    // ========================================================================
    // APPLY SOURCE/SINK TERMS (DILUTION/CONCENTRATION EFFECTS)
    // ========================================================================
    void apply_source_sink_terms(Scalar current_time) {
        if (!source_sink_manager_) return;
        
        auto _scalar_concentration = scalar_concentration;
        auto _scalar_mass = scalar_mass;
        auto _volume = state.volume;
        auto _volume_old = state.volume_old;
        auto _area_top = state.area_top;
        auto _depth = state.depth;
        
        const auto& ss_terms = source_sink_manager_->get_source_sink_terms();
        Scalar dt_local = dt;
        Scalar default_area = domain.dx * domain.dy;
        
        // Apply each source/sink term using device-side parallel operations
        for (const auto& ss : ss_terms) {
            Scalar ss_value = ss.get_value(current_time);
            Ordinal n_ss_cells = ss.cell_indices.size();
            
            if (n_ss_cells == 0) continue;
            
            // Create device view for source/sink cell indices
            View1D<Ordinal> d_ss_indices("ss_indices", n_ss_cells);
            auto h_ss_indices = Kokkos::create_mirror_view(d_ss_indices);
            for (Ordinal i = 0; i < n_ss_cells; ++i) {
                h_ss_indices(i) = ss.cell_indices[i];
            }
            Kokkos::deep_copy(d_ss_indices, h_ss_indices);
            
            // Capture source/sink properties by value for device
            bool has_scalar_conc = ss.has_scalar_concentration;
            Scalar scalar_conc = ss.scalar_concentration;
            SourceSinkType ss_type = ss.type;
            
            Kokkos::parallel_for(RangePolicy(0, n_ss_cells),
                KOKKOS_LAMBDA (const Ordinal i) {
                    Ordinal cell_idx = d_ss_indices(i);
                    
                    if (_depth(cell_idx) <= 0.0) return;  // Skip dry cells
                    
                    Scalar cell_area = _area_top(cell_idx);
                    if (cell_area <= 0.0) cell_area = default_area;
                    
                    Scalar V_old = _volume_old(cell_idx);
                    Scalar V_new = _volume(cell_idx);
                    
                    if (ss_type == SourceSinkType::VOLUME_FLUX) {
                        // Volume flux: adds/removes water, affects concentration through dilution/concentration
                        Scalar volume_change = ss_value * dt_local;
                        Scalar V_after_ss = V_new + volume_change;
                        
                        if (V_after_ss > 0.0) {
                            if (has_scalar_conc && volume_change > 0.0) {
                                // Source with concentration: add scalar mass
                                Scalar mass_added = scalar_conc * volume_change;
                                _scalar_mass(cell_idx) += mass_added;
                                // Update concentration from new mass and volume
                                _scalar_concentration(cell_idx) = _scalar_mass(cell_idx) / V_after_ss;
                            } else {
                                // Sink or source without concentration: dilution/concentration effect
                                if (V_old > 0.0) {
                                    _scalar_concentration(cell_idx) = _scalar_concentration(cell_idx) * V_old / V_after_ss;
                                }
                            }
                        }
                        
                    } else if (ss_type == SourceSinkType::DEPTH_RATE) {
                        // Depth rate: convert to volume change
                        Scalar depth_change = ss_value * dt_local;
                        Scalar volume_change = depth_change * cell_area;
                        Scalar V_after_ss = V_new + volume_change;
                        
                        if (V_after_ss > 0.0) {
                            if (has_scalar_conc && volume_change > 0.0) {
                                // Source with concentration
                                Scalar mass_added = scalar_conc * volume_change;
                                _scalar_mass(cell_idx) += mass_added;
                                _scalar_concentration(cell_idx) = _scalar_mass(cell_idx) / V_after_ss;
                            } else {
                                // Dilution/concentration
                                if (V_old > 0.0) {
                                    _scalar_concentration(cell_idx) = _scalar_concentration(cell_idx) * V_old / V_after_ss;
                                }
                            }
                        }
                        
                    } else if (ss_type == SourceSinkType::MASS_FLUX) {
                        // Mass flux: directly add scalar mass
                        Scalar mass_added = ss_value * dt_local;
                        _scalar_mass(cell_idx) += mass_added;
                        
                        // Update concentration from new mass
                        if (V_new > 0.0) {
                            _scalar_concentration(cell_idx) = _scalar_mass(cell_idx) / V_new;
                        }
                    }
                });
        }
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
            KOKKOS_LAMBDA (const Ordinal i) {
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
            KOKKOS_LAMBDA (const Ordinal i) {
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
            KOKKOS_LAMBDA (const Ordinal i) {
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
            KOKKOS_LAMBDA (const Ordinal i) {
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
            KOKKOS_LAMBDA (const Ordinal i) {
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
            KOKKOS_LAMBDA (const Ordinal i) {
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
    
    // --- Boundary Conditions ---
    GwScalarBoundaryConditionManager* bc_manager_;  // Pointer to boundary condition manager
    Ordinal scalar_index_;                          // Index of this scalar species
    
    // --- Boundary Concentrations for Advection ---
    // Store boundary concentrations for each face direction
    // bc_conc_xm/xp = boundary concentration for x-/x+ faces (NaN if no BC)
    View1D<Scalar> bc_conc_xm, bc_conc_xp;  // X-direction boundary concentrations
    View1D<Scalar> bc_conc_ym, bc_conc_yp;  // Y-direction boundary concentrations
    View1D<Scalar> bc_conc_zm, bc_conc_zp;  // Z-direction boundary concentrations
    
    // Constructor
    GwScalarTransportSolver(GwDomain& _domain,
                           ActiveCellMesh& _active_mesh,
                           GwStateVariables& _state,
                           Scalar _dt,
                           Ordinal _scalar_index = 0,
                           Scalar _diff_x = 0.0,
                           Scalar _diff_y = 0.0,
                           Scalar _diff_z = 0.0,
                           Scalar _disp_longitudinal = 0.0,
                           Scalar _disp_transverse = 0.0,
                           Scalar _s_lim_hi = 200.0,
                           Scalar _s_lim_lo = 0.0,
                           bool _use_tvd = true,
                           GwScalarBoundaryConditionManager* _bc_manager = nullptr)
        : BaseScalarTransportSolver(_dt, _s_lim_hi, _s_lim_lo, _use_tvd),
          domain(_domain), active_mesh(_active_mesh), state(_state),
          diff_x(_diff_x), diff_y(_diff_y), diff_z(_diff_z),
          disp_longitudinal(_disp_longitudinal), disp_transverse(_disp_transverse),
          bc_manager_(_bc_manager), scalar_index_(_scalar_index) {
        
        // Allocate scalar arrays
        Ordinal num_cells_3d = domain.num_cells_3d_total;
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
        
        // Allocate boundary concentration arrays and initialize to NaN (no BC)
        bc_conc_xm = View1D<Scalar>("bc_conc_xm", num_cells_3d);
        bc_conc_xp = View1D<Scalar>("bc_conc_xp", num_cells_3d);
        bc_conc_ym = View1D<Scalar>("bc_conc_ym", num_cells_3d);
        bc_conc_yp = View1D<Scalar>("bc_conc_yp", num_cells_3d);
        bc_conc_zm = View1D<Scalar>("bc_conc_zm", num_cells_3d);
        bc_conc_zp = View1D<Scalar>("bc_conc_zp", num_cells_3d);
        Kokkos::deep_copy(bc_conc_xm, std::numeric_limits<Scalar>::quiet_NaN());
        Kokkos::deep_copy(bc_conc_xp, std::numeric_limits<Scalar>::quiet_NaN());
        Kokkos::deep_copy(bc_conc_ym, std::numeric_limits<Scalar>::quiet_NaN());
        Kokkos::deep_copy(bc_conc_yp, std::numeric_limits<Scalar>::quiet_NaN());
        Kokkos::deep_copy(bc_conc_zm, std::numeric_limits<Scalar>::quiet_NaN());
        Kokkos::deep_copy(bc_conc_zp, std::numeric_limits<Scalar>::quiet_NaN());
    }
    
    // Set boundary condition manager
    void set_boundary_condition_manager(GwScalarBoundaryConditionManager* bc_manager) {
        bc_manager_ = bc_manager;
    }
    
    // Main solver function
    void solve(Scalar current_time = 0.0) {
        // Save old values
        Kokkos::deep_copy(scalar_old, scalar_concentration);
        
        // Populate boundary concentrations for advection
        populate_boundary_concentrations(current_time);
        
        // Compute scalar mass from concentration
        compute_scalar_mass();
        
        // Advective transport (uses boundary concentrations for inflow)
        compute_advection();
        
        // Compute dispersion tensor
        compute_dispersion_tensor();
        
        // Dispersive transport
        compute_dispersion();
        
        // Update concentration from mass
        update_concentration();
        
        // Exchange scalar with surface (if coupled)
        // This will be called from ModelDriver after surface solve
        
        // Apply limiters (to be implemented later)
        // apply_scalar_limiters();
        
        // Enforce boundary conditions (set prescribed values on boundary cells)
        enforce_boundary_conditions(current_time);
    }
    
    // Populate boundary concentrations for advection
    // This sets bc_conc_* arrays with boundary values for cells on domain boundaries
    void populate_boundary_concentrations(Scalar current_time) {
        if (!bc_manager_) return;
        
        // Reset to NaN (no BC)
        Kokkos::deep_copy(bc_conc_xm, std::numeric_limits<Scalar>::quiet_NaN());
        Kokkos::deep_copy(bc_conc_xp, std::numeric_limits<Scalar>::quiet_NaN());
        Kokkos::deep_copy(bc_conc_ym, std::numeric_limits<Scalar>::quiet_NaN());
        Kokkos::deep_copy(bc_conc_yp, std::numeric_limits<Scalar>::quiet_NaN());
        Kokkos::deep_copy(bc_conc_zm, std::numeric_limits<Scalar>::quiet_NaN());
        Kokkos::deep_copy(bc_conc_zp, std::numeric_limits<Scalar>::quiet_NaN());
        
        const auto& bcs = bc_manager_->get_boundary_conditions(scalar_index_);
        
        // Create host mirrors
        auto h_bc_xm = Kokkos::create_mirror_view(bc_conc_xm);
        auto h_bc_xp = Kokkos::create_mirror_view(bc_conc_xp);
        auto h_bc_ym = Kokkos::create_mirror_view(bc_conc_ym);
        auto h_bc_yp = Kokkos::create_mirror_view(bc_conc_yp);
        auto h_bc_zm = Kokkos::create_mirror_view(bc_conc_zm);
        auto h_bc_zp = Kokkos::create_mirror_view(bc_conc_zp);
        
        Kokkos::deep_copy(h_bc_xm, bc_conc_xm);
        Kokkos::deep_copy(h_bc_xp, bc_conc_xp);
        Kokkos::deep_copy(h_bc_ym, bc_conc_ym);
        Kokkos::deep_copy(h_bc_yp, bc_conc_yp);
        Kokkos::deep_copy(h_bc_zm, bc_conc_zm);
        Kokkos::deep_copy(h_bc_zp, bc_conc_zp);
        
        for (const auto& bc : bcs) {
            if (bc.type != ScalarBcType::PRESCRIBED_CONCENTRATION &&
                bc.type != ScalarBcType::CAUCHY) continue;
            
            Scalar bc_value = bc.get_value(current_time);
            
            // Assign to appropriate face based on boundary face
            for (Ordinal cell_idx : bc.cell_indices) {
                if (bc.face == "x-") {
                    h_bc_xm(cell_idx) = bc_value;
                } else if (bc.face == "x+") {
                    h_bc_xp(cell_idx) = bc_value;
                } else if (bc.face == "y-") {
                    h_bc_ym(cell_idx) = bc_value;
                } else if (bc.face == "y+") {
                    h_bc_yp(cell_idx) = bc_value;
                } else if (bc.face == "z-") {
                    h_bc_zm(cell_idx) = bc_value;
                } else if (bc.face == "z+") {
                    h_bc_zp(cell_idx) = bc_value;
                }
            }
        }
        
        // Copy back to device
        Kokkos::deep_copy(bc_conc_xm, h_bc_xm);
        Kokkos::deep_copy(bc_conc_xp, h_bc_xp);
        Kokkos::deep_copy(bc_conc_ym, h_bc_ym);
        Kokkos::deep_copy(bc_conc_yp, h_bc_yp);
        Kokkos::deep_copy(bc_conc_zm, h_bc_zm);
        Kokkos::deep_copy(bc_conc_zp, h_bc_zp);
    }
    
    // ========================================================================
    // EXCHANGE SCALAR WITH SURFACE
    // ========================================================================
    // Removes scalar mass from subsurface when it seeps to surface
    // sw_scalar_solver: Surface water scalar transport solver
    // sw_state: Surface water state variables (for accessing seepage rate)
    void exchange_scalar_with_surface(SwScalarTransportSolver* sw_scalar_solver,
                                     const SwStateVariables& sw_state,
                                     const SwDomain& sw_domain) {
        if (!sw_scalar_solver) return;
        
        auto _scalar_concentration = scalar_concentration;
        auto _scalar_mass = scalar_mass;
        auto _volume = state.volume;
        
        // Get surface seepage rate and scalar
        auto _sw_seepage_rate = sw_state.seepage_rate;
        auto _sw_scalar_concentration = sw_scalar_solver->scalar_concentration;
        auto _sw_depth = sw_state.depth;
        auto _sw_area_top = sw_state.area_top;
        
        // Mesh coupling maps
        auto _gw_to_sw = domain.mesh->gw_to_sw_idx;
        auto _sw_to_gw = sw_domain.mesh->sw_to_gw_idx;
        
        // Create host mirrors
        auto h_scalar = Kokkos::create_mirror_view(_scalar_concentration);
        auto h_scalar_mass = Kokkos::create_mirror_view(_scalar_mass);
        auto h_volume = Kokkos::create_mirror_view(_volume);
        auto h_sw_seepage_rate = Kokkos::create_mirror_view(_sw_seepage_rate);
        auto h_sw_scalar = Kokkos::create_mirror_view(_sw_scalar_concentration);
        auto h_sw_depth = Kokkos::create_mirror_view(_sw_depth);
        auto h_sw_area_top = Kokkos::create_mirror_view(_sw_area_top);
        auto h_gw_to_sw = Kokkos::create_mirror_view(_gw_to_sw);
        
        Kokkos::deep_copy(h_scalar, _scalar_concentration);
        Kokkos::deep_copy(h_scalar_mass, _scalar_mass);
        Kokkos::deep_copy(h_volume, _volume);
        Kokkos::deep_copy(h_sw_seepage_rate, _sw_seepage_rate);
        Kokkos::deep_copy(h_sw_scalar, _sw_scalar_concentration);
        Kokkos::deep_copy(h_sw_depth, _sw_depth);
        Kokkos::deep_copy(h_sw_area_top, _sw_area_top);
        Kokkos::deep_copy(h_gw_to_sw, _gw_to_sw);
        
        // Get diffusion coefficient (Dzz) - need to get from domain or state
        Scalar Dzz = 0.0; // Will need to be passed as parameter or stored in state
        
        // Get active mesh views
        auto h_active_to_domain = Kokkos::create_mirror_view(active_mesh.active_to_domain);
        auto h_coord_k = Kokkos::create_mirror_view(active_mesh.coord_k);
        auto h_neighbor_top = Kokkos::create_mirror_view(active_mesh.neighbor_top);
        Kokkos::deep_copy(h_active_to_domain, active_mesh.active_to_domain);
        Kokkos::deep_copy(h_coord_k, active_mesh.coord_k);
        Kokkos::deep_copy(h_neighbor_top, active_mesh.neighbor_top);
        
        // Process each active groundwater cell
        for (Ordinal i = 0; i < active_mesh.num_active; ++i) {
            Ordinal domain_idx = h_active_to_domain(i);
            if (domain.active_mask_3d(domain_idx) == 0) continue;
            
            // Check if this is a top cell
            Ordinal k = h_coord_k(i);
            Ordinal active_top = h_neighbor_top(i);
            bool is_top_cell = (active_top < 0) || (k == domain.nz - 1);
            
            if (is_top_cell) {
                // Get corresponding surface cell
                Ordinal sw_idx = h_gw_to_sw(domain_idx);
                if (sw_idx < 0 || sw_idx >= sw_domain.num_cells_total) continue;
                
                Scalar qss = h_sw_seepage_rate(sw_idx);
                Scalar dept = h_sw_depth(sw_idx);
                Scalar Az = domain.dx * domain.dy;
                if (h_sw_area_top(sw_idx) > 0.0) {
                    Az = h_sw_area_top(sw_idx);
                }
                
                // Get scalar concentrations
                Scalar s_subs = h_scalar(domain_idx); // Subsurface scalar
                Scalar s_surf = h_sw_scalar(sw_idx); // Surface scalar
                
                Scalar sseepage = 0.0;
                
                // Seepage (qss > 0): scalar leaves subsurface
                if (qss > 0.0) {
                    // Scalar doesn't leave subsurface if surface is dry
                    if (dept > 0.0) {
                        // Advective + diffusive flux
                        sseepage = qss * s_subs;
                        // TODO: Add diffusive term when Dzz is available
                        // Scalar dz = ...; // Get layer thickness
                        // sseepage = qss * s_subs + (2.0 * Dzz / dz) * (s_subs - s_surf);
                    } else {
                        sseepage = 0.0;
                    }
                }
                // Infiltration (qss < 0): scalar enters subsurface
                else if (qss < 0.0) {
                    // Advective + diffusive flux
                    sseepage = qss * s_surf;
                    // TODO: Add diffusive term when Dzz is available
                    // sseepage = qss * s_surf + (2.0 * Dzz / dz) * (s_surf - s_subs);
                }
                // No exchange (qss == 0)
                
                // Remove scalar mass from subsurface
                h_scalar_mass(domain_idx) -= dt * Az * sseepage;
            }
        }
        
        // Copy back to device
        Kokkos::deep_copy(_scalar_mass, h_scalar_mass);
    }
    
    // ========================================================================
    // ENFORCE BOUNDARY CONDITIONS
    // ========================================================================
    void enforce_boundary_conditions(Scalar current_time) {
        if (!bc_manager_) return;
        
        auto _scalar_concentration = scalar_concentration;
        const auto& bcs = bc_manager_->get_boundary_conditions(scalar_index_);
        
        // Process each boundary condition using device-side parallel operations
        for (const auto& bc : bcs) {
            Scalar bc_value = bc.get_value(current_time);
            Ordinal n_bc_cells = bc.cell_indices.size();
            
            if (n_bc_cells == 0) continue;
            
            // Create device view for BC cell indices
            View1D<Ordinal> d_bc_indices("bc_indices", n_bc_cells);
            auto h_bc_indices = Kokkos::create_mirror_view(d_bc_indices);
            for (Ordinal i = 0; i < n_bc_cells; ++i) {
                h_bc_indices(i) = bc.cell_indices[i];
            }
            Kokkos::deep_copy(d_bc_indices, h_bc_indices);
            
            if (bc.type == ScalarBcType::PRESCRIBED_CONCENTRATION ||
                bc.type == ScalarBcType::CAUCHY) {
                // Dirichlet or Cauchy BC: prescribed concentration
                Kokkos::parallel_for(RangePolicy(0, n_bc_cells),
                    KOKKOS_LAMBDA (const Ordinal i) {
                        Ordinal cell_idx = d_bc_indices(i);
                        _scalar_concentration(cell_idx) = bc_value;
                    });
            }
        }
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
            KOKKOS_LAMBDA (const Ordinal i) {
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
        auto _dz_layers = domain.layer_thickness;
        
        // Boundary concentrations for inflow handling
        auto _bc_xm = bc_conc_xm;
        auto _bc_xp = bc_conc_xp;
        auto _bc_ym = bc_conc_ym;
        auto _bc_yp = bc_conc_yp;
        auto _bc_zm = bc_conc_zm;
        auto _bc_zp = bc_conc_zp;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _neighbor_bottom = active_mesh.neighbor_bottom;
        auto _neighbor_top = active_mesh.neighbor_top;
        auto _active_to_domain = active_mesh.active_to_domain;
        auto _coord_k = active_mesh.coord_k;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            KOKKOS_LAMBDA (const Ordinal i) {
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
                // Flux convention: positive = flow INTO cell from x+ direction (from right)
                Scalar sip = 0.0;  // Scalar at i+1/2 face
                Scalar sim = 0.0;  // Scalar at i-1/2 face
                
                // i+1/2 face
                Scalar flux_x_here = _flux_x(domain_idx);
                if (flux_x_here > 0.0) {
                    // Flow from right to left (inflow from x+ direction)
                    // Use concentration from right neighbor (upwind)
                    if (i_right >= 0) {
                        sip = _scalar_concentration(i_right);
                        if (use_tvd && i_left >= 0) {
                            sip = tvd_superbee(_scalar_concentration(domain_idx),
                                              _scalar_concentration(i_right),
                                              _scalar_concentration(i_left),
                                              flux_x_here / Ax,
                                              domain.dx, dt);
                        }
                    } else {
                        // No right neighbor - check for boundary inflow concentration
                        Scalar bc_val = _bc_xp(domain_idx);
                        if (!Kokkos::isnan(bc_val)) {
                            sip = bc_val;  // Use boundary concentration for inflow
                        } else {
                            sip = _scalar_concentration(domain_idx);  // Zero-gradient fallback
                        }
                    }
                } else if (flux_x_here < 0.0) {
                    // Flow from left to right (outflow to x+ direction)
                    // Use concentration from current cell (upwind)
                    sip = _scalar_concentration(domain_idx);
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
                            sip = tvd_superbee(_scalar_concentration(i_right),
                                              _scalar_concentration(domain_idx),
                                              _scalar_concentration(i_right_right),
                                              flux_x_here / Ax,
                                              domain.dx, dt);
                        } else {
                            sip = _scalar_concentration(domain_idx);
                        }
                    }
                }
                
                // i-1/2 face (flux_x(i_left) = flux at x+ face of left cell)
                // This flux represents flow between left cell and current cell
                // flux_x(i_left) > 0 means flow INTO left cell from current cell (outflow from current)
                // flux_x(i_left) < 0 means flow INTO current cell from left cell (inflow to current)
                if (i_left >= 0) {
                    Scalar flux_x_left = _flux_x(i_left);
                    if (flux_x_left > 0.0) {
                        // Flow from current cell into left cell (outflow from current at i-1/2)
                        sim = _scalar_concentration(domain_idx);
                        if (use_tvd && i_right >= 0) {
                            sim = tvd_superbee(_scalar_concentration(i_left),
                                              _scalar_concentration(domain_idx),
                                              _scalar_concentration(i_right),
                                              flux_x_left / Ax,
                                              domain.dx, dt);
                        }
                    } else if (flux_x_left < 0.0) {
                        // Flow from left cell into current cell (inflow to current at i-1/2)
                        sim = _scalar_concentration(i_left);
                        if (use_tvd) {
                            // Get i-2 neighbor for TVD
                            Ordinal i_left_left = -1;
                            if (active_left >= 0) {
                                Ordinal active_left_left = _neighbor_left(active_left);
                                if (active_left_left >= 0) {
                                    i_left_left = _active_to_domain(active_left_left);
                                }
                            }
                            if (i_left_left >= 0) {
                                sim = tvd_superbee(_scalar_concentration(domain_idx),
                                                  _scalar_concentration(i_left),
                                                  _scalar_concentration(i_left_left),
                                                  flux_x_left / Ax,
                                                  domain.dx, dt);
                            } else {
                                sim = _scalar_concentration(i_left);
                            }
                        }
                    }
                } else {
                    // No left neighbor - this cell is on x- boundary
                    // For x- boundary with prescribed concentration (e.g., freshwater flux BC)
                    Scalar bc_val = _bc_xm(domain_idx);
                    if (!Kokkos::isnan(bc_val)) {
                        sim = bc_val;
                    } else {
                        sim = _scalar_concentration(domain_idx);  // Zero-gradient fallback
                    }
                }
                
                // Y-direction advection
                Scalar sjp = 0.0;  // Scalar at j+1/2 face
                Scalar sjm = 0.0;  // Scalar at j-1/2 face
                
                // j+1/2 face (flux_y > 0 means inflow from y+ direction)
                Scalar flux_y_here = _flux_y(domain_idx);
                if (flux_y_here > 0.0) {
                    // Inflow from front neighbor
                    if (i_front >= 0) {
                        sjp = _scalar_concentration(i_front);
                        if (use_tvd && i_back >= 0) {
                            sjp = tvd_superbee(_scalar_concentration(domain_idx),
                                              _scalar_concentration(i_front),
                                              _scalar_concentration(i_back),
                                              flux_y_here / Ay,
                                              domain.dy, dt);
                        }
                    } else {
                        // No front neighbor - boundary inflow from y+
                        Scalar bc_val = _bc_yp(domain_idx);
                        if (!Kokkos::isnan(bc_val)) {
                            sjp = bc_val;
                        } else {
                            sjp = _scalar_concentration(domain_idx);
                        }
                    }
                } else if (flux_y_here < 0.0) {
                    // Outflow to front neighbor
                    sjp = _scalar_concentration(domain_idx);
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
                            sjp = tvd_superbee(_scalar_concentration(i_front),
                                              _scalar_concentration(domain_idx),
                                              _scalar_concentration(i_front_front),
                                              flux_y_here / Ay,
                                              domain.dy, dt);
                        } else {
                            sjp = _scalar_concentration(domain_idx);
                        }
                    }
                }
                
                // j-1/2 face (flux_y(i_back) = flux at y+ face of back cell)
                // flux_y(i_back) > 0 means inflow to back cell from current (outflow from current)
                // flux_y(i_back) < 0 means inflow to current from back cell
                if (i_back >= 0) {
                    Scalar flux_y_back = _flux_y(i_back);
                    if (flux_y_back > 0.0) {
                        // Outflow from current to back cell
                        sjm = _scalar_concentration(domain_idx);
                        if (use_tvd && i_front >= 0) {
                            sjm = tvd_superbee(_scalar_concentration(i_back),
                                              _scalar_concentration(domain_idx),
                                              _scalar_concentration(i_front),
                                              flux_y_back / Ay,
                                              domain.dy, dt);
                        }
                    } else if (flux_y_back < 0.0) {
                        // Inflow from back cell to current
                        sjm = _scalar_concentration(i_back);
                        if (use_tvd) {
                            // Get j-2 neighbor for TVD
                            Ordinal i_back_back = -1;
                            if (active_back >= 0) {
                                Ordinal active_back_back = _neighbor_back(active_back);
                                if (active_back_back >= 0) {
                                    i_back_back = _active_to_domain(active_back_back);
                                }
                            }
                            if (i_back_back >= 0) {
                                sjm = tvd_superbee(_scalar_concentration(domain_idx),
                                                  _scalar_concentration(i_back),
                                                  _scalar_concentration(i_back_back),
                                                  flux_y_back / Ay,
                                                  domain.dy, dt);
                            } else {
                                sjm = _scalar_concentration(i_back);
                            }
                        }
                    }
                } else {
                    // No back neighbor - boundary on y- side
                    Scalar bc_val = _bc_ym(domain_idx);
                    if (!Kokkos::isnan(bc_val)) {
                        sjm = bc_val;
                    } else {
                        sjm = _scalar_concentration(domain_idx);
                    }
                }
                
                // Z-direction advection
                Scalar skp = 0.0;  // Scalar at k+1/2 face
                Scalar skm = 0.0;  // Scalar at k-1/2 face
                
                // k+1/2 face (flux_z > 0 means inflow from top, i.e., downward flow)
                Scalar flux_z_here = _flux_z(domain_idx);
                if (flux_z_here > 0.0) {
                    // Flow downward (from top into current cell)
                    if (i_top >= 0) {
                        skp = _scalar_concentration(i_top);
                        if (use_tvd && i_bottom >= 0) {
                            skp = tvd_superbee(_scalar_concentration(domain_idx),
                                              _scalar_concentration(i_top),
                                              _scalar_concentration(i_bottom),
                                              flux_z_here / Az,
                                              dz, dt);
                        }
                    } else {
                        // No top neighbor - boundary on z+ side
                        Scalar bc_val = _bc_zp(domain_idx);
                        if (!Kokkos::isnan(bc_val)) {
                            skp = bc_val;
                        } else {
                            skp = _scalar_concentration(domain_idx);
                        }
                    }
                } else if (flux_z_here < 0.0) {
                    // Flow upward (from current cell to top)
                    skp = _scalar_concentration(domain_idx);
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
                            skp = tvd_superbee(_scalar_concentration(i_top),
                                              _scalar_concentration(domain_idx),
                                              _scalar_concentration(i_top_top),
                                              flux_z_here / Az,
                                              dz, dt);
                        } else {
                            skp = _scalar_concentration(domain_idx);
                        }
                    }
                }
                
                // k-1/2 face (flux_z(i_bottom) = flux at z+ face of bottom cell)
                // flux_z(i_bottom) > 0 means inflow to bottom cell from current (downward from current)
                // flux_z(i_bottom) < 0 means outflow from bottom cell to current (upward from bottom)
                if (i_bottom >= 0) {
                    Scalar flux_z_bottom = _flux_z(i_bottom);
                    if (flux_z_bottom > 0.0) {
                        // Downward flow from current to bottom
                        skm = _scalar_concentration(domain_idx);
                        if (use_tvd && i_top >= 0) {
                            skm = tvd_superbee(_scalar_concentration(i_bottom),
                                              _scalar_concentration(domain_idx),
                                              _scalar_concentration(i_top),
                                              flux_z_bottom / Az,
                                              dz, dt);
                        }
                    } else if (flux_z_bottom < 0.0) {
                        // Upward flow from bottom to current
                        skm = _scalar_concentration(i_bottom);
                        if (use_tvd) {
                            // Get k-2 neighbor for TVD
                            Ordinal i_bottom_bottom = -1;
                            if (active_bottom >= 0) {
                                Ordinal active_bottom_bottom = _neighbor_bottom(active_bottom);
                                if (active_bottom_bottom >= 0) {
                                    i_bottom_bottom = _active_to_domain(active_bottom_bottom);
                                }
                            }
                            if (i_bottom_bottom >= 0) {
                                skm = tvd_superbee(_scalar_concentration(domain_idx),
                                                  _scalar_concentration(i_bottom),
                                                  _scalar_concentration(i_bottom_bottom),
                                                  flux_z_bottom / Az,
                                                  dz, dt);
                            } else {
                                skm = _scalar_concentration(i_bottom);
                            }
                        }
                    }
                } else {
                    // No bottom neighbor - boundary on z- side
                    Scalar bc_val = _bc_zm(domain_idx);
                    if (!Kokkos::isnan(bc_val)) {
                        skm = bc_val;
                    } else {
                        skm = _scalar_concentration(domain_idx);
                    }
                }
                
                // Update scalar mass from advection (legacy convention)
                // flux_x(i) > 0 means flow INTO cell i from x+ direction
                // Mass change = flux_x(i)*sip - flux_x(i-1)*sim
                //   - flux_x(i)*sip: if positive (inflow from right), adds mass with concentration sip
                //   - flux_x(i-1)*sim: outflow to left (=inflow to left cell), removes mass with concentration sim
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
            KOKKOS_LAMBDA (const Ordinal i) {
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
        auto _dz_layers = domain.layer_thickness;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _neighbor_bottom = active_mesh.neighbor_bottom;
        auto _neighbor_top = active_mesh.neighbor_top;
        auto _active_to_domain = active_mesh.active_to_domain;
        auto _coord_k = active_mesh.coord_k;
        
        // Boundary concentrations for dispersive flux at boundaries
        auto _bc_xm = bc_conc_xm;
        auto _bc_xp = bc_conc_xp;
        auto _bc_ym = bc_conc_ym;
        auto _bc_yp = bc_conc_yp;
        auto _bc_zm = bc_conc_zm;
        auto _bc_zp = bc_conc_zp;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            KOKKOS_LAMBDA (const Ordinal i) {
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
                if (_conductivity_x(domain_idx) > 0.0) {
                    if (i_right >= 0 && _active_mask_3d(i_right) > 0) {
                        // Interior face - flux with neighbor
                        Scalar Dxx_face = 0.5 * (_Dxx(domain_idx) + _Dxx(i_right));
                        jip_fx = Dxx_face * Ax * 
                                 (_scalar_concentration(i_right) - _scalar_concentration(domain_idx)) / dist_x;
                    } else {
                        // Boundary face - use boundary concentration if available
                        Scalar bc_val = _bc_xp(domain_idx);
                        if (!Kokkos::isnan(bc_val)) {
                            // Dispersive flux from boundary
                            jip_fx = _Dxx(domain_idx) * Ax * 
                                     (bc_val - _scalar_concentration(domain_idx)) / (0.5 * dist_x);
                        }
                    }
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
                if (_conductivity_x(domain_idx) > 0.0) {
                    if (i_left >= 0 && _active_mask_3d(i_left) > 0) {
                        // Interior face - flux from left
                        jim_fx = _Dxx(i_left) * Ax * 
                                 (_scalar_concentration(domain_idx) - _scalar_concentration(i_left)) / dist_x;
                    } else {
                        // Boundary face - use boundary concentration if available
                        Scalar bc_val = _bc_xm(domain_idx);
                        if (!Kokkos::isnan(bc_val)) {
                            // Dispersive flux from left boundary
                            jim_fx = _Dxx(domain_idx) * Ax * 
                                     (_scalar_concentration(domain_idx) - bc_val) / (0.5 * dist_x);
                        }
                    }
                }
                // Cross-term fy for left boundary (inside the i_left check originally)
                if (i_left >= 0 && _conductivity_x(i_left) > 0.0 && _active_mask_3d(i_left) > 0) {
                    
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
                
                // Update scalar mass from dispersive flux
                // Flux convention: j = D * A * (S_neighbor - S_here) / dx
                //   - If S_neighbor > S_here: j > 0, flux INTO cell from neighbor direction
                // Mass balance: dM/dt = sum of fluxes into cell
                //   - jip is flux at i+1/2 face: positive means flux INTO cell from right
                //   - jim is flux at i-1/2 face (from left cell's perspective): positive means flux OUT of left cell
                // Correct form: mass_change = (flux_in_from_left - flux_out_to_right) = jim - jip
                // But actually, let's reconsider:
                //   - jip = D * (S_right - S_here) / dx at i+1/2 face
                //   - jim = D * (S_here - S_left) / dx at i-1/2 face (computed from left cell)
                // The diffusion equation: dS/dt = D * d²S/dx² = D * (S_{i+1} - 2*S_i + S_{i-1}) / dx²
                // In FV: dM/dt = D*A/dx * (S_{i+1} - S_i) - D*A/dx * (S_i - S_{i-1})
                //              = D*A/dx * (S_{i+1} - 2*S_i + S_{i-1})
                // So: mass += dt * (jip - jim) where jip = D*A*(S_{i+1} - S_i)/dx, jim = D*A*(S_i - S_{i-1})/dx
                // This is correct if jim is computed as (S_here - S_left), but looking at code...
                // jim_fx = D * (S_here - S_left) / dx, so jip - jim = D*A/dx*(S_right - 2*S_here + S_left) ✓
                _scalar_mass(domain_idx) += dt * ((jip - jim) + (jjp - jjm) + (jkp - jkm));
            });
    }
    
    // Update concentration from mass with scalar limiter
    // Prevents numerical overshoots by clamping to neighbor min/max values
    void update_concentration() {
        auto _scalar_concentration = scalar_concentration;
        auto _scalar_mass = scalar_mass;
        auto _scalar_old = scalar_old;  // Old concentration for limiter bounds
        auto _volume = state.volume;
        auto _active_mask_3d = domain.active_mask_3d;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _neighbor_bottom = active_mesh.neighbor_bottom;
        auto _neighbor_top = active_mesh.neighbor_top;
        
        // Boundary concentrations for limiter bounds
        auto _bc_xm = bc_conc_xm;
        auto _bc_xp = bc_conc_xp;
        auto _bc_ym = bc_conc_ym;
        auto _bc_yp = bc_conc_yp;
        auto _bc_zm = bc_conc_zm;
        auto _bc_zp = bc_conc_zp;
        
        // Scalar limits for extreme value detection
        constexpr Scalar s_lim_hi = 200.0;  // Max reasonable concentration
        constexpr Scalar s_lim_lo = 0.0;    // Min reasonable concentration
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            KOKKOS_LAMBDA (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) <= 0) {
                    _scalar_concentration(domain_idx) = 0.0;
                    return;
                }
                
                // Compute new concentration from mass
                Scalar s_new = 0.0;
                if (_volume(domain_idx) > 0.0) {
                    s_new = _scalar_mass(domain_idx) / _volume(domain_idx);
                }
                
                // Compute min/max from neighbors and boundary values (using OLD concentrations)
                // Start with current cell's old concentration
                Scalar s_min = _scalar_old(domain_idx);
                Scalar s_max = _scalar_old(domain_idx);
                
                // Check all neighbors
                Ordinal active_left = _neighbor_left(i);
                Ordinal active_right = _neighbor_right(i);
                Ordinal active_back = _neighbor_back(i);
                Ordinal active_front = _neighbor_front(i);
                Ordinal active_bottom = _neighbor_bottom(i);
                Ordinal active_top = _neighbor_top(i);
                
                // X- direction
                if (active_left >= 0) {
                    Ordinal idx_left = _active_to_domain(active_left);
                    if (_active_mask_3d(idx_left) > 0) {
                        if (_scalar_old(idx_left) > s_max) s_max = _scalar_old(idx_left);
                        if (_scalar_old(idx_left) < s_min) s_min = _scalar_old(idx_left);
                    }
                } else {
                    // Check boundary
                    Scalar bc_val = _bc_xm(domain_idx);
                    if (!Kokkos::isnan(bc_val)) {
                        if (bc_val > s_max) s_max = bc_val;
                        if (bc_val < s_min) s_min = bc_val;
                    }
                }
                
                // X+ direction
                if (active_right >= 0) {
                    Ordinal idx_right = _active_to_domain(active_right);
                    if (_active_mask_3d(idx_right) > 0) {
                        if (_scalar_old(idx_right) > s_max) s_max = _scalar_old(idx_right);
                        if (_scalar_old(idx_right) < s_min) s_min = _scalar_old(idx_right);
                    }
                } else {
                    // Check boundary
                    Scalar bc_val = _bc_xp(domain_idx);
                    if (!Kokkos::isnan(bc_val)) {
                        if (bc_val > s_max) s_max = bc_val;
                        if (bc_val < s_min) s_min = bc_val;
                    }
                }
                
                // Y- direction
                if (active_back >= 0) {
                    Ordinal idx_back = _active_to_domain(active_back);
                    if (_active_mask_3d(idx_back) > 0) {
                        if (_scalar_old(idx_back) > s_max) s_max = _scalar_old(idx_back);
                        if (_scalar_old(idx_back) < s_min) s_min = _scalar_old(idx_back);
                    }
                } else {
                    Scalar bc_val = _bc_ym(domain_idx);
                    if (!Kokkos::isnan(bc_val)) {
                        if (bc_val > s_max) s_max = bc_val;
                        if (bc_val < s_min) s_min = bc_val;
                    }
                }
                
                // Y+ direction
                if (active_front >= 0) {
                    Ordinal idx_front = _active_to_domain(active_front);
                    if (_active_mask_3d(idx_front) > 0) {
                        if (_scalar_old(idx_front) > s_max) s_max = _scalar_old(idx_front);
                        if (_scalar_old(idx_front) < s_min) s_min = _scalar_old(idx_front);
                    }
                } else {
                    Scalar bc_val = _bc_yp(domain_idx);
                    if (!Kokkos::isnan(bc_val)) {
                        if (bc_val > s_max) s_max = bc_val;
                        if (bc_val < s_min) s_min = bc_val;
                    }
                }
                
                // Z- direction
                if (active_bottom >= 0) {
                    Ordinal idx_bottom = _active_to_domain(active_bottom);
                    if (_active_mask_3d(idx_bottom) > 0) {
                        if (_scalar_old(idx_bottom) > s_max) s_max = _scalar_old(idx_bottom);
                        if (_scalar_old(idx_bottom) < s_min) s_min = _scalar_old(idx_bottom);
                    }
                } else {
                    Scalar bc_val = _bc_zm(domain_idx);
                    if (!Kokkos::isnan(bc_val)) {
                        if (bc_val > s_max) s_max = bc_val;
                        if (bc_val < s_min) s_min = bc_val;
                    }
                }
                
                // Z+ direction
                if (active_top >= 0) {
                    Ordinal idx_top = _active_to_domain(active_top);
                    if (_active_mask_3d(idx_top) > 0) {
                        if (_scalar_old(idx_top) > s_max) s_max = _scalar_old(idx_top);
                        if (_scalar_old(idx_top) < s_min) s_min = _scalar_old(idx_top);
                    }
                } else {
                    Scalar bc_val = _bc_zp(domain_idx);
                    if (!Kokkos::isnan(bc_val)) {
                        if (bc_val > s_max) s_max = bc_val;
                        if (bc_val < s_min) s_min = bc_val;
                    }
                }
                
                // Apply scalar limiter: clamp to [s_min, s_max]
                // Only apply if bounds are reasonable (not at global limits)
                if (s_new > s_max && s_max < s_lim_hi) {
                    s_new = s_max;
                } else if (s_new < s_min && s_min > s_lim_lo) {
                    s_new = s_min;
                }
                
                // Apply hard limits for extreme values
                if (s_new > s_lim_hi) {
                    s_new = s_lim_hi;
                } else if (s_new < 0.0) {
                    s_new = 0.0;
                }
                
                _scalar_concentration(domain_idx) = s_new;
            });
    }
};

// ============================================================================
//                      DEFERRED IMPLEMENTATION
// ============================================================================
// Implement SwScalarTransportSolver::exchange_scalar_with_subsurface after
// GwScalarTransportSolver is fully defined

inline void SwScalarTransportSolver::exchange_scalar_with_subsurface(
    GwScalarTransportSolver* gw_scalar_solver,
    const GwStateVariables& gw_state,
    const GwDomain& gw_domain) {
    if (!gw_scalar_solver) return;
    
    exchange_scalar_with_subsurface_impl(
        gw_scalar_solver->scalar_concentration,
        gw_scalar_solver->scalar_mass,
        gw_state, gw_domain);
}

#endif // FREHG_SCALAR_TRANSPORT_SOLVER_HPP

