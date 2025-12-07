#ifndef FREHG_GROUNDWATER_SOLVER_HPP
#define FREHG_GROUNDWATER_SOLVER_HPP

#include "define.hpp"
#include "Domain.hpp"
#include "StateVariables.hpp"
#include "ActiveCellMesh.hpp"
#include "PCGSolver.hpp"
#include "BoundaryConditions.hpp"
#include <memory>
#include <cmath>
#include <Kokkos_Core.hpp>

// Use classes from Frehg namespace
using Frehg::GwBoundaryConditionManager;
using Frehg::GwBcType;

// Forward declarations
class SwStateVariables;
class SwDomain;

// ============================================================================
//                      GROUNDWATER EQUATION SOLVER
// ============================================================================
// Solves the 3D Richards equation for variably saturated groundwater flow
// Uses Kokkos for parallelization and PCG for linear system solution
// Implements predictor-corrector scheme (no Newton iteration)

class GroundwaterSolver {
public:
    // --- Domain and Mesh ---
    GwDomain& domain;
    ActiveCellMesh& active_mesh;
    
    // --- State Variables ---
    GwStateVariables& state;
    
    // --- Physical Parameters ---
    Scalar dt;                  // Time step
    Scalar Ss;                  // Specific storage
    Scalar min_depth;           // Minimum depth threshold
    bool use_corrector;         // Use corrector step
    bool follow_terrain;        // Follow terrain (slope effects)
    bool baroclinic;            // Baroclinic effects (density/viscosity)
    
    // --- Linear Solver ---
    std::unique_ptr<PCGSolver3D> pcg_solver;
    PreconditionerType precond_type;
    Scalar solver_tolerance;
    Ordinal max_solver_iterations;
    
    // --- Solution Vector (active cells only) ---
    View1D<Scalar> solution;    // Solution vector for PCG (h_new in active space)
    
    // --- Boundary Conditions ---
    GwBoundaryConditionManager* bc_manager_;  // Pointer to boundary condition manager
    
    // --- Surface-Subsurface Coupling ---
    bool coupled_with_sw_;          // Is coupled with shallow water
    SwStateVariables* sw_state_;    // Pointer to surface water state (optional)
    const SwDomain* sw_domain_;     // Pointer to surface water domain (optional)
    
    // Constructor
    GroundwaterSolver(GwDomain& _domain,
                      ActiveCellMesh& _active_mesh,
                      GwStateVariables& _state,
                      Scalar _dt,
                      Scalar _Ss = 1.0e-5,
                      Scalar _min_depth = 1.0e-6,
                      bool _use_corrector = true,
                      bool _follow_terrain = false,
                      bool _baroclinic = false,
                      PreconditionerType _precond_type = PreconditionerType::JACOBI,
                      Scalar _solver_tolerance = 1e-8,
                      Ordinal _max_solver_iterations = 10000,
                      GwBoundaryConditionManager* _bc_manager = nullptr)
        : domain(_domain), active_mesh(_active_mesh), state(_state),
          dt(_dt), Ss(_Ss), min_depth(_min_depth),
          use_corrector(_use_corrector), follow_terrain(_follow_terrain),
          baroclinic(_baroclinic),
          precond_type(_precond_type), solver_tolerance(_solver_tolerance),
          max_solver_iterations(_max_solver_iterations),
          bc_manager_(_bc_manager),
          coupled_with_sw_(false), sw_state_(nullptr), sw_domain_(nullptr) {
        
        // Initialize PCG solver
        pcg_solver = std::make_unique<PCGSolver3D>(
            active_mesh.num_active, precond_type, solver_tolerance, max_solver_iterations);
        
        // Allocate solution vector for active cells
        solution = View1D<Scalar>("solution", active_mesh.num_active);
    }
    
    // Set boundary condition manager
    void set_boundary_condition_manager(GwBoundaryConditionManager* bc_manager) {
        bc_manager_ = bc_manager;
    }
    
    // Set surface water coupling
    void set_surface_coupling(SwStateVariables* sw_state, const SwDomain* sw_domain) {
        coupled_with_sw_ = (sw_state != nullptr && sw_domain != nullptr);
        sw_state_ = sw_state;
        sw_domain_ = sw_domain;
    }
    
    // Set time step (for adaptive time stepping)
    void set_time_step(Scalar new_dt) {
        dt = new_dt;
    }
    
    // Get current time step
    Scalar get_time_step() const {
        return dt;
    }
    
    // ========================================================================
    // UPDATE DENSITY/VISCOSITY FROM SCALAR CONCENTRATION
    // ========================================================================
    // Updates density and viscosity ratios based on scalar (e.g., salinity) concentration
    // For seawater: r_rho = 1.0 + S * 0.000744 (where S is salinity in g/kg)
    //               r_visc = 1.0 / (1.0 + S * 0.0022)
    // This should be called before solve() when baroclinic effects are enabled
    void update_density_viscosity_from_scalar(const View1D<Scalar>& scalar_concentration,
                                              Scalar density_coeff = 0.000744,
                                              Scalar viscosity_coeff = 0.0022) {
        if (!baroclinic) return;
        
        auto _density_ratio = state.density_ratio;
        auto _viscosity_ratio = state.viscosity_ratio;
        auto _active_mask_3d = domain.active_mask_3d;
        auto _active_to_domain = active_mesh.active_to_domain;
        auto _scalar = scalar_concentration;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            KOKKOS_LAMBDA (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) {
                    _density_ratio(domain_idx) = 1.0;
                    _viscosity_ratio(domain_idx) = 1.0;
                    return;
                }
                
                Scalar s = _scalar(domain_idx);
                
                // Density ratio: rho/rho_0 = 1 + coeff * S
                // Typical for seawater: coeff ≈ 0.000744 per g/kg salinity
                _density_ratio(domain_idx) = 1.0 + s * density_coeff;
                
                // Viscosity ratio: mu_0/mu = 1 / (1 + coeff * S)
                // Typical for seawater: coeff ≈ 0.0022 per g/kg salinity
                _viscosity_ratio(domain_idx) = 1.0 / (1.0 + s * viscosity_coeff);
            });
    }

    // Main solver function (predictor-corrector scheme)
    void solve(Scalar current_time = 0.0) {
        // Save old values (copies pressure->pressure_old, pressure_old->head_prev_prev, etc.)
        state.update_old_values();
        
        // Compute water content from head and specific moisture capacity
        compute_water_content_from_head();
        compute_specific_moisture_capacity();
        
        // >>> PREDICTOR STEP <<<
        // 1. Compute hydraulic conductivity on faces
        compute_K_face();
        
        // 2. Compute baroclinic effects (density/viscosity ratios)
        if (baroclinic) {
            baroclinic_face();
        }
        
        // 3. Assemble linear system
        assemble_linear_system();
        
        // 4. Add surface-subsurface coupling contribution to RHS
        if (coupled_with_sw_) {
            add_surface_coupling_to_rhs();
        }
        
        // 5. Apply boundary conditions to matrix
        apply_boundary_conditions_to_matrix(current_time);
        
        // 6. Solve for new head using PCG
        solve_linear_system();
        
        // 7. Enforce boundary conditions on solution
        enforce_boundary_conditions_on_solution(current_time);
        
        // >>> CORRECTOR STEP <<<
        if (use_corrector) {
            // 1. Recompute hydraulic conductivity
            compute_K_face();
            
            // 2. Compute groundwater fluxes
            compute_groundwater_flux();
            
            // 3. Check available room (pore space)
            check_room();
            
            // 4. Update water content from fluxes
            update_water_content();
        } else {
            // Direct update from head
            compute_water_content_from_head();
        }
        
        // Final update of derived quantities
        compute_head_from_water_content();
        compute_water_content_from_head();
    }

private:
    // ========================================================================
    // HELPER FUNCTIONS (Mualem-van Genuchten Soil Property Functions)
    // ========================================================================
    // Implements the Mualem-van Genuchten model for variably saturated flow
    // References:
    // - van Genuchten (1980): A closed-form equation for predicting the hydraulic
    //   conductivity of unsaturated soils
    // - Mualem (1976): A new model for predicting the hydraulic conductivity
    
    // Compute effective saturation from head (Se = (θ - θr)/(θs - θr))
    // Modified van Genuchten model with air entry value (ha)
    // If ha = 0, returns to original van Genuchten model
    KOKKOS_INLINE_FUNCTION
    Scalar compute_Se_from_h(Scalar h, Ordinal idx) const {
        Scalar ha = state.vg_ha(idx);
        
        // Modified van Genuchten: if h >= -ha, then Se = 1.0 (saturated)
        // This improves convergence by maintaining saturation until air entry pressure
        if (h >= -ha) {
            return 1.0; // Saturated until air entry value
        }
        
        Scalar alpha = state.vg_alpha(idx);
        Scalar n = state.vg_n(idx);
        Scalar m = state.vg_m(idx);
        
        // Modified van Genuchten water retention curve
        // Se = [1 + (α(|h| - ha))^n]^(-m) for h < -ha
        // If ha = 0, this reduces to original: Se = [1 + (α|h|)^n]^(-m)
        Scalar abs_h = std::abs(h);
        Scalar h_effective = abs_h - ha;  // Effective head below air entry
        
        // Ensure h_effective is non-negative
        if (h_effective < 0.0) {
            return 1.0; // Shouldn't happen if h >= -ha check above, but safety check
        }
        
        Scalar term = 1.0 + std::pow(alpha * h_effective, n);
        Scalar Se = std::pow(term, -m);
        
        // Clamp to [0, 1]
        if (Se > 1.0) Se = 1.0;
        if (Se < 0.0) Se = 0.0;
        
        return Se;
    }
    
    // Compute water content from head using modified van Genuchten model
    KOKKOS_INLINE_FUNCTION
    Scalar compute_wch(Scalar h, Ordinal idx) const {
        Scalar wcs = state.water_content_sat(idx);
        Scalar wcr = state.water_content_res(idx);
        Scalar ha = state.vg_ha(idx);
        
        // Modified van Genuchten: if h >= -ha, then saturated
        if (h >= -ha) {
            return wcs; // Saturated until air entry value
        }
        
        Scalar Se = compute_Se_from_h(h, idx);
        return wcr + (wcs - wcr) * Se;
    }
    
    // Compute head from water content (inverse modified van Genuchten)
    KOKKOS_INLINE_FUNCTION
    Scalar compute_hwc(Ordinal idx) const {
        Scalar wc = state.water_content(idx);
        Scalar wcs = state.water_content_sat(idx);
        Scalar wcr = state.water_content_res(idx);
        Scalar ha = state.vg_ha(idx);
        
        if (wc >= wcs) {
            return -ha; // At air entry value when saturated
        }
        if (wc <= wcr) {
            return -1.0e6; // Very dry (large negative head)
        }
        
        // Compute effective saturation
        Scalar Se = (wc - wcr) / (wcs - wcr);
        if (Se >= 1.0) return -ha; // At air entry value
        if (Se <= 0.0) return -1.0e6;
        
        Scalar alpha = state.vg_alpha(idx);
        Scalar n = state.vg_n(idx);
        Scalar m = state.vg_m(idx);
        
        // Inverse modified van Genuchten: 
        // h = -ha - (1/α) * [(Se^(-1/m) - 1)^(1/n)]
        // If ha = 0, this reduces to original: h = -(1/α) * [(Se^(-1/m) - 1)^(1/n)]
        Scalar Se_inv_m = std::pow(Se, -1.0 / m);
        if (Se_inv_m <= 1.0) {
            return -ha; // Shouldn't happen, but handle edge case
        }
        Scalar term = std::pow(Se_inv_m - 1.0, 1.0 / n);
        Scalar h_effective = term / alpha;  // Effective head below air entry
        Scalar h = -ha - h_effective;
        
        return h;
    }
    
    // Compute relative permeability using Mualem model
    // Kr(Se) = Se^0.5 * [1 - (1 - Se^(1/m))^m]^2
    KOKKOS_INLINE_FUNCTION
    Scalar compute_Kr_Mualem(Scalar Se, Ordinal idx) const {
        if (Se >= 1.0) {
            return 1.0; // Fully saturated
        }
        if (Se <= 0.0) {
            return 0.0; // Completely dry
        }
        
        Scalar m = state.vg_m(idx);
        
        // Mualem relative permeability model
        Scalar Se_sqrt = std::sqrt(Se);
        Scalar Se_1_m = std::pow(Se, 1.0 / m);
        Scalar term = 1.0 - std::pow(1.0 - Se_1_m, m);
        Scalar Kr = Se_sqrt * term * term;
        
        // Clamp to [0, 1]
        if (Kr > 1.0) Kr = 1.0;
        if (Kr < 0.0) Kr = 0.0;
        
        return Kr;
    }
    
    // Compute hydraulic conductivity from head and saturated conductivity
    // K = Ks * Kr(h), where Kr is relative permeability (Mualem model)
    // Uses modified van Genuchten model with air entry value
    KOKKOS_INLINE_FUNCTION
    Scalar compute_K(Scalar h, Scalar Ks, Ordinal idx) const {
        Scalar ha = state.vg_ha(idx);
        
        // Modified van Genuchten: if h >= -ha, then saturated
        if (h >= -ha) {
            return Ks; // Saturated until air entry value
        }
        
        // Compute effective saturation from head
        Scalar Se = compute_Se_from_h(h, idx);
        
        // Compute relative permeability using Mualem model
        Scalar Kr = compute_Kr_Mualem(Se, idx);
        
        return Ks * Kr;
    }
    
    // Compute specific moisture capacity (dwc/dh)
    // C(h) = (θs - θr) * dSe/dh
    // Modified van Genuchten with air entry value
    KOKKOS_INLINE_FUNCTION
    Scalar compute_ch(Ordinal idx) const {
        Scalar h = state.pressure(idx);
        Scalar ha = state.vg_ha(idx);
        
        // Modified van Genuchten: if h >= -ha, then dSe/dh = 0
        if (h >= -ha) {
            return 0.0; // Saturated until air entry value (no change)
        }
        
        Scalar wcs = state.water_content_sat(idx);
        Scalar wcr = state.water_content_res(idx);
        Scalar alpha = state.vg_alpha(idx);
        Scalar n = state.vg_n(idx);
        Scalar m = state.vg_m(idx);
        
        // Derivative of modified van Genuchten water retention curve
        // For h < -ha: Se = [1 + (α(|h| - ha))^n]^(-m)
        // dSe/dh = dSe/d|h| * d|h|/dh
        // Since |h| = -h for h < 0, d|h|/dh = -1
        // dSe/d|h| = -m * n * α^n * (|h| - ha)^(n-1) * [1 + (α(|h| - ha))^n]^(-m-1)
        // Therefore: dSe/dh = m * n * α^n * (|h| - ha)^(n-1) * [1 + (α(|h| - ha))^n]^(-m-1)
        // If ha = 0, this reduces to original formula
        Scalar abs_h = std::abs(h);
        Scalar h_effective = abs_h - ha;  // Effective head below air entry
        
        if (h_effective <= 0.0) {
            return 0.0; // Safety check
        }
        
        Scalar alpha_h_eff_n = std::pow(alpha * h_effective, n);
        Scalar term1 = 1.0 + alpha_h_eff_n;
        Scalar term2 = std::pow(term1, -m - 1.0);
        Scalar h_eff_pow = (n > 1.0) ? std::pow(h_effective, n - 1.0) : 1.0;
        // For h < -ha, derivative is positive (increasing h increases Se)
        Scalar dSe_dh = m * n * std::pow(alpha, n) * h_eff_pow * term2;
        
        // Specific moisture capacity
        Scalar ch = (wcs - wcr) * dSe_dh;
        
        // Ensure non-negative
        if (ch < 0.0) ch = 0.0;
        
        return ch;
    }
    
    // ========================================================================
    // COMPUTE WATER CONTENT FROM HEAD
    // ========================================================================
    void compute_water_content_from_head() {
        auto _pressure = state.pressure;
        auto _water_content = state.water_content;
        auto _water_content_from_head = state.water_content_from_head;
        auto _active_mask_3d = domain.active_mask_3d;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            KOKKOS_LAMBDA (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) > 0) {
                    Scalar h = _pressure(domain_idx);
                    _water_content_from_head(domain_idx) = compute_wch(h, domain_idx);
                }
            });
    }
    
    // ========================================================================
    // COMPUTE HEAD FROM WATER CONTENT
    // ========================================================================
    void compute_head_from_water_content() {
        auto _head_from_water_content = state.head_from_water_content;
        auto _active_mask_3d = domain.active_mask_3d;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            KOKKOS_LAMBDA (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) > 0) {
                    _head_from_water_content(domain_idx) = compute_hwc(domain_idx);
                }
            });
    }
    
    // ========================================================================
    // COMPUTE SPECIFIC MOISTURE CAPACITY
    // ========================================================================
    void compute_specific_moisture_capacity() {
        // This is computed on-the-fly in matrix assembly
        // Placeholder function for now
    }
    
    // ========================================================================
    // COMPUTE HYDRAULIC CONDUCTIVITY ON FACES
    // ========================================================================
    void compute_K_face() {
        auto _pressure = state.pressure;
        auto _conductivity_x = state.conductivity_x;
        auto _conductivity_y = state.conductivity_y;
        auto _conductivity_z = state.conductivity_z;
        auto _conductivity_sat_x = state.conductivity_sat_x;
        auto _conductivity_sat_y = state.conductivity_sat_y;
        auto _conductivity_sat_z = state.conductivity_sat_z;
        auto _active_mask_3d = domain.active_mask_3d;
        
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _neighbor_top = active_mesh.neighbor_top;
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _neighbor_bottom = active_mesh.neighbor_bottom;
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            KOKKOS_LAMBDA (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) {
                    _conductivity_x(domain_idx) = 0.0;
                    _conductivity_y(domain_idx) = 0.0;
                    _conductivity_z(domain_idx) = 0.0;
                    return;
                }
                
                // X-direction (face between i and i+1)
                Ordinal active_right = _neighbor_right(i);
                if (active_right >= 0) {
                    Ordinal idx_right = _active_to_domain(active_right);
                    Scalar h_here = _pressure(domain_idx);
                    Scalar h_right = _pressure(idx_right);
                    Scalar Ks_here = _conductivity_sat_x(domain_idx);
                    Scalar Ks_right = _conductivity_sat_x(idx_right);
                    Scalar K_here = compute_K(h_here, Ks_here, domain_idx);
                    Scalar K_right = compute_K(h_right, Ks_right, idx_right);
                    _conductivity_x(domain_idx) = 0.5 * (K_here + K_right);
                } else {
                    _conductivity_x(domain_idx) = compute_K(_pressure(domain_idx), 
                                                             _conductivity_sat_x(domain_idx), domain_idx);
                }
                
                // Y-direction (face between j and j+1)
                Ordinal active_front = _neighbor_front(i);
                if (active_front >= 0) {
                    Ordinal idx_front = _active_to_domain(active_front);
                    Scalar h_here = _pressure(domain_idx);
                    Scalar h_front = _pressure(idx_front);
                    Scalar Ks_here = _conductivity_sat_y(domain_idx);
                    Scalar Ks_front = _conductivity_sat_y(idx_front);
                    Scalar K_here = compute_K(h_here, Ks_here, domain_idx);
                    Scalar K_front = compute_K(h_front, Ks_front, idx_front);
                    _conductivity_y(domain_idx) = 0.5 * (K_here + K_front);
                } else {
                    _conductivity_y(domain_idx) = compute_K(_pressure(domain_idx), 
                                                             _conductivity_sat_y(domain_idx), domain_idx);
                }
                
                // Z-direction (face between k and k+1)
                Ordinal active_top = _neighbor_top(i);
                if (active_top >= 0) {
                    Ordinal idx_top = _active_to_domain(active_top);
                    Scalar h_here = _pressure(domain_idx);
                    Scalar h_top = _pressure(idx_top);
                    Scalar Ks_here = _conductivity_sat_z(domain_idx);
                    Scalar Ks_top = _conductivity_sat_z(idx_top);
                    Scalar K_here = compute_K(h_here, Ks_here, domain_idx);
                    Scalar K_top = compute_K(h_top, Ks_top, idx_top);
                    _conductivity_z(domain_idx) = 0.5 * (K_here + K_top);
                } else {
                    // This is a top cell (no top neighbor in groundwater domain)
                    _conductivity_z(domain_idx) = compute_K(_pressure(domain_idx), 
                                                             _conductivity_sat_z(domain_idx), domain_idx);
                }
            });
        
        // Apply surface-subsurface interface conductivity (for coupled simulations)
        if (coupled_with_sw_ && sw_state_ != nullptr) {
            apply_surface_interface_conductivity();
        }
    }
    
    // ========================================================================
    // APPLY SURFACE-SUBSURFACE INTERFACE CONDUCTIVITY
    // ========================================================================
    // Adjusts Z-conductivity at the top boundary based on surface water conditions
    // - If surface has water: use average of surface K (saturated) and subsurface K
    // - If surface is dry: use subsurface K or zero based on BC type
    void apply_surface_interface_conductivity() {
        if (!coupled_with_sw_ || sw_state_ == nullptr) return;
        
        auto _conductivity_z = state.conductivity_z;
        auto _conductivity_sat_z = state.conductivity_sat_z;
        auto _pressure = state.pressure;
        auto _active_mask_3d = domain.active_mask_3d;
        
        auto _neighbor_top = active_mesh.neighbor_top;
        auto _active_to_domain = active_mesh.active_to_domain;
        auto _coord_k = active_mesh.coord_k;
        auto _gw_to_sw = domain.mesh->gw_to_sw_idx;
        auto _depth = sw_state_->depth;
        Scalar min_d = min_depth;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            KOKKOS_LAMBDA (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) return;
                
                Ordinal active_top = _neighbor_top(i);
                Ordinal k = _coord_k(i);
                
                // Check if this is a top cell (no top neighbor in groundwater domain)
                bool is_top_cell = (active_top < 0) || (k == domain.nz - 1);
                
                if (is_top_cell) {
                    Ordinal sw_idx = _gw_to_sw(domain_idx);
                    if (sw_idx >= 0) {
                        Scalar surface_depth = _depth(sw_idx);
                        Scalar Ks = _conductivity_sat_z(domain_idx);
                        Scalar K_cell = compute_K(_pressure(domain_idx), Ks, domain_idx);
                        
                        if (surface_depth > min_d) {
                            // Surface has water: use average of saturated K and cell K
                            _conductivity_z(domain_idx) = 0.5 * (Ks + K_cell);
                        } else {
                            // Surface is dry: use cell K (or could be zero for no-flux BC)
                            _conductivity_z(domain_idx) = K_cell;
                        }
                    }
                }
            });
    }
    
    // ========================================================================
    // COMPUTE BAROCLINIC EFFECTS (Density/Viscosity Ratios)
    // ========================================================================
    void baroclinic_face() {
        auto _density_ratio = state.density_ratio;
        auto _viscosity_ratio = state.viscosity_ratio;
        auto _density_ratio_xp = state.density_ratio_xp;
        auto _density_ratio_yp = state.density_ratio_yp;
        auto _density_ratio_zp = state.density_ratio_zp;
        auto _viscosity_ratio_xp = state.viscosity_ratio_xp;
        auto _viscosity_ratio_yp = state.viscosity_ratio_yp;
        auto _viscosity_ratio_zp = state.viscosity_ratio_zp;
        auto _active_mask_3d = domain.active_mask_3d;
        
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _neighbor_top = active_mesh.neighbor_top;
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            KOKKOS_LAMBDA (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) {
                    _density_ratio_xp(domain_idx) = 1.0;
                    _density_ratio_yp(domain_idx) = 1.0;
                    _density_ratio_zp(domain_idx) = 1.0;
                    _viscosity_ratio_xp(domain_idx) = 1.0;
                    _viscosity_ratio_yp(domain_idx) = 1.0;
                    _viscosity_ratio_zp(domain_idx) = 1.0;
                    return;
                }
                
                // X-direction
                Ordinal active_right = _neighbor_right(i);
                if (active_right >= 0) {
                    Ordinal idx_right = _active_to_domain(active_right);
                    _density_ratio_xp(domain_idx) = 0.5 * (_density_ratio(domain_idx) + _density_ratio(idx_right));
                    _viscosity_ratio_xp(domain_idx) = 0.5 * (_viscosity_ratio(domain_idx) + _viscosity_ratio(idx_right));
                } else {
                    _density_ratio_xp(domain_idx) = _density_ratio(domain_idx);
                    _viscosity_ratio_xp(domain_idx) = _viscosity_ratio(domain_idx);
                }
                
                // Y-direction
                Ordinal active_front = _neighbor_front(i);
                if (active_front >= 0) {
                    Ordinal idx_front = _active_to_domain(active_front);
                    _density_ratio_yp(domain_idx) = 0.5 * (_density_ratio(domain_idx) + _density_ratio(idx_front));
                    _viscosity_ratio_yp(domain_idx) = 0.5 * (_viscosity_ratio(domain_idx) + _viscosity_ratio(idx_front));
                } else {
                    _density_ratio_yp(domain_idx) = _density_ratio(domain_idx);
                    _viscosity_ratio_yp(domain_idx) = _viscosity_ratio(domain_idx);
                }
                
                // Z-direction
                Ordinal active_top = _neighbor_top(i);
                if (active_top >= 0) {
                    Ordinal idx_top = _active_to_domain(active_top);
                    _density_ratio_zp(domain_idx) = 0.5 * (_density_ratio(domain_idx) + _density_ratio(idx_top));
                    _viscosity_ratio_zp(domain_idx) = 0.5 * (_viscosity_ratio(domain_idx) + _viscosity_ratio(idx_top));
                } else {
                    _density_ratio_zp(domain_idx) = _density_ratio(domain_idx);
                    _viscosity_ratio_zp(domain_idx) = _viscosity_ratio(domain_idx);
                }
            });
    }
    
    // ========================================================================
    // ASSEMBLE LINEAR SYSTEM (Matrix Coefficients and RHS)
    // ========================================================================
    void assemble_linear_system() {
        compute_matrix_coefficients();
        compute_rhs();
    }
    
    // Compute matrix coefficients (Gxp, Gxm, Gyp, Gym, Gzp, Gzm, Gct)
    void compute_matrix_coefficients() {
        auto _matrix_diag = state.matrix_diag;
        auto _matrix_xp = state.matrix_xp;
        auto _matrix_xm = state.matrix_xm;
        auto _matrix_yp = state.matrix_yp;
        auto _matrix_ym = state.matrix_ym;
        auto _matrix_zp = state.matrix_zp;
        auto _matrix_zm = state.matrix_zm;
        auto _conductivity_x = state.conductivity_x;
        auto _conductivity_y = state.conductivity_y;
        auto _conductivity_z = state.conductivity_z;
        auto _density_ratio_xp = state.density_ratio_xp;
        auto _density_ratio_yp = state.density_ratio_yp;
        auto _density_ratio_zp = state.density_ratio_zp;
        auto _viscosity_ratio_xp = state.viscosity_ratio_xp;
        auto _viscosity_ratio_yp = state.viscosity_ratio_yp;
        auto _viscosity_ratio_zp = state.viscosity_ratio_zp;
        auto _water_content = state.water_content;
        auto _water_content_old = state.water_content_old;
        auto _water_content_sat = state.water_content_sat;
        auto _density_ratio = state.density_ratio;
        auto _active_mask_3d = domain.active_mask_3d;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _neighbor_bottom = active_mesh.neighbor_bottom;
        auto _neighbor_top = active_mesh.neighbor_top;
        auto _active_to_domain = active_mesh.active_to_domain;
        auto _coord_k = active_mesh.coord_k;
        
        // Get layer thicknesses from domain
        auto _dz_layers = domain.layer_thickness;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            KOKKOS_LAMBDA (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) {
                    _matrix_diag(domain_idx) = 1.0;
                    _matrix_xp(domain_idx) = 0.0; _matrix_xm(domain_idx) = 0.0;
                    _matrix_yp(domain_idx) = 0.0; _matrix_ym(domain_idx) = 0.0;
                    _matrix_zp(domain_idx) = 0.0; _matrix_zm(domain_idx) = 0.0;
                    return;
                }
                
                // Get neighbors
                Ordinal active_left = _neighbor_left(i);
                Ordinal active_right = _neighbor_right(i);
                Ordinal active_back = _neighbor_back(i);
                Ordinal active_front = _neighbor_front(i);
                Ordinal active_bottom = _neighbor_bottom(i);
                Ordinal active_top = _neighbor_top(i);
                
                // Get layer index (k)
                Ordinal k = _coord_k(i);
                Scalar dz = _dz_layers(k);
                Scalar dz_avg = dz; // Simplified: use current layer thickness
                
                // Compute specific moisture capacity
                Scalar ch = compute_ch(domain_idx);
                
                // Cell volume
                Scalar V = domain.dx * domain.dy * dz;
                
                // X-direction coefficients
                _matrix_xp(domain_idx) = 0.0;
                if (active_right >= 0) {
                    Scalar Kx = _conductivity_x(domain_idx);
                    Scalar r_rho = _density_ratio_xp(domain_idx);
                    Scalar r_visc = _viscosity_ratio_xp(domain_idx);
                    Scalar Ax = domain.dy * dz; // Face area
                    _matrix_xp(domain_idx) = -Kx * dt * r_rho * r_visc * Ax / (domain.dx * domain.dx);
                }
                
                _matrix_xm(domain_idx) = 0.0;
                if (active_left >= 0) {
                    Ordinal idx_left = _active_to_domain(active_left);
                    Scalar Kx = _conductivity_x(idx_left);
                    Scalar r_rho = _density_ratio_xp(idx_left);
                    Scalar r_visc = _viscosity_ratio_xp(idx_left);
                    Scalar Ax = domain.dy * dz;
                    _matrix_xm(domain_idx) = -Kx * dt * r_rho * r_visc * Ax / (domain.dx * domain.dx);
                }
                
                // Y-direction coefficients
                _matrix_yp(domain_idx) = 0.0;
                if (active_front >= 0) {
                    Scalar Ky = _conductivity_y(domain_idx);
                    Scalar r_rho = _density_ratio_yp(domain_idx);
                    Scalar r_visc = _viscosity_ratio_yp(domain_idx);
                    Scalar Ay = domain.dx * dz;
                    _matrix_yp(domain_idx) = -Ky * dt * r_rho * r_visc * Ay / (domain.dy * domain.dy);
                }
                
                _matrix_ym(domain_idx) = 0.0;
                if (active_back >= 0) {
                    Ordinal idx_back = _active_to_domain(active_back);
                    Scalar Ky = _conductivity_y(idx_back);
                    Scalar r_rho = _density_ratio_yp(idx_back);
                    Scalar r_visc = _viscosity_ratio_yp(idx_back);
                    Scalar Ay = domain.dx * dz;
                    _matrix_ym(domain_idx) = -Ky * dt * r_rho * r_visc * Ay / (domain.dy * domain.dy);
                }
                
                // Z-direction coefficients
                _matrix_zp(domain_idx) = 0.0;
                if (active_top >= 0) {
                    Scalar Kz = _conductivity_z(domain_idx);
                    Scalar r_rho = _density_ratio_zp(domain_idx);
                    Scalar r_visc = _viscosity_ratio_zp(domain_idx);
                    Scalar Az = domain.dx * domain.dy;
                    Scalar dz_top = (k < domain.nz - 1) ? _dz_layers(k + 1) : dz;
                    Scalar dz_avg_z = 0.5 * (dz + dz_top);
                    _matrix_zp(domain_idx) = -Kz * dt * r_rho * r_visc * V / (dz * dz_avg_z);
                }
                
                _matrix_zm(domain_idx) = 0.0;
                if (active_bottom >= 0) {
                    Ordinal idx_bottom = _active_to_domain(active_bottom);
                    Scalar Kz = _conductivity_z(idx_bottom);
                    Scalar r_rho = _density_ratio_zp(idx_bottom);
                    Scalar r_visc = _viscosity_ratio_zp(idx_bottom);
                    Scalar Az = domain.dx * domain.dy;
                    Ordinal k_bottom = _coord_k(active_bottom);
                    Scalar dz_bottom = _dz_layers(k_bottom);
                    Scalar dz_avg_z = 0.5 * (dz + dz_bottom);
                    _matrix_zm(domain_idx) = -Kz * dt * r_rho * r_visc * V / (dz * dz_avg_z);
                }
                
                // Diagonal coefficient
                Scalar ch_term = (ch + Ss * _water_content_old(domain_idx) / _water_content_sat(domain_idx)) * 
                                 _density_ratio(domain_idx);
                _matrix_diag(domain_idx) = ch_term * V;
                _matrix_diag(domain_idx) -= (_matrix_xp(domain_idx) + _matrix_xm(domain_idx) + 
                                            _matrix_yp(domain_idx) + _matrix_ym(domain_idx) +
                                            _matrix_zp(domain_idx) + _matrix_zm(domain_idx));
            });
    }
    
    // Compute right-hand side vector
    void compute_rhs() {
        auto _matrix_rhs = state.matrix_rhs;
        auto _pressure = state.pressure;
        auto _pressure_old = state.pressure_old;
        auto _water_content = state.water_content;
        auto _water_content_old = state.water_content_old;
        auto _water_content_sat = state.water_content_sat;
        auto _density_ratio = state.density_ratio;
        auto _conductivity_z = state.conductivity_z;
        auto _density_ratio_zp = state.density_ratio_zp;
        auto _viscosity_ratio_zp = state.viscosity_ratio_zp;
        auto _active_mask_3d = domain.active_mask_3d;
        
        auto _neighbor_bottom = active_mesh.neighbor_bottom;
        auto _active_to_domain = active_mesh.active_to_domain;
        auto _coord_k = active_mesh.coord_k;
        auto _dz_layers = domain.layer_thickness;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            KOKKOS_LAMBDA (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) {
                    _matrix_rhs(domain_idx) = _pressure_old(domain_idx);
                    return;
                }
                
                // Get layer index
                Ordinal k = _coord_k(i);
                Scalar dz = _dz_layers(k);
                Scalar V = domain.dx * domain.dy * dz;
                
                // Compute specific moisture capacity
                Scalar ch = compute_ch(domain_idx);
                
                // Base term: storage term
                Scalar ch_term = (ch + Ss * _water_content_old(domain_idx) / _water_content_sat(domain_idx)) * 
                                _density_ratio(domain_idx);
                _matrix_rhs(domain_idx) = ch_term * _pressure_old(domain_idx) * V;
                
                // Density-driven flow term (vertical)
                _matrix_rhs(domain_idx) -= dt * V * _conductivity_z(domain_idx) * 
                                          _density_ratio_zp(domain_idx) * _density_ratio_zp(domain_idx) * 
                                          _viscosity_ratio_zp(domain_idx) / dz;
                
                Ordinal active_bottom = _neighbor_bottom(i);
                if (active_bottom >= 0) {
                    Ordinal idx_bottom = _active_to_domain(active_bottom);
                    Ordinal k_bottom = _coord_k(active_bottom);
                    Scalar dz_bottom = _dz_layers(k_bottom);
                    _matrix_rhs(domain_idx) += dt * V * _conductivity_z(idx_bottom) * 
                                              _density_ratio_zp(idx_bottom) * _density_ratio_zp(idx_bottom) * 
                                              _viscosity_ratio_zp(idx_bottom) / dz_bottom;
                }
                
                // Boundary conditions and source/sink terms will be added here later
            });
    }
    
    // Add surface water depth contribution to RHS (for coupled mode)
    // When there's surface water above, the head at top boundary = bottom + depth
    // This is added as: RHS -= Gzm * (surface_depth)
    void add_surface_coupling_to_rhs() {
        if (!coupled_with_sw_ || sw_state_ == nullptr) return;
        
        auto _matrix_rhs = state.matrix_rhs;
        auto _matrix_zm = state.matrix_zm;
        auto _neighbor_top = active_mesh.neighbor_top;
        auto _active_to_domain = active_mesh.active_to_domain;
        auto _gw_to_sw = domain.mesh->gw_to_sw_idx;
        auto _depth = sw_state_->depth;
        Scalar min_d = min_depth;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            KOKKOS_LAMBDA (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                
                // Check if this is a top cell (no neighbor above)
                Ordinal active_top = _neighbor_top(i);
                if (active_top < 0) {
                    // This is a top cell - check for surface water
                    Ordinal sw_idx = _gw_to_sw(domain_idx);
                    if (sw_idx >= 0) {
                        Scalar surf_depth = _depth(sw_idx);
                        if (surf_depth > min_d) {
                            // Surface water present - add depth as boundary condition
                            // The surface water elevation acts as fixed head
                            _matrix_rhs(domain_idx) -= _matrix_zm(domain_idx) * surf_depth;
                        }
                    }
                }
            });
    }
    
    // Add terrain slope effects to flux computation (for follow_terrain mode)
    // When follow_terrain=true, horizontal flux includes terrain slope contribution
    // Legacy: q = K * visc * ((dh/dx) * cos + sin)  where sin/cos are terrain slopes
    void add_terrain_slope_effects() {
        if (!follow_terrain) return;
        
        // Terrain slope effects are already included in compute_groundwater_flux()
        // through the face area calculations and cos/sin terms
        // For inclined domains, the flux formula becomes:
        // qx = K * r_visc * A * (dh/dx * cos_x + sign * sin_x)
        // This is a simplified implementation - full implementation would need
        // terrain slope angles (sin_x, cos_x, sin_y, cos_y) stored in domain
        // For now, this is a placeholder for future enhancement
    }
    
    // ========================================================================
    // SOLVE LINEAR SYSTEM
    // ========================================================================
    void solve_linear_system() {
        // Initialize solution vector with current head as initial guess
        auto _pressure = state.pressure;
        auto _active_to_domain = active_mesh.active_to_domain;
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            KOKKOS_LAMBDA (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                solution(i) = _pressure(domain_idx);
            });
        
        pcg_solver->solve_3d(active_mesh,
                            state.matrix_diag,
                            state.matrix_xp,
                            state.matrix_xm,
                            state.matrix_yp,
                            state.matrix_ym,
                            state.matrix_zp,
                            state.matrix_zm,
                            state.matrix_rhs,
                            solution);
        
        // Copy solution back to state.pressure (h_new)
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            KOKKOS_LAMBDA (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                _pressure(domain_idx) = solution(i);
            });
    }
    
    // ========================================================================
    // APPLY BOUNDARY CONDITIONS TO MATRIX
    // ========================================================================
    void apply_boundary_conditions_to_matrix(Scalar current_time) {
        if (!bc_manager_) return;
        
        auto _matrix_diag = state.matrix_diag;
        auto _matrix_xp = state.matrix_xp;
        auto _matrix_xm = state.matrix_xm;
        auto _matrix_yp = state.matrix_yp;
        auto _matrix_ym = state.matrix_ym;
        auto _matrix_zp = state.matrix_zp;
        auto _matrix_zm = state.matrix_zm;
        auto _matrix_rhs = state.matrix_rhs;
        
        const auto& bcs = bc_manager_->get_boundary_conditions();
        Scalar dt_local = dt;
        
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
            
            if (bc.type == GwBcType::FIXED_HEAD) {
                // Dirichlet BC: prescribed hydraulic head
                Kokkos::parallel_for(RangePolicy(0, n_bc_cells),
                    KOKKOS_LAMBDA (const Ordinal i) {
                        Ordinal cell_idx = d_bc_indices(i);
                        _matrix_diag(cell_idx) = 1.0;
                        _matrix_xp(cell_idx) = 0.0;
                        _matrix_xm(cell_idx) = 0.0;
                        _matrix_yp(cell_idx) = 0.0;
                        _matrix_ym(cell_idx) = 0.0;
                        _matrix_zp(cell_idx) = 0.0;
                        _matrix_zm(cell_idx) = 0.0;
                        _matrix_rhs(cell_idx) = bc_value;
                    });
                    
            } else if (bc.type == GwBcType::FIXED_FLUX) {
                // Neumann BC: prescribed flux (add to RHS)
                Kokkos::parallel_for(RangePolicy(0, n_bc_cells),
                    KOKKOS_LAMBDA (const Ordinal i) {
                        Ordinal cell_idx = d_bc_indices(i);
                        Scalar flux_volume = bc_value * dt_local;
                        _matrix_rhs(cell_idx) += flux_volume;
                    });
            }
        }
    }
    
    // ========================================================================
    // ENFORCE BOUNDARY CONDITIONS ON SOLUTION
    // ========================================================================
    void enforce_boundary_conditions_on_solution(Scalar current_time) {
        if (!bc_manager_) return;
        
        auto _pressure = state.pressure;
        const auto& bcs = bc_manager_->get_boundary_conditions();
        
        // Process each boundary condition using device-side parallel operations
        for (const auto& bc : bcs) {
            if (bc.type != GwBcType::FIXED_HEAD) continue;
            
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
            
            Kokkos::parallel_for(RangePolicy(0, n_bc_cells),
                KOKKOS_LAMBDA (const Ordinal i) {
                    Ordinal cell_idx = d_bc_indices(i);
                    _pressure(cell_idx) = bc_value;
                });
        }
    }
    
    // ========================================================================
    // COMPUTE GROUNDWATER FLUX
    // ========================================================================
    void compute_groundwater_flux() {
        // Compute Darcy fluxes in x, y, z directions
        // This is a simplified version - full implementation would use darcy_flux function
        auto _pressure = state.pressure;
        auto _conductivity_x = state.conductivity_x;
        auto _conductivity_y = state.conductivity_y;
        auto _conductivity_z = state.conductivity_z;
        auto _density_ratio_xp = state.density_ratio_xp;
        auto _density_ratio_yp = state.density_ratio_yp;
        auto _density_ratio_zp = state.density_ratio_zp;
        auto _viscosity_ratio_xp = state.viscosity_ratio_xp;
        auto _viscosity_ratio_yp = state.viscosity_ratio_yp;
        auto _viscosity_ratio_zp = state.viscosity_ratio_zp;
        auto _flux_x = state.flux_x;
        auto _flux_y = state.flux_y;
        auto _flux_z = state.flux_z;
        auto _active_mask_3d = domain.active_mask_3d;
        
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _neighbor_top = active_mesh.neighbor_top;
        auto _active_to_domain = active_mesh.active_to_domain;
        auto _coord_k = active_mesh.coord_k;
        auto _dz_layers = domain.layer_thickness;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            KOKKOS_LAMBDA (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) {
                    _flux_x(domain_idx) = 0.0;
                    _flux_y(domain_idx) = 0.0;
                    _flux_z(domain_idx) = 0.0;
                    return;
                }
                
                Ordinal k = _coord_k(i);
                Scalar dz = _dz_layers(k);
                
                // X-direction flux (Darcy's law)
                // Convention: q = K * r_visc * A * (h_neighbor - h_here) / dx
                // Positive flux = flow toward current cell (inflow from neighbor)
                Ordinal active_right = _neighbor_right(i);
                if (active_right >= 0) {
                    Ordinal idx_right = _active_to_domain(active_right);
                    Scalar dh = _pressure(idx_right) - _pressure(domain_idx);
                    Scalar Kx = _conductivity_x(domain_idx);
                    Scalar r_visc = _viscosity_ratio_xp(domain_idx);
                    Scalar Ax = domain.dy * dz;
                    _flux_x(domain_idx) = Kx * r_visc * Ax * dh / domain.dx;
                } else {
                    _flux_x(domain_idx) = 0.0;
                }
                
                // Y-direction flux
                Ordinal active_front = _neighbor_front(i);
                if (active_front >= 0) {
                    Ordinal idx_front = _active_to_domain(active_front);
                    Scalar dh = _pressure(idx_front) - _pressure(domain_idx);
                    Scalar Ky = _conductivity_y(domain_idx);
                    Scalar r_visc = _viscosity_ratio_yp(domain_idx);
                    Scalar Ay = domain.dx * dz;
                    _flux_y(domain_idx) = Ky * r_visc * Ay * dh / domain.dy;
                } else {
                    _flux_y(domain_idx) = 0.0;
                }
                
                // Z-direction flux (includes buoyancy/density term)
                // q = K * r_visc * A * (dh/dz - r_rho)
                // The -r_rho term accounts for density-driven flow (buoyancy)
                Ordinal active_top = _neighbor_top(i);
                if (active_top >= 0) {
                    Ordinal idx_top = _active_to_domain(active_top);
                    Scalar dh = _pressure(idx_top) - _pressure(domain_idx);
                    Scalar Kz = _conductivity_z(domain_idx);
                    Scalar r_rho = _density_ratio_zp(domain_idx);
                    Scalar r_visc = _viscosity_ratio_zp(domain_idx);
                    Scalar Az = domain.dx * domain.dy;
                    Scalar dz_top = (k < domain.nz - 1) ? _dz_layers(k + 1) : dz;
                    Scalar dz_avg = 0.5 * (dz + dz_top);
                    _flux_z(domain_idx) = Kz * r_visc * Az * (dh / dz_avg - r_rho);
                } else {
                    _flux_z(domain_idx) = 0.0;
                }
            });
    }
    
    // ========================================================================
    // CHECK AVAILABLE ROOM (Pore Space)
    // ========================================================================
    void check_room() {
        auto _water_content = state.water_content;
        auto _water_content_sat = state.water_content_sat;
        auto _room = state.room;
        auto _active_mask_3d = domain.active_mask_3d;
        auto _dz_layers = domain.layer_thickness;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        auto _coord_k = active_mesh.coord_k;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            KOKKOS_LAMBDA (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) {
                    _room(domain_idx) = 0.0;
                    return;
                }
                
                Ordinal k = _coord_k(i);
                Scalar dz = _dz_layers(k);
                Scalar V = domain.dx * domain.dy * dz;
                _room(domain_idx) = (_water_content_sat(domain_idx) - _water_content(domain_idx)) * V;
            });
    }
    
    // ========================================================================
public:
    // COMPUTE SEEPAGE FLUX FROM TOP BOUNDARY
    // ========================================================================
    // Computes seepage rate (qss) from top boundary flux for surface-subsurface coupling
    // This should be called after compute_groundwater_flux()
    // sw_state: Surface water state variables (to store seepage_rate)
    // sw_domain: Surface water domain (for mapping)
    void compute_seepage_flux(SwStateVariables& sw_state, const SwDomain& sw_domain) {
        auto _flux_z = state.flux_z;
        auto _water_content = state.water_content;
        auto _water_content_sat = state.water_content_sat;
        auto _water_content_res = state.water_content_res;
        auto _seepage_top = state.seepage_top;
        auto _active_mask_3d = domain.active_mask_3d;
        auto _dz_layers = domain.layer_thickness;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        auto _coord_k = active_mesh.coord_k;
        auto _neighbor_top = active_mesh.neighbor_top;
        
        // Surface water state views
        auto _seepage = sw_state.seepage;
        auto _seepage_old = sw_state.seepage_old;
        auto _seepage_rate = sw_state.seepage_rate;
        auto _reset_seepage = sw_state.reset_seepage;
        auto _depth = sw_state.depth;
        auto _area_top = sw_state.area_top;
        
        // Mesh coupling maps
        auto _gw_to_sw = domain.mesh->gw_to_sw_idx;
        auto _sw_to_gw = sw_domain.mesh->sw_to_gw_idx;
        
        // Create host mirrors for surface water (need to update seepage)
        auto h_seepage = Kokkos::create_mirror_view(_seepage);
        auto h_seepage_old = Kokkos::create_mirror_view(_seepage_old);
        auto h_seepage_rate = Kokkos::create_mirror_view(_seepage_rate);
        auto h_reset_seepage = Kokkos::create_mirror_view(_reset_seepage);
        auto h_depth = Kokkos::create_mirror_view(_depth);
        auto h_area_top = Kokkos::create_mirror_view(_area_top);
        
        Kokkos::deep_copy(h_seepage_old, _seepage_old);
        Kokkos::deep_copy(h_seepage, _seepage);
        Kokkos::deep_copy(h_reset_seepage, _reset_seepage);
        Kokkos::deep_copy(h_depth, _depth);
        Kokkos::deep_copy(h_area_top, _area_top);
        
        // Device views for groundwater (read-only)
        auto h_flux_z = Kokkos::create_mirror_view(_flux_z);
        auto h_water_content = Kokkos::create_mirror_view(_water_content);
        auto h_water_content_sat = Kokkos::create_mirror_view(_water_content_sat);
        auto h_water_content_res = Kokkos::create_mirror_view(_water_content_res);
        auto h_seepage_top = Kokkos::create_mirror_view(_seepage_top);
        auto h_active_mask_3d = Kokkos::create_mirror_view(_active_mask_3d);
        auto h_coord_k = Kokkos::create_mirror_view(_coord_k);
        auto h_neighbor_top = Kokkos::create_mirror_view(_neighbor_top);
        auto h_active_to_domain = Kokkos::create_mirror_view(_active_to_domain);
        auto h_gw_to_sw = Kokkos::create_mirror_view(_gw_to_sw);
        
        Kokkos::deep_copy(h_flux_z, _flux_z);
        Kokkos::deep_copy(h_water_content, _water_content);
        Kokkos::deep_copy(h_water_content_sat, _water_content_sat);
        Kokkos::deep_copy(h_water_content_res, _water_content_res);
        Kokkos::deep_copy(h_seepage_top, _seepage_top);
        Kokkos::deep_copy(h_active_mask_3d, _active_mask_3d);
        Kokkos::deep_copy(h_coord_k, _coord_k);
        Kokkos::deep_copy(h_neighbor_top, _neighbor_top);
        Kokkos::deep_copy(h_active_to_domain, _active_to_domain);
        Kokkos::deep_copy(h_gw_to_sw, _gw_to_sw);
        
        // Process each active groundwater cell
        for (Ordinal i = 0; i < active_mesh.num_active; ++i) {
            Ordinal domain_idx = h_active_to_domain(i);
            if (h_active_mask_3d(domain_idx) == 0) continue;
            
            Ordinal k = h_coord_k(i);
            Ordinal active_top = h_neighbor_top(i);
            
            // Check if this is a top cell (no top neighbor or k == nz-1)
            bool is_top_cell = (active_top < 0) || (k == domain.nz - 1);
            
            if (is_top_cell) {
                // Get corresponding surface cell index
                Ordinal sw_idx = h_gw_to_sw(domain_idx);
                if (sw_idx < 0 || sw_idx >= sw_domain.num_cells_total) continue;
                
                // Get top boundary flux (flux_z at the top face)
                // For top cells, the flux_z represents flux out of the top
                Scalar qz_top = h_flux_z(domain_idx);
                
                // Get cell area (Az)
                Scalar Az = domain.dx * domain.dy;
                if (h_area_top(sw_idx) > 0.0) {
                    Az = h_area_top(sw_idx);
                }
                
                // Handle reset flag
                if (h_reset_seepage(sw_idx) == 1) {
                    h_seepage(sw_idx) = 0.0;
                    h_reset_seepage(sw_idx) = 0;
                }
                
                // Accumulate seepage: qseepage += qz / Az
                // qz is already a flux (volume/time), so dividing by area gives velocity
                h_seepage(sw_idx) += qz_top / Az;
                
                // Handle special case: if infiltration (qz < 0) and surface is dry
                // and there's not enough pore space, add water to surface
                if (qz_top < 0.0) {
                    Scalar wc = h_water_content(domain_idx);
                    Scalar wcs = h_water_content_sat(domain_idx);
                    auto h_dz_layers = Kokkos::create_mirror_view(domain.layer_thickness);
                    Kokkos::deep_copy(h_dz_layers, domain.layer_thickness);
                    Scalar dz_layer = h_dz_layers(k);
                    Scalar pore_space = (wcs - wc) * domain.dx * domain.dy * dz_layer;
                    Scalar vseep = std::abs(qz_top) * dt * wcs;
                    
                    if (vseep > pore_space && h_depth(sw_idx) <= min_depth) {
                        // Add water to surface (this is handled in legacy code)
                        // For now, we'll let the surface solver handle it
                    }
                }
                
                // Store seepage rate
                h_seepage_rate(sw_idx) = h_seepage(sw_idx);
            }
        }
        
        // Copy back to device
        Kokkos::deep_copy(_seepage, h_seepage);
        Kokkos::deep_copy(_seepage_rate, h_seepage_rate);
        Kokkos::deep_copy(_reset_seepage, h_reset_seepage);
    }
    
    // ========================================================================
    // UPDATE WATER CONTENT (From Fluxes)
    // ========================================================================
    void update_water_content() {
        auto _water_content = state.water_content;
        auto _water_content_old = state.water_content_old;
        auto _water_content_sat = state.water_content_sat;
        auto _flux_x = state.flux_x;
        auto _flux_y = state.flux_y;
        auto _flux_z = state.flux_z;
        auto _density_ratio = state.density_ratio;
        auto _density_ratio_xp = state.density_ratio_xp;
        auto _density_ratio_yp = state.density_ratio_yp;
        auto _density_ratio_zp = state.density_ratio_zp;
        auto _pressure = state.pressure;
        auto _pressure_old = state.pressure_old;
        auto _active_mask_3d = domain.active_mask_3d;
        auto _dz_layers = domain.layer_thickness;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _neighbor_bottom = active_mesh.neighbor_bottom;
        auto _active_to_domain = active_mesh.active_to_domain;
        auto _coord_k = active_mesh.coord_k;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            KOKKOS_LAMBDA (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) {
                    return;
                }
                
                Ordinal k = _coord_k(i);
                Scalar dz = _dz_layers(k);
                Scalar V = domain.dx * domain.dy * dz;
                
                // Compute specific moisture capacity
                Scalar ch = compute_ch(domain_idx);
                
                // Storage coefficient
                Scalar coeff = V * (_density_ratio(domain_idx) + 
                                   _density_ratio(domain_idx) * Ss * 
                                   (_pressure(domain_idx) - _pressure_old(domain_idx)) / 
                                   _water_content_sat(domain_idx));
                
                // Flux divergence
                Scalar dqx = 0.0;
                Ordinal active_left = _neighbor_left(i);
                if (active_left >= 0) {
                    Ordinal idx_left = _active_to_domain(active_left);
                    dqx = dt * (_flux_x(domain_idx) * _density_ratio_xp(domain_idx) - 
                               _flux_x(idx_left) * _density_ratio_xp(idx_left));
                } else {
                    dqx = dt * _flux_x(domain_idx) * _density_ratio_xp(domain_idx);
                }
                
                Scalar dqy = 0.0;
                Ordinal active_back = _neighbor_back(i);
                if (active_back >= 0) {
                    Ordinal idx_back = _active_to_domain(active_back);
                    dqy = dt * (_flux_y(domain_idx) * _density_ratio_yp(domain_idx) - 
                               _flux_y(idx_back) * _density_ratio_yp(idx_back));
                } else {
                    dqy = dt * _flux_y(domain_idx) * _density_ratio_yp(domain_idx);
                }
                
                Scalar dqz = 0.0;
                Ordinal active_bottom = _neighbor_bottom(i);
                if (active_bottom >= 0) {
                    Ordinal idx_bottom = _active_to_domain(active_bottom);
                    dqz = dt * (_flux_z(domain_idx) * _density_ratio_zp(domain_idx) - 
                               _flux_z(idx_bottom) * _density_ratio_zp(idx_bottom));
                } else {
                    dqz = dt * _flux_z(domain_idx) * _density_ratio_zp(domain_idx);
                }
                
                // Update water content
                _water_content(domain_idx) = (_water_content_old(domain_idx) * _density_ratio(domain_idx) * V + 
                                             dqx + dqy + dqz) / coeff;
                
                // Clamp to valid range
                if (_water_content(domain_idx) > _water_content_sat(domain_idx)) {
                    _water_content(domain_idx) = _water_content_sat(domain_idx);
                }
                if (_water_content(domain_idx) < state.water_content_res(domain_idx)) {
                    _water_content(domain_idx) = state.water_content_res(domain_idx);
                }
            });
    }
    
public:
    // ========================================================================
    // ADAPTIVE TIME STEPPING
    // ========================================================================
    // Adjusts time step based on flux divergence errors and Courant number
    // Following legacy implementation in groundwater.c:adaptive_time_step()
    // Returns the new recommended time step
    Scalar adaptive_time_step(Scalar dt_min, Scalar dt_max, Scalar Co_max, 
                               bool sync_coupling, Scalar surface_dt,
                               const SwStateVariables* sw_state = nullptr,
                               const SwDomain* sw_domain = nullptr) {
        // Parameters from legacy code
        Scalar r_red = 0.75;      // Reduction factor
        Scalar r_inc = 1.25;      // Increase factor
        Scalar q_target = 5e-4;   // Target flux for async coupling
        
        Scalar dq_max = 0.0;
        Scalar dt_Comin = 1e8;
        Scalar err_max = 0.0;
        Scalar dt_new = dt;
        Scalar dt_sync = dt_max;
        
        // Get state variables
        auto _pressure = state.pressure;
        auto _pressure_old = state.pressure_old;
        auto _head_prev_prev = state.head_prev_prev;
        auto _water_content = state.water_content;
        auto _water_content_sat = state.water_content_sat;
        auto _flux_x = state.flux_x;
        auto _flux_y = state.flux_y;
        auto _flux_z = state.flux_z;
        auto _conductivity_z = state.conductivity_z;
        auto _active_mask_3d = domain.active_mask_3d;
        auto _dz_layers = domain.layer_thickness;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _neighbor_bottom = active_mesh.neighbor_bottom;
        auto _active_to_domain = active_mesh.active_to_domain;
        auto _coord_k = active_mesh.coord_k;
        
        // Get face areas from domain (need to compute or store)
        // For now, we'll compute them on the fly
        Scalar dt_prev = dt;  // Previous time step (dtn in legacy)
        
        // Compute errors and flux divergence using separate reductions
        // Reduction 1: Maximum flux divergence error
        Kokkos::parallel_reduce(
            RangePolicy(0, active_mesh.num_active),
            KOKKOS_LAMBDA (const Ordinal i, Scalar& local_dq_max) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) return;
                
                // Flux divergence (only for interior cells, not top cells)
                Ordinal k = _coord_k(i);
                if (k > 0) {  // Not top cell
                    Scalar dz = _dz_layers(k);
                    Scalar Ax = domain.dy * dz;
                    Scalar Ay = domain.dx * dz;
                    Scalar Az = domain.dx * domain.dy;
                    
                    // Compute flux divergence
                    Scalar qin = 0.0, qou = 0.0;
                    
                    // X-direction
                    qin += _flux_x(domain_idx) / Ax;
                    Ordinal active_left = _neighbor_left(i);
                    if (active_left >= 0) {
                        Ordinal idx_left = _active_to_domain(active_left);
                        qou += _flux_x(idx_left) / Ax;
                    }
                    
                    // Y-direction
                    qin += _flux_y(domain_idx) / Ay;
                    Ordinal active_back = _neighbor_back(i);
                    if (active_back >= 0) {
                        Ordinal idx_back = _active_to_domain(active_back);
                        qou += _flux_y(idx_back) / Ay;
                    }
                    
                    // Z-direction
                    qin += _flux_z(domain_idx) / Az;
                    Ordinal active_bottom = _neighbor_bottom(i);
                    if (active_bottom >= 0) {
                        Ordinal idx_bottom = _active_to_domain(active_bottom);
                        qou += _flux_z(idx_bottom) / Az;
                    }
                    
                    Scalar dq = std::abs(qin - qou) * dt / dz;
                    if (dq > local_dq_max) local_dq_max = dq;
                }
            },
            Kokkos::Max<Scalar>(dq_max)
        );
        
        // Reduction 2: Minimum Courant-limited time step
        auto _water_content_res = state.water_content_res;
        auto _conductivity_sat_z = state.conductivity_sat_z;
        auto _vg_alpha = state.vg_alpha;
        auto _vg_n = state.vg_n;
        auto _vg_m = state.vg_m;
        auto _vg_ha = state.vg_ha;
        
        Kokkos::parallel_reduce(
            RangePolicy(0, active_mesh.num_active),
            KOKKOS_LAMBDA (const Ordinal i, Scalar& local_dt_Comin) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) return;
                
                Ordinal k = _coord_k(i);
                if (k > 0) {  // Not top cell
                    Scalar dz = _dz_layers(k);
                    Scalar wc = _water_content(domain_idx);
                    Scalar wcs = _water_content_sat(domain_idx);
                    
                    // Courant number constraint for unsaturated flow
                    if (wc < wcs) {
                        // Compute dK/dwc inline (can't call member function from lambda)
                        Scalar wcr = _water_content_res(domain_idx);
                        Scalar Ks = _conductivity_sat_z(domain_idx);
                        Scalar alpha = _vg_alpha(domain_idx);
                        Scalar n = _vg_n(domain_idx);
                        Scalar m = _vg_m(domain_idx);
                        Scalar ha = _vg_ha(domain_idx);
                        
                        // Lambda limiter
                        Scalar wc_lim = wc;
                        if (wc > 0.9999 * wcs && wc < wcs) {
                            wc_lim = 0.9999 * wcs;
                        }
                        
                        Scalar s = (wc_lim - wcr) / (wcs - wcr);
                        Scalar dKdwc = 0.0;
                        
                        if (s > 0.0 && s < 1.0) {
                            Scalar term0, term1, term2;
                            
                            if (ha > 0.0) {
                                // Modified van Genuchten
                                Scalar wcm = wcr + (wcs - wcr) * std::pow((1.0 + std::pow(ha * alpha, n)), m);
                                Scalar c2 = (wcs - wcr) / (wcm - wcr);
                                Scalar c1 = 1.0 / std::pow(1.0 - std::pow((1.0 - std::pow(c2, 1.0/m)), m), 2.0);
                                
                                term0 = std::pow(1.0 - std::pow(c2 * s, 1.0/m), m);
                                term1 = 0.5 * c1 * Ks * std::pow(s, -0.5) * (1.0 - term0) * (1.0 - term0);
                                term2 = 2.0 * c1 * Ks * c2 * std::pow(s, 0.5) * 
                                        std::pow(c2 * s, 1.0/m - 1.0) * (1.0 - term0) * 
                                        std::pow(1.0 - std::pow(c2 * s, 1.0/m), m - 1.0);
                            } else {
                                // Standard van Genuchten
                                term0 = std::pow(1.0 - std::pow(s, 1.0/m), m);
                                term1 = 0.5 * Ks * std::pow(s, -0.5) * (1.0 - term0) * (1.0 - term0);
                                term2 = 2.0 * Ks * std::pow(s, (2.0 - m) / (2.0 * m)) * (1.0 - term0) * 
                                        std::pow(1.0 - std::pow(s, 1.0/m), m - 1.0);
                            }
                            
                            dKdwc = (term1 + term2) / (wcs - wcr);
                            if (dKdwc < 0.0) dKdwc = 0.0;
                        }
                        
                        if (dKdwc > 0.0) {
                            Scalar dt_Co = Co_max * dz / dKdwc;
                            if (dt_Co < local_dt_Comin) local_dt_Comin = dt_Co;
                        }
                    }
                }
            },
            Kokkos::Min<Scalar>(dt_Comin)
        );
        
        // Reduction 3: Maximum error estimate
        Kokkos::parallel_reduce(
            RangePolicy(0, active_mesh.num_active),
            KOKKOS_LAMBDA (const Ordinal i, Scalar& local_err_max) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) return;
                
                // Error estimate (second-order time derivative)
                Scalar h = _pressure(domain_idx);
                Scalar hn = _pressure_old(domain_idx);
                Scalar hnm = _head_prev_prev(domain_idx);
                Scalar err = 0.5 * dt * std::abs((h - hn) / dt - (hn - hnm) / dt_prev);
                if (err > local_err_max) local_err_max = err;
            },
            Kokkos::Max<Scalar>(err_max)
        );
        
        // Adjust dt due to async coupling (if not synchronous)
        if (!sync_coupling && sw_state && sw_domain) {
            auto _depth = sw_state->depth;
            auto _pressure_sw = sw_state->pressure;
            auto _pressure_sw_old = sw_state->pressure_old;
            auto _gw_to_sw = domain.mesh->gw_to_sw_idx;
            
            // Find minimum dt_eta from surface changes
            Kokkos::parallel_reduce(
                RangePolicy(0, active_mesh.num_active),
                KOKKOS_LAMBDA (const Ordinal i, Scalar& local_dt_sync) {
                    Ordinal domain_idx = _active_to_domain(i);
                    if (_active_mask_3d(domain_idx) == 0) return;
                    
                    Ordinal k = _coord_k(i);
                    if (k == 0) {  // Top cell
                        Ordinal sw_idx = _gw_to_sw(domain_idx);
                        if (sw_idx >= 0 && sw_idx < sw_domain->num_cells_total) {
                            Scalar depth = _depth(sw_idx);
                            Scalar eta = _pressure_sw(sw_idx);
                            Scalar etan = _pressure_sw_old(sw_idx);
                            
                            if (depth > 0.0 && std::abs(eta - etan) > 1e-10) {
                                Scalar dz = _dz_layers(k);
                                Scalar Kz = _conductivity_z(domain_idx);
                                if (Kz > 0.0) {
                                    Scalar dt_eta = std::abs(q_target * dz * surface_dt / Kz / (eta - etan));
                                    if (dt_eta < local_dt_sync) local_dt_sync = dt_eta;
                                }
                            }
                        }
                    }
                },
                Kokkos::Min<Scalar>(dt_sync)
            );
        }
        
        // Adjust dt based on flux divergence
        if (dq_max > 0.02) {
            dt_new = dt * r_red;
        } else if (dq_max < 0.01) {
            dt_new = dt * r_inc;
        }
        
        // Apply Courant number constraint
        if (dt_new > dt_Comin) {
            dt_new = dt_Comin;
        }
        
        // Apply async coupling constraint (if not synchronous)
        if (!sync_coupling && dt_new > dt_sync) {
            dt_new = dt_sync;
        }
        
        // Clamp to bounds
        if (dt_new > dt_max) dt_new = dt_max;
        if (dt_new < dt_min) dt_new = dt_min;
        
        return dt_new;
    }
    
private:
    // Helper function to compute dK/dwc for adaptive time stepping
    // Following legacy implementation in subroutines.c:compute_dKdwc()
    KOKKOS_INLINE_FUNCTION
    Scalar compute_dKdwc(Ordinal idx) const {
        Scalar wc = state.water_content(idx);
        Scalar wcs = state.water_content_sat(idx);
        Scalar wcr = state.water_content_res(idx);
        Scalar Ks = state.conductivity_sat_z(idx);
        Scalar alpha = state.vg_alpha(idx);
        Scalar n = state.vg_n(idx);
        Scalar m = state.vg_m(idx);
        Scalar ha = state.vg_ha(idx);
        
        // Lambda limiter (from legacy code)
        if (wc > 0.9999 * wcs && wc < wcs) {
            wc = 0.9999 * wcs;
        }
        
        Scalar s = (wc - wcr) / (wcs - wcr);
        if (s <= 0.0 || s >= 1.0) {
            return 0.0;  // At bounds, derivative is zero
        }
        
        Scalar term0, term1, term2, dKdwc;
        
        if (ha > 0.0) {
            // Modified van Genuchten
            Scalar wcm = wcr + (wcs - wcr) * std::pow((1.0 + std::pow(ha * alpha, n)), m);
            Scalar c2 = (wcs - wcr) / (wcm - wcr);
            Scalar c1 = 1.0 / std::pow(1.0 - std::pow((1.0 - std::pow(c2, 1.0/m)), m), 2.0);
            
            term0 = std::pow(1.0 - std::pow(c2 * s, 1.0/m), m);
            term1 = 0.5 * c1 * Ks * std::pow(s, -0.5) * (1.0 - term0) * (1.0 - term0);
            term2 = 2.0 * c1 * Ks * c2 * std::pow(s, 0.5) * 
                    std::pow(c2 * s, 1.0/m - 1.0) * (1.0 - term0) * 
                    std::pow(1.0 - std::pow(c2 * s, 1.0/m), m - 1.0);
        } else {
            // Standard van Genuchten
            term0 = std::pow(1.0 - std::pow(s, 1.0/m), m);
            term1 = 0.5 * Ks * std::pow(s, -0.5) * (1.0 - term0) * (1.0 - term0);
            term2 = 2.0 * Ks * std::pow(s, (2.0 - m) / (2.0 * m)) * (1.0 - term0) * 
                    std::pow(1.0 - std::pow(s, 1.0/m), m - 1.0);
        }
        
        dKdwc = (term1 + term2) / (wcs - wcr);
        
        // Ensure non-negative
        if (dKdwc < 0.0) dKdwc = 0.0;
        
        return dKdwc;
    }
};

#endif // FREHG_GROUNDWATER_SOLVER_HPP

