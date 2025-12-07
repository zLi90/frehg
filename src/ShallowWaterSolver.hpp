#ifndef FREHG_SHALLOW_WATER_SOLVER_HPP
#define FREHG_SHALLOW_WATER_SOLVER_HPP

#include "define.hpp"
#include "Domain.hpp"
#include "StateVariables.hpp"
#include "ActiveCellMesh.hpp"
#include "PCGSolver.hpp"
#include "BoundaryConditions.hpp"
#include "SourceSinkTerms.hpp"
#include <memory>
#include <cmath>

// ============================================================================
//                      SHALLOW WATER EQUATION SOLVER
// ============================================================================
// Solves the 2D depth-integrated Navier-Stokes equations (shallow water)
// Uses Kokkos for parallelization and PCG for linear system solution

class ShallowWaterSolver {
public:
    // --- Domain and Mesh ---
    SwDomain& domain;
    ActiveCellMesh& active_mesh;
    
    // --- State Variables ---
    SwStateVariables& state;
    
    // --- Physical Parameters ---
    Scalar dt;                  // Time step
    Scalar grav;                // Gravitational acceleration
    Scalar manning_n;           // Manning's roughness coefficient
    Scalar visc_x;              // Horizontal eddy viscosity (x-direction)
    Scalar visc_y;              // Horizontal eddy viscosity (y-direction)
    Scalar min_depth;           // Minimum depth threshold
    Scalar water_threshold;     // Water depth threshold for wetting (wtfh)
    
    // --- Wind Forcing Parameters ---
    bool sim_wind_;             // Enable wind forcing
    Scalar wind_speed_;         // Current wind speed (m/s)
    Scalar wind_direction_;     // Current wind direction from north (degrees)
    Scalar north_angle_;        // Angle of domain north from true north (degrees)
    Scalar rho_air_;            // Air density (kg/m³), default 1.225
    Scalar rho_water_;          // Water density (kg/m³), default 1000
    Scalar Cw_;                 // Wind drag coefficient, default 0.0013
    Scalar CwT_;                // Thin layer model coefficient, default 1.0
    Scalar hD_;                 // Reference depth for thin layer model (m), default 1.0
    
    // --- Linear Solver ---
    std::unique_ptr<PCGSolver2D> pcg_solver;
    PreconditionerType precond_type;
    Scalar solver_tolerance;
    Ordinal max_solver_iterations;
    
    // --- Solution Vector (active cells only) ---
    View1D<Scalar> solution;    // Solution vector for PCG (eta_new in active space)
    
    // --- Boundary Conditions ---
    SwBoundaryConditionManager* bc_manager_;  // Pointer to boundary condition manager
    
    // --- Source/Sink Terms ---
    SwSourceSinkManager* source_sink_manager_;  // Pointer to source/sink manager
    
    // Constructor
    ShallowWaterSolver(SwDomain& _domain,
                      ActiveCellMesh& _active_mesh,
                      SwStateVariables& _state,
                      Scalar _dt,
                      Scalar _grav = Constants::g,
                      Scalar _manning_n = 0.03,
                      Scalar _visc_x = 0.0,
                      Scalar _visc_y = 0.0,
                      Scalar _min_depth = 1.0e-6,
                      Scalar _water_threshold = 0.01,
                      PreconditionerType _precond = PreconditionerType::JACOBI,
                      Scalar _solver_tol = 1.0e-8,
                      Ordinal _max_iter = 10000,
                      SwBoundaryConditionManager* _bc_manager = nullptr,
                      SwSourceSinkManager* _source_sink_manager = nullptr)
        : domain(_domain),
          active_mesh(_active_mesh),
          state(_state),
          dt(_dt),
          grav(_grav),
          manning_n(_manning_n),
          visc_x(_visc_x),
          visc_y(_visc_y),
          min_depth(_min_depth),
          water_threshold(_water_threshold),
          sim_wind_(false),
          wind_speed_(0.0),
          wind_direction_(0.0),
          north_angle_(0.0),
          rho_air_(1.225),
          rho_water_(1000.0),
          Cw_(0.0013),
          CwT_(1.0),
          hD_(1.0),
          precond_type(_precond),
          solver_tolerance(_solver_tol),
          max_solver_iterations(_max_iter),
          bc_manager_(_bc_manager),
          source_sink_manager_(_source_sink_manager) {
        
        // Create PCG solver
        pcg_solver = std::make_unique<PCGSolver2D>(
            active_mesh.num_active, precond_type, solver_tolerance, max_solver_iterations);
        
        // Allocate solution vector (active cells only)
        solution = View1D<Scalar>("solution", active_mesh.num_active);
    }
    
    // Set boundary condition manager
    void set_boundary_condition_manager(SwBoundaryConditionManager* bc_manager) {
        bc_manager_ = bc_manager;
    }
    
    // Set source/sink manager
    void set_source_sink_manager(SwSourceSinkManager* source_sink_manager) {
        source_sink_manager_ = source_sink_manager;
    }
    
    // ========================================================================
    // WIND FORCING SETUP
    // ========================================================================
    // Enable/disable wind forcing
    void set_wind_enabled(bool enabled) {
        sim_wind_ = enabled;
    }
    
    // Set wind parameters
    void set_wind_parameters(Scalar rho_air, Scalar rho_water, Scalar Cw, 
                             Scalar CwT = 1.0, Scalar hD = 1.0, Scalar north_angle = 0.0) {
        rho_air_ = rho_air;
        rho_water_ = rho_water;
        Cw_ = Cw;
        CwT_ = CwT;
        hD_ = hD;
        north_angle_ = north_angle;
    }
    
    // Update wind conditions (call each time step if wind varies)
    void set_wind_conditions(Scalar wind_speed, Scalar wind_direction) {
        wind_speed_ = wind_speed;
        wind_direction_ = wind_direction;
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
    // ADAPTIVE TIME STEPPING (Based on CFL Condition)
    // ========================================================================
    // Computes recommended time step based on CFL number
    // Returns the new recommended time step
    Scalar adaptive_time_step(Scalar dt_min, Scalar dt_max, Scalar cfl_max = 0.7) {
        Scalar max_cfl = 0.0;
        Scalar dt_new = dt;
        
        auto _velocity_x = state.velocity_x;
        auto _velocity_y = state.velocity_y;
        
        // Compute CFL numbers from current velocities
        Kokkos::parallel_reduce(
            RangePolicy(0, domain.num_cells_total),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i, Scalar& local_max_cfl) {
                Scalar cfl_x_val = std::abs(_velocity_x(i) * dt / domain.dx);
                Scalar cfl_y_val = std::abs(_velocity_y(i) * dt / domain.dy);
                Scalar max_local = (cfl_x_val > cfl_y_val) ? cfl_x_val : cfl_y_val;
                if (max_local > local_max_cfl) local_max_cfl = max_local;
            },
            Kokkos::Max<Scalar>(max_cfl)
        );
        
        // Adjust time step based on CFL
        if (max_cfl > cfl_max && max_cfl > 0.0) {
            // Reduce time step to bring CFL below threshold
            dt_new = dt * cfl_max / max_cfl;
        } else if (max_cfl < 0.5 * cfl_max && max_cfl > 0.0) {
            // Increase time step (but be conservative)
            dt_new = dt * 1.1;  // 10% increase
        }
        
        // Clamp to bounds
        if (dt_new > dt_max) dt_new = dt_max;
        if (dt_new < dt_min) dt_new = dt_min;
        
        return dt_new;
    }
    
    // ========================================================================
    // MAIN SOLVER: Solve for surface elevation (eta)
    // ========================================================================
    void solve(Scalar current_time = 0.0) {
        // Step 1: Save old values
        state.update_old_values();
        
        // Step 2: Compute momentum source terms (Ex, Ey, Dx, Dy)
        compute_momentum_source();
        
        // Step 3: Compute matrix RHS
        compute_rhs();
        
        // Step 4: Add source/sink terms (inflow, rainfall, etc.) to RHS
        add_source_sink_to_rhs(current_time);
        
        // Step 5: Compute matrix coefficients
        compute_matrix_coefficients();
        
        // Step 6: Apply boundary conditions to matrix
        apply_boundary_conditions_to_matrix(current_time);
        
        // Step 7: Solve linear system
        solve_linear_system();
        
        // Step 8: Map solution back to domain
        map_solution_to_domain();
        
        // Step 9: Enforce boundary conditions on solution
        enforce_boundary_conditions_on_solution(current_time);
        
        // Step 10: CFL limiter
        apply_cfl_limiter();
        
        // Step 11: Apply source/sink terms (for post-processing effects like dilution)
        apply_source_sink_terms(current_time);
        
        // Step 12: Update depth
        update_depth();
    }
    
    // ========================================================================
    // VELOCITY UPDATE: Update velocities from pressure gradient
    // ========================================================================
    void update_velocity() {
        // Step 1: Update depth (if needed)
        update_depth();
        
        // Step 2: Update geometry (volumes and areas)
        update_geometry();
        
        // Step 3: Update drag coefficients
        update_drag_coefficients();
        
        // Step 4: Update velocity from pressure gradient
        compute_velocity_from_pressure();
        
        // Step 5: Apply velocity limiters
        apply_velocity_limiters();
        
        // Step 6: Update fluxes and CFL numbers
        update_fluxes_and_cfl();
        
        // Step 7: Interpolate velocities to cell centers
        interpolate_velocities();
    }
    
private:
    // ========================================================================
    // MOMENTUM SOURCE TERMS
    // ========================================================================
    void compute_momentum_source() {
        auto _velocity_x = state.velocity_x;
        auto _velocity_y = state.velocity_y;
        auto _velocity_x_at_y = state.velocity_x_at_y;
        auto _velocity_y_at_x = state.velocity_y_at_x;
        auto _momentum_x = state.momentum_x;
        auto _momentum_y = state.momentum_y;
        auto _drag_factor_x = state.drag_factor_x;
        auto _drag_factor_y = state.drag_factor_y;
        auto _drag_coef_x = state.drag_coef_x;
        auto _drag_coef_y = state.drag_coef_y;
        auto _volume_x = state.volume_x;
        auto _volume_y = state.volume_y;
        auto _area_x = state.area_x;
        auto _area_y = state.area_y;
        auto _area_top_x = state.area_top_x;
        auto _area_top_y = state.area_top_y;
        auto _cfl_x = state.cfl_x;
        auto _cfl_y = state.cfl_y;
        auto _depth = state.depth;
        auto _depth_x = state.depth_x;
        auto _depth_y = state.depth_y;
        auto _pressure = state.pressure;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                // Get domain index
                Ordinal domain_idx = _active_to_domain(i);
                
                // Get neighbor active indices
                Ordinal active_left = _neighbor_left(i);
                Ordinal active_right = _neighbor_right(i);
                Ordinal active_back = _neighbor_back(i);
                Ordinal active_front = _neighbor_front(i);
                
                // Get neighbor domain indices (neighbors are stored as active indices)
                Ordinal i_left = (active_left >= 0) ? _active_to_domain(active_left) : -1;
                Ordinal i_right = (active_right >= 0) ? _active_to_domain(active_right) : -1;
                Ordinal i_back = (active_back >= 0) ? _active_to_domain(active_back) : -1;
                Ordinal i_front = (active_front >= 0) ? _active_to_domain(active_front) : -1;
                
                // Advection terms (upwind scheme)
                Scalar adv_x = 0.0;
                Scalar adv_y = 0.0;
                
                // X-direction advection
                if (active_left >= 0 && active_right >= 0) {
                    Scalar u_here = _velocity_x(domain_idx);
                    Scalar u_left = (i_left >= 0) ? _velocity_x(i_left) : u_here;
                    Scalar u_right = (i_right >= 0) ? _velocity_x(i_right) : u_here;
                    Scalar v_interp = _velocity_y_at_x(domain_idx);
                    Scalar u_back = (i_back >= 0) ? _velocity_x(i_back) : u_here;
                    Scalar u_front = (i_front >= 0) ? _velocity_x(i_front) : u_here;
                    
                    adv_x = (0.5 / domain.dx) * 
                        ((u_here + std::abs(u_here)) * (u_here - u_left) +
                         (u_here - std::abs(u_here)) * (u_right - u_here)) +
                        (0.5 / domain.dy) *
                        ((v_interp + std::abs(v_interp)) * (u_here - u_back) +
                         (v_interp - std::abs(v_interp)) * (u_front - u_here));
                }
                
                // Y-direction advection
                if (active_back >= 0 && active_front >= 0) {
                    Scalar v_here = _velocity_y(domain_idx);
                    Scalar v_back = (i_back >= 0) ? _velocity_y(i_back) : v_here;
                    Scalar v_front = (i_front >= 0) ? _velocity_y(i_front) : v_here;
                    Scalar u_interp = _velocity_x_at_y(domain_idx);
                    Scalar v_left = (i_left >= 0) ? _velocity_y(i_left) : v_here;
                    Scalar v_right = (i_right >= 0) ? _velocity_y(i_right) : v_here;
                    
                    adv_y = (0.5 / domain.dx) *
                        ((u_interp + std::abs(u_interp)) * (v_here - v_left) +
                         (u_interp - std::abs(u_interp)) * (v_right - v_here)) +
                        (0.5 / domain.dy) *
                        ((v_here + std::abs(v_here)) * (v_here - v_back) +
                         (v_here - std::abs(v_here)) * (v_front - v_here));
                }
                
                // CFL-based advection limiter
                if (_velocity_x(domain_idx) == 0.0 || _cfl_x(domain_idx) > 0.7) {
                    adv_x = 0.0;
                } else if (_cfl_x(domain_idx) > 0.5) {
                    adv_x *= (0.7 - _cfl_x(domain_idx)) / (0.7 - 0.5);
                }
                
                if (_velocity_y(domain_idx) == 0.0 || _cfl_y(domain_idx) > 0.7) {
                    adv_y = 0.0;
                } else if (_cfl_y(domain_idx) > 0.5) {
                    adv_y *= (0.7 - _cfl_y(domain_idx)) / (0.7 - 0.5);
                }
                
                // Diffusion terms
                Scalar diff_x = 0.0;
                Scalar diff_y = 0.0;
                
                if (_volume_x(domain_idx) > 0.0 && visc_x > 0.0) {
                    if (active_left >= 0 && active_right >= 0) {
                        diff_x = (visc_x / _volume_x(domain_idx) / domain.dx) *
                            (_area_x(domain_idx) * (_velocity_x(i_right) - _velocity_x(domain_idx)) -
                             _area_x(domain_idx) * (_velocity_x(domain_idx) - _velocity_x(i_left)));
                    }
                    if (active_back >= 0 && active_front >= 0) {
                        Scalar area_y_back = (i_back >= 0) ? _area_y(i_back) : 0.0;
                        diff_x += (visc_y / _volume_x(domain_idx) / domain.dy) *
                            (_area_y(domain_idx) * (_velocity_x(i_front) - _velocity_x(domain_idx)) -
                             area_y_back * (_velocity_x(domain_idx) - _velocity_x(i_back)));
                    }
                }
                
                if (_volume_y(domain_idx) > 0.0 && visc_y > 0.0) {
                    if (active_left >= 0 && active_right >= 0) {
                        Scalar area_x_left = (i_left >= 0) ? _area_x(i_left) : 0.0;
                        diff_y = (visc_x / _volume_y(domain_idx) / domain.dx) *
                            (_area_x(domain_idx) * (_velocity_y(i_right) - _velocity_y(domain_idx)) -
                             area_x_left * (_velocity_y(domain_idx) - _velocity_y(i_left)));
                    }
                    if (active_back >= 0 && active_front >= 0) {
                        diff_y += (visc_y / _volume_y(domain_idx) / domain.dy) *
                            (_area_y(domain_idx) * (_velocity_y(i_front) - _velocity_y(domain_idx)) -
                             _area_y(domain_idx) * (_velocity_y(domain_idx) - _velocity_y(i_back)));
                    }
                }
                
                // Drag terms
                Scalar fac_dx = 0.0;
                Scalar fac_dy = 0.0;
                Scalar vel_mag_x = std::sqrt(_velocity_x(domain_idx) * _velocity_x(domain_idx) + 
                                            _velocity_y_at_x(domain_idx) * _velocity_y_at_x(domain_idx));
                Scalar vel_mag_y = std::sqrt(_velocity_x_at_y(domain_idx) * _velocity_x_at_y(domain_idx) + 
                                            _velocity_y(domain_idx) * _velocity_y(domain_idx));
                
                if (_volume_x(domain_idx) > 0.0) {
                    fac_dx = _area_top_x(domain_idx) / _volume_x(domain_idx);
                }
                if (_volume_y(domain_idx) > 0.0) {
                    fac_dy = _area_top_y(domain_idx) / _volume_y(domain_idx);
                }
                
                // Standard momentum source
                _drag_factor_x(domain_idx) = 1.0 / (0.5 * dt * _drag_coef_x(domain_idx) * vel_mag_x * fac_dx + 1.0);
                _drag_factor_y(domain_idx) = 1.0 / (0.5 * dt * _drag_coef_y(domain_idx) * vel_mag_y * fac_dy + 1.0);
                
                _momentum_x(domain_idx) = _velocity_x(domain_idx) + dt * (diff_x - adv_x);
                _momentum_y(domain_idx) = _velocity_y(domain_idx) + dt * (diff_y - adv_y);
                
                // Apply drag factor
                _momentum_x(domain_idx) *= _drag_factor_x(domain_idx);
                _momentum_y(domain_idx) *= _drag_factor_y(domain_idx);
            });
        
        // Add wind forcing if enabled
        if (sim_wind_) {
            apply_wind_forcing();
        }
    }
    
    // ========================================================================
    // WIND FORCING
    // ========================================================================
    // Applies wind stress to momentum source terms
    // Based on: τ = ρ_air * Cw * (U_wind - u·cos(ω) - v·sin(ω))²
    // With thin-layer model for shallow water
    void apply_wind_forcing() {
        auto _momentum_x = state.momentum_x;
        auto _momentum_y = state.momentum_y;
        auto _velocity_x = state.velocity_x;
        auto _velocity_y = state.velocity_y;
        auto _depth_x = state.depth_x;
        auto _depth_y = state.depth_y;
        auto _drag_factor_x = state.drag_factor_x;
        auto _drag_factor_y = state.drag_factor_y;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        
        // Wind parameters (captured by value for device)
        const Scalar wind_spd = wind_speed_;
        const Scalar wind_dir = wind_direction_;
        const Scalar north_ang = north_angle_;
        const Scalar rho_a = rho_air_;
        const Scalar rho_w = rho_water_;
        const Scalar cw = Cw_;
        const Scalar cwt = CwT_;
        const Scalar hd = hD_;
        const Scalar dt_local = dt;
        const Scalar pi = 3.14159265358979323846;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                
                // Wind direction: phi from north, omega in radians from positive x
                Scalar phi = wind_dir + north_ang;
                Scalar omega = phi * pi / 180.0;
                
                Scalar u = _velocity_x(domain_idx);
                Scalar v = _velocity_y(domain_idx);
                
                // Relative wind speed (wind relative to water surface velocity)
                Scalar u_rel = wind_spd - u * std::cos(omega) - v * std::sin(omega);
                
                // Total wind stress magnitude: τ = ρ_air * Cw * u_rel²
                Scalar tau = rho_a * cw * u_rel * u_rel;
                
                // X-direction wind stress with thin-layer model
                Scalar depth_x = _depth_x(domain_idx);
                Scalar tau_x = tau;
                if (depth_x < hd) {
                    // Thin layer model: reduce wind stress exponentially for shallow water
                    tau_x = tau * std::exp(cwt * (depth_x - hd) / hd);
                    if (depth_x < 0.5 * hd) {
                        tau_x = 0.0;  // No wind stress for very shallow water
                    }
                }
                
                // Y-direction wind stress with thin-layer model
                Scalar depth_y = _depth_y(domain_idx);
                Scalar tau_y = tau;
                if (depth_y < hd) {
                    tau_y = tau * std::exp(cwt * (depth_y - hd) / hd);
                    if (depth_y < 0.5 * hd) {
                        tau_y = 0.0;
                    }
                }
                
                // Add wind stress to momentum
                // Wind stress contribution: dv/dt = τ * cos(ω) / (depth * ρ_water)
                if (depth_x > 0.0) {
                    Scalar wind_mom_x = dt_local * tau_x * std::cos(omega) / (depth_x * rho_w);
                    _momentum_x(domain_idx) += wind_mom_x * _drag_factor_x(domain_idx);
                }
                
                if (depth_y > 0.0) {
                    Scalar wind_mom_y = dt_local * tau_y * std::sin(omega) / (depth_y * rho_w);
                    _momentum_y(domain_idx) += wind_mom_y * _drag_factor_y(domain_idx);
                }
            });
    }
    
    // ========================================================================
    // COMPUTE RHS
    // ========================================================================
    void compute_rhs() {
        auto _pressure = state.pressure;
        auto _momentum_x = state.momentum_x;
        auto _momentum_y = state.momentum_y;
        auto _area_x = state.area_x;
        auto _area_y = state.area_y;
        auto _area_top = state.area_top;
        auto _matrix_rhs = state.matrix_rhs;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                
                // Get neighbors
                Ordinal active_left = _neighbor_left(i);
                Ordinal active_back = _neighbor_back(i);
                Ordinal i_left = (active_left >= 0) ? _active_to_domain(active_left) : -1;
                Ordinal i_back = (active_back >= 0) ? _active_to_domain(active_back) : -1;
                
                // RHS = eta * Asz - dt * div(Ex, Ey)
                Scalar div_flux = _area_x(domain_idx) * _momentum_x(domain_idx);
                if (i_left >= 0) {
                    div_flux -= _area_x(i_left) * _momentum_x(i_left);
                }
                
                div_flux += _area_y(domain_idx) * _momentum_y(domain_idx);
                if (i_back >= 0) {
                    div_flux -= _area_y(i_back) * _momentum_y(i_back);
                }
                
                _matrix_rhs(domain_idx) = _pressure(domain_idx) * _area_top(domain_idx) - dt * div_flux;
            });
    }
    
    // ========================================================================
    // ADD SOURCE/SINK TERMS TO RHS (Inflow, Rainfall, Evaporation, etc.)
    // ========================================================================
    // Source/sink terms are added to the RHS as: RHS += source_value * dt
    // - Positive values: add water (inflow, rainfall)
    // - Negative values: remove water (outflow, evaporation)
    void add_source_sink_to_rhs(Scalar current_time) {
        if (!source_sink_manager_) return;
        
        auto _matrix_rhs = state.matrix_rhs;
        auto _area_top = state.area_top;
        
        const auto& ss_terms = source_sink_manager_->get_source_sink_terms();
        
        // Host mirrors for source/sink term application
        auto h_matrix_rhs = Kokkos::create_mirror_view(_matrix_rhs);
        auto h_area_top = Kokkos::create_mirror_view(_area_top);
        
        Kokkos::deep_copy(h_matrix_rhs, _matrix_rhs);
        Kokkos::deep_copy(h_area_top, _area_top);
        
        for (const auto& ss : ss_terms) {
            Scalar ss_value = ss.get_value(current_time);
            Ordinal n_cells = ss.cell_indices.size();
            
            if (n_cells == 0) continue;
            
            for (Ordinal cell_idx : ss.cell_indices) {
                if (ss.type == SourceSinkType::VOLUME_FLUX) {
                    // Volume flux (m³/s) - distribute among cells
                    h_matrix_rhs(cell_idx) += ss_value * dt / static_cast<Scalar>(n_cells);
                    
                } else if (ss.type == SourceSinkType::DEPTH_RATE) {
                    // Depth rate (m/s) - convert to volume using cell area
                    h_matrix_rhs(cell_idx) += ss_value * dt * h_area_top(cell_idx);
                }
            }
        }
        
        Kokkos::deep_copy(_matrix_rhs, h_matrix_rhs);
    }
    
    // ========================================================================
    // COMPUTE MATRIX COEFFICIENTS
    // ========================================================================
    void compute_matrix_coefficients() {
        auto _matrix_diag = state.matrix_diag;
        auto _matrix_xp = state.matrix_xp;
        auto _matrix_xm = state.matrix_xm;
        auto _matrix_yp = state.matrix_yp;
        auto _matrix_ym = state.matrix_ym;
        auto _matrix_rhs = state.matrix_rhs;
        auto _drag_factor_x = state.drag_factor_x;
        auto _drag_factor_y = state.drag_factor_y;
        auto _area_x = state.area_x;
        auto _area_y = state.area_y;
        auto _area_top = state.area_top;
        auto _volume_x = state.volume_x;
        auto _volume_y = state.volume_y;
        auto _depth = state.depth;
        auto _velocity_x = state.velocity_x;
        auto _velocity_y = state.velocity_y;
        auto _pressure = state.pressure;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_back = active_mesh.neighbor_back;
        
        Scalar coef = grav * dt * dt;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                
                Ordinal active_left = _neighbor_left(i);
                Ordinal active_back = _neighbor_back(i);
                Ordinal i_left = (active_left >= 0) ? _active_to_domain(active_left) : -1;
                Ordinal i_back = (active_back >= 0) ? _active_to_domain(active_back) : -1;
                
                // Initialize coefficients
                _matrix_xp(domain_idx) = 0.0;
                _matrix_xm(domain_idx) = 0.0;
                _matrix_yp(domain_idx) = 0.0;
                _matrix_ym(domain_idx) = 0.0;
                
                // X+ coefficient
                if (_volume_x(domain_idx) > 0.0) {
                    _matrix_xp(domain_idx) = coef * _area_x(domain_idx) * _area_x(domain_idx) * 
                                            _drag_factor_x(domain_idx) / _volume_x(domain_idx);
                }
                
                // X- coefficient (from left neighbor)
                if (i_left >= 0 && _volume_x(i_left) > 0.0) {
                    _matrix_xm(domain_idx) = coef * _area_x(i_left) * _area_x(i_left) * 
                                            _drag_factor_x(i_left) / _volume_x(i_left);
                }
                
                // Y+ coefficient
                if (_volume_y(domain_idx) > 0.0) {
                    _matrix_yp(domain_idx) = coef * _area_y(domain_idx) * _area_y(domain_idx) * 
                                            _drag_factor_y(domain_idx) / _volume_y(domain_idx);
                }
                
                // Y- coefficient (from back neighbor)
                if (i_back >= 0 && _volume_y(i_back) > 0.0) {
                    _matrix_ym(domain_idx) = coef * _area_y(i_back) * _area_y(i_back) * 
                                            _drag_factor_y(i_back) / _volume_y(i_back);
                }
                
                // Diagonal
                _matrix_diag(domain_idx) = _area_top(domain_idx) + 
                                          _matrix_xp(domain_idx) + _matrix_xm(domain_idx) +
                                          _matrix_yp(domain_idx) + _matrix_ym(domain_idx);
                
                // Avoid singularity
                if (_depth(domain_idx) == 0.0) {
                    _matrix_diag(domain_idx) = domain.dx * domain.dy;
                    _matrix_rhs(domain_idx) = _pressure(domain_idx) * domain.dx * domain.dy;
                    
                    if (_velocity_x(domain_idx) == 0.0 && 
                        (i_left < 0 || _velocity_x(i_left) == 0.0) &&
                        _velocity_y(domain_idx) == 0.0 &&
                        (i_back < 0 || _velocity_y(i_back) == 0.0)) {
                        _matrix_xp(domain_idx) = 0.0;
                        _matrix_xm(domain_idx) = 0.0;
                        _matrix_yp(domain_idx) = 0.0;
                        _matrix_ym(domain_idx) = 0.0;
                    }
                } else if (_matrix_diag(domain_idx) == 0.0) {
                    _matrix_diag(domain_idx) = domain.dx * domain.dy;
                    _matrix_rhs(domain_idx) = _pressure(domain_idx) * domain.dx * domain.dy;
                }
                
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
        auto _matrix_rhs = state.matrix_rhs;
        auto _bottom = state.bottom;
        
        const auto& bcs = bc_manager_->get_boundary_conditions();
        Scalar dt_local = dt;
        
        // Process each boundary condition using device-side parallel operations
        for (const auto& bc : bcs) {
            Scalar bc_value = bc.get_value(current_time);
            Ordinal n_bc_cells = bc.cell_indices.size();
            
            if (n_bc_cells == 0) continue;
            
            // Create device view for BC cell indices (only for this BC)
            View1D<Ordinal> d_bc_indices("bc_indices", n_bc_cells);
            auto h_bc_indices = Kokkos::create_mirror_view(d_bc_indices);
            for (Ordinal i = 0; i < n_bc_cells; ++i) {
                h_bc_indices(i) = bc.cell_indices[i];
            }
            Kokkos::deep_copy(d_bc_indices, h_bc_indices);
            
            if (bc.type == SwBcType::FREE_SURFACE_ELEVATION) {
                // Dirichlet BC: prescribed surface elevation
                Kokkos::parallel_for(RangePolicy(0, n_bc_cells),
                    [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                        Ordinal cell_idx = d_bc_indices(i);
                        _matrix_diag(cell_idx) = 1.0;
                        _matrix_xp(cell_idx) = 0.0;
                        _matrix_xm(cell_idx) = 0.0;
                        _matrix_yp(cell_idx) = 0.0;
                        _matrix_ym(cell_idx) = 0.0;
                        _matrix_rhs(cell_idx) = bc_value;
                    });
                    
            } else if (bc.type == SwBcType::WATER_DEPTH) {
                // Dirichlet BC: prescribed water depth
                Kokkos::parallel_for(RangePolicy(0, n_bc_cells),
                    [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                        Ordinal cell_idx = d_bc_indices(i);
                        Scalar eta_value = _bottom(cell_idx) + bc_value;
                        _matrix_diag(cell_idx) = 1.0;
                        _matrix_xp(cell_idx) = 0.0;
                        _matrix_xm(cell_idx) = 0.0;
                        _matrix_yp(cell_idx) = 0.0;
                        _matrix_ym(cell_idx) = 0.0;
                        _matrix_rhs(cell_idx) = eta_value;
                    });
                    
            } else if (bc.type == SwBcType::FLOW_RATE) {
                // Flow rate BC: add to RHS as source term
                Kokkos::parallel_for(RangePolicy(0, n_bc_cells),
                    [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                        Ordinal cell_idx = d_bc_indices(i);
                        Scalar source_flux = bc_value * dt_local;
                        _matrix_rhs(cell_idx) += source_flux;
                    });
                    
            } else if (bc.type == SwBcType::FREE_OUTFLOW) {
                // Free outflow: zero gradient (Neumann BC)
                // Natural outflow - no modification needed
            }
        }
    }
    
    // ========================================================================
    // ENFORCE BOUNDARY CONDITIONS ON SOLUTION
    // ========================================================================
    void enforce_boundary_conditions_on_solution(Scalar current_time) {
        if (!bc_manager_) return;
        
        auto _pressure = state.pressure;
        auto _bottom = state.bottom;
        
        const auto& bcs = bc_manager_->get_boundary_conditions();
        
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
            
            if (bc.type == SwBcType::FREE_SURFACE_ELEVATION) {
                Kokkos::parallel_for(RangePolicy(0, n_bc_cells),
                    [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                        Ordinal cell_idx = d_bc_indices(i);
                        _pressure(cell_idx) = bc_value;
                    });
                    
            } else if (bc.type == SwBcType::WATER_DEPTH) {
                Kokkos::parallel_for(RangePolicy(0, n_bc_cells),
                    [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                        Ordinal cell_idx = d_bc_indices(i);
                        _pressure(cell_idx) = _bottom(cell_idx) + bc_value;
                    });
                    
            } else if (bc.type == SwBcType::FREE_OUTFLOW) {
                Kokkos::parallel_for(RangePolicy(0, n_bc_cells),
                    [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                        Ordinal cell_idx = d_bc_indices(i);
                        if (_pressure(cell_idx) < _bottom(cell_idx)) {
                            _pressure(cell_idx) = _bottom(cell_idx);
                        }
                    });
            }
        }
    }
    
    // ========================================================================
    // APPLY SOURCE/SINK TERMS
    // ========================================================================
    void apply_source_sink_terms(Scalar current_time) {
        if (!source_sink_manager_) return;
        
        auto _pressure = state.pressure;
        auto _area_top = state.area_top;
        
        const auto& ss_terms = source_sink_manager_->get_source_sink_terms();
        Scalar dt_local = dt;
        Scalar dx_local = domain.dx;
        Scalar dy_local = domain.dy;
        
        // Process each source/sink term using device-side parallel operations
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
            
            if (ss.type == SourceSinkType::VOLUME_FLUX) {
                // Volume flux (m³/s): add/remove volume, convert to depth change
                Kokkos::parallel_for(RangePolicy(0, n_ss_cells),
                    [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                        Ordinal cell_idx = d_ss_indices(i);
                        Scalar cell_area = _area_top(cell_idx);
                        if (cell_area <= 0.0) cell_area = dx_local * dy_local;
                        Scalar volume_change = ss_value * dt_local;
                        Scalar depth_change = volume_change / cell_area;
                        _pressure(cell_idx) += depth_change;
                    });
                    
            } else if (ss.type == SourceSinkType::DEPTH_RATE) {
                // Depth rate (m/s): directly add/remove depth
                Kokkos::parallel_for(RangePolicy(0, n_ss_cells),
                    [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                        Ordinal cell_idx = d_ss_indices(i);
                        Scalar depth_change = ss_value * dt_local;
                        _pressure(cell_idx) += depth_change;
                    });
                    
            } else if (ss.type == SourceSinkType::MASS_FLUX) {
                // Mass flux: for scalar transport only, not applied here
                // This will be handled in scalar transport solver
            }
        }
        
        // Ensure eta >= bottom for all cells (device-side)
        auto _bottom = state.bottom;
        auto _active_to_domain = active_mesh.active_to_domain;
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_pressure(domain_idx) < _bottom(domain_idx)) {
                    _pressure(domain_idx) = _bottom(domain_idx);
                }
            });
    }
    
    // ========================================================================
    // SOLVE LINEAR SYSTEM
    // ========================================================================
    void solve_linear_system() {
        // Create views for active cells only
        View1D<Scalar> diag_active("diag_active", active_mesh.num_active);
        View1D<Scalar> xp_active("xp_active", active_mesh.num_active);
        View1D<Scalar> xm_active("xm_active", active_mesh.num_active);
        View1D<Scalar> yp_active("yp_active", active_mesh.num_active);
        View1D<Scalar> ym_active("ym_active", active_mesh.num_active);
        View1D<Scalar> rhs_active("rhs_active", active_mesh.num_active);
        
        // Map domain variables to active cell space
        auto _matrix_diag = state.matrix_diag;
        auto _matrix_xp = state.matrix_xp;
        auto _matrix_xm = state.matrix_xm;
        auto _matrix_yp = state.matrix_yp;
        auto _matrix_ym = state.matrix_ym;
        auto _matrix_rhs = state.matrix_rhs;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                diag_active(i) = _matrix_diag(domain_idx);
                xp_active(i) = _matrix_xp(domain_idx);
                xm_active(i) = _matrix_xm(domain_idx);
                yp_active(i) = _matrix_yp(domain_idx);
                ym_active(i) = _matrix_ym(domain_idx);
                rhs_active(i) = _matrix_rhs(domain_idx);
            });
        
        // Solve using PCG
        pcg_solver->solve_2d(active_mesh,
                            diag_active, xp_active, xm_active, yp_active, ym_active,
                            rhs_active, solution);
    }
    
    // ========================================================================
    // MAP SOLUTION BACK TO DOMAIN
    // ========================================================================
    void map_solution_to_domain() {
        auto _pressure = state.pressure;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                _pressure(domain_idx) = solution(i);
            });
    }
    
    // ========================================================================
    // APPLY SEEPAGE TO SURFACE ELEVATION
    // ========================================================================
    // Applies seepage/infiltration from subsurface to surface water
    // This should be called after groundwater solver computes seepage_rate
    void apply_seepage() {
        auto _pressure = state.pressure;
        auto _bottom = state.bottom;
        auto _depth = state.depth;
        auto _seepage_rate = state.seepage_rate;
        auto _reset_seepage = state.reset_seepage;
        
        Kokkos::parallel_for(RangePolicy(0, domain.num_cells_total),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Scalar qss = _seepage_rate(i);
                Scalar diff = _pressure(i) - _bottom(i);
                
                // Infiltration (qss < 0): water goes from surface to subsurface
                if (qss < 0.0) {
                    _pressure(i) += qss * dt;
                    _reset_seepage(i) = 1;
                    // Ensure pressure doesn't go below bottom
                    if (_pressure(i) < _bottom(i)) {
                        _pressure(i) = _bottom(i);
                    }
                }
                // Seepage (qss > 0): water comes from subsurface to surface
                else if (qss > 0.0) {
                    // If surface is dry (depth <= min_depth)
                    if (diff <= min_depth) {
                        // Only apply if seepage is significant enough
                        if (qss * dt > min_depth) {
                            _pressure(i) += qss * dt;
                            _reset_seepage(i) = 1;
                        } else {
                            _reset_seepage(i) = 0;
                        }
                    } else {
                        // Surface is wet: always apply seepage
                        _pressure(i) += qss * dt;
                        _reset_seepage(i) = 1;
                    }
                }
                // No seepage (qss == 0): do nothing
                
                // Update depth
                Scalar new_depth = _pressure(i) - _bottom(i);
                if (new_depth < 0.0) {
                    _pressure(i) = _bottom(i);
                    _depth(i) = 0.0;
                } else if (new_depth < min_depth) {
                    _pressure(i) = _bottom(i);
                    _depth(i) = 0.0;
                } else {
                    _depth(i) = new_depth;
                }
            });
    }
    
    // ========================================================================
    // UPDATE DEPTH
    // ========================================================================
    void update_depth() {
        auto _pressure = state.pressure;
        auto _bottom = state.bottom;
        auto _depth = state.depth;
        auto _depth_x = state.depth_x;
        auto _depth_y = state.depth_y;
        
        Kokkos::parallel_for(RangePolicy(0, domain.num_cells_total),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Scalar h = _pressure(i) - _bottom(i);
                _depth(i) = (h > 0.0) ? h : 0.0;
                
                // Update face depths (simplified - would use neighbor info in full implementation)
                _depth_x(i) = _depth(i);
                _depth_y(i) = _depth(i);
            });
    }
    
    // ========================================================================
    // UPDATE GEOMETRY (Volumes and Areas)
    // ========================================================================
    void update_geometry() {
        auto _depth = state.depth;
        auto _area_top = state.area_top;
        auto _area_x = state.area_x;
        auto _area_y = state.area_y;
        auto _volume = state.volume;
        auto _volume_x = state.volume_x;
        auto _volume_y = state.volume_y;
        
        Kokkos::parallel_for(RangePolicy(0, domain.num_cells_total),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                // Top area
                if (_depth(i) > 0.0) {
                    _area_top(i) = domain.dx * domain.dy;
                } else {
                    _area_top(i) = 0.0;
                }
                
                // Face areas
                _area_x(i) = _depth(i) * domain.dy;
                _area_y(i) = _depth(i) * domain.dx;
                
                // Volumes
                _volume(i) = _depth(i) * domain.dx * domain.dy;
                
                // Staggered volumes (simplified - would average with neighbors)
                _volume_x(i) = _volume(i);
                _volume_y(i) = _volume(i);
            });
    }
    
    // ========================================================================
    // UPDATE DRAG COEFFICIENTS
    // ========================================================================
    void update_drag_coefficients() {
        auto _depth = state.depth;
        auto _volume = state.volume;
        auto _drag_coef_x = state.drag_coef_x;
        auto _drag_coef_y = state.drag_coef_y;
        
        Scalar coef = grav * manning_n * manning_n;
        
        Kokkos::parallel_for(RangePolicy(0, domain.num_cells_total),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                if (_volume(i) > 0.0) {
                    Scalar eff_depth = _volume(i) / (domain.dx * domain.dy);
                    Scalar expo = (eff_depth < 0.1) ? (2.0 / 3.0) : (1.0 / 3.0);  // Thin layer model
                    Scalar cd = coef / std::pow(eff_depth, expo);
                    _drag_coef_x(i) = cd;
                    _drag_coef_y(i) = cd;
                } else {
                    _drag_coef_x(i) = 0.0;
                    _drag_coef_y(i) = 0.0;
                }
            });
    }
    
    // ========================================================================
    // COMPUTE VELOCITY FROM PRESSURE GRADIENT
    // ========================================================================
    void compute_velocity_from_pressure() {
        auto _pressure = state.pressure;
        auto _momentum_x = state.momentum_x;
        auto _momentum_y = state.momentum_y;
        auto _drag_factor_x = state.drag_factor_x;
        auto _drag_factor_y = state.drag_factor_y;
        auto _velocity_x = state.velocity_x;
        auto _velocity_y = state.velocity_y;
        auto _volume_x = state.volume_x;
        auto _volume_y = state.volume_y;
        auto _area_x = state.area_x;
        auto _area_y = state.area_y;
        auto _depth = state.depth;
        
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_front = active_mesh.neighbor_front;
        
        Scalar coef = grav * dt;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                
                Ordinal active_right = _neighbor_right(i);
                Ordinal active_front = _neighbor_front(i);
                Ordinal i_right = (active_right >= 0) ? _active_to_domain(active_right) : -1;
                Ordinal i_front = (active_front >= 0) ? _active_to_domain(active_front) : -1;
                
                Scalar eff_hx = 0.0;
                Scalar eff_hy = 0.0;
                
                if (_volume_x(domain_idx) > 0.0) {
                    eff_hx = _area_x(domain_idx) / _volume_x(domain_idx);
                }
                if (_volume_y(domain_idx) > 0.0) {
                    eff_hy = _area_y(domain_idx) / _volume_y(domain_idx);
                }
                
                // Standard velocity update
                if (i_right >= 0) {
                    _velocity_x(domain_idx) = (_momentum_x(domain_idx) - 
                                               coef * eff_hx * (_pressure(i_right) - _pressure(domain_idx))) *
                                              _drag_factor_x(domain_idx);
                } else {
                    _velocity_x(domain_idx) = _momentum_x(domain_idx) * _drag_factor_x(domain_idx);
                }
                
                if (i_front >= 0) {
                    _velocity_y(domain_idx) = (_momentum_y(domain_idx) - 
                                              coef * eff_hy * (_pressure(i_front) - _pressure(domain_idx))) *
                                             _drag_factor_y(domain_idx);
                } else {
                    _velocity_y(domain_idx) = _momentum_y(domain_idx) * _drag_factor_y(domain_idx);
                }
            });
    }
    
    // ========================================================================
    // APPLY VELOCITY LIMITERS
    // ========================================================================
    void apply_velocity_limiters() {
        auto _velocity_x = state.velocity_x;
        auto _velocity_y = state.velocity_y;
        auto _area_x = state.area_x;
        auto _area_y = state.area_y;
        auto _depth = state.depth;
        auto _cfl_active = state.cfl_active;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_back = active_mesh.neighbor_back;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                
                // Zero velocity when face area is too small
                if (_area_x(domain_idx) < water_threshold * domain.dy) {
                    _velocity_x(domain_idx) = 0.0;
                }
                if (_area_y(domain_idx) < water_threshold * domain.dx) {
                    _velocity_y(domain_idx) = 0.0;
                }
                
                // Zero velocity out of dry cells
                if (_depth(domain_idx) < water_threshold) {
                    if (_velocity_x(domain_idx) > 0.0) {
                        _velocity_x(domain_idx) = 0.0;
                    }
                    Ordinal active_left = _neighbor_left(i);
                    if (active_left >= 0) {
                        Ordinal i_left = _active_to_domain(active_left);
                        if (_velocity_x(i_left) < 0.0) {
                            _velocity_x(i_left) = 0.0;
                        }
                    }
                    
                    if (_velocity_y(domain_idx) > 0.0) {
                        _velocity_y(domain_idx) = 0.0;
                    }
                    Ordinal active_back = _neighbor_back(i);
                    if (active_back >= 0) {
                        Ordinal i_back = _active_to_domain(active_back);
                        if (_velocity_y(i_back) < 0.0) {
                            _velocity_y(i_back) = 0.0;
                        }
                    }
                }
                
                // CFL limiter
                if (_cfl_active(domain_idx) == 1) {
                    _velocity_x(domain_idx) = 0.0;
                    Ordinal active_left = _neighbor_left(i);
                    if (active_left >= 0) {
                        Ordinal i_left = _active_to_domain(active_left);
                        _velocity_x(i_left) = 0.0;
                    }
                    _velocity_y(domain_idx) = 0.0;
                    Ordinal active_back = _neighbor_back(i);
                    if (active_back >= 0) {
                        Ordinal i_back = _active_to_domain(active_back);
                        _velocity_y(i_back) = 0.0;
                    }
                    _cfl_active(domain_idx) = 0;
                }
            });
    }
    
    // ========================================================================
    // UPDATE FLUXES AND CFL NUMBERS
    // ========================================================================
    void update_fluxes_and_cfl() {
        auto _velocity_x = state.velocity_x;
        auto _velocity_y = state.velocity_y;
        auto _area_x = state.area_x;
        auto _area_y = state.area_y;
        auto _flux_x = state.flux_x;
        auto _flux_y = state.flux_y;
        auto _cfl_x = state.cfl_x;
        auto _cfl_y = state.cfl_y;
        
        Kokkos::parallel_for(RangePolicy(0, domain.num_cells_total),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                _flux_x(i) = _velocity_x(i) * _area_x(i);
                _flux_y(i) = _velocity_y(i) * _area_y(i);
                _cfl_x(i) = std::abs(_velocity_x(i) * dt / domain.dx);
                _cfl_y(i) = std::abs(_velocity_y(i) * dt / domain.dy);
            });
    }
    
    // ========================================================================
    // INTERPOLATE VELOCITIES
    // ========================================================================
    void interpolate_velocities() {
        auto _velocity_x = state.velocity_x;
        auto _velocity_y = state.velocity_y;
        auto _velocity_x_at_y = state.velocity_x_at_y;
        auto _velocity_y_at_x = state.velocity_y_at_x;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                
                Ordinal active_left = _neighbor_left(i);
                Ordinal active_back = _neighbor_back(i);
                Ordinal active_right = _neighbor_right(i);
                Ordinal active_front = _neighbor_front(i);
                Ordinal i_left = (active_left >= 0) ? _active_to_domain(active_left) : domain_idx;
                Ordinal i_back = (active_back >= 0) ? _active_to_domain(active_back) : domain_idx;
                Ordinal i_right = (active_right >= 0) ? _active_to_domain(active_right) : domain_idx;
                Ordinal i_front = (active_front >= 0) ? _active_to_domain(active_front) : domain_idx;
                
                // Interpolate u to y-face (uy)
                _velocity_x_at_y(domain_idx) = 0.25 * (_velocity_x(domain_idx) + _velocity_x(i_left) +
                                                       _velocity_x(i_front) + 
                                                       (active_front >= 0 && _neighbor_left(active_front) >= 0 ? 
                                                        _velocity_x(_active_to_domain(_neighbor_left(active_front))) : 
                                                        _velocity_x(i_front)));
                
                // Interpolate v to x-face (vx)
                _velocity_y_at_x(domain_idx) = 0.25 * (_velocity_y(domain_idx) + _velocity_y(i_back) +
                                                       _velocity_y(i_right) +
                                                       (active_right >= 0 && _neighbor_back(active_right) >= 0 ?
                                                        _velocity_y(_active_to_domain(_neighbor_back(active_right))) :
                                                        _velocity_y(i_right)));
            });
    }
    
    // ========================================================================
    // APPLY CFL LIMITER
    // ========================================================================
    void apply_cfl_limiter() {
        auto _pressure = state.pressure;
        auto _bottom = state.bottom;
        auto _depth = state.depth;
        auto _cfl_active = state.cfl_active;
        
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                
                Scalar diff = _pressure(domain_idx) - _bottom(domain_idx);
                
                // Check for isolated wet cells
                if (_depth(domain_idx) <= 0.0 && diff > 0.0) {
                    bool has_wet_neighbor = false;
                    
                    Ordinal active_right = _neighbor_right(i);
                    Ordinal active_left = _neighbor_left(i);
                    Ordinal active_front = _neighbor_front(i);
                    Ordinal active_back = _neighbor_back(i);
                    Ordinal i_right = (active_right >= 0) ? _active_to_domain(active_right) : -1;
                    Ordinal i_left = (active_left >= 0) ? _active_to_domain(active_left) : -1;
                    Ordinal i_front = (active_front >= 0) ? _active_to_domain(active_front) : -1;
                    Ordinal i_back = (active_back >= 0) ? _active_to_domain(active_back) : -1;
                    
                    if (i_right >= 0 && _depth(i_right) > 0.0) has_wet_neighbor = true;
                    if (i_left >= 0 && _depth(i_left) > 0.0) has_wet_neighbor = true;
                    if (i_front >= 0 && _depth(i_front) > 0.0) has_wet_neighbor = true;
                    if (i_back >= 0 && _depth(i_back) > 0.0) has_wet_neighbor = true;
                    
                    if (!has_wet_neighbor) {
                        _pressure(domain_idx) = _bottom(domain_idx);
                        _cfl_active(domain_idx) = 1;
                    }
                }
                
                // Remove small depth
                if (diff > 0.0 && diff < min_depth) {
                    _pressure(domain_idx) = _bottom(domain_idx);
                }
            });
    }
};

#endif // FREHG_SHALLOW_WATER_SOLVER_HPP

