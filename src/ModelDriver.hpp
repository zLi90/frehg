#ifndef FREHG_MODEL_DRIVER_HPP
#define FREHG_MODEL_DRIVER_HPP

#include "define.hpp"
#include "Domain.hpp"
#include "StateVariables.hpp"
#include "ActiveCellMesh.hpp"
#include "ShallowWaterSolver.hpp"
#include "GroundwaterSolver.hpp"
#include "ScalarTransportSolver.hpp"
#include "Initializer.hpp"
#include "OutputWriter.hpp"
#include <Kokkos_Core.hpp> // Required for Kokkos::Timer
#include <string>
#include <iostream>
#include <memory> 
#include <vector>
#include <cmath>

namespace Frehg {

class ModelDriver {
private:
    // Performance timer
    Kokkos::Timer timer_;
    
    // Note: Domains, state variables, meshes, and solvers are owned by initializer_
    // Access them via initializer_->get_*() methods
    
    // --- Simulation Parameters ---
    bool sim_shallowwater_;
    bool sim_groundwater_;
    bool sync_coupling_;          // Synchronous coupling (same dt) vs async
    bool use_adaptive_dt_;        // Use adaptive time stepping
    int n_scalar_;                // Number of scalar species
    bool baroclinic_;             // Baroclinic effects (density/viscosity from scalar)
    Scalar dt_;                   // Surface water time step
    Scalar dtg_;                  // Groundwater time step (may differ)
    Scalar dt_min_;               // Minimum time step
    Scalar dt_max_;               // Maximum time step
    Scalar Co_max_;               // Maximum Courant number for groundwater
    Scalar cfl_max_;              // Maximum CFL number for shallow water
    Scalar end_time_;             // Total simulation time
    Scalar dt_out_;               // Output interval
    int output_interval_;          // Output every N steps
    
    // --- Output Directory ---
    std::string output_dir_;
    
    // --- Initializer ---
    std::unique_ptr<Initializer> initializer_;
    
    // --- Output Writer ---
    std::unique_ptr<OutputWriter> output_writer_;

public:
    // ========================================================================
    // Constructor
    // ========================================================================
    ModelDriver() = default;

    // ========================================================================
    // Destructor
    // ========================================================================
    ~ModelDriver() = default;

    // ========================================================================
    // Initialization Phase
    // ========================================================================
    void initialize(const std::string& input_dir, const std::string& output_dir) {
        output_dir_ = output_dir;
        
        std::cout << "------------------------------------------------------------\n";
        std::cout << "[Driver] Phase 1: Initialization\n";
        std::cout << "[Driver] Reading configuration from: " << input_dir << "\n";
        std::cout << "[Driver] Writing results to        : " << output_dir << "\n";
        std::cout << "------------------------------------------------------------\n";

        // Create initializer and perform all initialization steps
        initializer_ = std::make_unique<Initializer>(input_dir, output_dir);
        initializer_->initialize_all();
        
        // Extract configuration from initializer
        const auto& config = initializer_->get_config();
        sim_shallowwater_ = config.sim_shallowwater;
        sim_groundwater_ = config.sim_groundwater;
        sync_coupling_ = config.sync_coupling;
        n_scalar_ = config.n_scalar;
        baroclinic_ = config.baroclinic;
        dt_ = config.dt;
        dtg_ = config.dt;
        
        // Adaptive time stepping parameters (with defaults)
        use_adaptive_dt_ = config.dt_adjust;  // Will be read from input
        dt_min_ = (config.dt_min > 0.0) ? config.dt_min : config.dt * 0.1;
        dt_max_ = (config.dt_max > 0.0) ? config.dt_max : config.dt;
        Co_max_ = (config.Co_max > 0.0) ? config.Co_max : 2.0;
        cfl_max_ = 0.7;  // Default CFL limit for shallow water
        
        end_time_ = config.Tend;
        dt_out_ = config.dt_out;
        output_interval_ = static_cast<int>(config.dt_out / config.dt);
        
        // Create OutputWriter
        output_writer_ = std::make_unique<OutputWriter>(
            output_dir_,
            dt_out_,
            initializer_->get_sw_domain(),
            initializer_->get_gw_domain(),
            initializer_->get_sw_state(),
            initializer_->get_gw_state(),
            initializer_->get_sw_active_mesh(),
            initializer_->get_gw_active_mesh()
        );
        
        // Set scalar solvers if available
        if (n_scalar_ > 0) {
            std::vector<SwScalarTransportSolver*> sw_scalars;
            std::vector<GwScalarTransportSolver*> gw_scalars;
            
            if (sim_shallowwater_) {
                const auto& sw_solvers = initializer_->get_sw_scalar_solvers();
                for (const auto& solver : sw_solvers) {
                    sw_scalars.push_back(solver.get());
                }
            }
            
            if (sim_groundwater_) {
                const auto& gw_solvers = initializer_->get_gw_scalar_solvers();
                for (const auto& solver : gw_solvers) {
                    gw_scalars.push_back(solver.get());
                }
            }
            
            output_writer_->set_scalar_solvers(sw_scalars, gw_scalars);
        }
        
        // Initialize OutputWriter
        output_writer_->initialize(dt_, n_scalar_);
        
        // IMPORTANT: Initialize volume and volume_old for groundwater
        // This is needed for scalar mass calculation (mass = concentration * volume)
        if (sim_groundwater_) {
            initialize_gw_volume();
        }
        
        std::cout << "[Driver] Initialization complete.\n" << std::endl;
    }
    
    // ========================================================================
    // Initialize groundwater volume from water content
    // ========================================================================
    void initialize_gw_volume() {
        auto* gw_state = initializer_->get_gw_state();
        auto* gw_domain = initializer_->get_gw_domain();
        if (!gw_state || !gw_domain) return;
        
        auto _volume = gw_state->volume;
        auto _volume_old = gw_state->volume_old;
        auto _water_content = gw_state->water_content;
        auto _active_mask_3d = gw_domain->active_mask_3d;
        auto _layer_thickness = gw_domain->layer_thickness;
        
        Scalar dx = gw_domain->dx;
        Scalar dy = gw_domain->dy;
        Ordinal nx = gw_domain->nx;
        Ordinal ny = gw_domain->ny;
        Ordinal nz = gw_domain->nz;
        
        using RangePolicy = Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>;
        Kokkos::parallel_for(RangePolicy(0, gw_domain->num_cells_3d_total),
            KOKKOS_LAMBDA (const Ordinal idx) {
                if (_active_mask_3d(idx) > 0) {
                    // Compute layer k from 3D index
                    Ordinal k = idx / (nx * ny);
                    Scalar dz = _layer_thickness(k);
                    Scalar cell_volume = dx * dy * dz;
                    
                    // Volume = water content * cell volume
                    _volume(idx) = _water_content(idx) * cell_volume;
                    _volume_old(idx) = _volume(idx);
                }
            });
        Kokkos::fence();
        
        std::cout << "[Driver] Initialized GW volume from water content" << std::endl;
    }

    // ========================================================================
    // Execution Phase (Time Loop)
    // ========================================================================
    void run() {
        std::cout << "------------------------------------------------------------\n";
        std::cout << "[Driver] Phase 2: Simulation Loop\n";
        std::cout << "------------------------------------------------------------\n";
        
        // Reset timer to measure pure computation time (excluding setup)
        timer_.reset();
        
        // Initialize time stepping variables
        Scalar current_time = 0.0;
        Scalar t_subsurface = 0.0;  // Subsurface time (for async coupling)
        int step = 0;
        
        // Set groundwater time step
        if (!sim_groundwater_ || sync_coupling_) {
            dtg_ = dt_;
        } else {
            dtg_ = dt_min_;
        }
        
        // Write initial conditions
        std::cout << "[Driver] Writing initial conditions...\n";
        if (output_writer_) {
            output_writer_->write_spatial_fields(current_time, step);
            output_writer_->write_time_series(current_time, step);
        }
        
        std::cout << "[Driver] Beginning time loop...\n";

        while (current_time < end_time_) {
            // Adaptive time stepping: adjust dt before the step
            Scalar dt_sw = dt_;
            Scalar dt_gw = dtg_;
            
            if (use_adaptive_dt_) {
                // 1. Compute adaptive time step for shallow water (if enabled)
                if (sim_shallowwater_ && initializer_) {
                    auto* sw_solver = initializer_->get_sw_solver();
                    if (sw_solver) {
                        // Compute adaptive time step based on current velocities
                        dt_sw = sw_solver->adaptive_time_step(dt_min_, dt_max_, cfl_max_);
                        sw_solver->set_time_step(dt_sw);
                    }
                }
                
                // 2. Compute adaptive time step for groundwater (if enabled)
                if (sim_groundwater_ && initializer_) {
                    auto* gw_solver = initializer_->get_gw_solver();
                    auto* sw_solver = (sim_shallowwater_) ? initializer_->get_sw_solver() : nullptr;
                    auto* sw_state = (sim_shallowwater_) ? initializer_->get_sw_state() : nullptr;
                    auto* sw_domain = (sim_shallowwater_) ? initializer_->get_sw_domain() : nullptr;
                    
                    if (gw_solver) {
                        dt_gw = gw_solver->adaptive_time_step(
                            dt_min_, dt_max_, Co_max_,
                            sync_coupling_, dt_sw,
                            sw_state, sw_domain);
                        gw_solver->set_time_step(dt_gw);
                    }
                }
                
                // 3. Handle synchronous coupling: use minimum of both time steps
                if (sync_coupling_ && sim_shallowwater_ && sim_groundwater_) {
                    dt_sw = std::min(dt_sw, dt_gw);
                    dt_gw = dt_sw;
                    if (sim_shallowwater_ && initializer_) {
                        auto* sw_solver = initializer_->get_sw_solver();
                        if (sw_solver) sw_solver->set_time_step(dt_sw);
                    }
                    if (sim_groundwater_ && initializer_) {
                        auto* gw_solver = initializer_->get_gw_solver();
                        if (gw_solver) gw_solver->set_time_step(dt_gw);
                    }
                }
            }
            
            // Update time
            current_time += dt_sw;
            if (!sim_groundwater_ || sync_coupling_) {
                t_subsurface = current_time;
            }
            step++;
            
            if (step % 10 == 0) {
                std::cout << "Step: " << step << " | Time: " << current_time 
                         << "s / " << end_time_ << "s | dt_sw: " << dt_sw 
                         << " | dt_gw: " << dt_gw << "\r" << std::flush;
            }
            
            // Note: Boundary conditions and source/sink terms are automatically
            // handled by the solvers when solve(current_time) is called.
            // The solvers access their respective BC managers and source/sink managers
            // which were set up during initialization.
            
            // 1. Solve Shallow Water (if enabled)
            if (sim_shallowwater_ && initializer_) {
                auto* solver = initializer_->get_sw_solver();
                if (solver) solver->solve(current_time);
            }
            
            // 2. Solve Groundwater (if enabled)
            if (sim_groundwater_ && initializer_) {
                auto* solver = initializer_->get_gw_solver();
                auto* sw_solver = (sim_shallowwater_) ? initializer_->get_sw_solver() : nullptr;
                auto* sw_state = (sim_shallowwater_) ? initializer_->get_sw_state() : nullptr;
                auto* sw_domain = (sim_shallowwater_) ? initializer_->get_sw_domain() : nullptr;
                
                if (solver) {
                    if (sync_coupling_) {
                        // Synchronous coupling: same time step
                        solver->solve(current_time);
                    } else {
                        // Asynchronous coupling: multiple groundwater steps per surface step
                        while (t_subsurface + dt_gw <= current_time) {
                            t_subsurface += dt_gw;
                            solver->solve(t_subsurface);
                            
                            // Recompute adaptive time step for next groundwater step
                            if (use_adaptive_dt_) {
                                dt_gw = solver->adaptive_time_step(
                                    dt_min_, dt_max_, Co_max_,
                                    sync_coupling_, dt_sw,
                                    sw_state, sw_domain);
                                solver->set_time_step(dt_gw);
                            }
                            
                            std::cout << "   >>> Subsurface executed with dt = " << dt_gw << "\n";
                        }
                    }
                    
                    // Compute seepage flux from groundwater to surface
                    if (sw_solver && sw_state && sw_domain) {
                        solver->compute_seepage_flux(*sw_state, *sw_domain);
                    }
                }
            }
            
            // 3. Apply seepage to surface water (if both domains are active)
            if (sim_shallowwater_ && sim_groundwater_ && initializer_) {
                auto* solver = initializer_->get_sw_solver();
                if (solver) {
                    solver->apply_seepage();
                }
            }
            
            // 4. Update Shallow Water Velocity (after both solvers)
            if (sim_shallowwater_ && initializer_) {
                auto* solver = initializer_->get_sw_solver();
                if (solver) solver->update_velocity();
            }
            
            // 5. Solve Scalar Transport (if enabled)
            if (n_scalar_ > 0 && initializer_) {
                // Surface water scalar transport
                if (sim_shallowwater_) {
                    const auto& sw_solvers = initializer_->get_sw_scalar_solvers();
                    const auto& gw_solvers = initializer_->get_gw_scalar_solvers();
                    auto* gw_state = (sim_groundwater_) ? initializer_->get_gw_state() : nullptr;
                    auto* gw_domain = (sim_groundwater_) ? initializer_->get_gw_domain() : nullptr;
                    auto* sw_state = initializer_->get_sw_state();
                    auto* sw_domain = initializer_->get_sw_domain();
                    
                    for (size_t kk = 0; kk < sw_solvers.size() && kk < static_cast<size_t>(n_scalar_); ++kk) {
                        if (sw_solvers[kk]) {
                            sw_solvers[kk]->solve(current_time);
                            
                            // Exchange scalar with subsurface
                            if (sim_groundwater_ && gw_state && gw_domain && 
                                kk < gw_solvers.size() && gw_solvers[kk]) {
                                sw_solvers[kk]->exchange_scalar_with_subsurface(
                                    gw_solvers[kk].get(), *gw_state, *gw_domain);
                            }
                        }
                    }
                }
                
                // Groundwater scalar transport
                if (sim_groundwater_) {
                    const auto& gw_solvers = initializer_->get_gw_scalar_solvers();
                    const auto& sw_solvers = initializer_->get_sw_scalar_solvers();
                    auto* sw_state = (sim_shallowwater_) ? initializer_->get_sw_state() : nullptr;
                    auto* sw_domain = (sim_shallowwater_) ? initializer_->get_sw_domain() : nullptr;
                    
                    for (size_t kk = 0; kk < gw_solvers.size() && kk < static_cast<size_t>(n_scalar_); ++kk) {
                        if (gw_solvers[kk]) {
                            gw_solvers[kk]->solve(current_time);
                            
                            // Exchange scalar with surface
                            if (sim_shallowwater_ && sw_state && sw_domain &&
                                kk < sw_solvers.size() && sw_solvers[kk]) {
                                gw_solvers[kk]->exchange_scalar_with_surface(
                                    sw_solvers[kk].get(), *sw_state, *sw_domain);
                            }
                        }
                    }
                    
                    // Update density/viscosity from scalar concentration (for baroclinic flow)
                    // Uses the first scalar species (typically salinity)
                    if (baroclinic_ && n_scalar_ > 0 && !gw_solvers.empty() && gw_solvers[0]) {
                        auto* gw_solver = initializer_->get_gw_solver();
                        if (gw_solver) {
                            // Default coefficients for seawater salinity:
                            // r_rho = 1.0 + S * 0.000744 (density ratio)
                            // r_visc = 1.0 / (1.0 + S * 0.0022) (viscosity ratio)
                            gw_solver->update_density_viscosity_from_scalar(
                                gw_solvers[0]->scalar_concentration, 0.000744, 0.0022);
                        }
                    }
                }
            }
            
            // 5. Write output at specified intervals
            if (output_writer_) {
                // Write spatial fields if it's time
                if (output_writer_->should_write_spatial(current_time, step)) {
                    output_writer_->write_spatial_fields(current_time, step);
                }
                
                // Write time-series data every step (continuous monitoring)
                output_writer_->write_time_series(current_time, step);
            }
            
            // 6. Report time step completion
            if (step % 100 == 0) {
                std::cout << "\n[Driver] Step " << step << " (" << current_time 
                         << "s of " << end_time_ << "s) completed, dt = " << dt_ << "\n";
            }
        }
        
        double total_time = timer_.seconds();
        std::cout << "\n[Driver] Simulation completed!\n";
        std::cout << "[Driver] Total steps: " << step << "\n";
        std::cout << "[Driver] Total wall time: " << total_time << " s\n";
        std::cout << "[Driver] Average time per step: " << total_time / step << " s\n";
    }
    
private:
    // Note: Boundary conditions, source/sink terms, and output writing are now
    // handled automatically by the solvers and OutputWriter classes.
    // No additional helper functions are needed here.
};

} // namespace Frehg

#endif // FREHG_MODEL_DRIVER_HPP
