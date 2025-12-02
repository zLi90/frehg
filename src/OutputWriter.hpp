#ifndef FREHG_OUTPUT_WRITER_HPP
#define FREHG_OUTPUT_WRITER_HPP

#include "define.hpp"
#include "Domain.hpp"
#include "StateVariables.hpp"
#include "ActiveCellMesh.hpp"
#include "ScalarTransportSolver.hpp"
#include <string>
#include <fstream>
#include <vector>
#include <memory>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <cmath>
#include <algorithm>

namespace Frehg {

// ============================================================================
//                      OUTPUT WRITER CLASS
// ============================================================================
// Handles writing model outputs:
// 1. Spatial fields in VTK format at user-defined intervals
// 2. Time-series data for global variables (continuous output)
// ============================================================================

class OutputWriter {
private:
    // --- Output Configuration ---
    std::string output_dir_;              // Output directory path
    Scalar dt_out_;                       // Output interval for spatial fields
    int output_interval_;                 // Output every N time steps
    bool write_sw_fields_;                // Write surface water fields
    bool write_gw_fields_;                // Write groundwater fields
    bool write_scalar_fields_;           // Write scalar concentration fields
    bool write_time_series_;             // Write time-series data
    
    // --- Output Counters ---
    int output_step_counter_;            // Counter for spatial field outputs
    int time_series_counter_;             // Counter for time-series outputs
    
    // --- Output Timing ---
    mutable Scalar next_output_time_;     // Next time to write output
    mutable Scalar output_time_tolerance_; // Tolerance for time comparison
    
    // --- Domain and State Pointers ---
    SwDomain* sw_domain_;
    GwDomain* gw_domain_;
    SwStateVariables* sw_state_;
    GwStateVariables* gw_state_;
    ActiveCellMesh* sw_active_mesh_;
    ActiveCellMesh* gw_active_mesh_;
    
    // --- Scalar Transport Solvers ---
    std::vector<SwScalarTransportSolver*> sw_scalar_solvers_;
    std::vector<GwScalarTransportSolver*> gw_scalar_solvers_;
    
    // --- Time-Series Output File ---
    std::ofstream time_series_file_;
    std::string time_series_filename_;
    
    // --- Output Timing ---
    mutable Scalar next_output_time_;     // Next time to write output
    mutable Scalar output_time_tolerance_; // Tolerance for time comparison
    
    // --- Boundary Condition Managers (for computing fluxes) ---
    // Note: These are optional - boundary fluxes can be computed if available
    // For now, we'll compute volumes and leave flux computation as placeholder
    
    // --- Flags for one-time output ---
    bool bathymetry_written_;             // Track if bathymetry has been written
    
public:
    // ========================================================================
    // Constructor
    // ========================================================================
    OutputWriter(const std::string& output_dir,
                 Scalar dt_out,
                 SwDomain* sw_domain = nullptr,
                 GwDomain* gw_domain = nullptr,
                 SwStateVariables* sw_state = nullptr,
                 GwStateVariables* gw_state = nullptr,
                 ActiveCellMesh* sw_active_mesh = nullptr,
                 ActiveCellMesh* gw_active_mesh = nullptr)
        : output_dir_(output_dir),
          dt_out_(dt_out),
          output_interval_(1),
          write_sw_fields_(sw_domain != nullptr),
          write_gw_fields_(gw_domain != nullptr),
          write_scalar_fields_(false),
          write_time_series_(true),
          output_step_counter_(0),
          time_series_counter_(0),
          next_output_time_(0.0),
          output_time_tolerance_(1.0e-10),
          bathymetry_written_(false),
          sw_domain_(sw_domain),
          gw_domain_(gw_domain),
          sw_state_(sw_state),
          gw_state_(gw_state),
          sw_active_mesh_(sw_active_mesh),
          gw_active_mesh_(gw_active_mesh) {
        
        // Calculate output interval based on dt_out
        // This will be set properly during initialization
    }
    
    // ========================================================================
    // Destructor
    // ========================================================================
    ~OutputWriter() {
        if (time_series_file_.is_open()) {
            time_series_file_.close();
        }
    }
    
    // ========================================================================
    // Initialize Output Writer
    // ========================================================================
    // Sets up output directories and files
    void initialize(Scalar dt, int n_scalar = 0) {
        // Create output directory if it doesn't exist
        struct stat info;
        if (stat(output_dir_.c_str(), &info) != 0) {
            // Directory doesn't exist, create it
            #ifdef _WIN32
                _mkdir(output_dir_.c_str());
            #else
                mkdir(output_dir_.c_str(), 0755);
            #endif
        }
        
        // Set output interval based on dt_out_ and dt
        if (dt_out_ > 0.0 && dt > 0.0) {
            output_interval_ = static_cast<int>(std::ceil(dt_out_ / dt));
            if (output_interval_ < 1) {
                output_interval_ = 1;
            }
        } else {
            output_interval_ = 1;  // Write every step if dt_out not specified
        }
        
        // Initialize next output time (start at t=0 to write initial conditions)
        next_output_time_ = 0.0;
        
        // Set time tolerance based on dt_out or a small fraction of dt
        if (dt_out_ > 0.0) {
            output_time_tolerance_ = dt_out_ * 1.0e-6;  // 0.0001% of output interval
        } else {
            output_time_tolerance_ = dt * 1.0e-6;
        }
        
        // Initialize time-series output file
        time_series_filename_ = output_dir_ + "/time_series.txt";
        time_series_file_.open(time_series_filename_, std::ios::out);
        
        if (!time_series_file_.is_open()) {
            throw std::runtime_error("Cannot open time-series output file: " + time_series_filename_);
        }
        
        // Write header for time-series file
        time_series_file_ << std::scientific << std::setprecision(10);
        time_series_file_ << "# Time-series output for FREHG model\n";
        time_series_file_ << "# Format: time step time";
        
        if (write_sw_fields_) {
            time_series_file_ << " sw_total_volume sw_inflow sw_outflow";
        }
        if (write_gw_fields_) {
            time_series_file_ << " gw_total_volume gw_inflow gw_outflow";
        }
        if (write_scalar_fields_) {
            for (int k = 0; k < n_scalar; ++k) {
                time_series_file_ << " sw_scalar_" << k << "_mass gw_scalar_" << k << "_mass";
            }
        }
        time_series_file_ << "\n";
        time_series_file_ << "# Units: time [s], volume [m³], flux [m³/s], mass [kg]\n";
        time_series_file_.flush();
    }
    
    // ========================================================================
    // Set Scalar Transport Solvers
    // ========================================================================
    void set_scalar_solvers(
        const std::vector<SwScalarTransportSolver*>& sw_scalars,
        const std::vector<GwScalarTransportSolver*>& gw_scalars) {
        sw_scalar_solvers_ = sw_scalars;
        gw_scalar_solvers_ = gw_scalars;
        write_scalar_fields_ = (!sw_scalars.empty() || !gw_scalars.empty());
    }
    
    // ========================================================================
    // Check if Spatial Output Should Be Written
    // ========================================================================
    bool should_write_spatial(Scalar current_time, int time_step) const {
        // Time-based check: more robust for variable time steps
        // Check if current_time has reached or passed the next output time
        
        if (current_time >= next_output_time_ - output_time_tolerance_) {
            // Update next_output_time_ for the next output
            // Use while loop to handle cases where multiple outputs were skipped
            while (next_output_time_ <= current_time + output_time_tolerance_) {
                next_output_time_ += dt_out_;
            }
            return true;
        }
        
        return false;
    }
    
    // ========================================================================
    // Write Spatial Fields (VTK Format)
    // ========================================================================
    // Writes all spatial fields at current time step
    void write_spatial_fields(Scalar current_time, int time_step) {
        // Write surface water fields if enabled (includes surface scalars)
        if (write_sw_fields_ && sw_domain_ && sw_state_) {
            write_sw_vtk(current_time, time_step);
        }
        
        // Write groundwater fields if enabled (includes groundwater scalars)
        if (write_gw_fields_ && gw_domain_ && gw_state_) {
            write_gw_vtk(current_time, time_step);
        }
        
        output_step_counter_++;
    }
    
    // ========================================================================
    // Write Time-Series Data
    // ========================================================================
    // Writes global variables (total volume, boundary fluxes, etc.)
    void write_time_series(Scalar current_time, int time_step) {
        if (!write_time_series_ || !time_series_file_.is_open()) {
            return;
        }
        
        // Write time step and time
        time_series_file_ << time_step << " " << current_time;
        
        // Compute boundary fluxes once (if needed)
        Scalar sw_inflow = 0.0, sw_outflow = 0.0;
        Scalar gw_inflow = 0.0, gw_outflow = 0.0;
        if (write_sw_fields_ || write_gw_fields_) {
            compute_boundary_fluxes(sw_inflow, sw_outflow, gw_inflow, gw_outflow);
        }
        
        // Compute and write surface water variables
        if (write_sw_fields_ && sw_domain_ && sw_state_ && sw_active_mesh_) {
            Scalar sw_volume = compute_total_sw_volume();
            time_series_file_ << " " << sw_volume << " " << sw_inflow << " " << sw_outflow;
        }
        
        // Compute and write groundwater variables
        if (write_gw_fields_ && gw_domain_ && gw_state_ && gw_active_mesh_) {
            Scalar gw_volume = compute_total_gw_volume();
            time_series_file_ << " " << gw_volume << " " << gw_inflow << " " << gw_outflow;
        }
        
        // Compute and write scalar mass
        if (write_scalar_fields_) {
            // Surface water scalar mass
            for (size_t k = 0; k < sw_scalar_solvers_.size(); ++k) {
                if (sw_scalar_solvers_[k] && sw_active_mesh_) {
                    Scalar scalar_mass = compute_total_sw_scalar_mass(k);
                    time_series_file_ << " " << scalar_mass;
                } else {
                    time_series_file_ << " 0.0";
                }
            }
            
            // Groundwater scalar mass
            for (size_t k = 0; k < gw_scalar_solvers_.size(); ++k) {
                if (gw_scalar_solvers_[k] && gw_active_mesh_) {
                    Scalar scalar_mass = compute_total_gw_scalar_mass(k);
                    time_series_file_ << " " << scalar_mass;
                } else {
                    time_series_file_ << " 0.0";
                }
            }
        }
        
        time_series_file_ << "\n";
        time_series_file_.flush();
        time_series_counter_++;
    }
    
private:
    // ========================================================================
    // Write Surface Water VTK File
    // ========================================================================
    void write_sw_vtk(Scalar current_time, int time_step) {
        std::string filename = generate_filename("sw", "vtk", time_step);
        std::ofstream vtk_file(filename, std::ios::out);
        
        if (!vtk_file.is_open()) {
            std::cerr << "Warning: Cannot open VTK file: " << filename << std::endl;
            return;
        }
        
        // Get domain dimensions
        Ordinal nx = sw_domain_->nx;
        Ordinal ny = sw_domain_->ny;
        Scalar dx = sw_domain_->dx;
        Scalar dy = sw_domain_->dy;
        Ordinal num_points = (nx + 1) * (ny + 1);
        Ordinal num_cells = nx * ny;
        
        // Create host mirrors for data access
        auto h_pressure = Kokkos::create_mirror_view(sw_state_->pressure);
        auto h_depth = Kokkos::create_mirror_view(sw_state_->depth);
        auto h_velocity_x = Kokkos::create_mirror_view(sw_state_->velocity_x);
        auto h_velocity_y = Kokkos::create_mirror_view(sw_state_->velocity_y);
        auto h_bottom = Kokkos::create_mirror_view(sw_state_->bottom);
        auto h_active_mask = Kokkos::create_mirror_view(sw_domain_->active_mask);
        
        Kokkos::deep_copy(h_pressure, sw_state_->pressure);
        Kokkos::deep_copy(h_depth, sw_state_->depth);
        Kokkos::deep_copy(h_velocity_x, sw_state_->velocity_x);
        Kokkos::deep_copy(h_velocity_y, sw_state_->velocity_y);
        Kokkos::deep_copy(h_bottom, sw_state_->bottom);
        Kokkos::deep_copy(h_active_mask, sw_domain_->active_mask);
        
        // Write VTK header (Legacy format, version 2.0)
        vtk_file << "# vtk DataFile Version 2.0\n";
        vtk_file << "FREHG Surface Water Output, Time = " << current_time << " s\n";
        vtk_file << "ASCII\n";
        vtk_file << "DATASET STRUCTURED_GRID\n";
        vtk_file << "DIMENSIONS " << (nx + 1) << " " << (ny + 1) << " 1\n";
        vtk_file << "POINTS " << num_points << " double\n";
        
        // Write grid points (vertices at cell corners)
        // Points are ordered: (i, j, k) where i varies fastest
        vtk_file << std::scientific << std::setprecision(10);
        for (Ordinal j = 0; j <= ny; ++j) {
            for (Ordinal i = 0; i <= nx; ++i) {
                Scalar x = i * dx;
                Scalar y = j * dy;
                Scalar z = 0.0;  // 2D, so z = 0
                vtk_file << x << " " << y << " " << z << "\n";
            }
        }
        
        // Write cell data
        vtk_file << "\nCELL_DATA " << num_cells << "\n";
        
        // Free surface elevation (eta)
        vtk_file << "SCALARS eta double 1\n";
        vtk_file << "LOOKUP_TABLE default\n";
        for (Ordinal j = 0; j < ny; ++j) {
            for (Ordinal i = 0; i < nx; ++i) {
                Ordinal cell_idx = i + j * nx;
                vtk_file << h_pressure(cell_idx) << "\n";
            }
        }
        
        // Water depth
        vtk_file << "\nSCALARS depth double 1\n";
        vtk_file << "LOOKUP_TABLE default\n";
        for (Ordinal j = 0; j < ny; ++j) {
            for (Ordinal i = 0; i < nx; ++i) {
                Ordinal cell_idx = i + j * nx;
                vtk_file << h_depth(cell_idx) << "\n";
            }
        }
        
        // Velocity (vector field)
        vtk_file << "\nVECTORS velocity double\n";
        for (Ordinal j = 0; j < ny; ++j) {
            for (Ordinal i = 0; i < nx; ++i) {
                Ordinal cell_idx = i + j * nx;
                vtk_file << h_velocity_x(cell_idx) << " " 
                        << h_velocity_y(cell_idx) << " 0.0\n";
            }
        }
        
        // Active mask (for visualization)
        vtk_file << "\nSCALARS active int 1\n";
        vtk_file << "LOOKUP_TABLE default\n";
        for (Ordinal j = 0; j < ny; ++j) {
            for (Ordinal i = 0; i < nx; ++i) {
                Ordinal cell_idx = i + j * nx;
                vtk_file << h_active_mask(cell_idx) << "\n";
            }
        }
        
        // Bathymetry (only write once if not already written)
        // Or write every time - user said "only save once" but it's useful to have in each file
        // for independent visualization. We'll write it but mark that it's time-invariant.
        vtk_file << "\nSCALARS bathymetry double 1\n";
        vtk_file << "LOOKUP_TABLE default\n";
        for (Ordinal j = 0; j < ny; ++j) {
            for (Ordinal i = 0; i < nx; ++i) {
                Ordinal cell_idx = i + j * nx;
                vtk_file << h_bottom(cell_idx) << "\n";
            }
        }
        
        // Write scalar concentrations if available
        if (write_scalar_fields_ && !sw_scalar_solvers_.empty()) {
            for (size_t k = 0; k < sw_scalar_solvers_.size(); ++k) {
                if (sw_scalar_solvers_[k]) {
                    auto _scalar_concentration = sw_scalar_solvers_[k]->scalar_concentration;
                    auto h_scalar = Kokkos::create_mirror_view(_scalar_concentration);
                    Kokkos::deep_copy(h_scalar, _scalar_concentration);
                    
                    vtk_file << "\nSCALARS scalar_" << k << " double 1\n";
                    vtk_file << "LOOKUP_TABLE default\n";
                    for (Ordinal j = 0; j < ny; ++j) {
                        for (Ordinal i = 0; i < nx; ++i) {
                            Ordinal cell_idx = i + j * nx;
                            vtk_file << h_scalar(cell_idx) << "\n";
                        }
                    }
                }
            }
        }
        
        vtk_file.close();
        
        // Mark bathymetry as written (for potential separate file in the future)
        if (!bathymetry_written_) {
            bathymetry_written_ = true;
            // Optionally write a separate bathymetry-only file
            write_bathymetry_vtk();
        }
    }
    
    // ========================================================================
    // Write Groundwater VTK File
    // ========================================================================
    void write_gw_vtk(Scalar current_time, int time_step) {
        std::string filename = generate_filename("gw", "vtk", time_step);
        std::ofstream vtk_file(filename, std::ios::out);
        
        if (!vtk_file.is_open()) {
            std::cerr << "Warning: Cannot open VTK file: " << filename << std::endl;
            return;
        }
        
        // Get domain dimensions
        Ordinal nx = gw_domain_->nx;
        Ordinal ny = gw_domain_->ny;
        Ordinal nz = gw_domain_->nz;
        Scalar dx = gw_domain_->dx;
        Scalar dy = gw_domain_->dy;
        Scalar dz = gw_domain_->dz;  // Assuming uniform dz for now
        Ordinal num_points = (nx + 1) * (ny + 1) * (nz + 1);
        Ordinal num_cells = nx * ny * nz;
        
        // Create host mirrors for data access
        auto h_pressure = Kokkos::create_mirror_view(gw_state_->pressure);
        auto h_water_content = Kokkos::create_mirror_view(gw_state_->water_content);
        auto h_velocity_x = Kokkos::create_mirror_view(gw_state_->velocity_x);
        auto h_velocity_y = Kokkos::create_mirror_view(gw_state_->velocity_y);
        auto h_flux_x = Kokkos::create_mirror_view(gw_state_->flux_x);
        auto h_flux_y = Kokkos::create_mirror_view(gw_state_->flux_y);
        auto h_flux_z = Kokkos::create_mirror_view(gw_state_->flux_z);
        auto h_active_mask_3d = Kokkos::create_mirror_view(gw_domain_->active_mask_3d);
        
        Kokkos::deep_copy(h_pressure, gw_state_->pressure);
        Kokkos::deep_copy(h_water_content, gw_state_->water_content);
        Kokkos::deep_copy(h_velocity_x, gw_state_->velocity_x);
        Kokkos::deep_copy(h_velocity_y, gw_state_->velocity_y);
        Kokkos::deep_copy(h_flux_x, gw_state_->flux_x);
        Kokkos::deep_copy(h_flux_y, gw_state_->flux_y);
        Kokkos::deep_copy(h_flux_z, gw_state_->flux_z);
        Kokkos::deep_copy(h_active_mask_3d, gw_domain_->active_mask_3d);
        
        // Write VTK header (Legacy format, version 2.0)
        vtk_file << "# vtk DataFile Version 2.0\n";
        vtk_file << "FREHG Groundwater Output, Time = " << current_time << " s\n";
        vtk_file << "ASCII\n";
        vtk_file << "DATASET STRUCTURED_GRID\n";
        vtk_file << "DIMENSIONS " << (nx + 1) << " " << (ny + 1) << " " << (nz + 1) << "\n";
        vtk_file << "POINTS " << num_points << " double\n";
        
        // Write grid points (vertices at cell corners)
        // Points are ordered: (i, j, k) where i varies fastest
        vtk_file << std::scientific << std::setprecision(10);
        for (Ordinal k = 0; k <= nz; ++k) {
            for (Ordinal j = 0; j <= ny; ++j) {
                for (Ordinal i = 0; i <= nx; ++i) {
                    Scalar x = i * dx;
                    Scalar y = j * dy;
                    Scalar z = k * dz;  // Depth coordinate (0 at top, increases downward)
                    vtk_file << x << " " << y << " " << z << "\n";
                }
            }
        }
        
        // Write cell data
        vtk_file << "\nCELL_DATA " << num_cells << "\n";
        
        // Pressure head (h)
        vtk_file << "SCALARS pressure_head double 1\n";
        vtk_file << "LOOKUP_TABLE default\n";
        for (Ordinal k = 0; k < nz; ++k) {
            for (Ordinal j = 0; j < ny; ++j) {
                for (Ordinal i = 0; i < nx; ++i) {
                    Ordinal cell_idx = i + j * nx + k * nx * ny;
                    vtk_file << h_pressure(cell_idx) << "\n";
                }
            }
        }
        
        // Water content
        vtk_file << "\nSCALARS water_content double 1\n";
        vtk_file << "LOOKUP_TABLE default\n";
        for (Ordinal k = 0; k < nz; ++k) {
            for (Ordinal j = 0; j < ny; ++j) {
                for (Ordinal i = 0; i < nx; ++i) {
                    Ordinal cell_idx = i + j * nx + k * nx * ny;
                    vtk_file << h_water_content(cell_idx) << "\n";
                }
            }
        }
        
        // Velocity (3D vector field)
        // Note: Groundwater velocity can be very small, so we maintain full precision
        // Velocity can be computed from Darcy flux: v = q / (porosity * area)
        // For now, we'll use the stored velocity if available, or compute from flux
        vtk_file << "\nVECTORS velocity double\n";
        for (Ordinal k = 0; k < nz; ++k) {
            for (Ordinal j = 0; j < ny; ++j) {
                for (Ordinal i = 0; i < nx; ++i) {
                    Ordinal cell_idx = i + j * nx + k * nx * ny;
                    
                    // Use stored velocity if available (may be computed in solver)
                    // Otherwise, velocity is very small in groundwater
                    Scalar vx = h_velocity_x(cell_idx);
                    Scalar vy = h_velocity_y(cell_idx);
                    
                    // For z-component, we can use flux_z as an approximation
                    // In groundwater, velocity_z might not be explicitly stored
                    // We'll use flux_z directly or set to a small value
                    Scalar vz = 0.0;
                    if (h_water_content(cell_idx) > 0.0) {
                        // Approximate vertical velocity from vertical flux
                        // This is a simplification; actual implementation may vary
                        vz = h_flux_z(cell_idx) / (h_water_content(cell_idx) * dx * dy);
                    }
                    
                    vtk_file << vx << " " << vy << " " << vz << "\n";
                }
            }
        }
        
        // Darcy flux (alternative to velocity, useful for groundwater)
        vtk_file << "\nVECTORS darcy_flux double\n";
        for (Ordinal k = 0; k < nz; ++k) {
            for (Ordinal j = 0; j < ny; ++j) {
                for (Ordinal i = 0; i < nx; ++i) {
                    Ordinal cell_idx = i + j * nx + k * nx * ny;
                    vtk_file << h_flux_x(cell_idx) << " " 
                            << h_flux_y(cell_idx) << " " 
                            << h_flux_z(cell_idx) << "\n";
                }
            }
        }
        
        // Active mask (for visualization)
        vtk_file << "\nSCALARS active int 1\n";
        vtk_file << "LOOKUP_TABLE default\n";
        for (Ordinal k = 0; k < nz; ++k) {
            for (Ordinal j = 0; j < ny; ++j) {
                for (Ordinal i = 0; i < nx; ++i) {
                    Ordinal cell_idx = i + j * nx + k * nx * ny;
                    vtk_file << h_active_mask_3d(cell_idx) << "\n";
                }
            }
        }
        
        // Write scalar concentrations if available
        if (write_scalar_fields_ && !gw_scalar_solvers_.empty()) {
            for (size_t kk = 0; kk < gw_scalar_solvers_.size(); ++kk) {
                if (gw_scalar_solvers_[kk]) {
                    auto _scalar_concentration = gw_scalar_solvers_[kk]->scalar_concentration;
                    auto h_scalar = Kokkos::create_mirror_view(_scalar_concentration);
                    Kokkos::deep_copy(h_scalar, _scalar_concentration);
                    
                    vtk_file << "\nSCALARS scalar_" << kk << " double 1\n";
                    vtk_file << "LOOKUP_TABLE default\n";
                    for (Ordinal k = 0; k < nz; ++k) {
                        for (Ordinal j = 0; j < ny; ++j) {
                            for (Ordinal i = 0; i < nx; ++i) {
                                Ordinal cell_idx = i + j * nx + k * nx * ny;
                                vtk_file << h_scalar(cell_idx) << "\n";
                            }
                        }
                    }
                }
            }
        }
        
        vtk_file.close();
    }
    
    
    // ========================================================================
    // Helper: Write Bathymetry-Only VTK File (one-time)
    // ========================================================================
    void write_bathymetry_vtk() {
        std::string filename = output_dir_ + "/bathymetry.vtk";
        std::ofstream vtk_file(filename, std::ios::out);
        
        if (!vtk_file.is_open()) {
            std::cerr << "Warning: Cannot open bathymetry VTK file: " << filename << std::endl;
            return;
        }
        
        // Get domain dimensions
        Ordinal nx = sw_domain_->nx;
        Ordinal ny = sw_domain_->ny;
        Scalar dx = sw_domain_->dx;
        Scalar dy = sw_domain_->dy;
        Ordinal num_points = (nx + 1) * (ny + 1);
        Ordinal num_cells = nx * ny;
        
        // Create host mirror for bathymetry
        auto h_bottom = Kokkos::create_mirror_view(sw_state_->bottom);
        auto h_active_mask = Kokkos::create_mirror_view(sw_domain_->active_mask);
        Kokkos::deep_copy(h_bottom, sw_state_->bottom);
        Kokkos::deep_copy(h_active_mask, sw_domain_->active_mask);
        
        // Write VTK header
        vtk_file << "# vtk DataFile Version 2.0\n";
        vtk_file << "FREHG Bathymetry (Time-Invariant)\n";
        vtk_file << "ASCII\n";
        vtk_file << "DATASET STRUCTURED_GRID\n";
        vtk_file << "DIMENSIONS " << (nx + 1) << " " << (ny + 1) << " 1\n";
        vtk_file << "POINTS " << num_points << " double\n";
        
        // Write grid points
        vtk_file << std::scientific << std::setprecision(10);
        for (Ordinal j = 0; j <= ny; ++j) {
            for (Ordinal i = 0; i <= nx; ++i) {
                Scalar x = i * dx;
                Scalar y = j * dy;
                Scalar z = 0.0;
                vtk_file << x << " " << y << " " << z << "\n";
            }
        }
        
        // Write cell data
        vtk_file << "\nCELL_DATA " << num_cells << "\n";
        
        // Bathymetry
        vtk_file << "SCALARS bathymetry double 1\n";
        vtk_file << "LOOKUP_TABLE default\n";
        for (Ordinal j = 0; j < ny; ++j) {
            for (Ordinal i = 0; i < nx; ++i) {
                Ordinal cell_idx = i + j * nx;
                vtk_file << h_bottom(cell_idx) << "\n";
            }
        }
        
        // Active mask
        vtk_file << "\nSCALARS active int 1\n";
        vtk_file << "LOOKUP_TABLE default\n";
        for (Ordinal j = 0; j < ny; ++j) {
            for (Ordinal i = 0; i < nx; ++i) {
                Ordinal cell_idx = i + j * nx;
                vtk_file << h_active_mask(cell_idx) << "\n";
            }
        }
        
        vtk_file.close();
    }
    
    // ========================================================================
    // Helper: Generate Output Filename
    // ========================================================================
    std::string generate_filename(const std::string& prefix, 
                                  const std::string& suffix,
                                  int time_step) const {
        std::ostringstream oss;
        oss << output_dir_ << "/" << prefix << "_" 
            << std::setfill('0') << std::setw(6) << time_step 
            << "." << suffix;
        return oss.str();
    }
    
    // ========================================================================
    // Helper: Compute Global Surface Water Volume
    // ========================================================================
    Scalar compute_total_sw_volume() const {
        if (!sw_state_ || !sw_active_mesh_) {
            return 0.0;
        }
        
        auto _volume = sw_state_->volume;
        auto _active_mask = sw_domain_->active_mask;
        auto _active_to_domain = sw_active_mesh_->active_to_domain;
        
        // Create host mirror and compute sum
        auto h_volume = Kokkos::create_mirror_view(_volume);
        auto h_active_mask = Kokkos::create_mirror_view(_active_mask);
        Kokkos::deep_copy(h_volume, _volume);
        Kokkos::deep_copy(h_active_mask, _active_mask);
        
        Scalar total_volume = 0.0;
        for (Ordinal i = 0; i < sw_active_mesh_->num_active; ++i) {
            Ordinal domain_idx = _active_to_domain(i);
            if (h_active_mask(domain_idx) > 0) {
                total_volume += h_volume(domain_idx);
            }
        }
        
        return total_volume;
    }
    
    // ========================================================================
    // Helper: Compute Global Groundwater Volume
    // ========================================================================
    Scalar compute_total_gw_volume() const {
        if (!gw_state_ || !gw_active_mesh_) {
            return 0.0;
        }
        
        auto _water_content = gw_state_->water_content;
        auto _active_mask_3d = gw_domain_->active_mask_3d;
        auto _active_to_domain = gw_active_mesh_->active_to_domain;
        
        // Create host mirror and compute sum
        auto h_water_content = Kokkos::create_mirror_view(_water_content);
        auto h_active_mask_3d = Kokkos::create_mirror_view(_active_mask_3d);
        Kokkos::deep_copy(h_water_content, _water_content);
        Kokkos::deep_copy(h_active_mask_3d, _active_mask_3d);
        
        // For groundwater, we need cell volume (dx * dy * dz)
        // Assuming uniform grid for now
        Scalar cell_volume = gw_domain_->dx * gw_domain_->dy * gw_domain_->dz;
        
        Scalar total_volume = 0.0;
        for (Ordinal i = 0; i < gw_active_mesh_->num_active; ++i) {
            Ordinal domain_idx = _active_to_domain(i);
            if (h_active_mask_3d(domain_idx) > 0) {
                total_volume += h_water_content(domain_idx) * cell_volume;
            }
        }
        
        return total_volume;
    }
    
    // ========================================================================
    // Helper: Compute Total Surface Water Scalar Mass
    // ========================================================================
    Scalar compute_total_sw_scalar_mass(size_t scalar_index) const {
        if (scalar_index >= sw_scalar_solvers_.size() || !sw_scalar_solvers_[scalar_index]) {
            return 0.0;
        }
        
        auto* solver = sw_scalar_solvers_[scalar_index];
        auto _scalar_mass = solver->scalar_mass;
        auto _active_mask = sw_domain_->active_mask;
        auto _active_to_domain = sw_active_mesh_->active_to_domain;
        
        // Create host mirror and compute sum
        auto h_scalar_mass = Kokkos::create_mirror_view(_scalar_mass);
        auto h_active_mask = Kokkos::create_mirror_view(_active_mask);
        Kokkos::deep_copy(h_scalar_mass, _scalar_mass);
        Kokkos::deep_copy(h_active_mask, _active_mask);
        
        Scalar total_mass = 0.0;
        for (Ordinal i = 0; i < sw_active_mesh_->num_active; ++i) {
            Ordinal domain_idx = _active_to_domain(i);
            if (h_active_mask(domain_idx) > 0) {
                total_mass += h_scalar_mass(domain_idx);
            }
        }
        
        return total_mass;
    }
    
    // ========================================================================
    // Helper: Compute Total Groundwater Scalar Mass
    // ========================================================================
    Scalar compute_total_gw_scalar_mass(size_t scalar_index) const {
        if (scalar_index >= gw_scalar_solvers_.size() || !gw_scalar_solvers_[scalar_index]) {
            return 0.0;
        }
        
        auto* solver = gw_scalar_solvers_[scalar_index];
        auto _scalar_mass = solver->scalar_mass;
        auto _active_mask_3d = gw_domain_->active_mask_3d;
        auto _active_to_domain = gw_active_mesh_->active_to_domain;
        
        // Create host mirror and compute sum
        auto h_scalar_mass = Kokkos::create_mirror_view(_scalar_mass);
        auto h_active_mask_3d = Kokkos::create_mirror_view(_active_mask_3d);
        Kokkos::deep_copy(h_scalar_mass, _scalar_mass);
        Kokkos::deep_copy(h_active_mask_3d, _active_mask_3d);
        
        Scalar total_mass = 0.0;
        for (Ordinal i = 0; i < gw_active_mesh_->num_active; ++i) {
            Ordinal domain_idx = _active_to_domain(i);
            if (h_active_mask_3d(domain_idx) > 0) {
                total_mass += h_scalar_mass(domain_idx);
            }
        }
        
        return total_mass;
    }
    
    // ========================================================================
    // Helper: Compute Boundary Inflow/Outflow Rates
    // ========================================================================
    void compute_boundary_fluxes(Scalar& sw_inflow, Scalar& sw_outflow,
                                 Scalar& gw_inflow, Scalar& gw_outflow) const {
        // TODO: Compute total boundary fluxes from boundary conditions
        // This requires access to boundary condition managers and source/sink managers
        // For now, this is a placeholder - implementation would need:
        // 1. Access to SwBoundaryConditionManager and GwBoundaryConditionManager
        // 2. Access to SwSourceSinkManager
        // 3. Sum of flow rates from boundary conditions and source/sink terms
        
        // Placeholder: set to zero
        sw_inflow = 0.0;
        sw_outflow = 0.0;
        gw_inflow = 0.0;
        gw_outflow = 0.0;
    }
};

} // namespace Frehg

#endif // FREHG_OUTPUT_WRITER_HPP

