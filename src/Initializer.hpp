#ifndef FREHG_INITIALIZER_HPP
#define FREHG_INITIALIZER_HPP

#include "define.hpp"
#include "InputReader.hpp"
#include "Domain.hpp"
#include "StateVariables.hpp"
#include "ActiveCellMesh.hpp"
#include "ShallowWaterSolver.hpp"
#include "GroundwaterSolver.hpp"
#include "ScalarTransportSolver.hpp"
#include <memory>
#include <string>
#include <vector>
#include <limits>
#include <iostream>

namespace Frehg {

// ============================================================================
//                      INITIALIZATION CLASS
// ============================================================================
// Handles all initialization tasks:
// - Reading configuration from input file
// - Initializing domains
// - Initializing state variables
// - Building active cell meshes
// - Creating solvers
// This keeps initialization logic separate from the ModelDriver

class Initializer {
private:
    // Input file reader
    std::unique_ptr<InputReader> input_reader_;
    std::string input_dir_;
    std::string output_dir_;
    
    // Configuration parameters (will be populated from input file)
    struct Config {
        // Directories
        std::string finput;
        std::string foutput;
        std::string sim_id;
        
        // Domain geometry
        Ordinal NX, NY, NZ;
        Scalar dx, dy, dz;
        Scalar botZ;
        Scalar dz_incre;
        
        // MPI settings
        bool use_mpi;
        Ordinal mpi_nx, mpi_ny;
        Ordinal nthreads;
        
        // Time control
        Scalar dt;
        Scalar Tend;
        Scalar dt_out;
        Ordinal n_monitor;
        
        // Physics flags
        bool sim_shallowwater;
        bool sim_groundwater;
        bool sync_coupling;
        
        // Scalar transport
        Ordinal n_scalar;
        bool baroclinic;
        bool superbee;
        
        // Default constructor
        Config() : NX(0), NY(0), NZ(0), dx(0.0), dy(0.0), dz(0.0), botZ(0.0),
                   dz_incre(1.0), use_mpi(false), mpi_nx(1), mpi_ny(1),
                   nthreads(1), dt(1.0), Tend(100.0), dt_out(10.0), n_monitor(0),
                   sim_shallowwater(true), sim_groundwater(true), sync_coupling(true),
                   n_scalar(0), baroclinic(false), superbee(true) {}
    };
    
    Config config_;
    
    // Initialized objects (will be passed to ModelDriver)
    std::unique_ptr<SwDomain> sw_domain_;
    std::unique_ptr<GwDomain> gw_domain_;
    std::unique_ptr<SwStateVariables> sw_state_;
    std::unique_ptr<GwStateVariables> gw_state_;
    std::unique_ptr<ActiveCellMesh> sw_active_mesh_;
    std::unique_ptr<ActiveCellMesh> gw_active_mesh_;
    std::unique_ptr<ShallowWaterSolver> sw_solver_;
    std::unique_ptr<GroundwaterSolver> gw_solver_;
    std::vector<std::unique_ptr<SwScalarTransportSolver>> sw_scalar_solvers_;
    std::vector<std::unique_ptr<GwScalarTransportSolver>> gw_scalar_solvers_;

public:
    // ========================================================================
    // Constructor
    // ========================================================================
    Initializer(const std::string& input_dir, const std::string& output_dir)
        : input_dir_(input_dir), output_dir_(output_dir) {
        // Input file path (assuming config.txt in input directory)
        std::string config_file = input_dir + "/config.txt";
        input_reader_ = std::make_unique<InputReader>(config_file);
    }
    
    // ========================================================================
    // Destructor
    // ========================================================================
    ~Initializer() = default;
    
    // ========================================================================
    // READ CONFIGURATION FROM INPUT FILE
    // ========================================================================
    // This is the main function that reads all parameters from config.txt
    // For now, it's a skeleton that only reads a few essential fields
    // More fields will be added as needed
    void read_configuration() {
        if (!input_reader_) {
            throw std::runtime_error("InputReader not initialized");
        }
        
        // Read directory settings
        config_.finput = input_reader_->read_string("finput", input_dir_);
        config_.foutput = input_reader_->read_string("foutput", output_dir_);
        config_.sim_id = input_reader_->read_string("sim_id", "default");
        
        // Read domain geometry
        config_.NX = input_reader_->read_int("NX", 10);
        config_.NY = input_reader_->read_int("NY", 10);
        config_.NZ = input_reader_->read_int("NZ", 10);
        config_.dx = input_reader_->read_double("dx", 1.0);
        config_.dy = input_reader_->read_double("dy", 1.0);
        config_.dz = input_reader_->read_double("dz", 1.0);
        config_.botZ = input_reader_->read_double("botZ", 0.0);
        config_.dz_incre = input_reader_->read_double("dz_incre", 1.0);
        
        // Read MPI settings
        config_.use_mpi = input_reader_->read_bool("use_mpi", false);
        config_.mpi_nx = input_reader_->read_int("mpi_nx", 1);
        config_.mpi_ny = input_reader_->read_int("mpi_ny", 1);
        config_.nthreads = input_reader_->read_int("nthreads", 1);
        
        // Read time control
        config_.dt = input_reader_->read_double("dt", 1.0);
        config_.Tend = input_reader_->read_double("Tend", 100.0);
        config_.dt_out = input_reader_->read_double("dt_out", 10.0);
        config_.n_monitor = input_reader_->read_int("n_monitor", 0);
        
        // Read physics flags
        config_.sim_shallowwater = input_reader_->read_bool("sim_shallowwater", true);
        config_.sim_groundwater = input_reader_->read_bool("sim_groundwater", true);
        config_.sync_coupling = input_reader_->read_bool("sync_coupling", true);
        
        // Read scalar transport settings
        config_.n_scalar = input_reader_->read_int("n_scalar", 0);
        config_.baroclinic = input_reader_->read_bool("baroclinic", false);
        config_.superbee = input_reader_->read_bool("superbee", true);
        
        // TODO: Add more fields as needed:
        // - Bathymetry file settings
        // - Physical parameters (grav, visc, manning, etc.)
        // - Boundary conditions
        // - Initial conditions
        // - Wind model settings
        // - Evaporation/rainfall settings
        // - Inflow/tide settings
        // - Groundwater parameters (Ksx, Ksy, Ksz, Ss, wcs, wcr, etc.)
        // - Scalar transport parameters (diffusion, dispersion, etc.)
    }
    
    // ========================================================================
    // INITIALIZE DOMAINS
    // ========================================================================
    void initialize_domains() {
        // Read bathymetry file settings
        bool bath_file = input_reader_->read_bool("bath_file", false);
        
        std::vector<Scalar> bathymetry_data;
        std::vector<int> active_mask_data;
        Ordinal dem_ncols = 0, dem_nrows = 0;
        Scalar dem_cellsize = 0.0;
        Scalar dem_nodata = -9999.0;
        
        if (bath_file) {
            // Read DEM file path
            std::string dem_filename = input_reader_->read_string("bathymetry_file", "");
            if (dem_filename.empty()) {
                // Try default name in input directory
                dem_filename = config_.finput + "/bath.asc";
            }
            
            // Read structured data file (ESRI ASCII Grid format)
            StructuredGrid dem_grid = read_structured_data(dem_filename);
            
            dem_ncols = dem_grid.ncols;
            dem_nrows = dem_grid.nrows;
            dem_cellsize = dem_grid.cellsize;
            dem_nodata = dem_grid.nodata_value;
            bathymetry_data = dem_grid.data;
            // Active mask is automatically generated from nodata_value
            active_mask_data = dem_grid.active_mask;
            
            // Update domain dimensions from DEM if not specified in config
            if (config_.NX == 0) config_.NX = dem_ncols;
            if (config_.NY == 0) config_.NY = dem_nrows;
            if (config_.dx == 0.0) config_.dx = dem_cellsize;
            if (config_.dy == 0.0) config_.dy = dem_cellsize;
            
            // Validate dimensions match
            if (config_.NX != dem_ncols || config_.NY != dem_nrows) {
                throw std::runtime_error("Domain dimensions (NX=" + std::to_string(config_.NX) + 
                                       ", NY=" + std::to_string(config_.NY) + 
                                       ") do not match DEM dimensions (" + 
                                       std::to_string(dem_ncols) + "x" + std::to_string(dem_nrows) + 
                                       ") in file: " + dem_filename);
            }
            
            // Validate cell size matches (with tolerance)
            if (std::abs(config_.dx - dem_cellsize) > 1e-6 || 
                std::abs(config_.dy - dem_cellsize) > 1e-6) {
                std::cerr << "WARNING: Domain cell size (dx=" << config_.dx 
                         << ", dy=" << config_.dy 
                         << ") does not match DEM cellsize (" << dem_cellsize 
                         << ") in file: " << dem_filename << std::endl;
            }
            
            // Calculate bathymetry offset (make minimum elevation = 0)
            Scalar z_min = std::numeric_limits<Scalar>::max();
            for (const auto& z : bathymetry_data) {
                if (std::abs(z - dem_nodata) > 1e-10) {
                    if (z < z_min) z_min = z;
                }
            }
            Scalar offset = (z_min < 0.0) ? -z_min : 0.0;
            
            // Apply offset to bathymetry (only to non-nodata values)
            for (auto& z : bathymetry_data) {
                if (std::abs(z - dem_nodata) > 1e-10) {
                    z += offset;
                }
            }
            
        } else {
            // No bathymetry file - create flat domain
            bathymetry_data.resize(config_.NX * config_.NY, config_.botZ);
            active_mask_data.resize(config_.NX * config_.NY, 1);
        }
        
        // Initialize surface water domain
        if (config_.sim_shallowwater) {
            sw_domain_ = std::make_unique<SwDomain>(
                config_.NX, config_.NY, config_.dx, config_.dy);
            sw_domain_->initialize(active_mask_data, bathymetry_data);
        }
        
        // Initialize groundwater domain
        if (config_.sim_groundwater) {
            gw_domain_ = std::make_unique<GwDomain>(
                config_.NX, config_.NY, config_.NZ,
                config_.dx, config_.dy, config_.dz,
                config_.botZ, config_.dz_incre);
            
            // For groundwater, we pass 2D active mask and bathymetry
            // The domain will build 3D structure internally
            gw_domain_->initialize(active_mask_data, bathymetry_data);
        }
    }
    
    // ========================================================================
    // INITIALIZE STATE VARIABLES
    // ========================================================================
    void initialize_state_variables() {
        if (config_.sim_shallowwater && sw_domain_) {
            sw_state_ = std::make_unique<SwStateVariables>(
                sw_domain_->num_cells_total);
            
            // TODO: Set initial conditions from input file
            // - Initial surface elevation (eta)
            // - Initial velocities (u, v)
        }
        
        if (config_.sim_groundwater && gw_domain_) {
            gw_state_ = std::make_unique<GwStateVariables>(
                gw_domain_->num_cells_3d_total);
            
            // TODO: Set initial conditions from input file
            // - Initial head (h)
            // - Initial water content (wc)
            // - Hydraulic conductivity (Ksx, Ksy, Ksz)
            // - Soil properties (wcs, wcr, Ss, etc.)
        }
    }
    
    // ========================================================================
    // BUILD ACTIVE CELL MESHES
    // ========================================================================
    void build_active_meshes() {
        if (config_.sim_shallowwater && sw_domain_) {
            sw_active_mesh_ = std::make_unique<ActiveCellMesh>(
                sw_domain_->num_cells_total);
            sw_active_mesh_->build_from_2d(*sw_domain_);
        }
        
        if (config_.sim_groundwater && gw_domain_) {
            gw_active_mesh_ = std::make_unique<ActiveCellMesh>(
                gw_domain_->num_cells_3d_total);
            gw_active_mesh_->build_from_3d(*gw_domain_);
        }
    }
    
    // ========================================================================
    // CREATE SOLVERS
    // ========================================================================
    void create_solvers() {
        // Create shallow water solver
        if (config_.sim_shallowwater && sw_domain_ && sw_active_mesh_ && sw_state_) {
            // TODO: Read solver parameters from config
            Scalar grav = input_reader_->read_double("grav", 9.81);
            Scalar manning_n = input_reader_->read_double("manning", 0.03);
            Scalar visc_x = input_reader_->read_double("viscx", 0.0);
            Scalar visc_y = input_reader_->read_double("viscy", 0.0);
            
            sw_solver_ = std::make_unique<ShallowWaterSolver>(
                *sw_domain_, *sw_active_mesh_, *sw_state_,
                config_.dt, grav, manning_n, visc_x, visc_y);
        }
        
        // Create groundwater solver
        if (config_.sim_groundwater && gw_domain_ && gw_active_mesh_ && gw_state_) {
            // TODO: Read solver parameters from config
            Scalar Ss = input_reader_->read_double("Ss", 1.0e-5);
            bool use_corrector = input_reader_->read_bool("use_corrector", true);
            bool follow_terrain = input_reader_->read_bool("follow_terrain", false);
            bool baroclinic = config_.baroclinic;
            
            gw_solver_ = std::make_unique<GroundwaterSolver>(
                *gw_domain_, *gw_active_mesh_, *gw_state_,
                config_.dt, Ss, 1.0e-6, use_corrector, follow_terrain, baroclinic);
        }
        
        // Create scalar transport solvers
        if (config_.n_scalar > 0) {
            // Surface water scalar solvers
            if (config_.sim_shallowwater && sw_domain_ && sw_active_mesh_ && sw_state_) {
                // TODO: Read scalar transport parameters
                Scalar diff_x = input_reader_->read_double("difux", 0.0);
                Scalar diff_y = input_reader_->read_double("difuy", 0.0);
                
                for (Ordinal kk = 0; kk < config_.n_scalar; ++kk) {
                    sw_scalar_solvers_.push_back(
                        std::make_unique<SwScalarTransportSolver>(
                            *sw_domain_, *sw_active_mesh_, *sw_state_,
                            config_.dt, diff_x, diff_y, 200.0, 0.0, config_.superbee));
                }
            }
            
            // Groundwater scalar solvers
            if (config_.sim_groundwater && gw_domain_ && gw_active_mesh_ && gw_state_) {
                // TODO: Read scalar transport parameters
                Scalar diff_x = input_reader_->read_double("difux", 0.0);
                Scalar diff_y = input_reader_->read_double("difuy", 0.0);
                Scalar diff_z = input_reader_->read_double("difuz", 0.0);
                Scalar disp_long = input_reader_->read_double("disp_lon", 0.0);
                Scalar disp_trans = input_reader_->read_double("disp_lat", 0.0);
                
                for (Ordinal kk = 0; kk < config_.n_scalar; ++kk) {
                    gw_scalar_solvers_.push_back(
                        std::make_unique<GwScalarTransportSolver>(
                            *gw_domain_, *gw_active_mesh_, *gw_state_,
                            config_.dt, diff_x, diff_y, diff_z,
                            disp_long, disp_trans, 200.0, 0.0, config_.superbee));
                }
            }
        }
    }
    
    // ========================================================================
    // MAIN INITIALIZATION FUNCTION
    // ========================================================================
    // Call this to perform all initialization steps in order
    void initialize_all() {
        // Step 1: Read configuration from input file
        read_configuration();
        
        // Step 2: Initialize domains
        initialize_domains();
        
        // Step 3: Initialize state variables
        initialize_state_variables();
        
        // Step 4: Build active cell meshes
        build_active_meshes();
        
        // Step 5: Create solvers
        create_solvers();
    }
    
    // ========================================================================
    // GETTERS (for ModelDriver to access initialized objects)
    // ========================================================================
    const Config& get_config() const { return config_; }
    
    SwDomain* get_sw_domain() { return sw_domain_.get(); }
    GwDomain* get_gw_domain() { return gw_domain_.get(); }
    
    SwStateVariables* get_sw_state() { return sw_state_.get(); }
    GwStateVariables* get_gw_state() { return gw_state_.get(); }
    
    ActiveCellMesh* get_sw_active_mesh() { return sw_active_mesh_.get(); }
    ActiveCellMesh* get_gw_active_mesh() { return gw_active_mesh_.get(); }
    
    ShallowWaterSolver* get_sw_solver() { return sw_solver_.get(); }
    GroundwaterSolver* get_gw_solver() { return gw_solver_.get(); }
    
    const std::vector<std::unique_ptr<SwScalarTransportSolver>>& get_sw_scalar_solvers() const {
        return sw_scalar_solvers_;
    }
    
    const std::vector<std::unique_ptr<GwScalarTransportSolver>>& get_gw_scalar_solvers() const {
        return gw_scalar_solvers_;
    }
};

} // namespace Frehg

#endif // FREHG_INITIALIZER_HPP

