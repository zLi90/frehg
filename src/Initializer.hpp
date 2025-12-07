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
#include "BoundaryConditions.hpp"
#include "SourceSinkTerms.hpp"
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
        bool dt_adjust;            // Enable adaptive time stepping
        Scalar dt_min;             // Minimum time step
        Scalar dt_max;             // Maximum time step
        Scalar Co_max;             // Maximum Courant number for groundwater
        
        // Physics flags
        bool sim_shallowwater;
        bool sim_groundwater;
        bool sync_coupling;
        
        // Scalar transport
        Ordinal n_scalar;
        bool baroclinic;
        bool superbee;
        
        // Wind forcing parameters
        bool sim_wind;              // Enable wind forcing
        Scalar wind_speed;          // Initial wind speed (m/s)
        Scalar wind_direction;      // Initial wind direction from north (degrees)
        Scalar north_angle;         // Domain north angle from true north (degrees)
        Scalar rho_air;             // Air density (kg/m³)
        Scalar rho_water;           // Water density (kg/m³)
        Scalar Cw;                  // Wind drag coefficient
        Scalar CwT;                 // Thin layer model coefficient
        Scalar hD;                  // Reference depth for thin layer (m)
        
        // Default constructor
        Config() : NX(0), NY(0), NZ(0), dx(0.0), dy(0.0), dz(0.0), botZ(0.0),
                   dz_incre(1.0), use_mpi(false), mpi_nx(1), mpi_ny(1),
                   nthreads(1), dt(1.0), Tend(100.0), dt_out(10.0), n_monitor(0),
                   dt_adjust(false), dt_min(0.0), dt_max(0.0), Co_max(0.0),
                   sim_shallowwater(true), sim_groundwater(true), sync_coupling(true),
                   n_scalar(0), baroclinic(false), superbee(true),
                   sim_wind(false), wind_speed(0.0), wind_direction(0.0),
                   north_angle(0.0), rho_air(1.225), rho_water(1000.0),
                   Cw(0.0013), CwT(1.0), hD(1.0) {}
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
    std::unique_ptr<SwBoundaryConditionManager> sw_bc_manager_;
    std::unique_ptr<GwBoundaryConditionManager> gw_bc_manager_;
    std::unique_ptr<SwScalarBoundaryConditionManager> sw_scalar_bc_manager_;
    std::unique_ptr<GwScalarBoundaryConditionManager> gw_scalar_bc_manager_;
    std::unique_ptr<SwSourceSinkManager> sw_source_sink_manager_;

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
        config_.dt_adjust = input_reader_->read_bool("dt_adjust", false);
        config_.dt_min = input_reader_->read_double("dt_min", 0.0);
        config_.dt_max = input_reader_->read_double("dt_max", 0.0);
        config_.Co_max = input_reader_->read_double("Co_max", 0.0);
        
        // Read physics flags
        config_.sim_shallowwater = input_reader_->read_bool("sim_shallowwater", true);
        config_.sim_groundwater = input_reader_->read_bool("sim_groundwater", true);
        config_.sync_coupling = input_reader_->read_bool("sync_coupling", true);
        
        // Read scalar transport settings
        config_.n_scalar = input_reader_->read_int("n_scalar", 0);
        config_.baroclinic = input_reader_->read_bool("baroclinic", false);
        config_.superbee = input_reader_->read_bool("superbee", true);
        
        // Read wind forcing settings
        config_.sim_wind = input_reader_->read_bool("sim_wind", false);
        config_.wind_speed = input_reader_->read_double("init_windspd", 0.0);
        config_.wind_direction = input_reader_->read_double("init_winddir", 0.0);
        config_.north_angle = input_reader_->read_double("north_angle", 0.0);
        config_.rho_air = input_reader_->read_double("rhoa", 1.225);
        config_.rho_water = input_reader_->read_double("rhow", 1000.0);
        config_.Cw = input_reader_->read_double("Cw", 0.0013);
        config_.CwT = input_reader_->read_double("CwT", 1.0);
        config_.hD = input_reader_->read_double("hD", 1.0);
        
        // TODO: Add more fields as needed:
        // - Bathymetry file settings
        // - Physical parameters (grav, visc, manning, etc.)
        // - Boundary conditions
        // - Initial conditions
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
            
            // Copy bathymetry to state
            auto h_bottom = Kokkos::create_mirror_view(sw_state_->bottom);
            auto h_bath = Kokkos::create_mirror_view(sw_domain_->bathymetry);
            Kokkos::deep_copy(h_bath, sw_domain_->bathymetry);
            Kokkos::deep_copy(h_bottom, h_bath);
            Kokkos::deep_copy(sw_state_->bottom, h_bottom);
            
            // Initialize surface water state variables
            initialize_sw_state();
        }
        
        if (config_.sim_groundwater && gw_domain_) {
            gw_state_ = std::make_unique<GwStateVariables>(
                gw_domain_->num_cells_3d_total);
            
            // Read and set soil parameters (heterogeneous or homogeneous)
            initialize_soil_parameters();
            
            // Initialize groundwater state variables
            initialize_gw_state();
        }
    }
    
private:
    // ========================================================================
    // INITIALIZE SURFACE WATER STATE VARIABLES
    // ========================================================================
    void initialize_sw_state() {
        // Create host mirrors for initialization
        auto h_pressure = Kokkos::create_mirror_view(sw_state_->pressure);
        auto h_velocity_x = Kokkos::create_mirror_view(sw_state_->velocity_x);
        auto h_velocity_y = Kokkos::create_mirror_view(sw_state_->velocity_y);
        auto h_bottom = Kokkos::create_mirror_view(sw_state_->bottom);
        auto h_active = Kokkos::create_mirror_view(sw_domain_->active_mask);
        
        Kokkos::deep_copy(h_bottom, sw_state_->bottom);
        Kokkos::deep_copy(h_active, sw_domain_->active_mask);
        
        // Read initial condition type for surface elevation
        std::string ic_eta_type = input_reader_->read_string("ic_eta_type", "constant");
        std::vector<Scalar> eta_data;
        
        if (ic_eta_type == "constant") {
            // Constant initial surface elevation
            Scalar init_eta = input_reader_->read_double("init_eta", 0.0);
            eta_data.resize(sw_domain_->num_cells_total, init_eta);
        } else if (ic_eta_type == "file") {
            // Read from file
            std::string eta_filename = input_reader_->read_string("ic_eta_file", "");
            if (eta_filename.empty()) {
                eta_filename = config_.finput + "/surf_ic.asc";
            }
            StructuredGrid eta_grid = read_structured_data(eta_filename);
            
            // Validate dimensions
            if (eta_grid.ncols != sw_domain_->nx || eta_grid.nrows != sw_domain_->ny) {
                throw std::runtime_error("Initial surface elevation file dimensions (" +
                                       std::to_string(eta_grid.ncols) + "x" + 
                                       std::to_string(eta_grid.nrows) + 
                                       ") do not match domain dimensions (" +
                                       std::to_string(sw_domain_->nx) + "x" +
                                       std::to_string(sw_domain_->ny) + ")");
            }
            eta_data = eta_grid.data;
        } else {
            throw std::runtime_error("Unknown ic_eta_type: " + ic_eta_type + 
                                   " (must be 'constant' or 'file')");
        }
        
        // Read initial condition type for water depth (alternative to eta)
        std::string ic_depth_type = input_reader_->read_string("ic_depth_type", "none");
        if (ic_depth_type != "none") {
            std::vector<Scalar> depth_data;
            
            if (ic_depth_type == "constant") {
                Scalar init_depth = input_reader_->read_double("init_depth", 0.0);
                depth_data.resize(sw_domain_->num_cells_total, init_depth);
            } else if (ic_depth_type == "file") {
                std::string depth_filename = input_reader_->read_string("ic_depth_file", "");
                if (depth_filename.empty()) {
                    depth_filename = config_.finput + "/depth_ic.asc";
                }
                StructuredGrid depth_grid = read_structured_data(depth_filename);
                
                if (depth_grid.ncols != sw_domain_->nx || 
                    depth_grid.nrows != sw_domain_->ny) {
                    throw std::runtime_error("Initial depth file dimensions do not match domain");
                }
                depth_data = depth_grid.data;
            }
            
            // Convert depth to surface elevation (eta = bottom + depth)
            for (Ordinal i = 0; i < sw_domain_->num_cells_total; ++i) {
                if (h_active(i) > 0) {
                    eta_data[i] = h_bottom(i) + depth_data[i];
                }
            }
        }
        
        // Read initial velocities
        std::string ic_uv_type = input_reader_->read_string("ic_uv_type", "constant");
        std::vector<Scalar> u_data, v_data;
        
        if (ic_uv_type == "constant") {
            Scalar init_u = input_reader_->read_double("init_u", 0.0);
            Scalar init_v = input_reader_->read_double("init_v", 0.0);
            u_data.resize(sw_domain_->num_cells_total, init_u);
            v_data.resize(sw_domain_->num_cells_total, init_v);
        } else if (ic_uv_type == "file") {
            std::string u_filename = input_reader_->read_string("ic_u_file", "");
            std::string v_filename = input_reader_->read_string("ic_v_file", "");
            if (u_filename.empty()) u_filename = config_.finput + "/uu_ic.asc";
            if (v_filename.empty()) v_filename = config_.finput + "/vv_ic.asc";
            
            StructuredGrid u_grid = read_structured_data(u_filename);
            StructuredGrid v_grid = read_structured_data(v_filename);
            
            if (u_grid.ncols != sw_domain_->nx || u_grid.nrows != sw_domain_->ny ||
                v_grid.ncols != sw_domain_->nx || v_grid.nrows != sw_domain_->ny) {
                throw std::runtime_error("Initial velocity file dimensions do not match domain");
            }
            u_data = u_grid.data;
            v_data = v_grid.data;
        } else {
            throw std::runtime_error("Unknown ic_uv_type: " + ic_uv_type);
        }
        
        // Assign values to state variables
        for (Ordinal i = 0; i < sw_domain_->num_cells_total; ++i) {
            if (h_active(i) > 0) {
                // Ensure eta >= bottom
                h_pressure(i) = std::max(eta_data[i], h_bottom(i));
                h_velocity_x(i) = u_data[i];
                h_velocity_y(i) = v_data[i];
            } else {
                h_pressure(i) = h_bottom(i);
                h_velocity_x(i) = 0.0;
                h_velocity_y(i) = 0.0;
            }
        }
        
        // Copy to device
        Kokkos::deep_copy(sw_state_->pressure, h_pressure);
        Kokkos::deep_copy(sw_state_->velocity_x, h_velocity_x);
        Kokkos::deep_copy(sw_state_->velocity_y, h_velocity_y);
        Kokkos::deep_copy(sw_state_->pressure_old, h_pressure);
        Kokkos::deep_copy(sw_state_->velocity_x_old, h_velocity_x);
        Kokkos::deep_copy(sw_state_->velocity_y_old, h_velocity_y);
    }
    
    // ========================================================================
    // INITIALIZE SOIL PARAMETERS
    // ========================================================================
    // Reads soil parameters from files or sets homogeneous values
    void initialize_soil_parameters() {
        // Create host mirrors for soil parameters
        auto h_Ksx = Kokkos::create_mirror_view(gw_state_->conductivity_sat_x);
        auto h_Ksy = Kokkos::create_mirror_view(gw_state_->conductivity_sat_y);
        auto h_Ksz = Kokkos::create_mirror_view(gw_state_->conductivity_sat_z);
        auto h_porosity = Kokkos::create_mirror_view(gw_state_->water_content_sat);
        auto h_theta_r = Kokkos::create_mirror_view(gw_state_->water_content_res);
        auto h_alpha = Kokkos::create_mirror_view(gw_state_->vg_alpha);
        auto h_n = Kokkos::create_mirror_view(gw_state_->vg_n);
        auto h_m = Kokkos::create_mirror_view(gw_state_->vg_m);
        auto h_ha = Kokkos::create_mirror_view(gw_state_->vg_ha);
        auto h_active_3d = Kokkos::create_mirror_view(gw_domain_->active_mask_3d);
        
        Kokkos::deep_copy(h_active_3d, gw_domain_->active_mask_3d);
        
        // Check if using heterogeneous soil parameters
        bool use_heterogeneous = input_reader_->read_bool("use_heterogeneous_soil", false);
        
        if (use_heterogeneous) {
            // Read soil type map and soil parameter file
            std::string soil_type_map_file = input_reader_->read_string("soil_type_map_file", "");
            std::string soil_param_file = input_reader_->read_string("soil_param_file", "");
            
            if (soil_type_map_file.empty() || soil_param_file.empty()) {
                throw std::runtime_error("Heterogeneous soil enabled but soil_type_map_file or "
                                       "soil_param_file not specified");
            }
            
            // Construct full paths
            std::string soil_type_map_path = config_.finput + "/" + soil_type_map_file;
            std::string soil_param_path = config_.finput + "/" + soil_param_file;
            
            // Check if binary format
            bool binary_soil_map = input_reader_->read_bool("soil_type_map_binary", false);
            
            // Read soil type map (3D)
            std::vector<int> soil_type_map = InputReader::read_3d_soil_type_map(
                soil_type_map_path, gw_domain_->nx, gw_domain_->ny, gw_domain_->nz, binary_soil_map);
            
            // Read soil parameters
            std::map<int, InputReader::SoilParameters> soil_params = 
                InputReader::read_soil_parameter_file(soil_param_path);
            
            // Assign parameters to each cell based on soil type
            for (Ordinal k = 0; k < gw_domain_->nz; ++k) {
                for (Ordinal j = 0; j < gw_domain_->ny; ++j) {
                    for (Ordinal i = 0; i < gw_domain_->nx; ++i) {
                        Ordinal idx = i + j * gw_domain_->nx + k * gw_domain_->nx * gw_domain_->ny;
                        
                        if (h_active_3d(idx) > 0) {
                            int soil_type = soil_type_map[idx];
                            
                            // Check if soil type exists in parameter map
                            if (soil_params.find(soil_type) == soil_params.end()) {
                                throw std::runtime_error("Soil type " + std::to_string(soil_type) + 
                                                       " not found in soil parameter file at "
                                                       "cell (i=" + std::to_string(i) + 
                                                       ", j=" + std::to_string(j) + 
                                                       ", k=" + std::to_string(k) + ")");
                            }
                            
                            const auto& params = soil_params[soil_type];
                            
                            // Set parameters
                            h_Ksx(idx) = params.Ksx;
                            h_Ksy(idx) = params.Ksy;
                            h_Ksz(idx) = params.Ksz;
                            h_porosity(idx) = params.porosity;
                            h_theta_r(idx) = params.theta_r;
                            h_alpha(idx) = params.alpha;
                            h_n(idx) = params.n;
                            h_m(idx) = 1.0 - 1.0 / params.n;  // m = 1 - 1/n
                            h_ha(idx) = params.ha;
                        }
                    }
                }
            }
        } else {
            // Homogeneous soil parameters
            Scalar Ksx = input_reader_->read_double("Ksx", 1.0e-4);
            Scalar Ksy = input_reader_->read_double("Ksy", 1.0e-4);
            Scalar Ksz = input_reader_->read_double("Ksz", 1.0e-4);
            Scalar porosity = input_reader_->read_double("porosity", 0.4);
            Scalar theta_r = input_reader_->read_double("theta_r", 0.05);
            Scalar alpha = input_reader_->read_double("vg_alpha", 0.01);
            Scalar n = input_reader_->read_double("vg_n", 1.5);
            Scalar ha = input_reader_->read_double("vg_ha", 0.0);
            Scalar m = 1.0 - 1.0 / n;
            
            // Set homogeneous values for all active cells
            for (Ordinal i = 0; i < gw_domain_->num_cells_3d_total; ++i) {
                if (h_active_3d(i) > 0) {
                    h_Ksx(i) = Ksx;
                    h_Ksy(i) = Ksy;
                    h_Ksz(i) = Ksz;
                    h_porosity(i) = porosity;
                    h_theta_r(i) = theta_r;
                    h_alpha(i) = alpha;
                    h_n(i) = n;
                    h_m(i) = m;
                    h_ha(i) = ha;
                }
            }
        }
        
        // Copy to device
        Kokkos::deep_copy(gw_state_->conductivity_sat_x, h_Ksx);
        Kokkos::deep_copy(gw_state_->conductivity_sat_y, h_Ksy);
        Kokkos::deep_copy(gw_state_->conductivity_sat_z, h_Ksz);
        Kokkos::deep_copy(gw_state_->water_content_sat, h_porosity);
        Kokkos::deep_copy(gw_state_->water_content_res, h_theta_r);
        Kokkos::deep_copy(gw_state_->vg_alpha, h_alpha);
        Kokkos::deep_copy(gw_state_->vg_n, h_n);
        Kokkos::deep_copy(gw_state_->vg_m, h_m);
        Kokkos::deep_copy(gw_state_->vg_ha, h_ha);
    }
    
    // ========================================================================
    // INITIALIZE GROUNDWATER STATE VARIABLES
    // ========================================================================
    void initialize_gw_state() {
        // Create host mirrors
        auto h_pressure = Kokkos::create_mirror_view(gw_state_->pressure);
        auto h_water_content = Kokkos::create_mirror_view(gw_state_->water_content);
        auto h_active_3d = Kokkos::create_mirror_view(gw_domain_->active_mask_3d);
        auto h_layer_bottoms = Kokkos::create_mirror_view(gw_domain_->layer_bottoms);
        auto h_layer_thickness = Kokkos::create_mirror_view(gw_domain_->layer_thickness);
        
        Kokkos::deep_copy(h_active_3d, gw_domain_->active_mask_3d);
        Kokkos::deep_copy(h_layer_bottoms, gw_domain_->layer_bottoms);
        Kokkos::deep_copy(h_layer_thickness, gw_domain_->layer_thickness);
        
        // Create bathymetry mirror only if surface domain exists
        View1D<Scalar> h_bath;
        if (sw_domain_) {
            h_bath = Kokkos::create_mirror_view(sw_domain_->bathymetry);
            Kokkos::deep_copy(h_bath, sw_domain_->bathymetry);
        }
        
        // Initialize all cells to zero first
        for (Ordinal i = 0; i < gw_domain_->num_cells_3d_total; ++i) {
            h_pressure(i) = 0.0;
            h_water_content(i) = 0.0;
        }
        
        // Read initial condition type for groundwater head
        std::string ic_h_type = input_reader_->read_string("ic_h_type", "constant");
        std::vector<Scalar> h_data;
        bool h_from_file = false;
        
        if (ic_h_type == "constant") {
            Scalar init_h = input_reader_->read_double("init_h", 0.0);
            h_data.resize(gw_domain_->num_cells_3d_total, init_h);
        } else if (ic_h_type == "file") {
            std::string h_filename = input_reader_->read_string("ic_h_file", "");
            if (h_filename.empty()) {
                h_filename = config_.finput + "/head_ic.asc";
            }
            // Note: For 3D, we might need a different format or multiple files
            // For now, assume a 2D file that applies to all layers
            StructuredGrid h_grid = read_structured_data(h_filename);
            if (h_grid.ncols != gw_domain_->nx || h_grid.nrows != gw_domain_->ny) {
                throw std::runtime_error("Initial head file dimensions do not match domain");
            }
            // Expand 2D to 3D (same value for all layers at each (i,j))
            h_data.resize(gw_domain_->num_cells_3d_total);
            for (Ordinal k = 0; k < gw_domain_->nz; ++k) {
                for (Ordinal j = 0; j < gw_domain_->ny; ++j) {
                    for (Ordinal i = 0; i < gw_domain_->nx; ++i) {
                        Ordinal idx_2d = i + j * gw_domain_->nx;
                        Ordinal idx_3d = i + j * gw_domain_->nx + k * gw_domain_->nx * gw_domain_->ny;
                        h_data[idx_3d] = h_grid.data[idx_2d];
                    }
                }
            }
            h_from_file = true;
        } else if (ic_h_type == "water_table") {
            // Initialize from water table elevation
            std::string ic_wt_type = input_reader_->read_string("ic_wt_type", "constant");
            h_data.resize(gw_domain_->num_cells_3d_total);
            
            if (ic_wt_type == "constant") {
                Scalar init_wt = input_reader_->read_double("init_wt", 0.0);
                bool wt_relative = input_reader_->read_bool("init_wt_relative", false);
                
                if (wt_relative && sw_domain_) {
                    // Water table relative to surface (bathymetry)
                    for (Ordinal j = 0; j < gw_domain_->ny; ++j) {
                        for (Ordinal i = 0; i < gw_domain_->nx; ++i) {
                            Ordinal sw_idx = i + j * gw_domain_->nx;
                            Scalar surface_elev = h_bath(sw_idx);
                            Scalar wt_elev = surface_elev - init_wt;
                            
                            // Assign to all layers
                            for (Ordinal k = 0; k < gw_domain_->nz; ++k) {
                                Ordinal gw_idx = i + j * gw_domain_->nx + 
                                               k * gw_domain_->nx * gw_domain_->ny;
                                Scalar cell_bottom = h_layer_bottoms(k);
                                Scalar cell_thickness = h_layer_thickness(k);
                                Scalar cell_center = cell_bottom + 0.5 * cell_thickness;
                                
                                // Compute head assuming hydrostatic: h = wt_elev - cell_center
                                h_data[gw_idx] = wt_elev - cell_center;
                            }
                        }
                    }
                } else {
                    // Fixed water table elevation
                    for (Ordinal k = 0; k < gw_domain_->nz; ++k) {
                        for (Ordinal j = 0; j < gw_domain_->ny; ++j) {
                            for (Ordinal i = 0; i < gw_domain_->nx; ++i) {
                                Ordinal gw_idx = i + j * gw_domain_->nx + 
                                               k * gw_domain_->nx * gw_domain_->ny;
                                Scalar cell_bottom = h_layer_bottoms(k);
                                Scalar cell_thickness = h_layer_thickness(k);
                                Scalar cell_center = cell_bottom + 0.5 * cell_thickness;
                                h_data[gw_idx] = init_wt - cell_center;
                            }
                        }
                    }
                }
            } else if (ic_wt_type == "file") {
                std::string wt_filename = input_reader_->read_string("ic_wt_file", "");
                if (wt_filename.empty()) {
                    wt_filename = config_.finput + "/wt_ic.asc";
                }
                StructuredGrid wt_grid = read_structured_data(wt_filename);
                if (wt_grid.ncols != gw_domain_->nx || wt_grid.nrows != gw_domain_->ny) {
                    throw std::runtime_error("Initial water table file dimensions do not match domain");
                }
                
                // Convert water table elevation to head for each cell
                for (Ordinal k = 0; k < gw_domain_->nz; ++k) {
                    for (Ordinal j = 0; j < gw_domain_->ny; ++j) {
                        for (Ordinal i = 0; i < gw_domain_->nx; ++i) {
                            Ordinal idx_2d = i + j * gw_domain_->nx;
                            Ordinal idx_3d = i + j * gw_domain_->nx + 
                                          k * gw_domain_->nx * gw_domain_->ny;
                            Scalar wt_elev = wt_grid.data[idx_2d];
                            Scalar cell_bottom = h_layer_bottoms(k);
                            Scalar cell_thickness = h_layer_thickness(k);
                            Scalar cell_center = cell_bottom + 0.5 * cell_thickness;
                            h_data[idx_3d] = wt_elev - cell_center;
                        }
                    }
                }
            } else {
                throw std::runtime_error("Unknown ic_wt_type: " + ic_wt_type);
            }
        } else {
            throw std::runtime_error("Unknown ic_h_type: " + ic_h_type);
        }
        
        // Read initial condition type for water content
        std::string ic_wc_type = input_reader_->read_string("ic_wc_type", "from_head");
        std::vector<Scalar> wc_data;
        bool wc_from_file = false;
        
        if (ic_wc_type == "constant") {
            Scalar init_wc = input_reader_->read_double("init_wc", 0.0);
            wc_data.resize(gw_domain_->num_cells_3d_total, init_wc);
            wc_from_file = true;
        } else if (ic_wc_type == "file") {
            std::string wc_filename = input_reader_->read_string("ic_wc_file", "");
            if (wc_filename.empty()) {
                wc_filename = config_.finput + "/moisture_ic.asc";
            }
            StructuredGrid wc_grid = read_structured_data(wc_filename);
            if (wc_grid.ncols != gw_domain_->nx || wc_grid.nrows != gw_domain_->ny) {
                throw std::runtime_error("Initial water content file dimensions do not match domain");
            }
            // Expand 2D to 3D
            wc_data.resize(gw_domain_->num_cells_3d_total);
            for (Ordinal k = 0; k < gw_domain_->nz; ++k) {
                for (Ordinal j = 0; j < gw_domain_->ny; ++j) {
                    for (Ordinal i = 0; i < gw_domain_->nx; ++i) {
                        Ordinal idx_2d = i + j * gw_domain_->nx;
                        Ordinal idx_3d = i + j * gw_domain_->nx + 
                                       k * gw_domain_->nx * gw_domain_->ny;
                        wc_data[idx_3d] = wc_grid.data[idx_2d];
                    }
                }
            }
            wc_from_file = true;
        } else if (ic_wc_type == "from_head") {
            // Will compute from head after head is set
            wc_from_file = false;
        } else {
            throw std::runtime_error("Unknown ic_wc_type: " + ic_wc_type);
        }
        
        // Assign head values to active cells
        for (Ordinal i = 0; i < gw_domain_->num_cells_3d_total; ++i) {
            if (h_active_3d(i) > 0) {
                h_pressure(i) = h_data[i];
            }
        }
        
        // Assign water content values
        if (wc_from_file) {
            for (Ordinal i = 0; i < gw_domain_->num_cells_3d_total; ++i) {
                if (h_active_3d(i) > 0) {
                    h_water_content(i) = wc_data[i];
                }
            }
        } else {
            // TODO: Compute water content from head using soil properties
            // For now, set to a default value
            // This will be implemented when soil property functions are available
            Scalar default_wc = input_reader_->read_double("init_wc_default", 0.1);
            for (Ordinal i = 0; i < gw_domain_->num_cells_3d_total; ++i) {
                if (h_active_3d(i) > 0) {
                    h_water_content(i) = default_wc;
                }
            }
        }
        
        // Copy to device
        Kokkos::deep_copy(gw_state_->pressure, h_pressure);
        Kokkos::deep_copy(gw_state_->water_content, h_water_content);
        Kokkos::deep_copy(gw_state_->pressure_old, h_pressure);
        Kokkos::deep_copy(gw_state_->water_content_old, h_water_content);
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
            
            Scalar min_depth = input_reader_->read_double("min_dept", 1.0e-6);
            Scalar water_threshold = input_reader_->read_double("wtfh", 0.01);
            
            sw_solver_ = std::make_unique<ShallowWaterSolver>(
                *sw_domain_, *sw_active_mesh_, *sw_state_,
                config_.dt, grav, manning_n, visc_x, visc_y,
                min_depth, water_threshold,
                PreconditionerType::JACOBI, 1.0e-8, 10000,
                sw_bc_manager_.get(), sw_source_sink_manager_.get());
            
            // Configure wind forcing
            if (config_.sim_wind) {
                sw_solver_->set_wind_enabled(true);
                sw_solver_->set_wind_parameters(
                    config_.rho_air, config_.rho_water, config_.Cw,
                    config_.CwT, config_.hD, config_.north_angle);
                sw_solver_->set_wind_conditions(config_.wind_speed, config_.wind_direction);
            }
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
                config_.dt, Ss, 1.0e-6, use_corrector, follow_terrain, baroclinic,
                PreconditionerType::JACOBI, 1.0e-8, 10000,
                gw_bc_manager_.get());
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
                            config_.dt, kk, diff_x, diff_y, 200.0, 0.0, config_.superbee,
                            sw_scalar_bc_manager_.get(), sw_source_sink_manager_.get()));
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
                            config_.dt, kk, diff_x, diff_y, diff_z,
                            disp_long, disp_trans, 200.0, 0.0, config_.superbee,
                            gw_scalar_bc_manager_.get()));
                }
            }
        }
    }
    
    // ========================================================================
    // INITIALIZE BOUNDARY CONDITIONS
    // ========================================================================
    void initialize_boundary_conditions() {
        if (config_.sim_shallowwater && sw_domain_) {
            sw_bc_manager_ = std::make_unique<SwBoundaryConditionManager>(sw_domain_.get());
            sw_bc_manager_->read_from_input(input_reader_.get(), config_.finput);
            
            // Initialize source/sink terms
            sw_source_sink_manager_ = std::make_unique<SwSourceSinkManager>(sw_domain_.get());
            sw_source_sink_manager_->read_from_input(input_reader_.get(), config_.finput);
        }
        
        if (config_.sim_groundwater && gw_domain_) {
            gw_bc_manager_ = std::make_unique<GwBoundaryConditionManager>(gw_domain_.get());
            gw_bc_manager_->read_from_input(input_reader_.get(), config_.finput);
        }
        
        // Initialize scalar transport boundary conditions
        if (config_.n_scalar > 0) {
            if (config_.sim_shallowwater && sw_domain_) {
                sw_scalar_bc_manager_ = std::make_unique<SwScalarBoundaryConditionManager>(
                    sw_domain_.get(), config_.n_scalar);
                sw_scalar_bc_manager_->read_from_input(input_reader_.get(), config_.finput);
            }
            
            if (config_.sim_groundwater && gw_domain_) {
                gw_scalar_bc_manager_ = std::make_unique<GwScalarBoundaryConditionManager>(
                    gw_domain_.get(), config_.n_scalar);
                gw_scalar_bc_manager_->read_from_input(input_reader_.get(), config_.finput);
            }
        }
    }
    
public:
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
        
        // Step 5: Initialize boundary conditions
        initialize_boundary_conditions();
        
        // Step 6: Create solvers
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
    
    SwBoundaryConditionManager* get_sw_bc_manager() { return sw_bc_manager_.get(); }
    GwBoundaryConditionManager* get_gw_bc_manager() { return gw_bc_manager_.get(); }
    SwScalarBoundaryConditionManager* get_sw_scalar_bc_manager() { return sw_scalar_bc_manager_.get(); }
    GwScalarBoundaryConditionManager* get_gw_scalar_bc_manager() { return gw_scalar_bc_manager_.get(); }
    SwSourceSinkManager* get_sw_source_sink_manager() { return sw_source_sink_manager_.get(); }
};

} // namespace Frehg

#endif // FREHG_INITIALIZER_HPP

