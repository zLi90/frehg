#ifndef FREHG_BOUNDARY_CONDITIONS_HPP
#define FREHG_BOUNDARY_CONDITIONS_HPP

#include "define.hpp"
#include "Domain.hpp"
#include "InputReader.hpp"
#include <vector>
#include <string>
#include <memory>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <numeric>

namespace Frehg {

// ============================================================================
//                      BOUNDARY CONDITION TYPES
// ============================================================================
enum class SwBcType {
    FREE_SURFACE_ELEVATION,  // Prescribed surface elevation (eta)
    WATER_DEPTH,             // Prescribed water depth
    FLOW_RATE,               // Prescribed flow rate (discharge)
    FREE_OUTFLOW             // Free outflow (zero gradient)
};

// ============================================================================
//                      BOUNDARY CONDITION VALUE TYPE
// ============================================================================
enum class BcValueType {
    CONSTANT,                // Constant value
    TIME_SERIES              // Time-varying from file
};

// ============================================================================
//                      SURFACE WATER BOUNDARY CONDITION
// ============================================================================
class SwBoundaryCondition {
public:
    SwBcType type;                    // Type of boundary condition
    BcValueType value_type;            // Constant or time series
    
    // Location range (in grid coordinates: i, j indices)
    Ordinal i_min, i_max;             // X-direction range
    Ordinal j_min, j_max;             // Y-direction range
    
    // Cell indices where BC is applied (computed from range)
    std::vector<Ordinal> cell_indices;
    
    // Constant value (if value_type == CONSTANT)
    Scalar constant_value;
    
    // Time series data (if value_type == TIME_SERIES)
    std::vector<Scalar> time_values;   // Time points
    std::vector<Scalar> data_values;   // BC values at each time point
    std::string time_series_file;      // File path for time series
    
    // Constructor
    SwBoundaryCondition() 
        : type(SwBcType::FREE_SURFACE_ELEVATION),
          value_type(BcValueType::CONSTANT),
          i_min(0), i_max(0), j_min(0), j_max(0),
          constant_value(0.0) {}
    
    // Get current BC value at given time (interpolates if time series)
    Scalar get_value(Scalar current_time) const {
        if (value_type == BcValueType::CONSTANT) {
            return constant_value;
        } else {
            // Linear interpolation for time series
            if (time_values.empty()) {
                return constant_value;
            }
            
            // Clamp to bounds
            if (current_time <= time_values.front()) {
                return data_values.front();
            }
            if (current_time >= time_values.back()) {
                return data_values.back();
            }
            
            // Find interpolation interval
            for (size_t i = 0; i < time_values.size() - 1; ++i) {
                if (time_values[i] <= current_time && current_time <= time_values[i+1]) {
                    Scalar t0 = time_values[i];
                    Scalar t1 = time_values[i+1];
                    Scalar v0 = data_values[i];
                    Scalar v1 = data_values[i+1];
                    Scalar alpha = (current_time - t0) / (t1 - t0);
                    return v0 + alpha * (v1 - v0);
                }
            }
            
            return data_values.back();
        }
    }
    
    // Load time series from file
    void load_time_series(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open boundary condition time series file: " + filename);
        }
        
        time_values.clear();
        data_values.clear();
        
        std::string line;
        while (std::getline(file, line)) {
            // Skip empty lines and comments
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            if (line.empty() || line[0] == '#') continue;
            
            std::istringstream iss(line);
            Scalar t, v;
            if (iss >> t >> v) {
                time_values.push_back(t);
                data_values.push_back(v);
            }
        }
        
        file.close();
        
        if (time_values.empty()) {
            throw std::runtime_error("No valid data found in boundary condition file: " + filename);
        }
        
        // Sort by time (in case file is not sorted)
        std::vector<size_t> indices(time_values.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), 
                 [&](size_t a, size_t b) { return time_values[a] < time_values[b]; });
        
        std::vector<Scalar> sorted_times, sorted_values;
        for (size_t idx : indices) {
            sorted_times.push_back(time_values[idx]);
            sorted_values.push_back(data_values[idx]);
        }
        time_values = sorted_times;
        data_values = sorted_values;
    }
};

// ============================================================================
//                      BOUNDARY CONDITION MANAGER
// ============================================================================
class SwBoundaryConditionManager {
private:
    SwDomain* domain_;
    std::vector<SwBoundaryCondition> boundary_conditions_;
    
public:
    SwBoundaryConditionManager(SwDomain* domain) : domain_(domain) {}
    
    // Add a boundary condition
    void add_boundary_condition(const SwBoundaryCondition& bc) {
        boundary_conditions_.push_back(bc);
        // Identify cells for this BC
        identify_boundary_cells(boundary_conditions_.back());
    }
    
    // Read boundary conditions from input file
    void read_from_input(InputReader* input_reader, const std::string& input_dir) {
        // Read number of boundary conditions
        Ordinal n_bc = input_reader->read_int("n_sw_bc", 0);
        
        for (Ordinal i = 0; i < n_bc; ++i) {
            SwBoundaryCondition bc;
            
            // Read BC type
            std::string type_str = input_reader->read_string("sw_bc_type_" + std::to_string(i), "eta");
            if (type_str == "eta" || type_str == "surface_elevation") {
                bc.type = SwBcType::FREE_SURFACE_ELEVATION;
            } else if (type_str == "depth") {
                bc.type = SwBcType::WATER_DEPTH;
            } else if (type_str == "flow_rate" || type_str == "discharge") {
                bc.type = SwBcType::FLOW_RATE;
            } else if (type_str == "free_outflow" || type_str == "outflow") {
                bc.type = SwBcType::FREE_OUTFLOW;
            } else {
                throw std::runtime_error("Unknown surface water BC type: " + type_str);
            }
            
            // Read location range (in grid coordinates)
            std::vector<int> loc_range = input_reader->read_int_array("sw_bc_loc_" + std::to_string(i), 
                                                                      std::vector<int>{0, 0, 0, 0});
            if (loc_range.size() != 4) {
                throw std::runtime_error("sw_bc_loc must have 4 values: [i_min, i_max, j_min, j_max]");
            }
            bc.i_min = static_cast<Ordinal>(loc_range[0]);
            bc.i_max = static_cast<Ordinal>(loc_range[1]);
            bc.j_min = static_cast<Ordinal>(loc_range[2]);
            bc.j_max = static_cast<Ordinal>(loc_range[3]);
            
            // Validate range
            if (bc.i_min > bc.i_max || bc.j_min > bc.j_max) {
                throw std::runtime_error("Invalid BC location range: i_min > i_max or j_min > j_max");
            }
            if (bc.i_max >= domain_->nx || bc.j_max >= domain_->ny) {
                throw std::runtime_error("BC location range exceeds domain dimensions");
            }
            
            // Read value type
            std::string value_type_str = input_reader->read_string("sw_bc_value_type_" + std::to_string(i), "constant");
            if (value_type_str == "constant") {
                bc.value_type = BcValueType::CONSTANT;
                bc.constant_value = input_reader->read_double("sw_bc_value_" + std::to_string(i), 0.0);
            } else if (value_type_str == "time_series" || value_type_str == "file") {
                bc.value_type = BcValueType::TIME_SERIES;
                bc.time_series_file = input_reader->read_string("sw_bc_file_" + std::to_string(i), "");
                if (bc.time_series_file.empty()) {
                    bc.time_series_file = input_dir + "/sw_bc_" + std::to_string(i) + ".txt";
                }
                bc.load_time_series(bc.time_series_file);
            } else {
                throw std::runtime_error("Unknown BC value type: " + value_type_str);
            }
            
            // Identify boundary cells
            identify_boundary_cells(bc);
            
            // Add to list
            boundary_conditions_.push_back(bc);
        }
    }
    
    // Get all boundary conditions
    const std::vector<SwBoundaryCondition>& get_boundary_conditions() const {
        return boundary_conditions_;
    }
    
    // Get boundary conditions for a specific cell
    std::vector<const SwBoundaryCondition*> get_bc_for_cell(Ordinal cell_idx) const {
        std::vector<const SwBoundaryCondition*> result;
        for (const auto& bc : boundary_conditions_) {
            if (std::find(bc.cell_indices.begin(), bc.cell_indices.end(), cell_idx) 
                != bc.cell_indices.end()) {
                result.push_back(&bc);
            }
        }
        return result;
    }
    
    // Check if a cell has a boundary condition
    bool has_boundary_condition(Ordinal cell_idx) const {
        for (const auto& bc : boundary_conditions_) {
            if (std::find(bc.cell_indices.begin(), bc.cell_indices.end(), cell_idx) 
                != bc.cell_indices.end()) {
                return true;
            }
        }
        return false;
    }
    
private:
    // Identify boundary cells from location range
    // A cell is considered a boundary cell if it is:
    // 1. Within the specified range
    // 2. Active
    // 3. At the domain boundary OR has at least one inactive neighbor
    // This handles both regular rectangular domains and irregular domains
    void identify_boundary_cells(SwBoundaryCondition& bc) {
        bc.cell_indices.clear();
        
        // Get active mask
        auto h_active = Kokkos::create_mirror_view(domain_->active_mask);
        Kokkos::deep_copy(h_active, domain_->active_mask);
        
        // Iterate through all cells in the range
        for (Ordinal j = bc.j_min; j <= bc.j_max; ++j) {
            for (Ordinal i = bc.i_min; i <= bc.i_max; ++i) {
                Ordinal cell_idx = i + j * domain_->nx;
                
                // Check if cell is active
                if (h_active(cell_idx) <= 0) continue;
                
                // Check if cell is at domain boundary or has inactive neighbor
                bool is_boundary = false;
                
                // Check if at domain edge (physical boundary)
                if (i == 0 || i == domain_->nx - 1 || 
                    j == 0 || j == domain_->ny - 1) {
                    is_boundary = true;
                } else {
                    // Check if any neighbor is inactive (boundary in irregular domain)
                    Ordinal left = (i - 1) + j * domain_->nx;
                    Ordinal right = (i + 1) + j * domain_->nx;
                    Ordinal bottom = i + (j - 1) * domain_->nx;
                    Ordinal top = i + (j + 1) * domain_->nx;
                    
                    if (h_active(left) <= 0 || h_active(right) <= 0 ||
                        h_active(bottom) <= 0 || h_active(top) <= 0) {
                        is_boundary = true;
                    }
                }
                
                if (is_boundary) {
                    bc.cell_indices.push_back(cell_idx);
                }
            }
        }
    }
};

// ============================================================================
//                      GROUNDWATER BOUNDARY CONDITION TYPES
// ============================================================================
enum class GwBcType {
    FIXED_HEAD,          // Prescribed hydraulic head (Dirichlet)
    FIXED_FLUX,          // Prescribed flux (Neumann)
    NO_FLOW,             // Zero flux (no flow boundary)
    FREE_OUTFLOW         // Free outflow (zero gradient)
};

// ============================================================================
//                      GROUNDWATER BOUNDARY CONDITION
// ============================================================================
class GwBoundaryCondition {
public:
    GwBcType type;                    // Type of boundary condition
    BcValueType value_type;            // Constant or time series
    
    // Location range (in grid coordinates: i, j, k indices)
    Ordinal i_min, i_max;             // X-direction range
    Ordinal j_min, j_max;             // Y-direction range
    Ordinal k_min, k_max;             // Z-direction range (optional, -1 means all layers)
    
    // Boundary face (if specified, applies to entire face)
    // If not specified, applies to cells within range
    std::string face;                  // "x+", "x-", "y+", "y-", "z+", "z-" or "cell"
    
    // Cell indices where BC is applied (computed from range)
    std::vector<Ordinal> cell_indices;
    
    // Constant value (if value_type == CONSTANT)
    Scalar constant_value;
    
    // Time series data (if value_type == TIME_SERIES)
    std::vector<Scalar> time_values;   // Time points
    std::vector<Scalar> data_values;   // BC values at each time point
    std::string time_series_file;      // File path for time series
    
    // Constructor
    GwBoundaryCondition() 
        : type(GwBcType::NO_FLOW),
          value_type(BcValueType::CONSTANT),
          i_min(0), i_max(0), j_min(0), j_max(0), k_min(-1), k_max(-1),
          face("cell"), constant_value(0.0) {}
    
    // Get current BC value at given time (interpolates if time series)
    Scalar get_value(Scalar current_time) const {
        if (value_type == BcValueType::CONSTANT) {
            return constant_value;
        } else {
            // Linear interpolation for time series (same as SwBoundaryCondition)
            if (time_values.empty()) {
                return constant_value;
            }
            
            if (current_time <= time_values.front()) {
                return data_values.front();
            }
            if (current_time >= time_values.back()) {
                return data_values.back();
            }
            
            for (size_t i = 0; i < time_values.size() - 1; ++i) {
                if (time_values[i] <= current_time && current_time <= time_values[i+1]) {
                    Scalar t0 = time_values[i];
                    Scalar t1 = time_values[i+1];
                    Scalar v0 = data_values[i];
                    Scalar v1 = data_values[i+1];
                    Scalar alpha = (current_time - t0) / (t1 - t0);
                    return v0 + alpha * (v1 - v0);
                }
            }
            
            return data_values.back();
        }
    }
    
    // Load time series from file (same as SwBoundaryCondition)
    void load_time_series(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open boundary condition time series file: " + filename);
        }
        
        time_values.clear();
        data_values.clear();
        
        std::string line;
        while (std::getline(file, line)) {
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            if (line.empty() || line[0] == '#') continue;
            
            std::istringstream iss(line);
            Scalar t, v;
            if (iss >> t >> v) {
                time_values.push_back(t);
                data_values.push_back(v);
            }
        }
        
        file.close();
        
        if (time_values.empty()) {
            throw std::runtime_error("No valid data found in boundary condition file: " + filename);
        }
        
        // Sort by time
        std::vector<size_t> indices(time_values.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), 
                 [&](size_t a, size_t b) { return time_values[a] < time_values[b]; });
        
        std::vector<Scalar> sorted_times, sorted_values;
        for (size_t idx : indices) {
            sorted_times.push_back(time_values[idx]);
            sorted_values.push_back(data_values[idx]);
        }
        time_values = sorted_times;
        data_values = sorted_values;
    }
};

// ============================================================================
//                      GROUNDWATER BOUNDARY CONDITION MANAGER
// ============================================================================
class GwBoundaryConditionManager {
private:
    GwDomain* domain_;
    std::vector<GwBoundaryCondition> boundary_conditions_;
    
public:
    GwBoundaryConditionManager(GwDomain* domain) : domain_(domain) {}
    
    // Add a boundary condition
    void add_boundary_condition(const GwBoundaryCondition& bc) {
        boundary_conditions_.push_back(bc);
        identify_boundary_cells(boundary_conditions_.back());
    }
    
    // Read boundary conditions from input file
    void read_from_input(InputReader* input_reader, const std::string& input_dir) {
        // Read number of boundary conditions
        Ordinal n_bc = input_reader->read_int("n_gw_bc", 0);
        
        for (Ordinal i = 0; i < n_bc; ++i) {
            GwBoundaryCondition bc;
            
            // Read BC type
            std::string type_str = input_reader->read_string("gw_bc_type_" + std::to_string(i), "no_flow");
            if (type_str == "fixed_head" || type_str == "head") {
                bc.type = GwBcType::FIXED_HEAD;
            } else if (type_str == "fixed_flux" || type_str == "flux") {
                bc.type = GwBcType::FIXED_FLUX;
            } else if (type_str == "no_flow") {
                bc.type = GwBcType::NO_FLOW;
            } else if (type_str == "free_outflow") {
                bc.type = GwBcType::FREE_OUTFLOW;
            } else {
                throw std::runtime_error("Unknown groundwater BC type: " + type_str);
            }
            
            // Read boundary face or location range
            bc.face = input_reader->read_string("gw_bc_face_" + std::to_string(i), "cell");
            
            if (bc.face == "cell") {
                // Cell-by-cell BC: read location range
                std::vector<int> loc_range = input_reader->read_int_array("gw_bc_loc_" + std::to_string(i), 
                                                                          std::vector<int>{0, 0, 0, 0, -1, -1});
                if (loc_range.size() < 4) {
                    throw std::runtime_error("gw_bc_loc must have at least 4 values: [i_min, i_max, j_min, j_max] or 6 values with [k_min, k_max]");
                }
                bc.i_min = static_cast<Ordinal>(loc_range[0]);
                bc.i_max = static_cast<Ordinal>(loc_range[1]);
                bc.j_min = static_cast<Ordinal>(loc_range[2]);
                bc.j_max = static_cast<Ordinal>(loc_range[3]);
                if (loc_range.size() >= 6) {
                    bc.k_min = static_cast<Ordinal>(loc_range[4]);
                    bc.k_max = static_cast<Ordinal>(loc_range[5]);
                } else {
                    bc.k_min = -1;  // All layers
                    bc.k_max = -1;
                }
                
                // Validate range
                if (bc.i_min > bc.i_max || bc.j_min > bc.j_max) {
                    throw std::runtime_error("Invalid BC location range");
                }
                if (bc.i_max >= domain_->nx || bc.j_max >= domain_->ny) {
                    throw std::runtime_error("BC location range exceeds domain dimensions");
                }
            } else {
                // Face BC: validate face name
                if (bc.face != "x+" && bc.face != "x-" && 
                    bc.face != "y+" && bc.face != "y-" &&
                    bc.face != "z+" && bc.face != "z-") {
                    throw std::runtime_error("Invalid boundary face: " + bc.face);
                }
            }
            
            // Read value type
            std::string value_type_str = input_reader->read_string("gw_bc_value_type_" + std::to_string(i), "constant");
            if (value_type_str == "constant") {
                bc.value_type = BcValueType::CONSTANT;
                bc.constant_value = input_reader->read_double("gw_bc_value_" + std::to_string(i), 0.0);
            } else if (value_type_str == "time_series" || value_type_str == "file") {
                bc.value_type = BcValueType::TIME_SERIES;
                bc.time_series_file = input_reader->read_string("gw_bc_file_" + std::to_string(i), "");
                if (bc.time_series_file.empty()) {
                    bc.time_series_file = input_dir + "/gw_bc_" + std::to_string(i) + ".txt";
                }
                bc.load_time_series(bc.time_series_file);
            } else {
                throw std::runtime_error("Unknown BC value type: " + value_type_str);
            }
            
            // Identify boundary cells
            identify_boundary_cells(bc);
            
            // Add to list
            boundary_conditions_.push_back(bc);
        }
    }
    
    // Get all boundary conditions
    const std::vector<GwBoundaryCondition>& get_boundary_conditions() const {
        return boundary_conditions_;
    }
    
    // Get boundary conditions for a specific cell
    std::vector<const GwBoundaryCondition*> get_bc_for_cell(Ordinal cell_idx) const {
        std::vector<const GwBoundaryCondition*> result;
        for (const auto& bc : boundary_conditions_) {
            if (std::find(bc.cell_indices.begin(), bc.cell_indices.end(), cell_idx) 
                != bc.cell_indices.end()) {
                result.push_back(&bc);
            }
        }
        return result;
    }
    
    // Check if a cell has a boundary condition
    bool has_boundary_condition(Ordinal cell_idx) const {
        for (const auto& bc : boundary_conditions_) {
            if (std::find(bc.cell_indices.begin(), bc.cell_indices.end(), cell_idx) 
                != bc.cell_indices.end()) {
                return true;
            }
        }
        return false;
    }
    
private:
    // Identify boundary cells from location range or face
    void identify_boundary_cells(GwBoundaryCondition& bc) {
        bc.cell_indices.clear();
        
        // Get active mask
        auto h_active_3d = Kokkos::create_mirror_view(domain_->active_mask_3d);
        Kokkos::deep_copy(h_active_3d, domain_->active_mask_3d);
        
        if (bc.face != "cell") {
            // Face boundary: identify all cells on that face
            if (bc.face == "x+") {
                // Right face: i = nx - 1
                for (Ordinal k = 0; k < domain_->nz; ++k) {
                    for (Ordinal j = 0; j < domain_->ny; ++j) {
                        Ordinal cell_idx = (domain_->nx - 1) + j * domain_->nx + k * domain_->nx * domain_->ny;
                        if (h_active_3d(cell_idx) > 0) {
                            bc.cell_indices.push_back(cell_idx);
                        }
                    }
                }
            } else if (bc.face == "x-") {
                // Left face: i = 0
                for (Ordinal k = 0; k < domain_->nz; ++k) {
                    for (Ordinal j = 0; j < domain_->ny; ++j) {
                        Ordinal cell_idx = 0 + j * domain_->nx + k * domain_->nx * domain_->ny;
                        if (h_active_3d(cell_idx) > 0) {
                            bc.cell_indices.push_back(cell_idx);
                        }
                    }
                }
            } else if (bc.face == "y+") {
                // Top face: j = ny - 1
                for (Ordinal k = 0; k < domain_->nz; ++k) {
                    for (Ordinal i = 0; i < domain_->nx; ++i) {
                        Ordinal cell_idx = i + (domain_->ny - 1) * domain_->nx + k * domain_->nx * domain_->ny;
                        if (h_active_3d(cell_idx) > 0) {
                            bc.cell_indices.push_back(cell_idx);
                        }
                    }
                }
            } else if (bc.face == "y-") {
                // Bottom face: j = 0
                for (Ordinal k = 0; k < domain_->nz; ++k) {
                    for (Ordinal i = 0; i < domain_->nx; ++i) {
                        Ordinal cell_idx = i + 0 * domain_->nx + k * domain_->nx * domain_->ny;
                        if (h_active_3d(cell_idx) > 0) {
                            bc.cell_indices.push_back(cell_idx);
                        }
                    }
                }
            } else if (bc.face == "z+") {
                // Bottom layer: k = nz - 1
                for (Ordinal j = 0; j < domain_->ny; ++j) {
                    for (Ordinal i = 0; i < domain_->nx; ++i) {
                        Ordinal cell_idx = i + j * domain_->nx + (domain_->nz - 1) * domain_->nx * domain_->ny;
                        if (h_active_3d(cell_idx) > 0) {
                            bc.cell_indices.push_back(cell_idx);
                        }
                    }
                }
            } else if (bc.face == "z-") {
                // Top layer: k = 0
                for (Ordinal j = 0; j < domain_->ny; ++j) {
                    for (Ordinal i = 0; i < domain_->nx; ++i) {
                        Ordinal cell_idx = i + j * domain_->nx + 0 * domain_->nx * domain_->ny;
                        if (h_active_3d(cell_idx) > 0) {
                            bc.cell_indices.push_back(cell_idx);
                        }
                    }
                }
            }
        } else {
            // Cell-by-cell: iterate through range
            Ordinal k_start = (bc.k_min >= 0) ? bc.k_min : 0;
            Ordinal k_end = (bc.k_max >= 0) ? bc.k_max : (domain_->nz - 1);
            
            for (Ordinal k = k_start; k <= k_end && k < domain_->nz; ++k) {
                for (Ordinal j = bc.j_min; j <= bc.j_max; ++j) {
                    for (Ordinal i = bc.i_min; i <= bc.i_max; ++i) {
                        Ordinal cell_idx = i + j * domain_->nx + k * domain_->nx * domain_->ny;
                        
                        // Check if cell is active
                        if (h_active_3d(cell_idx) <= 0) continue;
                        
                        // Check if cell is at domain boundary or has inactive neighbor
                        bool is_boundary = false;
                        
                        // Check if at domain edge
                        if (i == 0 || i == domain_->nx - 1 || 
                            j == 0 || j == domain_->ny - 1 ||
                            k == 0 || k == domain_->nz - 1) {
                            is_boundary = true;
                        } else {
                            // Check neighbors
                            Ordinal left = (i - 1) + j * domain_->nx + k * domain_->nx * domain_->ny;
                            Ordinal right = (i + 1) + j * domain_->nx + k * domain_->nx * domain_->ny;
                            Ordinal back = i + (j - 1) * domain_->nx + k * domain_->nx * domain_->ny;
                            Ordinal front = i + (j + 1) * domain_->nx + k * domain_->nx * domain_->ny;
                            Ordinal bottom = i + j * domain_->nx + (k - 1) * domain_->nx * domain_->ny;
                            Ordinal top = i + j * domain_->nx + (k + 1) * domain_->nx * domain_->ny;
                            
                            if (h_active_3d(left) <= 0 || h_active_3d(right) <= 0 ||
                                h_active_3d(back) <= 0 || h_active_3d(front) <= 0 ||
                                h_active_3d(bottom) <= 0 || h_active_3d(top) <= 0) {
                                is_boundary = true;
                            }
                        }
                        
                        if (is_boundary) {
                            bc.cell_indices.push_back(cell_idx);
                        }
                    }
                }
            }
        }
    }
};

// ============================================================================
//                      SCALAR TRANSPORT BOUNDARY CONDITION TYPES
// ============================================================================
enum class ScalarBcType {
    PRESCRIBED_CONCENTRATION,  // Prescribed scalar concentration (Dirichlet)
    CAUCHY                     // Cauchy type (mixed) for interface exchange
};

// ============================================================================
//                      SCALAR TRANSPORT BOUNDARY CONDITION
// ============================================================================
class ScalarBoundaryCondition {
public:
    ScalarBcType type;                    // Type of boundary condition
    BcValueType value_type;                // Constant or time series
    
    // Scalar index (which scalar species this BC applies to)
    Ordinal scalar_index;
    
    // Location range (in grid coordinates)
    Ordinal i_min, i_max;                 // X-direction range
    Ordinal j_min, j_max;                 // Y-direction range
    Ordinal k_min, k_max;                 // Z-direction range (for 3D, -1 means all layers)
    
    // Boundary face (if specified, applies to entire face)
    std::string face;                      // "x+", "x-", "y+", "y-", "z+", "z-" or "cell"
    
    // For Cauchy type: exchange coefficient
    Scalar exchange_coefficient;           // α in: flux = α * (s_external - s_internal)
    
    // Cell indices where BC is applied (computed from range)
    std::vector<Ordinal> cell_indices;
    
    // Constant value (if value_type == CONSTANT)
    Scalar constant_value;
    
    // Time series data (if value_type == TIME_SERIES)
    std::vector<Scalar> time_values;       // Time points
    std::vector<Scalar> data_values;       // BC values at each time point
    std::string time_series_file;          // File path for time series
    
    // Constructor
    ScalarBoundaryCondition() 
        : type(ScalarBcType::PRESCRIBED_CONCENTRATION),
          value_type(BcValueType::CONSTANT),
          scalar_index(0),
          i_min(0), i_max(0), j_min(0), j_max(0), k_min(-1), k_max(-1),
          face("cell"), exchange_coefficient(1.0), constant_value(0.0) {}
    
    // Get current BC value at given time
    Scalar get_value(Scalar current_time) const {
        if (value_type == BcValueType::CONSTANT) {
            return constant_value;
        } else {
            // Linear interpolation for time series
            if (time_values.empty()) {
                return constant_value;
            }
            
            if (current_time <= time_values.front()) {
                return data_values.front();
            }
            if (current_time >= time_values.back()) {
                return data_values.back();
            }
            
            for (size_t i = 0; i < time_values.size() - 1; ++i) {
                if (time_values[i] <= current_time && current_time <= time_values[i+1]) {
                    Scalar t0 = time_values[i];
                    Scalar t1 = time_values[i+1];
                    Scalar v0 = data_values[i];
                    Scalar v1 = data_values[i+1];
                    Scalar alpha = (current_time - t0) / (t1 - t0);
                    return v0 + alpha * (v1 - v0);
                }
            }
            
            return data_values.back();
        }
    }
    
    // Load time series from file
    void load_time_series(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open scalar boundary condition time series file: " + filename);
        }
        
        time_values.clear();
        data_values.clear();
        
        std::string line;
        while (std::getline(file, line)) {
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            if (line.empty() || line[0] == '#') continue;
            
            std::istringstream iss(line);
            Scalar t, v;
            if (iss >> t >> v) {
                time_values.push_back(t);
                data_values.push_back(v);
            }
        }
        
        file.close();
        
        if (time_values.empty()) {
            throw std::runtime_error("No valid data found in scalar boundary condition file: " + filename);
        }
        
        // Sort by time
        std::vector<size_t> indices(time_values.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), 
                 [&](size_t a, size_t b) { return time_values[a] < time_values[b]; });
        
        std::vector<Scalar> sorted_times, sorted_values;
        for (size_t idx : indices) {
            sorted_times.push_back(time_values[idx]);
            sorted_values.push_back(data_values[idx]);
        }
        time_values = sorted_times;
        data_values = sorted_values;
    }
};

// ============================================================================
//                      SURFACE WATER SCALAR BOUNDARY CONDITION MANAGER
// ============================================================================
class SwScalarBoundaryConditionManager {
private:
    SwDomain* domain_;
    Ordinal num_scalars_;
    std::vector<std::vector<ScalarBoundaryCondition>> boundary_conditions_;  // [scalar_index][bc_index]
    
public:
    SwScalarBoundaryConditionManager(SwDomain* domain, Ordinal num_scalars) 
        : domain_(domain), num_scalars_(num_scalars) {
        boundary_conditions_.resize(num_scalars_);
    }
    
    // Read boundary conditions from input file
    void read_from_input(InputReader* input_reader, const std::string& input_dir) {
        for (Ordinal kk = 0; kk < num_scalars_; ++kk) {
            // Read number of boundary conditions for this scalar
            Ordinal n_bc = input_reader->read_int("n_sw_scalar_bc_" + std::to_string(kk), 0);
            
            for (Ordinal i = 0; i < n_bc; ++i) {
                ScalarBoundaryCondition bc;
                bc.scalar_index = kk;
                
                // Read BC type
                std::string type_str = input_reader->read_string("sw_scalar_bc_type_" + std::to_string(kk) + "_" + std::to_string(i), "prescribed");
                if (type_str == "prescribed" || type_str == "concentration") {
                    bc.type = ScalarBcType::PRESCRIBED_CONCENTRATION;
                } else if (type_str == "cauchy") {
                    bc.type = ScalarBcType::CAUCHY;
                    bc.exchange_coefficient = input_reader->read_double("sw_scalar_bc_exchange_" + std::to_string(kk) + "_" + std::to_string(i), 1.0);
                } else {
                    throw std::runtime_error("Unknown scalar BC type: " + type_str);
                }
                
                // Read location range
                std::vector<int> loc_range = input_reader->read_int_array("sw_scalar_bc_loc_" + std::to_string(kk) + "_" + std::to_string(i), 
                                                                          std::vector<int>{0, 0, 0, 0});
                if (loc_range.size() != 4) {
                    throw std::runtime_error("sw_scalar_bc_loc must have 4 values: [i_min, i_max, j_min, j_max]");
                }
                bc.i_min = static_cast<Ordinal>(loc_range[0]);
                bc.i_max = static_cast<Ordinal>(loc_range[1]);
                bc.j_min = static_cast<Ordinal>(loc_range[2]);
                bc.j_max = static_cast<Ordinal>(loc_range[3]);
                
                // Read value type
                std::string value_type_str = input_reader->read_string("sw_scalar_bc_value_type_" + std::to_string(kk) + "_" + std::to_string(i), "constant");
                if (value_type_str == "constant") {
                    bc.value_type = BcValueType::CONSTANT;
                    bc.constant_value = input_reader->read_double("sw_scalar_bc_value_" + std::to_string(kk) + "_" + std::to_string(i), 0.0);
                } else if (value_type_str == "time_series" || value_type_str == "file") {
                    bc.value_type = BcValueType::TIME_SERIES;
                    bc.time_series_file = input_reader->read_string("sw_scalar_bc_file_" + std::to_string(kk) + "_" + std::to_string(i), "");
                    if (bc.time_series_file.empty()) {
                        bc.time_series_file = input_dir + "/sw_scalar_" + std::to_string(kk) + "_bc_" + std::to_string(i) + ".txt";
                    }
                    bc.load_time_series(bc.time_series_file);
                }
                
                // Identify boundary cells (similar to SwBoundaryConditionManager)
                identify_boundary_cells(bc);
                
                boundary_conditions_[kk].push_back(bc);
            }
        }
    }
    
    // Get boundary conditions for a specific scalar
    const std::vector<ScalarBoundaryCondition>& get_boundary_conditions(Ordinal scalar_index) const {
        if (scalar_index >= num_scalars_) {
            throw std::runtime_error("Scalar index out of range");
        }
        return boundary_conditions_[scalar_index];
    }
    
    // Get boundary conditions for a specific cell and scalar
    std::vector<const ScalarBoundaryCondition*> get_bc_for_cell(Ordinal scalar_index, Ordinal cell_idx) const {
        std::vector<const ScalarBoundaryCondition*> result;
        if (scalar_index >= num_scalars_) return result;
        
        for (const auto& bc : boundary_conditions_[scalar_index]) {
            if (std::find(bc.cell_indices.begin(), bc.cell_indices.end(), cell_idx) 
                != bc.cell_indices.end()) {
                result.push_back(&bc);
            }
        }
        return result;
    }
    
private:
    // Identify boundary cells (same logic as SwBoundaryConditionManager)
    void identify_boundary_cells(ScalarBoundaryCondition& bc) {
        bc.cell_indices.clear();
        
        auto h_active = Kokkos::create_mirror_view(domain_->active_mask);
        Kokkos::deep_copy(h_active, domain_->active_mask);
        
        for (Ordinal j = bc.j_min; j <= bc.j_max; ++j) {
            for (Ordinal i = bc.i_min; i <= bc.i_max; ++i) {
                Ordinal cell_idx = i + j * domain_->nx;
                
                if (h_active(cell_idx) <= 0) continue;
                
                bool is_boundary = false;
                
                if (i == 0 || i == domain_->nx - 1 || 
                    j == 0 || j == domain_->ny - 1) {
                    is_boundary = true;
                } else {
                    Ordinal left = (i - 1) + j * domain_->nx;
                    Ordinal right = (i + 1) + j * domain_->nx;
                    Ordinal bottom = i + (j - 1) * domain_->nx;
                    Ordinal top = i + (j + 1) * domain_->nx;
                    
                    if (h_active(left) <= 0 || h_active(right) <= 0 ||
                        h_active(bottom) <= 0 || h_active(top) <= 0) {
                        is_boundary = true;
                    }
                }
                
                if (is_boundary) {
                    bc.cell_indices.push_back(cell_idx);
                }
            }
        }
    }
};

// ============================================================================
//                      GROUNDWATER SCALAR BOUNDARY CONDITION MANAGER
// ============================================================================
class GwScalarBoundaryConditionManager {
private:
    GwDomain* domain_;
    Ordinal num_scalars_;
    std::vector<std::vector<ScalarBoundaryCondition>> boundary_conditions_;  // [scalar_index][bc_index]
    
public:
    GwScalarBoundaryConditionManager(GwDomain* domain, Ordinal num_scalars) 
        : domain_(domain), num_scalars_(num_scalars) {
        boundary_conditions_.resize(num_scalars_);
    }
    
    // Read boundary conditions from input file
    void read_from_input(InputReader* input_reader, const std::string& input_dir) {
        for (Ordinal kk = 0; kk < num_scalars_; ++kk) {
            Ordinal n_bc = input_reader->read_int("n_gw_scalar_bc_" + std::to_string(kk), 0);
            
            for (Ordinal i = 0; i < n_bc; ++i) {
                ScalarBoundaryCondition bc;
                bc.scalar_index = kk;
                
                std::string type_str = input_reader->read_string("gw_scalar_bc_type_" + std::to_string(kk) + "_" + std::to_string(i), "prescribed");
                if (type_str == "prescribed" || type_str == "concentration") {
                    bc.type = ScalarBcType::PRESCRIBED_CONCENTRATION;
                } else if (type_str == "cauchy") {
                    bc.type = ScalarBcType::CAUCHY;
                    bc.exchange_coefficient = input_reader->read_double("gw_scalar_bc_exchange_" + std::to_string(kk) + "_" + std::to_string(i), 1.0);
                } else {
                    throw std::runtime_error("Unknown scalar BC type: " + type_str);
                }
                
                bc.face = input_reader->read_string("gw_scalar_bc_face_" + std::to_string(kk) + "_" + std::to_string(i), "cell");
                
                if (bc.face == "cell") {
                    std::vector<int> loc_range = input_reader->read_int_array("gw_scalar_bc_loc_" + std::to_string(kk) + "_" + std::to_string(i), 
                                                                              std::vector<int>{0, 0, 0, 0, -1, -1});
                    if (loc_range.size() < 4) {
                        throw std::runtime_error("gw_scalar_bc_loc must have at least 4 values");
                    }
                    bc.i_min = static_cast<Ordinal>(loc_range[0]);
                    bc.i_max = static_cast<Ordinal>(loc_range[1]);
                    bc.j_min = static_cast<Ordinal>(loc_range[2]);
                    bc.j_max = static_cast<Ordinal>(loc_range[3]);
                    if (loc_range.size() >= 6) {
                        bc.k_min = static_cast<Ordinal>(loc_range[4]);
                        bc.k_max = static_cast<Ordinal>(loc_range[5]);
                    }
                }
                
                std::string value_type_str = input_reader->read_string("gw_scalar_bc_value_type_" + std::to_string(kk) + "_" + std::to_string(i), "constant");
                if (value_type_str == "constant") {
                    bc.value_type = BcValueType::CONSTANT;
                    bc.constant_value = input_reader->read_double("gw_scalar_bc_value_" + std::to_string(kk) + "_" + std::to_string(i), 0.0);
                } else if (value_type_str == "time_series" || value_type_str == "file") {
                    bc.value_type = BcValueType::TIME_SERIES;
                    bc.time_series_file = input_reader->read_string("gw_scalar_bc_file_" + std::to_string(kk) + "_" + std::to_string(i), "");
                    if (bc.time_series_file.empty()) {
                        bc.time_series_file = input_dir + "/gw_scalar_" + std::to_string(kk) + "_bc_" + std::to_string(i) + ".txt";
                    }
                    bc.load_time_series(bc.time_series_file);
                }
                
                identify_boundary_cells(bc);
                boundary_conditions_[kk].push_back(bc);
            }
        }
    }
    
    const std::vector<ScalarBoundaryCondition>& get_boundary_conditions(Ordinal scalar_index) const {
        if (scalar_index >= num_scalars_) {
            throw std::runtime_error("Scalar index out of range");
        }
        return boundary_conditions_[scalar_index];
    }
    
    std::vector<const ScalarBoundaryCondition*> get_bc_for_cell(Ordinal scalar_index, Ordinal cell_idx) const {
        std::vector<const ScalarBoundaryCondition*> result;
        if (scalar_index >= num_scalars_) return result;
        
        for (const auto& bc : boundary_conditions_[scalar_index]) {
            if (std::find(bc.cell_indices.begin(), bc.cell_indices.end(), cell_idx) 
                != bc.cell_indices.end()) {
                result.push_back(&bc);
            }
        }
        return result;
    }
    
private:
    // Identify boundary cells (similar to GwBoundaryConditionManager)
    void identify_boundary_cells(ScalarBoundaryCondition& bc) {
        bc.cell_indices.clear();
        
        auto h_active_3d = Kokkos::create_mirror_view(domain_->active_mask_3d);
        Kokkos::deep_copy(h_active_3d, domain_->active_mask_3d);
        
        if (bc.face != "cell") {
            // Face boundary (same logic as GwBoundaryConditionManager)
            if (bc.face == "x+") {
                for (Ordinal k = 0; k < domain_->nz; ++k) {
                    for (Ordinal j = 0; j < domain_->ny; ++j) {
                        Ordinal cell_idx = (domain_->nx - 1) + j * domain_->nx + k * domain_->nx * domain_->ny;
                        if (h_active_3d(cell_idx) > 0) {
                            bc.cell_indices.push_back(cell_idx);
                        }
                    }
                }
            } else if (bc.face == "x-") {
                for (Ordinal k = 0; k < domain_->nz; ++k) {
                    for (Ordinal j = 0; j < domain_->ny; ++j) {
                        Ordinal cell_idx = 0 + j * domain_->nx + k * domain_->nx * domain_->ny;
                        if (h_active_3d(cell_idx) > 0) {
                            bc.cell_indices.push_back(cell_idx);
                        }
                    }
                }
            } else if (bc.face == "y+") {
                for (Ordinal k = 0; k < domain_->nz; ++k) {
                    for (Ordinal i = 0; i < domain_->nx; ++i) {
                        Ordinal cell_idx = i + (domain_->ny - 1) * domain_->nx + k * domain_->nx * domain_->ny;
                        if (h_active_3d(cell_idx) > 0) {
                            bc.cell_indices.push_back(cell_idx);
                        }
                    }
                }
            } else if (bc.face == "y-") {
                for (Ordinal k = 0; k < domain_->nz; ++k) {
                    for (Ordinal i = 0; i < domain_->nx; ++i) {
                        Ordinal cell_idx = i + 0 * domain_->nx + k * domain_->nx * domain_->ny;
                        if (h_active_3d(cell_idx) > 0) {
                            bc.cell_indices.push_back(cell_idx);
                        }
                    }
                }
            } else if (bc.face == "z+") {
                for (Ordinal j = 0; j < domain_->ny; ++j) {
                    for (Ordinal i = 0; i < domain_->nx; ++i) {
                        Ordinal cell_idx = i + j * domain_->nx + (domain_->nz - 1) * domain_->nx * domain_->ny;
                        if (h_active_3d(cell_idx) > 0) {
                            bc.cell_indices.push_back(cell_idx);
                        }
                    }
                }
            } else if (bc.face == "z-") {
                for (Ordinal j = 0; j < domain_->ny; ++j) {
                    for (Ordinal i = 0; i < domain_->nx; ++i) {
                        Ordinal cell_idx = i + j * domain_->nx + 0 * domain_->nx * domain_->ny;
                        if (h_active_3d(cell_idx) > 0) {
                            bc.cell_indices.push_back(cell_idx);
                        }
                    }
                }
            }
        } else {
            // Cell-by-cell (same logic as GwBoundaryConditionManager)
            Ordinal k_start = (bc.k_min >= 0) ? bc.k_min : 0;
            Ordinal k_end = (bc.k_max >= 0) ? bc.k_max : (domain_->nz - 1);
            
            for (Ordinal k = k_start; k <= k_end && k < domain_->nz; ++k) {
                for (Ordinal j = bc.j_min; j <= bc.j_max; ++j) {
                    for (Ordinal i = bc.i_min; i <= bc.i_max; ++i) {
                        Ordinal cell_idx = i + j * domain_->nx + k * domain_->nx * domain_->ny;
                        
                        if (h_active_3d(cell_idx) <= 0) continue;
                        
                        bool is_boundary = false;
                        
                        if (i == 0 || i == domain_->nx - 1 || 
                            j == 0 || j == domain_->ny - 1 ||
                            k == 0 || k == domain_->nz - 1) {
                            is_boundary = true;
                        } else {
                            Ordinal left = (i - 1) + j * domain_->nx + k * domain_->nx * domain_->ny;
                            Ordinal right = (i + 1) + j * domain_->nx + k * domain_->nx * domain_->ny;
                            Ordinal back = i + (j - 1) * domain_->nx + k * domain_->nx * domain_->ny;
                            Ordinal front = i + (j + 1) * domain_->nx + k * domain_->nx * domain_->ny;
                            Ordinal bottom = i + j * domain_->nx + (k - 1) * domain_->nx * domain_->ny;
                            Ordinal top = i + j * domain_->nx + (k + 1) * domain_->nx * domain_->ny;
                            
                            if (h_active_3d(left) <= 0 || h_active_3d(right) <= 0 ||
                                h_active_3d(back) <= 0 || h_active_3d(front) <= 0 ||
                                h_active_3d(bottom) <= 0 || h_active_3d(top) <= 0) {
                                is_boundary = true;
                            }
                        }
                        
                        if (is_boundary) {
                            bc.cell_indices.push_back(cell_idx);
                        }
                    }
                }
            }
        }
    }
};

} // namespace Frehg

#endif // FREHG_BOUNDARY_CONDITIONS_HPP

