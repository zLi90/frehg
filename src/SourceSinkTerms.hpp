#ifndef FREHG_SOURCE_SINK_TERMS_HPP
#define FREHG_SOURCE_SINK_TERMS_HPP

#include "define.hpp"
#include "Domain.hpp"
#include "InputReader.hpp"
#include "BoundaryConditions.hpp"  // For BcValueType enum
#include <vector>
#include <string>
#include <memory>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <numeric>

namespace Frehg {

// ============================================================================
//                      SOURCE/SINK TERM TYPES
// ============================================================================
enum class SourceSinkType {
    VOLUME_FLUX,           // Volume flux (m³/s) - adds/removes water
    DEPTH_RATE,            // Depth rate (m/s) - adds/removes water depth
    MASS_FLUX              // Mass flux (kg/s) - for scalar transport
};

// ============================================================================
//                      SURFACE WATER SOURCE/SINK TERM
// ============================================================================
class SwSourceSinkTerm {
public:
    SourceSinkType type;                   // Type of source/sink term
    BcValueType value_type;                // Constant or time series
    
    // Location range (in grid coordinates: i, j indices)
    Ordinal i_min, i_max;                  // X-direction range
    Ordinal j_min, j_max;                  // Y-direction range
    
    // Cell indices where source/sink is applied (computed from range)
    std::vector<Ordinal> cell_indices;
    
    // Constant value (if value_type == CONSTANT)
    Scalar constant_value;
    
    // Time series data (if value_type == TIME_SERIES)
    std::vector<Scalar> time_values;       // Time points
    std::vector<Scalar> data_values;       // Source/sink values at each time point
    std::string time_series_file;          // File path for time series
    
    // For scalar transport: concentration of source (if applicable)
    Scalar scalar_concentration;           // Concentration in source water
    bool has_scalar_concentration;         // Whether scalar concentration is specified
    
    // Constructor
    SwSourceSinkTerm() 
        : type(SourceSinkType::VOLUME_FLUX),
          value_type(BcValueType::CONSTANT),
          i_min(0), i_max(0), j_min(0), j_max(0),
          constant_value(0.0),
          scalar_concentration(0.0),
          has_scalar_concentration(false) {}
    
    // Get current value at given time
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
            throw std::runtime_error("Cannot open source/sink time series file: " + filename);
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
            throw std::runtime_error("No valid data found in source/sink file: " + filename);
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
//                      SURFACE WATER SOURCE/SINK MANAGER
// ============================================================================
class SwSourceSinkManager {
private:
    SwDomain* domain_;
    std::vector<SwSourceSinkTerm> source_sink_terms_;
    
public:
    SwSourceSinkManager(SwDomain* domain) : domain_(domain) {}
    
    // Read source/sink terms from input file
    void read_from_input(InputReader* input_reader, const std::string& input_dir) {
        // Read number of source/sink terms
        Ordinal n_ss = input_reader->read_int("n_sw_source_sink", 0);
        
        for (Ordinal i = 0; i < n_ss; ++i) {
            SwSourceSinkTerm ss;
            
            // Read source/sink type
            std::string type_str = input_reader->read_string("sw_ss_type_" + std::to_string(i), "volume_flux");
            if (type_str == "volume_flux" || type_str == "flux") {
                ss.type = SourceSinkType::VOLUME_FLUX;
            } else if (type_str == "depth_rate" || type_str == "depth") {
                ss.type = SourceSinkType::DEPTH_RATE;
            } else if (type_str == "mass_flux") {
                ss.type = SourceSinkType::MASS_FLUX;
            } else {
                throw std::runtime_error("Unknown source/sink type: " + type_str);
            }
            
            // Read location range
            std::vector<int> loc_range = input_reader->read_int_array("sw_ss_loc_" + std::to_string(i), 
                                                                      std::vector<int>{0, 0, 0, 0});
            if (loc_range.size() != 4) {
                throw std::runtime_error("sw_ss_loc must have 4 values: [i_min, i_max, j_min, j_max]");
            }
            ss.i_min = static_cast<Ordinal>(loc_range[0]);
            ss.i_max = static_cast<Ordinal>(loc_range[1]);
            ss.j_min = static_cast<Ordinal>(loc_range[2]);
            ss.j_max = static_cast<Ordinal>(loc_range[3]);
            
            // Validate range
            if (ss.i_min > ss.i_max || ss.j_min > ss.j_max) {
                throw std::runtime_error("Invalid source/sink location range");
            }
            if (ss.i_max >= domain_->nx || ss.j_max >= domain_->ny) {
                throw std::runtime_error("Source/sink location range exceeds domain dimensions");
            }
            
            // Read value type
            std::string value_type_str = input_reader->read_string("sw_ss_value_type_" + std::to_string(i), "constant");
            if (value_type_str == "constant") {
                ss.value_type = BcValueType::CONSTANT;
                ss.constant_value = input_reader->read_double("sw_ss_value_" + std::to_string(i), 0.0);
            } else if (value_type_str == "time_series" || value_type_str == "file") {
                ss.value_type = BcValueType::TIME_SERIES;
                ss.time_series_file = input_reader->read_string("sw_ss_file_" + std::to_string(i), "");
                if (ss.time_series_file.empty()) {
                    ss.time_series_file = input_dir + "/sw_ss_" + std::to_string(i) + ".txt";
                }
                ss.load_time_series(ss.time_series_file);
            } else {
                throw std::runtime_error("Unknown source/sink value type: " + value_type_str);
            }
            
            // Read scalar concentration (optional, for scalar transport)
            ss.has_scalar_concentration = input_reader->read_bool("sw_ss_has_scalar_" + std::to_string(i), false);
            if (ss.has_scalar_concentration) {
                ss.scalar_concentration = input_reader->read_double("sw_ss_scalar_" + std::to_string(i), 0.0);
            }
            
            // Identify cells within range
            identify_cells(ss);
            
            source_sink_terms_.push_back(ss);
        }
    }
    
    // Get all source/sink terms
    const std::vector<SwSourceSinkTerm>& get_source_sink_terms() const {
        return source_sink_terms_;
    }
    
    // Get source/sink terms for a specific cell
    std::vector<const SwSourceSinkTerm*> get_ss_for_cell(Ordinal cell_idx) const {
        std::vector<const SwSourceSinkTerm*> result;
        for (const auto& ss : source_sink_terms_) {
            if (std::find(ss.cell_indices.begin(), ss.cell_indices.end(), cell_idx) 
                != ss.cell_indices.end()) {
                result.push_back(&ss);
            }
        }
        return result;
    }
    
private:
    // Identify cells within range (all active cells in range, not just boundaries)
    void identify_cells(SwSourceSinkTerm& ss) {
        ss.cell_indices.clear();
        
        auto h_active = Kokkos::create_mirror_view(domain_->active_mask);
        Kokkos::deep_copy(h_active, domain_->active_mask);
        
        for (Ordinal j = ss.j_min; j <= ss.j_max; ++j) {
            for (Ordinal i = ss.i_min; i <= ss.i_max; ++i) {
                Ordinal cell_idx = i + j * domain_->nx;
                
                // Include all active cells in range
                if (h_active(cell_idx) > 0) {
                    ss.cell_indices.push_back(cell_idx);
                }
            }
        }
    }
};

} // namespace Frehg

#endif // FREHG_SOURCE_SINK_TERMS_HPP

