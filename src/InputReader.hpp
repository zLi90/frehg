#ifndef FREHG_INPUT_READER_HPP
#define FREHG_INPUT_READER_HPP

#include "define.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>

namespace Frehg {

// ============================================================================
//                      STRUCTURED GRID DATA STRUCTURE
// ============================================================================
// Structure to hold structured grid file information and data
// Supports ESRI ASCII Grid format
struct StructuredGrid {
    Ordinal ncols;
    Ordinal nrows;
    Scalar xllcorner;      // or xllcenter (lower-left x coordinate)
    Scalar yllcorner;      // or yllcenter (lower-left y coordinate)
    Scalar cellsize;
    Scalar nodata_value;
    bool use_corner;       // true if xllcorner/yllcorner, false if xllcenter/yllcenter
    
    std::vector<Scalar> data;  // Grid data (stored as i + j * ncols, where j=0 is bottom row)
    std::vector<int> active_mask;  // Active mask (1=active, 0=inactive based on nodata)
    
    // Default constructor
    StructuredGrid() : ncols(0), nrows(0), xllcorner(0.0), yllcorner(0.0),
                      cellsize(0.0), nodata_value(-9999.0), use_corner(true) {}
};

// ============================================================================
//                      READ STRUCTURED DATA FILE
// ============================================================================
// Reads a structured grid data file (ESRI ASCII Grid format)
// Format:
//   ncols <value>
//   nrows <value>
//   xllcorner <value>  (or xllcenter)
//   yllcorner <value>  (or yllcenter)
//   cellsize <value>
//   nodata_value <value>
//   <data values> (space or newline separated, stored top-to-bottom, left-to-right)
//
// Returns: StructuredGrid structure with data converted to bottom-to-top indexing
//          (our internal format: index = i + j * ncols, where j=0 is bottom row)
//          Active mask is automatically generated based on nodata_value
inline StructuredGrid read_structured_data(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open structured data file: " + filename);
    }
    
    StructuredGrid grid;
    std::string line;
    
    // Helper function to trim whitespace
    auto trim = [](const std::string& str) -> std::string {
        size_t first = str.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) return "";
        size_t last = str.find_last_not_of(" \t\r\n");
        return str.substr(first, (last - first + 1));
    };
    
    // Read headers
    bool found_data_start = false;
    while (std::getline(file, line) && !found_data_start) {
        line = trim(line);
        if (line.empty()) continue;
        
        // Convert to lowercase for case-insensitive matching
        std::string lower = line;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        
        // Parse header fields
        if (lower.find("ncols") == 0) {
            std::istringstream iss(line);
            std::string key;
            iss >> key >> grid.ncols;
        } else if (lower.find("nrows") == 0) {
            std::istringstream iss(line);
            std::string key;
            iss >> key >> grid.nrows;
        } else if (lower.find("xllcorner") == 0) {
            std::istringstream iss(line);
            std::string key;
            iss >> key >> grid.xllcorner;
            grid.use_corner = true;
        } else if (lower.find("xllcenter") == 0) {
            std::istringstream iss(line);
            std::string key;
            iss >> key >> grid.xllcorner;
            grid.use_corner = false;
        } else if (lower.find("yllcorner") == 0) {
            std::istringstream iss(line);
            std::string key;
            iss >> key >> grid.yllcorner;
            grid.use_corner = true;
        } else if (lower.find("yllcenter") == 0) {
            std::istringstream iss(line);
            std::string key;
            iss >> key >> grid.yllcorner;
            grid.use_corner = false;
        } else if (lower.find("cellsize") == 0) {
            std::istringstream iss(line);
            std::string key;
            iss >> key >> grid.cellsize;
        } else if (lower.find("nodata_value") == 0 || lower.find("nodata") == 0) {
            std::istringstream iss(line);
            std::string key;
            iss >> key >> grid.nodata_value;
        } else {
            // First non-header line - this is the start of data
            // Push this line back so we can read it as data
            file.seekg(-static_cast<long>(line.length() + 1), std::ios::cur);
            found_data_start = true;
        }
    }
    
    // Validate headers
    if (grid.ncols == 0 || grid.nrows == 0) {
        throw std::runtime_error("Invalid structured data file: ncols or nrows not specified in " + filename);
    }
    
    // Allocate data arrays
    grid.data.resize(grid.ncols * grid.nrows);
    grid.active_mask.resize(grid.ncols * grid.nrows, 0);
    
    // Read elevation data
    // ESRI ASCII Grid stores data from top-left, row by row (top to bottom)
    // We need to convert to our format: bottom-left, row by row (bottom to top)
    // Our format: index = i + j * ncols, where i is column, j is row (j=0 is bottom)
    
    for (Ordinal j = 0; j < grid.nrows; ++j) {
        // Read one row from file (this is row from top: j_file = grid.nrows - 1 - j)
        Ordinal row_idx = grid.nrows - 1 - j;  // Convert to bottom-up indexing
        
        for (Ordinal i = 0; i < grid.ncols; ++i) {
            Scalar value;
            if (!(file >> value)) {
                throw std::runtime_error("Unexpected end of file while reading grid data in " + filename + 
                                       " at row " + std::to_string(j) + ", col " + std::to_string(i));
            }
            
            // Calculate index in our format (i, j from bottom)
            Ordinal idx = i + row_idx * grid.ncols;
            grid.data[idx] = value;
            
            // Determine if cell is active (not nodata)
            if (std::abs(value - grid.nodata_value) > 1e-10) {
                grid.active_mask[idx] = 1;
            } else {
                grid.active_mask[idx] = 0;
            }
        }
    }
    
    file.close();
    return grid;
}

// ============================================================================
//                      INPUT FILE READER
// ============================================================================
// Flexible input file reader that supports:
// - Comments (lines starting with # or //)
// - Empty lines
// - Format: field_name = value
// - String values
// - Numeric values (int, double)
// - Arrays (comma-separated values)
// - Easy to add/remove fields without code changes

class InputReader {
private:
    std::string filename_;
    std::vector<std::string> lines_;  // Cached lines for multiple reads
    
    // Read and cache all lines from file
    void read_file() {
        std::ifstream file(filename_);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open input file: " + filename_);
        }
        
        lines_.clear();
        std::string line;
        while (std::getline(file, line)) {
            lines_.push_back(line);
        }
        file.close();
    }
    
    // Remove leading and trailing whitespace
    static std::string trim(const std::string& str) {
        size_t first = str.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) {
            return "";
        }
        size_t last = str.find_last_not_of(" \t\r\n");
        return str.substr(first, (last - first + 1));
    }
    
    // Check if line is a comment or empty
    static bool is_comment_or_empty(const std::string& line) {
        std::string trimmed = trim(line);
        return trimmed.empty() || 
               trimmed[0] == '#' || 
               (trimmed.length() >= 2 && trimmed.substr(0, 2) == "//");
    }
    
    // Find the line containing the field
    std::string find_field_line(const std::string& field_name) const {
        for (const auto& line : lines_) {
            // Skip comments and empty lines
            if (is_comment_or_empty(line)) {
                continue;
            }
            
            // Check if this line contains the field
            size_t eq_pos = line.find('=');
            if (eq_pos == std::string::npos) {
                continue;  // No '=' found, skip this line
            }
            
            std::string field_part = trim(line.substr(0, eq_pos));
            if (field_part == field_name) {
                return line;
            }
        }
        return "";  // Field not found
    }

public:
    // Constructor
    InputReader(const std::string& filename) : filename_(filename) {
        read_file();
    }
    
    // Reload file (useful if file is modified)
    void reload() {
        read_file();
    }
    
    // ========================================================================
    // READ ONE FIELD (Generic - returns string value)
    // ========================================================================
    // This is the base function that all other read functions use
    // Returns the value as a string, or throws if field not found
    std::string read_field(const std::string& field_name, 
                          const std::string& default_value = "") const {
        std::string line = find_field_line(field_name);
        
        if (line.empty()) {
            if (!default_value.empty()) {
                return default_value;
            }
            throw std::runtime_error("Field '" + field_name + "' not found in input file: " + filename_);
        }
        
        // Extract value part (after '=')
        size_t eq_pos = line.find('=');
        if (eq_pos == std::string::npos) {
            throw std::runtime_error("Invalid format for field '" + field_name + "' in input file: " + filename_);
        }
        
        std::string value = trim(line.substr(eq_pos + 1));
        
        // Remove trailing comment if present
        size_t comment_pos = value.find('#');
        if (comment_pos != std::string::npos) {
            value = trim(value.substr(0, comment_pos));
        }
        comment_pos = value.find("//");
        if (comment_pos != std::string::npos) {
            value = trim(value.substr(0, comment_pos));
        }
        
        return value;
    }
    
    // ========================================================================
    // READ ONE FIELD (String)
    // ========================================================================
    std::string read_string(const std::string& field_name, 
                           const std::string& default_value = "") const {
        return read_field(field_name, default_value);
    }
    
    // ========================================================================
    // READ ONE FIELD (Double)
    // ========================================================================
    double read_double(const std::string& field_name, 
                      double default_value = 0.0) const {
        try {
            std::string value_str = read_field(field_name);
            if (value_str.empty()) {
                return default_value;
            }
            return std::stod(value_str);
        } catch (const std::exception& e) {
            throw std::runtime_error("Cannot convert field '" + field_name + 
                                   "' to double: " + e.what());
        }
    }
    
    // ========================================================================
    // READ ONE FIELD (Integer)
    // ========================================================================
    int read_int(const std::string& field_name, 
                int default_value = 0) const {
        try {
            std::string value_str = read_field(field_name);
            if (value_str.empty()) {
                return default_value;
            }
            return std::stoi(value_str);
        } catch (const std::exception& e) {
            throw std::runtime_error("Cannot convert field '" + field_name + 
                                   "' to int: " + e.what());
        }
    }
    
    // ========================================================================
    // READ ONE FIELD (Boolean)
    // ========================================================================
    bool read_bool(const std::string& field_name, 
                  bool default_value = false) const {
        try {
            std::string value_str = read_field(field_name);
            if (value_str.empty()) {
                return default_value;
            }
            
            // Convert to lowercase for comparison
            std::string lower = value_str;
            std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
            
            if (lower == "1" || lower == "true" || lower == "yes" || lower == "on") {
                return true;
            } else if (lower == "0" || lower == "false" || lower == "no" || lower == "off") {
                return false;
            } else {
                // Try to parse as integer
                int val = std::stoi(value_str);
                return (val != 0);
            }
        } catch (const std::exception& e) {
            throw std::runtime_error("Cannot convert field '" + field_name + 
                                   "' to bool: " + e.what());
        }
    }
    
    // ========================================================================
    // READ ONE FIELD (Array of doubles - comma-separated)
    // ========================================================================
    std::vector<double> read_double_array(const std::string& field_name) const {
        std::string value_str = read_field(field_name);
        
        std::vector<double> result;
        if (value_str.empty()) {
            return result;
        }
        
        // Split by comma
        std::istringstream iss(value_str);
        std::string token;
        while (std::getline(iss, token, ',')) {
            token = trim(token);
            if (!token.empty()) {
                try {
                    result.push_back(std::stod(token));
                } catch (const std::exception& e) {
                    throw std::runtime_error("Cannot convert array element '" + token + 
                                           "' to double in field '" + field_name + "'");
                }
            }
        }
        
        return result;
    }
    
    // ========================================================================
    // READ ONE FIELD (Array of integers - comma-separated)
    // ========================================================================
    std::vector<int> read_int_array(const std::string& field_name) const {
        std::string value_str = read_field(field_name);
        
        std::vector<int> result;
        if (value_str.empty()) {
            return result;
        }
        
        // Split by comma
        std::istringstream iss(value_str);
        std::string token;
        while (std::getline(iss, token, ',')) {
            token = trim(token);
            if (!token.empty()) {
                try {
                    result.push_back(std::stoi(token));
                } catch (const std::exception& e) {
                    throw std::runtime_error("Cannot convert array element '" + token + 
                                           "' to int in field '" + field_name + "'");
                }
            }
        }
        
        return result;
    }
    
    // ========================================================================
    // CHECK IF FIELD EXISTS
    // ========================================================================
    bool has_field(const std::string& field_name) const {
        return !find_field_line(field_name).empty();
    }
    
    // ========================================================================
    // GET FILENAME
    // ========================================================================
    std::string get_filename() const {
        return filename_;
    }
    
    // ========================================================================
    // READ 3D SOIL TYPE MAP
    // ========================================================================
    // Reads a 3D soil type index map from file
    // Format: Binary file (int32) or ASCII file (space/newline separated)
    // Returns: 3D vector (i + j*nx + k*nx*ny indexing)
    static std::vector<int> read_3d_soil_type_map(const std::string& filename,
                                                  Ordinal nx, Ordinal ny, Ordinal nz,
                                                  bool binary = false) {
        std::vector<int> soil_type_map(nx * ny * nz);
        
        if (binary) {
            // Read binary file (int32 per value)
            std::ifstream file(filename, std::ios::binary);
            if (!file.is_open()) {
                throw std::runtime_error("Cannot open binary soil type map file: " + filename);
            }
            
            for (Ordinal k = 0; k < nz; ++k) {
                for (Ordinal j = 0; j < ny; ++j) {
                    for (Ordinal i = 0; i < nx; ++i) {
                        Ordinal idx = i + j * nx + k * nx * ny;
                        int32_t value;
                        file.read(reinterpret_cast<char*>(&value), sizeof(int32_t));
                        if (!file) {
                            throw std::runtime_error("Unexpected end of file in soil type map at "
                                                   "i=" + std::to_string(i) + ", j=" + std::to_string(j) +
                                                   ", k=" + std::to_string(k));
                        }
                        soil_type_map[idx] = static_cast<int>(value);
                    }
                }
            }
            file.close();
        } else {
            // Read ASCII file (space or newline separated)
            std::ifstream file(filename);
            if (!file.is_open()) {
                throw std::runtime_error("Cannot open ASCII soil type map file: " + filename);
            }
            
            for (Ordinal k = 0; k < nz; ++k) {
                for (Ordinal j = 0; j < ny; ++j) {
                    for (Ordinal i = 0; i < nx; ++i) {
                        Ordinal idx = i + j * nx + k * nx * ny;
                        int value;
                        if (!(file >> value)) {
                            throw std::runtime_error("Unexpected end of file in soil type map at "
                                                   "i=" + std::to_string(i) + ", j=" + std::to_string(j) +
                                                   ", k=" + std::to_string(k));
                        }
                        soil_type_map[idx] = value;
                    }
                }
            }
            file.close();
        }
        
        return soil_type_map;
    }
    
    // ========================================================================
    // READ SOIL PARAMETER FILE
    // ========================================================================
    // Reads soil parameters for each soil type from a file
    // Format: CSV or space-separated, one row per soil type
    // Columns: soil_type_id, Ksx, Ksy, Ksz, porosity, theta_r, alpha, n, ha
    // Returns: Map from soil_type_id to soil parameters
    struct SoilParameters {
        Scalar Ksx;      // Saturated hydraulic conductivity in x-direction
        Scalar Ksy;      // Saturated hydraulic conductivity in y-direction
        Scalar Ksz;      // Saturated hydraulic conductivity in z-direction
        Scalar porosity; // Porosity (saturated water content)
        Scalar theta_r;  // Residual water content
        Scalar alpha;    // Van Genuchten alpha parameter
        Scalar n;        // Van Genuchten n parameter
        Scalar ha;       // Air entry value
        
        SoilParameters() : Ksx(1.0e-4), Ksy(1.0e-4), Ksz(1.0e-4),
                          porosity(0.4), theta_r(0.05), alpha(0.01), n(1.5), ha(0.0) {}
    };
    
    static std::map<int, SoilParameters> read_soil_parameter_file(const std::string& filename) {
        std::map<int, SoilParameters> soil_params;
        
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open soil parameter file: " + filename);
        }
        
        std::string line;
        bool first_line = true;
        
        while (std::getline(file, line)) {
            line = trim(line);
            
            // Skip empty lines and comments
            if (line.empty() || line[0] == '#' || 
                (line.length() >= 2 && line.substr(0, 2) == "//")) {
                continue;
            }
            
            // Skip header line if present
            if (first_line) {
                std::string lower = line;
                std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
                if (lower.find("soil_type") != std::string::npos || 
                    lower.find("type_id") != std::string::npos) {
                    first_line = false;
                    continue;
                }
            }
            first_line = false;
            
            // Parse line (CSV or space-separated)
            std::istringstream iss(line);
            std::string token;
            std::vector<std::string> tokens;
            
            // Check if CSV (comma-separated) or space-separated
            bool is_csv = (line.find(',') != std::string::npos);
            
            if (is_csv) {
                while (std::getline(iss, token, ',')) {
                    tokens.push_back(trim(token));
                }
            } else {
                while (iss >> token) {
                    tokens.push_back(token);
                }
            }
            
            if (tokens.size() < 8) {
                throw std::runtime_error("Invalid soil parameter line (need at least 8 values): " + line);
            }
            
            try {
                int soil_type_id = std::stoi(tokens[0]);
                SoilParameters params;
                params.Ksx = std::stod(tokens[1]);
                params.Ksy = std::stod(tokens[2]);
                params.Ksz = std::stod(tokens[3]);
                params.porosity = std::stod(tokens[4]);
                params.theta_r = std::stod(tokens[5]);
                params.alpha = std::stod(tokens[6]);
                params.n = std::stod(tokens[7]);
                
                // Air entry value (optional, default 0.0)
                if (tokens.size() >= 9) {
                    params.ha = std::stod(tokens[8]);
                } else {
                    params.ha = 0.0;
                }
                
                // Compute m = 1 - 1/n (will be computed when setting parameters)
                
                soil_params[soil_type_id] = params;
            } catch (const std::exception& e) {
                throw std::runtime_error("Error parsing soil parameter line: " + line + 
                                       " - " + e.what());
            }
        }
        
        file.close();
        return soil_params;
    }
};

} // namespace Frehg

#endif // FREHG_INPUT_READER_HPP

