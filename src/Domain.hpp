#ifndef FREHG_DOMAIN_HPP
#define FREHG_DOMAIN_HPP

#include "define.hpp"
#include "mesh.hpp"
#include <vector>
#include <string>
#include <memory>
#include <algorithm>

// ============================================================================
//                           BASE DOMAIN CLASS
// ============================================================================
// Shared properties for both surface water and groundwater domains

class BaseDomain {
public:
    // --- Grid Dimensions ---
    Ordinal nx, ny;              // Horizontal grid dimensions
    Scalar dx, dy;               // Horizontal grid spacing
    
    // --- Domain Extent ---
    Scalar x_min, x_max;         // Domain bounds in x-direction
    Scalar y_min, y_max;         // Domain bounds in y-direction
    
    // --- Total Cell Counts (including inactive) ---
    Ordinal num_cells_total;     // Total cells in rectangular grid
    
    // --- Active Cell Information ---
    Ordinal num_cells_active;    // Number of active cells
    View1D<int> active_mask;     // Active mask (1=active, 0=inactive) on device
    
    // --- Bathymetry/Elevation Data ---
    View1D<Scalar> bathymetry;   // Bathymetry/elevation data (DEM)
    
    // --- Mesh System ---
    std::shared_ptr<MeshSystem> mesh;  // Shared mesh system
    
    // Constructor
    BaseDomain(Ordinal _nx, Ordinal _ny, Scalar _dx, Scalar _dy)
        : nx(_nx), ny(_ny), dx(_dx), dy(_dy),
          x_min(0.0), x_max(_nx * _dx),
          y_min(0.0), y_max(_ny * _dy),
          num_cells_total(_nx * _ny),
          num_cells_active(0) {
        
        // Allocate views
        active_mask = View1D<int>("active_mask", num_cells_total);
        bathymetry = View1D<Scalar>("bathymetry", num_cells_total);
    }
    
    // Virtual destructor for polymorphism
    virtual ~BaseDomain() = default;
    
    // Initialize from DEM data
    // active_mask_input: 1D array of size nx*ny (1=active, 0=inactive)
    // bathymetry_input: 1D array of size nx*ny (elevation values)
    virtual void initialize(const std::vector<int>& active_mask_input,
                          const std::vector<Scalar>& bathymetry_input) {
        
        // Create host mirrors
        auto h_active = Kokkos::create_mirror_view(active_mask);
        auto h_bath = Kokkos::create_mirror_view(bathymetry);
        
        // Count active cells and copy data
        num_cells_active = 0;
        for (Ordinal i = 0; i < num_cells_total; ++i) {
            h_active(i) = active_mask_input[i];
            h_bath(i) = bathymetry_input[i];
            if (active_mask_input[i] > 0) {
                num_cells_active++;
            }
        }
        
        // Transfer to device
        Kokkos::deep_copy(active_mask, h_active);
        Kokkos::deep_copy(bathymetry, h_bath);
    }
    
    // Get 2D index from flattened index
    KOKKOS_INLINE_FUNCTION
    void get_2d_index(Ordinal flat_idx, Ordinal& i, Ordinal& j) const {
        j = flat_idx / nx;
        i = flat_idx - j * nx;
    }
    
    // Get flattened index from 2D coordinates
    KOKKOS_INLINE_FUNCTION
    Ordinal get_flat_index(Ordinal i, Ordinal j) const {
        return i + j * nx;
    }
    
    // Check if cell is active
    KOKKOS_INLINE_FUNCTION
    bool is_active(Ordinal flat_idx) const {
        return active_mask(flat_idx) > 0;
    }
};

// ============================================================================
//                      SURFACE WATER DOMAIN (2D)
// ============================================================================

class SwDomain : public BaseDomain {
public:
    // Additional 2D-specific properties can be added here
    
    SwDomain(Ordinal _nx, Ordinal _ny, Scalar _dx, Scalar _dy)
        : BaseDomain(_nx, _ny, _dx, _dy) {
        
        // Create mesh system (2D only, nz=1 for surface)
        mesh = std::make_shared<MeshSystem>(nx, ny, 1, dx, dy);
    }
    
    // Initialize surface water domain
    void initialize(const std::vector<int>& active_mask_input,
                   const std::vector<Scalar>& bathymetry_input) override {
        
        // Initialize base domain
        BaseDomain::initialize(active_mask_input, bathymetry_input);
        
        // Initialize mesh with surface water data
        // For 2D surface, we use a single layer (nz=1)
        std::vector<Scalar> dz_layers = {1.0}; // Dummy layer thickness
        mesh->initialize(active_mask_input, dz_layers);
    }
    
    // Get domain name
    std::string get_name() const { return "SurfaceWater"; }
};

// ============================================================================
//                      GROUNDWATER DOMAIN (3D)
// ============================================================================

class GwDomain : public BaseDomain {
public:
    // --- 3D-Specific Properties ---
    Ordinal nz;                  // Number of vertical layers
    Scalar dz;                   // Base vertical spacing (may vary)
    Scalar dz_increment;         // Layer thickness increment factor
    Scalar bottom_z;             // Bottom elevation of domain
    bool follow_terrain;         // Whether to follow terrain or use regular grid
    
    // --- Vertical Layer Information ---
    View1D<Scalar> layer_bottoms;  // Bottom elevation of each layer
    View1D<Scalar> layer_thickness; // Thickness of each layer
    
    // --- 3D Active Mask ---
    View1D<int> active_mask_3d;    // 3D active mask (nx*ny*nz)
    Ordinal num_cells_3d_total;     // Total 3D cells
    Ordinal num_cells_3d_active;     // Active 3D cells
    
    GwDomain(Ordinal _nx, Ordinal _ny, Ordinal _nz, 
             Scalar _dx, Scalar _dy, Scalar _dz,
             Scalar _bottom_z = 0.0, Scalar _dz_increment = 1.0,
             bool _follow_terrain = false)
        : BaseDomain(_nx, _ny, _dx, _dy),
          nz(_nz), dz(_dz), dz_increment(_dz_increment),
          bottom_z(_bottom_z), follow_terrain(_follow_terrain),
          num_cells_3d_total(_nx * _ny * _nz),
          num_cells_3d_active(0) {
        
        // Create mesh system (3D)
        mesh = std::make_shared<MeshSystem>(nx, ny, nz, dx, dy);
        
        // Allocate 3D views
        active_mask_3d = View1D<int>("active_mask_3d", num_cells_3d_total);
        layer_bottoms = View1D<Scalar>("layer_bottoms", nz);
        layer_thickness = View1D<Scalar>("layer_thickness", nz);
    }
    
    // Initialize groundwater domain
    void initialize(const std::vector<int>& active_mask_2d,
                   const std::vector<Scalar>& bathymetry_input) override {
        
        // Initialize base domain (2D surface)
        BaseDomain::initialize(active_mask_2d, bathymetry_input);
        
        // Build 3D layers
        build_3d_layers(bathymetry_input);
        
        // Build 3D active mask
        build_3d_active_mask(active_mask_2d, bathymetry_input);
        
        // Initialize mesh with 3D data
        std::vector<Scalar> dz_layers(nz);
        auto h_layer_thickness = Kokkos::create_mirror_view(layer_thickness);
        Kokkos::deep_copy(h_layer_thickness, layer_thickness);
        for (Ordinal k = 0; k < nz; ++k) {
            dz_layers[k] = h_layer_thickness(k);
        }
        mesh->initialize(active_mask_2d, dz_layers);
    }
    
    // Get 3D index from flattened index
    KOKKOS_INLINE_FUNCTION
    void get_3d_index(Ordinal flat_idx, Ordinal& i, Ordinal& j, Ordinal& k) const {
        Ordinal cells_per_layer = nx * ny;
        k = flat_idx / cells_per_layer;
        Ordinal layer_idx = flat_idx - k * cells_per_layer;
        j = layer_idx / nx;
        i = layer_idx - j * nx;
    }
    
    // Get flattened index from 3D coordinates
    KOKKOS_INLINE_FUNCTION
    Ordinal get_flat_index_3d(Ordinal i, Ordinal j, Ordinal k) const {
        return i + j * nx + k * nx * ny;
    }
    
    // Check if 3D cell is active
    KOKKOS_INLINE_FUNCTION
    bool is_active_3d(Ordinal flat_idx) const {
        return active_mask_3d(flat_idx) > 0;
    }
    
    // Get domain name
    std::string get_name() const { return "Groundwater"; }
    
private:
    // Build 3D layer structure
    void build_3d_layers(const std::vector<Scalar>& bathymetry) {
        auto h_layer_bottoms = Kokkos::create_mirror_view(layer_bottoms);
        auto h_layer_thickness = Kokkos::create_mirror_view(layer_thickness);
        
        // Find max and min bathymetry
        Scalar max_bath = *std::max_element(bathymetry.begin(), bathymetry.end());
        Scalar min_bath = *std::min_element(bathymetry.begin(), bathymetry.end());
        
        if (follow_terrain) {
            // Terrain-following: layers adapt to local bathymetry
            // For now, use uniform layers - can be enhanced later
            Scalar total_depth = max_bath - bottom_z;
            Scalar avg_dz = total_depth / nz;
            
            for (Ordinal k = 0; k < nz; ++k) {
                h_layer_thickness(k) = avg_dz;
                h_layer_bottoms(k) = max_bath - (k + 1) * avg_dz;
            }
        } else {
            // Regular Cartesian grid
            if (dz_increment == 1.0) {
                // Uniform layers
                for (Ordinal k = 0; k < nz; ++k) {
                    h_layer_thickness(k) = dz;
                    h_layer_bottoms(k) = max_bath - (k + 1) * dz;
                }
            } else {
                // Variable thickness layers
                Scalar current_dz = dz;
                Scalar current_bottom = max_bath - dz;
                h_layer_bottoms(0) = current_bottom;
                h_layer_thickness(0) = dz;
                
                for (Ordinal k = 1; k < nz; ++k) {
                    current_dz *= dz_increment;
                    current_bottom -= current_dz;
                    h_layer_bottoms(k) = current_bottom;
                    h_layer_thickness(k) = current_dz;
                }
            }
        }
        
        Kokkos::deep_copy(layer_bottoms, h_layer_bottoms);
        Kokkos::deep_copy(layer_thickness, h_layer_thickness);
    }
    
    // Build 3D active mask based on bathymetry
    void build_3d_active_mask(const std::vector<int>& active_2d,
                              const std::vector<Scalar>& bathymetry) {
        auto h_active_3d = Kokkos::create_mirror_view(active_mask_3d);
        auto h_layer_bottoms = Kokkos::create_mirror_view(layer_bottoms);
        auto h_layer_thickness = Kokkos::create_mirror_view(layer_thickness);
        
        Kokkos::deep_copy(h_layer_bottoms, layer_bottoms);
        Kokkos::deep_copy(h_layer_thickness, layer_thickness);
        
        num_cells_3d_active = 0;
        
        for (Ordinal k = 0; k < nz; ++k) {
            Scalar layer_bottom = h_layer_bottoms(k);
            Scalar layer_top = layer_bottom + h_layer_thickness(k);
            
            for (Ordinal j = 0; j < ny; ++j) {
                for (Ordinal i = 0; i < nx; ++i) {
                    Ordinal sw_idx = i + j * nx;
                    Ordinal gw_idx = i + j * nx + k * nx * ny;
                    
                    // Cell is active if:
                    // 1. Surface cell is active
                    // 2. Layer bottom is below bathymetry
                    if (active_2d[sw_idx] > 0 && layer_bottom < bathymetry[sw_idx]) {
                        h_active_3d(gw_idx) = 1;
                        num_cells_3d_active++;
                    } else {
                        h_active_3d(gw_idx) = 0;
                    }
                }
            }
        }
        
        Kokkos::deep_copy(active_mask_3d, h_active_3d);
    }
};

#endif // FREHG_DOMAIN_HPP

