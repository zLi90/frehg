#ifndef FREHG_MESH_HPP
#define FREHG_MESH_HPP

#include "define.hpp"
#include <vector>
#include <iostream>

// ============================================================================
//                                MESH SYSTEM
// ============================================================================

class MeshSystem {
public:
    // --- Dimensions ---
    Ordinal nx, ny, nz;           // Grid dimensions
    Ordinal num_cells_2d;         // Total surface cells (nx * ny)
    Ordinal num_cells_3d;         // Total subsurface cells (nx * ny * nz)

    // --- Geometry ---
    Scalar dx, dy;                // Horizontal grid spacing (constant)
    
    // Vertical grid spacing can vary by layer, stored on Device
    View1D<Scalar> dz;            // Size: num_cells_3d (allows fully variable dz)
    View1D<Scalar> z_coords;      // Elevation of cell centers

    // --- Active/Inactive Masks (1 = Active, 0 = Inactive) ---
    View1D<int> sw_active;        // Surface Water active mask
    View1D<int> gw_active;        // Groundwater active mask

    // --- Coordinates / Indices mapping ---
    // Mapping from flattened index 'id' to (i, j, k)
    // Useful for reconstructing position inside a flattened kernel
    View1D<Ordinal> sw_ix, sw_iy;           // Surface Water (i, j)
    View1D<Ordinal> gw_ix, gw_iy, gw_iz;    // Groundwater (i, j, k)

    // --- Connectivity / Topology ---
    // Stores the flattened index of the neighbor.
    // IMPT: Boundary Handling Strategy
    // To avoid "if" statements in GPU kernels, boundary neighbors point to:
    // 1. The cell itself (Clamped/Neumann-0 style default)
    // 2. OR a specific logic handled by the solver.
    // Here we map to the neighbor index if valid, or the current index if boundary.
    
    // Surface Water Connectivity (2D)
    View1D<Ordinal> sw_left;   // i-1
    View1D<Ordinal> sw_right;  // i+1
    View1D<Ordinal> sw_back;   // j-1
    View1D<Ordinal> sw_front;  // j+1

    // Groundwater Connectivity (3D)
    View1D<Ordinal> gw_left;   // i-1
    View1D<Ordinal> gw_right;  // i+1
    View1D<Ordinal> gw_back;   // j-1
    View1D<Ordinal> gw_front;  // j+1
    View1D<Ordinal> gw_bot;    // k-1 (deeper)
    View1D<Ordinal> gw_top;    // k+1 (shallower)

    // --- Coupling Maps ---
    // Links surface 2D cell to the top-most active 3D cell in that column
    View1D<Ordinal> sw_to_gw_idx; 
    // Links 3D cell to surface 2D cell (usually just i + j*nx)
    View1D<Ordinal> gw_to_sw_idx; 

    // Constructor
    MeshSystem(Ordinal _nx, Ordinal _ny, Ordinal _nz, Scalar _dx, Scalar _dy) 
        : nx(_nx), ny(_ny), nz(_nz), dx(_dx), dy(_dy) {
        
        num_cells_2d = nx * ny;
        num_cells_3d = nx * ny * nz;

        // Allocation
        dz = View1D<Scalar>("dz", num_cells_3d);
        z_coords = View1D<Scalar>("z_coords", num_cells_3d);

        sw_active = View1D<int>("sw_active", num_cells_2d);
        gw_active = View1D<int>("gw_active", num_cells_3d);

        sw_ix = View1D<Ordinal>("sw_ix", num_cells_2d);
        sw_iy = View1D<Ordinal>("sw_iy", num_cells_2d);
        
        gw_ix = View1D<Ordinal>("gw_ix", num_cells_3d);
        gw_iy = View1D<Ordinal>("gw_iy", num_cells_3d);
        gw_iz = View1D<Ordinal>("gw_iz", num_cells_3d);

        sw_left = View1D<Ordinal>("sw_left", num_cells_2d);
        sw_right = View1D<Ordinal>("sw_right", num_cells_2d);
        sw_back = View1D<Ordinal>("sw_back", num_cells_2d);
        sw_front = View1D<Ordinal>("sw_front", num_cells_2d);

        gw_left = View1D<Ordinal>("gw_left", num_cells_3d);
        gw_right = View1D<Ordinal>("gw_right", num_cells_3d);
        gw_back = View1D<Ordinal>("gw_back", num_cells_3d);
        gw_front = View1D<Ordinal>("gw_front", num_cells_3d);
        gw_bot = View1D<Ordinal>("gw_bot", num_cells_3d);
        gw_top = View1D<Ordinal>("gw_top", num_cells_3d);

        sw_to_gw_idx = View1D<Ordinal>("sw_to_gw", num_cells_2d);
        gw_to_sw_idx = View1D<Ordinal>("gw_to_sw", num_cells_3d);
    }

    // Initialization Function
    // This runs on the HOST to build the topology, then copies to DEVICE.
    // 'active_mask_input' is a 1D array of size nx*ny (0 or 1).
    // 'dz_layers' is a vector of size nz (thickness of each layer).
    void initialize(const std::vector<int>& active_mask_2d, const std::vector<Scalar>& dz_layers) {
        
        // 1. Create Host Mirrors to populate data locally
        auto h_sw_active = Kokkos::create_mirror_view(sw_active);
        auto h_gw_active = Kokkos::create_mirror_view(gw_active);
        
        auto h_sw_ix = Kokkos::create_mirror_view(sw_ix);
        auto h_sw_iy = Kokkos::create_mirror_view(sw_iy);
        
        auto h_gw_ix = Kokkos::create_mirror_view(gw_ix);
        auto h_gw_iy = Kokkos::create_mirror_view(gw_iy);
        auto h_gw_iz = Kokkos::create_mirror_view(gw_iz);

        auto h_sw_left  = Kokkos::create_mirror_view(sw_left);
        auto h_sw_right = Kokkos::create_mirror_view(sw_right);
        auto h_sw_back  = Kokkos::create_mirror_view(sw_back);
        auto h_sw_front = Kokkos::create_mirror_view(sw_front);

        auto h_gw_left  = Kokkos::create_mirror_view(gw_left);
        auto h_gw_right = Kokkos::create_mirror_view(gw_right);
        auto h_gw_back  = Kokkos::create_mirror_view(gw_back);
        auto h_gw_front = Kokkos::create_mirror_view(gw_front);
        auto h_gw_bot   = Kokkos::create_mirror_view(gw_bot);
        auto h_gw_top   = Kokkos::create_mirror_view(gw_top);

        auto h_dz = Kokkos::create_mirror_view(dz);
        auto h_z_coords = Kokkos::create_mirror_view(z_coords);
        auto h_sw_to_gw = Kokkos::create_mirror_view(sw_to_gw_idx);
        auto h_gw_to_sw = Kokkos::create_mirror_view(gw_to_sw_idx);

        // 2. Build Surface Mesh (2D)
        for (Ordinal j = 0; j < ny; ++j) {
            for (Ordinal i = 0; i < nx; ++i) {
                Ordinal id = i + j * nx; // Flattened index

                // Coordinates
                h_sw_ix(id) = i;
                h_sw_iy(id) = j;

                // Active Flag
                h_sw_active(id) = active_mask_2d[id];

                // Connectivity (Clamped)
                // If at boundary, point to self. Solver handles physics BC.
                h_sw_left(id)  = (i > 0)      ? (id - 1)  : id;
                h_sw_right(id) = (i < nx - 1) ? (id + 1)  : id;
                h_sw_back(id)  = (j > 0)      ? (id - nx) : id;
                h_sw_front(id) = (j < ny - 1) ? (id + nx) : id;
                
                // Initialize coupling map default
                h_sw_to_gw(id) = -1; 
            }
        }

        // 3. Build Subsurface Mesh (3D)
        // k=0 is the bottom layer, k=nz-1 is the top layer
        Scalar current_z_bottom = 0.0; // Assume flat bottom for now
        
        for (Ordinal k = 0; k < nz; ++k) {
            Scalar layer_thickness = dz_layers[k];
            
            for (Ordinal j = 0; j < ny; ++j) {
                for (Ordinal i = 0; i < nx; ++i) {
                    Ordinal sw_id = i + j * nx;             // 2D index
                    Ordinal gw_id = i + j * nx + k * num_cells_2d; // 3D index

                    // Coordinates
                    h_gw_ix(gw_id) = i;
                    h_gw_iy(gw_id) = j;
                    h_gw_iz(gw_id) = k;

                    // Active Flag (Columnar based on surface)
                    h_gw_active(gw_id) = h_sw_active(sw_id);

                    // Vertical Geometry
                    h_dz(gw_id) = layer_thickness;
                    // Simple Z calculation (needs refinement for complex bathymetry)
                    h_z_coords(gw_id) = current_z_bottom + 0.5 * layer_thickness;

                    // Connectivity (Horizontal - Clamped)
                    h_gw_left(gw_id)  = (i > 0)      ? (gw_id - 1)  : gw_id;
                    h_gw_right(gw_id) = (i < nx - 1) ? (gw_id + 1)  : gw_id;
                    h_gw_back(gw_id)  = (j > 0)      ? (gw_id - nx) : gw_id;
                    h_gw_front(gw_id) = (j < ny - 1) ? (gw_id + nx) : gw_id;

                    // Connectivity (Vertical - Clamped)
                    h_gw_bot(gw_id) = (k > 0)      ? (gw_id - num_cells_2d) : gw_id;
                    h_gw_top(gw_id) = (k < nz - 1) ? (gw_id + num_cells_2d) : gw_id;

                    // Coupling Maps
                    h_gw_to_sw(gw_id) = sw_id;
                    if (k == nz - 1) {
                        // This is the top layer
                        h_sw_to_gw(sw_id) = gw_id;
                    }
                }
            }
            current_z_bottom += layer_thickness;
        }

        // 4. Transfer to Device
        Kokkos::deep_copy(sw_active, h_sw_active);
        Kokkos::deep_copy(gw_active, h_gw_active);
        
        Kokkos::deep_copy(sw_ix, h_sw_ix);
        Kokkos::deep_copy(sw_iy, h_sw_iy);
        
        Kokkos::deep_copy(gw_ix, h_gw_ix);
        Kokkos::deep_copy(gw_iy, h_gw_iy);
        Kokkos::deep_copy(gw_iz, h_gw_iz);

        Kokkos::deep_copy(sw_left, h_sw_left);
        Kokkos::deep_copy(sw_right, h_sw_right);
        Kokkos::deep_copy(sw_back, h_sw_back);
        Kokkos::deep_copy(sw_front, h_sw_front);

        Kokkos::deep_copy(gw_left, h_gw_left);
        Kokkos::deep_copy(gw_right, h_gw_right);
        Kokkos::deep_copy(gw_back, h_gw_back);
        Kokkos::deep_copy(gw_front, h_gw_front);
        Kokkos::deep_copy(gw_bot, h_gw_bot);
        Kokkos::deep_copy(gw_top, h_gw_top);

        Kokkos::deep_copy(dz, h_dz);
        Kokkos::deep_copy(z_coords, h_z_coords);
        Kokkos::deep_copy(sw_to_gw_idx, h_sw_to_gw);
        Kokkos::deep_copy(gw_to_sw_idx, h_gw_to_sw);
    }
};

#endif // FREHG_MESH_HPP