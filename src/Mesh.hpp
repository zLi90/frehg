#ifndef MESH_HPP
#define MESH_HPP

#include "Types.hpp"
#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm> // for std::max

// -----------------------------------------------------------------------------
// Direction Enums for Readability
// -----------------------------------------------------------------------------
enum { EAST=0, WEST=1, NORTH=2, SOUTH=3, TOP=4, BOTTOM=5, NUM_DIRS=6 };

class UnstructuredMesh {
public:
    std::string name;
    int spatial_dim;      // 2 for Surface, 3 for Subsurface
    int num_active_cells; // The "Packed" size (N)

    // Original Structured Dimensions (for reference/IO)
    int nx_global, ny_global, nz_global;

    // -------------------------------------------------------------------------
    // 1. Connectivity & Geometry (Device Views)
    // -------------------------------------------------------------------------
    // Neighbors: Stores the ID of the neighbor. -1 if boundary/inactive.
    // Dimensions: [num_active_cells][6] (Even for 2D, we keep 6 for consistency, Top/Bot are -1)
    View2D<int> neighbors;

    // Geometry: Structure of Arrays
    View1D<Real> center_x;
    View1D<Real> center_y;
    View1D<Real> center_z;
    View1D<Real> dz;      // Vertical thickness
    View1D<Real> volume;  // Volume (3D) or Area (2D)

    // -------------------------------------------------------------------------
    // 2. Coupling & Boundaries
    // -------------------------------------------------------------------------
    // For Surface Mesh: Stores the ID of the top-most active Subsurface cell below it.
    // For Subsurface Mesh: Stores the ID of the Surface cell above it (if it is a top cell).
    // Value is -1 if no connection exists.
    View1D<int> coupling_index;

    // Boundary Flags (useful for quick logic in solvers)
    View1D<int> is_boundary; 

    // -------------------------------------------------------------------------
    // 3. Constructor
    // -------------------------------------------------------------------------
    UnstructuredMesh(std::string _name, int n_cells, int dim) 
        : name(_name), num_active_cells(n_cells), spatial_dim(dim) 
    {
        // Allocate Device Views
        neighbors      = View2D<int>(name + "_neighbors", n_cells, NUM_DIRS);
        center_x       = View1D<Real>(name + "_x", n_cells);
        center_y       = View1D<Real>(name + "_y", n_cells);
        center_z       = View1D<Real>(name + "_z", n_cells);
        dz             = View1D<Real>(name + "_dz", n_cells);
        volume         = View1D<Real>(name + "_vol", n_cells);
        coupling_index = View1D<int>(name + "_couple", n_cells);
        is_boundary    = View1D<int>(name + "_is_bnd", n_cells);
    }

    // -------------------------------------------------------------------------
    // 4. Builder: Surface Domain (2D)
    // -------------------------------------------------------------------------
    // Logic adapted from 'build_surf_map'
    void init_surface(int nx, int ny, double dx, double dy, 
                      const std::vector<int>& active_mask) 
    {
        nx_global = nx; ny_global = ny; nz_global = 1;

        // A. Host Mirrors for Initialization
        auto h_neighbors = Kokkos::create_mirror_view(neighbors);
        auto h_x         = Kokkos::create_mirror_view(center_x);
        auto h_y         = Kokkos::create_mirror_view(center_y);
        auto h_z         = Kokkos::create_mirror_view(center_z); // Typically 0 or surface elevation
        auto h_vol       = Kokkos::create_mirror_view(volume);
        auto h_couple    = Kokkos::create_mirror_view(coupling_index);

        // B. Mapping: Global (i,j) -> Local (cell_id)
        std::vector<int> g2l(nx * ny, -1);
        std::vector<int> l2g(num_active_cells);
        
        int count = 0;
        for (int i = 0; i < nx * ny; i++) {
            if (active_mask[i] > 0) {
                if (count >= num_active_cells) break; // Safety
                g2l[i] = count;
                l2g[count] = i;
                count++;
            }
        }

        // C. Build Connectivity
        for (int id = 0; id < num_active_cells; id++) {
            int global = l2g[id];
            int j = global / nx;
            int i = global % nx;

            // Geometry
            h_x(id) = (i + 0.5) * dx;
            h_y(id) = (j + 0.5) * dy;
            h_z(id) = 0.0;     // Placeholder, update if you have elevation data
            h_vol(id) = dx * dy;
            h_couple(id) = -1; // Will be filled by coupling step later

            // Neighbors (Standard 5-point stencil logic)
            h_neighbors(id, EAST)  = (i < nx - 1) ? g2l[global + 1] : -1;
            h_neighbors(id, WEST)  = (i > 0)      ? g2l[global - 1] : -1;
            h_neighbors(id, NORTH) = (j < ny - 1) ? g2l[global + nx] : -1;
            h_neighbors(id, SOUTH) = (j > 0)      ? g2l[global - nx] : -1;
            h_neighbors(id, TOP)    = -1; // No vertical neighbors in 2D
            h_neighbors(id, BOTTOM) = -1;
        }

        // D. deep_copy to Device
        Kokkos::deep_copy(neighbors, h_neighbors);
        Kokkos::deep_copy(center_x, h_x);
        Kokkos::deep_copy(center_y, h_y);
        Kokkos::deep_copy(center_z, h_z);
        Kokkos::deep_copy(volume, h_vol);
        Kokkos::deep_copy(coupling_index, h_couple);
        
        std::cout << ">> Surface Mesh Initialized: " << num_active_cells << " active cells." << std::endl;
    }

    // -------------------------------------------------------------------------
    // 5. Builder: Subsurface Domain (3D)
    // -------------------------------------------------------------------------
    // Logic adapted from 'build_subsurf_map' including terrain following
    void init_subsurface(int nx, int ny, int nz, 
                         double dx, double dy, double botZ, 
                         const std::vector<double>& bathymetry,
                         const std::vector<int>& surf_mask) 
    {
        nx_global = nx; ny_global = ny; nz_global = nz;

        // A. Host Mirrors
        auto h_neighbors = Kokkos::create_mirror_view(neighbors);
        auto h_x         = Kokkos::create_mirror_view(center_x);
        auto h_y         = Kokkos::create_mirror_view(center_y);
        auto h_z         = Kokkos::create_mirror_view(center_z);
        auto h_dz        = Kokkos::create_mirror_view(dz);
        auto h_vol       = Kokkos::create_mirror_view(volume);
        auto h_couple    = Kokkos::create_mirror_view(coupling_index);
        auto h_bnd       = Kokkos::create_mirror_view(is_boundary);

        // B. Geometry & Mapping Pre-Calculation
        // We must replicate the 'map.c' logic where inactive cells (above bathymetry) are skipped.
        int total_cells = nx * ny * nz;
        std::vector<int> g2l(total_cells, -1);
        std::vector<int> l2g(num_active_cells);
        
        // Temporary storage for calculated geometry before packing
        std::vector<double> tmp_z(total_cells, 0.0);
        std::vector<double> tmp_dz(total_cells, 0.0);
        
        int count = 0;

        // Loop over columns
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int surf_idx = j * nx + i;
                
                // If surface is inactive, subsurface column is likely inactive too (check map.c logic)
                if (surf_mask[surf_idx] <= 0) continue;

                double surf_elev = bathymetry[surf_idx];
                double total_depth = surf_elev - botZ;
                
                // Simple Sigma Stretch (Uniform dz distribution) - Matches map.c 'follow_terrain'
                double col_dz = total_depth / nz;

                // Loop layers (k=0 is bottom in map.c usually, check config. Here assuming k=0 is bottom)
                for (int k = 0; k < nz; k++) {
                    int global_idx = k * (nx * ny) + j * nx + i;

                    // Calculate Geometry
                    double cell_bot = botZ + k * col_dz;
                    tmp_dz[global_idx] = col_dz;
                    tmp_z[global_idx]  = cell_bot + 0.5 * col_dz;

                    // Register Active Cell
                    if (col_dz > 1e-6) {
                        if (count < num_active_cells) {
                            g2l[global_idx] = count;
                            l2g[count] = global_idx;
                            count++;
                        }
                    }
                }
            }
        }

        // C. Build Connectivity (Packed)
        for (int id = 0; id < num_active_cells; id++) {
            int global = l2g[id];
            
            // Decode i, j, k
            int k = global / (nx * ny);
            int rem = global % (nx * ny);
            int j = rem / nx;
            int i = rem % nx;

            h_x(id)  = (i + 0.5) * dx;
            h_y(id)  = (j + 0.5) * dy;
            h_z(id)  = tmp_z[global];
            h_dz(id) = tmp_dz[global];
            h_vol(id)= dx * dy * h_dz(id);
            h_couple(id) = -1; // Default
            h_bnd(id) = 0;

            // --- Horizontal Neighbors ---
            h_neighbors(id, EAST)  = (i < nx - 1) ? g2l[global + 1] : -1;
            h_neighbors(id, WEST)  = (i > 0)      ? g2l[global - 1] : -1;
            h_neighbors(id, NORTH) = (j < ny - 1) ? g2l[global + nx] : -1;
            h_neighbors(id, SOUTH) = (j > 0)      ? g2l[global - nx] : -1;

            // --- Vertical Neighbors ---
            // TOP (k+1)
            if (k < nz - 1) {
                h_neighbors(id, TOP) = g2l[global + (nx * ny)];
            } else {
                h_neighbors(id, TOP) = -1; 
                h_bnd(id) = 1; // Mark as Top Boundary Cell (Coupling Interface)
            }

            // BOTTOM (k-1)
            if (k > 0) {
                h_neighbors(id, BOTTOM) = g2l[global - (nx * ny)];
            } else {
                h_neighbors(id, BOTTOM) = -1;
            }
        }

        // D. Upload
        Kokkos::deep_copy(neighbors, h_neighbors);
        Kokkos::deep_copy(center_x, h_x);
        Kokkos::deep_copy(center_y, h_y);
        Kokkos::deep_copy(center_z, h_z);
        Kokkos::deep_copy(dz, h_dz);
        Kokkos::deep_copy(volume, h_vol);
        Kokkos::deep_copy(coupling_index, h_couple);
        Kokkos::deep_copy(is_boundary, h_bnd);

        std::cout << ">> Subsurface Mesh Initialized: " << num_active_cells << " active cells." << std::endl;
    }
};

#endif // MESH_HPP