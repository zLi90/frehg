#ifndef FREHG_ACTIVE_CELL_MESH_HPP
#define FREHG_ACTIVE_CELL_MESH_HPP

#include "define.hpp"
#include "Domain.hpp"
#include <vector>
#include <unordered_map>
#include <algorithm>

// ============================================================================
//                      ACTIVE CELL MESH SYSTEM
// ============================================================================
// Efficient 1D storage for only active cells with connectivity patterns
// This is essential for sparse matrix assembly in irregular domains

class ActiveCellMesh {
public:
    // --- Dimensions ---
    Ordinal num_active;          // Number of active cells
    
    // --- Index Mapping ---
    // Map from domain index (2D or 3D flattened) to active cell index
    View1D<Ordinal> domain_to_active;  // Size: num_domain_total, -1 if inactive
    
    // Map from active cell index to domain index
    View1D<Ordinal> active_to_domain;  // Size: num_active
    
    // --- Connectivity in Active Space ---
    // Neighbor indices in active cell space (1D)
    // -1 indicates boundary or inactive neighbor
    View1D<Ordinal> neighbor_left;    // i-1 neighbor in active space
    View1D<Ordinal> neighbor_right;   // i+1 neighbor in active space
    View1D<Ordinal> neighbor_back;    // j-1 neighbor in active space
    View1D<Ordinal> neighbor_front;   // j+1 neighbor in active space
    
    // For 3D only:
    View1D<Ordinal> neighbor_bottom;  // k-1 neighbor (deeper)
    View1D<Ordinal> neighbor_top;     // k+1 neighbor (shallower)
    
    // --- Stencil Offsets for Matrix Assembly ---
    // Distance in active space to neighbors (for sparse matrix)
    // Positive = forward, negative = backward, -1 = no neighbor
    View1D<Ordinal> dist_left;        // Offset to left neighbor
    View1D<Ordinal> dist_right;        // Offset to right neighbor
    View1D<Ordinal> dist_back;         // Offset to back neighbor
    View1D<Ordinal> dist_front;        // Offset to front neighbor
    View1D<Ordinal> dist_bottom;       // Offset to bottom neighbor (3D)
    View1D<Ordinal> dist_top;          // Offset to top neighbor (3D)
    
    // --- Original Coordinates ---
    // Store original (i,j) or (i,j,k) coordinates for each active cell
    View1D<Ordinal> coord_i;          // i-coordinate
    View1D<Ordinal> coord_j;          // j-coordinate
    View1D<Ordinal> coord_k;          // k-coordinate (3D only)
    
    // Constructor
    ActiveCellMesh(Ordinal num_domain_total)
        : num_active(0) {
        
        // Allocate domain mapping (will be resized after counting active cells)
        domain_to_active = View1D<Ordinal>("domain_to_active", num_domain_total);
    }
    
    // Build active cell mesh from 2D domain
    void build_from_2d(const SwDomain& domain) {
        build_mapping_2d(domain);
        build_connectivity_2d(domain);
    }
    
    // Build active cell mesh from 3D domain
    void build_from_3d(const GwDomain& domain) {
        build_mapping_3d(domain);
        build_connectivity_3d(domain);
    }
    
private:
    // Build index mapping for 2D domain
    void build_mapping_2d(const SwDomain& domain) {
        auto h_active_mask = Kokkos::create_mirror_view(domain.active_mask);
        Kokkos::deep_copy(h_active_mask, domain.active_mask);
        
        auto h_domain_to_active = Kokkos::create_mirror_view(domain_to_active);
        
        // First pass: count active cells and build domain_to_active map
        num_active = 0;
        for (Ordinal idx = 0; idx < domain.num_cells_total; ++idx) {
            if (h_active_mask(idx) > 0) {
                h_domain_to_active(idx) = num_active;
                num_active++;
            } else {
                h_domain_to_active(idx) = -1;
            }
        }
        
        // Allocate active arrays
        active_to_domain = View1D<Ordinal>("active_to_domain", num_active);
        neighbor_left = View1D<Ordinal>("neighbor_left", num_active);
        neighbor_right = View1D<Ordinal>("neighbor_right", num_active);
        neighbor_back = View1D<Ordinal>("neighbor_back", num_active);
        neighbor_front = View1D<Ordinal>("neighbor_front", num_active);
        dist_left = View1D<Ordinal>("dist_left", num_active);
        dist_right = View1D<Ordinal>("dist_right", num_active);
        dist_back = View1D<Ordinal>("dist_back", num_active);
        dist_front = View1D<Ordinal>("dist_front", num_active);
        coord_i = View1D<Ordinal>("coord_i", num_active);
        coord_j = View1D<Ordinal>("coord_j", num_active);
        
        auto h_active_to_domain = Kokkos::create_mirror_view(active_to_domain);
        auto h_coord_i = Kokkos::create_mirror_view(coord_i);
        auto h_coord_j = Kokkos::create_mirror_view(coord_j);
        
        // Second pass: build active_to_domain map and store coordinates
        Ordinal active_idx = 0;
        for (Ordinal idx = 0; idx < domain.num_cells_total; ++idx) {
            if (h_active_mask(idx) > 0) {
                h_active_to_domain(active_idx) = idx;
                domain.get_2d_index(idx, h_coord_i(active_idx), h_coord_j(active_idx));
                active_idx++;
            }
        }
        
        Kokkos::deep_copy(domain_to_active, h_domain_to_active);
        Kokkos::deep_copy(active_to_domain, h_active_to_domain);
        Kokkos::deep_copy(coord_i, h_coord_i);
        Kokkos::deep_copy(coord_j, h_coord_j);
    }
    
    // Build connectivity for 2D domain
    void build_connectivity_2d(const SwDomain& domain) {
        auto h_neighbor_left = Kokkos::create_mirror_view(neighbor_left);
        auto h_neighbor_right = Kokkos::create_mirror_view(neighbor_right);
        auto h_neighbor_back = Kokkos::create_mirror_view(neighbor_back);
        auto h_neighbor_front = Kokkos::create_mirror_view(neighbor_front);
        auto h_dist_left = Kokkos::create_mirror_view(dist_left);
        auto h_dist_right = Kokkos::create_mirror_view(dist_right);
        auto h_dist_back = Kokkos::create_mirror_view(dist_back);
        auto h_dist_front = Kokkos::create_mirror_view(dist_front);
        auto h_domain_to_active = Kokkos::create_mirror_view(domain_to_active);
        auto h_active_to_domain = Kokkos::create_mirror_view(active_to_domain);
        auto h_active_mask = Kokkos::create_mirror_view(domain.active_mask);
        
        Kokkos::deep_copy(h_domain_to_active, domain_to_active);
        Kokkos::deep_copy(h_active_to_domain, active_to_domain);
        Kokkos::deep_copy(h_active_mask, domain.active_mask);
        
        // Use mesh connectivity if available
        auto h_sw_left = Kokkos::create_mirror_view(domain.mesh->sw_left);
        auto h_sw_right = Kokkos::create_mirror_view(domain.mesh->sw_right);
        auto h_sw_back = Kokkos::create_mirror_view(domain.mesh->sw_back);
        auto h_sw_front = Kokkos::create_mirror_view(domain.mesh->sw_front);
        
        Kokkos::deep_copy(h_sw_left, domain.mesh->sw_left);
        Kokkos::deep_copy(h_sw_right, domain.mesh->sw_right);
        Kokkos::deep_copy(h_sw_back, domain.mesh->sw_back);
        Kokkos::deep_copy(h_sw_front, domain.mesh->sw_front);
        
        for (Ordinal active_idx = 0; active_idx < num_active; ++active_idx) {
            Ordinal domain_idx = h_active_to_domain(active_idx);
            
            // Left neighbor (i-1)
            Ordinal left_domain = h_sw_left(domain_idx);
            if (h_active_mask(left_domain) > 0) {
                Ordinal left_active = h_domain_to_active(left_domain);
                h_neighbor_left(active_idx) = left_active;
                h_dist_left(active_idx) = active_idx - left_active;
            } else {
                h_neighbor_left(active_idx) = -1;
                h_dist_left(active_idx) = -1;
            }
            
            // Right neighbor (i+1)
            Ordinal right_domain = h_sw_right(domain_idx);
            if (h_active_mask(right_domain) > 0) {
                Ordinal right_active = h_domain_to_active(right_domain);
                h_neighbor_right(active_idx) = right_active;
                h_dist_right(active_idx) = right_active - active_idx;
            } else {
                h_neighbor_right(active_idx) = -1;
                h_dist_right(active_idx) = -1;
            }
            
            // Back neighbor (j-1)
            Ordinal back_domain = h_sw_back(domain_idx);
            if (h_active_mask(back_domain) > 0) {
                Ordinal back_active = h_domain_to_active(back_domain);
                h_neighbor_back(active_idx) = back_active;
                h_dist_back(active_idx) = active_idx - back_active;
            } else {
                h_neighbor_back(active_idx) = -1;
                h_dist_back(active_idx) = -1;
            }
            
            // Front neighbor (j+1)
            Ordinal front_domain = h_sw_front(domain_idx);
            if (h_active_mask(front_domain) > 0) {
                Ordinal front_active = h_domain_to_active(front_domain);
                h_neighbor_front(active_idx) = front_active;
                h_dist_front(active_idx) = front_active - active_idx;
            } else {
                h_neighbor_front(active_idx) = -1;
                h_dist_front(active_idx) = -1;
            }
        }
        
        Kokkos::deep_copy(neighbor_left, h_neighbor_left);
        Kokkos::deep_copy(neighbor_right, h_neighbor_right);
        Kokkos::deep_copy(neighbor_back, h_neighbor_back);
        Kokkos::deep_copy(neighbor_front, h_neighbor_front);
        Kokkos::deep_copy(dist_left, h_dist_left);
        Kokkos::deep_copy(dist_right, h_dist_right);
        Kokkos::deep_copy(dist_back, h_dist_back);
        Kokkos::deep_copy(dist_front, h_dist_front);
    }
    
    // Build index mapping for 3D domain
    void build_mapping_3d(const GwDomain& domain) {
        auto h_active_mask_3d = Kokkos::create_mirror_view(domain.active_mask_3d);
        Kokkos::deep_copy(h_active_mask_3d, domain.active_mask_3d);
        
        auto h_domain_to_active = Kokkos::create_mirror_view(domain_to_active);
        
        // First pass: count active cells
        num_active = 0;
        for (Ordinal idx = 0; idx < domain.num_cells_3d_total; ++idx) {
            if (h_active_mask_3d(idx) > 0) {
                h_domain_to_active(idx) = num_active;
                num_active++;
            } else {
                h_domain_to_active(idx) = -1;
            }
        }
        
        // Allocate active arrays
        active_to_domain = View1D<Ordinal>("active_to_domain", num_active);
        neighbor_left = View1D<Ordinal>("neighbor_left", num_active);
        neighbor_right = View1D<Ordinal>("neighbor_right", num_active);
        neighbor_back = View1D<Ordinal>("neighbor_back", num_active);
        neighbor_front = View1D<Ordinal>("neighbor_front", num_active);
        neighbor_bottom = View1D<Ordinal>("neighbor_bottom", num_active);
        neighbor_top = View1D<Ordinal>("neighbor_top", num_active);
        dist_left = View1D<Ordinal>("dist_left", num_active);
        dist_right = View1D<Ordinal>("dist_right", num_active);
        dist_back = View1D<Ordinal>("dist_back", num_active);
        dist_front = View1D<Ordinal>("dist_front", num_active);
        dist_bottom = View1D<Ordinal>("dist_bottom", num_active);
        dist_top = View1D<Ordinal>("dist_top", num_active);
        coord_i = View1D<Ordinal>("coord_i", num_active);
        coord_j = View1D<Ordinal>("coord_j", num_active);
        coord_k = View1D<Ordinal>("coord_k", num_active);
        
        auto h_active_to_domain = Kokkos::create_mirror_view(active_to_domain);
        auto h_coord_i = Kokkos::create_mirror_view(coord_i);
        auto h_coord_j = Kokkos::create_mirror_view(coord_j);
        auto h_coord_k = Kokkos::create_mirror_view(coord_k);
        
        // Second pass: build active_to_domain map
        Ordinal active_idx = 0;
        for (Ordinal idx = 0; idx < domain.num_cells_3d_total; ++idx) {
            if (h_active_mask_3d(idx) > 0) {
                h_active_to_domain(active_idx) = idx;
                domain.get_3d_index(idx, h_coord_i(active_idx), 
                                   h_coord_j(active_idx), h_coord_k(active_idx));
                active_idx++;
            }
        }
        
        Kokkos::deep_copy(domain_to_active, h_domain_to_active);
        Kokkos::deep_copy(active_to_domain, h_active_to_domain);
        Kokkos::deep_copy(coord_i, h_coord_i);
        Kokkos::deep_copy(coord_j, h_coord_j);
        Kokkos::deep_copy(coord_k, h_coord_k);
    }
    
    // Build connectivity for 3D domain
    void build_connectivity_3d(const GwDomain& domain) {
        auto h_neighbor_left = Kokkos::create_mirror_view(neighbor_left);
        auto h_neighbor_right = Kokkos::create_mirror_view(neighbor_right);
        auto h_neighbor_back = Kokkos::create_mirror_view(neighbor_back);
        auto h_neighbor_front = Kokkos::create_mirror_view(neighbor_front);
        auto h_neighbor_bottom = Kokkos::create_mirror_view(neighbor_bottom);
        auto h_neighbor_top = Kokkos::create_mirror_view(neighbor_top);
        auto h_dist_left = Kokkos::create_mirror_view(dist_left);
        auto h_dist_right = Kokkos::create_mirror_view(dist_right);
        auto h_dist_back = Kokkos::create_mirror_view(dist_back);
        auto h_dist_front = Kokkos::create_mirror_view(dist_front);
        auto h_dist_bottom = Kokkos::create_mirror_view(dist_bottom);
        auto h_dist_top = Kokkos::create_mirror_view(dist_top);
        auto h_domain_to_active = Kokkos::create_mirror_view(domain_to_active);
        auto h_active_to_domain = Kokkos::create_mirror_view(active_to_domain);
        auto h_active_mask_3d = Kokkos::create_mirror_view(domain.active_mask_3d);
        
        Kokkos::deep_copy(h_domain_to_active, domain_to_active);
        Kokkos::deep_copy(h_active_to_domain, active_to_domain);
        Kokkos::deep_copy(h_active_mask_3d, domain.active_mask_3d);
        
        // Use mesh connectivity
        auto h_gw_left = Kokkos::create_mirror_view(domain.mesh->gw_left);
        auto h_gw_right = Kokkos::create_mirror_view(domain.mesh->gw_right);
        auto h_gw_back = Kokkos::create_mirror_view(domain.mesh->gw_back);
        auto h_gw_front = Kokkos::create_mirror_view(domain.mesh->gw_front);
        auto h_gw_bot = Kokkos::create_mirror_view(domain.mesh->gw_bot);
        auto h_gw_top = Kokkos::create_mirror_view(domain.mesh->gw_top);
        
        Kokkos::deep_copy(h_gw_left, domain.mesh->gw_left);
        Kokkos::deep_copy(h_gw_right, domain.mesh->gw_right);
        Kokkos::deep_copy(h_gw_back, domain.mesh->gw_back);
        Kokkos::deep_copy(h_gw_front, domain.mesh->gw_front);
        Kokkos::deep_copy(h_gw_bot, domain.mesh->gw_bot);
        Kokkos::deep_copy(h_gw_top, domain.mesh->gw_top);
        
        for (Ordinal active_idx = 0; active_idx < num_active; ++active_idx) {
            Ordinal domain_idx = h_active_to_domain(active_idx);
            
            // Horizontal neighbors (same logic as 2D)
            Ordinal left_domain = h_gw_left(domain_idx);
            if (h_active_mask_3d(left_domain) > 0) {
                Ordinal left_active = h_domain_to_active(left_domain);
                h_neighbor_left(active_idx) = left_active;
                h_dist_left(active_idx) = active_idx - left_active;
            } else {
                h_neighbor_left(active_idx) = -1;
                h_dist_left(active_idx) = -1;
            }
            
            Ordinal right_domain = h_gw_right(domain_idx);
            if (h_active_mask_3d(right_domain) > 0) {
                Ordinal right_active = h_domain_to_active(right_domain);
                h_neighbor_right(active_idx) = right_active;
                h_dist_right(active_idx) = right_active - active_idx;
            } else {
                h_neighbor_right(active_idx) = -1;
                h_dist_right(active_idx) = -1;
            }
            
            Ordinal back_domain = h_gw_back(domain_idx);
            if (h_active_mask_3d(back_domain) > 0) {
                Ordinal back_active = h_domain_to_active(back_domain);
                h_neighbor_back(active_idx) = back_active;
                h_dist_back(active_idx) = active_idx - back_active;
            } else {
                h_neighbor_back(active_idx) = -1;
                h_dist_back(active_idx) = -1;
            }
            
            Ordinal front_domain = h_gw_front(domain_idx);
            if (h_active_mask_3d(front_domain) > 0) {
                Ordinal front_active = h_domain_to_active(front_domain);
                h_neighbor_front(active_idx) = front_active;
                h_dist_front(active_idx) = front_active - active_idx;
            } else {
                h_neighbor_front(active_idx) = -1;
                h_dist_front(active_idx) = -1;
            }
            
            // Vertical neighbors
            Ordinal bottom_domain = h_gw_bot(domain_idx);
            if (h_active_mask_3d(bottom_domain) > 0) {
                Ordinal bottom_active = h_domain_to_active(bottom_domain);
                h_neighbor_bottom(active_idx) = bottom_active;
                h_dist_bottom(active_idx) = active_idx - bottom_active;
            } else {
                h_neighbor_bottom(active_idx) = -1;
                h_dist_bottom(active_idx) = -1;
            }
            
            Ordinal top_domain = h_gw_top(domain_idx);
            if (h_active_mask_3d(top_domain) > 0) {
                Ordinal top_active = h_domain_to_active(top_domain);
                h_neighbor_top(active_idx) = top_active;
                h_dist_top(active_idx) = top_active - active_idx;
            } else {
                h_neighbor_top(active_idx) = -1;
                h_dist_top(active_idx) = -1;
            }
        }
        
        Kokkos::deep_copy(neighbor_left, h_neighbor_left);
        Kokkos::deep_copy(neighbor_right, h_neighbor_right);
        Kokkos::deep_copy(neighbor_back, h_neighbor_back);
        Kokkos::deep_copy(neighbor_front, h_neighbor_front);
        Kokkos::deep_copy(neighbor_bottom, h_neighbor_bottom);
        Kokkos::deep_copy(neighbor_top, h_neighbor_top);
        Kokkos::deep_copy(dist_left, h_dist_left);
        Kokkos::deep_copy(dist_right, h_dist_right);
        Kokkos::deep_copy(dist_back, h_dist_back);
        Kokkos::deep_copy(dist_front, h_dist_front);
        Kokkos::deep_copy(dist_bottom, h_dist_bottom);
        Kokkos::deep_copy(dist_top, h_dist_top);
    }
};

#endif // FREHG_ACTIVE_CELL_MESH_HPP

