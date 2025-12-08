#ifndef FREHG_PCG_SOLVER_HPP
#define FREHG_PCG_SOLVER_HPP

#include "define.hpp"
#include "ActiveCellMesh.hpp"
#include <cmath>
#include <iostream>

// ============================================================================
//                      PRECONDITIONED CONJUGATE GRADIENT SOLVER
// ============================================================================
// Kokkos-based PCG solver for sparse linear systems
// Works with stencil-based matrices (no need for full sparse matrix storage)

enum class PreconditionerType {
    NONE,      // No preconditioner (identity)
    JACOBI,    // Diagonal/Jacobi preconditioner
    SSOR       // Symmetric Successive Over-Relaxation
};

template <int Dim>
class PCGSolver {
private:
    Ordinal num_active;
    PreconditionerType precond_type;
    Scalar tolerance;
    Ordinal max_iterations;
    Scalar omega;  // Relaxation parameter for SSOR (typically 1.0-1.5)
    
    // Working vectors
    View1D<Scalar> r;      // Residual
    View1D<Scalar> z;      // Preconditioned residual
    View1D<Scalar> p;      // Search direction
    View1D<Scalar> Ap;      // Matrix-vector product A*p
    
    // Preconditioner storage
    View1D<Scalar> M_inv;  // Inverse diagonal (for Jacobi/SSOR)
    
    // Statistics
    Ordinal iterations;
    Scalar final_residual;
    bool converged;
    
public:
    PCGSolver(Ordinal num_active_cells, 
              PreconditionerType prec = PreconditionerType::JACOBI,
              Scalar tol = 1.0e-8,
              Ordinal max_iter = 10000,
              Scalar ssor_omega = 1.0)
        : num_active(num_active_cells),
          precond_type(prec),
          tolerance(tol),
          max_iterations(max_iter),
          omega(ssor_omega),
          iterations(0),
          final_residual(0.0),
          converged(false) {
        
        // Allocate working vectors
        r = View1D<Scalar>("r", num_active);
        z = View1D<Scalar>("z", num_active);
        p = View1D<Scalar>("p", num_active);
        Ap = View1D<Scalar>("Ap", num_active);
        M_inv = View1D<Scalar>("M_inv", num_active);
    }
    
    // Solve Ax = b for 2D system (surface water)
    // Matrix coefficients: diag, xp, xm, yp, ym
    // ActiveCellMesh provides neighbor connectivity
    // x: Solution vector (input: initial guess, output: solution)
    // use_initial_guess: If true, uses x as initial guess; if false, starts from zero
    void solve_2d(const ActiveCellMesh& mesh,
                  const View1D<Scalar>& diag,      // Diagonal (Sct)
                  const View1D<Scalar>& xp,        // X+ coefficient (Sxp)
                  const View1D<Scalar>& xm,        // X- coefficient (Sxm)
                  const View1D<Scalar>& yp,        // Y+ coefficient (Syp)
                  const View1D<Scalar>& ym,        // Y- coefficient (Sym)
                  const View1D<Scalar>& rhs,       // Right-hand side (Srhs)
                  View1D<Scalar>& x,               // Solution (input/output)
                  bool use_initial_guess = false) { // Use x as initial guess (default: no for 2D)
        
        // Build preconditioner
        build_preconditioner_2d(mesh, diag, xp, xm, yp, ym);
        
        // Compute initial residual: r = b - A*x
        if (use_initial_guess) {
            // Compute r = b - A*x using the provided initial guess
            matrix_vector_product_2d(mesh, diag, xp, xm, yp, ym, x, r);
            auto _r = r;
            auto _rhs = rhs;
            Kokkos::parallel_for(RangePolicy(0, num_active),
                KOKKOS_LAMBDA (const Ordinal i) {
                    _r(i) = _rhs(i) - _r(i);
                });
        } else {
            // Zero initial guess: r = b - A*0 = b
            Kokkos::deep_copy(x, 0.0);
            Kokkos::deep_copy(r, rhs);
        }
        
        // Compute initial residual norm
        Scalar rnorm0 = compute_norm(r);
        if (rnorm0 < tolerance) {
            converged = true;
            iterations = 0;
            final_residual = rnorm0;
            return;
        }
        
        // Apply preconditioner: z = M^{-1} * r
        apply_preconditioner_2d(mesh, diag, xp, xm, yp, ym, r, z);
        
        // Initialize search direction: p = z
        Kokkos::deep_copy(p, z);
        
        // Compute initial dot product: rz = r^T * z
        Scalar rz_old = dot_product(r, z);
        
        // Main PCG iteration
        converged = false;
        for (iterations = 0; iterations < max_iterations; ++iterations) {
            // Compute Ap = A * p
            matrix_vector_product_2d(mesh, diag, xp, xm, yp, ym, p, Ap);
            
            // Compute alpha = (r^T * z) / (p^T * Ap)
            Scalar pAp = dot_product(p, Ap);
            if (std::abs(pAp) < 1.0e-20) {
                // Breakdown - p is in null space
                break;
            }
            Scalar alpha = rz_old / pAp;
            
            // Update solution: x = x + alpha * p
            update_vector(x, p, alpha);
            
            // Update residual: r = r - alpha * Ap
            update_vector(r, Ap, -alpha);
            
            // Check convergence
            Scalar rnorm = compute_norm(r);
            final_residual = rnorm;
            
            if (rnorm < tolerance * rnorm0) {
                converged = true;
                break;
            }
            
            // Apply preconditioner: z = M^{-1} * r
            apply_preconditioner_2d(mesh, diag, xp, xm, yp, ym, r, z);
            
            // Compute beta = (r_new^T * z_new) / (r_old^T * z_old)
            Scalar rz_new = dot_product(r, z);
            Scalar beta = rz_new / rz_old;
            
            // Update search direction: p = z + beta * p
            update_vector(p, z, 1.0, beta);
            
            rz_old = rz_new;
        }
    }
    
    // Solve Ax = b for 3D system (groundwater)
    // Matrix coefficients: diag, xp, xm, yp, ym, zp, zm
    // x: Solution vector (input: initial guess, output: solution)
    // use_initial_guess: If true, uses x as initial guess; if false, starts from zero
    void solve_3d(const ActiveCellMesh& mesh,
                  const View1D<Scalar>& diag,      // Diagonal (Gct) - indexed by domain
                  const View1D<Scalar>& xp,        // X+ coefficient (Gxp) - indexed by domain
                  const View1D<Scalar>& xm,        // X- coefficient (Gxm) - indexed by domain
                  const View1D<Scalar>& yp,        // Y+ coefficient (Gyp) - indexed by domain
                  const View1D<Scalar>& ym,        // Y- coefficient (Gym) - indexed by domain
                  const View1D<Scalar>& zp,        // Z+ coefficient (Gzp) - indexed by domain
                  const View1D<Scalar>& zm,        // Z- coefficient (Gzm) - indexed by domain
                  const View1D<Scalar>& rhs,       // Right-hand side (Grhs) - indexed by domain
                  View1D<Scalar>& x,               // Solution (input/output) - indexed by active
                  bool use_initial_guess = true) { // Use x as initial guess
        
        // Build preconditioner
        build_preconditioner_3d(mesh, diag, xp, xm, yp, ym, zp, zm);
        
        auto _active_to_domain = mesh.active_to_domain;
        
        // Compute initial residual: r = b - A*x
        if (use_initial_guess) {
            // Compute r = b - A*x using the provided initial guess
            matrix_vector_product_3d(mesh, diag, xp, xm, yp, ym, zp, zm, x, r);
            auto _r = r;
            auto _rhs = rhs;
            Kokkos::parallel_for(RangePolicy(0, num_active),
                KOKKOS_LAMBDA (const Ordinal i) {
                    Ordinal domain_idx = _active_to_domain(i);
                    _r(i) = _rhs(domain_idx) - _r(i);
                });
        } else {
            // Zero initial guess: r = b - A*0 = b
            // Copy RHS with index mapping (rhs is indexed by domain, r by active)
            Kokkos::deep_copy(x, 0.0);
            auto _r = r;
            auto _rhs = rhs;
            Kokkos::parallel_for(RangePolicy(0, num_active),
                KOKKOS_LAMBDA (const Ordinal i) {
                    Ordinal domain_idx = _active_to_domain(i);
                    _r(i) = _rhs(domain_idx);
                });
        }
        
        // Compute initial residual norm
        Scalar rnorm0 = compute_norm(r);
        if (rnorm0 < tolerance) {
            converged = true;
            iterations = 0;
            final_residual = rnorm0;
            return;
        }
        
        // Apply preconditioner: z = M^{-1} * r
        apply_preconditioner_3d(mesh, diag, xp, xm, yp, ym, zp, zm, r, z);
        
        // Initialize search direction: p = z
        Kokkos::deep_copy(p, z);
        
        // Compute initial dot product: rz = r^T * z
        Scalar rz_old = dot_product(r, z);
        
        // Main PCG iteration
        converged = false;
        for (iterations = 0; iterations < max_iterations; ++iterations) {
            // Compute Ap = A * p
            matrix_vector_product_3d(mesh, diag, xp, xm, yp, ym, zp, zm, p, Ap);
            
            // Compute alpha = (r^T * z) / (p^T * Ap)
            Scalar pAp = dot_product(p, Ap);
            if (std::abs(pAp) < 1.0e-20) {
                // Breakdown - p is in null space
                break;
            }
            Scalar alpha = rz_old / pAp;
            
            // Update solution: x = x + alpha * p
            update_vector(x, p, alpha);
            
            // Update residual: r = r - alpha * Ap
            update_vector(r, Ap, -alpha);
            
            // Check convergence
            Scalar rnorm = compute_norm(r);
            final_residual = rnorm;
            
            if (rnorm < tolerance * rnorm0) {
                converged = true;
                break;
            }
            
            // Check for NaN
            if (std::isnan(rnorm)) {
                break;
            }
            
            // Apply preconditioner: z = M^{-1} * r
            apply_preconditioner_3d(mesh, diag, xp, xm, yp, ym, zp, zm, r, z);
            
            // Compute beta = (r_new^T * z_new) / (r_old^T * z_old)
            Scalar rz_new = dot_product(r, z);
            Scalar beta = rz_new / rz_old;
            
            // Update search direction: p = z + beta * p
            update_vector(p, z, 1.0, beta);
            
            rz_old = rz_new;
        }
    }
    
    // Get solver statistics
    Ordinal get_iterations() const { return iterations; }
    Scalar get_final_residual() const { return final_residual; }
    bool is_converged() const { return converged; }
    
private:
    // Build preconditioner (diagonal inverse for Jacobi/SSOR)
    void build_preconditioner_2d(const ActiveCellMesh& mesh,
                                const View1D<Scalar>& diag,
                                const View1D<Scalar>& xp,
                                const View1D<Scalar>& xm,
                                const View1D<Scalar>& yp,
                                const View1D<Scalar>& ym) {
        
        if (precond_type == PreconditionerType::NONE) {
            Kokkos::deep_copy(M_inv, 1.0);
            return;
        }
        
        // For Jacobi: M = diag, M_inv = 1/diag
        // For SSOR: M = (D/omega) * (D + omega*L) * D^{-1} * (D + omega*U)
        // Simplified: use diagonal approximation M ≈ diag
        auto _diag = diag;
        auto _M_inv = M_inv;
        
        Kokkos::parallel_for(RangePolicy(0, num_active), 
            KOKKOS_LAMBDA (const Ordinal i) {
                if (std::abs(_diag(i)) > 1.0e-20) {
                    _M_inv(i) = 1.0 / _diag(i);
                } else {
                    _M_inv(i) = 1.0;  // Avoid division by zero
                }
            });
    }
    
    void build_preconditioner_3d(const ActiveCellMesh& mesh,
                                 const View1D<Scalar>& diag,
                                 const View1D<Scalar>& xp,
                                 const View1D<Scalar>& xm,
                                 const View1D<Scalar>& yp,
                                 const View1D<Scalar>& ym,
                                 const View1D<Scalar>& zp,
                                 const View1D<Scalar>& zm) {
        
        if (precond_type == PreconditionerType::NONE) {
            Kokkos::deep_copy(M_inv, 1.0);
            return;
        }
        
        auto _diag = diag;
        auto _M_inv = M_inv;
        auto _active_to_domain = mesh.active_to_domain;
        
        // Note: diag is indexed by domain index, M_inv is indexed by active index
        Kokkos::parallel_for(RangePolicy(0, num_active),
            KOKKOS_LAMBDA (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (Kokkos::abs(_diag(domain_idx)) > 1.0e-20) {
                    _M_inv(i) = 1.0 / _diag(domain_idx);
                } else {
                    _M_inv(i) = 1.0;
                }
            });
    }
    
    // Apply preconditioner: z = M^{-1} * r
    void apply_preconditioner_2d(const ActiveCellMesh& mesh,
                                 const View1D<Scalar>& diag,
                                 const View1D<Scalar>& xp,
                                 const View1D<Scalar>& xm,
                                 const View1D<Scalar>& yp,
                                 const View1D<Scalar>& ym,
                                 const View1D<Scalar>& r,
                                 View1D<Scalar>& z) {
        
        if (precond_type == PreconditionerType::NONE) {
            Kokkos::deep_copy(z, r);
            return;
        }
        
        if (precond_type == PreconditionerType::JACOBI) {
            // Simple diagonal scaling: z = M_inv * r
            auto _M_inv = M_inv;
            auto _r = r;
            auto _z = z;
            
            Kokkos::parallel_for(RangePolicy(0, num_active),
                KOKKOS_LAMBDA (const Ordinal i) {
                    _z(i) = _M_inv(i) * _r(i);
                });
        } else if (precond_type == PreconditionerType::SSOR) {
            // SSOR: Forward and backward sweeps
            // Forward: (D/omega + L) * z1 = r
            // Backward: (D/omega + U) * z = (D/omega) * z1
            // Simplified: use diagonal approximation
            apply_preconditioner_2d(mesh, diag, xp, xm, yp, ym, r, z);
            // For now, use Jacobi as approximation
            auto _M_inv = M_inv;
            auto _r = r;
            auto _z = z;
            
            Kokkos::parallel_for(RangePolicy(0, num_active),
                KOKKOS_LAMBDA (const Ordinal i) {
                    _z(i) = _M_inv(i) * _r(i);
                });
        }
    }
    
    void apply_preconditioner_3d(const ActiveCellMesh& mesh,
                                 const View1D<Scalar>& diag,
                                 const View1D<Scalar>& xp,
                                 const View1D<Scalar>& xm,
                                 const View1D<Scalar>& yp,
                                 const View1D<Scalar>& ym,
                                 const View1D<Scalar>& zp,
                                 const View1D<Scalar>& zm,
                                 const View1D<Scalar>& r,
                                 View1D<Scalar>& z) {
        
        if (precond_type == PreconditionerType::NONE) {
            Kokkos::deep_copy(z, r);
            return;
        }
        
        if (precond_type == PreconditionerType::JACOBI) {
            auto _M_inv = M_inv;
            auto _r = r;
            auto _z = z;
            
            Kokkos::parallel_for(RangePolicy(0, num_active),
                KOKKOS_LAMBDA (const Ordinal i) {
                    _z(i) = _M_inv(i) * _r(i);
                });
        } else if (precond_type == PreconditionerType::SSOR) {
            // Simplified SSOR (use Jacobi for now)
            auto _M_inv = M_inv;
            auto _r = r;
            auto _z = z;
            
            Kokkos::parallel_for(RangePolicy(0, num_active),
                KOKKOS_LAMBDA (const Ordinal i) {
                    _z(i) = _M_inv(i) * _r(i);
                });
        }
    }
    
    // Matrix-vector product: Ap = A * p (2D)
    void matrix_vector_product_2d(const ActiveCellMesh& mesh,
                                 const View1D<Scalar>& diag,
                                 const View1D<Scalar>& xp,
                                 const View1D<Scalar>& xm,
                                 const View1D<Scalar>& yp,
                                 const View1D<Scalar>& ym,
                                 const View1D<Scalar>& p,
                                 View1D<Scalar>& Ap) {
        
        auto _diag = diag;
        auto _xp = xp;
        auto _xm = xm;
        auto _yp = yp;
        auto _ym = ym;
        auto _p = p;
        auto _Ap = Ap;
        
        auto _neighbor_left = mesh.neighbor_left;
        auto _neighbor_right = mesh.neighbor_right;
        auto _neighbor_back = mesh.neighbor_back;
        auto _neighbor_front = mesh.neighbor_front;
        
        // Matrix-vector product for SPD matrix (2D):
        // A*p = diag*p + xp*p[right] + xm*p[left] + yp*p[front] + ym*p[back]
        Kokkos::parallel_for(RangePolicy(0, num_active),
            KOKKOS_LAMBDA (const Ordinal i) {
                // Diagonal contribution
                Scalar result = _diag(i) * _p(i);
                
                // X+ neighbor (right)
                Ordinal n_right = _neighbor_right(i);
                if (n_right >= 0 && n_right < num_active) {
                    result += _xp(i) * _p(n_right);
                }
                
                // X- neighbor (left)
                Ordinal n_left = _neighbor_left(i);
                if (n_left >= 0 && n_left < num_active) {
                    result += _xm(i) * _p(n_left);
                }
                
                // Y+ neighbor (front)
                Ordinal n_front = _neighbor_front(i);
                if (n_front >= 0 && n_front < num_active) {
                    result += _yp(i) * _p(n_front);
                }
                
                // Y- neighbor (back)
                Ordinal n_back = _neighbor_back(i);
                if (n_back >= 0 && n_back < num_active) {
                    result += _ym(i) * _p(n_back);
                }
                
                _Ap(i) = result;
            });
    }
    
    // Matrix-vector product: Ap = A * p (3D)
    // Note: Matrix coefficients (diag, xp, etc.) are indexed by DOMAIN index
    //       Solution vectors (p, Ap) are indexed by ACTIVE index
    void matrix_vector_product_3d(const ActiveCellMesh& mesh,
                                 const View1D<Scalar>& diag,
                                 const View1D<Scalar>& xp,
                                 const View1D<Scalar>& xm,
                                 const View1D<Scalar>& yp,
                                 const View1D<Scalar>& ym,
                                 const View1D<Scalar>& zp,
                                 const View1D<Scalar>& zm,
                                 const View1D<Scalar>& p,
                                 View1D<Scalar>& Ap) {
        
        auto _diag = diag;
        auto _xp = xp;
        auto _xm = xm;
        auto _yp = yp;
        auto _ym = ym;
        auto _zp = zp;
        auto _zm = zm;
        auto _p = p;
        auto _Ap = Ap;
        
        auto _neighbor_left = mesh.neighbor_left;
        auto _neighbor_right = mesh.neighbor_right;
        auto _neighbor_back = mesh.neighbor_back;
        auto _neighbor_front = mesh.neighbor_front;
        auto _neighbor_bottom = mesh.neighbor_bottom;
        auto _neighbor_top = mesh.neighbor_top;
        auto _active_to_domain = mesh.active_to_domain;
        
        // Matrix-vector product for SPD matrix:
        // A*p = diag*p + xp*p[right] + xm*p[left] + yp*p[front] + ym*p[back] + zp*p[top] + zm*p[bottom]
        // Off-diagonal coefficients (xp, xm, etc.) are NEGATIVE for SPD matrix
        Kokkos::parallel_for(RangePolicy(0, num_active),
            KOKKOS_LAMBDA (const Ordinal i) {
                // Get domain index for matrix coefficient access
                Ordinal domain_idx = _active_to_domain(i);
                
                // Diagonal contribution
                Scalar result = _diag(domain_idx) * _p(i);
                
                // X+ neighbor (off-diag is stored as negative value)
                Ordinal n_right = _neighbor_right(i);
                if (n_right >= 0 && n_right < num_active) {
                    result += _xp(domain_idx) * _p(n_right);
                }
                
                // X- neighbor
                Ordinal n_left = _neighbor_left(i);
                if (n_left >= 0 && n_left < num_active) {
                    result += _xm(domain_idx) * _p(n_left);
                }
                
                // Y+ neighbor
                Ordinal n_front = _neighbor_front(i);
                if (n_front >= 0 && n_front < num_active) {
                    result += _yp(domain_idx) * _p(n_front);
                }
                
                // Y- neighbor
                Ordinal n_back = _neighbor_back(i);
                if (n_back >= 0 && n_back < num_active) {
                    result += _ym(domain_idx) * _p(n_back);
                }
                
                // Z+ neighbor (top/shallow)
                Ordinal n_top = _neighbor_top(i);
                if (n_top >= 0 && n_top < num_active) {
                    result += _zp(domain_idx) * _p(n_top);
                }
                
                // Z- neighbor (bottom/deep)
                Ordinal n_bottom = _neighbor_bottom(i);
                if (n_bottom >= 0 && n_bottom < num_active) {
                    result += _zm(domain_idx) * _p(n_bottom);
                }
                
                _Ap(i) = result;
            });
    }
    
    // Vector operations
    Scalar compute_norm(const View1D<Scalar>& v) {
        Scalar norm_sq;
        Kokkos::parallel_reduce(RangePolicy(0, num_active),
            KOKKOS_LAMBDA (const Ordinal i, Scalar& sum) {
                sum += v(i) * v(i);
            }, norm_sq);
        return std::sqrt(norm_sq);
    }
    
    Scalar dot_product(const View1D<Scalar>& a, const View1D<Scalar>& b) {
        Scalar dot;
        auto _a = a;
        auto _b = b;
        Kokkos::parallel_reduce(RangePolicy(0, num_active),
            KOKKOS_LAMBDA (const Ordinal i, Scalar& sum) {
                sum += _a(i) * _b(i);
            }, dot);
        return dot;
    }
    
    // x = x + alpha*y (default), or x = beta*x + alpha*y (if beta specified)
    void update_vector(View1D<Scalar>& x, const View1D<Scalar>& y, Scalar alpha, Scalar beta = 1.0) {
        auto _x = x;
        auto _y = y;
        Kokkos::parallel_for(RangePolicy(0, num_active),
            KOKKOS_LAMBDA (const Ordinal i) {
                _x(i) = beta * _x(i) + alpha * _y(i);
            });
    }
};

// Type aliases for convenience
using PCGSolver2D = PCGSolver<2>;
using PCGSolver3D = PCGSolver<3>;

#endif // FREHG_PCG_SOLVER_HPP

