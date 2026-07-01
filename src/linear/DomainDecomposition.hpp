// Model-owned, backend-agnostic structured parallel layout (P2.4).
//
// DomainDecomposition replaces the PETSc DMDA. It describes which global cells each rank
// owns, the canonical global cell numbering used for matrix rows (contiguous within a
// rank, ranks ordered by the process grid), the STAR stencil nonzero structure, and the
// ghost exchange (built on MpiComm). It contains NO PETSc/Trilinos/solver-library types,
// so it links against Kokkos + MPI only (verified by tests/linear/test_decomp_no_petsc).
//
// Decomp2D is used by the SWE solver (P4); Decomp3D by the RE solver (P5). The vertical
// (z) direction is on-rank only (no MPI in z), matching MpiComm::haloExchange3D.
#ifndef FREHG2_LINEAR_DOMAIN_DECOMPOSITION_HPP
#define FREHG2_LINEAR_DOMAIN_DECOMPOSITION_HPP

#include <array>
#include <tuple>

#include "core/MpiComm.hpp"
#include "frehg2/core/define.hpp"

namespace frehg2 {

// Maximum STAR-stencil columns (7-point 3D includes the center).
inline constexpr int kMaxStencil = 7;

// Polymorphic base consumed by SparseSystem::preallocate and the assembly kernels. The
// "owned rows" are indexed 0..ownedRowCount()-1 in the same order as the global numbering
// within this rank (so owned local index L maps to global row ownershipRange().first + L).
class DecompBase {
 public:
  explicit DecompBase(const MpiComm& comm) : comm_(comm) {}
  virtual ~DecompBase() = default;

  const MpiComm& comm() const { return comm_; }

  // Total global rows across all ranks.
  virtual int globalRowCount() const = 0;
  // Half-open [row_start, row_end) global-row range owned by this rank.
  virtual std::pair<int, int> ownershipRange() const = 0;
  int ownedRowCount() const {
    auto r = ownershipRange();
    return r.second - r.first;
  }

  // For owned local row L: write its global row index and its STAR stencil columns
  // (global indices, in-domain only). cols must have capacity kMaxStencil.
  virtual void ownedRowStencil(int local_owned, int& global_row, int* cols,
                               int& ncols) const = 0;

  // Storage size of a halo-padded local field for this decomposition.
  virtual int localHaloFieldSize() const = 0;
  // Fill ghost cells of a halo-padded local field via MpiComm.
  virtual void haloExchange(const RealArr1D& field) const = 0;

  // Convert between halo-padded local storage and the contiguous owned vector used for
  // RHS/solution (owned[L] corresponds to global row ownershipRange().first + L).
  virtual void packOwned(const RealArr1D& field_halo, const RealArr1D& owned) const = 0;
  virtual void unpackOwned(const RealArr1D& owned, const RealArr1D& field_halo) const = 0;

 protected:
  const MpiComm& comm_;
};

// 2D structured decomposition (SWE). 5-point STAR.
class Decomp2D : public DecompBase {
 public:
  explicit Decomp2D(const MpiComm& comm);

  std::pair<int, int> globalSize() const { return {nx_, ny_}; }
  // (i0, j0, ni, nj) owned interior range for this rank.
  std::tuple<int, int, int, int> localCorners() const;

  // Global matrix row for ANY in-domain cell (gi,gj), resolving the owner rank.
  int globalRow(int gi, int gj) const;

  // STAR columns at (gi,gj): center + in-domain {W,E,S,N}. Up to 5.
  void stencilColumns(int gi, int gj, int* cols, int& ncols) const;

  int globalRowCount() const override { return nx_ * ny_; }
  std::pair<int, int> ownershipRange() const override;
  void ownedRowStencil(int local_owned, int& global_row, int* cols,
                       int& ncols) const override;
  int localHaloFieldSize() const override {
    return (comm_.localNx() + 2) * (comm_.localNy() + 2);
  }
  void haloExchange(const RealArr1D& field) const override;
  void packOwned(const RealArr1D& field_halo, const RealArr1D& owned) const override;
  void unpackOwned(const RealArr1D& owned, const RealArr1D& field_halo) const override;

 private:
  int rowStartOfRank(int px, int py) const;
  int nx_;
  int ny_;
};

// 3D structured decomposition (RE). 7-point STAR; z on-rank.
class Decomp3D : public DecompBase {
 public:
  Decomp3D(const MpiComm& comm, int nz);

  std::tuple<int, int, int> globalSize() const { return {nx_, ny_, nz_}; }
  std::tuple<int, int, int, int, int, int> localCorners() const;  // i0,j0,k0,ni,nj,nk

  int globalRow(int gi, int gj, int k) const;
  void stencilColumns(int gi, int gj, int k, int* cols, int& ncols) const;

  int globalRowCount() const override { return nx_ * ny_ * nz_; }
  std::pair<int, int> ownershipRange() const override;
  void ownedRowStencil(int local_owned, int& global_row, int* cols,
                       int& ncols) const override;
  int localHaloFieldSize() const override {
    return (comm_.localNx() + 2) * (comm_.localNy() + 2) * (nz_ + 2);
  }
  void haloExchange(const RealArr1D& field) const override;
  void packOwned(const RealArr1D& field_halo, const RealArr1D& owned) const override;
  void unpackOwned(const RealArr1D& owned, const RealArr1D& field_halo) const override;

 private:
  int rowStartOfRank(int px, int py) const;
  int nx_;
  int ny_;
  int nz_;
};

}  // namespace frehg2

#endif  // FREHG2_LINEAR_DOMAIN_DECOMPOSITION_HPP
