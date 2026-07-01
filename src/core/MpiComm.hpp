// Model-owned MPI communication primitives (P2.3).
//
// MpiComm owns the structured 2D process grid (mpi_nx x mpi_ny) over the global domain
// (nx x ny) and provides block decomposition queries, hand-written MPI_Sendrecv halo
// exchange for halo-padded fields, and a gather-to-rank-0. This is the model's ONLY
// halo-exchange mechanism; DomainDecomposition (P2.4) builds on it. There is no PETSc
// DMDA scatter. The communicator is MPI_COMM_WORLD (== PETSC_COMM_WORLD by default), so
// this header carries NO PETSc dependency.
#ifndef FREHG2_CORE_MPI_COMM_HPP
#define FREHG2_CORE_MPI_COMM_HPP

#include <mpi.h>

#include <tuple>

#include "frehg2/core/define.hpp"

namespace frehg2 {

// MPI message tags for the four halo directions.
struct HaloTags {
  int left = 1001;
  int right = 1002;
  int down = 1003;
  int up = 1004;
};

class MpiComm {
 public:
  // Build a process grid over the global domain. Requires mpi_nx*mpi_ny == comm size.
  MpiComm(int global_nx, int global_ny, int mpi_nx, int mpi_ny,
          MPI_Comm comm = MPI_COMM_WORLD);

  int rank() const { return rank_; }
  int size() const { return size_; }
  MPI_Comm comm() const { return comm_; }
  int mpiNx() const { return mpi_nx_; }
  int mpiNy() const { return mpi_ny_; }
  int globalNx() const { return global_nx_; }
  int globalNy() const { return global_ny_; }

  // This rank's position in the process grid.
  int px() const { return rank_ % mpi_nx_; }
  int py() const { return rank_ / mpi_nx_; }

  // This rank's owned interior block.
  int localNx() const { return partSize(px(), global_nx_, mpi_nx_); }
  int localNy() const { return partSize(py(), global_ny_, mpi_ny_); }
  int i0() const { return partStart(px(), global_nx_, mpi_nx_); }
  int j0() const { return partStart(py(), global_ny_, mpi_ny_); }

  // Block-distribution helpers (also used by DomainDecomposition). part index in [0,P).
  static int partSize(int part, int n, int P);
  static int partStart(int part, int n, int P);
  static int partOf(int global_index, int n, int P);

  // Global (i,j) -> (local_i, local_j, owner_rank).
  std::tuple<int, int, int> globalToLocal(int gi, int gj) const;
  // This rank's (local_i, local_j) -> (global_i, global_j).
  std::tuple<int, int> localToGlobal(int li, int lj) const;

  // Neighbor ranks (MPI_PROC_NULL at domain edges).
  int leftRank() const { return px() > 0 ? rank_ - 1 : MPI_PROC_NULL; }
  int rightRank() const { return px() < mpi_nx_ - 1 ? rank_ + 1 : MPI_PROC_NULL; }
  int downRank() const { return py() > 0 ? rank_ - mpi_nx_ : MPI_PROC_NULL; }
  int upRank() const { return py() < mpi_ny_ - 1 ? rank_ + mpi_nx_ : MPI_PROC_NULL; }

  // Fill the one-cell ghost ring of a halo-padded 2D field. Storage is column-major-free
  // flat: index = c + r*nx_with_halo, with interior c in [1,localNx], r in [1,localNy].
  // nx_with_halo must equal localNx()+2 and ny_with_halo localNy()+2.
  void haloExchange2D(const RealArr1D& field, int nx_with_halo, int ny_with_halo,
                      HaloTags tags = HaloTags()) const;

  // Same, applied independently to each of nz layers. The vertical (z) direction is
  // on-rank only: z-adjacent cells are directly addressable in storage (k +/- 1), so no
  // MPI message is exchanged in z.
  void haloExchange3D(const RealArr1D& field, int nx_with_halo, int ny_with_halo, int nz,
                      HaloTags tags = HaloTags()) const;

  // Gather all ranks' interior surface cells onto rank 0 in global row-major order
  // (gi + gj*global_nx). `local` is the halo-padded local field; `global` is sized
  // global_nx*global_ny and only meaningful on rank 0.
  void gatherToRank0(const RealArr1D& local, int nx_with_halo, int ny_with_halo,
                     const RealArr1DHost& global) const;

 private:
  MPI_Comm comm_;
  int rank_ = 0;
  int size_ = 1;
  int global_nx_ = 1;
  int global_ny_ = 1;
  int mpi_nx_ = 1;
  int mpi_ny_ = 1;
};

}  // namespace frehg2

#endif  // FREHG2_CORE_MPI_COMM_HPP
