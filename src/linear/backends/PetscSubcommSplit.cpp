// P10 Task 10.3.7 — implementation of the GPU async solver-group split. Compiled ONLY on a
// Kokkos GPU backend; an empty object file otherwise (the CPU/OpenMP async path in
// src/driver/AsyncPipeline is the production reference). See PetscSubcommSplit.hpp and
// docs/research_notes/async_gpu.md.
#include "linear/backends/PetscSubcommSplit.hpp"

#if FREHG2_GPU_ASSEMBLY

#include <stdexcept>

namespace frehg2 {

namespace {
void checkSub(PetscErrorCode ierr, const char* what) {
  if (ierr != 0) {
    throw std::runtime_error(std::string("PetscSubcommSplit: ") + what + " failed");
  }
}
}  // namespace

AsyncSolverGroups::AsyncSolverGroups() {
  // Two contiguous color groups over PETSC_COMM_WORLD: {surface, groundwater}.
  checkSub(PetscSubcommCreate(PETSC_COMM_WORLD, &subcomm_), "PetscSubcommCreate");
  checkSub(PetscSubcommSetNumber(subcomm_, 2), "PetscSubcommSetNumber");
  checkSub(PetscSubcommSetType(subcomm_, PETSC_SUBCOMM_CONTIGUOUS), "PetscSubcommSetType");

  // One stream-backed device exec-space instance per domain so the GPU overlaps the two
  // domains' kernels and solves. partition_space hands back independent stream instances.
  auto streams = Kokkos::Experimental::partition_space(DeviceExec(), 1, 1);
  stream_surface_ = streams[0];
  stream_groundwater_ = streams[1];
}

AsyncSolverGroups::~AsyncSolverGroups() { PetscSubcommDestroy(&subcomm_); }

MPI_Comm AsyncSolverGroups::comm(Group g) const {
  // The child communicator is the same handle on every rank; the rank's color (which group it
  // actually participates in) is set by PetscSubcomm. Callers build that domain's solver on it.
  (void)g;
  return PetscSubcommChild(subcomm_);
}

const DeviceExec& AsyncSolverGroups::stream(Group g) const {
  return (g == Group::Surface) ? stream_surface_ : stream_groundwater_;
}

void AsyncSolverGroups::fenceAll() const {
  stream_surface_.fence();
  stream_groundwater_.fence();
}

}  // namespace frehg2

#endif  // FREHG2_GPU_ASSEMBLY
