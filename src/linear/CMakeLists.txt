# Backend-agnostic parallel layout (P2.4). Links Kokkos + MPI only (via frehg2_core) and
# carries NO PETSc symbols, so DomainDecomposition is unit-testable without PETSc
# (tests/linear/test_decomp_no_petsc).
add_library(frehg2_decomp STATIC
  DomainDecomposition.cpp)
add_library(Frehg2::decomp ALIAS frehg2_decomp)

target_link_libraries(frehg2_decomp PUBLIC Frehg2::core)
apply_strict_warnings(frehg2_decomp)

# Solver backends (P2.5/P2.6): the PETSc backend is isolated here so physics and the
# decomposition never link PETSc directly.
add_subdirectory(backends)
