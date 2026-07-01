# PETSc solver backend (P2.5.4 + P2.6). This is the ONLY library that links PETSc; the
# decomposition and physics modules never do. PETSc Mat/Vec/KSP types are permitted here.
add_library(frehg2_petsc_backend STATIC
  PetscLinearSolver.cpp
  KokkosPetscBridge.cpp
  KokkosPetscBridgeDevice.cpp   # Option B (GPU-native COO); empty TU on non-GPU builds (P10)
  PetscSubcommSplit.cpp         # GPU async solver-group comm split; empty on non-GPU (P10.3.7)
  LinearSolverFactory.cpp)
add_library(Frehg2::petsc_backend ALIAS frehg2_petsc_backend)

target_link_libraries(frehg2_petsc_backend PUBLIC Frehg2::decomp PETSc::PETSc)
apply_strict_warnings(frehg2_petsc_backend)
