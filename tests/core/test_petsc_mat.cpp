// P1.6 acceptance: verify the local PETSc install supports the design's matrix path.
//
// The Frehg2 solver always builds MATMPIAIJ (serial == one-rank MPIAIJ; no
// MatCreateSeqAIJ, per .cursorrules). PETSc's *default* "aij" type resolves to seqaij on
// a one-rank communicator, which is exactly why the solver selects MATMPIAIJ explicitly.
// This test (1) reports the default type for information and (2) asserts that an
// explicitly-typed MATMPIAIJ matrix builds and reports type "mpiaij".
#include <petscmat.h>

#include <cstdio>

int main(int argc, char** argv) {
  if (PetscInitialize(&argc, &argv, nullptr, nullptr) != 0) {
    std::fprintf(stderr, "PetscInitialize failed\n");
    return 1;
  }

  PetscMPIInt comm_size = 1;
  MPI_Comm_size(PETSC_COMM_WORLD, &comm_size);

  // (1) Informational: PETSc default type via MatSetFromOptions.
  {
    Mat d;
    MatCreate(PETSC_COMM_WORLD, &d);
    MatSetSizes(d, PETSC_DECIDE, PETSC_DECIDE, 4, 4);
    MatSetFromOptions(d);
    MatSetUp(d);
    MatType default_type;
    MatGetType(d, &default_type);
    PetscPrintf(PETSC_COMM_WORLD, "Default mat type (comm size %d): %s\n",
                static_cast<int>(comm_size), default_type);
    MatDestroy(&d);
  }

  // (2) Design path: explicitly MATMPIAIJ, preallocated from (here trivial) ownership.
  Mat A;
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, 4, 4);
  MatSetType(A, MATMPIAIJ);
  MatMPIAIJSetPreallocation(A, 1, nullptr, 0, nullptr);
  MatSetUp(A);

  PetscInt rstart = 0;
  PetscInt rend = 0;
  MatGetOwnershipRange(A, &rstart, &rend);
  for (PetscInt r = rstart; r < rend; ++r) {
    PetscScalar v = 2.0;
    MatSetValues(A, 1, &r, 1, &r, &v, INSERT_VALUES);
  }
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  MatType solver_type;
  MatGetType(A, &solver_type);
  PetscPrintf(PETSC_COMM_WORLD, "Solver mat type: %s\n", solver_type);

  PetscBool is_mpiaij = PETSC_FALSE;
  PetscStrcmp(solver_type, MATMPIAIJ, &is_mpiaij);

  MatDestroy(&A);
  PetscFinalize();

  if (is_mpiaij != PETSC_TRUE) {
    std::fprintf(stderr, "FAIL: expected mpiaij, got a different type\n");
    return 1;
  }
  return 0;
}
