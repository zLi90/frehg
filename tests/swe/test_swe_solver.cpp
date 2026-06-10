#include "swe/SweSolver.hpp"

#ifdef USE_PETSC
#include <petsc.h>
#endif

#include <cmath>

int main(int argc, char** argv)
{
#ifdef USE_PETSC
    PetscInitialize(&argc, &argv, nullptr, nullptr);
#else
    (void)argc;
    (void)argv;
#endif

    int result = 0;
    try {
        frehg2::GridSpec spec;
        spec.nx = 1;
        spec.ny = 1;
        spec.nz = 1;

        const frehg2::SweSolver solver(frehg2::Grid(spec), frehg2::SweParameters{});
        frehg2::SweLinearSystem system;
        system.n = 2;
        system.rhs = {1.0, 2.0};
        system.entries = {
            {0, 0, 2.0},
            {1, 1, 4.0},
        };

        const auto solution = solver.solveLinearSystem(system);
        if (solution.size() != 2 || std::abs(solution[0] - 0.5) > 1.0e-10 ||
            std::abs(solution[1] - 0.5) > 1.0e-10) {
            result = 1;
        }
    } catch (...) {
        result = 1;
    }

#ifdef USE_PETSC
    PetscFinalize();
#endif
    return result;
}
