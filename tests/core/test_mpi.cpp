#include "core/MpiComm.hpp"

#ifdef USE_MPI
#include <mpi.h>
#endif

int main(int argc, char** argv)
{
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
#else
    (void)argc;
    (void)argv;
#endif

    int result = 0;

    try {
        frehg2::MpiComm comm(2, 2);

        if (comm.size() != 4) {
            result = 1;
        }
        if (comm.rankX() != comm.rank() % 2) {
            result = 1;
        }
        if (comm.rankY() != comm.rank() / 2) {
            result = 1;
        }

#ifdef USE_MPI
        if (comm.rank() == 0) {
            if (comm.neighbor(frehg2::Direction::XM) != MPI_PROC_NULL ||
                comm.neighbor(frehg2::Direction::YM) != MPI_PROC_NULL ||
                comm.neighbor(frehg2::Direction::XP) != 1 ||
                comm.neighbor(frehg2::Direction::YP) != 2) {
                result = 1;
            }
        }
#endif
    } catch (...) {
        result = 1;
    }

#ifdef USE_MPI
    int global_result = 0;
    MPI_Allreduce(&result, &global_result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Finalize();
    return global_result;
#else
    return result;
#endif
}
