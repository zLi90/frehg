// ============================================================================
//                           TOP LEVEL CODE OF FREHG
//                              Zhi Li 2020-05-19
// ============================================================================
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<mpi.h>
// #include<omp.h>

#include"configuration.h"
#include"initialize.h"
#include"map.h"
#include"solve.h"
#include"utility.h"

int main(int argc, char *argv[])
{
    Config *param;
    Data *data;
    Map *smap;
    Map *gmap;
    int irank = 0,  nrank = 1;
    float t0, t1, t_total;

    read_input(&param);
    // omp_set_num_threads(param->nthreads);
    if (param->use_mpi == 1)
    {
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &irank);
        MPI_Comm_size(MPI_COMM_WORLD, &nrank);
    }
    if (irank == 0) {t0 = clock();}

    mpi_print("\n\n>>>>>>  Starting FREHG simulation  <<<<<< ",irank);
    mpi_print("   >>>  Author: Zhi Li(LBNL, 2020) <<<  \n\n",irank);
    init(&data, &smap, &gmap, &param, irank, nrank);
    solve(&data, smap, gmap, param, irank, nrank);


    if (param->use_mpi == 1)    {MPI_Finalize();}

    if (irank == 0) {
        t1 = clock();
        t_total = (t1 - t0)/(float)CLOCKS_PER_SEC;
        printf("  >>>> Total time = %.4f sec <<<<  \n", t_total);
    }
    mpi_print(">>>>>>   Ending FREHG simulation  <<<<<< ",irank);

    return 0;
}
