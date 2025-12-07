// Head file for mpifunctions.c

#include "configuration.h"
#include "map.h"

void mpi_bcast_int(int *y, int n, int root);
void mpi_bcast_double(double *y, int n, int root);
void mpi_gather_int(int *y_root, int *y, int n, int root);
void mpi_gather_double(double *y_root, double *y, int n, int root);
void mpi_exchange_surf(double *y, Map *smap, int data_type, Config *param, int irank, int nrank);
void mpi_exchange_subsurf(double *y, Map *gmap, int data_type, Config *param, int irank, int nrank);
