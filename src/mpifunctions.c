// Functions used for massage passing
#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>

#include "configuration.h"
#include "map.h"

void mpi_bcast_int(int *y, int n, int root);
void mpi_bcast_double(double *y, int n, int root);
void mpi_gather_int(int *y_root, int *y, int n, int root);
void mpi_gather_double(double *y_root, double *y, int n, int root);
void mpi_exchange_surf(double *y, Map *smap, int data_type, Config *param, int irank, int nrank);
void mpi_exchange_subsurf(double *y, Map *gmap, int data_type, Config *param, int irank, int nrank);


// >>>>> MPI Broadcast <<<<<
void mpi_bcast_int(int *y, int n, int root)
{
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&y[0], n, MPI_INT, root, MPI_COMM_WORLD);
}

// >>>>> MPI Broadcast for double <<<<<
void mpi_bcast_double(double *y, int n, int root)
{
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&y[0], n, MPI_DOUBLE, root, MPI_COMM_WORLD);
}

// >>>>> MPI Gather <<<<<
void mpi_gather_int(int *y_root, int *y, int n, int root)
{MPI_Gather(&y[0], n, MPI_INT, &y_root[0], n, MPI_INT, root, MPI_COMM_WORLD);}

// >>>>> MPI Gather for double <<<<<
void mpi_gather_double(double *y_root, double *y, int n, int root)
{MPI_Gather(&y[0], n, MPI_DOUBLE, &y_root[0], n, MPI_DOUBLE, root, MPI_COMM_WORLD);}

// >>>>> MPI exchange for surface domain <<<<<
void mpi_exchange_surf(double *y, Map *smap, int data_type, Config *param, int irank, int nrank)
{
    // get the rank number
    MPI_Status status;
    int idown, iup, ileft, iright;
    if (irank < param->mpi_nx)  {ileft = MPI_PROC_NULL;}
    else    {ileft = irank - param->mpi_nx;}
    if (irank >= param->mpi_nx*(param->mpi_ny-1))   {iright = MPI_PROC_NULL;}
    else    {iright = irank + param->mpi_nx;}
    if (irank % param->mpi_nx == 0) {iup = MPI_PROC_NULL;}
    else    {iup = irank - 1;}
    if ((irank+1) % param->mpi_nx == 0) {idown = MPI_PROC_NULL;}
    else    {idown = irank + 1;}
    // start index of send/recv
    int sendleft = smap->jMin[0], recvleft = smap->jPou[0];
    int sendright = smap->jPin[0], recvright = smap->jMou[0];
    int sendup = smap->iMin[0], recvup = smap->iPou[0];
    int senddown = smap->iPin[0], recvdown = smap->iMou[0];
    // define mpi datatype for left-right exchange
    MPI_Datatype mpi_column;
    if (data_type == 1)
    {MPI_Type_contiguous(param->nx, MPI_INT, &mpi_column);}
    else if (data_type == 2)
    {MPI_Type_contiguous(param->nx, MPI_DOUBLE, &mpi_column);}
    MPI_Type_commit(&mpi_column);
    // define mpi datatype for up-down exchange
    MPI_Datatype mpi_vecsend;
    MPI_Datatype mpi_vecrecv;
    if (data_type == 1)
    {
        MPI_Type_vector(param->ny, 1, param->nx, MPI_INT, &mpi_vecsend);
        MPI_Type_contiguous(param->ny, MPI_INT, &mpi_vecrecv);
    }
    else if (data_type == 2)
    {
        MPI_Type_vector(param->ny, 1, param->nx, MPI_DOUBLE, &mpi_vecsend);
        MPI_Type_contiguous(param->ny, MPI_DOUBLE, &mpi_vecrecv);
    }
    MPI_Type_commit(&mpi_vecsend);
    MPI_Type_commit(&mpi_vecrecv);
    // perform data exchange
    MPI_Sendrecv(&y[sendleft], 1, mpi_column, ileft, 9, &y[recvleft], 1, \
        mpi_column, iright, 9, MPI_COMM_WORLD, &status);
    MPI_Sendrecv(&y[sendright], 1, mpi_column, iright, 9, &y[recvright], 1, \
        mpi_column, ileft, 9, MPI_COMM_WORLD, &status);
    MPI_Sendrecv(&y[sendup], 1, mpi_vecsend, iup, 9, &y[recvup], 1, \
        mpi_vecrecv, idown, 9, MPI_COMM_WORLD, &status);
    MPI_Sendrecv(&y[senddown], 1, mpi_vecsend, idown, 9, &y[recvdown], 1, \
        mpi_vecrecv, iup, 9, MPI_COMM_WORLD, &status);
    MPI_Type_free(&mpi_column);
    MPI_Type_free(&mpi_vecsend);
    MPI_Type_free(&mpi_vecrecv);
}


// >>>>> MPI exchange for subsurface domain <<<<<
void mpi_exchange_subsurf(double *y, Map *gmap, int data_type, Config *param, int irank, int nrank)
{
    // get the rank number
    MPI_Status status;
    int idown, iup, ileft, iright;
    if (irank < param->mpi_nx)  {ileft = MPI_PROC_NULL;}
    else    {ileft = irank - param->mpi_nx;}
    if (irank >= param->mpi_nx*(param->mpi_ny-1))   {iright = MPI_PROC_NULL;}
    else    {iright = irank + param->mpi_nx;}
    if (irank % param->mpi_nx == 0) {iup = MPI_PROC_NULL;}
    else    {iup = irank - 1;}
    if ((irank+1) % param->mpi_nx == 0) {idown = MPI_PROC_NULL;}
    else    {idown = irank + 1;}
    // start index of send/recv
    int sendleft = gmap->jMin[0], recvleft = gmap->jPou[0];
    int sendright = gmap->jPin[0], recvright = gmap->jMou[0];
    int sendup = gmap->iMin[0], recvup = gmap->iPou[0];
    int senddown = gmap->iPin[0], recvdown = gmap->iMou[0];
    // define mpi datatype for left-right exchange
    MPI_Datatype mpi_column;
    if (data_type == 1)
    {MPI_Type_contiguous(param->nx*param->nz, MPI_INT, &mpi_column);}
    else if (data_type == 2)
    {MPI_Type_contiguous(param->nx*param->nz, MPI_DOUBLE, &mpi_column);}
    MPI_Type_commit(&mpi_column);
    // define mpi datatype for up-down exchange
    MPI_Datatype mpi_vecsend;
    MPI_Datatype mpi_vecrecv;
    if (data_type == 1)
    {
        MPI_Type_vector(param->ny, param->nz, param->nx*param->nz, MPI_INT, &mpi_vecsend);
        MPI_Type_contiguous(param->ny*param->nz, MPI_INT, &mpi_vecrecv);
    }
    else if (data_type == 2)
    {
        MPI_Type_vector(param->ny, param->nz, param->nx*param->nz, MPI_DOUBLE, &mpi_vecsend);
        MPI_Type_contiguous(param->ny*param->nz, MPI_DOUBLE, &mpi_vecrecv);
    }
    MPI_Type_commit(&mpi_vecsend);
    MPI_Type_commit(&mpi_vecrecv);
    // perform data exchange
    MPI_Sendrecv(&y[sendleft], 1, mpi_column, ileft, 9, &y[recvleft], 1, \
        mpi_column, iright, 9, MPI_COMM_WORLD, &status);
    MPI_Sendrecv(&y[sendright], 1, mpi_column, iright, 9, &y[recvright], 1, \
        mpi_column, ileft, 9, MPI_COMM_WORLD, &status);
    MPI_Sendrecv(&y[sendup], 1, mpi_vecsend, iup, 9, &y[recvup], 1, \
        mpi_vecrecv, idown, 9, MPI_COMM_WORLD, &status);
    MPI_Sendrecv(&y[senddown], 1, mpi_vecsend, idown, 9, &y[recvdown], 1, \
        mpi_vecrecv, iup, 9, MPI_COMM_WORLD, &status);
    MPI_Type_free(&mpi_column);
    MPI_Type_free(&mpi_vecsend);
    MPI_Type_free(&mpi_vecrecv);
}
