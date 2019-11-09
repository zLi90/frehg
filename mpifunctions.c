#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<mpi.h>

// -----------------------------------------------------------------------------
// Functions for MPI use
// - Zhi Li 2017-05-10 -
// -----------------------------------------------------------------------------

#include "configuration.h"
#include "map.h"

void combineAllRanks(double *ally, double *y, Config *setting, int root);
void combineAllRanksGround(double *ally, double *y, Config *setting, int root);
void mpiexchange(double *S, Maps *map, Config *setting, int irank, int nrank);
void mpiexchangeInt(int *S, Maps *map, Config *setting, int irank, int nrank);
void mpiexchangeGround(double *S, Gmaps *gmap, Config *setting, int irank, int nrank);

// =============== Combine all ranks into one ===============
void combineAllRanks(double *ally, double *y, Config *setting, int root)
{
  MPI_Gather(&y[0], setting->N2ci, MPI_DOUBLE, \
    &ally[0], setting->N2ci, MPI_DOUBLE, root, MPI_COMM_WORLD);
}

// =============== Combine all ranks into one ===============
void combineAllRanksGround(double *ally, double *y, Config *setting, int root)
{
  MPI_Gather(&y[0], setting->N3ci, MPI_DOUBLE, \
    &ally[0], setting->N3ci, MPI_DOUBLE, root, MPI_COMM_WORLD);
}

// =============== Message Passing via MPI ===============
void mpiexchange(double *S, Maps *map, Config *setting, int irank, int nrank)
{
  // set send and receive ranks
  MPI_Status status;
  int ileft = irank - 1;
  int iright = irank + 1;
  if (irank == 0) {ileft = MPI_PROC_NULL; }
  if (irank == nrank-1) {iright = MPI_PROC_NULL; }
   // starting indices
  int sendleft = 0;
  int recvleft = map->jPgt[0];
  int sendright = setting->nx*(setting->ny-1);
  int recvright = map->jMgt[0];
   // define column datatype
  MPI_Datatype mpi_column;
  MPI_Type_contiguous(setting->nx, MPI_DOUBLE, &mpi_column);
  MPI_Type_commit(&mpi_column);
  // send left
  MPI_Sendrecv(&S[sendleft], 1, mpi_column, ileft, 9, &S[recvleft], 1, \
    mpi_column, iright, 9, MPI_COMM_WORLD, &status);
  // send right
  MPI_Sendrecv(&S[sendright], 1, mpi_column, iright, 9, &S[recvright], 1, \
    mpi_column, ileft, 9, MPI_COMM_WORLD, &status);
  MPI_Type_free(&mpi_column);
}

// =============== Message Passing via MPI, the int version ===============
void mpiexchangeInt(int *S, Maps *map, Config *setting, int irank, int nrank)
{
  // set send and receive ranks
  MPI_Status status;
  int ileft = irank - 1;
  int iright = irank + 1;
  if (irank == 0) {ileft = MPI_PROC_NULL; }
  if (irank == nrank-1) {iright = MPI_PROC_NULL; }
   // starting indices
  int sendleft = 0;
  int recvleft = map->jPgt[0];
  int sendright = setting->nx*(setting->ny-1);
  int recvright = map->jMgt[0];
   // define column datatype
  MPI_Datatype mpi_column;
  MPI_Type_contiguous(setting->nx, MPI_INT, &mpi_column);
  MPI_Type_commit(&mpi_column);
  // send left
  MPI_Sendrecv(&S[sendleft], 1, mpi_column, ileft, 9, &S[recvleft], 1, \
    mpi_column, iright, 9, MPI_COMM_WORLD, &status);
  // send right
  MPI_Sendrecv(&S[sendright], 1, mpi_column, iright, 9, &S[recvright], 1, \
    mpi_column, ileft, 9, MPI_COMM_WORLD, &status);
  MPI_Type_free(&mpi_column);
}

// =============== Message Passing via MPI, groundwater ===============
void mpiexchangeGround(double *S, Gmaps *gmap, Config *setting, int irank, int nrank)
{
	// set send and receive ranks
	MPI_Status status;
	int ileft = irank - 1;
	int iright = irank + 1;
	if (irank == 0) {ileft = MPI_PROC_NULL; }
	if (irank == nrank-1) {iright = MPI_PROC_NULL; }
	// starting indices
	int sendleft = 0;
	int recvleft = setting->N3ci;
	int sendright = setting->nx * (setting->ny-1) * gmap->maxLay;
	int recvright = setting->N3ci + setting->nx * gmap->maxLay;
	// define column datatype
	MPI_Datatype mpi_column;
	MPI_Type_contiguous(setting->nx * gmap->maxLay, MPI_DOUBLE, &mpi_column);
	MPI_Type_commit(&mpi_column);
	// send left
	MPI_Sendrecv(&S[sendleft], 1, mpi_column, ileft, 9, &S[recvleft], 1, \
	mpi_column, iright, 9, MPI_COMM_WORLD, &status);
	// send right
	MPI_Sendrecv(&S[sendright], 1, mpi_column, iright, 9, &S[recvright], 1, \
	mpi_column, ileft, 9, MPI_COMM_WORLD, &status);
	MPI_Type_free(&mpi_column);
}


