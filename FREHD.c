#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<mpi.h>

// -----------------------------------------------------------------------------
// This is the top level script of the SubFREHD-C model.
// - Zhi Li 2017-05-10 -
// -----------------------------------------------------------------------------

#include "bathymetry.h"
#include "configuration.h"
#include "fileio.h"
#include "groundwater.h"
#include "initialize.h"
#include "map.h"
#include "mpifunctions.h"
#include "nsfunctions.h"
#include "nssolve.h"
#include "subgrid.h"
#include "utilities.h"

// number of config items to be read as doubles
#define ND 120
// number of total lines to be read
#define NL 150

void ReadInputFile(char filename[], double *value, char *date, char *inputFolder, char *outputFolder, char *subgridFolder, int N)
{
    int ii = 0;
    char *ptr, *token, arr[N];
    double ret;
    FILE *fid;
    fid = fopen(filename, "r");
    if (fid == NULL)
    {printf("WARNING: Unable to open the bathymetry file! \n");}
    while (fgets(arr, N, fid) != NULL)
    {
        // Ignore lines start with '#'
        char *line = arr;
        if (*line == '#')
        {continue;}
        // Ignore contents after ';'
        token = strtok(arr, ";");
        if (ii == 0)
        {strcpy(date, token);}
        else if (ii == 1)
        {strcpy(inputFolder, token);}
        else if (ii == 2)
        {strcpy(outputFolder, token);}
        else if (ii == 3)
        {strcpy(subgridFolder, token);}
        else
        {
            ret = strtod(token, &ptr);
            value[ii] = ret;
        }
        ii++;
    }
    fclose(fid);
}

// ==================== Main Function ====================
int main(int argc, char *argv[])
{
  Config *setting;
  Bath *bath;
  Maps *map;
  Gmaps *gmap;
  Data *data;
  IC *ic;
  BC *bc;
  Sub *sub;
  Ground *ground;
  int irank = 0, nrank = 1;
  writeText("===== FrehdC Model Invoked ... =====", irank);
  double ts, te;
  double *value = malloc(ND*sizeof(double));
  char *date = malloc(10*sizeof(char));
  char *inputFolder = malloc(100*sizeof(char));
  char *outputFolder = malloc(100*sizeof(char));
  char *subgridFolder = malloc(100*sizeof(char));
  writeText("----- STATUS ----- Data structure initialized...", irank);
  // timer
  if (irank == 0) {ts = clock();}
  // read user settings
  writeText("----- STATUS ----- Reading user settings...", irank);
  ReadInputFile("data.dat", value, date, inputFolder, outputFolder, subgridFolder, NL);
  ReadUserSettings(&setting, value, date, inputFolder, outputFolder, subgridFolder);
  // start MPI
  if (setting->useMPI == 1)
  {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
  // read bathymetry
  writeText("----- STATUS ----- Reading bathymetry...", irank);
  ReadBathymetry(&bath, setting, irank);
  InitBathymetry(bath, setting, irank);
  writeText("----- STATUS ----- Creating maps...", irank);
  createMaps(&map, setting);
  if (setting->useSubsurface == 1)
  {createGmaps(&gmap, bath, map, setting, 0, irank, nrank);}
  // Initialize the model
  writeText("----- STATUS ----- Initializing...", irank);
  Init(&data, &map, &ground, &gmap, &ic, &bc, &sub, bath, setting, irank, nrank);
  // Solve the Navier Stokes equations
  writeText("----- STATUS ----- Begin to solve...", irank);
  SolveAll(&data, &sub, &ground, map, gmap, bc, bath, setting, irank, nrank);
  // end MPI
  if (setting->useMPI == 1)
  { MPI_Finalize();}
  writeText("----- STATUS ----- Completing...", irank);
  // timer
  if (irank == 0) {te = clock(); printf("Total time = %lf min\n",(te-ts)/60/(float)CLOCKS_PER_SEC);}
  return 0;
}
