#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

// -----------------------------------------------------------------------------
// Functions to create the bathymetry. It reads the bathymetry data from input
// file and convert them to the format can be used by the solver
// - Zhi Li 2017-05-04 -
// -----------------------------------------------------------------------------

#include "configuration.h"
#include "bathymetry.h"
#include "fileio.h"
#include "subgrid.h"
#include "utilities.h"

// ========== Read bathymetry data from file ==========
void ReadBathymetry(Bath **bath, Config *setting, int irank)
{
  // allocate memory for the bathymetry
  int ii = 0, jj, N = 1;
  double Zmin, nextvalue;
  *bath = malloc(sizeof(Bath));
  (*bath)->allZ = malloc(setting->N2CI * sizeof(double));
  (*bath)->allEdgeX = malloc(setting->N2CI * sizeof(double));
  (*bath)->allEdgeY = malloc(setting->N2CI * sizeof(double));
  (*bath)->bottomZ = malloc(setting->N2ct * sizeof(double));
  (*bath)->edgeX = malloc(setting->N2ci * sizeof(double));
  (*bath)->edgeY = malloc(setting->N2ci * sizeof(double));
  (*bath)->bottomXP = malloc(setting->N2ct * sizeof(double));
  (*bath)->bottomYP = malloc(setting->N2ct * sizeof(double));
  (*bath)->offset = malloc(N*sizeof(double));
  // read bathymetry data from file
  char fullname[100];
  strcpy(fullname, setting->inputFolder);
  if (setting->useSubgrid == 0)
  {
    strcat(fullname, "bath.dat");
    ReadFile((*bath)->allZ, fullname, setting->N2CI);
  }
  else
  {
    strcat(fullname, setting->subgridFolder);
    strcat(fullname, "bottom.dat");
    ReadFile((*bath)->allZ, fullname, setting->N2CI);
  }
  if (setting->useCellEdge == 1)
  {
    char edgexname[100], edgeyname[100];
    strcpy(edgexname, setting->inputFolder);
    strcat(edgexname, "edgex.dat");
    strcpy(edgeyname, setting->inputFolder);
    strcat(edgeyname, "edgey.dat");
    ReadFile((*bath)->allEdgeX, edgexname, setting->N2CI);
    ReadFile((*bath)->allEdgeY, edgeyname, setting->N2CI);
  }
  else if (setting->useSubgrid != 0)
  {
    char edgexname[100], edgeyname[100];
    strcpy(edgexname, setting->inputFolder);
    strcat(edgexname, setting->subgridFolder);
    strcat(edgexname, "bottomXP.dat");
    strcpy(edgeyname, setting->inputFolder);
    strcat(edgeyname, setting->subgridFolder);
    strcat(edgeyname, "bottomYP.dat");
    ReadFile((*bath)->allEdgeX, edgexname, setting->N2CI);
    ReadFile((*bath)->allEdgeY, edgeyname, setting->N2CI);
  }
}

// ========== Init bathymetry ==========
void InitBathymetry(Bath *bath, Config *setting, int irank)
{
  int jj, N = 1;
  double Zmin;
  // compute the offset
  Zmin = getMin(bath->allZ, setting->N2CI);
  if (Zmin >= 0)
  {bath->offset[0] = 0;}
  else
  {bath->offset[0] = -Zmin;}
  //assign bathymetry to each subdomain
  for (jj = 0; jj < setting->N2ct; jj++)
  {bath->bottomZ[jj] = 0;}
  for (jj = 0; jj < setting->N2ci; jj++)
  {
    bath->bottomZ[jj] = bath->allZ[jj+irank*setting->N2ci];
    bath->bottomZ[jj] = bath->bottomZ[jj] + bath->offset[0];
  }
  if (setting->useSubsurface == 1)
    {setting->zbot += bath->offset[0];}

  //assign cell edges to each subdomain
  if (setting->useCellEdge == 1 | setting->useSubgrid != 0)
  {
    for (jj = 0; jj < setting->N2ci; jj++)
    {
      bath->edgeX[jj] = bath->allEdgeX[jj+irank*setting->N2ci];
      bath->edgeX[jj] = bath->edgeX[jj] + bath->offset[0];
      bath->edgeY[jj] = bath->allEdgeY[jj+irank*setting->N2ci];
      bath->edgeY[jj] = bath->edgeY[jj] + bath->offset[0];
    }
  }
}
