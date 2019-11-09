#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

// -----------------------------------------------------------------------------
// Functions to read data from input files and to write model results into
// output files.
// - Zhi Li 2017-05-04 -
// -----------------------------------------------------------------------------

#include "fileio.h"
#include "configuration.h"
#include "groundwater.h"
#include "initialize.h"
#include "mpifunctions.h"
#include "scalar.h"

void ReadFile(double *arr, char filename[], int N);
void ReadFileInt(int *arr, char filename[], int N);
void DataOutput(Data **data, Ground **ground, Bath *bath, Config *setting, int tt, int root, int irank);
void SubOutput(Sub **sub, Config *setting, int tt, int root, int irank);
void WriteFile(double *ally, char *filename, Config *setting, int tt, int N);
void combineSerial(double *ally, double *y, Config *setting, int N);
void writeRestartFile(Data **data, Bath *bath, Config *setting, int tt, int root, int irank);

// ==================== read data from input files ====================
void ReadFile(double *arr, char filename[], int N)
{
  int ii = 0, jj;
  FILE *fid;
  fid = fopen(filename, "r");
  if (fid == NULL)
  {
    printf("WARNING: Unable to open the bathymetry file! \n");
  }
  for (ii = 0; ii < N; ii++)
  {
    fscanf(fid, "%lf,", &arr[ii]);
  }
  fclose(fid);
}

// ==================== read data from input files (int) ====================
void ReadFileInt(int *arr, char filename[], int N)
{
  int ii = 0, jj;
  FILE *fid;
  fid = fopen(filename, "r");
  if (fid == NULL)
  {
    printf("WARNING: Unable to open the bathymetry file! \n");
  }
  for (ii = 0; ii < N; ii++)
  {
    fscanf(fid, "%d,", &arr[ii]);
  }
  fclose(fid);
}

// =============== Top level code for data output ===============
void DataOutput(Data **data, Ground **ground, Bath *bath, Config *setting, int tt, int root, int irank)
{
  int ii;
  for (ii = 0; ii < setting->N2ci; ii++)
  {(*data)->surf[ii] = (*data)->surf[ii] - bath->offset[0];}
  if (setting->useMPI == 1)
  {
    combineAllRanks((*data)->alluuXP, (*data)->uuXP, setting, root);
    combineAllRanks((*data)->allvvYP, (*data)->vvYP, setting, root);
    combineAllRanks((*data)->allsurf, (*data)->surf, setting, root);
    combineAllRanks((*data)->alldepth, (*data)->depth, setting, root);
    //if (setting->CDnotN == 0)
    //{
      combineAllRanks((*data)->allCDXP, (*data)->CDXP, setting, root);
      combineAllRanks((*data)->allCDYP, (*data)->CDYP, setting, root);
    //}
    if (setting->useScalar == 1)
    {combineAllRanks((*data)->allS, (*data)->S, setting, root);}
    if (setting->checkConservation == 1)
    {combineAllRanks((*data)->allVloss, (*data)->Vloss, setting, root);}
	if (setting->useSubsurface == 1)
	{combineAllRanksGround((*ground)->allh, (*ground)->h, setting, root);}
    if (setting->useSubscalar == 1)
	{combineAllRanksGround((*ground)->allS, (*ground)->S, setting, root);}
  }
  else
  {
    combineSerial((*data)->alluuXP, (*data)->uuXP, setting, setting->N2ci);
    combineSerial((*data)->allvvYP, (*data)->vvYP, setting, setting->N2ci);
    combineSerial((*data)->allsurf, (*data)->surf, setting, setting->N2ci);
    combineSerial((*data)->alldepth, (*data)->depth, setting, setting->N2ci);
    //if (setting->CDnotN == 0)
    //{
      combineSerial((*data)->allCDXP, (*data)->CDXP, setting, setting->N2ci);
      combineSerial((*data)->allCDYP, (*data)->CDYP, setting, setting->N2ci);
    //}
    if (setting->useScalar == 1)
    {combineSerial((*data)->allS, (*data)->S, setting, setting->N2ci);}
    if (setting->checkConservation == 1)
    {combineSerial((*data)->allVloss, (*data)->Vloss, setting, setting->N2ci);}
    if (setting->useSubsurface == 1)
	{
    combineSerial((*ground)->allh, (*ground)->h, setting, setting->N3ci);
    combineSerial((*ground)->allSw, (*ground)->Sw, setting, setting->N3ci);
  }
    if (setting->useSubscalar == 1)
	{combineSerial((*ground)->allS, (*ground)->S, setting, setting->N3ci);}
  }
  if (irank == 0)
  {
    if (setting->saveuu == 1)
    {WriteFile((*data)->alluuXP, "uVelocity", setting, tt, setting->N2CI);}
    if (setting->savevv == 1)
    {WriteFile((*data)->allvvYP, "vVelocity", setting, tt, setting->N2CI);}
    if (setting->savesurface == 1)
    {WriteFile((*data)->allsurf, "surfaceZ", setting, tt, setting->N2CI);}
    if (setting->savedepth == 1)
    {WriteFile((*data)->alldepth, "depth", setting, tt, setting->N2CI);}
    //if (setting->CDnotN == 0 & setting->saveCD == 1)
        if (setting->saveCD == 1)
    {
      WriteFile((*data)->allCDXP, "cdxp", setting, tt, setting->N2CI);
      WriteFile((*data)->allCDYP, "cdyp", setting, tt, setting->N2CI);
    }
    if (setting->useScalar == 1 & setting->savescalar == 1)
    {WriteFile((*data)->allS, "scalar", setting, tt, setting->N2CI);}
    if (setting->checkConservation == 1)
    {WriteFile((*data)->allVloss, "vloss", setting, tt, setting->N2CI);}
    // groundwater
    if (setting->useSubsurface == 1)
    {
        WriteFile((*ground)->allh, "head", setting, tt, setting->N3CI);
        WriteFile((*ground)->allSw, "satu", setting, tt, setting->N3CI);
        //WriteFile((*ground)->Quu, "Quu", setting, tt, setting->N3ci);
        //WriteFile((*ground)->Qvv, "Qvv", setting, tt, setting->N3ci);
        //WriteFile((*ground)->Qww, "Qww", setting, tt, setting->N3ci);
        WriteFile((*ground)->nlay, "nlay", setting, tt, setting->N2ci);
		if (setting->useSubscalar == 1)
		{WriteFile((*ground)->allS, "subscalar", setting, tt, setting->N3CI);}
    }
  }
  for (ii = 0; ii < setting->N2ci; ii++)
  {(*data)->surf[ii] = (*data)->surf[ii] + bath->offset[0];}
}

// =============== Output the subgrid variables ===============
void SubOutput(Sub **sub, Config *setting, int tt, int root, int irank)
{
    int ii;
    if (setting->useMPI == 1)
    {
        combineAllRanks((*sub)->allNx, (*sub)->Nx, setting, root);
        combineAllRanks((*sub)->allOy, (*sub)->Oy, setting, root);
        combineAllRanks((*sub)->allVx, (*sub)->Vx, setting, root);
        combineAllRanks((*sub)->allVy, (*sub)->Vy, setting, root);
        combineAllRanks((*sub)->allVo, (*sub)->V, setting, root);
        combineAllRanks((*sub)->allZo, (*sub)->Z, setting, root);
    }
    else
    {
        combineSerial((*sub)->allNx, (*sub)->Nx, setting, setting->N2ci);
        combineSerial((*sub)->allOy, (*sub)->Oy, setting, setting->N2ci);
        combineSerial((*sub)->allVx, (*sub)->Vx, setting, setting->N2ci);
        combineSerial((*sub)->allVy, (*sub)->Vy, setting, setting->N2ci);
        combineSerial((*sub)->allVo, (*sub)->V, setting, setting->N2ci);
        combineSerial((*sub)->allZo, (*sub)->Z, setting, setting->N2ci);
    }
    if (irank == 0 & setting->savesub == 1)
    {
        WriteFile((*sub)->allNx, "SubNx", setting, tt, setting->N2CI);
        WriteFile((*sub)->allOy, "SubOy", setting, tt, setting->N2CI);
        WriteFile((*sub)->allVx, "SubVx", setting, tt, setting->N2CI);
        WriteFile((*sub)->allVy, "SubVy", setting, tt, setting->N2CI);
        WriteFile((*sub)->allVo, "SubV", setting, tt, setting->N2CI);
        WriteFile((*sub)->allZo, "SubZ", setting, tt, setting->N2CI);
    }
}

// =============== Write output file ===============
void WriteFile(double *ally, char *filename, Config *setting, int tt, int N)
{
  int ii;
  FILE *fp;
  // convert time step into string
  char tstr[10];
  sprintf(tstr, "%d", tt);
  // create filename for saving
  char *fullname = malloc(100);
  strcpy(fullname, setting->saveFolder);
  //strcat(fullname, "/");
  strcat(fullname, filename);
  strcat(fullname, "_");
  strcat(fullname, tstr);
  strcat(fullname, ".dat");
  // open and write to file
  fp = fopen(fullname, "w");
  for (ii = 0; ii < N; ii++)
  {
    fprintf(fp, "%6.6f \n", ally[ii]);
  }
  fclose(fp);
  free(fullname);
}

// =============== If serial, copy fields to all fields ===============
void combineSerial(double *ally, double *y, Config *setting, int N)
{
  int ii;
  for (ii = 0; ii < N; ii++)
  {ally[ii] = y[ii];}
}

// =============== Write the restart file ===============
void writeRestartFile(Data **data, Bath *bath, Config *setting, int tt, int root, int irank)
{
  int ii, jj = 0, N, n = 3;
  double *allRestart;
  for (ii = 0; ii < setting->N2ci; ii++)
  {(*data)->surf[ii] = (*data)->surf[ii] - bath->offset[0];}
  // combine variables to rank 0
  if (setting->useMPI == 1)
  {
    combineAllRanks((*data)->alluuXP, (*data)->uuXP, setting, root);
    combineAllRanks((*data)->allvvYP, (*data)->vvYP, setting, root);
    combineAllRanks((*data)->allsurf, (*data)->surf, setting, root);
    if (setting->useScalar == 1)
    {combineAllRanks((*data)->allS, (*data)->S, setting, root);}
  }
  else
  {
    combineSerial((*data)->alluuXP, (*data)->uuXP, setting, setting->N2ci);
    combineSerial((*data)->allvvYP, (*data)->vvYP, setting, setting->N2ci);
    combineSerial((*data)->allsurf, (*data)->surf, setting, setting->N2ci);
    if (setting->useScalar == 1)
    {combineSerial((*data)->allS, (*data)->S, setting, setting->N2ci);}
  }
  if (irank == root)
  {
    // initilize the restart file structure
    // [tNStart, surf, uu, vv, scalar]
    N = setting->N2CI*4 + n;
    allRestart = malloc(N*sizeof(double));
    // compute the new starting time
    allRestart[0] = setting->tNEnd;
    allRestart[1] = setting->dx;
    allRestart[2] = setting->dy;
    // copy all other variables
    for (ii = n; ii < setting->N2CI+n; ii++)
    {allRestart[ii] = (*data)->allsurf[jj]; jj++;}
    jj = 0;
    for (ii = setting->N2CI+n; ii < 2*setting->N2CI + n; ii++)
    {allRestart[ii] = (*data)->alluuXP[jj]; jj++;}
    jj = 0;
    for (ii = 2*setting->N2CI + n; ii < 3*setting->N2CI+n; ii++)
    {allRestart[ii] = (*data)->allvvYP[jj]; jj++;}
    jj = 0;
    for (ii = 3*setting->N2CI + n; ii < 4*setting->N2CI+n; ii++)
    {allRestart[ii] = (*data)->allS[jj]; jj++;}
    // write to the restart file
    WriteFile(allRestart, "restart", setting, tt, N);
  }
  for (ii = 0; ii < setting->N2ci; ii++)
  {(*data)->surf[ii] = (*data)->surf[ii] + bath->offset[0];}
}
