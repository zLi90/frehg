#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<sys/stat.h>

// -----------------------------------------------------------------------------
// Functions to initialize the model, which includes reading the bathymetry,
// creating the grids and maps, reading the initial and boundary conditions and
// set up all other parameters for time stepping
// - Zhi Li 2017-05-04 -
// -----------------------------------------------------------------------------

#include "bathymetry.h"
#include "configuration.h"
#include "fileio.h"
#include "groundwater.h"
#include "initialize.h"
#include "map.h"
#include "nsfunctions.h"
#include "mpifunctions.h"
#include "scalar.h"
#include "subgrid.h"
#include "utilities.h"

void Init(Data **data, Maps **map, Ground **ground, Gmaps **gmap, IC **ic, BC **bc, Sub **sub, Bath *bath, Config *setting, int irank, int nrank);
void enforceBathBC(Bath *bath, Maps *map, Config *setting, int irank, int nrank);
void initDataArrays(Data **data, Config *setting);
void initFieldValues(Data **data, Bath *bath, IC *ic, Config *setting, int irank);
void initICArrays(IC **ic, Config *setting);
void initBCArrays(BC **bc, Config *setting);
void initCD(Data **data, Config *setting);
void readIC(IC **ic, Bath *bath, Config *setting);
void readBC(BC **bc, Config *setting, int irank);
//void readRestartFile(IC **ic, Config **setting);
void readOneBC(double *arr, char filename[], Config *setting, int N);
void initGroundArrays(Ground **ground, Gmaps *gmap, Bath *bath, Config *setting);


// =============== Top level code for initialization ===============
void Init(Data **data, Maps **map, Ground **ground, Gmaps **gmap, IC **ic, BC **bc, Sub **sub, Bath *bath, Config *setting, int irank, int nrank)
{
  enforceBathBC(bath, *map, setting, irank, nrank);
  updateAllFaceElev(bath, *map, setting);
  if (setting->useMPI == 1)
  {
    mpiexchange(bath->bottomZ, *map, setting, irank, nrank);
    mpiexchange(bath->bottomXP, *map, setting, irank, nrank);
    mpiexchange(bath->bottomYP, *map, setting, irank, nrank);
  }
  // allocate memory for the data arrays
  if (irank == 0)
  {printf("Initializing data arrays ...\n");}
  initDataArrays(data, setting);
  // allocate memory for groundwater model
  if (setting->useSubsurface == 1)
  {initGroundArrays(ground, *gmap, bath, setting);}
  // read initial conditions
  if (irank == 0)
  {printf("Setting up initial conditions ...\n");}
  initICArrays(ic, setting);
  readIC(ic, bath, setting);
  // read boundary conditions
  if (irank == 0)
  {printf("Setting up boundary conditions ...\n");}
  initBCArrays(bc, setting);
  readBC(bc, setting, irank);
  // initialize bottom CD
  if (irank == 0)
  {printf("Initializing data arrays ...\n");}
  initCD(data, setting);
  // initialize U, V, surf with IC for internal cells
  initFieldValues(data, bath, *ic, setting, irank);
  // adjust free surface at dry regions
  zeroNegSurf(data, bath, setting, 0);
  // assign BC on ghost cells
  enforceFreeSurfaceBC(data, *map, setting);
  enforceVelocityBC(data, *map, setting);
  if (setting->useMPI == 1)
  {
    mpiexchange((*data)->uuXP, *map, setting, irank, nrank);
    mpiexchange((*data)->vvYP, *map, setting, irank, nrank);
  }
  // enforce tidal BC at tidal boundary
  enforceTidalBC(data, bath, *bc, setting, 0, irank, nrank);
  if (setting->useMPI == 1)
  {
    mpiexchange((*data)->surf, *map, setting, irank, nrank);
  }
  // initialize for cell center and edge depth
  updateDepth(data, *map, bath, setting);
  updateAllFaceDepth(data, bath, *map, setting, irank, nrank);
  if (setting->useMPI == 1)
  {
    mpiexchange((*data)->depthXP, *map, setting, irank, nrank);
    mpiexchange((*data)->depthYP, *map, setting, irank, nrank);
  }
  // initialize subgrid
  if (setting->useSubgrid == 1 || setting->CDnotN == 2)
  {initAllSubVar(sub, bath, setting, irank);}
//  initSubArea(sub, bath, setting, irank);
//  readSubArea(sub, bath, setting, irank);
//  splitSubArea(sub, bath, setting, irank);

  if (setting->useSubgrid != 0 || setting->CDnotN == 2)
  {
    writeText("----- STATUS ----- Subgrid model is invoked...", irank);
    computeSubArea(sub, *data, bath, *map, setting);
    if (setting->useMPI == 1)
    {
      mpiexchangeInt((*sub)->ind, *map, setting, irank, nrank);
      mpiexchange((*sub)->Op, *map, setting, irank, nrank);
      mpiexchange((*sub)->Om, *map, setting, irank, nrank);
      mpiexchange((*sub)->Vyp, *map, setting, irank, nrank);
      mpiexchange((*sub)->Vym, *map, setting, irank, nrank);
//      mpiexchange((*sub)->Zyp, *map, setting, irank, nrank);
//      mpiexchange((*sub)->Zym, *map, setting, irank, nrank);
    }
    combineSubArea(sub, *data, bath, *map, setting, irank, nrank);
    if (setting->useMPI == 1)
    {
      mpiexchange((*sub)->Nx, *map, setting, irank, nrank);
      mpiexchange((*sub)->Oy, *map, setting, irank, nrank);
//      mpiexchange((*sub)->N, *map, setting, irank, nrank);
//      mpiexchange((*sub)->O, *map, setting, irank, nrank);
      mpiexchange((*sub)->V, *map, setting, irank, nrank);
      mpiexchange((*sub)->Z, *map, setting, irank, nrank);
      mpiexchange((*sub)->Vx, *map, setting, irank, nrank);
      mpiexchange((*sub)->Vy, *map, setting, irank, nrank);
//      mpiexchange((*sub)->Zx, *map, setting, irank, nrank);
//      mpiexchange((*sub)->Zy, *map, setting, irank, nrank);
    }
  }
  // compute cell volumes
  computeVolume(data, setting);
  // initialize for the uuYP and vvXP velocities
  velocityInterp(data, setting);
  // correct the drag coefficient
  thinLayerDrag(data, setting);
  // convert the inflow location into the transposed coordinate
  findInflowLocation(data, *map, setting);
  // initialize scalar mass
  if (setting->useScalar == 1)
  {
    int ii;
    for (ii = 0; ii < setting->N2ci; ii++)
    {if ((*data)->depth[ii] <= 0) {(*data)->S[ii] = 0;}}
    scalarMass(data, setting);
    if (setting->useMPI == 1)
    {mpiexchange((*data)->S, *map, setting, irank, nrank);}
    //enforceScalarBC(data, *bc, *map, setting, 0);
  }
}

// =============== message passing for bathymetry ===============
void enforceBathBC(Bath *bath, Maps *map, Config *setting, int irank, int nrank)
{
  int ii;
  // bathymetry values at ghost cells
  for (ii = 0; ii < setting->nx; ii++)
  {
    bath->bottomZ[map->jMgt[ii]] = bath->bottomZ[map->jMbd[ii]];
    bath->bottomZ[map->jPgt[ii]] = bath->bottomZ[map->jPbd[ii]];
  }
  for (ii = 0; ii < setting->ny; ii++)
  {
    bath->bottomZ[map->iMgt[ii]] = bath->bottomZ[map->iMbd[ii]];
    bath->bottomZ[map->iPgt[ii]] = bath->bottomZ[map->iPbd[ii]];
  }
  // message passing
  if (setting->useMPI == 1)
  {mpiexchange(bath->bottomZ, map, setting, irank, nrank);}
}

// =============== allocate memory for data drrays ===============
void initDataArrays(Data **data, Config *setting)
{
  int ii;
  // initialize the data arrays
  *data = malloc(sizeof(Data));
  (*data)->uuXP = malloc(setting->N2ct*sizeof(double));
  (*data)->uuYP = malloc(setting->N2ci*sizeof(double));
  (*data)->vvYP = malloc(setting->N2ct*sizeof(double));
  (*data)->vvXP = malloc(setting->N2ci*sizeof(double));
  (*data)->Fuu = malloc(setting->N2ct*sizeof(double));
  (*data)->Fvv = malloc(setting->N2ct*sizeof(double));
  (*data)->surf = malloc(setting->N2ct*sizeof(double));
  (*data)->surfOld = malloc(setting->N2ct*sizeof(double));
  (*data)->Vloss = malloc(setting->N2ci*sizeof(double));
  (*data)->depth = malloc(setting->N2ci*sizeof(double));
  (*data)->depthXP = malloc(setting->N2ct*sizeof(double));
  (*data)->depthYP = malloc(setting->N2ct*sizeof(double));
  (*data)->EnXP = malloc(setting->N2ct*sizeof(double));
  (*data)->EnYP = malloc(setting->N2ct*sizeof(double));
  (*data)->DragXP = malloc(setting->N2ct*sizeof(double));
  (*data)->DragYP = malloc(setting->N2ct*sizeof(double));
  (*data)->CDXP = malloc(setting->N2ct*sizeof(double));
  (*data)->CDYP = malloc(setting->N2ct*sizeof(double));
  (*data)->GnXP = malloc(setting->N2ci*sizeof(double));
  (*data)->GnXM = malloc(setting->N2ci*sizeof(double));
  (*data)->GnYP = malloc(setting->N2ci*sizeof(double));
  (*data)->GnYM = malloc(setting->N2ci*sizeof(double));
  (*data)->GnCt = malloc(setting->N2ci*sizeof(double));
  (*data)->alluuXP = malloc(setting->N2CI*sizeof(double));
  (*data)->allvvYP = malloc(setting->N2CI*sizeof(double));
  (*data)->allsurf = malloc(setting->N2CI*sizeof(double));
  (*data)->alldepth = malloc(setting->N2CI*sizeof(double));
  (*data)->allCDXP = malloc(setting->N2CI*sizeof(double));
  (*data)->allCDYP = malloc(setting->N2CI*sizeof(double));
  (*data)->allVloss = malloc(setting->N2CI*sizeof(double));
  (*data)->inflowLoc = malloc(setting->inflowLocLength*sizeof(int));
  (*data)->z = malloc(setting->N2ci*sizeof(double));
  (*data)->wtfXP = malloc(setting->N2ct*sizeof(double));
  (*data)->wtfYP = malloc(setting->N2ct*sizeof(double));
  (*data)->cellV = malloc(setting->N2ci*sizeof(double));
  (*data)->S = malloc(setting->N2ct*sizeof(double));
  (*data)->Sm = malloc(setting->N2ci*sizeof(double));
  (*data)->allS = malloc(setting->N2CI*sizeof(double));
  (*data)->rmvCFL = malloc(setting->N2ci*sizeof(double));
  (*data)->advX = malloc(setting->N2ct*sizeof(double));
  (*data)->advY = malloc(setting->N2ct*sizeof(double));
  (*data)->CFLx = malloc(setting->N2ct*sizeof(double));
  (*data)->CFLy = malloc(setting->N2ct*sizeof(double));
  (*data)->Vloss1 = malloc(1*sizeof(double));
  (*data)->Vloss2 = malloc(1*sizeof(double));
  (*data)->Vloss3 = malloc(1*sizeof(double));
  (*data)->Vloss4 = malloc(1*sizeof(double));
  (*data)->Vloss1[0] = 0;
  (*data)->Vloss2[0] = 0;
  (*data)->Vloss3[0] = 0;
  (*data)->Vloss4[0] = 0;
  (*data)->Qseep = malloc(setting->N2ci*sizeof(double));
  (*data)->Ssub = malloc(setting->N2ci*sizeof(double));
  (*data)->dz = malloc(setting->N2ci*sizeof(double));
  (*data)->rainC = malloc(setting->N2ci*sizeof(double));
//    (*data)->rainSum = malloc(1.0*sizeof(double));
//    (*data)->rainSum[0] = 0.0;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    (*data)->GnCt[ii] = 0;
    (*data)->GnXP[ii] = 0;
    (*data)->GnXM[ii] = 0;
    (*data)->GnYP[ii] = 0;
    (*data)->GnYM[ii] = 0;
    (*data)->z[ii] = 0;
    (*data)->cellV[ii] = 0;
    (*data)->Sm[ii] = 0;
    (*data)->rmvCFL[ii] = 0;
    (*data)->Vloss[ii] = 0;
    (*data)->Qseep[ii] = 0;
	(*data)->dz[ii] = 0;
	(*data)->rainC[ii] = 0;
	if (setting->useSubscalar == 1)
	{(*data)->Ssub[ii] = setting->subS0;}
  }
  for (ii = 0; ii < setting->N2ct; ii++)
  {
    (*data)->Fuu[ii] = 0;
    (*data)->Fvv[ii] = 0;
    (*data)->EnXP[ii] = 0;
    (*data)->EnYP[ii] = 0;
    (*data)->DragXP[ii] = 0;
    (*data)->DragYP[ii] = 0;
    (*data)->wtfXP[ii] = 0;
    (*data)->wtfYP[ii] = 0;
    (*data)->S[ii] = 0;
  }
  for (ii = 0; ii < setting->N2CI; ii++)
  {
    (*data)->alluuXP[ii] = 0;
    (*data)->allvvYP[ii] = 0;
    (*data)->allsurf[ii] = 0;
    (*data)->alldepth[ii] = 0;
    (*data)->allCDXP[ii] = 0;
    (*data)->allCDYP[ii] = 0;
    (*data)->allS[ii] = 0;
    (*data)->allVloss[ii] = 0;
  }
}

// =============== initialize data array for IC ===============
void initICArrays(IC **ic, Config *setting)
{
  *ic = malloc(sizeof(IC));
}

// =============== read initial condition from user setting ===============
void readIC(IC **ic, Bath *bath, Config *setting)
{
  // initial surface elevation
  if (setting->useConstSurf0 == 1)
  {(*ic)->surf0 = setting->initSurf;}
  else
  {
    (*ic)->allSurf0 = malloc(setting->N2CI * sizeof(double));
    char sicname[100];
    strcpy(sicname, setting->inputFolder);
    strcat(sicname, "surfic.dat");
    ReadFile((*ic)->allSurf0, sicname, setting->N2CI);
  }
    // initial velocity
    if (setting->useConstU0 == 1)
    {(*ic)->U0 = setting->initU;}
    else
    {
        (*ic)->allU0 = malloc(setting->N2CI * sizeof(double));
        char sicname[100];
        strcpy(sicname, setting->inputFolder);
        strcat(sicname, "uuic.dat");
        ReadFile((*ic)->allU0, sicname, setting->N2CI);
    }
    if (setting->useConstV0 == 1)
    {(*ic)->V0 = setting->initV;}
    else
    {
        (*ic)->allV0 = malloc(setting->N2CI * sizeof(double));
        char sicname[100];
        strcpy(sicname, setting->inputFolder);
        strcat(sicname, "vvic.dat");
        ReadFile((*ic)->allV0, sicname, setting->N2CI);
    }
    // initial scalar
  if (setting->useScalar == 1)
  {
    if (setting->useConstInitS == 1)
    {(*ic)->S0 = setting->initS;}
    else
    {
      (*ic)->allS0 = malloc(setting->N2CI * sizeof(double));
      char sicname[100];
      strcpy(sicname, setting->inputFolder);
      strcat(sicname, "scalaric.dat");
      ReadFile((*ic)->allS0, sicname, setting->N2CI);
    }
  }
}

// =============== initialize data array for BC ===============
void initBCArrays(BC **bc, Config *setting)
{
  *bc = malloc(sizeof(BC));
  (*bc)->tideP = malloc(setting->Nt*sizeof(double));
  (*bc)->tideM = malloc(setting->Nt*sizeof(double));
  (*bc)->inflow = malloc(setting->Nt*sizeof(double));
  (*bc)->windspd = malloc(setting->Nt*sizeof(double));
  (*bc)->winddir = malloc(setting->Nt*sizeof(double));
  (*bc)->tidalPS = malloc(setting->Nt*sizeof(double));
  (*bc)->tidalMS = malloc(setting->Nt*sizeof(double));
  (*bc)->evap = malloc(setting->Nt*sizeof(double));
  (*bc)->rain = malloc(setting->Nt*sizeof(double));
}

// =============== read boundary conditions from data files ===============
void readBC(BC **bc, Config *setting, int irank)
{
  if (setting->bcType != 2)
  {
    writeText("Loading tide boundary condition ...", irank);
    readOneBC((*bc)->tideP, "tideP.dat", setting, setting->tideNP);
  }
  if (setting->bcType != 1)
  {
    writeText("Loading tide boundary condition ...", irank);
    readOneBC((*bc)->tideM, "tideM.dat", setting, setting->tideNM);
  }
  if (setting->bcType != 3)
  {
    writeText("Loading inflow boundary condition ...", irank);
    readOneBC((*bc)->inflow, "inflow.dat", setting, setting->inflowN);
  }
  if (setting->useWind == 1)
  {
    writeText("Loading wind boundary condition ...", irank);
    readOneBC((*bc)->windspd, "windspd.dat", setting, setting->windspdN);
    readOneBC((*bc)->winddir, "winddir.dat", setting, setting->winddirN);
  }
  if (setting->useConstTidePS == 0)
  {
    writeText("Loading tidal salinity boundary condition ...", irank);
    readOneBC((*bc)->tidalPS, "salinityBC.dat", setting, setting->tidalPSN);
  }
  if (setting->useConstTideMS == 0)
  {
    writeText("Loading tidal salinity boundary condition ...", irank);
    readOneBC((*bc)->tidalMS, "salinityBC.dat", setting, setting->tidalMSN);
  }
  if (setting->useEvap == 1)
  {
    writeText("Loading evaporation boundary condition ...", irank);
    readOneBC((*bc)->evap, "evap.dat", setting, setting->evapN);
  }
  if (setting->useRain == 1)
  {
    writeText("Loading rainfall boundary condition ...", irank);
    readOneBC((*bc)->rain, "rain.dat", setting, setting->rainN);
  }
}

// =============== initialize drag coefficient ================
void initCD(Data **data, Config *setting)
{
  int ii;
  for (ii = 0; ii < setting->N2ct; ii++)
  {
    (*data)->CDXP[ii] = setting->CDx;
    (*data)->CDYP[ii] = setting->CDy;
  }
}

// =============== set initial conditions ===============
void initFieldValues(Data **data, Bath *bath, IC *ic, Config *setting, int irank)
{
  int ii, col;
    double dsurf;
    for (ii = 0; ii < setting->N2ci; ii++)
    {
        // surface
        if (setting->useConstSurf0 == 1)
        {(*data)->surf[ii] = ic->surf0 + bath->offset[0];}
        else
        {(*data)->surf[ii] = ic->allSurf0[ii+irank*setting->N2ci] + bath->offset[0];}
        // velocity
        if (setting->useConstU0 == 1)
        {(*data)->uuXP[ii] = ic->U0;}
        else
        {(*data)->uuXP[ii] = ic->allU0[ii+irank*setting->N2ci];}
        if (setting->useConstV0 == 1)
        {(*data)->vvYP[ii] = ic->V0;}
        else
        {(*data)->vvYP[ii] = ic->allV0[ii+irank*setting->N2ci];}
        // scalar
      if (setting->useScalar == 1)
      {
        if (setting->useConstInitS == 1)
        {(*data)->S[ii] = ic->S0;}
        else
        {(*data)->S[ii] = ic->allS0[ii+irank*setting->N2ci];}
      }
    }
    // spatially-varied initial surface, ZhiLi20180808
//    if (setting->bcType == 3 && setting->useConstSurf0 == 1 && bc->tideP[0] != bc->tideM[0])
//    {
//        dsurf = (bc->tideP[0] - bc->tideM[0]) / setting->NY;
//        for (ii = 0; ii < setting->N2ci; ii ++)
//        {
//            col = floor(ii / setting->NX);
//            (*data)->surf[ii] = bc->tideM[0] + (irank*setting->ny + col)*dsurf + bath->offset[0];
//        }
//    }
}

// =============== read the restart file ===============
//void readRestartFile(IC **ic, Config **setting)
//{
//  int N, n = 3, ii, jj, r = 1, n2ci;
//  double *oldConfig = malloc(n*sizeof(double));
//  ReadFile(oldConfig, (*setting)->restartFile, n);
//  // compute the new model time
//  (*setting)->tNStart = oldConfig[0];
//  (*setting)->tNEnd = (*setting)->tNStart + (*setting)->dt*(*setting)->Nt;
//  // check if the grid size changes, if not then interpolate
//  if (oldConfig[1] != (*setting)->dx)
//  {
//    r = oldConfig[1] / (*setting)->dx;
//    n2ci = (*setting)->N2CI / r / r;
//    N = n2ci * 4 + n;
//    double *allRestartIC = malloc(N*sizeof(double));
//    ReadFile(allRestartIC, (*setting)->restartFile, N);
//    // interpolate
//    interpRestart(allRestartIC, (*ic)->restartSurf0, n, n2ci+n, r, (*setting)->NX);
//    interpRestart(allRestartIC, (*ic)->restartU0, n2ci+n, 2*n2ci+n, r, (*setting)->NX);
//    interpRestart(allRestartIC, (*ic)->restartV0, 2*n2ci+n, 3*n2ci+n, r, (*setting)->NX);
//    interpRestart(allRestartIC, (*ic)->restartS0, 3*n2ci+n, 4*n2ci+n, r, (*setting)->NX);
//  }
//  else
//  {
//    N = (*setting)->N2CI * 4 + n;
//    double *allRestartIC = malloc(N*sizeof(double));
//    ReadFile(allRestartIC, (*setting)->restartFile, N);
//    jj = 0;
//    for (ii = n; ii < (*setting)->N2CI+n; ii++)
//    {(*ic)->restartSurf0[jj] = allRestartIC[ii]; jj++;}
//    jj = 0;
//    for (ii = (*setting)->N2CI+n; ii < 2*(*setting)->N2CI+n; ii++)
//    {(*ic)->restartU0[jj] = allRestartIC[ii]; jj++;}
//    jj = 0;
//    for (ii = 2*(*setting)->N2CI+n; ii < 3*(*setting)->N2CI+n; ii++)
//    {(*ic)->restartV0[jj] = allRestartIC[ii]; jj++;}
//    jj = 0;
//    for (ii = 3*(*setting)->N2CI+n; ii < 4*(*setting)->N2CI+n; ii++)
//    {(*ic)->restartS0[jj] = allRestartIC[ii]; jj++;}
//  }
//
//}

// =============== read one specific BC field ===============
void readOneBC(double *arr, char filename[], Config *setting, int N)
{
  int ind1, ind2, itval, ii;
  char *fullname = malloc(50);
  strcpy(fullname, setting->inputFolder);
  strcat(fullname, filename);
  // load all data for this bc
  double *allBC = malloc(2*N*sizeof(double));
  double *t = malloc(N*sizeof(double));
  double *value = malloc(N*sizeof(double));
  ReadFile(allBC, fullname, 2*N);
  // split time and value
  dataSplit(t, value, allBC, N);
  // find the index of the starting and ending elements
  itval = t[1] - t[0];
  ind1 = searchInd(t, setting->tNStart, itval, N);
  ind2 = searchInd(t, setting->tNEnd, itval, N);
  if (ind1 == ind2)
  {printf("WARNING: Run time is shorter than the BC time intervals...\n");}
  // truncate the bc for the given range
  double *trunkBC = malloc((ind2-ind1)*sizeof(double));
  for (ii = 0; ii < ind2-ind1; ii++)
  {trunkBC[ii] = value[ii+ind1];}
  // interpolate to the required data interval
  dataInterp(arr, t, value, ind1, setting->Nt, setting->dt, setting->tNStart);
  // free temporary arrays
  free(allBC);
  free(t);
  free(value);
  free(trunkBC);
}

// =============== allocate memory for groundwater ===============
void initGroundArrays(Ground **ground, Gmaps *gmap, Bath *bath, Config *setting)
{
    double Sres, a, n, htop;
    int ii;
    *ground = malloc(sizeof(Ground));
    (*ground)->allh = malloc(setting->N3CI*sizeof(double));
    (*ground)->allSw = malloc(setting->N3CI*sizeof(double));
    (*ground)->h = malloc(setting->N3cf*sizeof(double));
    (*ground)->hOld = malloc(setting->N3cf*sizeof(double));
    (*ground)->Sw = malloc(setting->N3cf*sizeof(double));
    (*ground)->Cx = malloc(setting->N3ct*sizeof(double));
    (*ground)->Cy = malloc(setting->N3ct*sizeof(double));
    (*ground)->Cz = malloc(setting->N3ci*sizeof(double));
    (*ground)->Kx = malloc(setting->N3ct*sizeof(double));
    (*ground)->Ky = malloc(setting->N3ct*sizeof(double));
    (*ground)->Kz = malloc(setting->N3ci*sizeof(double));
    (*ground)->V = malloc(setting->N3ci*sizeof(double));
    (*ground)->GnCt = malloc(setting->N3ci*sizeof(double));
    (*ground)->GnXP = malloc(setting->N3ci*sizeof(double));
    (*ground)->GnXM = malloc(setting->N3ci*sizeof(double));
    (*ground)->GnYP = malloc(setting->N3ci*sizeof(double));
    (*ground)->GnYM = malloc(setting->N3ci*sizeof(double));
    (*ground)->GnZP = malloc(setting->N3ci*sizeof(double));
    (*ground)->GnZM = malloc(setting->N3ci*sizeof(double));
    (*ground)->B = malloc(setting->N3ci*sizeof(double));
    (*ground)->nlay = malloc(setting->N2ci*sizeof(double));
    (*ground)->Quu = malloc(setting->N3ct*sizeof(double));
    (*ground)->Qvv = malloc(setting->N3ct*sizeof(double));
    (*ground)->Qww = malloc(setting->N3ct*sizeof(double));
    (*ground)->S = malloc(setting->N3cf*sizeof(double));
    (*ground)->Sm = malloc(setting->N3cf*sizeof(double));
    (*ground)->allS = malloc(setting->N3CI*sizeof(double));
    // get initial head from water table
    for (ii = 0; ii < setting->N3cf; ii++)
    {
      (*ground)->h[ii] = 0.0;
      (*ground)->Sw[ii] = updateSaturation((*ground)->h[ii], setting);
    }
    n = setting->a1;
    Sres = setting->Sres;
    a = setting->a2;
    for (ii = 0; ii < setting->N3ci; ii++)
    {
        if (gmap->actv[ii] == 1 & gmap->bot3d[ii] >= setting->H0)
        {
            htop = bath->bottomZ[gmap->top2D[ii]];
            (*ground)->Sw[ii] = 1.0 - (1.0-Sres)*(gmap->bot3d[ii]-setting->H0)/(htop-setting->H0);
            if ((*ground)->Sw[ii] > 0.999999)  {(*ground)->Sw[ii] = 1.0;}
            if ((*ground)->Sw[ii] < Sres+0.005)  {(*ground)->Sw[ii] = Sres+0.005;}
            (*ground)->h[ii] = -(1.0/a) * pow((pow((1-Sres)/((*ground)->Sw[ii]-Sres),(n/(n-1.0)))-1.0),(1.0/n));
        }
    }
    // initialize
    for (ii = 0; ii < setting->N3ct; ii++)
    {
        (*ground)->Cx[ii] = 0.0;
        (*ground)->Cy[ii] = 0.0;
        // Kx
        if (gmap->iMjckc[ii] >= setting->N3ci)
        {(*ground)->Kx[gmap->iMjckc[ii]] = 0.0;}
        if (gmap->iPjckc[ii] >= setting->N3ci)
        {(*ground)->Kx[ii] = 0.0;}
        else if (gmap->actv[ii] == 0 | gmap->actv[gmap->iPjckc[ii]] == 0)
        {(*ground)->Kx[ii] = 0.0;}
        else
        {(*ground)->Kx[ii] = setting->Kxx;}
        // Ky
        if (gmap->icjMkc[ii] >= setting->N3ci)
        {(*ground)->Ky[gmap->icjMkc[ii]] = 0.0;}
        if (gmap->icjPkc[ii] >= setting->N3ci)
        {(*ground)->Ky[ii] = 0.0;}
        else if (gmap->actv[ii] == 0 | gmap->actv[gmap->icjPkc[ii]] == 0)
        {(*ground)->Ky[ii] = 0.0;}
        else
        {(*ground)->Ky[ii] = setting->Kyy;}
        (*ground)->Quu[ii] = 0.0;
        (*ground)->Qvv[ii] = 0.0;
        (*ground)->Qww[ii] = 0.0;
    }
    for (ii = 0; ii < setting->N3ci; ii++)
    {
        (*ground)->Cz[ii] = 0.0;
        // NOTE: Unlike Kx and Ky, Kz[ii] is its kM face (upward face)
        if (gmap->actv[ii] == 1)
        {(*ground)->Kz[ii] = setting->Kzz;}
        else
        {(*ground)->Kz[ii] = 0.0;}
        (*ground)->V[ii] = 0.0;
        (*ground)->B[ii] = 0.0;
    }
	for (ii = 0; ii < setting->N3cf; ii++)
	{
		// (*ground)->h[ii] = setting->H0;
		(*ground)->hOld[ii] = (*ground)->h[ii];
		(*ground)->S[ii] = setting->subS0;
        (*ground)->Sm[ii] = 0.0;
	}
}
