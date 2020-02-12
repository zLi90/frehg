#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

#include "configuration.h"
#include "utilities.h"

// -----------------------------------------------------------------------------
// This file contains all user-defined parameters for the FREHD model. This
// should be the only file user needs to modify. A complete description of
// all the fields in setting could be found in configuration.h
// - Zhi Li 2017-05-05 -
// -----------------------------------------------------------------------------

// ========= Assign setting values to its fields =========
void ReadUserSettings(Config **setting, double *configArr, char *date, char *inputFolder, char *outputFolder, char *subgridFolder)
{
    int ii, jj, kk = 0;
    *setting = malloc(sizeof(Config));
    // -------------------- GRID SETTINGS ----------------------------------------
    (*setting)->dx = configArr[4];    (*setting)->dy = configArr[5];
    (*setting)->NX = configArr[6];  (*setting)->NY = configArr[7];
    // -------------------- OPERATION SETTINGS -----------------------------------
    (*setting)->dt = configArr[8];    (*setting)->Nt = configArr[9]*86400/(*setting)->dt;
    (*setting)->OutItvl = configArr[11]*3600/(*setting)->dt;
    strcpy((*setting)->tStart, date);
    strcpy((*setting)->saveFolder, outputFolder);
    strcpy((*setting)->inputFolder, inputFolder);
    (*setting)->useCellEdge = configArr[10];
    (*setting)->savesurface = configArr[12];
    (*setting)->saveuu = configArr[13];
    (*setting)->savevv = configArr[14];
    (*setting)->savedepth = configArr[15];
    (*setting)->savescalar = configArr[16];
    (*setting)->saveCD = configArr[17];
    (*setting)->savesub = configArr[18];
    // -------------------- BOUNDARY CONDITIONS ----------------------------------
    (*setting)->bcType = configArr[19];
    // --- Tide P ---
    (*setting)->tideLocLengthP = configArr[20];;
    int tideLocP[(*setting)->tideLocLengthP];
    for (ii = 0; ii < (*setting)->tideLocLengthP; ii++)
    {tideLocP[ii] = configArr[21] + ii;}
    (*setting)->tideNP = configArr[22];
    // --- Tide M ---
    (*setting)->tideLocLengthM = configArr[23];;
    int tideLocM[(*setting)->tideLocLengthM];
    for (ii = 0; ii < (*setting)->tideLocLengthM; ii++)
    {tideLocM[ii] = configArr[24] + ii;}
    (*setting)->tideNM = configArr[25];
    // --- Inflow ---
    (*setting)->inflowLocLength = configArr[26];
    int inflowLoc[(*setting)->inflowLocLength];
    for (ii = configArr[27]; ii < configArr[28]; ii++)
    {
        for (jj = configArr[29]; jj < configArr[30]; jj++)
        {inflowLoc[kk] = jj * (*setting)->NX + ii - 1;      kk += 1;}
    }
    (*setting)->inflowN = configArr[31];
    // --- Wind ---
    (*setting)->useWind = configArr[32];
    (*setting)->northAngle = configArr[33];
    (*setting)->windspdN = configArr[34];
    (*setting)->winddirN = configArr[35];
    // --- Scalar ---
    (*setting)->useScalar = configArr[36];
    (*setting)->scalarAdv = configArr[37];
    (*setting)->useConstInflowS = configArr[38];
    (*setting)->inflowS = configArr[39];
    (*setting)->useConstTidePS = configArr[40];
    (*setting)->tidePS = configArr[41];
    (*setting)->tidalPSN = configArr[42];
    (*setting)->useConstTideMS = configArr[43];
    (*setting)->tideMS = configArr[44];
    (*setting)->tidalMSN = configArr[45];
    // -------------------- INITIAL CONDITIONS -----------------------------------
    (*setting)->initU = configArr[46];
    (*setting)->initV = configArr[47];
    (*setting)->initSurf = configArr[48];
    (*setting)->initS = configArr[49];
    (*setting)->useConstSurf0 = configArr[50];
    (*setting)->useConstU0 = configArr[51];
    (*setting)->useConstV0 = configArr[52];
    (*setting)->useConstInitS = configArr[53];
    // -------------------- SUBGRID SETTINGS -------------------------------------
    (*setting)->useSubgrid = configArr[54];
    (*setting)->useSubDrag = configArr[55];
    strcpy((*setting)->subgridFolder, subgridFolder);
    (*setting)->dxf = configArr[56];
    (*setting)->dyf = configArr[57];
    (*setting)->dA = (*setting)->dxf * (*setting)->dyf;
    (*setting)->subR = (*setting)->dx / (*setting)->dxf;
    (*setting)->surfmax = configArr[58];
    (*setting)->surfmin = configArr[59];
    (*setting)->dsurf = configArr[60];
    (*setting)->staggeredV = configArr[61];
    (*setting)->useminA = configArr[62];
    (*setting)->beta = configArr[63];
    // -------------------- MPI SETTINGS -----------------------------------------
    (*setting)->useMPI = configArr[64];
    (*setting)->np = configArr[65];
    // -------------------- PHYSICAL PROPERTIES ----------------------------------
    (*setting)->g = configArr[66];
    (*setting)->NUx = configArr[67];
    (*setting)->NUy = configArr[68];
    (*setting)->CDnotN = configArr[69];
    (*setting)->CDx = configArr[70];
    (*setting)->CDy = configArr[71];
    (*setting)->manningN = configArr[72];
    (*setting)->useThinLayer = configArr[73];
    (*setting)->z0 = configArr[74];
    (*setting)->hD = (*setting)->z0;
    (*setting)->CDmax = configArr[75];
    (*setting)->CwT = configArr[76];
    (*setting)->rhoa = configArr[77];
    (*setting)->Cw = configArr[78];
    (*setting)->Kx = configArr[79];
    (*setting)->Ky = configArr[80];
    (*setting)->minDepth = configArr[81];
    (*setting)->wtfh = configArr[82];
    (*setting)->CFLl = configArr[83];
    (*setting)->CFLh = configArr[84];
    (*setting)->eps = configArr[85];
    (*setting)->maxIter = configArr[86];
  // -------------------- GROUNDWATER ----------------------------------
    (*setting)->useSubsurface = configArr[87];
    (*setting)->zbot = configArr[88];
    (*setting)->layZ = configArr[89];
    (*setting)->dtg = configArr[90];
    (*setting)->porosity = configArr[91];
    (*setting)->compS = configArr[92];
    (*setting)->compW = configArr[93];
    (*setting)->Kxx = configArr[94];
    (*setting)->Kyy = configArr[95];
    (*setting)->Kzz = configArr[96];
    (*setting)->H0 = configArr[97];
    (*setting)->useSubscalar = configArr[98];
    (*setting)->subS0 = configArr[99];
	(*setting)->Kmx = configArr[100];
	(*setting)->Kmy = configArr[101];
	(*setting)->Kmz = configArr[102];
	(*setting)->a1 = configArr[103];
	(*setting)->a2 = configArr[104];
	(*setting)->Sres = configArr[105];
  (*setting)->useUnSat = configArr[106];
	// ------------------------ EVAP/RAIN --------------------------------
	(*setting)->useRain = configArr[107];
    (*setting)->qRain = configArr[108];
    (*setting)->rTstart = configArr[109];
    (*setting)->rTend = configArr[120];
    (*setting)->rainN = configArr[121];
	(*setting)->useEvap = configArr[122];
	(*setting)->qEvap = configArr[123];
    (*setting)->eTstart = configArr[124];
    (*setting)->eTend = configArr[125];
	(*setting)->evapN = configArr[126];
    (*setting)->kkext = configArr[111];
  // ---------------------------------------------------------------------------
  // -------------------- NO USER CHANGE BELOW ---------------------------------
  // ---------------------------------------------------------------------------
  // ---------- Derived Values ----------
  // sizes of data arrays
  (*setting)->N2CI = (*setting)->NX*(*setting)->NY;
  (*setting)->N2CT = ((*setting)->NX+2)*((*setting)->NY+2);
  if ((*setting)->NY % (*setting)->np != 0)
  {printf("WARNING: The computation domain cannot be evenly divided by nranks!\n");}
  else
  {
    (*setting)->nx = (*setting)->NX;    (*setting)->ny = round((*setting)->NY/(*setting)->np);
    (*setting)->N2ci = (*setting)->nx*(*setting)->ny;
    (*setting)->N2ct = ((*setting)->nx+2)*((*setting)->ny+2);
  }
  // start and end time represented by date number
  (*setting)->tNStart = dateNum((*setting)->tStart);
  (*setting)->tNEnd = (*setting)->tNStart + (*setting)->dt*(*setting)->Nt;
  // locations of boundary conditions
  for (ii = 0; ii < (*setting)->tideLocLengthP; ii++)
  {(*setting)->tideLocP[ii] = tideLocP[ii] + (*setting)->nx*((*setting)->ny-1);}
  for (ii = 0; ii < (*setting)->tideLocLengthM; ii++)
  {(*setting)->tideLocM[ii] = tideLocM[ii];}
  for (ii = 0; ii < (*setting)->inflowLocLength; ii++)
  {(*setting)->inflowLoc[ii] = inflowLoc[ii];}// + (*setting)->nx;}
  // check if subgrid bathymetry can be created
  if ((*setting)->useSubgrid == 1)
  {
    if (fmod((*setting)->dx,(*setting)->dxf) != 0 | fmod((*setting)->dy,(*setting)->dyf) != 0)
    {printf("WARNING: The coarse grid cannot be divided by the fine grid!\n");}
  }
}
