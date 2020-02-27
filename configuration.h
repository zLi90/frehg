#ifndef CONFIGURATION_H
#define CONFIGURATION_H

typedef struct Config
{
  // -------------------- grid properties --------------------------------------
  // nx, ny = x and y dimensions of computation domain [m]
  int nx, ny;
  // NX, NY = x and y dimensions plus the ghost cells [m]
  int NX, NY;
  // N2ci, N2ct = number of cells without/with ghost cells for each subdomain
  int N2ci, N2ct;
  // N2CI, N2CT = number of cells without/with ghost cells for the entire
  // domain. If MPI is not used, N2CI = N2ci, N2CT = N2CT
  int N2CI, N2CT;
  // dx, dy = grid size in x and y directions [m]
  double dx, dy;
  // -------------------- operation settings -----------------------------------
  // Nt = number of time steps will be simulated
  int Nt;
  // OutItvl = time interval to output the model results [steps]
  int OutItvl;
  // dt = size of the time step [s]
  double dt;
  // tNStart, tNEnd = start and end time of the model represented in C time number
  double tNStart, tNEnd;
  // tStart = start date of model written as a string
  char tStart[11];
  // saveFolder = new folder created to save the output files
  char saveFolder[100], inputFolder[100];
  // use the cell edge file
  int useCellEdge;
  // the variables to be saved
  int savesurface, saveuu, savevv, savedepth, savescalar, saveCD, savesub;
  // whether or not check conservation and save volume loss
  int checkConservation;
  // -------------------- physical properties ----------------------------------
  // g = 9.81
  double g;
  // NUx, NUy = eddy viscosity in x and y directions
  double NUx, NUy;
  // CDx, CDy = bottom drag coefficient in x and y directions
  double CDx, CDy, manningN;
  // use CD, not Manning's n
  int CDnotN;
  // bottom roughness height
  double z0;
  // air density and wind drag coefficient for the wind model
  double rhoa, Cw;
  // -------------------- boundary conditions ----------------------------------
  // note: here we assume the tide is added on the YP boundary, inflow is added
  // on the YM boundary.
  // bcType = 1 means tideP/inflow bc,
  // bcType = 2 means tideM/inflow bc (requires further debugging)
  // bcType = 3 means tideP/tideM bc
  int bcType;
  // tideN = number of tide boundary data in the input file
  int tideNP, tideNM;
  // tideLoc = a vector storing the location where tide bc is added
  int tideLocP[1000], tideLocM[1000];
  // whether to use wind, evaporation and rain model or not
  int useWind;
  // north angle for wind directions
  // definition: angle from negative x clockwise to the north in the model
  double northAngle;
  // same for the inflow bc and wind bc
  int inflowN, inflowLoc[1000], windspdN, winddirN;
  // tide/inflowLocLength = number of cells the tide/inflow bc is applied
  int tideLocLengthP, tideLocLengthM, inflowLocLength;
  // -------------------- initial conditions -----------------------------------
  // initU, initV, initSurf = initial condition for U, V, and free surface
  double initU, initV, initSurf;
  int useConstSurf0, useConstU0, useConstV0;
  // -------------------- scalar settings --------------------------------------
  // use constant or time dependent salinity boundary conditions
  int useScalar, useConstTidePS, useConstTideMS, useConstInitS, useConstInflowS;
  int tidalPSN, tidalMSN, scalarAdv, useSubStep;
  double initS, tidePS, tideMS, inflowS, Kx, Ky;
  // -------------------- subgrid settings -------------------------------------
  // use subgrid model or not
  int useSubgrid;
  // use the subgrid drag model or not
  int useSubDrag;
  // coefficient of the subgrid drag method by ZhiLi20180329
  double lambda1, lambda2, beta;
  // use p-c method to correct the nonlinear subgrid areas
  int useCorrector;
  char subgridFolder[50];
  // the fine grid size and area
  double dxf, dyf, dA, subR;
  // the range and resolution of the surface elevation vector
  double surfmax, surfmin, dsurf;
  // use staggered subgrid volume or not
  int staggeredV;
  // use min face area
  int useminA, phiSurface, phiNonlinear;
  // -------------------- MPI settings -----------------------------------------
  int useMPI, np;
  // -------------------- setting for performance ------------------------------
  // residual for pcg
  double eps;
  // maximum iteration for pcg
  int maxIter;
  // minimum cell depth
  double minDepth;
  // thin layer model
  int useThinLayer;
  double CDmax, hD, CwT;
  // waterfall model
  double wtfh;
  // CFL limiter for the advection term
  double CFLl, CFLh;
  // -------------------- setting for groundwater ------------------------------
    int useSubsurface, useSubscalar, N3ci, N3ct, N3cf, N3ca, N3CI, useUnSat;
    int topBC, botBC;
    double zbot, layZ, porosity, Kxx, Kyy, Kzz, H0, dtg;
    double subS0, Kmx, Kmy, Kmz;
	double a1, a2, Sres;
  // -------------------- setting for evap/rain -------------------------------
	int useEvap, useRain, evapN, rainN, eTstart, eTend, rTstart, rTend;
    double qRain, qEvap, kkext;
}Config;

#endif

//void ReadUserSettings(Config **setting);
void ReadUserSettings(Config **setting, double *configArr, char *date, char *inputFolder, char *outputFolder, char *subgridFolder);
