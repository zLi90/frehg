
#include "bathymetry.h"
#include "configuration.h"
#include "map.h"
#include "subgrid.h"

#ifndef INITIALIZE_H
#define INITIALIZE_H

typedef struct Data
{
  double *uuXP, *uuYP, *vvXP, *vvYP, *surf, *surfOld, *Fuu, *Fvv, *cellV;
  double *alluuXP, *allvvYP, *allsurf, *alldepth, *allS, *allCDXP, *allCDYP;
  double *depth, *depthXP, *depthYP;
  double *EnXP, *EnYP, *DragXP, *DragYP, *CDXP, *CDYP;
  double *GnCt, *GnXP, *GnYP, *GnXM, *GnYM;
  double *z, *wtfXP, *wtfYP, *S, *Sm, *rmvCFL;
  double *advX, *advY, *CFLx, *CFLy;
  double *Vloss, *allVloss, *Vloss1, *Vloss2, *Vloss3, *Vloss4, *Qseep, *Ssub, *dz, *rainC;
  int *inflowLoc;
}Data;

typedef struct IC
{
  double surf0, U0, V0, S0, *allS0, *allSurf0, *allU0, *allV0, *restartU0, *restartV0, *restartSurf0, *restartS0;
}IC;

typedef struct BC
{
  double *tideP, *tideM, *inflow, *windspd, *winddir, *tidalPS, *tidalMS;
  double *evap, *rain;
}BC;

typedef struct Ground
{
    double H0, SS, *allh0, *allh, *h, *hOld, *Cx, *Cy, *Cz, *Kx, *Ky, *Kz, *V, *B, *Sw, *allSw;
    double *GnCt, *GnXP, *GnXM, *GnYP, *GnYM, *GnZP, *GnZM;
    double *nlay, *Quu, *Qvv, *Qww;
    double *S, *Sm, *allS;
}Ground;


#endif

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


//#endif
