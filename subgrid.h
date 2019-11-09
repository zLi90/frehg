//#include "configuration.h"
//#include "initialize.h"
//#include "map.h"

#ifndef SUBGRID_H
#define SUBGRID_H

typedef struct Sub
{
  double *bath, *bathXP, *bathYP, *allSurf;
  double *allV, *allN, *allO, *allZ, *allVxp, *allVyp, *allNp, *allOp, *allZxp, *allZyp;
  double *allCvX, *allCvY, *allYh, *allredCD;
  double *V, *N, *O, *Z, *Vx, *Vy, *Zx, *Zy, *Nx, *Oy;
  double *allVxm, *allVym, *allNm, *allOm, *allZxm, *allZym;
  double *Vxp, *Vyp, *Vxm, *Vym, *Zxp, *Zyp, *Zxm, *Zym, *Np, *Nm, *Op, *Om, *CvX, *CvY, *Yh, *redCD;
  double *wdNp, *wdOp, *wdNm, *wdOm, *wdNpAR, *wdOpAR, *wdNmAR, *wdOmAR;
  int *ind, *biX1, *biX2, *biY1, *biY2, *Exflag, *Eyflag;
  int Ne;
  double *allNx, *allOy, *allVx, *allVy;
  double *bathAR, *bathXPAR, *bathYPAR;
  double *allVAR, *allNAR, *allOAR, *allZAR, *allVxpAR, *allVypAR, *allNpAR, *allOpAR, *allZxpAR, *allZypAR;
  double *allVxmAR, *allVymAR, *allNmAR, *allOmAR, *allZxmAR, *allZymAR, *allCvXAR, *allCvYAR;
  double *biX1AR, *biX2AR, *biY1AR, *biY2AR, *ExflagAR, *EyflagAR, *allYhAR, *allredCDAR;
  double *CdsX, *CdsY, *allCdsX, *allCdsY, *allCdsXAR, *allCdsYAR;
  double *allCnAR, *allCn, *Cn;
    double *allvar, *allvarAR, *allVo, *allZo, *allNmin, *allOmin, *Nmin, *Omin;
    double *allphiX, *allphiY, *phiX, *phiY, *allZ0, *Z0;
}Sub;

typedef struct Data Data;

typedef struct Maps Maps;

typedef struct Config Config;

typedef struct BC BC;

typedef struct IC IC;

typedef struct Bath bath;

#endif

void initOneSubVar(Sub **sub, char vname[], Config *setting, int irank);
void initAllSubVar(Sub **sub, Bath *bath, Config *setting, int irank);
void initSubArea(Sub **sub, Bath *bath, Config *setting, int irank);
void readSubArea(Sub **sub, Bath *bath, Config *setting, int irank);
void splitSubArea(Sub **sub, Bath *bath, Config *setting, int irank);
void computeSubArea(Sub **sub, Data *data, Bath *bath, Maps *map, Config *setting);
void combineSubArea(Sub **sub, Data *data, Bath *bath, Maps *map, Config *setting, int irank, int nrank);
void updateAllFaceDepthSub(Data **data, Sub *sub, Maps *map, Config *setting, int irank, int nrank);
void advectionTermSub(Data **data, Maps *map, Sub *sub, Config *setting);
void diffusionTermSub(Data **data, Maps *map, Sub *sub, Config *setting);
void windTermSub(Data **data, BC *bc, Sub *sub, Config *setting, int tt);
void dragTermSub(Data **data, Maps *map, Sub *sub, Config *setting);
void matrixSourceTermSub(Data **data, Maps *map, BC *bc, Sub *sub, Config *setting, int tt, int irank);
void matrixCoeffSub(Data **data, Maps *map, Sub *sub, Config *setting);
void updateVelocitySub(Data **data, Maps *map, Sub *sub, Config *setting);
void computeFaceFlowRateSub(Data **data, Maps *map, Sub *sub, Config *setting);
void computeVolumeSub(Data **data, Sub *sub, Config *setting);
void volumeByFluxSub(Data **data, Sub *sub, Maps *map, BC *bc, Config *setting, int tt, int irank);
void velocityInterpSub(Data **data, Sub *sub, Config *setting);
void adjustTidalVelocitySub(Data **data, Maps *map, Sub *sub, Config *setting, int irank, int nrank);
void detectWaterfallLocationSub(Data **data, Sub *sub, Bath *bath, Maps *map, Config *setting);
void waterfallCorrectionSub(Data **data, Sub *sub, Bath *bath, Maps *map, Config *setting);
void updateCDSub(Data **data, Sub *sub, Maps *map, Config *setting);
void reduceCDSub(Data **data, Sub *sub, Config *setting);
void scalarDiffusionSub(Data **data, Maps *map, Sub *sub, Config *setting);
void updateScalarSub(Data **data, BC *bc, Bath *bath, Maps *map, Sub *sub, Config *setting, int tt, int irank, int nrank);
void extendSurface(Data *data, Sub *sub, double *eqsurf, int *j1, int *j2, int ii);
