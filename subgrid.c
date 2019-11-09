

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<sys/stat.h>

// -----------------------------------------------------------------------------
// This file contains functions required for the subgrid bathymetry model.
// - Zhi Li 2017-07-09 -
// -----------------------------------------------------------------------------

#include "bathymetry.h"
#include "configuration.h"
#include "fileio.h"
#include "initialize.h"
#include "map.h"
#include "nsfunctions.h"
#include "mpifunctions.h"
#include "subgrid.h"
#include "utilities.h"

void initOneSubVar(Sub **sub, char vname[], Config *setting, int irank);
void initAllSubVar(Sub **sub, Bath *bath, Config *setting, int irank);
void initSubArea(Sub **sub, Bath *bath, Config *setting, int irank);
void readSubArea(Sub **sub, Bath *bath, Config *setting, int irank);
void splitSubArea(Sub **sub, Bath *bath, Config *setting, int irank);
void computeSubArea(Sub **sub, Data *data, Bath *bath, Maps *map, Config *setting);
void combineSubArea(Sub **sub, Data *data, Bath *, Maps *map, Config *setting, int irank, int nrank);
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


// ----- extract the subgrid data for current surface elevation -----
void computeSubArea(Sub **sub, Data *data, Bath *bath, Maps *map, Config *setting)
{
  /*
    Illustration of this function:
    This function is invoked when the surface elevations are updated. It updates
    the subgrid variables based on the surface elevations. Note that only the
    subgrid variables for the half cells (or the 'staggered' cells) are updated.
    The next function ('combineSubArea') is used to combine the subgrid variables
    of each half cell, generating the final subgrid variables required in the
    discretization scheme.
    Definition of variables:
    allSurf - This is a vector that stores all reference surface elevations. The
              elevations values are spaced by 'dsurf'.
    ind - Every cell has an 'ind', which is the index in allSurf that is closest
          to the true surface elevation of this cell.
    j1, j2 - Since the true surface elevation often lies between two reference
             elevations in allSurf, j1 and j2 are the index of the reference
             elevations that above and below the true elevation. j1 = ind.
    k1, k2 - A subgrid variable for all possible surface elevations are stored
            in one 1D vector. k1 and k2 are the corresponding index of j1 and j2
            in that 1D vector. To compute k1 and k2, one needs to know j1, j2 as
            well as the index of the cell in the computation domain ('ii').
    bottomXP, YP - Bottom elevations of the cell edges. These values are computed
                   only using the two neighboring rows at the cell edges.
    wdNp, Nm, Op, Om - Bottom elevations of the cell interiors. There exists
            situations where the cell edge is inundated by the interior is
            higher thus remains dry. These variables are used to prohibit
            inundation of the entire cell when such situation happens.
  */
  // find the index for the surface elevation
  int ii, j1, j2, k1, k2;
  double wdNp, wdNm, wdOp, wdOm;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    // for cell ii, find the index in 'allSurf' corresponds to its free surface
    // elevation 'surf[ii]', save this index as 'ind[ii]'.
    while (fabs((*sub)->allSurf[(*sub)->ind[ii]]-data->surf[ii]) > \
        0.5*setting->dsurf + setting->dsurf/100)
    {
      // searching for the index
      if (data->surf[ii] > (*sub)->allSurf[(*sub)->ind[ii]])
      {(*sub)->ind[ii] += 1;}
      else
      {(*sub)->ind[ii] -= 1;}
      // stop when the index is out of bound
      if ((*sub)->ind[ii] < 0)
      {(*sub)->ind[ii] = 0; break;}
      if ((*sub)->ind[ii] > (*sub)->Ne - 1)
      {(*sub)->ind[ii] = (*sub)->Ne - 1; break;}
    }
    // interpolate for the subgrid areas and volumes
    j1 = (*sub)->ind[ii];
    if (j1 == 0)
    {j2 = 0;}
    else if (j1 == (*sub)->Ne - 1)
    {j2 = (*sub)->Ne - 1;}
    else
    {
      if (data->surf[ii] > (*sub)->allSurf[(*sub)->ind[ii]])
      {j2 = j1 + 1;}
      else
      {j2 = j1 - 1;}
    }
    k1 = setting->N2ci*j1 + ii;
    k2 = setting->N2ci*j2 + ii;
    if (setting->phiSurface == 1 | setting->phiNonlinear == 1)
    {
        (*sub)->phiX[ii] = (*sub)->allphiX[k1];
        (*sub)->phiY[ii] = (*sub)->allphiY[k1];
    }
    if (setting->useSubDrag == 2)
    {(*sub)->Yh[ii] = (*sub)->allYh[k1];}
    if (setting->useSubDrag == 3)
    {
      (*sub)->CvX[ii] = (*sub)->allCvX[k1];
      (*sub)->CvY[ii] = (*sub)->allCvY[k1];
      (*sub)->CdsX[ii] = (*sub)->allCdsX[k1];
      (*sub)->CdsY[ii] = (*sub)->allCdsY[k1];
    }
    if (setting->useSubDrag == 4)
    {(*sub)->Cn[ii] = (*sub)->allCn[k1];}
    // (*sub)->redCD[ii] = (*sub)->allredCD[k1];
    // variables on cell centers
    (*sub)->V[ii] = interpSurf(data->surf[ii],(*sub)->V[ii],(*sub)->allSurf[j1], \
      (*sub)->allV[k1],(*sub)->allSurf[j2],(*sub)->allV[k2]);
    // (*sub)->N[ii] = interpSurf(data->surf[ii],(*sub)->N[ii],(*sub)->allSurf[j1], \
    //   (*sub)->allN[k1],(*sub)->allSurf[j2],(*sub)->allN[k2]);
    // (*sub)->O[ii] = interpSurf(data->surf[ii],(*sub)->O[ii],(*sub)->allSurf[j1], \
    //   (*sub)->allO[k1],(*sub)->allSurf[j2],(*sub)->allO[k2]);
    (*sub)->Z[ii] = interpSurf(data->surf[ii],(*sub)->Z[ii],(*sub)->allSurf[j1], \
      (*sub)->allZ[k1],(*sub)->allSurf[j2],(*sub)->allZ[k2]);
    // volumes on cell faces
    (*sub)->Vxp[ii] = interpSurf(data->surf[ii],(*sub)->Vxp[ii],(*sub)->allSurf[j1], \
      (*sub)->allVxp[k1],(*sub)->allSurf[j2],(*sub)->allVxp[k2]);
    (*sub)->Vyp[ii] = interpSurf(data->surf[ii],(*sub)->Vyp[ii],(*sub)->allSurf[j1], \
      (*sub)->allVyp[k1],(*sub)->allSurf[j2],(*sub)->allVyp[k2]);
    (*sub)->Vxm[ii] = interpSurf(data->surf[ii],(*sub)->Vxm[ii],(*sub)->allSurf[j1], \
      (*sub)->allVxm[k1],(*sub)->allSurf[j2],(*sub)->allVxm[k2]);
    (*sub)->Vym[ii] = interpSurf(data->surf[ii],(*sub)->Vym[ii],(*sub)->allSurf[j1], \
      (*sub)->allVym[k1],(*sub)->allSurf[j2],(*sub)->allVym[k2]);
    // surface areas on cell faces
    // (*sub)->Zxp[ii] = interpSurf(data->surf[ii],(*sub)->Zxp[ii],(*sub)->allSurf[j1], \
    //   (*sub)->allZxp[k1],(*sub)->allSurf[j2],(*sub)->allZxp[k2]);
    // (*sub)->Zyp[ii] = interpSurf(data->surf[ii],(*sub)->Zyp[ii],(*sub)->allSurf[j1], \
    //   (*sub)->allZyp[k1],(*sub)->allSurf[j2],(*sub)->allZyp[k2]);
    // (*sub)->Zxm[ii] = interpSurf(data->surf[ii],(*sub)->Zxm[ii],(*sub)->allSurf[j1], \
    //   (*sub)->allZxm[k1],(*sub)->allSurf[j2],(*sub)->allZxm[k2]);
    // (*sub)->Zym[ii] = interpSurf(data->surf[ii],(*sub)->Zym[ii],(*sub)->allSurf[j1], \
    //   (*sub)->allZym[k1],(*sub)->allSurf[j2],(*sub)->allZym[k2]);
    // face areas on cell faces
    /*wdNp = (*sub)->wdNp[ii] + bath->offset[0];
    wdNm = (*sub)->wdNm[ii] + bath->offset[0];
    wdOp = (*sub)->wdOp[ii] + bath->offset[0];
    wdOm = (*sub)->wdOm[ii] + bath->offset[0];*/
    wdNp = bath->edgeX[ii];
    wdNm = bath->edgeX[map->iMjc[ii]];
    wdOp = bath->edgeY[ii];
    wdOm = bath->edgeY[map->icjM[ii]];
    // Np
    if (wdNp > (*sub)->allSurf[j1] & wdNp < (*sub)->allSurf[j2])
    {
      if (wdNp > data->surf[ii])
      {(*sub)->Np[ii] = 0;}
      else
      {(*sub)->Np[ii] = interpSurf(data->surf[ii],(*sub)->Np[ii],wdNp, \
        0.0,(*sub)->allSurf[j2],(*sub)->allNp[k2]);}
    }
    else
    {
      (*sub)->Np[ii] = interpSurf(data->surf[ii],(*sub)->Np[ii],(*sub)->allSurf[j1], \
        (*sub)->allNp[k1],(*sub)->allSurf[j2],(*sub)->allNp[k2]);
    }
    // Nm
    if (wdNm > (*sub)->allSurf[j1] & wdNm < (*sub)->allSurf[j2])
    {
      if (wdNm > data->surf[ii])
      {(*sub)->Nm[ii] = 0;}
      else
      {(*sub)->Nm[ii] = interpSurf(data->surf[ii],(*sub)->Nm[ii],wdNm, \
        0.0,(*sub)->allSurf[j2],(*sub)->allNm[k2]);}
    }
    else
    {
      (*sub)->Nm[ii] = interpSurf(data->surf[ii],(*sub)->Nm[ii],(*sub)->allSurf[j1], \
        (*sub)->allNm[k1],(*sub)->allSurf[j2],(*sub)->allNm[k2]);
    }
    // Op
    if (wdOp > (*sub)->allSurf[j1] & wdOp < (*sub)->allSurf[j2])
    {
      if (wdOp > data->surf[ii])
      {(*sub)->Op[ii] = 0;}
      else
      {(*sub)->Op[ii] = interpSurf(data->surf[ii],(*sub)->Op[ii],wdOp, \
        0.0,(*sub)->allSurf[j2],(*sub)->allOp[k2]);}
    }
    else
    {
      (*sub)->Op[ii] = interpSurf(data->surf[ii],(*sub)->Op[ii],(*sub)->allSurf[j1], \
        (*sub)->allOp[k1],(*sub)->allSurf[j2],(*sub)->allOp[k2]);
    }
    // Om
    if (wdOm > (*sub)->allSurf[j1] & wdOm < (*sub)->allSurf[j2])
    {
      if (wdOm > data->surf[ii])
      {(*sub)->Om[ii] = 0;}
      else
      {(*sub)->Om[ii] = interpSurf(data->surf[ii],(*sub)->Om[ii],wdOm, \
        0.0,(*sub)->allSurf[j2],(*sub)->allOm[k2]);}
    }
    else
    {
      (*sub)->Om[ii] = interpSurf(data->surf[ii],(*sub)->Om[ii],(*sub)->allSurf[j1], \
        (*sub)->allOm[k1],(*sub)->allSurf[j2],(*sub)->allOm[k2]);
    }
    // Nmin and Omin, added by ZhiLi20181004
      (*sub)->Nmin[ii] = interpSurf(data->surf[ii],(*sub)->Nmin[ii],(*sub)->allSurf[j1], \
                                    (*sub)->allNmin[k1],(*sub)->allSurf[j2],(*sub)->allNmin[k2]);
      (*sub)->Omin[ii] = interpSurf(data->surf[ii],(*sub)->Omin[ii],(*sub)->allSurf[j1], \
                                  (*sub)->allOmin[k1],(*sub)->allSurf[j2],(*sub)->allOmin[k2]);


    /*(*sub)->Np[ii] = interpSurf(data->surf[ii],(*sub)->Np[ii],(*sub)->allSurf[j1], \
      (*sub)->allNp[k1],(*sub)->allSurf[j2],(*sub)->allNp[k2]);
    (*sub)->Op[ii] = interpSurf(data->surf[ii],(*sub)->Op[ii],(*sub)->allSurf[j1], \
      (*sub)->allOp[k1],(*sub)->allSurf[j2],(*sub)->allOp[k2]);
    (*sub)->Nm[ii] = interpSurf(data->surf[ii],(*sub)->Nm[ii],(*sub)->allSurf[j1], \
      (*sub)->allNm[k1],(*sub)->allSurf[j2],(*sub)->allNm[k2]);
    (*sub)->Om[ii] = interpSurf(data->surf[ii],(*sub)->Om[ii],(*sub)->allSurf[j1], \
      (*sub)->allOm[k1],(*sub)->allSurf[j2],(*sub)->allOm[k2]);*/
    // zero out areas for dry cell
    if (data->depth[ii] == 0)
    {
      (*sub)->Vxp[ii] = 0;
      (*sub)->Vyp[ii] = 0;
      (*sub)->Vxm[ii] = 0;
      (*sub)->Vym[ii] = 0;
      // (*sub)->Zxp[ii] = 0;
      // (*sub)->Zyp[ii] = 0;
      // (*sub)->Zxm[ii] = 0;
      // (*sub)->Zym[ii] = 0;
    }
    if (data->surf[ii] - bath->bottomXP[ii] <= 0)
    {(*sub)->Np[ii] = 0;}
    if (data->surf[ii] - bath->bottomYP[ii] <= 0)
    {(*sub)->Op[ii] = 0;}
    if (map->iMjc[ii] < setting->N2ci)
    {if (data->surf[ii] - bath->bottomXP[map->iMjc[ii]] <= 0) {(*sub)->Nm[ii] = 0;}}
    if (map->icjM[ii] < setting->N2ci)
    {if (data->surf[ii] - bath->bottomYP[map->icjM[ii]] <= 0) {(*sub)->Om[ii] = 0;}}
    /*if (data->surf[ii] <= (*sub)->wdNp[ii] + bath->offset[0])
    {(*sub)->Np[ii] = 0;}
    if (data->surf[ii] <= (*sub)->wdNm[ii] + bath->offset[0])
    {(*sub)->Nm[ii] = 0;}
    if (data->surf[ii] <= (*sub)->wdOp[ii] + bath->offset[0])
    {(*sub)->Op[ii] = 0;}
    if (data->surf[ii] <= (*sub)->wdOm[ii] + bath->offset[0])
    {(*sub)->Om[ii] = 0;}*/
  }
}

// =============== Combine subgrid area of the half cells ===============
void combineSubArea(Sub **sub, Data *data, Bath *bath, Maps *map, Config *setting, int irank, int nrank)
{
  int ii, j1, j2, k1, k2, v1, v2, eqind, flag = 0;
  //double *eqsurf = malloc(sizeof(double));
  double eqN, realN, eqV;
  double Vmin = 0.002;
  double Amin = 0.002;
  double Azmin = 0.002;
  // combine the half-cell areas and volumes
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    // --- for areas (N and O), choose the larger between the two
    // in x direction
    if (map->iPjc[ii] < setting->N2ci)
    {
      if ((*sub)->Np[ii] >= (*sub)->Nm[map->iPjc[ii]])
      {(*sub)->Nx[ii] = (*sub)->Np[ii];}
      else
      {(*sub)->Nx[ii] = (*sub)->Nm[map->iPjc[ii]];}
    }
    else
    {
      (*sub)->Nx[ii] = (*sub)->Np[ii];
    }
    // in y direction
    if (map->icjP[ii] < setting->N2ci)
    {
      if ((*sub)->Op[ii] >= (*sub)->Om[map->icjP[ii]])
      {(*sub)->Oy[ii] = (*sub)->Op[ii];}
      else
      {(*sub)->Oy[ii] = (*sub)->Om[map->icjP[ii]];}
    }
    else
    {
      (*sub)->Oy[ii] = (*sub)->Op[ii];
    }

    // --- volumes V and surface areas Z are simply added ---
    // for x direction
    if (map->iPjc[ii] < setting->N2ci)
    {
      if ((*sub)->Vxp[ii] > 0 & (*sub)->Vxm[map->iPjc[ii]] > 0)
      {(*sub)->Vx[ii] = (*sub)->Vxp[ii] + (*sub)->Vxm[map->iPjc[ii]];}
      else
      {(*sub)->Vx[ii] = 2.0*((*sub)->Vxp[ii] + (*sub)->Vxm[map->iPjc[ii]]);}
      //(*sub)->Zx[ii] = (*sub)->Zxp[ii] + (*sub)->Zxm[map->iPjc[ii]];
      //(*sub)->Vx[ii] = (*sub)->V[ii];
      (*sub)->Zx[ii] = (*sub)->Z[ii];
      if ((*sub)->Vx[ii] > 0 & (*sub)->Zx[ii] <= 0)
      {(*sub)->Zx[ii] = setting->dx * setting->dy;}
      // if waterfall condition is met, only use half of the volumes and areas
      if (data->surf[ii] < bath->bottomXP[ii])
      {
        (*sub)->Vx[ii] = (*sub)->Vxm[map->iPjc[ii]];
        (*sub)->Zx[ii] = 0.5 * (*sub)->Z[ii];//(*sub)->Zxm[map->iPjc[ii]];
      }
      if (data->surf[map->iPjc[ii]] < bath->bottomXP[ii])
      {
        (*sub)->Vx[ii] = (*sub)->Vxp[ii];
        (*sub)->Zx[ii] = 0.5 * (*sub)->Z[ii];//(*sub)->Zxp[ii];
      }
    }
    else
    {
      (*sub)->Vx[ii] = 2 * (*sub)->Vxp[ii];
      (*sub)->Zx[ii] = 2 * (*sub)->Z[ii];//(*sub)->Zxp[ii];
    }

    // for y direction
    if (map->icjP[ii] < setting->N2ci)
    {
      if ((*sub)->Vyp[ii] > 0 & (*sub)->Vym[map->icjP[ii]] > 0)
      {(*sub)->Vy[ii] = (*sub)->Vyp[ii] + (*sub)->Vym[map->icjP[ii]];}
      else
      {(*sub)->Vy[ii] = 2.0*((*sub)->Vyp[ii] + (*sub)->Vym[map->icjP[ii]]);}
      //(*sub)->Zy[ii] = (*sub)->Zyp[ii] + (*sub)->Zym[map->icjP[ii]];
      //(*sub)->Vy[ii] = (*sub)->V[ii];
      (*sub)->Zy[ii] = (*sub)->Z[ii];
      if ((*sub)->Vy[ii] > 0 & (*sub)->Zy[ii] <= 0)
      {(*sub)->Zy[ii] = setting->dx * setting->dy;}
      // if waterfall condition is met, only use half of the volumes and areas
      if (data->surf[ii] < bath->bottomYP[ii])
      {
        (*sub)->Vy[ii] = (*sub)->Vym[map->icjP[ii]];
        (*sub)->Zy[ii] = 0.5 * (*sub)->Z[ii];//(*sub)->Zym[map->icjP[ii]];
      }
      if (data->surf[map->icjP[ii]] < bath->bottomYP[ii])
      {
        (*sub)->Vy[ii] = (*sub)->Vyp[ii];
        (*sub)->Zy[ii] = 0.5 * (*sub)->Z[ii];//(*sub)->Zyp[ii];
      }
    }
    else
    {
      (*sub)->Vy[ii] = 2 * (*sub)->Vyp[ii];
      (*sub)->Zy[ii] = 2 * (*sub)->Z[ii];//(*sub)->Zyp[ii];
    }
      // ZhiLi20181004
      // adjust the face area with Omin and Nmin, added by ZhiLi20181004
      if (setting->useminA == 1)
      {
          if ((*sub)->Nx[ii] > (*sub)->Nmin[ii] & (*sub)->Nmin[ii] > 0)
          {(*sub)->Nx[ii] = (*sub)->Nx[ii] * ((*sub)->Nmin[ii]/(*sub)->Nx[ii]);}
          if ((*sub)->Oy[ii] > (*sub)->Omin[ii] & (*sub)->Omin[ii] > 0)
          {(*sub)->Oy[ii] = (*sub)->Oy[ii] * ((*sub)->Omin[ii]/(*sub)->Oy[ii]);}
      }
      // ZhiLi20181018
      if (setting->phiSurface == 1)
      {
          if ((*sub)->Vx[ii] > 0) {(*sub)->Vx[ii] = (*sub)->Vx[ii] * (*sub)->phiX[ii];}
          if ((*sub)->Vy[ii] > 0) {(*sub)->Vy[ii] = (*sub)->Vy[ii] * (*sub)->phiY[ii];}
      }
  }
  // ========== Reduce cell volume for topographical dissipation, ZHILI20180914 ==========
//    double beta;
//    beta = 0.2;
//    for (ii = 0; ii < setting->N2ci; ii++)
//    {
//        if ((*sub)->Nx[ii] == 0)
//        {
//            if ((*sub)->Oy[ii] == 0 & (*sub)->Oy[map->icjM[ii]] > 0)
//            {
//                (*sub)->Vx[map->iMjc[ii]] = beta * (*sub)->Vx[map->iMjc[ii]];
//                (*sub)->Vy[map->icjM[ii]] = beta * (*sub)->Vy[map->icjM[ii]];
//            }
//            else if ((*sub)->Oy[ii] > 0 & (*sub)->Oy[map->icjM[ii]] == 0)
//            {
//                (*sub)->Vx[map->iMjc[ii]] = beta * (*sub)->Vx[map->iMjc[ii]];
//                (*sub)->Vy[ii] = beta * (*sub)->Vy[ii];
//            }
//        }
//        else if ((*sub)->Nx[map->iMjc[ii]] == 0)
//        {
//            if ((*sub)->Oy[ii] == 0 & (*sub)->Oy[map->icjM[ii]] > 0)
//            {
//                (*sub)->Vx[ii] = beta * (*sub)->Vx[ii];
//                (*sub)->Vy[map->icjM[ii]] = beta * (*sub)->Vy[map->icjM[ii]];
//            }
//            else if ((*sub)->Oy[ii] > 0 & (*sub)->Oy[map->icjM[ii]] == 0)
//            {
//                (*sub)->Vx[ii] = beta * (*sub)->Vx[ii];
//                (*sub)->Vy[ii] = beta * (*sub)->Vy[ii];
//            }
//        }
//    }
  // ============ REMOVE THE STAGGERED VOLUME FOR DEBUGGING, ZHILI20170722 ============
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    /*(*sub)->Vx[ii] = (*sub)->V[ii];
    (*sub)->Vy[ii] = (*sub)->V[ii];
    (*sub)->Zx[ii] = (*sub)->Z[ii];
    (*sub)->Zy[ii] = (*sub)->Z[ii];*/
    // remove small areas and volumes
    if ((*sub)->Vx[ii] < Vmin) {(*sub)->Vx[ii] = 0.0;}
    if ((*sub)->Vy[ii] < Vmin) {(*sub)->Vy[ii] = 0.0;}
    if ((*sub)->V[ii] < Vmin) {(*sub)->V[ii] = 0.0;}
    if ((*sub)->Z[ii] < Azmin) {(*sub)->Z[ii] = 1.0;}
    if ((*sub)->Zx[ii] < Amin) {(*sub)->Zx[ii] = 0.0;}
    if ((*sub)->Zy[ii] < Amin) {(*sub)->Zy[ii] = 0.0;}
    if ((*sub)->Nx[ii] < Amin) {(*sub)->Nx[ii] = 0.0;}
    if ((*sub)->Oy[ii] < Amin) {(*sub)->Oy[ii] = 0.0;}
    // if ((*sub)->N[ii] < Amin) {(*sub)->N[ii] = 0.0;}
    // if ((*sub)->O[ii] < Amin) {(*sub)->O[ii] = 0.0;}
  }
  // adjust for boundaries
  for (ii = 0; ii < setting->ny; ii++)
  {
    (*sub)->Vx[map->iMgt[ii]] = (*sub)->Vx[map->iMbd[ii]];
    (*sub)->Vy[map->iMgt[ii]] = (*sub)->Vy[map->iMbd[ii]];
    (*sub)->Zx[map->iMgt[ii]] = (*sub)->Zx[map->iMbd[ii]];
    (*sub)->Zy[map->iMgt[ii]] = (*sub)->Zy[map->iMbd[ii]];
    (*sub)->Nx[map->iMgt[ii]] = (*sub)->Nx[map->iMbd[ii]];
    (*sub)->Oy[map->iMgt[ii]] = (*sub)->Oy[map->iMbd[ii]];
    (*sub)->V[map->iMgt[ii]] = (*sub)->V[map->iMbd[ii]];
    // (*sub)->N[map->iMgt[ii]] = (*sub)->N[map->iMbd[ii]];
    // (*sub)->O[map->iMgt[ii]] = (*sub)->O[map->iMbd[ii]];
    (*sub)->Z[map->iMgt[ii]] = (*sub)->Z[map->iMbd[ii]];
    (*sub)->Vx[map->iPgt[ii]] = (*sub)->Vx[map->iPbd[ii]];
    (*sub)->Vy[map->iPgt[ii]] = (*sub)->Vy[map->iPbd[ii]];
    (*sub)->Zx[map->iPgt[ii]] = (*sub)->Zx[map->iPbd[ii]];
    (*sub)->Zy[map->iPgt[ii]] = (*sub)->Zy[map->iPbd[ii]];
    (*sub)->Nx[map->iPgt[ii]] = (*sub)->Nx[map->iPbd[ii]];
    (*sub)->Oy[map->iPgt[ii]] = (*sub)->Oy[map->iPbd[ii]];
    (*sub)->V[map->iPgt[ii]] = (*sub)->V[map->iPbd[ii]];
    // (*sub)->N[map->iPgt[ii]] = (*sub)->N[map->iPbd[ii]];
    // (*sub)->O[map->iPgt[ii]] = (*sub)->O[map->iPbd[ii]];
    (*sub)->Z[map->iPgt[ii]] = (*sub)->Z[map->iPbd[ii]];
    if (setting->useSubDrag == 2)
    {(*sub)->Yh[map->iPgt[ii]] = (*sub)->Yh[map->iPbd[ii]];}
    if (setting->useSubDrag == 3)
    {
      (*sub)->CvX[map->iPgt[ii]] = (*sub)->CvX[map->iPbd[ii]];
      (*sub)->CvY[map->iPgt[ii]] = (*sub)->CvY[map->iPbd[ii]];
      (*sub)->CvX[map->iMgt[ii]] = (*sub)->CvX[map->iMbd[ii]];
      (*sub)->CvY[map->iMgt[ii]] = (*sub)->CvY[map->iMbd[ii]];
      (*sub)->CdsX[map->iPgt[ii]] = (*sub)->CdsX[map->iPbd[ii]];
      (*sub)->CdsY[map->iPgt[ii]] = (*sub)->CdsY[map->iPbd[ii]];
    }
    if (setting->useSubDrag == 4)
    {(*sub)->Cn[map->iPgt[ii]] = (*sub)->Cn[map->iPbd[ii]];}
    // (*sub)->redCD[map->iPgt[ii]] = (*sub)->redCD[map->iPbd[ii]];
    // following FREHD?
    (*sub)->Nx[map->iMgt[ii]] = 0;
    (*sub)->Nx[map->iPgt[ii]] = 0;
    (*sub)->Nx[map->iPbd[ii]] = 0;
  }
  for (ii = 0; ii < setting->nx; ii++)
  {
    (*sub)->Vx[map->jMgt[ii]] = (*sub)->Vx[map->jMbd[ii]];
    (*sub)->Vy[map->jMgt[ii]] = (*sub)->Vy[map->jMbd[ii]];
    (*sub)->Zx[map->jMgt[ii]] = (*sub)->Zx[map->jMbd[ii]];
    (*sub)->Zy[map->jMgt[ii]] = (*sub)->Zy[map->jMbd[ii]];
    (*sub)->Nx[map->jMgt[ii]] = (*sub)->Nx[map->jMbd[ii]];
    (*sub)->Oy[map->jMgt[ii]] = (*sub)->Oy[map->jMbd[ii]];
    (*sub)->V[map->jMgt[ii]] = (*sub)->V[map->jMbd[ii]];
    // (*sub)->N[map->jMgt[ii]] = (*sub)->N[map->jMbd[ii]];
    // (*sub)->O[map->jMgt[ii]] = (*sub)->O[map->jMbd[ii]];
    (*sub)->Z[map->jMgt[ii]] = (*sub)->Z[map->jMbd[ii]];
    (*sub)->Vx[map->jPgt[ii]] = (*sub)->Vx[map->jPbd[ii]];
    (*sub)->Vy[map->jPgt[ii]] = (*sub)->Vy[map->jPbd[ii]];
    (*sub)->Zx[map->jPgt[ii]] = (*sub)->Zx[map->jPbd[ii]];
    (*sub)->Zy[map->jPgt[ii]] = (*sub)->Zy[map->jPbd[ii]];
    (*sub)->Nx[map->jPgt[ii]] = (*sub)->Nx[map->jPbd[ii]];
    (*sub)->Oy[map->jPgt[ii]] = (*sub)->Oy[map->jPbd[ii]];
    (*sub)->V[map->jPgt[ii]] = (*sub)->V[map->jPbd[ii]];
    // (*sub)->N[map->jPgt[ii]] = (*sub)->N[map->jPbd[ii]];
    // (*sub)->O[map->jPgt[ii]] = (*sub)->O[map->jPbd[ii]];
    (*sub)->Z[map->jPgt[ii]] = (*sub)->Z[map->jPbd[ii]];
    if (setting->useSubDrag == 2)
    {(*sub)->Yh[map->jPgt[ii]] = (*sub)->Yh[map->jPbd[ii]];}
    if (setting->useSubDrag == 3)
    {
      (*sub)->CvX[map->jPgt[ii]] = (*sub)->CvX[map->jPbd[ii]];
      (*sub)->CvY[map->jPgt[ii]] = (*sub)->CvY[map->jPbd[ii]];
      (*sub)->CvX[map->jMgt[ii]] = (*sub)->CvX[map->jMbd[ii]];
      (*sub)->CvY[map->jMgt[ii]] = (*sub)->CvY[map->jMbd[ii]];
      (*sub)->CdsX[map->jPgt[ii]] = (*sub)->CdsX[map->jPbd[ii]];
      (*sub)->CdsY[map->jPgt[ii]] = (*sub)->CdsY[map->jPbd[ii]];
    }
    if (setting->useSubDrag == 4)
    {(*sub)->Cn[map->jPgt[ii]] = (*sub)->Cn[map->jPbd[ii]];}
    // (*sub)->redCD[map->jPgt[ii]] = (*sub)->redCD[map->jPbd[ii]];
    // following FREHD?
    (*sub)->Oy[map->jMgt[ii]] = 0;
    //if (irank == 0 & setting->bcType == 1) {(*sub)->Oy[map->jMgt[ii]] = 0;}
    //if (irank == nrank-1 & setting->bcType == 2) {(*sub)->Oy[map->jPgt[ii]] = 0;}

  }
    // zero the face area along YP face except where tidal BC is added, ZhiLi20180617
    if (irank == nrank - 1 & setting->bcType != 2)
    {
        int isTideBC[setting->nx], jj;
        for (ii = 0; ii < setting->nx; ii++) {isTideBC[ii] = 0;}
        for (ii = 0; ii < setting->tideLocLengthP; ii++)
        {jj = setting->tideLocP[ii] - setting->nx*(setting->ny-1); isTideBC[jj] = 1;}
        for (ii = 0; ii < setting->nx; ii++)
        {
            if (isTideBC[ii] == 0)
            {
                (*sub)->Oy[map->jPbd[ii]] = 0;
                (*sub)->Oy[map->jPgt[ii]] = 0;
            }
        }
    }
}
// ==================== Update subgrid face depth ====================
void updateAllFaceDepthSub(Data **data, Sub *sub, Maps *map, Config *setting, int irank, int nrank)
{
  int ii;
  for (ii = 0; ii < setting->N2ct; ii++)
  {
    (*data)->depthXP[ii] = sub->Nx[ii] / setting->dy;
    (*data)->depthYP[ii] = sub->Oy[ii] / setting->dx;
  }
  // zero out the depth at y outer boundaries of the domain
  if (irank == nrank - 1)
  {
    for (ii = 0; ii < setting->nx; ii++)
    {
      (*data)->depthYP[map->jPbd[ii]] = 0;
      (*data)->depthYP[map->jPgt[ii]] = 0;
    }
  }
  if (irank == 0)
  {
    for (ii = 0; ii < setting->nx; ii++)
    {(*data)->depthYP[map->jMgt[ii]] = 0;}
  }
  // remove negative depth
  for (ii = 0; ii < setting->N2ct; ii++)
  {
    if ((*data)->depthXP[ii] < 0) {(*data)->depthXP[ii] = 0;}
    if ((*data)->depthYP[ii] < 0) {(*data)->depthYP[ii] = 0;}
  }

}
// ==================== Compute the advection term ====================
// only the 1st order upwind scheme is implemented so far
void advectionTermSub(Data **data, Maps *map, Sub *sub, Config *setting)
{
    int ii;
    double advX, advY;
    double CFLh = 0.7, CFLl = 0.5, CFLx, CFLy;
    for (ii = 0; ii < setting->N2ci; ii++)
    {
        (*data)->advX[ii] = (0.5/(setting->dx*sub->phiX[ii]))*((*data)->uuXP[ii] + fabs((*data)->uuXP[ii])) * \
        ((*data)->uuXP[ii] - (*data)->uuXP[map->iMjc[ii]]) + \
        (0.5/(setting->dx*sub->phiX[ii]))*((*data)->uuXP[ii] - fabs((*data)->uuXP[ii])) * \
        ((*data)->uuXP[map->iPjc[ii]] - (*data)->uuXP[ii]) + \
        (0.5/(setting->dy*sub->phiY[ii]))*((*data)->vvXP[ii] + fabs((*data)->vvXP[ii])) * \
        ((*data)->uuXP[ii] - (*data)->uuXP[map->icjM[ii]]) + \
        (0.5/(setting->dy*sub->phiY[ii]))*((*data)->vvXP[ii] - fabs((*data)->vvXP[ii])) * \
        ((*data)->uuXP[map->icjP[ii]] - (*data)->uuXP[ii]);

        (*data)->advY[ii] = (0.5/(setting->dx*sub->phiX[ii]))*((*data)->uuYP[ii] + fabs((*data)->uuYP[ii])) * \
        ((*data)->vvYP[ii] - (*data)->vvYP[map->iMjc[ii]]) + \
        (0.5/(setting->dx*sub->phiX[ii]))*((*data)->uuYP[ii] - fabs((*data)->uuYP[ii])) * \
        ((*data)->vvYP[map->iPjc[ii]] - (*data)->vvYP[ii]) + \
        (0.5/(setting->dy*sub->phiY[ii]))*((*data)->vvYP[ii] + fabs((*data)->vvYP[ii])) * \
        ((*data)->vvYP[ii] - (*data)->vvYP[map->icjM[ii]]) + \
        (0.5/(setting->dy*sub->phiY[ii]))*((*data)->vvYP[ii] - fabs((*data)->vvYP[ii])) * \
        ((*data)->vvYP[map->icjP[ii]] - (*data)->vvYP[ii]);
        if ((*data)->uuXP[ii] == 0) {(*data)->advX[ii] = 0;}
        if ((*data)->vvYP[ii] == 0) {(*data)->advY[ii] = 0;}
    }
    // Apply a CFL limiter, implemented by ZhiLi 20170922
    for (ii = 0; ii < setting->N2ci; ii++)
    {
        // adjust the nonlinear terms based on the CFL number
        if ((*data)->CFLx[ii] > setting->CFLl & (*data)->CFLx[ii] <= setting->CFLh)
        {
            (*data)->advX[ii] = (*data)->advX[ii] * (setting->CFLh - (*data)->CFLx[ii]) / \
            (setting->CFLh - setting->CFLl);
            (*data)->advX[map->iPjc[ii]] = (*data)->advX[map->iPjc[ii]] * \
            (setting->CFLh - (*data)->CFLx[map->iPjc[ii]]) / (setting->CFLh - setting->CFLl);
            (*data)->advX[map->iMjc[ii]] = (*data)->advX[map->iMjc[ii]] * \
            (setting->CFLh - (*data)->CFLx[map->iMjc[ii]]) / (setting->CFLh - setting->CFLl);
            (*data)->advX[map->icjP[ii]] = (*data)->advX[map->icjP[ii]] * \
            (setting->CFLh - (*data)->CFLx[map->icjP[ii]]) / (setting->CFLh - setting->CFLl);
            (*data)->advX[map->icjM[ii]] = (*data)->advX[map->icjM[ii]] * \
            (setting->CFLh - (*data)->CFLx[map->icjM[ii]]) / (setting->CFLh - setting->CFLl);
        }
        else if ((*data)->CFLx[ii] > setting->CFLh)
        {
            (*data)->advX[ii] = 0;
            (*data)->advX[map->iPjc[ii]] = 0;
            (*data)->advX[map->iMjc[ii]] = 0;
            (*data)->advX[map->icjP[ii]] = 0;
            (*data)->advX[map->icjM[ii]] = 0;
        }
        if ((*data)->CFLy[ii] > setting->CFLl & (*data)->CFLy[ii] < setting->CFLh)
        {
            (*data)->advY[ii] = (*data)->advY[ii] * (setting->CFLh - (*data)->CFLy[ii]) / \
            (setting->CFLh - setting->CFLl);
            (*data)->advY[map->iPjc[ii]] = (*data)->advY[map->iPjc[ii]] * \
            (setting->CFLh - (*data)->CFLy[map->iPjc[ii]]) / (setting->CFLh - setting->CFLl);
            (*data)->advY[map->iMjc[ii]] = (*data)->advY[map->iMjc[ii]] * \
            (setting->CFLh - (*data)->CFLy[map->iMjc[ii]]) / (setting->CFLh - setting->CFLl);
            (*data)->advY[map->icjP[ii]] = (*data)->advY[map->icjP[ii]] * \
            (setting->CFLh - (*data)->CFLy[map->icjP[ii]]) / (setting->CFLh - setting->CFLl);
            (*data)->advY[map->icjM[ii]] = (*data)->advY[map->icjM[ii]] * \
            (setting->CFLh - (*data)->CFLy[map->icjM[ii]]) / (setting->CFLh - setting->CFLl);
        }
        else if ((*data)->CFLx[ii] > setting->CFLh)
        {
            (*data)->advY[ii] = 0;
            (*data)->advY[map->iPjc[ii]] = 0;
            (*data)->advY[map->iMjc[ii]] = 0;
            (*data)->advY[map->icjP[ii]] = 0;
            (*data)->advY[map->icjM[ii]] = 0;
        }

        (*data)->EnXP[ii] = (*data)->EnXP[ii] - setting->dt * (*data)->advX[ii];
        (*data)->EnYP[ii] = (*data)->EnYP[ii] - setting->dt * (*data)->advY[ii];
    }
//  int ii;
//  double advX, advY;
//  double CFLh = 0.7, CFLl = 0.5, CFLx, CFLy;
//  for (ii = 0; ii < setting->N2ci; ii++)
//  {
//    advX = (0.5/setting->dx)*((*data)->uuXP[ii] + fabs((*data)->uuXP[ii])) * \
//      ((*data)->uuXP[ii] - (*data)->uuXP[map->iMjc[ii]]) * sub->Nx[ii] + \
//      (0.5/setting->dx)*((*data)->uuXP[ii] - fabs((*data)->uuXP[ii])) * \
//      ((*data)->uuXP[map->iPjc[ii]] - (*data)->uuXP[ii]) * sub->Nx[ii] + \
//      (0.5/setting->dy)*((*data)->vvXP[ii] + fabs((*data)->vvXP[ii])) * \
//      ((*data)->uuXP[ii] - (*data)->uuXP[map->icjM[ii]]) * sub->Oy[map->icjM[ii]] + \
//      (0.5/setting->dy)*((*data)->vvXP[ii] - fabs((*data)->vvXP[ii])) * \
//      ((*data)->uuXP[map->icjP[ii]] - (*data)->uuXP[ii]) * sub->Oy[ii];
//    advY = (0.5/setting->dx)*((*data)->uuYP[ii] + fabs((*data)->uuYP[ii])) * \
//      ((*data)->vvYP[ii] - (*data)->vvYP[map->iMjc[ii]]) * sub->Nx[map->iMjc[ii]] + \
//      (0.5/setting->dx)*((*data)->uuYP[ii] - fabs((*data)->uuYP[ii])) * \
//      ((*data)->vvYP[map->iPjc[ii]] - (*data)->vvYP[ii]) * sub->Nx[ii] + \
//      (0.5/setting->dy)*((*data)->vvYP[ii] + fabs((*data)->vvYP[ii])) * \
//      ((*data)->vvYP[ii] - (*data)->vvYP[map->icjM[ii]]) * sub->Oy[ii] + \
//      (0.5/setting->dy)*((*data)->vvYP[ii] - fabs((*data)->vvYP[ii])) * \
//      ((*data)->vvYP[map->icjP[ii]] - (*data)->vvYP[ii]) * sub->Oy[ii];
//    // A limiter to remove unexisted flow due to interpolation, ZhiLi 20170927
//    if ((*data)->uuXP[ii] == 0) {advX = 0;}
//    if ((*data)->vvYP[ii] == 0) {advY = 0;}
//    // A CFL limiter implemented by Zhi Li 20170922
//    CFLx = (*data)->uuXP[ii] * setting->dt / setting->dx;
//    CFLy = (*data)->vvYP[ii] * setting->dt / setting->dy;
//    if (CFLx > CFLl & CFLx < CFLh)
//    { advX = advX * (CFLh - CFLx) / (CFLh - CFLl);}
//    else if (CFLx >= CFLh)
//    { advX = 0;}
//    if (CFLy > CFLl & CFLy < CFLh)
//    { advY = advY * (CFLh - CFLy) / (CFLh - CFLl);}
//    else if (CFLy >= CFLh)
//    { advY = 0;}
//
//    if (sub->Vx[ii] != 0)
//    {(*data)->EnXP[ii] = (*data)->EnXP[ii] - setting->dt * advX / sub->Vx[ii];}
//    if (sub->Vy[ii] != 0)
//    {(*data)->EnYP[ii] = (*data)->EnYP[ii] - setting->dt * advY / sub->Vy[ii];}
//  }
}

// ==================== Compute the diffusion term ====================
void diffusionTermSub(Data **data, Maps *map, Sub *sub, Config *setting)
{
  int ii;
  double diffX, diffY;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    if (sub->Vx[ii] > 0)
    {
      diffX = (setting->NUx / sub->Vx[ii]) * \
        (sub->Nx[map->iPjc[ii]]*((*data)->uuXP[map->iPjc[ii]]-(*data)->uuXP[ii])/setting->dx - \
        sub->Nx[ii]*(-(*data)->uuXP[map->iMjc[ii]]+(*data)->uuXP[ii])/setting->dx) + \
        (setting->NUy / sub->Vx[ii]) * \
        (sub->Oy[ii]*((*data)->uuXP[map->icjP[ii]]-(*data)->uuXP[ii])/setting->dy - \
        sub->Oy[map->icjM[ii]]*(-(*data)->uuXP[map->icjM[ii]]+(*data)->uuXP[ii])/setting->dy);
    }
    else
    {diffX = 0;}
    if (sub->Vy[ii] > 0)
    {
      diffY = (setting->NUx / sub->Vy[ii]) * \
        (sub->Nx[ii]*((*data)->vvYP[map->iPjc[ii]]-(*data)->vvYP[ii])/setting->dx - \
        sub->Nx[map->iMjc[ii]]*(-(*data)->vvYP[map->iMjc[ii]]+(*data)->vvYP[ii])/setting->dx) + \
        (setting->NUy / sub->Vy[ii]) * \
        (sub->Oy[map->icjP[ii]]*((*data)->vvYP[map->icjP[ii]]-(*data)->vvYP[ii])/setting->dy - \
        sub->Oy[ii]*(-(*data)->vvYP[map->icjM[ii]]+(*data)->vvYP[ii])/setting->dy);
    }
    else
    {diffY = 0;}
    (*data)->EnXP[ii] = (*data)->EnXP[ii] + setting->dt * diffX;
    (*data)->EnYP[ii] = (*data)->EnYP[ii] + setting->dt * diffY;
  }
}

// ==================== Transverse momentum, ZhiLi20181014 ====================
//void transMomTermSub(Data **data, Maps *map, Sub *sub, Config *setting)
//{
//    int ii;
//    double transX, transY;
//    for (ii = 0; ii < setting->N2ci; ii++)
//    {
//        if (sub->Vx[ii] > 0)
//        {
//            transX = 0.5*((*data)->uuXP[map->iMjc[ii]] + fabs((*data)->uuXP[map->iMjc[ii]])) * \
//            sub->xiXM[ii] * sub->Nx[ii] / sub->Vx[ii] + \
//            0.5*((*data)->uuXP[map->iPjc[ii]] - fabs((*data)->uuXP[map->iPjc[ii]])) * \
//            sub->xiXP[ii] * sub->Nx[ii] / sub->Vx[ii];
//        }
//        else
//        {transX = 0;}
//        if (sub->Vy[ii] > 0)
//        {
//            transY = 0.5*((*data)->vvYP[map->icjM[ii]] + fabs((*data)->vvYP[map->icjM[ii]])) * \
//            sub->xiYM[ii] * sub->Oy[ii] / sub->Vy[ii] + \
//            0.5*((*data)->vvYP[map->icjP[ii]] + fabs((*data)->vvYP[map->icjP[ii]])) * \
//            sub->xiYP[ii] * sub->Oy[ii] / sub->Vy[ii];
//        }
//        else
//        {transY = 0;}
//        (*data)->EnXP[ii] = (*data)->EnXP[ii] + setting->dt * setting->dt * setting->g * transX;
//        (*data)->EnYP[ii] = (*data)->EnYP[ii] + setting->dt * setting->dt * setting->g * transY;
//    }
//}

// ==================== Compute the wind term ====================
void windTermSub(Data **data, BC *bc, Sub *sub, Config *setting, int tt)
{
  double phi, omega, tau, Cw, pi = 3.1415926, rho = 1000.0, tauXP, tauYP;
  int ii;
  if (setting->useWind == 1)
  {
    // phi is the wind direction from the north
    phi = bc->winddir[tt] + setting->northAngle;
    // omega is the wind direction in rad from positive x
    omega = phi * pi / 180.0;
    // apply wind term uniformly on each grid cell
    for (ii = 0; ii < setting->N2ci; ii++)
    {
      // tau is the total wind stress
      tau = setting->rhoa * setting->Cw * \
        (bc->windspd[tt] - (*data)->uuXP[ii]*cos(omega) - (*data)->vvYP[ii]*sin(omega)) * \
        (bc->windspd[tt] - (*data)->uuXP[ii]*cos(omega) - (*data)->vvYP[ii]*sin(omega));
      // apply the thin layer model when necessary
      if ((*data)->depthXP[ii] < setting->hD & setting->useThinLayer == 1)
      {tauXP = tau * exp(setting->CwT*((*data)->depthXP[ii]-setting->hD)/setting->hD);}
      else
      {tauXP = tau;}
      if ((*data)->depthYP[ii] < setting->hD & setting->useThinLayer == 1)
      {tauYP = tau * exp(setting->CwT*((*data)->depthYP[ii]-setting->hD)/setting->hD);}
      else
      {tauYP = tau;}
      // add wind drag to the source term
      if (sub->Vx[ii] > 0 && sub->Nx[ii] > 0)
      {
        (*data)->EnXP[ii] = (*data)->EnXP[ii] + setting->dt * tauXP * cos(omega) * sub->Zx[ii] / \
          (sub->Vx[ii] * rho);
      }
      if (sub->Vy[ii] > 0 && sub->Oy[ii] > 0)
      {
        (*data)->EnYP[ii] = (*data)->EnYP[ii] - setting->dt * tauYP * sin(omega) * sub->Zy[ii] / \
          (sub->Vy[ii] * rho);
      }
    }
  }
}

// ==================== Compute the drag term ====================
void dragTermSub(Data **data, Maps *map, Sub *sub, Config *setting)
{
  int ii;
  double velx, vely, effHx, effHy;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    velx = sqrt((*data)->uuXP[ii]*(*data)->uuXP[ii] + \
      (*data)->vvXP[ii]*(*data)->vvXP[ii]);
    vely = sqrt((*data)->uuYP[ii]*(*data)->uuYP[ii] + \
      (*data)->vvYP[ii]*(*data)->vvYP[ii]);
    // modified for non-staggered wetting and drying, ZhiLi20180222
    // in x direction
    if (setting->staggeredV == 1)
    {
      if (sub->Vx[ii] > 0) {effHx = sub->Zx[ii] / sub->Vx[ii];}
      else {effHx = 0;}
      if (sub->Vy[ii] > 0) {effHy = sub->Zy[ii] / sub->Vy[ii];}
      else {effHy = 0;}
      // use V instead of Vx, Vy, ZHILI20180916
//        if (sub->V[ii] > 0) {effHx = sub->Zx[ii] / sub->V[ii];}
//        else {effHx = 0;}
//        if (sub->V[ii] > 0) {effHy = sub->Zy[ii] / sub->V[ii];}
//        else {effHy = 0;}
    }
    else
    {
      if (sub->V[ii] > 0 & sub->Vx[ii] > 0) {effHx = sub->Zx[ii] / sub->Vx[ii];}
      else if ((*data)->depthXP[ii] > 0) {effHx = 1.0 / (*data)->depthXP[ii];}
      else {effHx = 0;}
      if (sub->V[ii] > 0 & sub->Vy[ii] > 0) {effHy = sub->Zy[ii] / sub->Vy[ii];}
      else if ((*data)->depthYP[ii] > 0) {effHy = 1.0 / (*data)->depthYP[ii];}
      else {effHy = 0;}
    }
    (*data)->DragXP[ii] = 0.5 * setting->dt * (*data)->CDXP[ii] * velx * effHx;
    (*data)->DragYP[ii] = 0.5 * setting->dt * (*data)->CDYP[ii] * vely * effHy;

    /*if (sub->V[ii] > 0 & sub->Vx[ii] > 0)
    {
      (*data)->DragXP[ii] = 0.5 * setting->dt * (*data)->CDXP[ii] * \
      velx * sub->Zx[ii] / sub->Vx[ii];
    }
    else if ((*data)->depthXP[ii] > 0)
    {
      (*data)->DragXP[ii] = 0.5 * setting->dt * (*data)->CDXP[ii] * \
      velx / (*data)->depthXP[ii];;
    }
    else
    {
      (*data)->DragXP[ii] = 0;
    }
    // in y direction
    if (sub->V[ii] > 0 & sub->Vy[ii] > 0)
    {
      (*data)->DragYP[ii] = 0.5 * setting->dt * (*data)->CDYP[ii] * \
        vely * sub->Zy[ii] / sub->Vy[ii];
    }
    else if ((*data)->depthYP[ii] > 0)
    {
      (*data)->DragYP[ii] = 0.5 * setting->dt * (*data)->CDYP[ii] * \
        vely / (*data)->depthYP[ii];
    }
    else
    {
      (*data)->DragYP[ii] = 0;
    }*/
  }
  for (ii = 0; ii < setting->nx; ii++)
  {(*data)->DragXP[map->jMgt[ii]] = (*data)->DragXP[map->jMbd[ii]];}
  for (ii = 0; ii < setting->ny; ii++)
  {(*data)->DragYP[map->iMgt[ii]] = (*data)->DragYP[map->iMbd[ii]];}
}

// ============== Compute the RHS of the matrix equation ===============
void matrixSourceTermSub(Data **data, Maps *map, BC *bc, Sub *sub, Config *setting, int tt, int irank)
{
  int ii, jj, loc[setting->inflowLocLength];
  double p1, p2, p3, p4;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    // the non-staggered volume approach, ZhiLi20180222
    if (sub->Z[map->trps[ii]] > 0 | setting->staggeredV == 1)
    {
      (*data)->z[ii] = (*data)->surf[map->trps[ii]] * sub->Z[map->trps[ii]] - \
         setting->dt * \
        ((*data)->EnXP[map->trps[ii]] * sub->Nx[map->trps[ii]] - \
         (*data)->EnXP[map->iMjc[map->trps[ii]]] * sub->Nx[map->iMjc[map->trps[ii]]]) - \
         setting->dt *\
        ((*data)->EnYP[map->trps[ii]] * sub->Oy[map->trps[ii]] - \
         (*data)->EnYP[map->icjM[map->trps[ii]]] * sub->Oy[map->icjM[map->trps[ii]]]);
    }
    else
    {
      (*data)->z[ii] = (*data)->surf[map->trps[ii]] * setting->dx * setting->dy - \
         setting->dt * \
        ((*data)->EnXP[map->trps[ii]] * (*data)->depthXP[map->trps[ii]] * setting->dy - \
         (*data)->EnXP[map->iMjc[map->trps[ii]]] * (*data)->depthXP[map->iMjc[map->trps[ii]]]) * setting->dy - \
         setting->dt *\
        ((*data)->EnYP[map->trps[ii]] * (*data)->depthYP[map->trps[ii]] * setting->dx - \
         (*data)->EnYP[map->icjM[map->trps[ii]]] * (*data)->depthYP[map->icjM[map->trps[ii]]]) * setting->dx;
    }
  }

  if (irank == 0 & setting->bcType != 3)
  {
    for (ii = 0; ii < setting->inflowLocLength; ii++)
    {
      loc[ii] = (setting->inflowLoc[ii]%setting->nx)*setting->ny + \
        floor(setting->inflowLoc[ii]/setting->nx);
      (*data)->z[loc[ii]] = (*data)->z[loc[ii]] + \
        setting->dt * (bc->inflow[tt] / setting->inflowLocLength);
    }
  }
  // to make sure the matrix is non-singular
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    if ((*data)->z[ii] == 0 | (*data)->z[ii] == (*data)->surf[map->trps[ii]]*setting->dx*setting->dy)
    {(*data)->z[ii] = (*data)->surf[map->trps[ii]]*setting->dx*setting->dy;}
  }
}

// =============== Compute the coefficients of the matrix A ===============
void matrixCoeffSub(Data **data, Maps *map, Sub *sub, Config *setting)
{
  int ii,jj,kk;
  double coefx, coefy, area, Gmin = 0.00000001;
  coefx = setting->g * (setting->dt*setting->dt);
  coefy = setting->g * (setting->dt*setting->dt);
  area = setting->dx * setting->dy;
  if (setting->useSubgrid == 1)
  {
    for (ii = 0; ii < setting->N2ci; ii++)
    {
      //note : Theoretically we should apply the drag inversion when computing
      // these matrix coefficients, but it caused inconsistent flux. ZhiLi20170521

      jj = map->iMjc[ii];
      kk = map->icjM[ii];
      (*data)->GnXP[ii] = 0;
      (*data)->GnXM[ii] = 0;
      (*data)->GnYP[ii] = 0;
      (*data)->GnYM[ii] = 0;
      /*if (sub->Vx[ii] > 0)
      {(*data)->GnXP[ii] = coefx * sub->Nx[ii] * sub->Nx[ii] / sub->Vx[ii] / (1 + (*data)->DragXP[ii]);;}
      if (sub->Vx[jj] > 0)
      {(*data)->GnXM[ii] = coefx * sub->Nx[jj] * sub->Nx[jj] / sub->Vx[jj] / (1 + (*data)->DragXP[jj]);;}
      if (sub->Vy[ii] > 0)
      {(*data)->GnYP[ii] = coefy * sub->Oy[ii] * sub->Oy[ii] / sub->Vy[ii] / (1 + (*data)->DragYP[ii]);;}
      if (sub->Vy[kk] > 0)
      {(*data)->GnYM[ii] = coefy * sub->Oy[kk] * sub->Oy[kk] / sub->Vy[kk] / (1 + (*data)->DragYP[kk]);;}
      (*data)->GnCt[ii] = sub->Z[ii] + (*data)->GnXP[ii] + (*data)->GnYP[ii] + (*data)->GnXM[ii] + \
        (*data)->GnYM[ii];*/
        // the non-staggered volume approach, ZhiLi20180222
      if (sub->V[ii] > 0 | setting->staggeredV == 1)
      {
        if (sub->Vx[ii] > 0)
        {(*data)->GnXP[ii] = coefx * sub->Nx[ii] * sub->Nx[ii] / sub->Vx[ii] / (1 + (*data)->DragXP[ii]);}
        if (sub->Vx[jj] > 0)
        {(*data)->GnXM[ii] = coefx * sub->Nx[jj] * sub->Nx[jj] / sub->Vx[jj] / (1 + (*data)->DragXP[jj]);}
        if (sub->Vy[ii] > 0)
        {(*data)->GnYP[ii] = coefy * sub->Oy[ii] * sub->Oy[ii] / sub->Vy[ii] / (1 + (*data)->DragYP[ii]);}
        if (sub->Vy[kk] > 0)
        {(*data)->GnYM[ii] = coefy * sub->Oy[kk] * sub->Oy[kk] / sub->Vy[kk] / (1 + (*data)->DragYP[kk]);}
        (*data)->GnCt[ii] = sub->Z[ii] + (*data)->GnXP[ii] + (*data)->GnYP[ii] + (*data)->GnXM[ii] + \
          (*data)->GnYM[ii];
      }
      else
      {
        (*data)->GnXP[ii] = coefx * area * (*data)->depthXP[ii] / (1 + (*data)->DragXP[ii]);
        (*data)->GnXM[ii] = coefx * area * (*data)->depthXP[map->iMjc[ii]] / \
          (1 + (*data)->DragXP[map->iMjc[ii]]);
        (*data)->GnYP[ii] = coefy * area * (*data)->depthYP[ii] / (1 + (*data)->DragYP[ii]);
        (*data)->GnYM[ii] = coefy * area * (*data)->depthYP[map->icjM[ii]] / \
          (1 + (*data)->DragYP[map->icjM[ii]]);
        (*data)->GnCt[ii] = area + (*data)->GnXP[ii] + \
          (*data)->GnYP[ii] + (*data)->GnXM[ii] + (*data)->GnYM[ii];
      }
    }
  }
  else if (setting->useSubgrid == 2)
  {
    for (ii = 0; ii < setting->N2ci; ii++)
    {
      jj = map->iMjc[ii];
      kk = map->icjM[ii];
      (*data)->GnXP[ii] = 0;
      (*data)->GnXM[ii] = 0;
      (*data)->GnYP[ii] = 0;
      (*data)->GnYM[ii] = 0;
      if (sub->Vx[ii] > 0)
      {(*data)->GnXP[ii] = coefx * sub->Nx[ii] / setting->dx / (1 + (*data)->DragXP[ii]);}
      if (sub->Vx[jj] > 0)
      {(*data)->GnXM[ii] = coefx * sub->Nx[jj] / setting->dx / \
        (1 + (*data)->DragXP[map->iMjc[ii]]);}
      if (sub->Vy[ii] > 0)
      {(*data)->GnYP[ii] = coefy * sub->Oy[ii] / setting->dx / (1 + (*data)->DragYP[ii]);}
      if (sub->Vy[kk] > 0)
      {(*data)->GnYM[ii] = coefy * sub->Oy[jj] / setting->dx / \
        (1 + (*data)->DragYP[map->icjM[ii]]);}
      (*data)->GnCt[ii] = sub->Z[ii] + (*data)->GnXP[ii] + \
        (*data)->GnYP[ii] + (*data)->GnXM[ii] + (*data)->GnYM[ii];
    }
  }
  // avoid non-singularity of the matrix
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    if ((*data)->GnCt[ii] == 0.0 | (*data)->GnCt[ii] == setting->dx*setting->dy)
    {(*data)->GnCt[ii] = setting->dx*setting->dy;}//1.0;}// (*data)->z[ii] = (*data)->surf[map->trps[ii]];}
    // avoid inundation of dry cells when face velocities are zero, ZhiLi20171002
    if ((*data)->depth[ii] == 0)
    {
      if ((*data)->uuXP[ii] == 0 & (*data)->uuXP[map->iMjc[ii]] == 0 & (*data)->vvYP[ii] == 0 & \
          (*data)->vvYP[map->icjM[ii]] == 0)
      {
          (*data)->GnCt[ii] = setting->dx*setting->dy;//1.0;
          (*data)->GnXM[ii] = 0.0;
          (*data)->GnXP[ii] = 0.0;
          (*data)->GnYM[ii] = 0.0;
          (*data)->GnYP[ii] = 0.0;
          (*data)->z[map->sprt[ii]] = (*data)->surf[ii]*setting->dx*setting->dy;
      }
    }
  }
}

// ==================== Update velocities ====================
void updateVelocitySub(Data **data, Maps *map, Sub *sub, Config *setting)
{
  int ii, count = 0;
  double effHx, effHy;
  if (setting->useSubgrid == 1)
  {
    for (ii = 0; ii < setting->N2ci; ii++)
    {
      (*data)->uuXP[ii] = 0;
      (*data)->vvYP[ii] = 0;
      if (setting->staggeredV == 1)
      {
        if (sub->Vx[ii] > 0)
        {
          effHx = sub->Nx[ii] / sub->Vx[ii];
//          // remove too small effH, ZhiLi20180415
//          if (setting->useAVcutoff == 1 & effHx < setting->effHmin)
//          {effHx = setting->effHmin;}
        }
        else {effHx = 0;}
        if (sub->Vy[ii] > 0)
        {
          effHy = sub->Oy[ii] / sub->Vy[ii];
//          if (setting->useAVcutoff == 1 & effHy < setting->effHmin)
//          {effHy = setting->effHmin;}
        }
        else {effHy = 0;}
      }
      else
      {
        if (((*data)->depthXP[ii] > 0) & (sub->V[ii] == 0 | sub->V[map->iPjc[ii]] == 0))
        {effHx = 1.0 / setting->dx;}
        else if (sub->Vx[ii] > 0)
        {effHx = sub->Nx[ii] / sub->Vx[ii];}
        else
        {effHx = 0;}
        if (((*data)->depthYP[ii] > 0) & (sub->V[ii] == 0 | sub->V[map->icjP[ii]] == 0))
        {effHy = 1.0 / setting->dy;}
        else if (sub->Vy[ii] > 0)
        {effHy = sub->Oy[ii] / sub->Vy[ii];}
        else
        {effHy = 0;}
      }
      (*data)->uuXP[ii] = (*data)->EnXP[ii] - (setting->g * setting->dt * effHx) \
        * ((*data)->surf[map->iPjc[ii]] - (*data)->surf[ii]);// / (1.0 + (*data)->DragXP[ii]) ;
      (*data)->vvYP[ii] = (*data)->EnYP[ii] - (setting->g * setting->dt * effHy) \
        * ((*data)->surf[map->icjP[ii]] - (*data)->surf[ii]);// / (1.0 + (*data)->DragYP[ii]) ;
    }
  }
  // if using subgrid for continuity only
  else if (setting->useSubgrid == 2)
  {
    for (ii = 0; ii < setting->N2ci; ii++)
    {
      (*data)->uuXP[ii] = (*data)->EnXP[ii] - (setting->g * setting->dt / setting->dx) * \
        ((*data)->surf[map->iPjc[ii]] - (*data)->surf[ii]);// / (1 + (*data)->DragXP[ii]);
      (*data)->vvYP[ii] = (*data)->EnYP[ii] - (setting->g * setting->dt / setting->dy) * \
        ((*data)->surf[map->icjP[ii]] - (*data)->surf[ii]);// / (1 + (*data)->DragYP[ii]);
    }
  }
  // remove velocity on dry face
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    if (sub->Nx[ii] == 0)
    {(*data)->uuXP[ii] = 0;}
    if (sub->Oy[ii] == 0)
    {(*data)->vvYP[ii] = 0;}
  }
  // remove velocity out of dry cells, ZhiLi 20170925
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    if ((*data)->depth[ii] == 0)
    {
      if ((*data)->uuXP[ii] > 0) {(*data)->uuXP[ii] = 0;}
      if ((*data)->uuXP[map->iMjc[ii]] < 0) {(*data)->uuXP[map->iMjc[ii]] = 0;}
      if ((*data)->vvYP[ii] > 0) {(*data)->vvYP[ii] = 0;}
      if ((*data)->vvYP[map->icjM[ii]] < 0) {(*data)->vvYP[map->icjM[ii]] = 0;}
    }
  }
  // apply the CFL limiter
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    if ((*data)->rmvCFL[ii] == 1)
    {
      (*data)->uuXP[ii] = 0;
      (*data)->vvYP[ii] = 0;
      (*data)->uuXP[map->iMjc[ii]] = 0;
      (*data)->vvYP[map->icjM[ii]] = 0;
      (*data)->rmvCFL[ii] = 0;
      count += 1;
    }
  }
  // if (count > 0)
  // {printf("Number of cells removed by CFL limiter = %d\n",count); count = 0;}
  // apply another CFL limiter
  /*for (ii = 0; ii < setting->N2ci; ii++)
  {
    if ((*data)->uuXP[ii] > 0.99*setting->dxf/setting->dt)
    {(*data)->uuXP[ii] = 0.99*setting->dxf/setting->dt;}
    if ((*data)->vvYP[ii] > 0.99*setting->dyf/setting->dt)
    {(*data)->vvYP[ii] = 0.99*setting->dyf/setting->dt;}
  }*/
}

// =============== Compute face flow rates ===============
void computeFaceFlowRateSub(Data **data, Maps *map, Sub *sub, Config *setting)
{
  int ii;
  for (ii = 0; ii < setting->N2ct; ii++)
  {
    /*if ((*data)->wtfYP[ii] > 0)
    {wtfAy = sub->Op[ii];}
    else if ((*data)->wtfYP[ii] < 0)
    {wtfAy = sub->Om[map->icjP[ii]];}
    else
    {wtfAy = sub->Oy[ii];}

    if ((*data)->wtfXP[ii] > 0)
    {wtfAx = sub->Np[ii];}
    else if ((*data)->wtfXP[ii] < 0)
    {wtfAx = sub->Nm[map->iPjc[ii]];}
    else
    {wtfAx = sub->Nx[ii];}*/

    (*data)->Fuu[ii] = (*data)->uuXP[ii] * sub->Nx[ii];
    (*data)->Fvv[ii] = (*data)->vvYP[ii] * sub->Oy[ii];
  }
}

// =============== compute cell volumes ===============
void computeVolumeSub(Data **data, Sub *sub, Config *setting)
{
  int ii;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    (*data)->cellV[ii] = sub->V[ii];
    if (setting->useScalar == 1)
    {(*data)->Sm[ii] = (*data)->S[ii] * (*data)->cellV[ii];}
  }
}

// =============== compute cell volume using flux ===============
void volumeByFluxSub(Data **data, Sub *sub, Maps *map, BC *bc, Config *setting, int tt, int irank)
{
  int ii, jj, ind1, ind2;
  double nflux;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    (*data)->cellV[ii] = (*data)->cellV[ii] + setting->dt * \
      ((*data)->Fuu[map->iMjc[ii]] - (*data)->Fuu[ii] + \
      (*data)->Fvv[map->icjM[ii]] - (*data)->Fvv[ii]);
    if (irank == 0 & setting->bcType != 3)
    {
      for (jj = 0; jj < setting->inflowLocLength; jj++)
      {if (ii == setting->inflowLoc[jj])
        {(*data)->cellV[ii] = (*data)->cellV[ii] + bc->inflow[tt]*setting->dt;}}
    }
    if (setting->useEvap == 1)
    {
      (*data)->cellV[ii] = (*data)->cellV[ii] - bc->evap[tt] * \
        setting->dt * sub->Z[ii];
    }
    if (setting->useRain == 1)
    {
      (*data)->cellV[ii] = (*data)->cellV[ii] + bc->rain[tt] * \
        setting->dt * sub->Z[ii];
    }
    // added by ZhiLi 20180416
    if (setting->checkConservation == 1)
    {(*data)->Vloss[ii] = sub->V[ii] - (*data)->cellV[ii];}
    if ((*data)->cellV[ii] < 0)
    {
      ind2 = floor(ii / setting->nx);
      ind1 = ii - ind2 * setting->nx;
      //printf("WARNING: Net flux exceeds cell volume at (%d,%d)...\n",ind1,ind2);
      (*data)->cellV[ii] = 0;
    }
  }
}

// =============== interpolate velocity at cell faces ===============
void velocityInterpSub(Data **data, Sub *sub, Config *setting)
{
  int ii;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    if (sub->Nx[ii] + sub->Nx[ii+setting->nx] != 0)
    {
      (*data)->uuYP[ii] = (0.25 * ((*data)->Fuu[ii] + (*data)->Fuu[ii-1] + \
        (*data)->Fuu[ii+setting->nx] + (*data)->Fuu[ii-1+setting->nx])) / \
        (0.5 * (sub->Nx[ii] + sub->Nx[ii+setting->nx]));
    }
    else
    {(*data)->uuYP[ii] = 0;}
    if (sub->Oy[ii] + sub->Oy[ii+1] != 0)
    {
      (*data)->vvXP[ii] = 0.25 * ((*data)->Fvv[ii] + (*data)->Fvv[ii+1] + \
        (*data)->Fvv[ii-setting->nx] + (*data)->Fvv[ii+1-setting->nx]) / \
        (0.5 * (sub->Oy[ii] + sub->Oy[ii+1]));
    }
    else
    {(*data)->vvXP[ii] = 0;}
  }
}

// ===================== Adjust tidal boundary velocity ====================
void adjustTidalVelocitySub(Data **data, Maps *map, Sub *sub, Config *setting, int irank, int nrank)
{
  int ii, jj;
  double allFlux = 0, Fww = 0;
  if (irank == nrank - 1 & setting->bcType != 2)
  {
    for (ii = 0; ii < setting->tideLocLengthP; ii++)
    {
      jj = setting->tideLocP[ii];
      // vertical flux due to change in tidal bc
      Fww = ((*data)->surf[jj] - (*data)->surfOld[jj]) * sub->Z[jj] / setting->dt;
      // flux on the yp face from continuity
      allFlux = (*data)->Fvv[map->icjM[jj]] + (*data)->Fuu[map->iMjc[jj]]  - \
        (*data)->Fuu[jj] - Fww;
      // velocity is flux / area
      if (allFlux != 0 & sub->Oy[jj] > 0)
      {(*data)->vvYP[jj] = allFlux / sub->Oy[jj];}
      else
      {(*data)->vvYP[jj] = 0;}
      Fww = 0;
      allFlux = 0;
    }
  }
  if (irank == 0 & setting->bcType != 1)
  {
    for (ii = 0; ii < setting->tideLocLengthM; ii++)
    {
      jj = setting->tideLocM[ii];
      // vertical flux due to change in tidal bc
      Fww = ((*data)->surf[jj] - (*data)->surfOld[jj]) * sub->Z[jj] / setting->dt;
      // flux on the yp face from continuity
      allFlux = -(*data)->Fvv[jj] + (*data)->Fuu[map->iMjc[jj]]  - \
        (*data)->Fuu[jj] - Fww;
      // velocity is flux / area
      if (allFlux != 0 & sub->Oy[jj] > 0)
      {(*data)->vvYP[map->icjM[jj]] = -allFlux / sub->Oy[jj];}
      else
      {(*data)->vvYP[map->icjM[jj]] = 0;}
      Fww = 0;
      allFlux = 0;
    }
  }
}

void detectWaterfallLocationSub(Data **data, Sub *sub, Bath *bath, Maps *map, Config *setting)
{
  int ii, count = 0;
  double surf, bot, area, dep, dmin = 0.0;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    (*data)->wtfXP[ii] = 0;
    (*data)->wtfYP[ii] = 0;
  }
  // find waterfall locations
  // waterfall on XP face
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    surf = (*data)->surf[ii];
    bot = bath->bottomXP[ii];
    area = sub->Nx[ii];
    dep = (*data)->surf[map->iPjc[ii]] - bot;
    if (surf < bot-dmin & dep >= setting->wtfh)// & (bot-surf) > Ck*dep)
    {(*data)->wtfXP[ii] = -1;   count += 1;}
    //else if (surf < bot-dmin & dep < setting->wtfh)// & dep > 0 & (bot-surf) > Ck*dep)
    //{(*data)->wtfXP[ii] = -2;}
    //if (ii == 7234)
    //{printf("surf,bot,botjP,botYP = %lf,%lf,%lf,%lf\n",surf,bot,dep,(*data)->wtfXP[ii]);}
  }
  // waterfall on XM face
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    surf = (*data)->surf[ii];
    bot = bath->bottomXP[map->iMjc[ii]];
    area = sub->Nx[map->iMjc[ii]];
    dep = (*data)->surf[map->iMjc[ii]] - bot;
    if (surf < bot-dmin & dep >= setting->wtfh)// & (bot-surf) > Ck*dep)
    {(*data)->wtfXP[map->iMjc[ii]] = 1;   count += 1;}
    //else if (surf < bot-dmin & dep < setting->wtfh)// & dep > 0 & (bot-surf) > Ck*dep)
    //{(*data)->wtfXP[map->iMjc[ii]] = 2;}
  }
  // waterfall on YP face
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    surf = (*data)->surf[ii];
    bot = bath->bottomYP[ii];
    area = sub->Oy[ii];
    dep = (*data)->surf[map->icjP[ii]] - bot;
    if (surf < bot-dmin & dep >= setting->wtfh)// & (bot-surf) > Ck*dep)
    {(*data)->wtfYP[ii] = -1;   count += 1;}
    //else if (surf < bot-dmin & dep < setting->wtfh)// & dep > 0 & (bot-surf) > Ck*dep)
    //{(*data)->wtfYP[ii] = -2;}
    //if (ii == 7268+150)
    //{printf("surf,bot,botjP,botYP = %lf,%lf,%lf,%lf\n",surf,bot,dep,(*data)->wtfYP[ii]);}
  }
  // waterfall on YM face
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    surf = (*data)->surf[ii];
    bot = bath->bottomYP[map->icjM[ii]];
    area = sub->Oy[map->icjM[ii]];
    dep = (*data)->surf[map->icjM[ii]] - bot;
    if (surf < bot-dmin & dep >= setting->wtfh)// & (bot-surf) > Ck*dep)
    {(*data)->wtfYP[map->icjM[ii]] = 1;   count += 1;}
    /*if (ii == 10365)
    {
      printf("wtf,surf,surfM,bot,dep,avgdep,N = %f,%f,%f,%f,%f,%f,%f\n",\
      (*data)->wtfXP[ii-1],surf,(*data)->surf[ii-1],bot,dep,(*data)->depthXP[ii-1],area);
    }*/
    //else if (surf < bot-dmin & dep < setting->wtfh)// & dep > 0 & (bot-surf) > Ck*dep)
    //{(*data)->wtfYP[map->icjM[ii]] = 2;}
    //if (ii == 7268+150)
    //{printf("surf,bot,dep,botYM = %lf,%lf,%lf,%lf\n",surf,bot,dep,(*data)->wtfYP[map->icjM[ii]]);}
  }
  // if (count > 0)
  // {printf("Number of cell faces with waterfall = %d\n",count);  count = 0;}
  /*for (ii = 0; ii < setting->N2ci; ii++)
  {
    if ((*data)->wtfXP[ii] == 2)
    {
      (*data)->surf[ii] = bath->bottomZ[ii];
      (*data)->cellV[ii] = 0;
      (*data)->wtfXP[ii] = 0;
    }
    else if ((*data)->wtfXP[ii] == -2)
    {
      (*data)->surf[map->iPjc[ii]] = bath->bottomZ[map->iPjc[ii]];
      (*data)->cellV[map->iPjc[ii]] = 0;
      (*data)->wtfXP[ii] = 0;
    }
    if ((*data)->wtfYP[ii] == 2)
    {
      (*data)->surf[ii] = bath->bottomZ[ii];
      (*data)->cellV[ii] = 0;
      (*data)->wtfYP[ii] = 0;
    }
    else if ((*data)->wtfYP[ii] == -2)
    {
      (*data)->surf[map->icjP[ii]] = bath->bottomZ[map->icjP[ii]];
      (*data)->cellV[map->icjP[ii]] = 0;
      (*data)->wtfYP[ii] = 0;
    }
  }*/
}


// ===================== The water fall model ====================
void waterfallCorrectionSub(Data **data, Sub *sub, Bath *bath, Maps *map, Config *setting)
{
  int ii;
  double Cw = 0.7, CFLmax = 0.24/setting->subR, wtfA, vel0;
  double allwtf = 0, netflux = 0, diff = 0;
  // correct waterfall velocities
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    // compute waterfall velocity
    // This method has bug! Because the subgrid face area at wtf locations is
    // always 0, so no wtf velocity exists! ZhiLi20170827
    if ((*data)->wtfXP[ii] != 0)
    {
      vel0 = (*data)->uuXP[ii];
      // if ((*data)->wtfXP[ii] > 0)
      // {wtfA = sub->Np[ii];}
      // else
      // {wtfA = sub->Nm[map->iPjc[ii]];}

      (*data)->uuXP[ii] = (*data)->wtfXP[ii] * \
      Cw * sqrt(2 * setting->g * sub->Nx[ii] / setting->dy);
      //(*data)->uuXP[ii] = (*data)->wtfXP[ii] * \
      Cw * sqrt(2 * setting->g * wtfA / setting->dy);
      if (fabs((*data)->uuXP[ii]) > CFLmax * setting->dx / setting->dt)
      {(*data)->uuXP[ii] = (*data)->wtfXP[ii] * CFLmax * setting->dx / setting->dt;}
      (*data)->Vloss4[0] -= (vel0 - (*data)->uuXP[ii]) * setting->dt * sub->Nx[ii];
    }
    if ((*data)->wtfYP[ii] != 0)
    {
      vel0 = (*data)->vvYP[ii];
      // if ((*data)->wtfYP[ii] > 0)
      // {wtfA = sub->Op[ii];}
      // else
      // {wtfA = sub->Om[map->icjP[ii]];}

      (*data)->vvYP[ii] = (*data)->wtfYP[ii] * \
      Cw * sqrt(2 * setting->g * sub->Oy[ii] / setting->dx);
      //(*data)->vvYP[ii] = (*data)->wtfYP[ii] * \
      Cw * sqrt(2 * setting->g * wtfA / setting->dx);
      if (fabs((*data)->vvYP[ii]) > CFLmax * setting->dy / setting->dt)
      {(*data)->vvYP[ii] = (*data)->wtfYP[ii] * CFLmax * setting->dy / setting->dt;}
      (*data)->Vloss4[0] -= (vel0 - (*data)->vvYP[ii]) * setting->dt * sub->Oy[ii];
      /*if (ii == 879)
      {
        int mm = map->icjM[ii];
        //printf("faceAm, wtfYPm, vvYM, depthYM= %lf,%lf,%lf,%lf\n",sub->Op[mm],(*data)->wtfYP[mm],(*data)->vvYP[mm],(*data)->depthYP[mm]);
        printf("faceA, wtfYP, vvYP, depth= %lf,%lf,%lf,%lf\n",sub->Oy[ii],(*data)->wtfYP[ii],(*data)->vvYP[ii],(*data)->depth[ii]);
      }*/
    }
  }
  // check if too much flow out of a cell
  /*for (ii = 0; ii < setting->N2ci; ii++)
  {
    // total flow rates out of a cell
    netflux = (*data)->uuXP[map->iMjc[ii]] * sub->Nx[map->iMjc[ii]] + \
      (*data)->vvYP[map->icjM[ii]] * sub->Oy[map->icjM[ii]] - \
      (*data)->uuXP[ii] * sub->Nx[ii] - (*data)->vvYP[ii] * sub->Oy[ii];
    if (netflux < -(*data)->cellV[ii])
    {
      allwtf = (*data)->wtfXP[ii] - (*data)->wtfXP[map->iMjc[ii]] + \
        (*data)->wtfYP[ii] - (*data)->wtfYP[map->icjM[ii]];
      if (allwtf != 0)
      {
        diff = (-netflux - (*data)->cellV[ii]) / allwtf;
        if ((*data)->wtfXP[ii] == 1)
        {(*data)->uuXP[ii] = (*data)->uuXP[ii] - (diff/sub->Nx[ii]);}
        if ((*data)->wtfXP[map->iMjc[ii]] == -1)
        {(*data)->uuXP[map->iMjc[ii]] = (*data)->uuXP[map->iMjc[ii]] + (diff/sub->Nx[map->iMjc[ii]]);}
        if ((*data)->wtfYP[ii] == 1)
        {(*data)->vvYP[ii] = (*data)->vvYP[ii] - (diff/sub->Oy[ii]);}
        if ((*data)->wtfYP[map->icjM[ii]] == -1)
        {(*data)->vvYP[map->icjM[ii]] = (*data)->vvYP[map->icjM[ii]] - (diff/sub->Oy[map->icjM[ii]]);}
      }
    }
    allwtf = 0;
    diff = 0;
    netflux = 0;
  }*/
  // reset waterfall locations
  /*for (ii = 0; ii < setting->N2ci; ii++)
  {
    (*data)->wtfXP[ii] = 0;
    (*data)->wtfYP[ii] = 0;
  }*/
}

// =============== Update the bottom drag coefficient ===============
void updateCDSub(Data **data, Sub *sub, Maps *map, Config *setting)
{
  // the von Karman constant
  double K = 0.41, h, n;
  int ii, jj;
  double phi = 50, dNx, dNy;
  double dN = 0.1;
  double coeff, coeffX, coeffY;
  // update CD only when the CD is depth dependent
  if (setting->CDnotN == 0 | setting->CDnotN == 2)
  {
    // adjust the bottom CD by curvature
    if (setting->useSubDrag == 1)
    {
      for (ii = 0; ii < setting->N2ci; ii++)
      {
        if ((*data)->depth[ii] > setting->z0)
        {
          if (phi * sub->CvX[ii] / (setting->subR*setting->subR) > dN)
          {dNx = dN;}
          else
          {dNx = phi * sub->CvX[ii] / (setting->subR*setting->subR);}
          if (phi * sub->CvY[ii] / (setting->subR*setting->subR) > dN)
          {dNy = dN;}
          else
          {dNy = phi * sub->CvY[ii] / (setting->subR*setting->subR);}
          coeffX = setting->manningN + dNx;
          coeffY = setting->manningN + dNy;
          (*data)->CDXP[ii] = setting->g*coeffX*coeffX/pow((*data)->depth[ii], 1.0/3.0);
          (*data)->CDYP[ii] = setting->g*coeffY*coeffY/pow((*data)->depth[ii], 1.0/3.0);
        }
      }
    }
    // adjust the bottom CD by Volp2013
    else if (setting->useSubDrag == 2)
    {
      coeff = setting->g*setting->manningN*setting->manningN;
      for (ii = 0; ii < setting->N2ci; ii++)
      {
        // No anisotropy???!!!
        if (sub->Yh[ii] > 0)
        {
          //(*data)->CDXP[ii] = coeff*sub->Yh[ii] * sub->CvX[ii] * phi;
          //(*data)->CDYP[ii] = coeff*sub->Yh[ii] * sub->CvY[ii] * phi;
          (*data)->CDXP[ii] = coeff*sub->Yh[ii];
          (*data)->CDYP[ii] = coeff*sub->Yh[ii];
          if ((*data)->CDXP[ii] > setting->CDmax)
          {(*data)->CDXP[ii] = setting->CDmax;}
          if ((*data)->CDYP[ii] > setting->CDmax)
          {(*data)->CDYP[ii] = setting->CDmax;}
        }
        else
        {
          (*data)->CDXP[ii] = 0;
          (*data)->CDYP[ii] = 0;
        }
      }
    }
    // the subgrid drag model by ZhiLi20180329
    else if (setting->useSubDrag == 3)
    {
      int im, ip;
      coeff = setting->g * setting->manningN * setting->manningN;
      for (ii = 0; ii < setting->N2ci; ii++)
      {
        if (sub->V[ii] > 0)
        {
          h = sub->V[ii] / (setting->dx*setting->dy);
//          h = (*data)->depth[ii];
          // in x direction
        (*data)->CDXP[ii] = coeff / pow(h, 1.0/3.0) + \
          setting->lambda1 * sub->CdsX[ii] - \
          setting->lambda2 * sub->CvX[ii];
          // in y direction
        (*data)->CDYP[ii] = coeff / pow(h, 1.0/3.0) + \
          setting->lambda1 * sub->CdsY[ii] + \
          setting->lambda2 * sub->CvY[ii];

          // Combine Li with Casas, ZhiLi20180914
//            coeff = setting->g * sub->Cn[ii] * sub->Cn[ii];
//            (*data)->CDXP[ii] = setting->lambda2*coeff/pow(h, 1.0/3.0) + setting->lambda1*sub->CdsX[ii];
//            (*data)->CDYP[ii] = setting->lambda2*coeff/pow(h, 1.0/3.0) + setting->lambda1*sub->CdsY[ii];

          // remove negative drag coefficients
          if ((*data)->CDXP[ii] < 0) {(*data)->CDXP[ii] = 0;}
          if ((*data)->CDYP[ii] < 0) {(*data)->CDYP[ii] = 0;}

        }
      }
    }
    // the subgrid drag model by Casas2010, ZhiLi20180506
    else if (setting->useSubDrag == 4)
    {
        for (ii = 0; ii < setting->N2ci; ii++)
        {
            if (sub->V[ii] > 0)
            {
                h = sub->V[ii] / (setting->dx*setting->dy);
                (*data)->CDXP[ii] = setting->g * sub->Cn[ii] * sub->Cn[ii] / pow(h, 1.0/3.0);
                (*data)->CDYP[ii] = setting->g * sub->Cn[ii] * sub->Cn[ii] / pow(h, 1.0/3.0);
            }
        }
    }
    // the subgrid model based on roughness height, ZhiLi20181028
    else if (setting->useSubDrag == 5)
    {
      for (ii = 0; ii < setting->N2ci; ii++)
      {
        if (sub->V[ii] > 0)
        {
          h = sub->V[ii] / (setting->dx*setting->dy);
          n = 0.06 * pow(sub->Z0[ii], 1.0/7.0);
          if (n > setting->manningN)
          {coeff = setting->g * n * n;}
          else
          {coeff = setting->g * setting->manningN * setting->manningN;}
          (*data)->CDXP[ii] = coeff / pow(h, 1.0/3.0);
          (*data)->CDYP[ii] = coeff / pow(h, 1.0/3.0);
        }
      }
    }
    else
    {
      coeff = setting->g * setting->manningN * setting->manningN;
      for (ii = 0; ii < setting->N2ci; ii++)
      {
        if (sub->V[ii] > 0)
        {
          h = sub->V[ii] / (setting->dx*setting->dy);
          (*data)->CDXP[ii] = coeff / pow(h, 1.0/3.0);
          (*data)->CDYP[ii] = coeff / pow(h, 1.0/3.0);
            // Reduce drag at channel corners, ZhiLi20180911
//            double CdN = -0.002;
//            if ((sub->Nx[ii] == 0 & sub->Oy[ii] == 0) | (sub->Nx[ii] == 0 & sub->Oy[map->icjM[ii]] == 0) | (sub->Nx[map->iMjc[ii]] == 0 & sub->Oy[ii] == 0) | (sub->Nx[map->iMjc[ii]] == 0 & sub->Oy[map->icjM[ii]] == 0))
//            {
//                (*data)->CDXP[ii] = CdN; (*data)->CDYP[ii] = CdN;
//                (*data)->CDXP[map->iMjc[ii]] = CdN; (*data)->CDYP[map->icjM[ii]] = CdN;
//            }
        }
      }
    }
  }
}

// ========== reduce bottom CD in narrow channels ==========
void reduceCDSub(Data **data, Sub *sub, Config *setting)
{
  int ii;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    (*data)->CDXP[ii] = (*data)->CDXP[ii] / (sub->redCD[ii] * sub->redCD[ii]);
    (*data)->CDYP[ii] = (*data)->CDYP[ii] / (sub->redCD[ii] * sub->redCD[ii]);
  }
}

// ========== Scalar diffusion ==========
void scalarDiffusionSub(Data **data, Maps *map, Sub *sub, Config *setting)
{
  int ii;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    (*data)->Sm[ii] = (*data)->Sm[ii] + setting->dt * \
      ((setting->Kx * sub->Nx[ii] / setting->dx) * \
      ((*data)->S[map->iPjc[ii]] - (*data)->S[ii]) - \
      (setting->Kx * sub->Nx[map->iMjc[ii]] / setting->dx) * \
      ((*data)->S[ii] - (*data)->S[map->iMjc[ii]]) + \
      (setting->Ky * sub->Oy[ii] / setting->dy) * \
      ((*data)->S[map->icjP[ii]] - (*data)->S[ii]) - \
      (setting->Ky * sub->Oy[map->icjM[ii]] / setting->dy) * \
      ((*data)->S[ii] - (*data)->S[map->icjM[ii]]));
  }
}

// ================ update scalar concentration ===============
void updateScalarSub(Data **data, BC *bc, Bath *bath, Maps *map, Sub *sub, Config *setting, int tt, int irank, int nrank)
{
  int ii, jj, *N, flag = 0;
  int iP = 0, iM = 0, jP = 0, jM = 0;
  double *Sarr, *Smax, *Smin;
  double dXPold, dXMold, dYPold, dYMold;
  // local scalar limiter
  N = malloc(setting->N2ci * sizeof(int));
  Smax = malloc(setting->N2ci * sizeof(double));
  Smin = malloc(setting->N2ci * sizeof(double));
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    N[ii] = 0, Smax[ii] = 0, Smin[ii] = 0, jj = 0;
    // check if the target cell is connected to its neighbors
    // comment out this section to disable the limiter, ZhiLi 20170607
    dXPold = (*data)->surfOld[map->iPjc[ii]] - bath->bottomZ[map->iPjc[ii]];
    //if (dXPold > 0) {N[ii] += 1; iP = 1; dXPold = 0;}
    if (sub->Nx[ii] > 0 & dXPold > 0) {N[ii] += 1; iP = 1; dXPold = 0;}

    dXMold = (*data)->surfOld[map->iMjc[ii]] - bath->bottomZ[map->iMjc[ii]];
    //if (dXMold > 0) {N[ii] += 1; iM = 1; dXMold = 0;}
    if (sub->Nx[map->iMjc[ii]] > 0 & dXMold > 0) {N[ii] += 1; iM = 1; dXMold = 0;}

    if (irank != nrank - 1)
    {
      dYPold = (*data)->surfOld[map->icjP[ii]] - bath->bottomZ[map->icjP[ii]];
      //if (dYPold > 0) {N[ii] += 1; jP = 1; dYPold = 0;}
      if (sub->Oy[ii] > 0 & dYPold > 0) {N[ii] += 1; jP = 1; dYPold = 0;}
    }
    else
    {
      if (map->icjP[ii] <= setting->N2ci)
      {
        dYPold = (*data)->surfOld[map->icjP[ii]] - bath->bottomZ[map->icjP[ii]];
        //if (dYPold > 0) {N[ii] += 1; jP = 1; dYPold = 0;}
        if (sub->Oy[ii] > 0 & dYPold > 0) {N[ii] += 1; jP = 1; dYPold = 0;}
      }
    }

    if (irank != 0)
    {
      dYMold = (*data)->surfOld[map->icjM[ii]] - bath->bottomZ[map->icjM[ii]];
      //if (dYMold > 0) {N[ii] += 1; jM = 1; dYMold = 0;}
      if (sub->Oy[map->icjM[ii]] > 0 & dYMold > 0) {N[ii] += 1; jM = 1; dYMold = 0;}
    }
    else
    {
      if (map->icjM[ii] <= setting->N2ci)
      {
        dYMold = (*data)->surfOld[map->icjM[ii]] - bath->bottomZ[map->icjM[ii]];
        //if (dYMold > 0) {N[ii] += 1; jM = 1; dYMold = 0;}
        if (sub->Oy[map->icjM[ii]] > 0 & dYMold > 0) {N[ii] += 1; jM = 1; dYMold = 0;}
      }
    }


    /*if ((*data)->depthXP[ii] > 0) {N[ii] += 1; iP = 1;}
    if ((*data)->depthXP[map->iMjc[ii]] > 0) {N[ii] += 1; iM = 1;}
    if ((*data)->depthYP[ii] > 0) {N[ii] += 1; jP = 1;}
    if ((*data)->depthYP[map->icjM[ii]] > 0) {N[ii] += 1; jM = 1;}*/
    // if yes, get the max and min values of its neighbors

    if (N[ii] != 0)
    {
      Sarr = malloc(N[ii] * sizeof(double));
      if (iP == 1) {Sarr[jj] = (*data)->S[map->iPjc[ii]]; jj++;}
      if (iM == 1) {Sarr[jj] = (*data)->S[map->iMjc[ii]]; jj++;}
      if (jP == 1) {Sarr[jj] = (*data)->S[map->icjP[ii]]; jj++;}
      if (jM == 1) {Sarr[jj] = (*data)->S[map->icjM[ii]]; jj++;}
      Smax[ii] = getMax(Sarr, N[ii]);
      Smin[ii] = getMin(Sarr, N[ii]);
      iP = 0;
      iM = 0;
      jP = 0;
      jM = 0;
      free(Sarr);
    }
  }
  // scalar from inflow
  if (irank == 0 & setting->bcType != 3)
  {
    for (ii = 0; ii < setting->inflowLocLength; ii++)
    {(*data)->Sm[setting->inflowLoc[ii]] = (*data)->Sm[setting->inflowLoc[ii]] + \
      setting->inflowS * bc->inflow[tt] * setting->dt / setting->inflowLocLength;}
  }
  // update scalar concentration
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    if ((*data)->cellV[ii] > 0)
    {(*data)->S[ii] = (*data)->Sm[ii] / (*data)->cellV[ii];}
    else
    {(*data)->S[ii] = 0;}
    // apply the limiter
    if (irank == 0 & setting->bcType != 3)
    {
      for (jj = 0; jj < setting->inflowLocLength; jj++)
      {if (ii == setting->inflowLoc[jj]) {flag = 1; break;}}
    }
    // Note: The scalar limiter was not compatible with MPI, so it was disabled
    // by ZhiLi 20170607
    // debugged to allow salinity transferred across jP and jM boundary
    // by ZhiLi 20170815

    if (N[ii] != 0 & (*data)->S[ii] > Smax[ii] & flag == 0)
    {
      /*if (map->icjP[ii] < setting->N2ci & map->icjM[ii] < setting->N2ci)
      {(*data)->S[ii] = Smax[ii];}
      else if (irank == nrank-1 & map->icjP[ii] >= setting->N2ci)
      {(*data)->S[ii] = Smax[ii];}
      else if (irank == 0 & map->icjM[ii] >= setting->N2ci)
      {(*data)->S[ii] = Smax[ii];}*/
      (*data)->S[ii] = Smax[ii];
    }
    else if (N[ii] != 0 & (*data)->S[ii] < Smin[ii] & flag == 0)
    {
      /*if (map->icjP[ii] < setting->N2ci & map->icjM[ii] < setting->N2ci)
      {(*data)->S[ii] = Smin[ii];}
      else if (irank == nrank-1 & map->icjP[ii] >= setting->N2ci)
      {(*data)->S[ii] = Smin[ii];}
      else if (irank == 0 & map->icjM[ii] >= setting->N2ci)
      {(*data)->S[ii] = Smin[ii];}*/
      (*data)->S[ii] = Smin[ii];
    }
    flag = 0;
  }

  free(N);
  free(Smax);
  free(Smin);
  // remove high scalar singularities
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    if ((*data)->S[ii] > 100 | (*data)->S[ii] < 0)
    {
      if (!(irank == 0 & map->icjM[ii] > setting->N2ci) & \
       !(irank == nrank - 1 & map->icjP[ii] > setting->N2ci))
      {
        (*data)->S[ii] = 0;
        (*data)->surf[ii] = (*data)->surf[ii] - (*data)->depth[ii];
        (*data)->depth[ii] = 0;
      }
    }
  }
}

// ========== extend the surface to the neighbor cell ==========
void extendSurface(Data *data, Sub *sub, double *eqsurf, int *j1, int *j2, int ii)
{
  *eqsurf = data->surf[ii];
  *j1 = sub->ind[ii];
  if (*j1 == 0)
  {*j2 = 0;}
  else if (*j1 == sub->Ne - 1)
  {*j2 = sub->Ne - 1;}
  else
  {
    if (*eqsurf > sub->allSurf[*j1])
    {*j2 = *j1 + 1;}
    else
    {*j2 = *j1 - 1;}
  }
}





// ----- read and initialize the subgrid dataset -----
void initSubArea(Sub **sub, Bath *bath, Config *setting, int irank)
{
  int ii;
  *sub = malloc(sizeof(Sub));
  if (setting->useSubgrid != 0)
  {
    // initialize subgrid variables for all ranks
    // AR stands for all ranks
    //if (irank == 0)
    {
      // Ne is the number of possible surface elevations
      (*sub)->Ne = round((setting->surfmax - setting->surfmin)/setting->dsurf);
      (*sub)->allSurf = malloc((*sub)->Ne*sizeof(double));
      for (ii = 0; ii < (*sub)->Ne; ii++)
      {(*sub)->allSurf[ii] = setting->surfmin + ii*setting->dsurf + bath->offset[0];}
      /*(*sub)->bathAR = malloc(setting->N2CI*sizeof(double));
      (*sub)->bathXPAR = malloc(setting->N2CI*sizeof(double));
      (*sub)->bathYPAR = malloc(setting->N2CI*sizeof(double));*/
      // (*sub)->wdNpAR = malloc(setting->N2CI*sizeof(double));
      // (*sub)->wdOpAR = malloc(setting->N2CI*sizeof(double));
      // (*sub)->wdNmAR = malloc(setting->N2CI*sizeof(double));
      // (*sub)->wdOmAR = malloc(setting->N2CI*sizeof(double));
      (*sub)->allVAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      // (*sub)->allNAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      // (*sub)->allOAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      (*sub)->allZAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      (*sub)->allVxpAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      (*sub)->allVypAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      (*sub)->allVxmAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      (*sub)->allVymAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      (*sub)->allNpAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      (*sub)->allOpAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      (*sub)->allNmAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      (*sub)->allOmAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      // (*sub)->allZxpAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      // (*sub)->allZypAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      // (*sub)->allZxmAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      // (*sub)->allZymAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      (*sub)->allCvXAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      (*sub)->allCvYAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      (*sub)->allYhAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      (*sub)->allCnAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      // (*sub)->allredCDAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      (*sub)->allCdsXAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
      (*sub)->allCdsYAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
    }
    // below are initialization of the current ranks
    // initialize new bottom and edge elevations
    /*(*sub)->bath = malloc(setting->N2ci*sizeof(double));
    (*sub)->bathXP = malloc(setting->N2ci*sizeof(double));
    (*sub)->bathYP = malloc(setting->N2ci*sizeof(double));*/
    // (*sub)->wdNp = malloc(setting->N2ci*sizeof(double));
    // (*sub)->wdOp = malloc(setting->N2ci*sizeof(double));
    // (*sub)->wdNm = malloc(setting->N2ci*sizeof(double));
    // (*sub)->wdOm = malloc(setting->N2ci*sizeof(double));
    // initialize all subgrid data
    (*sub)->allV = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    // (*sub)->allN = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    // (*sub)->allO = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    (*sub)->allZ = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    (*sub)->allVxp = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    (*sub)->allVyp = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    (*sub)->allVxm = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    (*sub)->allVym = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    (*sub)->allNp = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    (*sub)->allOp = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    (*sub)->allNm = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    (*sub)->allOm = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    // (*sub)->allZxp = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    // (*sub)->allZyp = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    // (*sub)->allZxm = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    // (*sub)->allZym = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    (*sub)->allCvX = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    (*sub)->allCvY = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    (*sub)->allYh = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    (*sub)->allCn = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    // (*sub)->allredCD = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    (*sub)->allCdsX = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    (*sub)->allCdsY = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    // initialize subgrid data for the current surface elevations
    (*sub)->V = malloc(setting->N2ct*sizeof(double));
    // (*sub)->N = malloc(setting->N2ct*sizeof(double));
    // (*sub)->O = malloc(setting->N2ct*sizeof(double));
    (*sub)->Z = malloc(setting->N2ct*sizeof(double));
    (*sub)->Vxp = malloc(setting->N2ct*sizeof(double));
    (*sub)->Vyp = malloc(setting->N2ct*sizeof(double));
    (*sub)->Vxm = malloc(setting->N2ct*sizeof(double));
    (*sub)->Vym = malloc(setting->N2ct*sizeof(double));
    (*sub)->Np = malloc(setting->N2ct*sizeof(double));
    (*sub)->Op = malloc(setting->N2ct*sizeof(double));
    (*sub)->Nm = malloc(setting->N2ct*sizeof(double));
    (*sub)->Om = malloc(setting->N2ct*sizeof(double));
    // (*sub)->Zxp = malloc(setting->N2ct*sizeof(double));
    // (*sub)->Zyp = malloc(setting->N2ct*sizeof(double));
    // (*sub)->Zxm = malloc(setting->N2ct*sizeof(double));
    // (*sub)->Zym = malloc(setting->N2ct*sizeof(double));
    (*sub)->Vx = malloc(setting->N2ct*sizeof(double));
    (*sub)->Vy = malloc(setting->N2ct*sizeof(double));
    (*sub)->Nx = malloc(setting->N2ct*sizeof(double));
    (*sub)->Oy = malloc(setting->N2ct*sizeof(double));
    (*sub)->Zx = malloc(setting->N2ct*sizeof(double));
    (*sub)->Zy = malloc(setting->N2ct*sizeof(double));
    (*sub)->CvX = malloc(setting->N2ct*sizeof(double));
    (*sub)->CvY = malloc(setting->N2ct*sizeof(double));
    (*sub)->Yh = malloc(setting->N2ct*sizeof(double));
    (*sub)->Cn = malloc(setting->N2ct*sizeof(double));
    // (*sub)->redCD = malloc(setting->N2ct*sizeof(double));
    (*sub)->CdsX = malloc(setting->N2ct*sizeof(double));
    (*sub)->CdsY = malloc(setting->N2ct*sizeof(double));
      (*sub)->allNx = malloc(setting->N2CI*sizeof(double));
      (*sub)->allOy = malloc(setting->N2CI*sizeof(double));
      (*sub)->allVx = malloc(setting->N2CI*sizeof(double));
      (*sub)->allVy = malloc(setting->N2CI*sizeof(double));
      (*sub)->allVo = malloc(setting->N2CI*sizeof(double));
      (*sub)->allZo = malloc(setting->N2CI*sizeof(double));
    // initialize the searching index
    (*sub)->ind = malloc(setting->N2ct*sizeof(int));
    for (ii = 0; ii < setting->N2ct; ii++)
    {(*sub)->ind[ii] = round((*sub)->Ne / 2);}
  }
}
void readSubArea(Sub **sub, Bath *bath, Config *setting, int irank)
{
    // read in data
    if (setting->useSubgrid != 0)
    {
      char filepath[100];
      strcpy(filepath, setting->inputFolder);
      strcat(filepath, setting->subgridFolder);
      char filename[100];
      // read bath
      /*strcpy(filename, filepath);
      strcat(filename, "bottom.dat");
      ReadFile((*sub)->bathAR, filename, setting->N2CI);
      memset(filename, 0, sizeof(filename));
      // read bath edges
      strcpy(filename, filepath);
      strcat(filename, "bottomXP.dat");
      ReadFile((*sub)->bathXPAR, filename, setting->N2CI);
      memset(filename, 0, sizeof(filename));
      strcpy(filename, filepath);
      strcat(filename, "bottomYP.dat");
      ReadFile((*sub)->bathYPAR, filename, setting->N2CI);
      memset(filename, 0, sizeof(filename));*/
      // read the bathymetry index
      /*strcpy(filename, filepath);
      strcat(filename, "botindX1.dat");
      ReadFile((*sub)->biX1AR, filename, setting->N2CI);
      memset(filename, 0, sizeof(filename));
      strcpy(filename, filepath);
      strcat(filename, "botindX2.dat");
      ReadFile((*sub)->biX2AR, filename, setting->N2CI);
      memset(filename, 0, sizeof(filename));
      strcpy(filename, filepath);
      strcat(filename, "botindY1.dat");
      ReadFile((*sub)->biY1AR, filename, setting->N2CI);
      memset(filename, 0, sizeof(filename));
      strcpy(filename, filepath);
      strcat(filename, "botindY2.dat");
      ReadFile((*sub)->biY2AR, filename, setting->N2CI);
      memset(filename, 0, sizeof(filename));
      strcpy(filename, filepath);
      strcat(filename, "Exflag.dat");
      ReadFile((*sub)->ExflagAR, filename, setting->N2CI);
      memset(filename, 0, sizeof(filename));
      strcpy(filename, filepath);
      strcat(filename, "Eyflag.dat");
      ReadFile((*sub)->EyflagAR, filename, setting->N2CI);
      memset(filename, 0, sizeof(filename));*/
      // read cutoff bottom elevation
      // strcpy(filename, filepath);
      // strcat(filename, "wdNp.dat");
      // ReadFile((*sub)->wdNpAR, filename, setting->N2CI);
      // memset(filename, 0, sizeof(filename));
      // strcpy(filename, filepath);
      // strcat(filename, "wdOp.dat");
      // ReadFile((*sub)->wdOpAR, filename, setting->N2CI);
      // memset(filename, 0, sizeof(filename));
      // strcpy(filename, filepath);
      // strcat(filename, "wdNm.dat");
      // ReadFile((*sub)->wdNmAR, filename, setting->N2CI);
      // memset(filename, 0, sizeof(filename));
      // strcpy(filename, filepath);
      // strcat(filename, "wdOm.dat");
      // ReadFile((*sub)->wdOmAR, filename, setting->N2CI);
      // memset(filename, 0, sizeof(filename));
      // read curvatures
      strcpy(filename, filepath);
      strcat(filename, "Yh.dat");
      ReadFile((*sub)->allYhAR, filename, setting->N2CI*(*sub)->Ne);
      memset(filename, 0, sizeof(filename));
      if (setting->useSubDrag == 3)
      {
          strcpy(filename, filepath);
          strcat(filename, "CvX.dat");
          ReadFile((*sub)->allCvXAR, filename, setting->N2CI*(*sub)->Ne);
          memset(filename, 0, sizeof(filename));
          strcpy(filename, filepath);
          strcat(filename, "CvY.dat");
          ReadFile((*sub)->allCvYAR, filename, setting->N2CI*(*sub)->Ne);
          memset(filename, 0, sizeof(filename));
        strcpy(filename, filepath);
        strcat(filename, "effCdX.dat");
        ReadFile((*sub)->allCdsXAR, filename, setting->N2CI*(*sub)->Ne);
        memset(filename, 0, sizeof(filename));
        strcpy(filename, filepath);
        strcat(filename, "effCdY.dat");
        ReadFile((*sub)->allCdsYAR, filename, setting->N2CI*(*sub)->Ne);
        memset(filename, 0, sizeof(filename));
      }
      if (setting->useSubDrag == 4 | setting->useSubDrag == 3)
      {
          strcpy(filename, filepath);
          strcat(filename, "Cn.dat");
          ReadFile((*sub)->allCnAR, filename, setting->N2CI*(*sub)->Ne);
          memset(filename, 0, sizeof(filename));
      }
      /*strcpy(filename, filepath);
      strcat(filename, "redCD.dat");
      ReadFile((*sub)->allredCDAR, filename, setting->N2CI*(*sub)->Ne);
      memset(filename, 0, sizeof(filename));*/
      // read cell center volumes and areas
      strcpy(filename, filepath);
      strcat(filename, "V.dat");
      ReadFile((*sub)->allVAR, filename, setting->N2CI*(*sub)->Ne);
      memset(filename, 0, sizeof(filename));
      strcpy(filename, filepath);
      strcat(filename, "Z.dat");
      ReadFile((*sub)->allZAR, filename, setting->N2CI*(*sub)->Ne);
      memset(filename, 0, sizeof(filename));
      // strcpy(filename, filepath);
      // strcat(filename, "N.dat");
      // ReadFile((*sub)->allNAR, filename, setting->N2CI*(*sub)->Ne);
      // memset(filename, 0, sizeof(filename));
      // strcpy(filename, filepath);
      // strcat(filename, "O.dat");
      // ReadFile((*sub)->allOAR, filename, setting->N2CI*(*sub)->Ne);
      // memset(filename, 0, sizeof(filename));
      // read cell face volumes and areas
      // volumes
      strcpy(filename, filepath);
      strcat(filename, "Vxp.dat");
      ReadFile((*sub)->allVxpAR, filename, setting->N2CI*(*sub)->Ne);
      memset(filename, 0, sizeof(filename));
      strcpy(filename, filepath);
      strcat(filename, "Vyp.dat");
      ReadFile((*sub)->allVypAR, filename, setting->N2CI*(*sub)->Ne);
      memset(filename, 0, sizeof(filename));
      strcpy(filename, filepath);
      strcat(filename, "Vxm.dat");
      ReadFile((*sub)->allVxmAR, filename, setting->N2CI*(*sub)->Ne);
      memset(filename, 0, sizeof(filename));
      strcpy(filename, filepath);
      strcat(filename, "Vym.dat");
      ReadFile((*sub)->allVymAR, filename, setting->N2CI*(*sub)->Ne);
      memset(filename, 0, sizeof(filename));
      // surface areas
      // strcpy(filename, filepath);
      // strcat(filename, "Zxp.dat");
      // ReadFile((*sub)->allZxpAR, filename, setting->N2CI*(*sub)->Ne);
      // memset(filename, 0, sizeof(filename));
      // strcpy(filename, filepath);
      // strcat(filename, "Zyp.dat");
      // ReadFile((*sub)->allZypAR, filename, setting->N2CI*(*sub)->Ne);
      // memset(filename, 0, sizeof(filename));
      // strcpy(filename, filepath);
      // strcat(filename, "Zxm.dat");
      // ReadFile((*sub)->allZxmAR, filename, setting->N2CI*(*sub)->Ne);
      // memset(filename, 0, sizeof(filename));
      // strcpy(filename, filepath);
      // strcat(filename, "Zym.dat");
      // ReadFile((*sub)->allZymAR, filename, setting->N2CI*(*sub)->Ne);
      // memset(filename, 0, sizeof(filename));
      //face areas
      strcpy(filename, filepath);
      strcat(filename, "Np.dat");
      ReadFile((*sub)->allNpAR, filename, setting->N2CI*(*sub)->Ne);
      memset(filename, 0, sizeof(filename));
      strcpy(filename, filepath);
      strcat(filename, "Op.dat");
      ReadFile((*sub)->allOpAR, filename, setting->N2CI*(*sub)->Ne);
      memset(filename, 0, sizeof(filename));
      strcpy(filename, filepath);
      strcat(filename, "Nm.dat");
      ReadFile((*sub)->allNmAR, filename, setting->N2CI*(*sub)->Ne);
      memset(filename, 0, sizeof(filename));
      strcpy(filename, filepath);
      strcat(filename, "Om.dat");
      ReadFile((*sub)->allOmAR, filename, setting->N2CI*(*sub)->Ne);
      memset(filename, 0, sizeof(filename));
    }
}

// assign subgrid areas to each rank
void splitSubArea(Sub **sub, Bath *bath, Config *setting, int irank)
{
  int ii, jj, kk;
  if (setting->useSubgrid != 0)
  {
    // split variables with free surface dependency
    for (kk = 0; kk < (*sub)->Ne; kk++)
    {
      jj = kk*setting->N2CI;
      for (ii = 0; ii < setting->N2ci; ii++)
      {
        (*sub)->allV[setting->N2ci*kk+ii] = (*sub)->allVAR[jj+irank*setting->N2ci+ii];
        (*sub)->allZ[setting->N2ci*kk+ii] = (*sub)->allZAR[jj+irank*setting->N2ci+ii];
        //(*sub)->allN[setting->N2ci*kk+ii] = (*sub)->allNAR[jj+irank*setting->N2ci+ii];
        //(*sub)->allO[setting->N2ci*kk+ii] = (*sub)->allOAR[jj+irank*setting->N2ci+ii];
        (*sub)->allVxp[setting->N2ci*kk+ii] = (*sub)->allVxpAR[jj+irank*setting->N2ci+ii];
        (*sub)->allVxm[setting->N2ci*kk+ii] = (*sub)->allVxmAR[jj+irank*setting->N2ci+ii];
        (*sub)->allVyp[setting->N2ci*kk+ii] = (*sub)->allVypAR[jj+irank*setting->N2ci+ii];
        (*sub)->allVym[setting->N2ci*kk+ii] = (*sub)->allVymAR[jj+irank*setting->N2ci+ii];
        (*sub)->allNp[setting->N2ci*kk+ii] = (*sub)->allNpAR[jj+irank*setting->N2ci+ii];
        (*sub)->allNm[setting->N2ci*kk+ii] = (*sub)->allNmAR[jj+irank*setting->N2ci+ii];
        (*sub)->allOp[setting->N2ci*kk+ii] = (*sub)->allOpAR[jj+irank*setting->N2ci+ii];
        (*sub)->allOm[setting->N2ci*kk+ii] = (*sub)->allOmAR[jj+irank*setting->N2ci+ii];
        // (*sub)->allZxp[setting->N2ci*kk+ii] = (*sub)->allZxpAR[jj+irank*setting->N2ci+ii];
        // (*sub)->allZxm[setting->N2ci*kk+ii] = (*sub)->allZxmAR[jj+irank*setting->N2ci+ii];
        // (*sub)->allZyp[setting->N2ci*kk+ii] = (*sub)->allZypAR[jj+irank*setting->N2ci+ii];
        // (*sub)->allZym[setting->N2ci*kk+ii] = (*sub)->allZymAR[jj+irank*setting->N2ci+ii];
        (*sub)->allYh[setting->N2ci*kk+ii] = (*sub)->allYhAR[jj+irank*setting->N2ci+ii];
        if (setting->useSubDrag == 3)
        {
            (*sub)->allCvX[setting->N2ci*kk+ii] = (*sub)->allCvXAR[jj+irank*setting->N2ci+ii];
            (*sub)->allCvY[setting->N2ci*kk+ii] = (*sub)->allCvYAR[jj+irank*setting->N2ci+ii];
          (*sub)->allCdsX[setting->N2ci*kk+ii] = (*sub)->allCdsXAR[jj+irank*setting->N2ci+ii];
          (*sub)->allCdsY[setting->N2ci*kk+ii] = (*sub)->allCdsYAR[jj+irank*setting->N2ci+ii];
        }
        if (setting->useSubDrag == 4 | setting->useSubDrag == 3)
        {(*sub)->allCn[setting->N2ci*kk+ii] = (*sub)->allCnAR[jj+irank*setting->N2ci+ii];}
        // (*sub)->allredCD[setting->N2ci*kk+ii] = (*sub)->allredCDAR[jj+irank*setting->N2ci+ii];
      }
    }
    // split variables without free surface dependency
    // for (ii = 0; ii < setting->N2ci; ii++)
    // {
    //   /*(*sub)->bath[ii] = (*sub)->bathAR[ii+irank*setting->N2ci];
    //   (*sub)->bathXP[ii] = (*sub)->bathXPAR[ii+irank*setting->N2ci];
    //   (*sub)->bathYP[ii] = (*sub)->bathYPAR[ii+irank*setting->N2ci];*/
    //   //(*sub)->CvX[ii] = (*sub)->CvXAR[ii+irank*setting->N2ci];
    //   //(*sub)->CvY[ii] = (*sub)->CvYAR[ii+irank*setting->N2ci];
    //   (*sub)->wdNp[ii] = (*sub)->wdNpAR[ii+irank*setting->N2ci] ;
    //   (*sub)->wdOp[ii] = (*sub)->wdOpAR[ii+irank*setting->N2ci] ;
    //   (*sub)->wdNm[ii] = (*sub)->wdNmAR[ii+irank*setting->N2ci] ;
    //   (*sub)->wdOm[ii] = (*sub)->wdOmAR[ii+irank*setting->N2ci] ;
    // }
  }
}

// read one subgrid variable, assign it to all ranks and clear it from the memory
void initOneSubVar(Sub **sub, char vname[], Config *setting, int irank)
{
    (*sub)->allvarAR = malloc(setting->N2CI*(*sub)->Ne*sizeof(double));
    (*sub)->allvar = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    int ii,jj,kk;
    char filepath[100];
    strcpy(filepath, setting->inputFolder);
    strcat(filepath, setting->subgridFolder);
    char filename[100];
    strcpy(filename, filepath);
    strcat(filename, vname);
    ReadFile((*sub)->allvarAR, filename, setting->N2CI*(*sub)->Ne);
    memset(filename, 0, sizeof(filename));
    for (kk = 0; kk < (*sub)->Ne; kk++)
    {
        jj = kk*setting->N2CI;
        for (ii = 0; ii < setting->N2ci; ii++)
        {(*sub)->allvar[setting->N2ci*kk+ii] = (*sub)->allvarAR[jj+irank*setting->N2ci+ii];}
    }
    free((*sub)->allvarAR);
}

// initialize all the subgrid variables
void initAllSubVar(Sub **sub, Bath *bath, Config *setting, int irank)
{
    int ii, N;
    *sub = malloc(sizeof(Sub));
    (*sub)->Ne = round((setting->surfmax - setting->surfmin)/setting->dsurf);
    (*sub)->allSurf = malloc((*sub)->Ne*sizeof(double));
    for (ii = 0; ii < (*sub)->Ne; ii++)
    {(*sub)->allSurf[ii] = setting->surfmin + ii*setting->dsurf + bath->offset[0];}
    N = (*sub)->Ne * setting->N2ci;
    // subgrid cell volume
    initOneSubVar(sub, "V.dat", setting, irank);
    (*sub)->allV = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    for (ii = 0; ii < N; ii++)
    {(*sub)->allV[ii] = (*sub)->allvar[ii];}
    // subgrid cell area Z
    initOneSubVar(sub, "Z.dat", setting, irank);
    (*sub)->allZ = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    for (ii = 0; ii < N; ii++)
    {(*sub)->allZ[ii] = (*sub)->allvar[ii];}
    // half cell volumes
    initOneSubVar(sub, "Vxp.dat", setting, irank);
    (*sub)->allVxp = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    for (ii = 0; ii < N; ii++)
    {(*sub)->allVxp[ii] = (*sub)->allvar[ii];}
    initOneSubVar(sub, "Vxm.dat", setting, irank);
    (*sub)->allVxm = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    for (ii = 0; ii < N; ii++)
    {(*sub)->allVxm[ii] = (*sub)->allvar[ii];}
    initOneSubVar(sub, "Vyp.dat", setting, irank);
    (*sub)->allVyp = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    for (ii = 0; ii < N; ii++)
    {(*sub)->allVyp[ii] = (*sub)->allvar[ii];}
    initOneSubVar(sub, "Vym.dat", setting, irank);
    (*sub)->allVym = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    for (ii = 0; ii < N; ii++)
    {(*sub)->allVym[ii] = (*sub)->allvar[ii];}
    // cell face areas
    initOneSubVar(sub, "Np.dat", setting, irank);
    (*sub)->allNp = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    for (ii = 0; ii < N; ii++)
    {(*sub)->allNp[ii] = (*sub)->allvar[ii];}
    initOneSubVar(sub, "Nm.dat", setting, irank);
    (*sub)->allNm = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    for (ii = 0; ii < N; ii++)
    {(*sub)->allNm[ii] = (*sub)->allvar[ii];}
    initOneSubVar(sub, "Op.dat", setting, irank);
    (*sub)->allOp = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    for (ii = 0; ii < N; ii++)
    {(*sub)->allOp[ii] = (*sub)->allvar[ii];}
    initOneSubVar(sub, "Om.dat", setting, irank);
    (*sub)->allOm = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    for (ii = 0; ii < N; ii++)
    {(*sub)->allOm[ii] = (*sub)->allvar[ii];}
    // min face area, added by ZhiLi20181004
    initOneSubVar(sub, "Nmin.dat", setting, irank);
    (*sub)->allNmin = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    for (ii = 0; ii < N; ii++)
    {(*sub)->allNmin[ii] = (*sub)->allvar[ii];}
    initOneSubVar(sub, "Omin.dat", setting, irank);
    (*sub)->allOmin = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
    for (ii = 0; ii < N; ii++)
    {(*sub)->allOmin[ii] = (*sub)->allvar[ii];}
    // transverse momentum model
    if (setting->phiSurface == 1 | setting->phiNonlinear == 1)
    {
        initOneSubVar(sub, "phiX.dat", setting, irank);
        (*sub)->allphiX = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
        for (ii = 0; ii < N; ii++)
        {(*sub)->allphiX[ii] = (*sub)->allvar[ii];}
        initOneSubVar(sub, "phiY.dat", setting, irank);
        (*sub)->allphiY = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
        for (ii = 0; ii < N; ii++)
        {(*sub)->allphiY[ii] = (*sub)->allvar[ii];}
    }
    // drag model by Volp
    if (setting->useSubDrag == 2)
    {
        initOneSubVar(sub, "Yh.dat", setting, irank);
        (*sub)->allYh = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
        for (ii = 0; ii < N; ii++)
        {(*sub)->allYh[ii] = (*sub)->allvar[ii];}
    }
    // drag model by Casas
    if (setting->useSubDrag == 4 | setting->useSubDrag == 3)
    {
        initOneSubVar(sub, "Cn.dat", setting, irank);
        (*sub)->allCn = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
        for (ii = 0; ii < N; ii++)
        {(*sub)->allCn[ii] = (*sub)->allvar[ii];}
    }
    // drag model by Li
    if (setting->useSubDrag == 3)
    {
        initOneSubVar(sub, "CvX.dat", setting, irank);
        (*sub)->allCvX = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
        for (ii = 0; ii < N; ii++)
        {(*sub)->allCvX[ii] = (*sub)->allvar[ii];}
        initOneSubVar(sub, "CvY.dat", setting, irank);
        (*sub)->allCvY = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
        for (ii = 0; ii < N; ii++)
        {(*sub)->allCvY[ii] = (*sub)->allvar[ii];}
        initOneSubVar(sub, "effCdX.dat", setting, irank);
        (*sub)->allCdsX = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
        for (ii = 0; ii < N; ii++)
        {(*sub)->allCdsX[ii] = (*sub)->allvar[ii];}
        initOneSubVar(sub, "effCdY.dat", setting, irank);
        (*sub)->allCdsY = malloc(setting->N2ci*(*sub)->Ne*sizeof(double));
        for (ii = 0; ii < N; ii++)
        {(*sub)->allCdsY[ii] = (*sub)->allvar[ii];}
    }
    // drag model based on roughness height
//    if (setting->useSubDrag == 5)
//    {
//      (*sub)->allZ0 = malloc(setting->N2CI*sizeof(double));
//      (*sub)->Z0 = malloc(setting->N2ci*sizeof(double));
//      ReadFile((*bath)->allZ0, "Z0.dat", setting->N2CI);
//      for (ii = 0; ii < setting->N2ci; ii++)
//      {(*sub)->Z0[ii] = (*sub)->allZ0[ii+irank*setting->N2ci];}
//    }
    // initialize the subgrid variables for one step one rank
    (*sub)->V = malloc(setting->N2ct*sizeof(double));
    (*sub)->Z = malloc(setting->N2ct*sizeof(double));
    (*sub)->Vxp = malloc(setting->N2ct*sizeof(double));
    (*sub)->Vyp = malloc(setting->N2ct*sizeof(double));
    (*sub)->Vxm = malloc(setting->N2ct*sizeof(double));
    (*sub)->Vym = malloc(setting->N2ct*sizeof(double));
    (*sub)->Np = malloc(setting->N2ct*sizeof(double));
    (*sub)->Op = malloc(setting->N2ct*sizeof(double));
    (*sub)->Nm = malloc(setting->N2ct*sizeof(double));
    (*sub)->Om = malloc(setting->N2ct*sizeof(double));
    (*sub)->Vx = malloc(setting->N2ct*sizeof(double));
    (*sub)->Vy = malloc(setting->N2ct*sizeof(double));
    (*sub)->Nx = malloc(setting->N2ct*sizeof(double));
    (*sub)->Oy = malloc(setting->N2ct*sizeof(double));
    (*sub)->Zx = malloc(setting->N2ct*sizeof(double));
    (*sub)->Zy = malloc(setting->N2ct*sizeof(double));
    // min area, ZhiLi20181004
    (*sub)->Nmin = malloc(setting->N2ct*sizeof(double));
    (*sub)->Omin = malloc(setting->N2ct*sizeof(double));
    if (setting->phiSurface == 1 | setting->phiNonlinear == 1)
    {
        (*sub)->phiX = malloc(setting->N2ct*sizeof(double));
        (*sub)->phiY = malloc(setting->N2ct*sizeof(double));
    }
    if (setting->useSubDrag == 2)
    {(*sub)->Yh = malloc(setting->N2ct*sizeof(double));}
    if (setting->useSubDrag == 4 | setting->useSubDrag == 3)
    {(*sub)->Cn = malloc(setting->N2ct*sizeof(double));}
    if (setting->useSubDrag == 3)
    {
        (*sub)->CvX = malloc(setting->N2ct*sizeof(double));
        (*sub)->CvY = malloc(setting->N2ct*sizeof(double));
        (*sub)->CdsX = malloc(setting->N2ct*sizeof(double));
        (*sub)->CdsY = malloc(setting->N2ct*sizeof(double));
    }
    (*sub)->allNx = malloc(setting->N2CI*sizeof(double));
    (*sub)->allOy = malloc(setting->N2CI*sizeof(double));
    (*sub)->allVx = malloc(setting->N2CI*sizeof(double));
    (*sub)->allVy = malloc(setting->N2CI*sizeof(double));
    (*sub)->allVo = malloc(setting->N2CI*sizeof(double));
    (*sub)->allZo = malloc(setting->N2CI*sizeof(double));
    // initialize the searching index
    (*sub)->ind = malloc(setting->N2ct*sizeof(int));
    for (ii = 0; ii < setting->N2ct; ii++)
    {(*sub)->ind[ii] = round((*sub)->Ne / 2);}
}
