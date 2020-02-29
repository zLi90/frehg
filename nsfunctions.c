#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>

// -----------------------------------------------------------------------------
// This file contains the functions used in solving the Navier Stokes equations
// other than the direct FD/FV solvers. These functions were mainly used for
// data processing purposes. However, the difference between this file and the
// utilities.c is that, functions defined in this file are all specifically
// designed for the FREHD data structures.
// - ZhiLi 2017-05-05 -
// -----------------------------------------------------------------------------

#include "configuration.h"
#include "bathymetry.h"
#include "initialize.h"
#include "map.h"
#include "fileio.h"
#include "utilities.h"

#include "laspack/errhandl.h"
#include "laspack/vector.h"
#include "laspack/matrix.h"
#include "laspack/qmatrix.h"
#include "laspack/operats.h"
#include "laspack/factor.h"
#include "laspack/precond.h"
#include "laspack/eigenval.h"
#include "laspack/rtc.h"
#include "laspack/itersolv.h"
#include "laspack/mlsolv.h"
#include "laspack/version.h"
#include "laspack/copyrght.h"

void updateAllFaceElev(Bath *bath, Maps *map, Config *setting);
void updateAllFaceDepth(Data **data, Bath *bath, Maps *map, Config *setting, int irank, int nrank);
void updateDepth(Data **data, Maps *map, Bath *bath, Config *setting);
void updateVelocity(Data **data, Maps *map, Config *setting);
void enforceVelocityBC(Data **data, Maps *map, Config *setting);
void enforceFreeSurfaceBC(Data **data, Maps *map, Config *setting);
void velocityInterp(Data **data, Config *setting);
void computeFaceFlowRate(Data **data, Maps *map, Config *setting);
void computeVolume(Data **data, Config *setting);
void volumeByFlux(Data **data, Maps *map, BC *bc, Config *setting, int tt, int irank);
void enforceTidalBC(Data **data, Bath *bath, BC *bc, Config *setting, int tt, int irank, int nrank);
void zeroNegSurf(Data **data, Bath *bath, Config *setting, int tt);
void limiterCFL(Data **data, Bath *bath, Maps *map, Config *setting);
void explicitVelocity(Data **data, Maps *map, Config *setting);
void advectionTerm(Data **data, Maps *map, Config *setting);
void diffusionTerm(Data **data, Maps *map, Config *setting);
void windTerm(Data **data, BC *bc, Config *setting, int tt);
void dragTerm(Data **data, Maps *map, Config *setting);
void thinLayerDrag(Data **data, Config *setting);
void dragInversion(Data **data, Maps *map, Config *setting);
void rainfall(Data **data, BC *bc, int tt, Config *setting);
void matrixSourceTerm(Data **data, Maps *map, BC *bc, Config *setting, int tt, int irank);
void findInflowLocation(Data **data, Maps *map, Config *setting);
void matrixCoeff(Data **data, Maps *map, Config *setting);
void setupMatrix(Data *data, Maps *map, Config *setting, QMatrix A);
void adjustMatrixBoundary(Data **data, BC *bc, Maps *map, Config *setting, int irank, int nrank);
void adjustTidalBoundary(Data **data, BC *bc, Bath *bath, Config *setting, int tt, int irank, int nrank);
void solveMatrix(Data *data, Config *setting, QMatrix A, Vector x, Vector z);
void getFreeSurface(Data **data, Maps *map, Config *setting, Vector x);
void adjustTidalVelocity(Data **data, Maps *map, Config *setting, int irank, int nrank);
void detectWaterfallLocation(Data **data, Bath *bath, Maps *map, Config *setting);
void waterfallCorrection(Data **data, Bath *bath, Maps *map, Config *setting);
void updateCD(Data **data, Config *setting);
void monitorCFL(Data **data, Bath *bath, int irank, Config *setting, int tt, int root);
void infiltration(Data **data, Bath *bath, Maps *map, Config *setting);


// =========== Compute cell edge elevation from center elevation ===========
void updateAllFaceElev(Bath *bath, Maps *map, Config *setting)
{
  int ii;
  for (ii = 0; ii < setting->N2ct; ii++)
  { bath->bottomXP[ii] = 0.0; bath->bottomYP[ii] = 0.0;}
  // update bottomXP
  // internal faces
  if (setting->useCellEdge == 0 & setting->useSubgrid == 0)
  {
    for (ii = 0; ii < setting->N2ci; ii++)
    {
      if (bath->bottomZ[ii] > bath->bottomZ[map->iPjc[ii]])
      {bath->bottomXP[ii] = bath->bottomZ[ii];}
      else
      {bath->bottomXP[ii] = bath->bottomZ[map->iPjc[ii]];}
    }
  }
  else
  {
    for (ii = 0; ii < setting->N2ci; ii++)
    {bath->bottomXP[ii] = bath->edgeX[ii];}
  }
  // boundary faces
  for (ii = 0; ii < setting->ny; ii++)
  {
    bath->bottomXP[map->iMgt[ii]] = bath->bottomZ[map->iMbd[ii]];
    bath->bottomXP[map->iPgt[ii]] = bath->bottomXP[map->iPbd[ii]];
  }
  for (ii = 0; ii < setting->nx; ii++)
  {
    bath->bottomXP[map->jMgt[ii]] = bath->bottomXP[map->jMbd[ii]];
    bath->bottomXP[map->jPgt[ii]] = bath->bottomXP[map->jPbd[ii]];
  }
  // update bottomYP
  // internal faces
  if (setting->useCellEdge == 0 & setting->useSubgrid == 0)
  {
    for (ii = 0; ii < setting->N2ci; ii++)
    {
      if (bath->bottomZ[ii] > bath->bottomZ[map->icjP[ii]])
      {bath->bottomYP[ii] = bath->bottomZ[ii];}
      else
      {bath->bottomYP[ii] = bath->bottomZ[map->icjP[ii]];}
    }
  }
  else
  {
    for (ii = 0; ii < setting->N2ci; ii++)
    {bath->bottomYP[ii] = bath->edgeY[ii];}
  }
  // boundary faces
  for (ii = 0; ii < setting->nx; ii++)
  {
    if (bath->bottomZ[map->jMgt[ii]] > bath->bottomZ[ii])
    {bath->bottomYP[map->jMgt[ii]] = bath->bottomZ[map->jMgt[ii]];}
    else
    {bath->bottomYP[map->jMgt[ii]] = bath->bottomZ[ii];}
    bath->bottomYP[map->jPgt[ii]] = bath->bottomYP[map->jPbd[ii]];
  }
  for (ii = 0; ii < setting->ny; ii++)
  {
    bath->bottomYP[map->iMgt[ii]] = bath->bottomYP[map->iMbd[ii]];
    bath->bottomYP[map->iPgt[ii]] = bath->bottomYP[map->iPbd[ii]];
  }
}

// ========== Compute cell edge depth from center depth ==========
void updateAllFaceDepth(Data **data, Bath *bath, Maps *map, Config *setting, int irank, int nrank)
{
  int ii;
  for (ii = 0; ii < setting->N2ct; ii++)
  { (*data)->depthXP[ii] = 0; (*data)->depthYP[ii] = 0;}
  // depthXP
  // internal faces
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    if ((*data)->surf[ii] > (*data)->surf[map->iPjc[ii]])
    {(*data)->depthXP[ii] = (*data)->surf[ii] - bath->bottomXP[ii];}
    else
    {(*data)->depthXP[ii] = (*data)->surf[map->iPjc[ii]] - bath->bottomXP[ii];}
  }
  // boundary faces
  for (ii = 0; ii < setting->ny; ii++)
  {
    (*data)->depthXP[map->iMgt[ii]] = (*data)->surf[map->iMbd[ii]] - \
      bath->bottomXP[map->iMgt[ii]];
    (*data)->depthXP[map->iPgt[ii]] = (*data)->surf[map->iPgt[ii]] - \
      bath->bottomXP[map->iPgt[ii]];
    // following FREHD?
    (*data)->depthXP[map->iPbd[ii]] = 0;
    (*data)->depthXP[map->iPgt[ii]] = 0;
    (*data)->depthXP[map->iMgt[ii]] = 0;
  }

  for (ii = 0; ii < setting->nx; ii++)
  {
    if ((*data)->surf[map->jMgt[ii]] > (*data)->surf[map->jMgt[ii]+1])
    {(*data)->depthXP[map->jMgt[ii]] = (*data)->surf[map->jMgt[ii]] - bath->bottomXP[map->jMgt[ii]];}
    else
    {(*data)->depthXP[map->jMgt[ii]] = (*data)->surf[map->jMgt[ii]+1] - bath->bottomXP[map->jMgt[ii]];}

    if ((*data)->surf[map->jPgt[ii]] > (*data)->surf[map->jPgt[ii]+1])
    {(*data)->depthXP[map->jPgt[ii]] = (*data)->surf[map->jPgt[ii]] - bath->bottomXP[map->jPgt[ii]];}
    else
    {(*data)->depthXP[map->jPgt[ii]] = (*data)->surf[map->jPgt[ii]+1] - bath->bottomXP[map->jPgt[ii]];}
  }
  // depthYP
  // internal faces
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    if ((*data)->surf[ii] > (*data)->surf[map->icjP[ii]])
    {(*data)->depthYP[ii] = (*data)->surf[ii] - bath->bottomYP[ii];}
    else
    {(*data)->depthYP[ii] = (*data)->surf[map->icjP[ii]] - bath->bottomYP[ii];}
  }
  // boundary faces
  for (ii = 0; ii < setting->nx; ii++)
  {
    if ((*data)->surf[map->jMgt[ii]] > (*data)->surf[map->jMbd[ii]])
    {
      (*data)->depthYP[map->jMgt[ii]] = (*data)->surf[map->jMgt[ii]] - \
        bath->bottomYP[map->jMgt[ii]];
    }
    else
    {
      (*data)->depthYP[map->jMgt[ii]] = (*data)->surf[map->jMbd[ii]] - \
        bath->bottomYP[map->jMgt[ii]];
    }
    (*data)->depthYP[map->jPgt[ii]] = (*data)->surf[map->jPgt[ii]] - \
      bath->bottomYP[map->jPgt[ii]];
    // following FREHD?
    //(*data)->depthYP[map->jPbd[ii]] = 0;
    //(*data)->depthYP[map->jPgt[ii]] = 0;
    (*data)->depthYP[map->jMgt[ii]] = 0;
  }
  for (ii = 0; ii < setting->ny; ii++)
  {
    (*data)->depthYP[map->iMgt[ii]] = (*data)->depthYP[map->iMbd[ii]];
    (*data)->depthYP[map->iPgt[ii]] = (*data)->depthYP[map->iPbd[ii]];
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

// ========== Compute depth from bathymetry and free surface ==========
void updateDepth(Data **data, Maps *map, Bath *bath, Config *setting)
{
  int ii;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    (*data)->depth[ii] = (*data)->surf[ii] - bath->bottomZ[ii];
    if ((*data)->depth[ii] < 0) {(*data)->depth[ii] = 0;}
  }
}

// ==================== Update velocities ====================
void updateVelocity(Data **data, Maps *map, Config *setting)
{
  int ii, count = 0;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    (*data)->uuXP[ii] = (*data)->EnXP[ii] - (setting->g * setting->dt / setting->dx) * \
      ((*data)->surf[map->iPjc[ii]] - (*data)->surf[ii]);// / (1 + (*data)->DragXP[ii]);
    (*data)->vvYP[ii] = (*data)->EnYP[ii] - (setting->g * setting->dt / setting->dy) * \
      ((*data)->surf[map->icjP[ii]] - (*data)->surf[ii]);// / (1 + (*data)->DragYP[ii]);
      // if (ii == 189*88+40-1)
      // {printf("En, Baro = %f, %f\n",(*data)->EnYP[ii], - (setting->g * setting->dt / setting->dy) * \
      //   ((*data)->surf[map->icjP[ii]] - (*data)->surf[ii]));}
  }
  // remove velocity on dry face
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    if ((*data)->depthXP[ii] == 0)
    {(*data)->uuXP[ii] = 0;}
    if ((*data)->depthYP[ii] == 0)
    {(*data)->vvYP[ii] = 0;}
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
}

// =============== Enforce velocity BC at boundaries ===============
void enforceVelocityBC(Data **data, Maps *map, Config *setting)
{
  int ii;
  for (ii = 0; ii < setting->nx; ii++)
  {
    (*data)->uuXP[map->jPgt[ii]] = (*data)->uuXP[map->jPbd[ii]];
    (*data)->uuXP[map->jMgt[ii]] = (*data)->uuXP[map->jMbd[ii]];
    (*data)->vvYP[map->jPgt[ii]] = (*data)->vvYP[map->jPbd[ii]];
    (*data)->vvYP[map->jMgt[ii]] = (*data)->vvYP[map->jMbd[ii]];
  }
  for (ii = 0; ii < setting->ny; ii++)
  {
    (*data)->uuXP[map->iPgt[ii]] = (*data)->uuXP[map->iPbd[ii]];
    (*data)->uuXP[map->iMgt[ii]] = (*data)->uuXP[map->iMbd[ii]];
    (*data)->vvYP[map->iPgt[ii]] = (*data)->vvYP[map->iPbd[ii]];
    (*data)->vvYP[map->iMgt[ii]] = (*data)->vvYP[map->iMbd[ii]];
  }
  (*data)->uuXP[map->iPjP[0]] = (*data)->uuXP[map->jPgt[setting->nx-1]];
  (*data)->uuXP[map->iPjM[0]] = (*data)->uuXP[map->jMgt[setting->nx-1]];
  (*data)->uuXP[map->iMjP[0]] = (*data)->uuXP[map->jPgt[0]];
  (*data)->uuXP[map->iMjM[0]] = (*data)->uuXP[map->jMgt[0]];
  (*data)->vvYP[map->iPjP[0]] = (*data)->vvYP[map->iPgt[setting->ny-1]];
  (*data)->vvYP[map->iPjM[0]] = (*data)->vvYP[map->iMgt[setting->ny-1]];
  (*data)->vvYP[map->iMjP[0]] = (*data)->vvYP[map->iPgt[0]];
  (*data)->vvYP[map->iMjM[0]] = (*data)->vvYP[map->iMgt[0]];
}

// =============== Enforce free surface BC at boundaries ===============
void enforceFreeSurfaceBC(Data **data, Maps *map, Config *setting)
{
  int ii;
  for (ii = 0; ii < setting->nx; ii++)
  {
    (*data)->surf[map->jPgt[ii]] = (*data)->surf[map->jPbd[ii]];
    (*data)->surf[map->jMgt[ii]] = (*data)->surf[map->jMbd[ii]];
  }
  for (ii = 0; ii < setting->ny; ii++)
  {
    (*data)->surf[map->iPgt[ii]] = (*data)->surf[map->iPbd[ii]];
    (*data)->surf[map->iMgt[ii]] = (*data)->surf[map->iMbd[ii]];
  }
  (*data)->surf[map->iPjP[0]] = (*data)->surf[map->iPgt[setting->ny-1]];
  (*data)->surf[map->iPjM[0]] = (*data)->surf[map->iMgt[setting->ny-1]];
  (*data)->surf[map->iMjP[0]] = (*data)->surf[map->iPgt[0]];
  (*data)->surf[map->iMjM[0]] = (*data)->surf[map->iMgt[0]];
}

// ========== Interpolate velocity from i face to j face ==========
void velocityInterp(Data **data, Config *setting)
{
  int ii;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    (*data)->uuYP[ii] = 0.25 * ((*data)->uuXP[ii] + (*data)->uuXP[ii-1] + \
      (*data)->uuXP[ii+setting->nx] + (*data)->uuXP[ii-1+setting->nx]);
    (*data)->vvXP[ii] = 0.25 * ((*data)->vvYP[ii] + (*data)->vvYP[ii+1] + \
      (*data)->vvYP[ii-setting->nx] + (*data)->vvYP[ii+1-setting->nx]);
  }
}

// =============== Compute face flow rates ===============
void computeFaceFlowRate(Data **data, Maps *map, Config *setting)
{
  int ii;
  for (ii = 0; ii < setting->N2ct; ii++)
  {
    (*data)->Fuu[ii] = (*data)->uuXP[ii] * setting->dy * (*data)->depthXP[ii];
    (*data)->Fvv[ii] = (*data)->vvYP[ii] * setting->dx * (*data)->depthYP[ii];
  }
}

// =============== compute cell volumes ===============
void computeVolume(Data **data, Config *setting)
{
  int ii;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    (*data)->cellV[ii] = setting->dx * setting->dy * (*data)->depth[ii];
    if (setting->useScalar == 1)
    {(*data)->Sm[ii] = (*data)->S[ii] * (*data)->cellV[ii];}
  }
}

// =============== compute cell volume from flow rates ===============
void volumeByFlux(Data **data, Maps *map, BC *bc, Config *setting, int tt, int irank)
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
    if (setting->useEvap > 0)
    {
      (*data)->cellV[ii] = (*data)->cellV[ii] - bc->evap[tt] * \
        setting->dt * setting->dx * setting->dy;
    }
    if (setting->useRain > 0)
    {
      if ((*data)->rainC[ii] > setting->minDepth)
  		{
  			(*data)->cellV[ii] = (*data)->cellV[ii] + (*data)->rainC[ii] * \
          		setting->dx * setting->dy;
  		}

    }
    // added by ZhiLi20180416
    if (setting->checkConservation == 1)
    {(*data)->Vloss[ii] = (*data)->depth[ii]*setting->dx*setting->dy - (*data)->cellV[ii];}
    if ((*data)->cellV[ii] < 0)
    {
      ind2 = floor(ii / setting->nx);
      ind1 = ii - ind2 * setting->nx;
      //printf("WARNING: Net flux exceeds cell volume at (%d,%d)...\n",ind1,ind2);
      (*data)->cellV[ii] = 0;
    }
  }
}

// =============== enforce boundary conditions at boundaries ===============
void enforceTidalBC(Data **data, Bath *bath, BC *bc, Config *setting, int tt, int irank, int nrank)
{
  int ii;
  // enforce tidal bc on YP boundary
  if (irank == nrank - 1 & setting->bcType != 2)
  {
    for (ii = 0; ii < setting->tideLocLengthP; ii++)
    {
      if (bc->tideP[tt] + bath->offset[0] > bath->bottomZ[setting->tideLocP[ii]])
      {(*data)->surf[setting->tideLocP[ii]] = bc->tideP[tt] + bath->offset[0];}
    }
  }
  // enforce tidal bc on YM boundary
  if (irank == 0 & setting->bcType != 1)
  {
    for (ii = 0; ii < setting->tideLocLengthM; ii++)
    {
      if (bc->tideM[tt] + bath->offset[0] > bath->bottomZ[setting->tideLocM[ii]])
      {(*data)->surf[setting->tideLocM[ii]] = bc->tideM[tt] + bath->offset[0];}
    }
  }
}

// =============== Remove surfaces that are lower than bottom ===============
void zeroNegSurf(Data **data, Bath *bath, Config *setting, int tt)
{
  int ii;
  double diff;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    if ((*data)->surf[ii] < bath->bottomZ[ii])
    {
      diff = bath->bottomZ[ii] - (*data)->surf[ii];
      (*data)->surf[ii] = bath->bottomZ[ii];
    }
  }
}

// =============== Restrict wetting within 1 cell ===============
void limiterCFL(Data **data, Bath *bath, Maps *map, Config *setting)
{
  int ii, count = 0;
  double diff, surf, bott;
  bool iP = 0, iM = 0, jP = 0, jM = 0;
  // restrict wetting within 1 cell
  if (setting->useRain == 0)
  {
      for (ii = 0; ii < setting->N2ci; ii++)
      {
        diff = (*data)->surf[ii] - bath->bottomZ[ii];
        if ((*data)->depth[ii] <= 0 & diff > 0)
        {
          surf = (*data)->surfOld[ii];
          bott = bath->bottomZ[ii];
          iP = (*data)->surfOld[map->iPjc[ii]] <= surf+setting->minDepth | (*data)->depth[map->iPjc[ii]] <= 0;
          iM = (*data)->surfOld[map->iMjc[ii]] <= surf+setting->minDepth | (*data)->depth[map->iMjc[ii]] <= 0;
          jP = (*data)->surfOld[map->icjP[ii]] <= surf+setting->minDepth | (*data)->depth[map->icjP[ii]] <= 0;
          jM = (*data)->surfOld[map->icjM[ii]] <= surf+setting->minDepth | (*data)->depth[map->icjM[ii]] <= 0;
          if (iP & iM & jP & jM)
          {
            // does not account for tidal boundary volume loss
            if (ii < setting->N2ci-setting->dx)
            {(*data)->Vloss1[0] += diff * setting->dx * setting->dy;}
            (*data)->surf[ii] = bott;
            (*data)->rmvCFL[ii] = 1;
          }
        }
      }
  }
  // remove shallow depth
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    diff = (*data)->surf[ii] - bath->bottomZ[ii];
    if (diff <= setting->minDepth & diff > 0)
    {
      (*data)->rmvCFL[ii] = 2;
      (*data)->Vloss2[0] += diff * setting->dx * setting->dy;
      (*data)->surf[ii] = bath->bottomZ[ii];
      (*data)->cellV[ii] = 0;
      (*data)->rmvCFL[ii] = 0;
      count += 1;
    }
  }
  // if (count > 0)
  // {printf("Number of shallow cells removed = %d\n",count);  count = 0;}
}

// ========== Initialize the source term by the explicit velocity ==========
void explicitVelocity(Data **data, Maps *map, Config *setting)
{
  int ii;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    (*data)->EnXP[ii] = (*data)->uuXP[ii];
    (*data)->EnYP[ii] = (*data)->vvYP[ii];
  }
}

// ==================== Compute the advection term ====================
// only the 1st order upwind scheme is implemented so far
void advectionTerm(Data **data, Maps *map, Config *setting)
{
  int ii;
  double advX, advY;
  double CFLh = 0.7, CFLl = 0.5, CFLx, CFLy;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    (*data)->advX[ii] = (0.5/setting->dx)*((*data)->uuXP[ii] + fabs((*data)->uuXP[ii])) * \
      ((*data)->uuXP[ii] - (*data)->uuXP[map->iMjc[ii]]) + \
      (0.5/setting->dx)*((*data)->uuXP[ii] - fabs((*data)->uuXP[ii])) * \
      ((*data)->uuXP[map->iPjc[ii]] - (*data)->uuXP[ii]) + \
      (0.5/setting->dy)*((*data)->vvXP[ii] + fabs((*data)->vvXP[ii])) * \
      ((*data)->uuXP[ii] - (*data)->uuXP[map->icjM[ii]]) + \
      (0.5/setting->dy)*((*data)->vvXP[ii] - fabs((*data)->vvXP[ii])) * \
      ((*data)->uuXP[map->icjP[ii]] - (*data)->uuXP[ii]);
    (*data)->advY[ii] = (0.5/setting->dx)*((*data)->uuYP[ii] + fabs((*data)->uuYP[ii])) * \
      ((*data)->vvYP[ii] - (*data)->vvYP[map->iMjc[ii]]) + \
      (0.5/setting->dx)*((*data)->uuYP[ii] - fabs((*data)->uuYP[ii])) * \
      ((*data)->vvYP[map->iPjc[ii]] - (*data)->vvYP[ii]) + \
      (0.5/setting->dy)*((*data)->vvYP[ii] + fabs((*data)->vvYP[ii])) * \
      ((*data)->vvYP[ii] - (*data)->vvYP[map->icjM[ii]]) + \
      (0.5/setting->dy)*((*data)->vvYP[ii] - fabs((*data)->vvYP[ii])) * \
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
}

// ==================== Compute the diffusion term ====================
void diffusionTerm(Data **data, Maps *map, Config *setting)
{
  int ii;
  double diffX, diffY;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    diffX = (setting->NUx/(setting->dx*setting->dx)) * \
      ((*data)->uuXP[map->iPjc[ii]] - 2*(*data)->uuXP[ii] + (*data)->uuXP[map->iMjc[ii]]) \
      + (setting->NUy/(setting->dy*setting->dy)) * \
      ((*data)->uuXP[map->icjM[ii]] - 2*(*data)->uuXP[ii] + (*data)->uuXP[map->icjP[ii]]);
    diffY = (setting->NUx/(setting->dx*setting->dx)) * \
      ((*data)->vvYP[map->iPjc[ii]] - 2*(*data)->vvYP[ii] + (*data)->vvYP[map->iMjc[ii]]) \
      + (setting->NUy/(setting->dy*setting->dy)) * \
      ((*data)->vvYP[map->icjM[ii]] - 2*(*data)->vvYP[ii] + (*data)->vvYP[map->icjP[ii]]);
    (*data)->EnXP[ii] = (*data)->EnXP[ii] + setting->dt * diffX;
    (*data)->EnYP[ii] = (*data)->EnYP[ii] + setting->dt * diffY;
  }
}

// ==================== Compute the wind term ====================
void windTerm(Data **data, BC *bc, Config *setting, int tt)
{
  double phi, omega, tau, Cw, pi = 3.1415926, rho = 1000, tauXP, tauYP;
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
      if ((*data)->depthXP[ii] > 0)
      {
        (*data)->EnXP[ii] = (*data)->EnXP[ii] + setting->dt * tauXP * cos(omega) / \
          ((*data)->depthXP[ii] * rho);
      }
      if ((*data)->depthYP[ii] > 0)
      {
        (*data)->EnYP[ii] = (*data)->EnYP[ii] - setting->dt * tauYP * sin(omega) / \
          ((*data)->depthYP[ii] * rho);
      }
    }
  }
}

// ==================== Compute the drag term ====================
void dragTerm(Data **data, Maps *map, Config *setting)
{
  int ii;
  double velx, vely;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    velx = sqrt((*data)->uuXP[ii]*(*data)->uuXP[ii] + \
      (*data)->vvXP[ii]*(*data)->vvXP[ii]);
    vely = sqrt((*data)->uuYP[ii]*(*data)->uuYP[ii] + \
      (*data)->vvYP[ii]*(*data)->vvYP[ii]);
    if ((*data)->depthXP[ii] > 0)
    {
      (*data)->DragXP[ii] = 0.5 * setting->dt * (*data)->CDXP[ii] * \
      velx / (*data)->depthXP[ii];
    }
    else
    {
      (*data)->DragXP[ii] = 0;
    }
    if ((*data)->depthYP[ii] > 0)
    {
      (*data)->DragYP[ii] = 0.5 * setting->dt * (*data)->CDYP[ii] * \
        vely / (*data)->depthYP[ii];
    }
    else
    {
      (*data)->DragYP[ii] = 0;
    }
  }
  for (ii = 0; ii < setting->nx; ii++)
  {
    (*data)->DragXP[map->jMgt[ii]] = (*data)->DragXP[map->jMbd[ii]];
  }
  for (ii = 0; ii < setting->ny; ii++)
  {
    (*data)->DragYP[map->iMgt[ii]] = (*data)->DragYP[map->iMbd[ii]];
  }
}

// ==================== Thin layer drag model =====================
void thinLayerDrag(Data **data, Config *setting)
{
  int ii;
  double coeff = 0;
  if (setting->CDnotN == 1)
  {
    for (ii = 0; ii < setting->N2ct; ii++)
    {
      // x direction
      if ((*data)->depthXP[ii] <= setting->hD)
      {(*data)->CDXP[ii] = setting->CDx + (setting->CDmax - setting->CDx) * \
        (setting->hD - (*data)->depthXP[ii]) / setting->hD;}
      //else
      //{(*data)->CDXP[ii] = setting->CDx;}
      // y direction
      if ((*data)->depthYP[ii] <= setting->hD)
      {(*data)->CDYP[ii] = setting->CDy + (setting->CDmax - setting->CDy) * \
        (setting->hD - (*data)->depthYP[ii]) / setting->hD;}
      //else
      //{(*data)->CDYP[ii] = setting->CDy;}
    }
  }
  else
  {
    coeff = setting->g * setting->manningN * setting->manningN;
    for (ii = 0; ii < setting->N2ci; ii++)
    {
      if ((*data)->depth[ii] > 0)
      {
        if ((*data)->depthXP[ii] <= setting->hD)
        {(*data)->CDXP[ii] = coeff / pow((*data)->depth[ii], 2.0/3.0);}
        if ((*data)->depthYP[ii] <= setting->hD)
        {(*data)->CDYP[ii] = coeff / pow((*data)->depth[ii], 2.0/3.0);}
      }
    }
  }
}

// ==================== Inversion for the implicit drag ====================
void dragInversion(Data **data, Maps *map, Config *setting)
{
  int ii;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    (*data)->EnXP[ii] = (*data)->EnXP[ii] / (1 + (*data)->DragXP[ii]);
    (*data)->EnYP[ii] = (*data)->EnYP[ii] / (1 + (*data)->DragYP[ii]);
  }
  for (ii = 0; ii < setting->nx; ii++)
  {
    (*data)->EnXP[map->jMgt[ii]] = (*data)->EnXP[map->jMbd[ii]];
  }
  for (ii = 0; ii < setting->ny; ii++)
  {
    (*data)->EnYP[map->iMgt[ii]] = (*data)->EnYP[map->iMbd[ii]];
  }
}

// =============== Compute rainfall ===============
void rainfall(Data **data, BC *bc, int tt, Config *setting)
{
    int ii;
//    double h = bc->rain[tt];
//    (*data)->rainSum[0] += h*setting->dt;
//    if ((*data)->rainSum[0] > setting->minDepth)
//    {
//        for (ii = 0; ii < setting->N2ci; ii++)
//        {(*data)->surf[ii] += (*data)->rainSum;}
//        (*data)->rainSum[0] = 0.0;
//    }
}

// ============== Compute the RHS of the matrix equation ===============
void matrixSourceTerm(Data **data, Maps *map, BC *bc, Config *setting, int tt, int irank)
{
  int ii, jj, loc[setting->inflowLocLength];
  double p1, p2, p3, p4;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    (*data)->z[ii] = (*data)->surf[map->trps[ii]] - \
       (setting->dt/setting->dx) * \
      ((*data)->EnXP[map->trps[ii]] * (*data)->depthXP[map->trps[ii]] - \
       (*data)->EnXP[map->iMjc[map->trps[ii]]] * (*data)->depthXP[map->iMjc[map->trps[ii]]]) - \
       (setting->dt/setting->dy) * \
      ((*data)->EnYP[map->trps[ii]] * (*data)->depthYP[map->trps[ii]] - \
       (*data)->EnYP[map->icjM[map->trps[ii]]] * (*data)->depthYP[map->icjM[map->trps[ii]]]);
  }
  // note: Here inflow is only added to rank0. Inflow on other ranks requires
  // further debugging. ZhiLi 20170820
  if (irank == 0 & setting->bcType != 3)
  {
    for (ii = 0; ii < setting->inflowLocLength; ii++)
    {
      loc[ii] = (setting->inflowLoc[ii]%setting->nx)*setting->ny + \
        floor(setting->inflowLoc[ii]/setting->nx);
      (*data)->z[loc[ii]] = (*data)->z[loc[ii]] + \
        (setting->dt/(setting->dx*setting->dy)) * (bc->inflow[tt] / setting->inflowLocLength);
    }
  }
  // add rainfall and evaporation as source terms
  if (setting->useRain > 0)
  {
      for (ii = 0; ii < setting->N2ci; ii++)
      {
          // if (map->ii2d[map->trps[ii]] == 1)
          {(*data)->z[ii] += bc->rain[tt] * setting->dt;}
      }
      // printf("Rainfall rate = %f\n",1e5*bc->rain[tt]);
  }
  // evaporation on wet regions only
  if (setting->useEvap > 0)
  {
      for (ii = 0; ii < setting->N2ci; ii++)
      {
          if ((*data)->depth[map->trps[ii]] > 0.0)
          {(*data)->z[ii] -= bc->evap[tt] * setting->dt;}
      }
      // printf("Evap rate = %f\n",1e5*bc->evap[tt]);
  }
}

// ===== find the location to add inflow in the transposed array =====
void findInflowLocation(Data **data, Maps *map, Config *setting)
{
  int ii, jj;
  for (ii = 0; ii < setting->inflowLocLength; ii++)
  {
    for (jj = 0; jj < setting->N2ci; jj++)
    {
      if (map->trps[jj] == setting->inflowLoc[ii])
      {
        (*data)->inflowLoc[ii] = jj;
        break;
      }
    }
  }
}

// =============== Compute the coefficients of the matrix A ===============
void matrixCoeff(Data **data, Maps *map, Config *setting)
{
  int ii,jj;
  double coefx, coefy;
  coefx = setting->g * (setting->dt*setting->dt/(setting->dx*setting->dx));
  coefy = setting->g * (setting->dt*setting->dt/(setting->dy*setting->dy));
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    //note : Theoretically we should apply the drag inversion when computing
    // these matrix coefficients, but it caused inconsistent flux. ZhiLi20170521

    (*data)->GnXP[ii] = coefx * (*data)->depthXP[ii] / (1 + (*data)->DragXP[ii]);
    (*data)->GnXM[ii] = coefx * (*data)->depthXP[map->iMjc[ii]] / \
      (1 + (*data)->DragXP[map->iMjc[ii]]);
    (*data)->GnYP[ii] = coefy * (*data)->depthYP[ii] / (1 + (*data)->DragYP[ii]);
    (*data)->GnYM[ii] = coefy * (*data)->depthYP[map->icjM[ii]] / \
      (1 + (*data)->DragYP[map->icjM[ii]]);
    (*data)->GnCt[ii] = 1 + (*data)->GnXP[ii] + (*data)->GnYP[ii] + (*data)->GnXM[ii] + \
      (*data)->GnYM[ii];

    /*(*data)->GnXP[ii] = coefx * (*data)->depthXP[ii];
    (*data)->GnXM[ii] = coefx * (*data)->depthXP[map->iMjc[ii]];
    (*data)->GnYP[ii] = coefy * (*data)->depthYP[ii];
    (*data)->GnYM[ii] = coefy * (*data)->depthYP[map->icjM[ii]];
    (*data)->GnCt[ii] = 1 + (*data)->GnXP[ii] + (*data)->GnYP[ii] + (*data)->GnXM[ii] + \
      (*data)->GnYM[ii];*/

  }
}

// ================ Adjust for boundary conditions ===============
void adjustMatrixBoundary(Data **data, BC *bc, Maps *map, Config *setting, int irank, int nrank)
{
  int ii, jj, kk;
  // adjust for Dirchlet boundaries
  // adjust for Neumann boundaries
  for (ii = 0; ii < setting->ny; ii++)
  {
    jj = ii*setting->nx;
    kk = ii*setting->nx + setting->nx - 1;
    (*data)->GnCt[jj] = (*data)->GnCt[jj] - (*data)->GnXM[jj];
    (*data)->GnXM[jj] = 0;
    (*data)->GnCt[kk] = (*data)->GnCt[kk] - (*data)->GnXP[kk];
    (*data)->GnXP[kk] = 0;
  }
  for (ii = 0; ii < setting->nx; ii++)
  {
    jj = ii;
    kk = setting->N2ci - setting->nx + ii;

    if (irank == 0)
    {(*data)->GnCt[jj] = (*data)->GnCt[jj] - (*data)->GnYM[jj];}
    else
    {(*data)->z[map->sprt[jj]] = (*data)->z[map->sprt[jj]] + \
      (*data)->GnYM[jj] * (*data)->surf[map->icjM[jj]];}
    (*data)->GnYM[jj] = 0;

    if (irank == nrank - 1)
    {(*data)->GnCt[kk] = (*data)->GnCt[kk] - (*data)->GnYP[kk];}
    else
    {(*data)->z[map->sprt[kk]] = (*data)->z[map->sprt[kk]] + \
      (*data)->GnYP[kk] * (*data)->surf[map->icjP[kk]];}
    (*data)->GnYP[kk] = 0;

    /*(*data)->GnCt[jj] = (*data)->GnCt[jj] - (*data)->GnYM[jj];
    (*data)->GnYM[jj] = 0;
    (*data)->GnCt[kk] = (*data)->GnCt[kk] - (*data)->GnYP[kk];
    (*data)->GnYP[kk] = 0;*/
  }
}

// =============== adjust coefficient at tidal boundary ===============
void adjustTidalBoundary(Data **data, BC *bc, Bath *bath, Config *setting, \
  int tt, int irank, int nrank)
{
  int ii, jj, kk;
  // adjust tide on YP boundary
  if (irank == nrank - 1 & setting->bcType != 2)
  {
    for (ii = 0; ii < setting->tideLocLengthP; ii++)
    {
      jj = setting->tideLocP[ii];
      kk = (jj%setting->nx)*setting->ny + floor(jj/setting->nx);
      (*data)->GnCt[jj] = 1;
      (*data)->GnXM[jj] = 0;
      (*data)->GnXP[jj] = 0;
      (*data)->GnYM[jj] = 0;
      (*data)->GnYP[jj] = 0;
      (*data)->z[kk] = bc->tideP[tt] + bath->offset[0];
    }
  }
  // adjust tide on YM boundary

  if (irank == 0 & setting->bcType != 1)
  {

    for (ii = 0; ii < setting->tideLocLengthM; ii++)
    {

      jj = setting->tideLocM[ii];
      kk = (jj%setting->nx)*setting->ny + floor(jj/setting->nx);
      (*data)->GnCt[jj] = 1;
      (*data)->GnXM[jj] = 0;
      (*data)->GnXP[jj] = 0;
      (*data)->GnYM[jj] = 0;
      (*data)->GnYP[jj] = 0;
      (*data)->z[kk] = bc->tideM[tt] + bath->offset[0];
    }
  }
}

// =============== Generate Matrix A and z ===============
void setupMatrix(Data *data, Maps *map, Config *setting, QMatrix A)
{
  size_t ii;
  int jj;
   // for (jj = 0; jj < setting->N2ci; jj++)
   // {
   //     printf("%d,%f,%f,%f,%f,%f,%f\n",jj,data->GnCt[jj], \
   //          data->GnXP[jj],data->GnXM[jj],data->GnYP[jj],data->GnYM[jj], \
   //          data->z[map->sprt[jj]]);
   // }
  Q_SetLen(&A, 1, 3);
  Q_SetEntry(&A, 1, 0, 1, data->GnCt[0]);
  Q_SetEntry(&A, 1, 1, 2, -data->GnYP[0]);
  Q_SetEntry(&A, 1, 2, setting->ny+1, -data->GnXP[0]);
  for (ii = 2; ii < setting->N2ci; ii++)
  {
    jj = map->trps[ii - 1];

    if (ii < setting->ny)
    {
      Q_SetLen(&A, ii, 4);
      Q_SetEntry(&A, ii, 0, ii-1, -data->GnYM[jj]);
      Q_SetEntry(&A, ii, 1, ii, data->GnCt[jj]);
      Q_SetEntry(&A, ii, 2, ii+1, -data->GnYP[jj]);
      Q_SetEntry(&A, ii, 3, ii+setting->ny, -data->GnXP[jj]);
    }
    else if (ii == setting->ny)
    {
      Q_SetLen(&A, ii, 3);
      Q_SetEntry(&A, ii, 0, ii-1, -data->GnYM[jj]);
      Q_SetEntry(&A, ii, 1, ii, data->GnCt[jj]);
      Q_SetEntry(&A, ii, 2, ii+setting->ny, -data->GnXP[jj]);
    }
    else if (ii > setting->N2ci-setting->ny+1)
    {
      Q_SetLen(&A, ii, 4);
      Q_SetEntry(&A, ii, 0, ii-setting->ny, -data->GnXM[jj]);
      Q_SetEntry(&A, ii, 1, ii-1, -data->GnYM[jj]);
      Q_SetEntry(&A, ii, 2, ii, data->GnCt[jj]);
      Q_SetEntry(&A, ii, 3, ii+1, -data->GnYP[jj]);
    }
    else if (ii == setting->N2ci-setting->ny+1)
    {
      Q_SetLen(&A, ii, 3);
      Q_SetEntry(&A, ii, 0, ii-setting->ny, -data->GnXM[jj]);
      Q_SetEntry(&A, ii, 1, ii, data->GnCt[jj]);
      Q_SetEntry(&A, ii, 2, ii+1, -data->GnYP[jj]);
    }
    else
    {
      if (ii%setting->ny == 0)
      {
        Q_SetLen(&A, ii, 4);
        Q_SetEntry(&A, ii, 0, ii-setting->ny, -data->GnXM[jj]);
        Q_SetEntry(&A, ii, 1, ii-1, -data->GnYM[jj]);
        Q_SetEntry(&A, ii, 2, ii, data->GnCt[jj]);
        Q_SetEntry(&A, ii, 3, ii+setting->ny, -data->GnXP[jj]);
      }
      else if (ii%setting->ny == 1)
      {
        Q_SetLen(&A, ii, 4);
        Q_SetEntry(&A, ii, 0, ii-setting->ny, -data->GnXM[jj]);
        Q_SetEntry(&A, ii, 1, ii, data->GnCt[jj]);
        Q_SetEntry(&A, ii, 2, ii+1, -data->GnYP[jj]);
        Q_SetEntry(&A, ii, 3, ii+setting->ny, -data->GnXP[jj]);
      }
      else
      {
        Q_SetLen(&A, ii, 5);
        Q_SetEntry(&A, ii, 0, ii-setting->ny, -data->GnXM[jj]);
        Q_SetEntry(&A, ii, 1, ii-1, -data->GnYM[jj]);
        Q_SetEntry(&A, ii, 2, ii, data->GnCt[jj]);
        Q_SetEntry(&A, ii, 3, ii+1, -data->GnYP[jj]);
        Q_SetEntry(&A, ii, 4, ii+setting->ny, -data->GnXP[jj]);
      }
    }
  }
  Q_SetLen(&A, setting->N2ci, 3);
  Q_SetEntry(&A, setting->N2ci, 0, setting->N2ci-setting->ny, -data->GnXM[setting->N2ci-1]);
  Q_SetEntry(&A, setting->N2ci, 1, setting->N2ci-1, -data->GnYM[setting->N2ci-1]);
  Q_SetEntry(&A, setting->N2ci, 2, setting->N2ci, data->GnCt[setting->N2ci-1]);

}

// =============== Solve the linear system Ax = z ===============
void solveMatrix(Data *data, Config *setting, QMatrix A, Vector x, Vector z)
{
  size_t ii;
  // create the z Vector
  for (ii = 1; ii <= setting->N2ci; ii++)
  {V_SetCmp(&z, ii, data->z[ii-1]);}
  // initialize x
  V_SetAllCmp(&x, 0.0);
  // set stopping criteria
  SetRTCAccuracy(setting->eps);
  // solve the linear system
  // JacobiIter(&A, &x, &z, setting->maxIter, NULL, 1.1);
  CGIter(&A, &x, &z, setting->maxIter, SSORPrecond, 1);
}

// =============== Assign x back to surface ===============
void getFreeSurface(Data **data, Maps *map, Config *setting, Vector x)
{
  int ii;
  size_t kk;
  for (ii = 0; ii < setting->N2ct; ii++)
  {(*data)->surfOld[ii] = (*data)->surf[ii];}
  for (kk = 0; kk < setting->N2ci; kk++)
  {
    ii = map->trps[kk];
    (*data)->surf[ii] = V_GetCmp(&x, kk+1);
  }
  // printf("281,Xm,Ym,Ct,Yp,Xp,z,surf=%f,%f,%f,%f,%f,%f,%f\n",(*data)->GnXM[281],(*data)->GnYM[281],(*data)->GnCt[281],(*data)->GnYP[281],(*data)->GnXP[281],(*data)->z[map->sprt[281]],(*data)->surf[281]);
  // printf("284,Xm,Ym,Ct,Yp,Xp,z,surf=%f,%f,%f,%f,%f,%f,%f\n",(*data)->GnXM[284],(*data)->GnYM[284],(*data)->GnCt[284],(*data)->GnYP[284],(*data)->GnXP[284],(*data)->z[map->sprt[284]],(*data)->surf[284]);
}

// ===================== Adjust tidal boundary velocity ====================
void adjustTidalVelocity(Data **data, Maps *map, Config *setting, int irank, int nrank)
{
  int ii, jj;
  double allFlux = 0, Fww = 0;
  // tidal velocity on YP boundary
  if (irank == nrank - 1 & setting->bcType != 2)
  {
    for (ii = 0; ii < setting->tideLocLengthP; ii++)
    {
      jj = setting->tideLocP[ii];
      // vertical flux due to change in tidal bc
      Fww = ((*data)->surf[jj] - (*data)->surfOld[jj]) * setting->dx * \
        setting->dy / setting->dt;
      // flux on the yp face from continuity
      allFlux = (*data)->Fvv[map->icjM[jj]] + (*data)->Fuu[map->iMjc[jj]]  - \
        (*data)->Fuu[jj] - Fww;
      // velocity is flux / area
      if (allFlux != 0 & (*data)->depth[jj] > 0)
      {(*data)->vvYP[jj] = allFlux / ((*data)->depth[jj] * setting->dx);}
      else
      {(*data)->vvYP[jj] = 0;}
      Fww = 0;
      allFlux = 0;
    }
  }
  // tidal velocity on YM boundary
  if (irank == 0 & setting->bcType != 1)
  {
    for (ii = 0; ii < setting->tideLocLengthM; ii++)
    {
      jj = setting->tideLocM[ii];
      // vertical flux due to change in tidal bc
      Fww = ((*data)->surf[jj] - (*data)->surfOld[jj]) * setting->dx * \
        setting->dy / setting->dt;
      //printf("jj,surf,surfOld,Fww = %d,%lf,%lf,%lf\n",jj,(*data)->surf[jj],(*data)->surfOld[jj],Fww);
      // flux on the ym face from continuity
      allFlux = -(*data)->Fvv[jj] + (*data)->Fuu[map->iMjc[jj]]  - \
        (*data)->Fuu[jj] - Fww;
      // velocity is flux / area
      if (allFlux != 0 & (*data)->depth[jj] > 0)
      {(*data)->vvYP[map->icjM[jj]] = -allFlux / ((*data)->depth[jj] * setting->dx);}
      Fww = 0;
      allFlux = 0;
    }
  }
}
void detectWaterfallLocation(Data **data, Bath *bath, Maps *map, Config *setting)
{
  int ii, count = 0;
  double surf, bot, dep, dmin = 0.0;
  // find waterfall locations
  // waterfall on XP face
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    surf = (*data)->surf[ii];
    bot = bath->bottomXP[ii];
    dep = (*data)->depthXP[ii];
    if (surf < bot-dmin & dep >= setting->wtfh)// & (bot-surf) > Ck*dep)
    {(*data)->wtfXP[ii] = -1;   count += 1;}
    //else if (surf < bot-dmin & dep < setting->wtfh)// & dep > 0 & (bot-surf) > Ck*dep)
    //{(*data)->wtfXP[ii] = -2;}
  }
  // waterfall on XM face
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    surf = (*data)->surf[ii];
    bot = bath->bottomXP[map->iMjc[ii]];
    dep = (*data)->depthXP[map->iMjc[ii]];
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
    dep = (*data)->depthYP[ii];
    if (surf < bot-dmin & dep >= setting->wtfh)// & (bot-surf) > Ck*dep)
    {(*data)->wtfYP[ii] = -1;   count += 1;}
    //else if (surf < bot-dmin & dep < setting->wtfh)// & dep > 0 & (bot-surf) > Ck*dep)
    //{(*data)->wtfYP[ii] = -2;}
  }
  // waterfall on YM face
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    surf = (*data)->surf[ii];
    bot = bath->bottomYP[map->icjM[ii]];
    dep = (*data)->depthYP[map->icjM[ii]];
    if (surf < bot-dmin & dep >= setting->wtfh)// & (bot-surf) > Ck*dep)
    {(*data)->wtfYP[map->icjM[ii]] = 1;   count += 1;}
    //else if (surf < bot-dmin & dep < setting->wtfh)// & dep > 0 & (bot-surf) > Ck*dep)
    //{(*data)->wtfYP[map->icjM[ii]] = 2;}
  }
  // if (count > 0)
  // {printf("Number of cell faces with waterfall = %d\n",count);  count = 0;}
}
// ===================== The water fall model ====================
void waterfallCorrection(Data **data, Bath *bath, Maps *map, Config *setting)
{
  int ii;
  double Cw = 0.7, Ck = 0.1, Cv = 0, CFLmax = 0.24, surf, bot, dep, vel0;
  // correct waterfall velocities
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    if ((*data)->wtfXP[ii] != 0)
    {
      vel0 = (*data)->uuXP[ii];
      (*data)->uuXP[ii] = (*data)->wtfXP[ii] * \
      Cw * sqrt(2 * setting->g * (*data)->depthXP[ii]);
      if ((*data)->uuXP[ii] > CFLmax * setting->dx / setting->dt)
      {(*data)->uuXP[ii] = CFLmax * setting->dx / setting->dt;}
      (*data)->Vloss4[0] -= (vel0 - (*data)->uuXP[ii]) * \
          setting->dt * (*data)->depthXP[ii] * setting->dy;
    }
    if ((*data)->wtfYP[ii] != 0)
    {
      vel0 = (*data)->vvYP[ii];
      (*data)->vvYP[ii] = (*data)->wtfYP[ii] * \
      Cw * sqrt(2 * setting->g * (*data)->depthYP[ii]);
      if ((*data)->vvYP[ii] > CFLmax * setting->dy / setting->dt)
      {(*data)->vvYP[ii] = CFLmax * setting->dy / setting->dt;}
      (*data)->Vloss4[0] -= (vel0 - (*data)->vvYP[ii]) * \
          setting->dt * (*data)->depthYP[ii] * setting->dx;
    }

    // if (map->jj2d[ii] < 5 & map->ii2d[ii] == 1)
    // {printf("jj, wtf, v, vel0 depth = %d, %f, %f, %f, %f\n",map->jj2d[ii], (*data)->wtfYP[ii], (*data)->vvYP[ii], vel0, (*data)->depth[ii]);}
  }
  // reset waterfall locations
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    (*data)->wtfXP[ii] = 0;
    (*data)->wtfYP[ii] = 0;
  }
}

// =============== Update the bottom drag coefficient ===============
void updateCD(Data **data, Config *setting)
{
  // the von Karman constant
  double K = 0.41, effh;
  int ii;
  double coeff = 0;
  // update CD only when the CD is depth dependent
  if (setting->CDnotN == 0)
  {
    coeff = setting->g * setting->manningN * setting->manningN;
    for (ii = 0; ii < setting->N2ci; ii++)
    {
      if ((*data)->depth[ii] > 0)
      {
        (*data)->CDXP[ii] = coeff / pow((*data)->depth[ii], 1.0/3.0);
        (*data)->CDYP[ii] = coeff / pow((*data)->depth[ii], 1.0/3.0);
      }
    }
  }
}

// =============== Monitoring the CFL number ===============
void monitorCFL(Data **data, Bath *bath, int irank, Config *setting, int tt, int root)
{
  int ii, locxX, locxY, locyX, locyY, locC;
  double CFLx0 = 0, CFLy0 = 0;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    (*data)->CFLx[ii] = fabs((*data)->uuXP[ii]) * setting->dt / setting->dx;
    (*data)->CFLy[ii] = fabs((*data)->vvYP[ii]) * setting->dt / setting->dy;

    locxY = floor(ii/setting->nx);
    locxX = ii - setting->nx*floor(ii/setting->nx);
    locyY = floor(ii/setting->nx);
    locyX = ii - setting->nx*floor(ii/setting->nx);
    locC = ii;

    if ((*data)->CFLx[ii] > 2)
    {
        printf("WARNING: Max CFL in x direction is %lf at %d, (%d,%d), \
          step %d, rank %d! Velocity is %lf \n", \
          (*data)->CFLx[ii], locC, locxX, locxY, tt, irank, (*data)->uuXP[ii]);
    }
    if ((*data)->CFLy[ii] > 2)
    {
        printf("WARNING: Max CFL in y direction is %lf at %d, (%d,%d), \
          step %d, rank %d! Velocity is %lf \n", \
          (*data)->CFLy[ii], locC, locyX, locyY, tt, irank, (*data)->vvYP[ii]);
    }
    if ((*data)->CFLx[ii] > CFLx0)
    {CFLx0 = (*data)->CFLx[ii];}
    if ((*data)->CFLy[ii] > CFLy0)
    {CFLy0 = (*data)->CFLy[ii];}
  }

  // if CFL number is too big, force the program to quit
  if (fabs(CFLx0) > 20 | fabs(CFLy0) > 20)
  {
    //DataOutput(data, bath, setting, tt, root, irank);
    printf("WARNING: Max CFL exceeds 20, program forced to quit!\n\n");
    exit(EXIT_FAILURE);
  }
}

// ===================== Surface-subsurface exchange ====================
void infiltration(Data **data, Bath *bath, Maps *map, Config *setting)
{
    // surface-subsurface exchange, ZhiLi 20190109
    int ii;
    if (setting->useSubsurface == 1)
    {
      for (ii = 0; ii < setting->N2ci; ii++)
      {
          if ((*data)->depth[ii] > 0 & (*data)->Qseep[ii] < 0)
          {
              (*data)->surf[ii] -= -(*data)->Qseep[ii] * setting->dt * (setting->dt/setting->dtg) / (setting->dx * setting->dy);
              if ((*data)->surf[ii] < bath->bottomZ[ii] + setting->minDepth)
              {(*data)->surf[ii] = bath->bottomZ[ii];}
              (*data)->depth[ii] = (*data)->surf[ii] - bath->bottomZ[ii];
          }
          else if ((*data)->Qseep[ii] > 0)
          {
              (*data)->surf[ii] += (*data)->Qseep[ii] * setting->dt * (setting->dt/setting->dtg) / (setting->dx * setting->dy);
              (*data)->depth[ii] = (*data)->surf[ii] - bath->bottomZ[ii];
          }
      }
    }
}
