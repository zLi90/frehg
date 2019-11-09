#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

// -----------------------------------------------------------------------------
// This file contains functions used for scalar computations
// - ZhiLi 2017-05-13 -
// -----------------------------------------------------------------------------

#include "configuration.h"
#include "bathymetry.h"
#include "map.h"
#include "fileio.h"
#include "mpifunctions.h"
#include "scalar.h"
#include "subgrid.h"
#include "utilities.h"

void scalarTransport(Data **data, Maps *map, Bath *bath, BC *bc, Sub *sub, Config *setting, int tt, int irank, int nrank);
void scalarMass(Data **data, Config *setting);
void scalarAdvection(Data **data, Maps *map, Config *setting);
void scalarDiffusion(Data **data, Maps *map, Config *setting);
void updateScalar(Data **data, BC *bc, Bath *bath, Maps *map, Config *setting, int tt, int irank, int nrank);
void enforceScalarBC(Data **data, BC *bc, Maps *map, Config *setting, int tt, int irank, int nrank);
void scalarSubStepping(Data **data, Bath *bath, BC *bc, Maps *map, Config *setting, int tt, int irank, int nrank);

// ==================== Scalar top level script ====================
void scalarTransport(Data **data, Maps *map, Bath *bath, BC *bc, Sub *sub, Config *setting, int tt, int irank, int nrank)
{
  if (setting->useSubStep == 1)
  {scalarSubStepping(data, bath, bc, map, setting, tt, irank, nrank);}
  else
  {
      scalarAdvection(data, map, setting);
      if (setting->useSubgrid != 0)
      {
        scalarDiffusionSub(data, map, sub, setting);
        updateScalarSub(data, bc, bath, map, sub, setting, tt, irank, nrank);
      }
      else
      {
        scalarDiffusion(data, map, setting);
        updateScalar(data, bc, bath, map, setting, tt, irank, nrank);
      }
  }
  enforceScalarBC(data, bc, map, setting, tt, irank, nrank);
  //scalarMass(data, setting);
  if (setting->useMPI == 1)
  {
    mpiexchange((*data)->S, map, setting, irank, nrank);
    mpiexchange((*data)->Sm, map, setting, irank, nrank);
  }
}

// =============== compute mass of scalar in each cell ===============
void scalarMass(Data **data, Config *setting)
{
  int ii;
  for (ii = 0; ii < setting->N2ci; ii++)
  {(*data)->Sm[ii] = (*data)->S[ii] * (*data)->cellV[ii];}
}

// ==================== advective scalar flux ====================
void scalarAdvection(Data **data, Maps *map, Config *setting)
{
  int ii;
  // ---------- 1st order upwind ----------
  if (setting->scalarAdv == 1)
  {
    for (ii = 0; ii < setting->N2ci; ii++)
    {
        double Su, Sd, Sr, Sl;
        if ((*data)->Fuu[ii] > 0) {Sd = (*data)->S[ii];} else {Sd = (*data)->S[map->iPjc[ii]];}
        if ((*data)->Fuu[map->iMjc[ii]] > 0) {Su = (*data)->S[map->iMjc[ii]];} else {Su = (*data)->S[ii];}
        if ((*data)->Fvv[ii] > 0) {Sr = (*data)->S[ii];} else {Sr = (*data)->S[map->icjP[ii]];}
        if ((*data)->Fvv[map->icjM[ii]] > 0) {Sl = (*data)->S[map->icjM[ii]];} else {Sl = (*data)->S[ii];}
      (*data)->Sm[ii] = (*data)->Sm[ii] + setting->dt * \
        (-(*data)->Fuu[ii] * Sd + (*data)->Fuu[map->iMjc[ii]] * Su - \
         (*data)->Fvv[ii] * Sr + (*data)->Fvv[map->icjM[ii]] * Sl);
		// subsurface scalar, ZhiLi20180112
		if (setting->useSubsurface == 1 & setting->useSubscalar == 1)
		{(*data)->Sm[ii] = (*data)->Sm[ii] + setting->dt * (*data)->Qseep[ii] * (*data)->Ssub[ii];}
    }
  }
  // ---------- TVD with superbee limiter, ZhiLi20180405 ----------
  else if (setting->scalarAdv == 2)
  {
    double r, phi, p1, p2, Sr, Sl, Su, Sd, Cu, dSmin = 0.00001;
    for (ii = 0; ii < setting->N2ci; ii++)
    {
      // advective flux at XP face
      Cu = (*data)->uuXP[ii] * setting->dt / setting->dx;
      if ((*data)->Fuu[ii] > 0)
      {
        // if the stencil extends to dry cells, switch back to 1st upwind
        if ((*data)->Fuu[map->iMjc[ii]] <= 0)
        {Sd = (*data)->S[ii];}
        else
        {
          if (fabs((*data)->S[map->iPjc[ii]] - (*data)->S[ii]) < dSmin)
          {r = 0;}
          else
          {r = ((*data)->S[ii] - (*data)->S[map->iMjc[ii]]) / ((*data)->S[map->iPjc[ii]] - (*data)->S[ii]);}
          if (2*r > 1) {p1 = 1;} else {p1 = 2*r;}
          if (r > 2) {p2 = 2;} else {p2 = r;}
          if (p1 > 0) {phi = p1;} else {phi = 0;}
          if (p2 > phi) {phi = p2;}
          Sd = (*data)->S[ii] + 0.5*phi*((*data)->S[map->iPjc[ii]] - (*data)->S[ii])*Cu*(1-Cu);
        }
      }
      else
      {
        // if the stencil extends to dry cells, switch back to 1st upwind
        if ((*data)->Fuu[map->iPjc[ii]] >= 0)
        {Sd = (*data)->S[map->iPjc[ii]];}
        else
        {
          if (fabs((*data)->S[map->iPjc[ii]] - (*data)->S[ii]) < dSmin)
          {r = 0;}
          else
          {r = ((*data)->S[map->iPjc[ii]] - (*data)->S[map->iPPjc[ii]])/((*data)->S[ii] - (*data)->S[map->iPjc[ii]]);}
          if (2*r > 1) {p1 = 1;} else {p1 = 2*r;}
          if (r > 2) {p2 = 2;} else {p2 = r;}
          if (p1 > 0) {phi = p1;} else {phi = 0;}
          if (p2 > phi) {phi = p2;}
          Sd = (*data)->S[map->iPjc[ii]] + 0.5*phi*((*data)->S[ii] - (*data)->S[map->iPjc[ii]])*Cu*(1-Cu);
        }
      }
      // advective flux at XM face
      Cu = (*data)->uuXP[map->iMjc[ii]] * setting->dt / setting->dx;
      if ((*data)->Fuu[map->iMjc[ii]] > 0)
      {
        if ((*data)->Fuu[map->iMMjc[ii]] <= 0)
        {Su = (*data)->S[map->iMjc[ii]];}
        else
        {
          if (fabs((*data)->S[ii] - (*data)->S[map->iMjc[ii]]) < dSmin)
          {r = 0;}
          else
          {r = ((*data)->S[map->iMjc[ii]] - (*data)->S[map->iMMjc[ii]])/((*data)->S[ii] - (*data)->S[map->iMjc[ii]]);}
          if (2*r > 1) {p1 = 1;} else {p1 = 2*r;}
          if (r > 2) {p2 = 2;} else {p2 = r;}
          if (p1 > 0) {phi = p1;} else {phi = 0;}
          if (p2 > phi) {phi = p2;}
          Su = (*data)->S[map->iMjc[ii]] + 0.5*phi*((*data)->S[ii] - (*data)->S[map->iMjc[ii]])*Cu*(1-Cu);
        }
      }
      else
      {
        if ((*data)->Fuu[ii] >= 0)
        {Su = (*data)->S[ii];}
        else
        {
          if (fabs((*data)->S[ii] - (*data)->S[map->iMjc[ii]]) < dSmin)
          {r = 0;}
          else
          {r = ((*data)->S[ii] - (*data)->S[map->iPjc[ii]])/((*data)->S[map->iMjc[ii]] - (*data)->S[ii]);}
          if (2*r > 1) {p1 = 1;} else {p1 = 2*r;}
          if (r > 2) {p2 = 2;} else {p2 = r;}
          if (p1 > 0) {phi = p1;} else {phi = 0;}
          if (p2 > phi) {phi = p2;}
          Su = (*data)->S[ii] + 0.5*phi*((*data)->S[map->iMjc[ii]] - (*data)->S[ii])*Cu*(1-Cu);
        }
      }
      // advective flux at YP face
      Cu = (*data)->vvYP[ii] * setting->dt / setting->dy;
      if ((*data)->Fvv[ii] > 0)
      {
        if ((*data)->Fvv[map->icjM[ii]] <= 0)
        {Sr = (*data)->S[ii];}
        else
        {
          if (fabs((*data)->S[map->icjP[ii]] - (*data)->S[ii]) < dSmin)
          {r = 0;}
          else
          {r = ((*data)->S[ii] - (*data)->S[map->icjM[ii]]) / ((*data)->S[map->icjP[ii]] - (*data)->S[ii]);}
          if (2*r > 1) {p1 = 1;} else {p1 = 2*r;}
          if (r > 2) {p2 = 2;} else {p2 = r;}
          if (p1 > 0) {phi = p1;} else {phi = 0;}
          if (p2 > phi) {phi = p2;}
          Sr = (*data)->S[ii] + 0.5*phi*((*data)->S[map->icjP[ii]] - (*data)->S[ii])*Cu*(1-Cu);
        }
      }
      else
      {
        if ((*data)->Fvv[map->icjP[ii]] >= 0)
        {Sr = (*data)->S[map->icjP[ii]];}
        else
        {
          if (fabs((*data)->S[map->icjP[ii]] - (*data)->S[ii]) < dSmin)
          {r = 0;}
          else
          {r = ((*data)->S[map->icjP[ii]] - (*data)->S[map->icjPP[ii]])/((*data)->S[ii] - (*data)->S[map->icjP[ii]]);}
          if (2*r > 1) {p1 = 1;} else {p1 = 2*r;}
          if (r > 2) {p2 = 2;} else {p2 = r;}
          if (p1 > 0) {phi = p1;} else {phi = 0;}
          if (p2 > phi) {phi = p2;}
          Sr = (*data)->S[map->icjP[ii]] + 0.5*phi*((*data)->S[ii] - (*data)->S[map->icjP[ii]])*Cu*(1-Cu);
        }
      }
      // advective flux at YM face
      Cu = (*data)->vvYP[map->icjM[ii]] * setting->dt / setting->dx;
      if ((*data)->Fvv[map->icjM[ii]] > 0)
      {
        if ((*data)->Fvv[map->icjMM[ii]] <= 0)
        {Sl = (*data)->S[map->icjM[ii]];}
        else
        {
          if (fabs((*data)->S[ii] - (*data)->S[map->icjM[ii]]) < dSmin)
          {r = 0;}
          else
          {r = ((*data)->S[map->icjM[ii]] - (*data)->S[map->icjMM[ii]])/((*data)->S[ii] - (*data)->S[map->icjM[ii]]);}
          if (2*r > 1) {p1 = 1;} else {p1 = 2*r;}
          if (r > 2) {p2 = 2;} else {p2 = r;}
          if (p1 > 0) {phi = p1;} else {phi = 0;}
          if (p2 > phi) {phi = p2;}
          Sl = (*data)->S[map->icjM[ii]] + 0.5*phi*((*data)->S[ii] - (*data)->S[map->icjM[ii]])*Cu*(1-Cu);
        }
      }
      else
      {
        if ((*data)->Fvv[ii] >= 0)
        {Sl = (*data)->S[ii];}
        else
        {
          if (fabs((*data)->S[ii] - (*data)->S[map->icjM[ii]]) < dSmin)
          {r = 0;}
          else
          {r = ((*data)->S[ii] - (*data)->S[map->icjP[ii]])/((*data)->S[map->icjM[ii]] - (*data)->S[ii]);}
          if (2*r > 1) {p1 = 1;} else {p1 = 2*r;}
          if (r > 2) {p2 = 2;} else {p2 = r;}
          if (p1 > 0) {phi = p1;} else {phi = 0;}
          if (p2 > phi) {phi = p2;}
          Sl = (*data)->S[ii] + 0.5*phi*((*data)->S[map->icjM[ii]] - (*data)->S[ii])*Cu*(1-Cu);
        }
      }
      // compute the scalar flux
      (*data)->Sm[ii] = (*data)->Sm[ii] + setting->dt * \
        (-(*data)->Fuu[ii] * Sd + (*data)->Fuu[map->iMjc[ii]] * Su - \
          (*data)->Fvv[ii] * Sr + (*data)->Fvv[map->icjM[ii]] * Sl);
    }
  }
}

// ==================== diffusive scalar flux ====================
void scalarDiffusion(Data **data, Maps *map, Config *setting)
{
  int ii;
	double tKA, dz;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    (*data)->Sm[ii] = (*data)->Sm[ii] + setting->dt * \
      ((setting->Kx * setting->dy * (*data)->depthXP[ii] / setting->dx) * \
      ((*data)->S[map->iPjc[ii]] - (*data)->S[ii]) - \
      (setting->Kx * setting->dy * (*data)->depthXP[map->iMjc[ii]] / setting->dx) * \
      ((*data)->S[ii] - (*data)->S[map->iMjc[ii]]) + \
      (setting->Ky * setting->dx * (*data)->depthYP[ii] / setting->dy) * \
      ((*data)->S[map->icjP[ii]] - (*data)->S[ii]) - \
      (setting->Ky * setting->dx * (*data)->depthYP[map->icjM[ii]] / setting->dy) * \
      ((*data)->S[ii] - (*data)->S[map->icjM[ii]]));
		// subsurface scalar, ZhiLi20180112
		if (setting->useSubsurface == 1 & setting->useSubscalar == 1)
		{
			if ((*data)->dz[ii] > 0)
			{
				tKA = setting->dt * setting->Kmz * setting->dx * setting->dy * setting->porosity;
				(*data)->Sm[ii] = (*data)->Sm[ii] + tKA * ((*data)->Ssub[ii] - (*data)->S[ii]) / (*data)->dz[ii] ;
			}
		}
  }
}

// ================ update scalar concentration ===============
void updateScalar(Data **data, BC *bc, Bath *bath, Maps *map, Config *setting, int tt, int irank, int nrank)
{
  int ii, jj, *N, flag = 0;
  int iP = 0, iM = 0, jP = 0, jM = 0, kM = 0;
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
    if ((*data)->depthXP[ii] > 0 & dXPold > 0) {N[ii] += 1; iP = 1; dXPold = 0;}

    dXMold = (*data)->surfOld[map->iMjc[ii]] - bath->bottomZ[map->iMjc[ii]];
    if ((*data)->depthXP[map->iMjc[ii]] > 0 & dXMold > 0) {N[ii] += 1; iM = 1; dXMold = 0;}

    if (irank != nrank - 1)
    {
      dYPold = (*data)->surfOld[map->icjP[ii]] - bath->bottomZ[map->icjP[ii]];
      if ((*data)->depthYP[ii] > 0 & dYPold > 0) {N[ii] += 1; jP = 1; dYPold = 0;}
    }
    else
    {
      if (map->icjP[ii] <= setting->N2ci)
      {
        dYPold = (*data)->surfOld[map->icjP[ii]] - bath->bottomZ[map->icjP[ii]];
        if ((*data)->depthYP[ii] > 0 & dYPold > 0) {N[ii] += 1; jP = 1; dYPold = 0;}
      }
    }

    if (irank != 0)
    {
      dYMold = (*data)->surfOld[map->icjM[ii]] - bath->bottomZ[map->icjM[ii]];
      if ((*data)->depthYP[map->icjM[ii]] > 0 & dYMold > 0) {N[ii] += 1; jM = 1; dYMold = 0;}
    }
    else
    {
      if (map->icjM[ii] <= setting->N2ci)
      {
        dYMold = (*data)->surfOld[map->icjM[ii]] - bath->bottomZ[map->icjM[ii]];
        if ((*data)->depthYP[map->icjM[ii]] > 0 & dYMold > 0) {N[ii] += 1; jM = 1; dYMold = 0;}
      }
    }
	// connection to subsurface, ZhiLi20180112
	if ((*data)->dz[ii] > 0)
	{N[ii] += 1;	kM = 1;}
	// find max and min scalar of neighbor cells
    if (N[ii] != 0)
    {
      Sarr = malloc(N[ii] * sizeof(double));
      if (iP == 1) {Sarr[jj] = (*data)->S[map->iPjc[ii]]; jj++;}
      if (iM == 1) {Sarr[jj] = (*data)->S[map->iMjc[ii]]; jj++;}
      if (jP == 1) {Sarr[jj] = (*data)->S[map->icjP[ii]]; jj++;}
      if (jM == 1) {Sarr[jj] = (*data)->S[map->icjM[ii]]; jj++;}
      if (kM == 1) {Sarr[jj] = (*data)->Ssub[ii]; jj++;}
      Smax[ii] = getMax(Sarr, N[ii]);
      Smin[ii] = getMin(Sarr, N[ii]);
      iP = 0;
      iM = 0;
      jP = 0;
      jM = 0;
	  kM = 0;
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
	(*data)->cellV[ii] = (*data)->cellV[ii] + (*data)->Qseep[ii];
    
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

    if (N[ii] != 0 & (*data)->S[ii] > Smax[ii] & flag == 0)
	{
		//if (setting->useEvap == 0 | bc->evap[tt] == 0.0)
    	{(*data)->S[ii] = Smax[ii];}
	}
    else if (N[ii] != 0 & (*data)->S[ii] < Smin[ii] & flag == 0)
	{
		if (setting->useRain == 0 | (*data)->rainC[ii] >= 0.0)
    	{(*data)->S[ii] = Smin[ii];}
//        else if ((*data)->rainC[ii] < 0.0)
//        {(*data)->rainC[ii] = 0.0;}
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

// =============== enforce boundary condition ===============
void enforceScalarBC(Data **data, BC *bc, Maps *map, Config *setting, int tt, \
  int irank, int nrank)
{
  int ii;
  // Neumann bc at all boundaries
  for (ii = 0; ii < setting->nx; ii++)
  {
    (*data)->S[map->jPgt[ii]] = (*data)->S[map->jPbd[ii]];
    (*data)->S[map->jMgt[ii]] = (*data)->S[map->jMbd[ii]];
  }
  for (ii = 0; ii < setting->ny; ii++)
  {
    (*data)->S[map->iPgt[ii]] = (*data)->S[map->iPbd[ii]];
    (*data)->S[map->iMgt[ii]] = (*data)->S[map->iMbd[ii]];
  }
  // Direchlet bc at YP tidal boundaries
  if (irank == nrank - 1 & setting->bcType != 2)
  {
    if (setting->useConstTidePS == 1)
    {
      for (ii = 0; ii < setting->tideLocLengthP; ii++)
      {(*data)->S[setting->tideLocP[ii]] = setting->tidePS;}
    }
    else
    {
      for (ii = 0; ii < setting->tideLocLengthP; ii++)
      {(*data)->S[setting->tideLocP[ii]] = bc->tidalPS[tt];}
    }
  }
  // Direchlet bc at YM tidal boundaries
  if (irank == 0 & setting->bcType != 1)
  {
    if (setting->useConstTideMS == 1)
    {
      for (ii = 0; ii < setting->tideLocLengthM; ii++)
      {(*data)->S[setting->tideLocM[ii]] = setting->tideMS;}
    }
    else
    {
      for (ii = 0; ii < setting->tideLocLengthM; ii++)
      {(*data)->S[setting->tideLocM[ii]] = bc->tidalMS[tt];}
    }
  }
  // Direchlet bc at inflow boundaries
  //for (ii = 0; ii < setting->inflowLocLength; ii++)
  //{(*data)->S[setting->inflowLoc[ii]] = setting->inflowS;}
  // zero scalar over dry lands
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    if ((*data)->depth[ii] <= 0)
    {(*data)->S[ii] = 0;}
  }
}

// ====================== Check if subtime stepping is required =======================
void scalarSubStepping(Data **data, Bath *bath, BC *bc, Maps *map, Config *setting, int tt, int irank, int nrank)
{
    int ii, jj, kk, Nsub;
    double CFLx0 = 0.0, CFLy0 = 0.0, dt = 0.0;
    // global subtime stepping
    for (ii = 0; ii < setting->N2ci; ii++)
    {
        if ((*data)->CFLx[ii] > CFLx0)
        {CFLx0 = (*data)->CFLx[ii];}
        if ((*data)->CFLy[ii] > CFLy0)
        {CFLy0 = (*data)->CFLy[ii];}
    }
    //printf("Max CFL in x, y directions are %f, %f\n",CFLx0,CFLy0);
    if (CFLx0 > CFLy0) {Nsub = ceil(CFLx0);}
    else {Nsub = ceil(CFLy0);}
    if (Nsub > 1.0)
    {
        dt = setting->dt * 1.0 / Nsub;
        printf("Subtime-stepping is invoked with %d substeps...\n", Nsub);
        for (kk = 1; kk <= Nsub; kk++)
        {
            for (ii = 0; ii < setting->N2ci; ii++)
            {
                // subtime advection
                double Su, Sd, Sr, Sl;
                if ((*data)->Fuu[ii] > 0) {Sd = (*data)->S[ii];} else {Sd = (*data)->S[map->iPjc[ii]];}
                if ((*data)->Fuu[map->iMjc[ii]] > 0) {Su = (*data)->S[map->iMjc[ii]];} else {Su = (*data)->S[ii];}
                if ((*data)->Fvv[ii] > 0) {Sr = (*data)->S[ii];} else {Sr = (*data)->S[map->icjP[ii]];}
                if ((*data)->Fvv[map->icjM[ii]] > 0) {Sl = (*data)->S[map->icjM[ii]];} else {Sl = (*data)->S[ii];}
                (*data)->Sm[ii] = (*data)->Sm[ii] + dt * \
                (-(*data)->Fuu[ii] * Sd + (*data)->Fuu[map->iMjc[ii]] * Su - \
                 (*data)->Fvv[ii] * Sr + (*data)->Fvv[map->icjM[ii]] * Sl);
                // subtime diffusion
                (*data)->Sm[ii] = (*data)->Sm[ii] + dt * \
                ((setting->Kx * setting->dy * (*data)->depthXP[ii] / setting->dx) * \
                 ((*data)->S[map->iPjc[ii]] - (*data)->S[ii]) - \
                 (setting->Kx * setting->dy * (*data)->depthXP[map->iMjc[ii]] / setting->dx) * \
                 ((*data)->S[ii] - (*data)->S[map->iMjc[ii]]) + \
                 (setting->Ky * setting->dx * (*data)->depthYP[ii] / setting->dy) * \
                 ((*data)->S[map->icjP[ii]] - (*data)->S[ii]) - \
                 (setting->Ky * setting->dx * (*data)->depthYP[map->icjM[ii]] / setting->dy) * \
                 ((*data)->S[ii] - (*data)->S[map->icjM[ii]]));
                // subtime inflow
                if (irank == 0 & setting->bcType != 3)
                {
                    for (jj = 0; jj < setting->inflowLocLength; jj++)
                    {(*data)->Sm[setting->inflowLoc[jj]] = (*data)->Sm[setting->inflowLoc[jj]] + \
                        setting->inflowS * bc->inflow[tt] * dt / setting->inflowLocLength;}
                }
                // update volume
                (*data)->cellV[ii] = (*data)->cellV[ii] + dt * \
                ((*data)->Fuu[map->iMjc[ii]] - (*data)->Fuu[ii] + \
                 (*data)->Fvv[map->icjM[ii]] - (*data)->Fvv[ii]);
                if (irank == 0 & setting->bcType != 3)
                {
                    for (jj = 0; jj < setting->inflowLocLength; jj++)
                    {if (ii == setting->inflowLoc[jj])
                    {(*data)->cellV[ii] = (*data)->cellV[ii] + bc->inflow[tt]*dt;}}
                }
                // update scalar concentration
                if ((*data)->cellV[ii] > 0)
                {(*data)->S[ii] = (*data)->Sm[ii] / (*data)->cellV[ii];}
            }
        }
    }
    else
    {
        scalarAdvection(data, map, setting);
        scalarDiffusion(data, map, setting);
        updateScalar(data, bc, bath, map, setting, tt, irank, nrank);
    }
    
}










