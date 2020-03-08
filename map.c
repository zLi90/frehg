#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

// -----------------------------------------------------------------------------
// Functions defining the grid maps of the model
// - 2017-05-04 by Zhi Li -
// -----------------------------------------------------------------------------

#include "configuration.h"
#include "map.h"

void createMaps(Maps **map, Config *setting)
{
  int ii, jj, kk, col;
  *map = malloc(sizeof(Maps));
  (*map)->cntr = malloc(setting->N2ci*sizeof(int));
  (*map)->trps = malloc(setting->N2ci*sizeof(int));
  (*map)->sprt = malloc(setting->N2ci*sizeof(int));
  (*map)->iPjc = malloc(setting->N2ci*sizeof(int));
  (*map)->iMjc = malloc(setting->N2ci*sizeof(int));
  (*map)->icjP = malloc(setting->N2ci*sizeof(int));
  (*map)->icjM = malloc(setting->N2ci*sizeof(int));
  (*map)->iPPjc = malloc(setting->N2ci*sizeof(int));
  (*map)->iMMjc = malloc(setting->N2ci*sizeof(int));
  (*map)->icjPP = malloc(setting->N2ci*sizeof(int));
  (*map)->icjMM = malloc(setting->N2ci*sizeof(int));
  (*map)->iPjP = malloc(1*sizeof(int));
  (*map)->iPjM = malloc(1*sizeof(int));
  (*map)->iMjP = malloc(1*sizeof(int));
  (*map)->iMjM = malloc(1*sizeof(int));
  (*map)->ii2d = malloc(setting->N2ci*sizeof(int));
  (*map)->jj2d = malloc(setting->N2ci*sizeof(int));
  (*map)->iPbd = malloc(setting->ny*sizeof(int));
  (*map)->iPgt = malloc(setting->ny*sizeof(int));
  (*map)->iMbd = malloc(setting->ny*sizeof(int));
  (*map)->iMgt = malloc(setting->ny*sizeof(int));
  (*map)->jPbd = malloc(setting->nx*sizeof(int));
  (*map)->jPgt = malloc(setting->nx*sizeof(int));
  (*map)->jMbd = malloc(setting->nx*sizeof(int));
  (*map)->jMgt = malloc(setting->nx*sizeof(int));
  // set index for the center map
  for (ii = 0; ii < setting->N2ci; ii++)
  {(*map)->cntr[ii] = ii;}
  // set index for the transposed map
  for (kk = 0; kk < setting->N2ci; kk++)
  {
    ii = floor(kk / setting->ny);
    jj = kk % setting->ny;
    (*map)->trps[kk] = jj*setting->nx + ii;
    ii = floor(kk / setting->nx);
    jj = kk % setting->nx;
    (*map)->sprt[kk] = jj*setting->ny + ii;
  }
  // set index for the jP and jM maps
  for (ii = 0; ii < setting->N2ci-setting->nx; ii++)
  {
    (*map)->icjP[ii] = ii + setting->nx;
    (*map)->icjM[ii+setting->nx] = ii;
  }
  for (ii = 0; ii < setting->nx; ii++)
  {
    (*map)->icjP[ii+setting->N2ci-setting->nx] = ii + setting->N2ci;
    (*map)->icjM[ii] = ii + setting->N2ci + setting->nx;
  }
  // set index for the jPP and jMM maps, ZhiLi20180504
  for (ii = 0; ii < setting->N2ci-setting->nx; ii++)
  {
    (*map)->icjPP[ii] = (*map)->icjP[(*map)->icjP[ii]];
    (*map)->icjMM[ii+2*setting->nx] = (*map)->icjM[(*map)->icjM[ii+2*setting->nx]];
  }
  for (ii = 0; ii < setting->nx; ii++)
  {
    (*map)->icjPP[ii+setting->N2ci-setting->nx] = (*map)->icjP[ii+setting->N2ci-setting->nx];
    (*map)->icjMM[ii] = (*map)->icjM[ii];
  }
  // set index for the iP and iM maps
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    col = floor(ii/setting->nx);
    if (ii % setting->nx == setting->nx-1)
    {(*map)->iPjc[ii] = setting->N2ci + 2*setting->nx + col; (*map)->iMjc[ii] = ii - 1;}
    else if (ii % setting->nx == 0)
    {(*map)->iMjc[ii] = setting->N2ci + 2*setting->nx + setting->ny + col; (*map)->iPjc[ii] = ii + 1;}
    else
    {(*map)->iPjc[ii] = ii + 1;   (*map)->iMjc[ii] = ii - 1;}
  }
  // set index for the iPP and iMM maps, ZhiLi20180504
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    if (ii % setting->nx == setting->nx-1)
    {(*map)->iPPjc[ii] = (*map)->iPjc[ii]; (*map)->iMMjc[ii] = (*map)->iMjc[(*map)->iMjc[ii]];}
    else if (ii % setting->nx == 0)
    {(*map)->iMMjc[ii] = (*map)->iMjc[ii]; (*map)->iPPjc[ii] = (*map)->iPjc[(*map)->iPjc[ii]];}
    else
    {(*map)->iPPjc[ii] = (*map)->iPjc[(*map)->iPjc[ii]];   (*map)->iMMjc[ii] = (*map)->iMjc[(*map)->iMjc[ii]];}
  }
  // set corner maps
  (*map)->iMjM[0] = setting->N2ct - 1;
  (*map)->iPjM[0] = setting->N2ct - 2;
  (*map)->iPjP[0] = setting->N2ct - 3;
  (*map)->iMjP[0] = setting->N2ct - 4;
  // set the 2d (ii,jj) maps
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    col = floor(ii/setting->nx);
    (*map)->ii2d[ii] = ii - col*setting->nx;
    (*map)->jj2d[ii] = col;
  }
  // set the maps for the jP and jM ghost cells
  for (ii = 0; ii < setting->nx; ii++)
  {
    (*map)->jPbd[ii] = (*map)->cntr[setting->N2ci-setting->nx+ii];
    (*map)->jMbd[ii] = (*map)->cntr[ii];
    (*map)->jPgt[ii] = (*map)->icjP[setting->N2ci-setting->nx+ii];
    (*map)->jMgt[ii] = (*map)->icjM[ii];
  }
  // set the maps for the iP and iM ghost cells
  for (ii = 0; ii < setting->ny; ii++)
  {
    (*map)->iPbd[ii] = (*map)->cntr[((ii+1)*setting->nx)-1];
    (*map)->iMbd[ii] = (*map)->cntr[ii*setting->nx];
    (*map)->iPgt[ii] = (*map)->iPjc[((ii+1)*setting->nx)-1];
    (*map)->iMgt[ii] = (*map)->iMjc[ii*setting->nx];
  }
//    for (ii = 0; ii < setting->N2ci; ii++)
//    {printf("map ii = %d\n",(*map)->cntr[ii]);}
}


// ========== Create Map for Subsurface ==========
void createGmaps(Gmaps **gmap, Bath *bath, Maps *map, Config *setting)
{
    int ii, jj, kk, ll, pp, ind, maxLay, N3ca = 0;
    double columnH, nnlay, htopmin = 0.01;
    int ngtiP = 0, ngtiM = 0, ngtjP = 0, ngtjM = 0;
    int indgtiP = 0, indgtiM = 0, indgtjP = 0, indgtjM = 0;
    // initialize gmap, number of layer in each column, and top-layer thickness
    *gmap = malloc(sizeof(Gmaps));
    (*gmap)->nlay = malloc(setting->N2ci*sizeof(int));
    (*gmap)->htop = malloc(setting->N2ci*sizeof(double));
    for (ii = 0; ii < setting->N2ci; ii++)
    {
        (*gmap)->nlay[ii] = 0;
        (*gmap)->htop[ii] = -1;
    }
    // compute total number of subsurface cells
    maxLay = 0;
    for (ll = 0; ll < setting->N2ci; ll++)
    {
        columnH = bath->bottomZ[ll] - setting->zbot;
        if (columnH < 0) {printf("WARNING: Column height is negative!\n");}

        if (floor(columnH/setting->layZ) == columnH/setting->layZ)
        {(*gmap)->nlay[ll] = floor(columnH / setting->layZ);}
        else
        {
            if (columnH - floor(columnH/setting->layZ)*setting->layZ > 0.5*setting->layZ)
            {(*gmap)->nlay[ll] = floor(columnH/setting->layZ) + 1;}
            else
            {(*gmap)->nlay[ll] = floor(columnH/setting->layZ);}
        }
        (*gmap)->htop[ll] = setting->layZ;
        if ((*gmap)->htop[ll] < htopmin) {(*gmap)->htop[ll] = htopmin;}
        if ((*gmap)->nlay[ll] > maxLay) {maxLay = (*gmap)->nlay[ll];}
        N3ca += (*gmap)->nlay[ll];
    }
	(*gmap)->maxLay = maxLay;
    printf("Max number of layers = %d\n",maxLay);
    (*gmap)->ngtjP = (*gmap)->maxLay * setting->nx;
	(*gmap)->ngtjM = (*gmap)->maxLay * setting->nx;
	(*gmap)->ngtiP = (*gmap)->maxLay * setting->ny;
	(*gmap)->ngtiM = (*gmap)->maxLay * setting->ny;
    // compute total number of subsurface cells
    //for (ii = 0; ii < setting->nx; ii++)
    //{ngtjP += (*gmap)->nlay[map->jPbd[ii]];    ngtjM += (*gmap)->nlay[map->jMbd[ii]];}
    //for (ii = 0; ii < setting->ny; ii++)
    //{ngtiP += (*gmap)->nlay[map->iPbd[ii]];    ngtiM += (*gmap)->nlay[map->iMbd[ii]];}
	// N3ci = all interior grids for subsurface domain
    setting->N3ci = (*gmap)->maxLay * setting->nx * setting->ny;
	// N3ct = all grids for subsurface domain (interior + ghost)
    setting->N3ct = setting->N3ci + (*gmap)->ngtjP + (*gmap)->ngtjM + (*gmap)->ngtiP + (*gmap)->ngtiM;
	// N3CI = all interior grids for all ranks
	setting->N3CI = setting->N3ci * setting->np;
	// N3cf = all subsurface + surface grids for the current rank
    setting->N3cf = setting->N3ct + setting->N2ci;
	// N3ca = all active interior subsurface grids
	setting->N3ca = N3ca;

    // initialize Gmap fields
    (*gmap)->cntr = malloc(setting->N3ci*sizeof(int));
    (*gmap)->actv = malloc(setting->N3ci*sizeof(int));
    (*gmap)->iPjckc = malloc(setting->N3ci*sizeof(int));
    (*gmap)->iMjckc = malloc(setting->N3ci*sizeof(int));
    (*gmap)->icjPkc = malloc(setting->N3ci*sizeof(int));
    (*gmap)->icjMkc = malloc(setting->N3ci*sizeof(int));
    (*gmap)->icjckP = malloc(setting->N3ci*sizeof(int));
    (*gmap)->icjckM = malloc(setting->N3ci*sizeof(int));
    (*gmap)->ii = malloc(setting->N3ci*sizeof(int));
    (*gmap)->jj = malloc(setting->N3ci*sizeof(int));
    (*gmap)->kk = malloc(setting->N3ci*sizeof(int));
    (*gmap)->istop = malloc(setting->N3ci*sizeof(int));
    (*gmap)->top2D = malloc(setting->N3ci*sizeof(int));
    (*gmap)->dz3d = malloc(setting->N3ci*sizeof(double));
    (*gmap)->bot3d = malloc(setting->N3ci*sizeof(double));
    (*gmap)->bot2d = malloc(setting->N2ci*sizeof(double));
    // outer loop over 2D domain
    ind = 0;
    for (ll = 0; ll < setting->N2ci; ll++)
    {
        jj = floor(ll / setting->nx);
        ii = ll - jj * setting->nx;
        (*gmap)->bot2d[ll] = bath->bottomZ[ll];
        // inner loop over each column
        for (kk = 0; kk < (*gmap)->maxLay; kk++)
        {
            (*gmap)->ii[ind] = ii;
            (*gmap)->jj[ind] = jj;
            (*gmap)->kk[ind] = kk;
            (*gmap)->cntr[ind] = ind;
			if (kk < (*gmap)->maxLay - (*gmap)->nlay[ll])
			{(*gmap)->actv[ind] = 0;}
			else
			{(*gmap)->actv[ind] = 1;}
            // iP
            if (ii != setting->nx-1)
            {(*gmap)->iPjckc[ind] = ind + (*gmap)->maxLay;}
            else
            {
				(*gmap)->iPjckc[ind] = setting->N3ci + (*gmap)->ngtjP + (*gmap)->ngtjM + indgtiP;
				indgtiP += 1;
			}
            // iM
            if (ii != 0)
            {(*gmap)->iMjckc[ind] = ind - (*gmap)->maxLay;}
            else
            {
				(*gmap)->iMjckc[ind] = setting->N3ci + (*gmap)->ngtjP + (*gmap)->ngtjM + (*gmap)->ngtiP + indgtiM;
				indgtiM += 1;
			}
            // jP
            if (jj != setting->ny-1)
            {(*gmap)->icjPkc[ind] = ind + setting->nx * (*gmap)->maxLay;}
            else
            {
				(*gmap)->icjPkc[ind] = setting->N3ci + indgtjP;
				indgtjP += 1;
			}
            // jM
            if (jj != 0)
            {(*gmap)->icjMkc[ind] = ind - setting->nx * (*gmap)->maxLay;}
            else
            {
				(*gmap)->icjMkc[ind] = setting->N3ci + (*gmap)->ngtjP + indgtjM;
				indgtjM += 1;
			}
            // kP
            if (kk != (*gmap)->maxLay-1 & kk >= (*gmap)->maxLay - (*gmap)->nlay[ll])
            {(*gmap)->icjckP[ind] = ind + 1;}
            else
            {(*gmap)->icjckP[ind] = -1;}
            // kM
            if (kk != 0 & kk > (*gmap)->maxLay - (*gmap)->nlay[ll])
            {(*gmap)->icjckM[ind] = ind - 1;}
			else if (kk == (*gmap)->maxLay - (*gmap)->nlay[ll])
			{(*gmap)->icjckM[ind] = setting->N3ct + ll;}
            else
            {(*gmap)->icjckM[ind] = -1;}

			// identify if a cell is on top layer
            if (kk == (*gmap)->maxLay - (*gmap)->nlay[ll])
            {(*gmap)->istop[ind] = 1;   (*gmap)->top2D[ind] = ll;}
            else if (kk == (*gmap)->maxLay - (*gmap)->nlay[ll] + 1)
            {(*gmap)->istop[ind] = 2;   (*gmap)->top2D[ind] = ll;}
            else
            {(*gmap)->istop[ind] = 0;   (*gmap)->top2D[ind] = ll;}
            // calculate dz vector
            if ((*gmap)->actv[ind] == 1)
            {
                (*gmap)->dz3d[ind] = setting->layZ;
                // if ((*gmap)->istop[ind] == 1)
                // {(*gmap)->dz3d[ind] = (*gmap)->htop[(*gmap)->top2D[ind]];}
                // else
                // {(*gmap)->dz3d[ind] = setting->layZ;}
            }
            else
            {(*gmap)->dz3d[ind] = -1;}
            // update bottom elevation of each active cell
            if ((*gmap)->actv[ind] == 1)
            {
                if ((*gmap)->istop[ind] == 1)
                {(*gmap)->bot3d[ind] = bath->bottomZ[(*gmap)->top2D[ind]] - (*gmap)->dz3d[ind];}
                else
                {(*gmap)->bot3d[ind] = (*gmap)->bot3d[(*gmap)->icjckM[ind]] - setting->layZ;}
            }
            else
            {(*gmap)->bot3d[ind] = 100;}
            // update the counter
            ind++;

        }
    }
}
