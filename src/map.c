// Build maps (connections) between grid cells
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

// -----------------------------------------------------------------------------
// Functions defining the grid maps of the model
// - 2017-05-04 by Zhi Li -
// -----------------------------------------------------------------------------

#include"configuration.h"
#include"initialize.h"
#include"map.h"
#include"mpifunctions.h"
#include"utility.h"

void build_surf_map(Map **map, Config *param);
void build_surf_mat_map(double *actv, Map **map, Config *param);
void build_subsurf_map(Map **map, Map *smap, double *bath, double *offset, Config *param, int irank);
void build_subs_mat_map(double *actv, Map **map, Config *param);

// >>>>> Build connections for surface domain <<<<<
void build_surf_map(Map **map, Config *param)
{
    int ii;
    *map = malloc(sizeof(Map));
    (*map)->cntr = malloc(param->n2ci*sizeof(int));
    (*map)->ii = malloc(param->n2ci*sizeof(int));
    (*map)->jj = malloc(param->n2ci*sizeof(int));
    (*map)->iPjc = malloc(param->n2ci*sizeof(int));
    (*map)->iMjc = malloc(param->n2ci*sizeof(int));
    (*map)->icjP = malloc(param->n2ci*sizeof(int));
    (*map)->icjM = malloc(param->n2ci*sizeof(int));
    (*map)->iPjP = malloc(param->n2ci*sizeof(int));
    (*map)->iMjM = malloc(param->n2ci*sizeof(int));
    (*map)->iPjM = malloc(param->n2ci*sizeof(int));
    (*map)->iMjP = malloc(param->n2ci*sizeof(int));
    (*map)->iPin = malloc(param->ny*sizeof(int));
    (*map)->iPou = malloc(param->ny*sizeof(int));
    (*map)->iMin = malloc(param->ny*sizeof(int));
    (*map)->iMou = malloc(param->ny*sizeof(int));
    (*map)->jPin = malloc(param->nx*sizeof(int));
    (*map)->jPou = malloc(param->nx*sizeof(int));
    (*map)->jMin = malloc(param->nx*sizeof(int));
    (*map)->jMou = malloc(param->nx*sizeof(int));
    (*map)->ibot = malloc(param->n2ci*sizeof(int));
    (*map)->dz = malloc(param->n2ci*sizeof(double));

    // set map indexes
    for (ii = 0; ii < param->n2ci; ii++)
    {
        (*map)->cntr[ii] = ii;
        (*map)->jj[ii] = floor(ii/param->nx);
        (*map)->ii[ii] = ii - (*map)->jj[ii] * param->nx;
        // iP map
        if ((*map)->ii[ii] == param->nx-1)
        {(*map)->iPjc[ii] = param->n2ci + 2*param->nx + (*map)->jj[ii];}
        else
        {(*map)->iPjc[ii] = ii + 1;}
        // iM map
        if ((*map)->ii[ii] == 0)
        {(*map)->iMjc[ii] = param->n2ci + 2*param->nx + param->ny + (*map)->jj[ii];}
        else
        {(*map)->iMjc[ii] = ii - 1;}
        // jP map
        if ((*map)->jj[ii] == param->ny-1)
        {(*map)->icjP[ii] = param->n2ci + (*map)->ii[ii];}
        else
        {(*map)->icjP[ii] = ii + param->nx;}
        // jM map
        if ((*map)->jj[ii] == 0)
        {(*map)->icjM[ii] = param->n2ci + param->nx + (*map)->ii[ii];}
        else
        {(*map)->icjM[ii] = ii - param->nx;}
    }

    // diagonal maps
    for (ii = 0; ii < param->n2ci; ii++)
    {
        // iPjP map
        if ((*map)->ii[ii] == param->nx-1)
        {
            if ((*map)->jj[ii] == param->ny-1)
            {(*map)->iPjP[ii] = param->n2ct - 3;}
            else
            {(*map)->iPjP[ii] = (*map)->iPjc[(*map)->icjP[ii]];}
        }
        else if ((*map)->jj[ii] == param->ny-1)
        {(*map)->iPjP[ii] = (*map)->icjP[ii] + 1;}
        else
        {(*map)->iPjP[ii] = (*map)->iPjc[(*map)->icjP[ii]];}
        // iMjM map
        if ((*map)->ii[ii] == 0)
        {
            if ((*map)->jj[ii] == 0)
            {(*map)->iMjM[ii] = param->n2ct - 1;}
            else
            {(*map)->iMjM[ii] = (*map)->iMjc[(*map)->icjM[ii]];}
        }
        else if ((*map)->jj[ii] == 0)
        {(*map)->iMjM[ii] = (*map)->icjM[ii] - 1;}
        else
        {(*map)->iMjM[ii] = (*map)->iMjc[(*map)->icjM[ii]];}
        // iPjM map
        if ((*map)->ii[ii] == param->nx-1)
        {
            if ((*map)->jj[ii] == 0)
            {(*map)->iPjM[ii] = param->n2ct - 2;}
            else
            {(*map)->iPjM[ii] = (*map)->iPjc[(*map)->icjM[ii]];}
        }
        else if ((*map)->jj[ii] == 0)
        {(*map)->iPjM[ii] = (*map)->icjM[ii] + 1;}
        else
        {(*map)->iPjM[ii] = (*map)->iPjc[(*map)->icjM[ii]];}
        // iMjP map
        if ((*map)->ii[ii] == 0)
        {
            if ((*map)->jj[ii] == param->ny-1)
            {(*map)->iMjP[ii] = param->n2ct - 4;}
            else
            {(*map)->iMjP[ii] = (*map)->iMjc[(*map)->icjP[ii]];}
        }
        else if ((*map)->jj[ii] == param->ny-1)
        {(*map)->iMjP[ii] = (*map)->icjP[ii] - 1;}
        else
        {(*map)->iMjP[ii] = (*map)->iMjc[(*map)->icjP[ii]];}
    }
    // boundary cells for mpi exchange
    for (ii = 0; ii < param->nx; ii++)
    {
        (*map)->jPin[ii] = param->nx*(param->ny-1) + ii;
        (*map)->jPou[ii] = (*map)->icjP[(*map)->jPin[ii]];
        (*map)->jMin[ii] = ii;
        (*map)->jMou[ii] = (*map)->icjM[(*map)->jMin[ii]];
    }
    for (ii = 0; ii < param->ny; ii++)
    {
        (*map)->iPin[ii] = (ii+1) * param->nx - 1;
        (*map)->iPou[ii] = (*map)->iPjc[(*map)->iPin[ii]];
        (*map)->iMin[ii] = ii * param->nx;
        (*map)->iMou[ii] = (*map)->iMjc[(*map)->iMin[ii]];
    }
}

// >>>>> Build connections for surface domain <<<<<
void build_surf_mat_map(double *actv, Map **map, Config *param)   {
    int ii, jj, ixm, ixp, iym, iyp, irow = 0;
    (*map)->dom2mat = malloc(param->n2ci*sizeof(int));
    (*map)->mat2dom = malloc(param->nactv*sizeof(int));
    (*map)->distxp = malloc(param->nactv*sizeof(int));
    (*map)->distxm = malloc(param->nactv*sizeof(int));
    (*map)->distyp = malloc(param->nactv*sizeof(int));
    (*map)->distym = malloc(param->nactv*sizeof(int));
    for (ii = 0; ii < param->n2ci; ii++)    {
        if (actv[ii] > 0.0)    {
            (*map)->dom2mat[ii] = irow;
            (*map)->mat2dom[irow] = ii;
            irow += 1;
        }
        else {(*map)->dom2mat[ii] = -1;}
    }
    // compute distance in the matrix system
    for (ii = 0; ii < param->nactv; ii++)   {
        jj = (*map)->mat2dom[ii];
        //xp
        if (actv[(*map)->iPjc[jj]] > 0)    {
            ixp = (*map)->dom2mat[(*map)->iPjc[jj]];
            (*map)->distxp[ii] = ixp - ii;
        }
        else {(*map)->distxp[ii] = -1;}
        //xm
        if (actv[(*map)->iMjc[jj]] > 0)    {
            ixm = (*map)->dom2mat[(*map)->iMjc[jj]];
            (*map)->distxm[ii] = ii - ixm;
        }
        else {(*map)->distxm[ii] = -1;}
        //yp
        if (actv[(*map)->icjP[jj]] > 0)    {
            iyp = (*map)->dom2mat[(*map)->icjP[jj]];
            (*map)->distyp[ii] = iyp - ii;
        }
        else {(*map)->distyp[ii] = -1;}
        //ym
        if (actv[(*map)->icjM[jj]] > 0)    {
            iym = (*map)->dom2mat[(*map)->icjM[jj]];
            (*map)->distym[ii] = ii - iym;
        }
        else {(*map)->distym[ii] = -1;}
    }
}

// >>>>> Build connections for subsurface domain <<<<<
void build_subsurf_map(Map **map, Map *smap, double *bath, double *offset, Config *param, int irank)
{
    int ii, jj, nz_upper, nz_lower;
    double *bath_min, *bath_max, *bath_max_global, *bath_max_arr, *bath_min_global, *bath_min_arr;
    double bot_new, dz_new, hdiff, dist;

    *map = malloc(sizeof(Map));
    bath_min = malloc(sizeof(double));
    bath_max = malloc(sizeof(double));
    bath_max_global = malloc(sizeof(double));
    bath_max_arr = malloc(param->mpi_ny*param->mpi_nx*sizeof(double));
    bath_min_global = malloc(sizeof(double));
    bath_min_arr = malloc(param->mpi_ny*param->mpi_nx*sizeof(double));
    // calculate number of layers
    bath_max[0] = getMax(bath, param->n2ci);
    bath_min[0] = getMin(bath, param->n2ci);
    if (param->use_mpi == 1)
    {
        mpi_gather_double(bath_max_arr, bath_max, 1, 0);
        bath_max_global[0] = getMax(bath_max_arr, param->mpi_ny*param->mpi_nx);
        mpi_bcast_double(bath_max_global, 1, 0);
        mpi_gather_double(bath_min_arr, bath_min, 1, 0);
        bath_min_global[0] = getMin(bath_min_arr, param->mpi_ny*param->mpi_nx);
        mpi_bcast_double(bath_min_global, 1, 0);
    }
    else
    {
        bath_max_global[0] = bath_max[0];
        bath_min_global[0] = bath_min[0];
    }
    param->botZ += offset[0];
    if (param->dz_incre == 1.0)
    {
        param->nz = ceil((bath_max_global[0] - param->botZ) / param->dz);
        // param->nz = ceil((bath_min_global[0] - param->botZ) / param->dz);
        (*map)->bot1d = malloc(param->nz*sizeof(double));
        for (ii = 0; ii < param->nz; ii++)
        // {(*map)->bot1d[ii] = bath_min_global[0] - (ii+1)*param->dz;}
        {(*map)->bot1d[ii] = bath_max_global[0] - (ii+1)*param->dz;}
        nz_upper = param->nz;
        nz_lower = 0;
    }
    else
    {
        bot_new = bath_max_global[0] - param->dz;
        dz_new = param->dz;
        param->nz = 1;
        while (bot_new > param->botZ)   {
            dz_new = param->dz_incre * dz_new;
            bot_new -= dz_new;
            param->nz += 1;
        }
        (*map)->bot1d = malloc(param->nz*sizeof(double));
        (*map)->bot1d[0] = bath_max_global[0] - param->dz;
        dz_new = param->dz;
        for (ii = 1; ii < param->nz; ii++)
        {
            dz_new = param->dz_incre * dz_new;
            (*map)->bot1d[ii] = (*map)->bot1d[ii-1] - dz_new;
        }
    }
    if (irank == 0) {printf("   >> Total number of subsurface layer = %d\n",param->nz);}

    param->N3CI = param->NX * param->NY * param->nz;

    if ((*map)->bot1d[param->nz-1] > bath_min_global[0])
    {mpi_print("WARNING: Bottom of subsurface domain > min bathymetry!",irank);}

    // update domain dimensions
    param->n3ci = param->n2ci * param->nz;
    param->n3ct = param->n2ct * (param->nz+2);

    if (param->sim_groundwater == 1)
    {
        // calculate 3D dz, actv map and cntr map
        (*map)->cntr = malloc(param->n3ci*sizeof(int));
        (*map)->zcntr = malloc(param->n3ci*sizeof(double));
        (*map)->zcntr_root = malloc(param->N3CI*sizeof(double));
        (*map)->zcntr_out = malloc(param->N3CI*sizeof(double));
        (*map)->ii = malloc(param->n3ci*sizeof(int));
        (*map)->jj = malloc(param->n3ci*sizeof(int));
        (*map)->kk = malloc(param->n3ci*sizeof(int));
        (*map)->actv = malloc(param->n3ct*sizeof(int));
        (*map)->bot3d = malloc(param->n3ct*sizeof(double));
        (*map)->dz3d = malloc(param->n3ct*sizeof(double));
        (*map)->istop = malloc(param->n3ci*sizeof(int));
        (*map)->top2d = malloc(param->n3ci*sizeof(int));
        (*map)->sinx = malloc(param->n3ct*sizeof(double));
        (*map)->siny = malloc(param->n3ct*sizeof(double));
        (*map)->cosx = malloc(param->n3ct*sizeof(double));
        (*map)->cosy = malloc(param->n3ct*sizeof(double));
        (*map)->sintx = malloc(param->n2ct*sizeof(double));
        (*map)->sinty = malloc(param->n2ct*sizeof(double));
        (*map)->costx = malloc(param->n2ct*sizeof(double));
        (*map)->costy = malloc(param->n2ct*sizeof(double));
        (*map)->Ax = malloc(param->n3ct*sizeof(double));
        (*map)->Ay = malloc(param->n3ct*sizeof(double));
        (*map)->Az = malloc(param->n3ct*sizeof(double));
        (*map)->V = malloc(param->n3ct*sizeof(double));

        for (ii = 0; ii < param->n3ct; ii++)    {(*map)->actv[ii] = 0;}
        for (ii = 0; ii < param->n3ci; ii++)
        {
            (*map)->cntr[ii] = ii;
            (*map)->top2d[ii] = floor(ii / param->nz);
            (*map)->ii[ii] = smap->ii[(*map)->top2d[ii]];
            (*map)->jj[ii] = smap->jj[(*map)->top2d[ii]];
            (*map)->kk[ii] = ii - (*map)->top2d[ii] * param->nz;
            (*map)->istop[ii] = 0;
        }
        for (ii = 0; ii < param->n3ci; ii++)
        {
            if (param->sim_shallowwater == 1)   {
                if ((*map)->kk[ii] == 0)    {
                    smap->ibot[(*map)->top2d[ii]] = ii;
                }
            }
        }

        // get dz
        if (param->follow_terrain == 1) {
            if (param->dz_incre == 1.0) {
                for (ii = 0; ii < param->n3ci; ii++)    {
                    (*map)->actv[ii] = 1;
                    // (*map)->dz3d[ii] = param->dz;
                    (*map)->dz3d[ii] = (bath[(*map)->top2d[ii]] - param->botZ) / param->nz;
                    (*map)->bot3d[ii] = bath[(*map)->top2d[ii]] - ((*map)->kk[ii]+1)*(*map)->dz3d[ii];

                    // if ((*map)->kk[ii] < param->nz-1)   {
                    //     (*map)->dz3d[ii] = param->dz;
                    //     // (*map)->dz3d[ii] = (bath[(*map)->top2d[ii]] - param->botZ) / param->nz;
                    //     (*map)->bot3d[ii] = bath[(*map)->top2d[ii]] - ((*map)->kk[ii]+1)*(*map)->dz3d[ii];
                    // }
                    // else {
                    //     (*map)->dz3d[ii] = (*map)->bot3d[ii-1] - param->botZ;
                    //     (*map)->bot3d[ii] = param->botZ;
                    // }
                    if ((*map)->kk[ii] == 0)    {(*map)->istop[ii] = 1;}
                }
            }
            else {
                double incre_sum = 0.0;
                for (ii = 0; ii < param->nz; ii++)  {
                    incre_sum += pow(param->dz_incre, ii);
                }
                for (ii = 0; ii < param->n3ci; ii++)    {
                    (*map)->actv[ii] = 1;
                    if ((*map)->kk[ii] == 0)    {
                        (*map)->istop[ii] = 1;
                        (*map)->dz3d[ii] = (bath[(*map)->top2d[ii]] - param->botZ) / incre_sum;
                        (*map)->bot3d[ii] = bath[(*map)->top2d[ii]] - (*map)->dz3d[ii];
                        for (jj = 1; jj < param->nz; jj++)  {
                            (*map)->dz3d[ii+jj] = (*map)->dz3d[ii+jj-1] * param->dz_incre;
                            (*map)->bot3d[ii+jj] = (*map)->bot3d[ii+jj-1] - (*map)->dz3d[ii+jj];
                        }
                    }
                }
            }

        }
        // regular cartesian mesh
        else    {
            for (ii = 0; ii < param->n3ci; ii++)    {
                (*map)->bot3d[ii] = (*map)->bot1d[(*map)->kk[ii]];
                if ((*map)->kk[ii] == 0)    {(*map)->dz3d[ii] = param->dz;}
                else    {(*map)->dz3d[ii] = (*map)->bot3d[ii-1] - (*map)->bot3d[ii];}
            }
            // adjust dz with respect to bathymetry
            for (ii = 0; ii < param->n3ci; ii++)
            {
                if ((*map)->bot3d[ii] >= bath[(*map)->top2d[ii]] - 1e-5)    {
                    (*map)->actv[ii] = 0;
                    if ((*map)->bot3d[ii] - bath[(*map)->top2d[ii]] < param->dz)
                    {(*map)->bot3d[ii] = bath[(*map)->top2d[ii]];}
                }
                else
                {
                    (*map)->actv[ii] = 1;
                    if (bath[(*map)->top2d[ii]] - (*map)->bot3d[ii] <= param->dz)   {
                        if (bath[(*map)->top2d[ii]] - (*map)->bot3d[ii] >= 0.25*param->dz)  {
                            (*map)->dz3d[ii] = bath[(*map)->top2d[ii]] - (*map)->bot3d[ii];
                        }
                        else    {
                            (*map)->actv[ii] = 0;
                            (*map)->dz3d[ii+1] += bath[(*map)->top2d[ii]] - (*map)->bot3d[ii];
                        }
                    }
                }
            }
            // get the top cell index
            for (ii = 0; ii < param->n3ci; ii++)
            {
                if (ii == 0)    {if (bath[(*map)->top2d[ii]] - (*map)->bot3d[ii] <= param->dz)   {(*map)->istop[ii] = 1;}}
                else    {
                    if ((*map)->actv[ii] == 1 & (*map)->actv[ii-1] == 0)    {(*map)->istop[ii] = 1;}
                    else if ((*map)->actv[ii] == 1 & (*map)->kk[ii] == 0)    {(*map)->istop[ii] = 1;}
                }
            }
        }

        // build terrain-following mesh
        for (ii = 0; ii < param->n3ct; ii++)
        {
            (*map)->sinx[ii] = 0.0; (*map)->siny[ii] = 0.0;
            (*map)->cosx[ii] = 1.0; (*map)->cosy[ii] = 1.0;
            if (ii < param->n2ct)   {
                (*map)->sintx[ii] = 0.0; (*map)->sinty[ii] = 0.0;
                (*map)->costx[ii] = 0.0; (*map)->costy[ii] = 0.0;
            }
            (*map)->V[ii] = param->dx * param->dy * (*map)->dz3d[ii];
            (*map)->Az[ii] = param->dx * param->dy;
            (*map)->Ax[ii] = param->dy * (*map)->dz3d[ii];
            (*map)->Ay[ii] = param->dx * (*map)->dz3d[ii];
        }

        // get z-coordinates of each subsurface cell
        (*map)->nactv = 0;
        for (ii = 0; ii < param->n3ci; ii++)
        {
            // if ((*map)->actv[ii] == 0)  {(*map)->zcntr[ii] = 999;}
            // else    {(*map)->zcntr[ii] = (*map)->bot3d[ii] + 0.5*(*map)->dz3d[ii] - offset[0];}
            (*map)->zcntr[ii] = (*map)->bot3d[ii] + 0.5*(*map)->dz3d[ii] - offset[0];
            if ((*map)->actv[ii] == 1)  {(*map)->nactv += 1;}
        }
        // calculate iP, iM, jP, jM, kP, kM maps
        (*map)->iPjckc = malloc(param->n3ci*sizeof(int));
        (*map)->iMjckc = malloc(param->n3ci*sizeof(int));
        (*map)->icjPkc = malloc(param->n3ci*sizeof(int));
        (*map)->icjMkc = malloc(param->n3ci*sizeof(int));
        (*map)->icjckP = malloc(param->n3ci*sizeof(int));
        (*map)->icjckM = malloc(param->n3ci*sizeof(int));
        for (ii = 0; ii < param->n3ci; ii++)
        {
            jj = (*map)->top2d[ii];
            // iP and iM maps
            if ((*map)->ii[ii] == param->nx-1)
            {(*map)->iPjckc[ii] = param->n3ci + 2*param->nx*param->nz + (*map)->jj[ii]*param->nz + (*map)->kk[ii]; }
            else
            {(*map)->iPjckc[ii] = (*map)->cntr[ii] + param->nz;}
            if ((*map)->ii[ii] == 0)
            {(*map)->iMjckc[ii] = param->n3ci + 2*param->nx*param->nz + param->ny*param->nz + (*map)->jj[ii]*param->nz + (*map)->kk[ii];}
            else
            {(*map)->iMjckc[ii] = (*map)->cntr[ii] - param->nz;}
            // jP and jM maps
            if ((*map)->jj[ii] == param->ny-1)
            {(*map)->icjPkc[ii] = param->n3ci + (*map)->ii[ii]*param->nz + (*map)->kk[ii];}
            else
            {(*map)->icjPkc[ii] = (*map)->cntr[ii] + param->nx*param->nz;}
            if ((*map)->jj[ii] == 0)
            {(*map)->icjMkc[ii] = param->n3ci + param->nx*param->nz + (*map)->ii[ii]*param->nz + (*map)->kk[ii];}
            else
            {(*map)->icjMkc[ii] = (*map)->cntr[ii] - param->nx*param->nz;}
            // kP and kM maps
            if ((*map)->kk[ii] == param->nz-1)
            {(*map)->icjckP[ii] = param->n2ct*param->nz + (*map)->top2d[ii];}
            else
            {(*map)->icjckP[ii] = (*map)->cntr[ii] + 1;}
            if ((*map)->kk[ii] == 0)
            {(*map)->icjckM[ii] = param->n2ct*(param->nz+1) + (*map)->top2d[ii];}
            else
            {(*map)->icjckM[ii] = (*map)->cntr[ii] - 1;}
        }
        // calculate boundary maps
        (*map)->iPin = malloc(param->ny*param->nz*sizeof(int));
        (*map)->iPou = malloc(param->ny*param->nz*sizeof(int));
        (*map)->iMin = malloc(param->ny*param->nz*sizeof(int));
        (*map)->iMou = malloc(param->ny*param->nz*sizeof(int));
        (*map)->jPin = malloc(param->nx*param->nz*sizeof(int));
        (*map)->jPou = malloc(param->nx*param->nz*sizeof(int));
        (*map)->jMin = malloc(param->nx*param->nz*sizeof(int));
        (*map)->jMou = malloc(param->nx*param->nz*sizeof(int));
        (*map)->kPin = malloc(param->nx*param->ny*sizeof(int));
        (*map)->kPou = malloc(param->nx*param->ny*sizeof(int));
        (*map)->kMin = malloc(param->nx*param->ny*sizeof(int));
        (*map)->kMou = malloc(param->nx*param->ny*sizeof(int));
        for (ii = 0; ii < param->n3ci; ii++)
        {
            if ((*map)->ii[ii] == param->nx-1)
            {
                (*map)->iPin[(*map)->jj[ii]*param->nz+(*map)->kk[ii]] = ii;
                (*map)->iPou[(*map)->jj[ii]*param->nz+(*map)->kk[ii]] = (*map)->iPjckc[ii];
            }
            if ((*map)->ii[ii] == 0)
            {
                (*map)->iMin[(*map)->jj[ii]*param->nz+(*map)->kk[ii]] = ii;
                (*map)->iMou[(*map)->jj[ii]*param->nz+(*map)->kk[ii]] = (*map)->iMjckc[ii];
            }
            if ((*map)->jj[ii] == param->ny-1)
            {
                (*map)->jPin[(*map)->ii[ii]*param->nz+(*map)->kk[ii]] = ii;
                (*map)->jPou[(*map)->ii[ii]*param->nz+(*map)->kk[ii]] = (*map)->icjPkc[ii];
            }
            if ((*map)->jj[ii] == 0)
            {
                (*map)->jMin[(*map)->ii[ii]*param->nz+(*map)->kk[ii]] = ii;
                (*map)->jMou[(*map)->ii[ii]*param->nz+(*map)->kk[ii]] = (*map)->icjMkc[ii];
            }
            if ((*map)->kk[ii] == param->nz-1)
            {
                (*map)->kPin[(*map)->top2d[ii]] = ii;
                (*map)->kPou[(*map)->top2d[ii]] = (*map)->icjckP[ii];
            }
            if ((*map)->kk[ii] == 0)
            {
                (*map)->kMin[(*map)->top2d[ii]] = ii;
                (*map)->kMou[(*map)->top2d[ii]] = (*map)->icjckM[ii];
            }
        }
        // boundary actv
        for (ii = 0; ii < param->ny*param->nz; ii++)
        {
            // (*map)->actv[(*map)->iPou[ii]] = (*map)->actv[(*map)->iPin[ii]];
            // (*map)->actv[(*map)->iMou[ii]] = (*map)->actv[(*map)->iMin[ii]];
            (*map)->actv[(*map)->iPou[ii]] = 0;
            (*map)->actv[(*map)->iMou[ii]] = 0;
            (*map)->bot3d[(*map)->iPou[ii]] = (*map)->bot3d[(*map)->iPin[ii]];
            (*map)->bot3d[(*map)->iMou[ii]] = (*map)->bot3d[(*map)->iMin[ii]];
            (*map)->dz3d[(*map)->iPou[ii]] = (*map)->dz3d[(*map)->iPin[ii]];
            (*map)->dz3d[(*map)->iMou[ii]] = (*map)->dz3d[(*map)->iMin[ii]];
        }
        // yp and ym boundaries
        for (ii = 0; ii < param->nx*param->nz; ii++)
        {
            // (*map)->actv[(*map)->jPou[ii]] = (*map)->actv[(*map)->jPin[ii]];
            // (*map)->actv[(*map)->jMou[ii]] = (*map)->actv[(*map)->jMin[ii]];
            (*map)->actv[(*map)->jPou[ii]] = 0;
            (*map)->actv[(*map)->jMou[ii]] = 0;
            (*map)->bot3d[(*map)->jPou[ii]] = (*map)->bot3d[(*map)->jPin[ii]];
            (*map)->bot3d[(*map)->jMou[ii]] = (*map)->bot3d[(*map)->jMin[ii]];
            (*map)->dz3d[(*map)->jPou[ii]] = (*map)->dz3d[(*map)->jPin[ii]];
            (*map)->dz3d[(*map)->jMou[ii]] = (*map)->dz3d[(*map)->jMin[ii]];
        }
        // calculate the angles for terrain-following mesh

        if (param->follow_terrain == 1) {
            for (ii = 0; ii < param->n3ci; ii++)    {
                if ((*map)->ii[ii] < param->nx-1)    {
                    hdiff = fabs(((*map)->bot3d[(*map)->iPjckc[ii]] + 0.5*(*map)->dz3d[(*map)->iPjckc[ii]])
                        - ((*map)->bot3d[ii] + 0.5*(*map)->dz3d[ii]));
                    dist = sqrt(pow(hdiff,2.0) + pow(param->dx,2.0));
                    (*map)->sinx[ii] = hdiff / dist;
                    (*map)->cosx[ii] = param->dx / dist;
                }
                if ((*map)->jj[ii] < param->ny-1)    {
                    hdiff = fabs(((*map)->bot3d[(*map)->icjPkc[ii]] + 0.5*(*map)->dz3d[(*map)->icjPkc[ii]])
                        - ((*map)->bot3d[ii] + 0.5*(*map)->dz3d[ii]));
                    dist = sqrt(pow(hdiff,2.0) + pow(param->dy,2.0));
                    (*map)->siny[ii] = hdiff / dist;
                    (*map)->cosy[ii] = param->dy / dist;
                }
                // top angles
                if ((*map)->istop[ii] == 1)   {
                    jj = (*map)->top2d[ii];
                    // x
                    if (bath[jj] < bath[smap->iMjc[jj]])    {
                        hdiff = fabs(bath[jj] - bath[smap->iMjc[jj]] + 0.5*(*map)->dz3d[(*map)->iMjckc[ii]]);
                        dist = sqrt(pow(hdiff,2.0) + pow(param->dx,2.0));
                        (*map)->sintx[jj] = hdiff / dist;
                        (*map)->costx[jj] = param->dx / dist;
                    }
                    else if (bath[jj] < bath[smap->iPjc[jj]])   {
                        hdiff = fabs(bath[jj] - bath[smap->iPjc[jj]] + 0.5*(*map)->dz3d[(*map)->iPjckc[ii]]);
                        dist = sqrt(pow(hdiff,2.0) + pow(param->dx,2.0));
                        (*map)->sintx[jj] = hdiff / dist;
                        (*map)->costx[jj] = param->dx / dist;
                    }
                    // y
                    if (bath[jj] < bath[smap->icjM[jj]])    {
                        hdiff = fabs(bath[jj] - bath[smap->icjM[jj]] + 0.5*(*map)->dz3d[(*map)->icjMkc[ii]]);
                        dist = sqrt(pow(hdiff,2.0) + pow(param->dy,2.0));
                        (*map)->sinty[jj] = hdiff / dist;
                        (*map)->costy[jj] = param->dy / dist;
                    }
                    else if (bath[jj] < bath[smap->icjP[jj]])   {
                        hdiff = fabs(bath[jj] - bath[smap->icjP[jj]] + 0.5*(*map)->dz3d[(*map)->icjPkc[ii]]);
                        dist = sqrt(pow(hdiff,2.0) + pow(param->dy,2.0));
                        (*map)->sinty[jj] = hdiff / dist;
                        (*map)->costy[jj] = param->dy / dist;
                    }
                }
                // cell face areas and volumes
                (*map)->Ax[ii] = 0.5 * ((*map)->dz3d[ii] + (*map)->dz3d[(*map)->iPjckc[ii]]) * (*map)->cosx[ii] * param->dy;
                (*map)->Ay[ii] = 0.5 * ((*map)->dz3d[ii] + (*map)->dz3d[(*map)->icjPkc[ii]]) * (*map)->cosy[ii] * param->dx;
                if ((*map)->ii[ii] == 0)    {(*map)->Ax[(*map)->iMjckc[ii]] = (*map)->dz3d[ii] * param->dy;}
                if ((*map)->jj[ii] == 0)    {(*map)->Ay[(*map)->icjMkc[ii]] = (*map)->dz3d[ii] * param->dx;}
            }
        }

        // ZhiLi20201229!!!
        // for (ii = 0; ii < param->n3ct; ii++)
        // {(*map)->dz3d[ii] = param->dz;}
        // output the z-coordinates
        if (param->use_mpi == 1)
        {
            int root = 0;
            mpi_gather_double((*map)->zcntr_root, (*map)->zcntr, param->n3ci, root);
            if (irank == root)
            {
                reorder_subsurf((*map)->zcntr_out, (*map)->zcntr_root, *map, param);
                write_one_file((*map)->zcntr_out, "zcell", param, 0, param->N3CI);
            }
        }
        else
        {
            for (ii = 0; ii < param->n3ci; ii++)    {(*map)->zcntr_out[ii] = (*map)->zcntr[ii];}
            write_one_file((*map)->zcntr_out, "zcell", param, 0, param->N3CI);
        }
    }

    free(bath_min);
    free(bath_max);
    free(bath_max_global);
    free(bath_max_arr);
    free(bath_min_global);
    free(bath_min_arr);

}



// >>>>> Build connections for subsurface domain <<<<<
void build_subs_mat_map(double *actv, Map **map, Config *param)   {
    int ii, jj, kk, ixm, ixp, iym, iyp, izp, izm, irow = 0, nactv = 0;
    // remove inactive zones
    if (param->actv_file == 1)  {
        for (ii = 0; ii < param->n3ci; ii++)    {
            if (actv[(*map)->top2d[ii]] <= 0)   {(*map)->actv[ii] = 0;}
        }
        for (ii = 0; ii < param->n3ci; ii++)    {if ((*map)->actv[ii] == 1)   {nactv += 1;}}
        (*map)->nactv = nactv;
    }
    else if (param->follow_terrain == 1)    {(*map)->nactv = param->n3ci;}
    (*map)->dom2mat = malloc(param->n3ci*sizeof(int));
    (*map)->mat2dom = malloc((*map)->nactv*sizeof(int));
    (*map)->distxp = malloc((*map)->nactv*sizeof(int));
    (*map)->distxm = malloc((*map)->nactv*sizeof(int));
    (*map)->distyp = malloc((*map)->nactv*sizeof(int));
    (*map)->distym = malloc((*map)->nactv*sizeof(int));
    (*map)->distzp = malloc((*map)->nactv*sizeof(int));
    (*map)->distzm = malloc((*map)->nactv*sizeof(int));
    for (ii = 0; ii < param->n3ci; ii++)    {
        if ((*map)->actv[ii] == 1)    {
            (*map)->dom2mat[ii] = irow;
            (*map)->mat2dom[irow] = ii;
            irow += 1;
        }
        else {(*map)->dom2mat[ii] = -1;}
    }
    // compute distance in the matrix system
    for (ii = 0; ii < (*map)->nactv; ii++)   {
        jj = (*map)->mat2dom[ii];
        //xp
        if ((*map)->ii[jj] != param->nx-1 & (*map)->actv[(*map)->iPjckc[jj]] == 1)    {
            ixp = (*map)->dom2mat[(*map)->iPjckc[jj]];
            (*map)->distxp[ii] = ixp - ii;
        }
        else {(*map)->distxp[ii] = -1;}
        //xm
        if ((*map)->ii[jj] != 0 & (*map)->actv[(*map)->iMjckc[jj]] == 1)    {
            ixm = (*map)->dom2mat[(*map)->iMjckc[jj]];
            (*map)->distxm[ii] = ii - ixm;
        }
        else {(*map)->distxm[ii] = -1;}
        //yp
        if ((*map)->jj[jj] != param->ny-1 & (*map)->actv[(*map)->icjPkc[jj]] == 1)    {
            iyp = (*map)->dom2mat[(*map)->icjPkc[jj]];
            (*map)->distyp[ii] = iyp - ii;
        }
        else {(*map)->distyp[ii] = -1;}
        //ym
        if ((*map)->jj[jj] != 0 & (*map)->actv[(*map)->icjMkc[jj]] == 1)    {
            iym = (*map)->dom2mat[(*map)->icjMkc[jj]];
            (*map)->distym[ii] = ii - iym;
        }
        else {(*map)->distym[ii] = -1;}
        //zp
        if ((*map)->kk[jj] != param->nz-1 & (*map)->actv[(*map)->icjckP[jj]] == 1)    {
            izp = (*map)->dom2mat[(*map)->icjckP[jj]];
            (*map)->distzp[ii] = izp - ii;
        }
        else {(*map)->distzp[ii] = -1;}
        //zm
        if ((*map)->kk[jj] != 0 & (*map)->istop[jj] == 0 & (*map)->actv[(*map)->icjckM[jj]] == 1)    {
            izm = (*map)->dom2mat[(*map)->icjckM[jj]];
            (*map)->distzm[ii] = ii - izm;
        }
        else {(*map)->distzm[ii] = -1;}
    }
}
