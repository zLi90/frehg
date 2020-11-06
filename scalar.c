// scalar transport
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

#include "configuration.h"
#include "initialize.h"
#include "map.h"
#include "mpifunctions.h"
#include "scalar.h"
#include "utility.h"

void scalar_shallowwater(Data **data, Map *smap, Config *param, int irank, int nrank, int kk);
void scalar_groundwater(Data **data, Map *gmap, Config *param, int irank, int nrank, int kk);

// >>>>> Scalar transport for shallow water
void scalar_shallowwater(Data **data, Map *smap, Config *param, int irank, int nrank, int kk)
{
    int ii, jj, ll;
    double sip, sim, sjp, sjm, s_lim_hi, s_lim_lo, s_rainevap;
    double *s_min, *s_max;
    s_lim_hi = 1000.0;
    s_lim_lo = 0.0;
    s_min = malloc(param->n2ci*sizeof(double));
    s_max = malloc(param->n2ci*sizeof(double));
    for (ii = 0; ii < param->n2ci; ii++)
    {s_min[ii] = s_lim_hi; s_max[ii] = s_lim_lo;}
    for (ii = 0; ii < param->n2ci; ii++)
    {
        // // scalar mass
        (*data)->sm_surf[kk][ii] = (*data)->s_surf[kk][ii] * (*data)->Vsn[ii];
        // scalar advection
        if ((*data)->Fu[ii] > 0)    {sip = (*data)->s_surf[kk][ii];}
        else    {sip = (*data)->s_surf[kk][smap->iPjc[ii]];}
        if ((*data)->Fu[smap->iMjc[ii]] > 0)    {sim = (*data)->s_surf[kk][smap->iMjc[ii]];}
        else    {sim = (*data)->s_surf[kk][ii];}
        if ((*data)->Fv[ii] > 0)    {sjp = (*data)->s_surf[kk][ii];}
        else    {sjp = (*data)->s_surf[kk][smap->icjP[ii]];}
        if ((*data)->Fv[smap->icjM[ii]] > 0)    {sjm = (*data)->s_surf[kk][smap->icjM[ii]];}
        else    {sjm = (*data)->s_surf[kk][ii];}
        (*data)->sm_surf[kk][ii] = (*data)->sm_surf[kk][ii] + param->dt *
            (-(*data)->Fu[ii]*sip + (*data)->Fu[smap->iMjc[ii]]*sim - (*data)->Fv[ii]*sjp + (*data)->Fv[smap->icjM[ii]]*sjm);
        // scalar diffusion
        (*data)->sm_surf[kk][ii] = (*data)->sm_surf[kk][ii] + param->dt *
            ((param->difux*(*data)->Asx[ii]/param->dx) * ((*data)->s_surf[kk][smap->iPjc[ii]]-(*data)->s_surf[kk][ii]) -
            (param->difux*(*data)->Asx[smap->iMjc[ii]]/param->dx) * ((*data)->s_surf[kk][ii]-(*data)->s_surf[kk][smap->iMjc[ii]]) +
            (param->difuy*(*data)->Asy[ii]/param->dy) * ((*data)->s_surf[kk][smap->icjP[ii]]-(*data)->s_surf[kk][ii]) -
            (param->difuy*(*data)->Asy[smap->icjM[ii]]/param->dy) * ((*data)->s_surf[kk][ii]-(*data)->s_surf[kk][smap->icjM[ii]]));
        // scalar from subsurface
        if (param->sim_groundwater == 1)
        {(*data)->sm_surf[kk][ii] = (*data)->sm_surf[kk][ii] + (*data)->sseepage[kk][ii] * param->dt;}
        // scalar limiter
        if ((*data)->Asx[ii] > 0 & (*data)->etan[smap->iPjc[ii]] > (*data)->bottom[smap->iPjc[ii]])
        {
            if ((*data)->s_surf[kk][smap->iPjc[ii]] > s_max[ii])    {s_max[ii] = (*data)->s_surf[kk][smap->iPjc[ii]];}
            if ((*data)->s_surf[kk][smap->iPjc[ii]] < s_min[ii])    {s_min[ii] = (*data)->s_surf[kk][smap->iPjc[ii]];}
        }
        if ((*data)->Asx[smap->iMjc[ii]] > 0 & (*data)->etan[smap->iMjc[ii]] > (*data)->bottom[smap->iMjc[ii]])
        {
            if ((*data)->s_surf[kk][smap->iMjc[ii]] > s_max[ii])    {s_max[ii] = (*data)->s_surf[kk][smap->iMjc[ii]];}
            if ((*data)->s_surf[kk][smap->iMjc[ii]] < s_min[ii])    {s_min[ii] = (*data)->s_surf[kk][smap->iMjc[ii]];}
        }
        if ((*data)->Asy[ii] > 0 & (*data)->etan[smap->icjP[ii]] > (*data)->bottom[smap->icjP[ii]])
        {
            if ((*data)->s_surf[kk][smap->icjP[ii]] > s_max[ii])    {s_max[ii] = (*data)->s_surf[kk][smap->icjP[ii]];}
            if ((*data)->s_surf[kk][smap->icjP[ii]] < s_min[ii])    {s_min[ii] = (*data)->s_surf[kk][smap->icjP[ii]];}
        }
        if ((*data)->Asy[smap->icjM[ii]] > 0 & (*data)->etan[smap->icjM[ii]] > (*data)->bottom[smap->icjM[ii]])
        {
            if ((*data)->s_surf[kk][smap->icjM[ii]] > s_max[ii])    {s_max[ii] = (*data)->s_surf[kk][smap->icjM[ii]];}
            if ((*data)->s_surf[kk][smap->icjM[ii]] < s_min[ii])    {s_min[ii] = (*data)->s_surf[kk][smap->icjM[ii]];}
        }
        if (param->sim_groundwater == 1 & param->Ksz > 0)
        {
            if ((*data)->s_surfkP[kk][ii] > s_max[ii])  {s_max[ii] = (*data)->s_surfkP[kk][ii];}
            if ((*data)->s_surfkP[kk][ii] < s_min[ii])  {s_min[ii] = (*data)->s_surfkP[kk][ii];}
        }
        // scalar from inflow
        if (param->n_inflow > 0)
        {
            for (ll = 0; ll < param->n_inflow; ll++)
            {
                for (jj = 0; jj < (*data)->inflowloc_len[ll]; jj++)
                {
                    if (ii == (*data)->inflowloc[ll][jj] & (*data)->current_inflow[ll] > 0)
                    {
                        (*data)->sm_surf[kk][ii] += (*data)->current_inflow[ll]*param->dt*(*data)->current_s_inflow[kk][ll]/(*data)->inflowloc_len[ll];
                        s_max[ii] = s_lim_hi;
                        s_min[ii] = -1;
                    }
                }
            }
        }
        // modify limiter when evap/rain exist
        if ((*data)->Vflux[ii] > 0 & (*data)->dept[ii] > 0)
        {
            s_rainevap = (*data)->sm_surf[kk][ii] / (*data)->Vflux[ii];
            // net evaporation
            if ((*data)->rain[0] - (*data)->evap[0] < 0)
            {
                if (s_max[ii] < s_rainevap)
                {s_max[ii] = s_max[ii] * (*data)->Vflux[ii] / ((*data)->Vflux[ii]+((*data)->rain[0]-(*data)->evap[0])*param->dt*(*data)->Asz[ii]);}
                else
                {
                    s_rainevap = s_rainevap * (*data)->Vflux[ii] / ((*data)->Vflux[ii]+((*data)->rain[0]-(*data)->evap[0])*param->dt*(*data)->Asz[ii]);
                    if (s_rainevap > s_max[ii]) {s_max[ii] = s_rainevap;}
                }
            }
            // net rainfall
            else if ((*data)->rain[0] - (*data)->evap[0] > 0)
            {
                if (s_min[ii] > s_rainevap)
                {s_min[ii] = s_min[ii] * (*data)->Vflux[ii] / ((*data)->Vflux[ii]+((*data)->rain[0]-(*data)->evap[0])*param->dt*(*data)->Asz[ii]);}
                else
                {
                    s_rainevap = s_rainevap * (*data)->Vflux[ii] / ((*data)->Vflux[ii]+((*data)->rain[0]-(*data)->evap[0])*param->dt*(*data)->Asz[ii]);
                    if (s_rainevap < s_min[ii]) {s_min[ii] = s_rainevap;}
                }
            }
        }

    }
    // update scalar
    for (ii = 0; ii < param->n2ci; ii++)
    {
        if ((*data)->Vflux[ii] > 0 & (*data)->dept[ii] > 0)
        {(*data)->s_surf[kk][ii] = (*data)->sm_surf[kk][ii] / (*data)->Vflux[ii];}
        else
        {(*data)->s_surf[kk][ii] = 0.0;}
        // if (smap->ii[ii]==30 & smap->jj[ii]==37)
        // {
        //     printf("depth, vflux, dm, scalar = %f, %f, %f, %f\n",(*data)->dept[ii]*1e4,(*data)->Vflux[ii],(*data)->sm_surf[kk][ii],(*data)->s_surf[kk][ii]);
        //     printf("---\n");
        // }
        // apply the limiter
        if ((*data)->Vflux[ii] > 0 & (*data)->dept[ii] > 0)
        {
            if ((*data)->s_surf[kk][ii] > s_max[ii] & s_max[ii] != s_lim_hi)
            {(*data)->s_surf[kk][ii] = s_max[ii];}
            else if ((*data)->s_surf[kk][ii] < s_min[ii])
            {(*data)->s_surf[kk][ii] = s_min[ii];}
            if ((*data)->s_surf[kk][ii] > s_lim_hi | (*data)->s_surf[kk][ii] < s_lim_lo)
            {
                mpi_print("WARNING: Scalar extremes detected! Enforce the problematic cell to dry!", irank);
                (*data)->s_surf[kk][ii] = 0.0;
                (*data)->eta[ii] -= (*data)->dept[ii];
                (*data)->dept[ii] = 0.0;
                // printf("(%d,%d)\n",smap->ii[ii],smap->jj[ii]);
            }
        }

        // if (smap->ii[ii]==45 & smap->jj[ii]==15)
        // {
        //     printf("s_surf, vflux, dept = %f, %f, %f\n",(*data)->s_surf[kk][ii],(*data)->Vflux[ii],(*data)->dept[ii]);
        // }


    }
    // enforce scalar boundary condition
    for (ii = 0; ii < param->nx; ii++)
    {
        (*data)->s_surf[kk][smap->jPou[ii]] = (*data)->s_surf[kk][smap->jPin[ii]];
        (*data)->s_surf[kk][smap->jMou[ii]] = (*data)->s_surf[kk][smap->jMin[ii]];
    }
    for (ii = 0; ii < param->ny; ii++)
    {
        (*data)->s_surf[kk][smap->iPou[ii]] = (*data)->s_surf[kk][smap->iPin[ii]];
        (*data)->s_surf[kk][smap->iMou[ii]] = (*data)->s_surf[kk][smap->iMin[ii]];
    }
    // scalar for tide
    for (jj = 0; jj < param->n_tide; jj++)
    {
        if ((*data)->tideloc[jj][0] != -1)
        {
            for (ii = 0; ii < (*data)->tideloc_len[jj]; ii++)
            {(*data)->s_surf[kk][(*data)->tideloc[jj][ii]] = (*data)->current_s_tide[kk][jj];}
        }
    }
    for (ii = 0; ii < param->n2ci; ii++)
    {
        if ((*data)->dept[ii] <= 0.0)   {(*data)->s_surf[kk][ii] = 0.0;}
        // // scalar mass
        // (*data)->sm_surf[kk][ii] = (*data)->s_surf[kk][ii] * (*data)->Vs[ii];
    }

    free(s_min);
    free(s_max);

    if (param->use_mpi == 1)
    {
        mpi_exchange_surf((*data)->s_surf[kk], smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->sm_surf[kk], smap, 2, param, irank, nrank);
    }


}



// >>>>> Scalar transport for groundwater
void scalar_groundwater(Data **data, Map *gmap, Config *param, int irank, int nrank, int kk)
{
    int ii, jj;
    double sip, sim, sjp, sjm, skp, skm, s_lim_hi, s_lim_lo, s_rainevap;
    double jip, jim, jjp, jjm, jkp, jkm;
    double *s_min, *s_max;
    s_lim_hi = 1000.0;
    s_lim_lo = 0.0;
    s_min = malloc(param->n3ci*sizeof(double));
    s_max = malloc(param->n3ci*sizeof(double));
    for (ii = 0; ii < param->n3ci; ii++)
    {s_min[ii] = s_lim_hi; s_max[ii] = s_lim_lo;}

    for (ii = 0; ii < param->n3ci; ii++)
    {
        // scalar mass
        (*data)->sm_subs[kk][ii] = (*data)->s_subs[kk][ii] * (*data)->Vgn[ii];
        // scalar of kM boundary
        if (gmap->istop[ii] == 1)
        {
            if (param->sim_shallowwater == 1)
            {
                (*data)->s_subs[kk][gmap->icjckM[ii]] = (*data)->s_surf[kk][gmap->top2d[ii]];
                (*data)->sm_subs[kk][gmap->icjckM[ii]] = (*data)->sm_subs[kk][gmap->top2d[ii]];
                (*data)->s_surfkP[kk][gmap->top2d[ii]] = (*data)->s_subs[kk][ii];
                // if (gmap->ii[ii]==45 & gmap->jj[ii]==31 & gmap->istop[ii]==1)
                // {
                //     printf("kk=%d, depth=%f, s_surf = %f\n",gmap->kk[ii], (*data)->dept[gmap->top2d[ii]], (*data)->s_subs[kk][gmap->icjckM[ii]]);
                //     printf("---\n");
                // }
            }
            else
            {
                (*data)->s_subs[kk][gmap->icjckM[ii]] = (*data)->s_subs[kk][ii];
                (*data)->sm_subs[kk][gmap->icjckM[ii]] = (*data)->sm_subs[kk][ii];
            }
        }
        //
        // scalar advection
        //
        // x direction
        if ((*data)->qx[ii] > 0)    {sip = (*data)->s_subs[kk][ii];}
        else if ((*data)->qx[ii] < 0)   {sip = (*data)->s_subs[kk][gmap->iPjckc[ii]];}
        else    {sip = 0.0;}
        if ((*data)->qx[gmap->iMjckc[ii]] > 0)    {sim = (*data)->s_subs[kk][gmap->iMjckc[ii]];}
        else if ((*data)->qx[gmap->iMjckc[ii]] < 0)   {sim = (*data)->s_subs[kk][ii];}
        else    {sim = 0.0;}
        // y direction
        if ((*data)->qy[ii] > 0)    {sjp = (*data)->s_subs[kk][ii];}
        else if ((*data)->qy[ii] < 0)   {sjp = (*data)->s_subs[kk][gmap->icjPkc[ii]];}
        else    {sjp = 0.0;}
        if ((*data)->qy[gmap->icjMkc[ii]] > 0)    {sjm = (*data)->s_subs[kk][gmap->icjMkc[ii]];}
        else if ((*data)->qy[gmap->icjMkc[ii]] < 0)   {sjm = (*data)->s_subs[kk][ii];}
        else    {sjm = 0.0;}
        // z direction
        if ((*data)->qz[ii] > 0)    {skp = (*data)->s_subs[kk][gmap->icjckP[ii]];}
        else if ((*data)->qz[ii] < 0)   {skp = (*data)->s_subs[kk][ii];}
        else    {skp = 0.0;}
        if ((*data)->qz[gmap->icjckM[ii]] > 0)    {skm = (*data)->s_subs[kk][ii];}
        else if ((*data)->qz[gmap->icjckM[ii]] < 0)   {skm = (*data)->s_subs[kk][gmap->icjckM[ii]];}
        else    {skm = 0.0;}
        // advective transport
        (*data)->sm_subs[kk][ii] = (*data)->sm_subs[kk][ii] + param->dt * (*data)->wcs[ii] *\
            ((-(*data)->qx[ii] * sip + (*data)->qx[gmap->iMjckc[ii]] * sim) * param->dy * gmap->dz3d[ii] + \
            (-(*data)->qy[ii] * sjp + (*data)->qy[gmap->icjMkc[ii]] * sjm) * param->dx * gmap->dz3d[ii] + \
            ((*data)->qz[ii] * skp - (*data)->qz[gmap->icjckM[ii]] * skm) * param->dx * param->dy);
        //
        // scalar diffusion
        //
        // x direction
        if ((*data)->Kx[ii] > 0 & gmap->actv[gmap->iPjckc[ii]] == 1)
        {jip = gmap->dz3d[ii] * param->dy * (*data)->wcs[ii] * ((*data)->s_subs[kk][gmap->iPjckc[ii]] - (*data)->s_subs[kk][ii]) / param->dx;}
        else    {jip = 0.0;}
        if ((*data)->Kx[gmap->iMjckc[ii]] > 0 & gmap->actv[gmap->iMjckc[ii]] == 1)
        {jim = gmap->dz3d[ii] * param->dy * (*data)->wcs[ii] * ((*data)->s_subs[kk][ii] - (*data)->s_subs[kk][gmap->iMjckc[ii]]) / param->dx;}
        else    {jim = 0.0;}
        // y direction
        if ((*data)->Ky[ii] > 0 & gmap->actv[gmap->icjPkc[ii]] == 1)
        {jjp = gmap->dz3d[ii] * param->dx * (*data)->wcs[ii] * ((*data)->s_subs[kk][gmap->icjPkc[ii]] - (*data)->s_subs[kk][ii]) / param->dy;}
        else    {jjp = 0.0;}
        if ((*data)->Ky[gmap->icjMkc[ii]] > 0 & gmap->actv[gmap->icjMkc[ii]] == 1)
        {jjm = gmap->dz3d[ii] * param->dx * (*data)->wcs[ii] * ((*data)->s_subs[kk][ii] - (*data)->s_subs[kk][gmap->icjMkc[ii]]) / param->dy;}
        else    {jjm = 0.0;}
        // z direction
        if ((*data)->Kz[ii] > 0 & gmap->dz3d[ii] > 0 & gmap->actv[gmap->icjckP[ii]] == 1)
        {jkp = param->dx * param->dy * (*data)->wcs[ii] * ((*data)->s_subs[kk][ii] - (*data)->s_subs[kk][gmap->icjckP[ii]]) / gmap->dz3d[ii];}
        else    {jkp = 0.0;}
        if ((*data)->Kz[gmap->icjckM[ii]] > 0 & gmap->dz3d[ii] > 0)
        {
            if (param->sim_shallowwater == 1 & gmap->istop[ii] == 1 & (*data)->dept[gmap->top2d[ii]] > 0)
            {jkm = param->dx * param->dy * (*data)->wcs[ii] * ((*data)->s_subs[kk][gmap->icjckM[ii]] - (*data)->s_subs[kk][ii]) / gmap->dz3d[ii];}
            else if (gmap->actv[gmap->icjckM[ii]] == 1)
            {jkm = param->dx * param->dy * (*data)->wcs[ii] * ((*data)->s_subs[kk][gmap->icjckM[ii]] - (*data)->s_subs[kk][ii]) / gmap->dz3d[ii];}
        }
        else    {jkm = 0.0;}
        // diffusive transport
        (*data)->sm_subs[kk][ii] = (*data)->sm_subs[kk][ii] + param->dt * \
            (param->difux * (jip - jim) + param->difuy * (jjp - jjm) + param->difuz * (jkm - jkp));

        // surface-subsurface exchange
        if (param->sim_shallowwater == 1 & gmap->istop[ii] == 1)
        {
            jj = gmap->top2d[ii];
            if ((*data)->qseepage[jj] > 0)
            {
                (*data)->sseepage[kk][jj] = param->dx * param->dy * param->wcs * ((*data)->qseepage[jj] * (*data)->s_subs[kk][ii] + \
                    (2.0 * param->difuz * (*data)->wcs[ii] / gmap->dz3d[ii]) * ((*data)->s_subs[kk][ii] - (*data)->s_subs[kk][gmap->icjckM[ii]]));
            }
            else
            {
                (*data)->sseepage[kk][jj] = param->dx * param->dy * param->wcs * ((*data)->qseepage[jj] * (*data)->s_subs[kk][gmap->icjckM[ii]] + \
                    (2.0 * param->difuz * (*data)->wcs[ii] / gmap->dz3d[ii]) * ((*data)->s_subs[kk][ii] - (*data)->s_subs[kk][gmap->icjckM[ii]]));
            }
        }
        //
        // scalar limiter
        //
        if ((*data)->Kx[ii] > 0 & gmap->actv[gmap->iPjckc[ii]] == 1)
        {
            if ((*data)->s_subs[kk][gmap->iPjckc[ii]] > s_max[ii])    {s_max[ii] = (*data)->s_subs[kk][gmap->iPjckc[ii]];}
            if ((*data)->s_subs[kk][gmap->iPjckc[ii]] < s_min[ii])    {s_min[ii] = (*data)->s_subs[kk][gmap->iPjckc[ii]];}
        }
        if ((*data)->Kx[gmap->iMjckc[ii]] > 0 & gmap->actv[gmap->iMjckc[ii]] == 1)
        {
            if ((*data)->s_subs[kk][gmap->iMjckc[ii]] > s_max[ii])    {s_max[ii] = (*data)->s_subs[kk][gmap->iMjckc[ii]];}
            if ((*data)->s_subs[kk][gmap->iMjckc[ii]] < s_min[ii])    {s_min[ii] = (*data)->s_subs[kk][gmap->iMjckc[ii]];}
        }
        if ((*data)->Ky[ii] > 0 & gmap->actv[gmap->icjPkc[ii]] == 1)
        {
            if ((*data)->s_subs[kk][gmap->icjPkc[ii]] > s_max[ii])    {s_max[ii] = (*data)->s_subs[kk][gmap->icjPkc[ii]];}
            if ((*data)->s_subs[kk][gmap->icjPkc[ii]] < s_min[ii])    {s_min[ii] = (*data)->s_subs[kk][gmap->icjPkc[ii]];}
        }
        if ((*data)->Ky[gmap->icjMkc[ii]] > 0 & gmap->actv[gmap->icjMkc[ii]] == 1)
        {
            if ((*data)->s_subs[kk][gmap->icjMkc[ii]] > s_max[ii])    {s_max[ii] = (*data)->s_subs[kk][gmap->icjMkc[ii]];}
            if ((*data)->s_subs[kk][gmap->icjMkc[ii]] < s_min[ii])    {s_min[ii] = (*data)->s_subs[kk][gmap->icjMkc[ii]];}
        }
        if ((*data)->Kz[ii] > 0 & gmap->actv[gmap->icjckP[ii]] == 1)
        {
            if ((*data)->s_subs[kk][gmap->icjckP[ii]] > s_max[ii])    {s_max[ii] = (*data)->s_subs[kk][gmap->icjckP[ii]];}
            if ((*data)->s_subs[kk][gmap->icjckP[ii]] < s_min[ii])    {s_min[ii] = (*data)->s_subs[kk][gmap->icjckP[ii]];}
        }
        if ((*data)->Kz[gmap->icjckM[ii]] > 0)
        {
            if (gmap->istop[ii] == 1)
            {
                if (param->sim_shallowwater == 1)
                {
                    if ((*data)->dept[gmap->top2d[ii]] > 0)
                    {
                        if ((*data)->s_subs[kk][gmap->icjckM[ii]] > s_max[ii])    {s_max[ii] = (*data)->s_subs[kk][gmap->icjckM[ii]];}
                        if ((*data)->s_subs[kk][gmap->icjckM[ii]] < s_min[ii])    {s_min[ii] = (*data)->s_subs[kk][gmap->icjckM[ii]];}
                    }
                    else
                    {
                        // nothing is needed
                        // if (s_max[ii] < 20 || s_min[ii] < 20)
                        // {
                        //     printf("s, smax, smin = %f, %f, %f\n",(*data)->sm_subs[kk][ii] / (*data)->Vg[ii],s_max[ii],s_min[ii]);
                        //     printf("s = %f,%f,%f,%f,%f\n",(*data)->s_subs[kk][gmap->iPjckc[ii]],(*data)->s_subs[kk][gmap->iMjckc[ii]], \
                        //         (*data)->s_subs[kk][gmap->icjPkc[ii]],(*data)->s_subs[kk][gmap->icjMkc[ii]],(*data)->s_subs[kk][gmap->icjckP[ii]]);
                        //     printf("----\n");
                        // }
                    }
                }
                else
                {
                    // NOTE: Need to implement limiter for Dirichlet BC when modeling groundwater alone, ZhiLi20201102
                }
            }
            else
            {
                if (gmap->actv[gmap->icjckM[ii]] == 1)
                {
                    if ((*data)->s_subs[kk][gmap->icjckM[ii]] > s_max[ii])    {s_max[ii] = (*data)->s_subs[kk][gmap->icjckM[ii]];}
                    if ((*data)->s_subs[kk][gmap->icjckM[ii]] < s_min[ii])    {s_min[ii] = (*data)->s_subs[kk][gmap->icjckM[ii]];}
                }
            }
        }
        // remove limiter when evap/rain exist
        if (gmap->istop[ii] == 1 & param->bctype_GW[5] == 2 & (*data)->dept[gmap->top2d[ii]] <= 0 & (*data)->Vg[ii] > 0)
        {
            s_rainevap = (*data)->sm_subs[kk][ii] / (*data)->Vgflux[ii];
            // net evaporation
            if ((*data)->qtop > 0)
            {
                if (s_max[ii] < s_rainevap)
                {s_max[ii] = s_max[ii] * (*data)->Vgflux[ii] /  \
                        ((*data)->Vgflux[ii] - (*data)->qtop*param->dt*param->dx*param->dy);}
                else
                {
                    s_rainevap = s_rainevap * (*data)->Vgflux[ii] /  \
                        ((*data)->Vgflux[ii] - (*data)->qtop*param->dt*param->dx*param->dy);
                    if (s_rainevap > s_max[ii]) {s_max[ii] = s_rainevap;}
                }
                // printf("smax = %f\n",s_max[ii]);
            }
        }
    }
    for (ii = 0; ii < param->n3ci; ii++)
    {

        // update scalar concentration
        if ((*data)->Vgflux[ii] > 0)
        {(*data)->s_subs[kk][ii] = (*data)->sm_subs[kk][ii] / (*data)->Vgflux[ii];}
        else
        {(*data)->s_subs[kk][ii] = 0.0;}
        // if (gmap->ii[ii]==45 & gmap->jj[ii]==15 & gmap->istop[ii]==1)
        // {
        //     printf("---\n");
        //     printf("depth, vgflux, scalar = %f, %f, %f\n",(*data)->dept[gmap->top2d[ii]]*1e4,(*data)->Vgflux[ii],(*data)->s_subs[kk][ii]);
        // }
        // apply scalar limiter
        if ((*data)->s_subs[kk][ii] > s_max[ii] & s_max[ii] < s_lim_hi)
        {(*data)->s_subs[kk][ii] = s_max[ii];}
        else if ((*data)->s_subs[kk][ii] < s_min[ii] & s_min[ii] > 0)
        {(*data)->s_subs[kk][ii] = s_min[ii];}
        // detect extreme values
        if ((*data)->s_subs[kk][ii] > s_lim_hi)
        {
            mpi_print("WARNING: Scalar extremes > 1000 detected for groundwater!", irank);
            (*data)->s_subs[kk][ii] = 0.0;
        }
        else if ((*data)->s_subs[kk][ii] < 0)
        {
            mpi_print("WARNING: Scalar extremes < 0  detected for groundwater!", irank);
            // printf("scalar (%d,%d,%d) = %f   actv=%d wc=%f\n",gmap->ii[ii],gmap->jj[ii],gmap->kk[ii],(*data)->s_subs[kk][ii],gmap->actv[ii],(*data)->wc[ii]);
            (*data)->s_subs[kk][ii] = 0.0;
        }
        // if (gmap->ii[ii]==45 & gmap->jj[ii]==15 & gmap->istop[ii]==1)
        // {
        //
        //     printf("depth, vgflux, scalar = %f, %f, %f\n",(*data)->dept[gmap->top2d[ii]]*1e4,(*data)->Vgflux[ii],(*data)->s_subs[kk][ii]);
        //     printf("---\n");
        // }
    }
    // enforce boundary conditions
    for (ii = 0; ii < param->n3ci; ii++)
    {
        if (gmap->ii[ii] == 0)
        {
            (*data)->s_subs[kk][gmap->iMjckc[ii]] = (*data)->s_subs[kk][ii];
            (*data)->sm_subs[kk][gmap->iMjckc[ii]] = (*data)->sm_subs[kk][ii];
        }
        if (gmap->ii[ii] == param->nx-1)
        {
            (*data)->s_subs[kk][gmap->iPjckc[ii]] = (*data)->s_subs[kk][ii];
            (*data)->sm_subs[kk][gmap->iPjckc[ii]] = (*data)->sm_subs[kk][ii];
        }
        if (gmap->jj[ii] == 0)
        {
            (*data)->s_subs[kk][gmap->icjMkc[ii]] = (*data)->s_subs[kk][ii];
            (*data)->sm_subs[kk][gmap->icjMkc[ii]] = (*data)->sm_subs[kk][ii];
        }
        if (gmap->jj[ii] == param->ny-1)
        {
            (*data)->s_subs[kk][gmap->icjPkc[ii]] = (*data)->s_subs[kk][ii];
            (*data)->sm_subs[kk][gmap->icjPkc[ii]] = (*data)->sm_subs[kk][ii];
        }
        if (gmap->istop[ii] == 1)
        {
            if (param->sim_shallowwater == 1)
            {
                (*data)->s_subs[kk][gmap->icjckM[ii]] = (*data)->s_surf[kk][gmap->top2d[ii]];
                (*data)->sm_subs[kk][gmap->icjckM[ii]] = (*data)->sm_subs[kk][gmap->top2d[ii]];
                (*data)->s_surfkP[kk][gmap->top2d[ii]] = (*data)->s_subs[kk][ii];
            }
            else
            {
                (*data)->s_subs[kk][gmap->icjckM[ii]] = (*data)->s_subs[kk][ii];
                (*data)->sm_subs[kk][gmap->icjckM[ii]] = (*data)->sm_subs[kk][ii];
            }
        }
        if (gmap->kk[ii] == param->nz-1)
        {
            (*data)->s_subs[kk][gmap->icjckP[ii]] = (*data)->s_subs[kk][ii];
            (*data)->sm_subs[kk][gmap->icjckP[ii]] = (*data)->sm_subs[kk][ii];
        }
    }
    free(s_min);
    free(s_max);
    if (param->use_mpi == 1)
    {
        mpi_exchange_subsurf((*data)->s_subs[kk], gmap, 2, param, irank, nrank);
        mpi_exchange_subsurf((*data)->sm_subs[kk], gmap, 2, param, irank, nrank);
    }
}
