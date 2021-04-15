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
void advective_flux(Data **data, Map *gmap, Config *param, int icell, int kk);
double dispersive_flux(Data **data, Map *gmap, Config *param, int icell, int kk, char* axis);
void enforce_scalar_bc(Data **data, Map *gmap, Config *param, int kk, int irank);
void update_rhovisc(Data **data, Map *gmap, Config *param, int irank);
void dispersion_tensor(Data **data, Map *gmap, Config *param);

// >>>>> Scalar transport for shallow water
void scalar_shallowwater(Data **data, Map *smap, Config *param, int irank, int nrank, int kk)
{
    int ii, jj, ll;
    double sip, sim, sjp, sjm, s_lim_hi, s_lim_lo, s_rainevap, Vre;
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
    }
    // update scalar
    for (ii = 0; ii < param->n2ci; ii++)
    {
        if ((*data)->Vflux[ii] > 0 & (*data)->dept[ii] > 0)
        {(*data)->s_surf[kk][ii] = (*data)->sm_surf[kk][ii] / (*data)->Vflux[ii];}
        else
        {(*data)->s_surf[kk][ii] = 0.0;}
        // apply the limiter
        if ((*data)->Vflux[ii] > 0 & (*data)->dept[ii] > 0)
        {
            if ((*data)->s_surf[kk][ii] > s_max[ii] & s_max[ii] != s_lim_hi & s_max[ii] != s_lim_lo)
            {(*data)->s_surf[kk][ii] = s_max[ii];}
            else if ((*data)->s_surf[kk][ii] < s_min[ii] & s_min[ii] != s_lim_hi)
            {(*data)->s_surf[kk][ii] = s_min[ii];}
            if ((*data)->s_surf[kk][ii] >= s_lim_hi | (*data)->s_surf[kk][ii] < s_lim_lo)
            {
                mpi_print("WARNING: Scalar extremes detected! Enforce the problematic cell to dry!", irank);
                printf("(ii,jj)=(%d,%d)  surf=%f, dept=%f, s=%f\n",smap->ii[ii],smap->jj[ii],(*data)->eta[ii],(*data)->dept[ii],(*data)->s_surf[kk][ii]);
                (*data)->s_surf[kk][ii] = 0.0;
                (*data)->eta[ii] -= (*data)->dept[ii];
                (*data)->dept[ii] = 0.0;
            }
        }
        // rainfall/evaporation
        if ((*data)->Vs[ii] > 0 & (*data)->dept[ii] > 0)
        {
            if ((*data)->rain_sum[0] > param->min_dept)
            {Vre = ((*data)->rain_sum[0] - (*data)->evap[0]) * (*data)->Asz[ii] * param->dt;}
            else
            {Vre = -(*data)->evap[0] * (*data)->Asz[ii] * param->dt;}
            (*data)->s_surf[kk][ii] = (*data)->s_surf[kk][ii] * (*data)->Vs[ii] / ((*data)->Vs[ii] + Vre);
        }
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
    double jip, jim, jjp, jjm, jkp, jkm, coeff;
    double *s_min, *s_max;
    s_lim_hi = 1000.0;
    s_lim_lo = 0.0;
    s_min = malloc(param->n3ci*sizeof(double));
    s_max = malloc(param->n3ci*sizeof(double));
    for (ii = 0; ii < param->n3ci; ii++)
    {s_min[ii] = s_lim_hi; s_max[ii] = s_lim_lo;}
    // update dispersion tensor
    dispersion_tensor(data, gmap, param);

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
            }
            else
            {
                (*data)->s_subs[kk][gmap->icjckM[ii]] = (*data)->s_subs[kk][ii];
                (*data)->sm_subs[kk][gmap->icjckM[ii]] = (*data)->sm_subs[kk][ii];
            }
        }
        // advective transport
        advective_flux(data, gmap, param, ii, kk);
        // diffusive transport
        jip = dispersive_flux(data, gmap, param, ii, kk, "x");
        jim = dispersive_flux(data, gmap, param, gmap->iMjckc[ii], kk, "x");
        jjp = dispersive_flux(data, gmap, param, ii, kk, "y");
        jjm = dispersive_flux(data, gmap, param, gmap->icjMkc[ii], kk, "y");
        jkp = dispersive_flux(data, gmap, param, ii, kk, "z");
        if (gmap->istop[ii] != 1)
        {jkm = dispersive_flux(data, gmap, param, gmap->icjckM[ii], kk, "z");}
        else
        {
            if (param->sim_shallowwater == 1 & (*data)->dept[gmap->top2d[ii]] <= 0)
            {jkm = 0.0;}
            else if (param->sim_shallowwater == 0)
            {jkm = 0.0;}
            else
            {jkm = dispersive_flux(data, gmap, param, gmap->icjckM[ii], kk, "z");}
        }
        (*data)->sm_subs[kk][ii] = (*data)->sm_subs[kk][ii] + param->dt * ((jip - jim) + (jjp - jjm) + (jkp - jkm));

        // surface-subsurface exchange
        if (param->sim_shallowwater == 1 & gmap->istop[ii] == 1)
        {
            jj = gmap->top2d[ii];
            if ((*data)->qseepage[jj] > 0)
            {
                (*data)->sseepage[kk][jj] = param->dx * param->dy * (*data)->wc[ii] * ((*data)->qseepage[jj] * (*data)->s_subs[kk][ii] + \
                    (2.0 * (*data)->Dzz[ii] * (*data)->wc[ii] / gmap->dz3d[ii]) * ((*data)->s_subs[kk][ii] - (*data)->s_subs[kk][gmap->icjckM[ii]]));
            }
            else
            {
                (*data)->sseepage[kk][jj] = param->dx * param->dy * (*data)->wc[ii] * ((*data)->qseepage[jj] * (*data)->s_subs[kk][gmap->icjckM[ii]] + \
                    (2.0 * (*data)->Dzz[ii] * (*data)->wc[ii] / gmap->dz3d[ii]) * ((*data)->s_subs[kk][ii] - (*data)->s_subs[kk][gmap->icjckM[ii]]));
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

    }
    for (ii = 0; ii < param->n3ci; ii++)
    {

        (*data)->Vgflux[ii] = (*data)->wcs[ii] * param->dx * param->dy * gmap->dz3d[ii];
        // update salinity
        if ((*data)->Vgflux[ii] > 0)
        {(*data)->s_subs[kk][ii] = (*data)->sm_subs[kk][ii] / (*data)->Vgflux[ii];}
        else
        {(*data)->s_subs[kk][ii] = 0.0;}

        // apply scalar limiter except at top layer
        if (gmap->istop[ii] == 1 & param->bctype_GW[5] != 0 & gmap->ii[ii] >= 10 & gmap->ii[ii] < 190)
        {
            // nothing happens
        }
        else if (gmap->istop[ii] == 1)
        {
            // remove limiter for the top layer
        }
        else if (gmap->jj[ii] == 0 & param->bctype_GW[3] == 2)
        {
            // remove limiter for side boundaries
        }
        else if (gmap->jj[ii] == param->ny-1 & param->bctype_GW[2] == 2)
        {
            // remove limiter for side boundaries
        }
        else
        {
            if ((*data)->s_subs[kk][ii] > s_max[ii] & s_max[ii] < s_lim_hi)
            {(*data)->s_subs[kk][ii] = s_max[ii];}
            else if ((*data)->s_subs[kk][ii] < s_min[ii] & s_min[ii] > 0)
            {(*data)->s_subs[kk][ii] = s_min[ii];}
        }

        // detect extreme values
        if ((*data)->s_subs[kk][ii] > s_lim_hi & gmap->actv[ii] == 1)
        {
            mpi_print("WARNING: Scalar extremes > 1000 detected for groundwater!", irank);
            (*data)->s_subs[kk][ii] = 0.0;
        }
        else if ((*data)->s_subs[kk][ii] < 0 & gmap->actv[ii] == 1)
        {
            mpi_print("WARNING: Scalar extremes < 0  detected for groundwater!", irank);
            // printf("scalar (%d,%d,%d) = %f   actv=%d wc=%f\n",gmap->ii[ii],gmap->jj[ii],gmap->kk[ii],(*data)->s_subs[kk][ii],gmap->actv[ii],(*data)->wc[ii]);
            (*data)->s_subs[kk][ii] = 0.0;
        }
        if (gmap->actv[ii] == 0)    {(*data)->s_subs[kk][ii] = 0.0;}

    }

    // enforce scalar boundary condition
    enforce_scalar_bc(data, gmap, param, kk, irank);

    free(s_min);
    free(s_max);
    if (param->use_mpi == 1)
    {
        mpi_exchange_subsurf((*data)->s_subs[kk], gmap, 2, param, irank, nrank);
        mpi_exchange_subsurf((*data)->sm_subs[kk], gmap, 2, param, irank, nrank);
    }

}

// >>>>> calculate advective flux across one face
void advective_flux(Data **data, Map *gmap, Config *param, int icell, int kk)
{
    double sip=0.0, sim=0.0, sjp=0.0, sjm=0.0, skp=0.0, skm=0.0;
    double fx=0.0, fy=0.0, fz=0.0;
    // x direction
    if ((*data)->qx[icell] < 0)    {sip = (*data)->s_subs[kk][icell];}
    else if ((*data)->qx[icell] > 0)   {sip = (*data)->s_subs[kk][gmap->iPjckc[icell]];}
    if ((*data)->qx[gmap->iMjckc[icell]] < 0)    {sim = (*data)->s_subs[kk][gmap->iMjckc[icell]];}
    else if ((*data)->qx[gmap->iMjckc[icell]] > 0)   {sim = (*data)->s_subs[kk][icell];}
    fx = ((*data)->qx[icell] * sip - (*data)->qx[gmap->iMjckc[icell]] * sim) * param->dy * gmap->dz3d[icell];

    // y direction
    if ((*data)->qy[icell] < 0)    {sjp = (*data)->s_subs[kk][icell];}
    else if ((*data)->qy[icell] > 0)   {sjp = (*data)->s_subs[kk][gmap->icjPkc[icell]];}
    else    {sjp = 0.0;}
    if ((*data)->qy[gmap->icjMkc[icell]] < 0)    {sjm = (*data)->s_subs[kk][gmap->icjMkc[icell]];}
    else if ((*data)->qy[gmap->icjMkc[icell]] > 0)   {sjm = (*data)->s_subs[kk][icell];}
    else    {sjm = 0.0;}
    fy = ((*data)->qy[icell] * sjp - (*data)->qy[gmap->icjMkc[icell]] * sjm) * param->dx * gmap->dz3d[icell];

    // z direction
    if ((*data)->qz[icell] > 0)    {skp = (*data)->s_subs[kk][gmap->icjckP[icell]];}
    else if ((*data)->qz[icell] < 0)   {skp = (*data)->s_subs[kk][icell];}
    else    {skp = 0.0;}
    if ((*data)->qz[gmap->icjckM[icell]] > 0)
    {
        skm = (*data)->s_subs[kk][icell];
        if (param->bctype_GW[5] == 2 & gmap->istop[icell] == 1)
        {
            // scalar doesn't leave subsurface domain if surface is dry
            if (param->sim_shallowwater == 1 & (*data)->dept[gmap->top2d[icell]] <= 0.0)   {skm = 0.0;}
            else if (param->sim_shallowwater == 0)  {skm = 0.0;}
        }
    }
    else if ((*data)->qz[gmap->icjckM[icell]] < 0)
    {
        skm = (*data)->s_subs[kk][gmap->icjckM[icell]];
        if (param->sim_shallowwater == 1 & (*data)->dept[gmap->top2d[icell]] <= 0.0 & gmap->istop[icell] == 1)   {skm = 0.0;}
        else if (param->sim_shallowwater == 0 & gmap->istop[icell] == 1)  {skm = 0.0;}
    }
    else    {skm = 0.0;}
    fz = ((*data)->qz[icell] * skp - (*data)->qz[gmap->icjckM[icell]] * skm) * param->dx * param->dy;

    (*data)->sm_subs[kk][icell] = (*data)->sm_subs[kk][icell] + param->dt * (fx + fy + fz);

}

// >>>>> calculate diffusive-dispersive flux across one face
double dispersive_flux(Data **data, Map *gmap, Config *param, int icell, int kk, char* axis)
{
    double coeff, sxp, sxm, syp, sym, szp, szm, fx, fy, fz;
    int iPjP, iPjM, iPkP, iPkM, iMjP, iMkP, jPkP, jPkM, jMkP;
    fx = 0.0;
    fy = 0.0;
    fz = 0.0;

    if (strcmp(axis, "x") == 0 & gmap->actv[icell] == 1)
    {
        coeff = gmap->dz3d[icell] * param->dy * param->dx;
        if ((*data)->Kx[icell] > 0 & gmap->actv[gmap->iPjckc[icell]] == 1)
        {
            fx = coeff * (*data)->Dxx[icell] * ((*data)->s_subs[kk][gmap->iPjckc[icell]] - (*data)->s_subs[kk][icell]) / param->dx / param->dx;
        }
        if ((*data)->Ky[icell] > 0 & gmap->actv[gmap->icjPkc[icell]] == 1)
        {
            iPjP = gmap->iPjckc[icell] + param->nz*param->nx;
            iPjM = gmap->iPjckc[icell] - param->nz*param->nx;
            syp = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->iPjckc[icell]] + (*data)->s_subs[kk][gmap->icjPkc[icell]] + (*data)->s_subs[kk][iPjP]);
            sym = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->iPjckc[icell]] + (*data)->s_subs[kk][gmap->icjMkc[icell]] + (*data)->s_subs[kk][iPjM]);
            fy = coeff * (*data)->Dxy[icell] * (syp - sym) / param->dy / param->dy;
        }
        if ((*data)->Kz[icell] > 0 & gmap->actv[gmap->icjckP[icell]] == 1)
        {
            iPkP = gmap->iPjckc[icell] + 1;
            iPkM = gmap->iPjckc[icell] - 1;
            szp = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->iPjckc[icell]] + (*data)->s_subs[kk][gmap->icjckP[icell]] + (*data)->s_subs[kk][iPkP]);
            szm = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->iPjckc[icell]] + (*data)->s_subs[kk][gmap->icjckM[icell]] + (*data)->s_subs[kk][iPkM]);
            fz = coeff * (*data)->Dxz[icell] * (szp - szm) / gmap->dz3d[icell] / gmap->dz3d[icell];
        }
    }
    else if (strcmp(axis, "y") == 0 & gmap->actv[icell] == 1)
    {
        coeff = gmap->dz3d[icell] * param->dy * param->dx;
        if ((*data)->Kx[icell] > 0 & gmap->actv[gmap->iPjckc[icell]] == 1)
        {
            iPjP = gmap->icjPkc[icell] + param->nz;
            iMjP = gmap->icjPkc[icell] - param->nz;
            sxp = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->iPjckc[icell]] + (*data)->s_subs[kk][gmap->icjPkc[icell]] + (*data)->s_subs[kk][iPjP]);
            sxm = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->iMjckc[icell]] + (*data)->s_subs[kk][gmap->icjPkc[icell]] + (*data)->s_subs[kk][iMjP]);
            fx = coeff * (*data)->Dyx[icell] * (sxp - sxm) / param->dx / param->dx;
        }
        if ((*data)->Ky[icell] > 0 & gmap->actv[gmap->icjPkc[icell]] == 1)
        {
            fy = coeff * (*data)->Dyy[icell] * ((*data)->s_subs[kk][gmap->icjPkc[icell]] - (*data)->s_subs[kk][icell]) / param->dy / param->dy;
            // if (gmap->jj[icell] == param->ny-1 & param->bctype_GW[2] == 1)  {fy = fy * 2.0;}
        }
        if ((*data)->Kz[icell] > 0 & gmap->actv[gmap->icjckP[icell]] == 1)
        {
            jPkP = gmap->icjPkc[icell] + 1;
            jPkM = gmap->icjPkc[icell] - 1;
            szp = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->icjPkc[icell]] + (*data)->s_subs[kk][gmap->icjckP[icell]] + (*data)->s_subs[kk][jPkP]);
            szm = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->icjPkc[icell]] + (*data)->s_subs[kk][gmap->icjckM[icell]] + (*data)->s_subs[kk][jPkM]);
            fz = coeff * (*data)->Dyz[icell] * (szp - szm) / gmap->dz3d[icell] / gmap->dz3d[icell];
        }
    }
    else if (strcmp(axis, "z") == 0 & gmap->actv[icell] == 1)
    {
        coeff = gmap->dz3d[icell] * param->dy * param->dx;
        if ((*data)->Kx[icell] > 0 & gmap->actv[gmap->iPjckc[icell]] == 1)
        {
            iPkP = gmap->icjckP[icell] + param->nz;
            iMkP = gmap->icjckP[icell] - param->nz;
            sxp = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->iPjckc[icell]] + (*data)->s_subs[kk][gmap->icjckP[icell]] + (*data)->s_subs[kk][iPkP]);
            sxm = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->iMjckc[icell]] + (*data)->s_subs[kk][gmap->icjckP[icell]] + (*data)->s_subs[kk][iMkP]);
            fx = coeff * (*data)->Dzx[icell] * (sxp - sxm) / param->dx / param->dx;
        }
        if ((*data)->Ky[icell] > 0 & gmap->actv[gmap->icjPkc[icell]] == 1)
        {
            jPkP = gmap->icjckP[icell] + param->nz*param->nx;
            jMkP = gmap->icjckP[icell] - param->nz*param->nx;
            syp = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->icjckP[icell]] + (*data)->s_subs[kk][gmap->icjPkc[icell]] + (*data)->s_subs[kk][jPkP]);
            sym = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->icjckP[icell]] + (*data)->s_subs[kk][gmap->icjMkc[icell]] + (*data)->s_subs[kk][jMkP]);
            fy = coeff * (*data)->Dzy[icell] * (syp - sym) / param->dy / param->dy;
        }
        if ((*data)->Kz[icell] > 0 & gmap->actv[gmap->icjckP[icell]] == 1)
        {
            fz = coeff * (*data)->Dzz[icell] * ((*data)->s_subs[kk][gmap->icjckP[icell]] - (*data)->s_subs[kk][icell]) / gmap->dz3d[icell] / gmap->dz3d[icell];
        }
    }

    return fx + fy + fz;
}


// >>>>> enforce scalar boundary condition
void enforce_scalar_bc(Data **data, Map *gmap, Config *param, int kk, int irank)
{
    int ii;
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
            if (param->bctype_GW[3] != 0 & irank < param->mpi_nx)
            {
                (*data)->s_subs[kk][gmap->icjMkc[ii]] = param->s_ym[kk];
            }
            else
            {
                (*data)->s_subs[kk][gmap->icjMkc[ii]] = (*data)->s_subs[kk][ii];
                (*data)->sm_subs[kk][gmap->icjMkc[ii]] = (*data)->sm_subs[kk][ii];
            }
        }
        if (gmap->jj[ii] == param->ny-1)
        {
            if (param->bctype_GW[2] != 0 & irank >= param->mpi_nx*(param->mpi_ny-1))
            {
                (*data)->s_subs[kk][gmap->icjPkc[ii]] = param->s_yp[kk];
            }
            else
            {
                (*data)->s_subs[kk][gmap->icjPkc[ii]] = (*data)->s_subs[kk][ii];
                (*data)->sm_subs[kk][gmap->icjPkc[ii]] = (*data)->sm_subs[kk][ii];
            }
        }
        if (gmap->istop[ii] == 1)
        {
            if (param->sim_shallowwater == 1)
            {
                (*data)->s_subs[kk][gmap->icjckM[ii]] = (*data)->s_surf[kk][gmap->top2d[ii]];
                (*data)->sm_subs[kk][gmap->icjckM[ii]] = (*data)->sm_surf[kk][gmap->top2d[ii]];
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
}


// >>>>> correct density and viscosity for scalars
void update_rhovisc(Data **data, Map *gmap, Config *param, int irank)
{
    int ii;
    if (param->n_scalar > 0)
    {
        for (ii = 0; ii < param->n3ct; ii++)
        {
            (*data)->r_rhon[ii] = (*data)->r_rho[ii];
            if (gmap->actv[ii] == 1)
            {
                (*data)->r_rho[ii] = 1.0 + (*data)->s_subs[0][ii] * 0.0007;
                (*data)->r_visc[ii] = 1.0 / (1.0 + (*data)->s_subs[0][ii] * 0.0022);
            }
            else
            {
                (*data)->r_rho[ii] = 1.0;
                (*data)->r_visc[ii] = 1.0;
            }

        }
    }
}

// >>>>> calculate diffusion-dispersion coefficients
void dispersion_tensor(Data **data, Map *gmap, Config *param)
{
    // for now, only considers homogeneous dispersion coefficients
    int ii;
    double q_abs;
    for (ii = 0; ii < param->n3ct; ii++)
    {
        (*data)->Dxx[ii] = param->difux;
        (*data)->Dxy[ii] = 0.0;
        (*data)->Dxz[ii] = 0.0;
        (*data)->Dyy[ii] = param->difuy;
        (*data)->Dyx[ii] = 0.0;
        (*data)->Dyz[ii] = 0.0;
        (*data)->Dzz[ii] = param->difuz;
        (*data)->Dzx[ii] = 0.0;
        (*data)->Dzy[ii] = 0.0;
    }
    for (ii = 0; ii < param->n3ct; ii++)
    {
        q_abs = pow(pow((*data)->qx[ii],2.0) + pow((*data)->qy[ii],2.0) + pow((*data)->qz[ii],2.0), 0.5);
        if (q_abs > 0.0)
        {
            (*data)->Dxx[ii] += param->disp_lon * pow((*data)->qx[ii],2.0) / q_abs + \
                param->disp_lat * pow((*data)->qy[ii],2.0) / q_abs + \
                param->disp_lat * pow((*data)->qz[ii],2.0) / q_abs;
            (*data)->Dxy[ii] += (param->disp_lon - param->disp_lat) * (*data)->qx[ii] * (*data)->qy[ii] / q_abs;
            (*data)->Dxz[ii] += (param->disp_lon - param->disp_lat) * (*data)->qx[ii] * (*data)->qz[ii] / q_abs;
            (*data)->Dyy[ii] += param->disp_lon * pow((*data)->qy[ii],2.0) / q_abs + \
                param->disp_lat * pow((*data)->qx[ii],2.0) / q_abs + \
                param->disp_lat * pow((*data)->qz[ii],2.0) / q_abs;
            (*data)->Dyx[ii] += (param->disp_lon - param->disp_lat) * (*data)->qx[ii] * (*data)->qy[ii] / q_abs;
            (*data)->Dyz[ii] += (param->disp_lon - param->disp_lat) * (*data)->qy[ii] * (*data)->qz[ii] / q_abs;
            (*data)->Dzz[ii] += param->disp_lon * pow((*data)->qz[ii],2.0) / q_abs + \
                param->disp_lat * pow((*data)->qy[ii],2.0) / q_abs + \
                param->disp_lat * pow((*data)->qx[ii],2.0) / q_abs;
            (*data)->Dzx[ii] += (param->disp_lon - param->disp_lat) * (*data)->qx[ii] * (*data)->qz[ii] / q_abs;
            (*data)->Dzy[ii] += (param->disp_lon - param->disp_lat) * (*data)->qy[ii] * (*data)->qz[ii] / q_abs;
        }
        if ((*data)->Dxy[ii] < 0.0) {(*data)->Dxy[ii] = 0.0;}
        if ((*data)->Dxz[ii] < 0.0) {(*data)->Dxz[ii] = 0.0;}
        if ((*data)->Dyx[ii] < 0.0) {(*data)->Dyx[ii] = 0.0;}
        if ((*data)->Dyz[ii] < 0.0) {(*data)->Dyz[ii] = 0.0;}
        if ((*data)->Dzx[ii] < 0.0) {(*data)->Dzx[ii] = 0.0;}
        if ((*data)->Dzy[ii] < 0.0) {(*data)->Dzy[ii] = 0.0;}
    }
}
