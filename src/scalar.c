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
#include "subroutines.h"
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
    int ii, jj, ll, aa;
    double sip, sim, sjp, sjm, s_lim_hi, s_lim_lo, s_rainevap, Vre;
    double *s_min, *s_max;
    s_lim_hi = 200.0;
    s_lim_lo = 0.0;
    s_min = malloc(param->n2ci*sizeof(double));
    s_max = malloc(param->n2ci*sizeof(double));
    for (ii = 0; ii < param->n2ci; ii++)
    // {s_min[ii] = s_lim_hi; s_max[ii] = s_lim_lo;}
    {s_min[ii] = (*data)->s_surf[kk][ii]; s_max[ii] = (*data)->s_surf[kk][ii];}
    // for (ii = 0; ii < param->n2ci; ii++)
    for (aa = 0; aa < param->nactv; aa++)
    {
        ii = smap->mat2dom[aa];
        // // scalar mass
        (*data)->sm_surf[kk][ii] = (*data)->s_surf[kk][ii] * (*data)->Vsn[ii];

        // scalar advection
        // x advection
        if ((*data)->Fu[ii] > 0)
        {
            sip = (*data)->s_surf[kk][ii];
            if (param->superbee == 1)   {
                sip = tvd_superbee((*data)->s_surf[kk][smap->iPjc[ii]], (*data)->s_surf[kk][ii],
                    (*data)->s_surf[kk][smap->iMjc[ii]], (*data)->uu[ii], param->dx, param->dt, param);
            }
        }
        else
        {
            sip = (*data)->s_surf[kk][smap->iPjc[ii]];
            if (param->superbee == 1 & smap->ii[ii] != param->nx-1) {
                sip = tvd_superbee((*data)->s_surf[kk][ii], (*data)->s_surf[kk][smap->iPjc[ii]],
                    (*data)->s_surf[kk][smap->iPjc[smap->iPjc[ii]]], (*data)->uu[ii], param->dx, param->dt, param);
            }
        }
        if ((*data)->Fu[smap->iMjc[ii]] > 0)
        {
            sim = (*data)->s_surf[kk][smap->iMjc[ii]];
            if (param->superbee == 1 & smap->ii[ii] != 0)   {
                sim = tvd_superbee((*data)->s_surf[kk][ii], (*data)->s_surf[kk][smap->iMjc[ii]],
                    (*data)->s_surf[kk][smap->iMjc[smap->iMjc[ii]]], (*data)->uu[smap->iMjc[ii]], param->dx, param->dt, param);
            }
        }
        else
        {
            sim = (*data)->s_surf[kk][ii];
            if (param->superbee == 1)   {
                sim = tvd_superbee((*data)->s_surf[kk][smap->iMjc[ii]], (*data)->s_surf[kk][ii],
                    (*data)->s_surf[kk][smap->iPjc[ii]], (*data)->uu[smap->iMjc[ii]], param->dx, param->dt, param);
            }
        }
        // y advetcion
        if ((*data)->Fv[ii] > 0)
        {
            sjp = (*data)->s_surf[kk][ii];
            if (param->superbee == 1)   {
                sjp = tvd_superbee((*data)->s_surf[kk][smap->icjP[ii]], (*data)->s_surf[kk][ii],
                    (*data)->s_surf[kk][smap->icjM[ii]], (*data)->vv[ii], param->dy, param->dt, param);
            }
        }
        else
        {
            sjp = (*data)->s_surf[kk][smap->icjP[ii]];
            if (param->superbee == 1 & smap->jj[ii] != param->ny-1) {
                sjp = tvd_superbee((*data)->s_surf[kk][ii], (*data)->s_surf[kk][smap->icjP[ii]],
                    (*data)->s_surf[kk][smap->icjP[smap->icjP[ii]]], (*data)->vv[ii], param->dy, param->dt, param);
            }
        }
        if ((*data)->Fv[smap->icjM[ii]] > 0)
        {
            sjm = (*data)->s_surf[kk][smap->icjM[ii]];
            if (param->superbee == 1 & smap->jj[ii] != 0)   {
                sjm = tvd_superbee((*data)->s_surf[kk][ii], (*data)->s_surf[kk][smap->icjM[ii]],
                    (*data)->s_surf[kk][smap->icjM[smap->icjM[ii]]], (*data)->vv[smap->icjM[ii]], param->dy, param->dt, param);
            }
        }
        else
        {
            sjm = (*data)->s_surf[kk][ii];
            if (param->superbee == 1)   {
                sjm = tvd_superbee((*data)->s_surf[kk][smap->icjM[ii]], (*data)->s_surf[kk][ii],
                    (*data)->s_surf[kk][smap->icjP[ii]], (*data)->vv[smap->icjM[ii]], param->dy, param->dt, param);
            }
        }
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
        {
            // (*data)->sm_surf[kk][ii] = (*data)->sm_surf[kk][ii] + (*data)->sseepage[kk][ii] * param->dt;

            int i, j;
            double Aini = param->dx * param->dy;;
            i = smap->ii[ii];
            j = smap->jj[ii];
            jj = param->nz*param->nx*i + param->nz*j;
            if ((*data)->Asz[ii] > 0.0) {Aini = (*data)->Asz[ii];}
            if ((*data)->qss[ii] > 0)
            {
                // scalar doesn't leave subsurface if surface is dry
                if ((*data)->dept[ii] > 0.0)    {
                    (*data)->sseepage[kk][ii] = (*data)->qss[ii] * (*data)->s_surfkP[kk][ii] + \
                        (2.0 * (*data)->Dzz[jj] / smap->dz[ii]) * ((*data)->s_surfkP[kk][ii] - (*data)->s_surf[kk][ii]);
                }
                else {(*data)->sseepage[kk][ii] = 0.0;}
            }
            else if ((*data)->qss[ii] == 0)  {(*data)->sseepage[kk][ii] = 0.0;}
            else
            {
                (*data)->sseepage[kk][ii] = (*data)->qss[ii] * (*data)->s_surf[kk][ii] + \
                    (2.0 * (*data)->Dzz[jj] / smap->dz[ii]) * ((*data)->s_surfkP[kk][ii] - (*data)->s_surf[kk][ii]);
            }

            (*data)->sm_surf[kk][ii] += (*data)->sseepage[kk][ii] * Aini * param->dt;
        }
        else {
            (*data)->sseepage[kk][ii] = 0.0;
        }

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
            if ((*data)->uu[ii] != 0.0 | (*data)->uu[smap->iMjc[ii]] != 0.0
                    | (*data)->vv[ii] != 0.0 | (*data)->vv[smap->icjM[ii]] != 0.0)
            {
                if ((*data)->s_surf[kk][ii] > s_max[ii] & s_max[ii] != s_lim_hi & s_max[ii] != s_lim_lo)
                {(*data)->s_surf[kk][ii] = s_max[ii];}
                else if ((*data)->s_surf[kk][ii] < s_min[ii] & s_min[ii] != s_lim_hi)
                {(*data)->s_surf[kk][ii] = s_min[ii];}
            }
            else if (param->sim_groundwater == 1 & param->Ksz > 0)  {
                if ((*data)->s_surf[kk][ii] > s_max[ii] & s_max[ii] != s_lim_hi & s_max[ii] != s_lim_lo)
                {(*data)->s_surf[kk][ii] = s_max[ii];}
                else if ((*data)->s_surf[kk][ii] < s_min[ii] & s_min[ii] != s_lim_hi)
                {(*data)->s_surf[kk][ii] = s_min[ii];}
            }

            if ((*data)->s_surf[kk][ii] >= s_lim_hi)
            {
                (*data)->s_surf[kk][ii] = s_lim_hi;
            }
            else if ((*data)->s_surf[kk][ii] < s_lim_lo)
            {
                (*data)->s_surf[kk][ii] = s_lim_lo;
            }
        }
        //
        // if (smap->ii[ii] == 10 & smap->jj[ii] == 25)   {
        //     int i, j;
        //     i = smap->ii[ii];
        //     j = smap->jj[ii];
        //     jj = param->nz*param->nx*i + param->nz*j;
        //     // printf(" SURF : dept=%f, qseep=%f, s=%f, ssubs=%f, smax=%f, smin=%f\n",
        //     //     (*data)->dept[ii],1e5*(*data)->qss[ii],(*data)->s_surf[kk][ii],(*data)->s_subs[kk][jj],s_max[ii],s_min[ii]);
        // }

        // rainfall/evaporation
        if ((*data)->Vflux[ii] > 0)
        {
            Vre = 0.0;
            // if ((*data)->rain_sum[0] > param->min_dept)
            // {Vre = ((*data)->rain_sum[0] - (*data)->evap[0]) * (*data)->Asz[ii] * param->dt;}
            if ((*data)->dept[ii] > param->min_dept & (*data)->rain[0] > 0.0)
            {Vre += (*data)->rain[0] * (*data)->Asz[ii] * param->dt;}
            else
            {Vre -= (*data)->evap[ii] * (*data)->Asz[ii] * param->dt;}
            (*data)->s_surf[kk][ii] = (*data)->s_surf[kk][ii] * (*data)->Vflux[ii] / ((*data)->Vflux[ii] + Vre);
        }


    }
    // Kuan 2019
    // for (ii = 0; ii < param->n2ci; ii++)
    // {
    //     if ((*data)->Vflux[ii] > 0 & (*data)->dept[ii] > 0) {(*data)->s_surf[kk][ii] = 35.0;}
    // }
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
    int ii, jj, aa;
    double sip, sim, sjp, sjm, skp, skm, s_lim_hi, s_lim_lo, s_rainevap;
    double jip, jim, jjp, jjm, jkp, jkm, coeff;
    double *s_min, *s_max;
    s_lim_hi = 200.0;
    s_lim_lo = 0.0;
    s_min = malloc(param->n3ci*sizeof(double));
    s_max = malloc(param->n3ci*sizeof(double));
    for (ii = 0; ii < param->n3ci; ii++)
    {s_min[ii] = s_lim_hi; s_max[ii] = s_lim_lo;}
    // update dispersion tensor
    dispersion_tensor(data, gmap, param);
    // for (ii = 0; ii < param->n3ci; ii++)
    for (aa = 0; aa < gmap->nactv; aa++)
    {
        ii = gmap->mat2dom[aa];
        // scalar mass
        (*data)->sm_subs[kk][ii] = (*data)->s_subs[kk][ii] * (*data)->Vgn[ii];
        // (*data)->sm_subs[kk][ii] = (*data)->s_subs[kk][ii];

        // (*data)->sm_subs[kk][ii] = (*data)->s_subs[kk][ii] * (*data)->Vgflux[ii];

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
        if (gmap->jj[ii] == param->ny-1 & irank >= param->mpi_nx*(param->mpi_ny-1)) {jjp = jjp * 2.0;}
        jjm = dispersive_flux(data, gmap, param, gmap->icjMkc[ii], kk, "y");
        if (gmap->jj[ii] == 0 & irank < param->mpi_nx) {jjm = jjm * 2.0;}
        jkp = dispersive_flux(data, gmap, param, ii, kk, "z");
        if (gmap->istop[ii] != 1)
        {jkm = dispersive_flux(data, gmap, param, gmap->icjckM[ii], kk, "z");}
        else
        {
            // jkm = 0.0;
            if (param->sim_shallowwater == 1 & (*data)->dept[gmap->top2d[ii]] <= 0)
            {jkm = 0.0;}
            else if (param->sim_shallowwater == 0)
            {jkm = 0.0;}
            else
            {jkm = dispersive_flux(data, gmap, param, gmap->icjckM[ii], kk, "z");}
        }

        (*data)->sm_subs[kk][ii] = (*data)->sm_subs[kk][ii] + param->dt * ((jip - jim) + (jjp - jjm) + (jkp - jkm));

        // (*data)->sm_subs[kk][ii] = (*data)->sm_subs[kk][ii] + param->dt * ((jip - jim)/(param->dx)
        //     + (jjp - jjm)/(param->dy) + (jkp - jkm)/gmap->dz3d[ii]);

        // surface-subsurface exchange
        if (param->sim_shallowwater == 1 & gmap->istop[ii] == 1)
        {
            // (*data)->sm_subs[kk][ii] -= param->dt * (*data)->sseepage[kk][gmap->top2d[ii]] / gmap->dz3d[ii];
            (*data)->sm_subs[kk][ii] -= param->dt * gmap->Az[ii] * (*data)->sseepage[kk][gmap->top2d[ii]];
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

        if ((*data)->Ky[ii] > 0)
        {
            if (gmap->actv[gmap->icjPkc[ii]] == 1)  {
                if ((*data)->s_subs[kk][gmap->icjPkc[ii]] > s_max[ii])    {s_max[ii] = (*data)->s_subs[kk][gmap->icjPkc[ii]];}
                if ((*data)->s_subs[kk][gmap->icjPkc[ii]] < s_min[ii])    {s_min[ii] = (*data)->s_subs[kk][gmap->icjPkc[ii]];}
            }
            else if (gmap->jj[ii] == param->ny-1 & param->bctype_GW[2] != 0)    {
                if ((*data)->s_subs[kk][gmap->icjPkc[ii]] > s_max[ii])    {s_max[ii] = (*data)->s_subs[kk][gmap->icjPkc[ii]];}
                if ((*data)->s_subs[kk][gmap->icjPkc[ii]] < s_min[ii])    {s_min[ii] = (*data)->s_subs[kk][gmap->icjPkc[ii]];}
            }
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
                if (param->sim_shallowwater == 1 & (*data)->dept[gmap->top2d[ii]] > 0)
                {
                    if ((*data)->s_subs[kk][gmap->icjckM[ii]] > s_max[ii])    {s_max[ii] = (*data)->s_subs[kk][gmap->icjckM[ii]];}
                    if ((*data)->s_subs[kk][gmap->icjckM[ii]] < s_min[ii])    {s_min[ii] = (*data)->s_subs[kk][gmap->icjckM[ii]];}
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

    for (aa = 0; aa < gmap->nactv; aa++)
    {
        ii = gmap->mat2dom[aa];
    // for (ii = 0; ii < param->n3ci; ii++)
    // {
        // update scalar concentration
        if ((*data)->Vgflux[ii] > 0)
        {(*data)->s_subs[kk][ii] = (*data)->sm_subs[kk][ii] / (*data)->Vgflux[ii];}
        // if ((*data)->Vg[ii] > 0)
        // {(*data)->s_subs[kk][ii] = (*data)->sm_subs[kk][ii] / (*data)->Vg[ii];}
        else
        {(*data)->s_subs[kk][ii] = 0.0;}

        // (*data)->s_subs[kk][ii] = (*data)->sm_subs[kk][ii];

        if (param->sim_shallowwater == 1)
        {
            if (gmap->istop[ii] == 1 & (*data)->qtop[gmap->top2d[ii]] > 0 & (*data)->dept[gmap->top2d[ii]] <= 0.0)
            {
                s_max[ii] += 0.01;
            }
        }
        // apply scalar limiter
        if ((*data)->s_subs[kk][ii] > s_max[ii] & s_max[ii] < s_lim_hi)
        {(*data)->s_subs[kk][ii] = s_max[ii];}
        else if ((*data)->s_subs[kk][ii] < s_min[ii] & s_min[ii] > 0)
        {(*data)->s_subs[kk][ii] = s_min[ii];}

        // detect extreme values
        if ((*data)->s_subs[kk][ii] > s_lim_hi & gmap->actv[ii] == 1)
        {
            // mpi_print("WARNING: Scalar extremes > 1000 detected for groundwater!", irank);
            // printf("(ii,jj,kk,irank)=(%d,%d,%d,%d)  head=%f, wc=%f, s=%f\n",
            //     gmap->ii[ii],gmap->jj[ii],gmap->kk[ii],irank,(*data)->h[ii],(*data)->wc[ii],(*data)->s_subs[kk][ii]);
            // (*data)->s_subs[kk][ii] = 0.0;
            (*data)->s_subs[kk][ii] = s_lim_hi;
        }
        else if ((*data)->s_subs[kk][ii] < 0 & gmap->actv[ii] == 1)
        {
            if ((*data)->s_subs[kk][ii] < -0.01) {
                mpi_print("WARNING: Scalar extremes < 0  detected for groundwater!", irank);
                printf("(ii,jj,kk,irank)=(%d,%d,%d,%d)  head=%f, wc=%f, s=%f\n",
                    gmap->ii[ii],gmap->jj[ii],gmap->kk[ii],irank,(*data)->h[ii],(*data)->wc[ii],(*data)->s_subs[kk][ii]);
            }
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
    double sip=0.0, sim=0.0, sjp=0.0, sjm=0.0, skp=0.0, skm=0.0, dz_eff;
    double fx=0.0, fy=0.0, fz=0.0, distp, distm;
    // x direction
    // xp
    if (param->superbee == 1)
    {
        if ((*data)->qx[icell] > 0)
        {
            if (gmap->ii[icell] == param->nx-1)
            {sip = (*data)->s_subs[kk][gmap->iPjckc[icell]];}
            else
            {
                sip = tvd_superbee((*data)->s_subs[kk][icell], (*data)->s_subs[kk][gmap->iPjckc[icell]],
                    (*data)->s_subs[kk][gmap->iPjckc[gmap->iPjckc[icell]]], (*data)->qx[icell]/gmap->Ax[icell], param->dx, param->dt, param);
            }
        }
        else if ((*data)->qx[icell] < 0)
        {
            sip = tvd_superbee((*data)->s_subs[kk][gmap->iPjckc[icell]], (*data)->s_subs[kk][icell],
                (*data)->s_subs[kk][gmap->iMjckc[icell]], (*data)->qx[icell]/gmap->Ax[icell], param->dx, param->dt, param);
        }
        else    {sip = 0.0;}
    }
    else
    {
        if ((*data)->qx[icell] < 0)    {sip = (*data)->s_subs[kk][icell];}
        else if ((*data)->qx[icell] > 0)   {sip = (*data)->s_subs[kk][gmap->iPjckc[icell]];}
    }
    // xm
    if (param->superbee == 1)
    {
        if ((*data)->qx[gmap->iMjckc[icell]] > 0)
        {
            sim = tvd_superbee((*data)->s_subs[kk][gmap->iMjckc[icell]], (*data)->s_subs[kk][icell],
                (*data)->s_subs[kk][gmap->iPjckc[icell]], (*data)->qx[gmap->iMjckc[icell]]/gmap->Ax[gmap->iMjckc[icell]], param->dx, param->dt, param);
        }
        else if ((*data)->qx[gmap->iMjckc[icell]] < 0)
        {
            if (gmap->ii[icell] == 0)
            {sim = (*data)->s_subs[kk][gmap->iMjckc[icell]];}
            else
            {
                sim = tvd_superbee((*data)->s_subs[kk][icell], (*data)->s_subs[kk][gmap->iMjckc[icell]],
                    (*data)->s_subs[kk][gmap->iMjckc[gmap->iMjckc[icell]]], (*data)->qx[gmap->iMjckc[icell]]/gmap->Ax[gmap->iMjckc[icell]], param->dx, param->dt, param);
            }
        }
        else    {sim = 0.0;}
    }
    else
    {
        if ((*data)->qx[gmap->iMjckc[icell]] < 0)    {sim = (*data)->s_subs[kk][gmap->iMjckc[icell]];}
        else if ((*data)->qx[gmap->iMjckc[icell]] > 0)   {sim = (*data)->s_subs[kk][icell];}
    }
    // consider lateral surface-subsurface exchange
    // distp = param->dx / gmap->cosx[icell];
    // distm = param->dx / gmap->cosx[gmap->iMjckc[icell]];
    distp = param->dx;   distm = param->dx;
    // fx = ((*data)->qx[icell] * sip - (*data)->qx[gmap->iMjckc[icell]] * sim) * param->dy * gmap->dz3d[icell];
    // fx = (*data)->qx[icell] * sip / distp - (*data)->qx[gmap->iMjckc[icell]] * sim / distm;
    fx = (*data)->qx[icell] * sip - (*data)->qx[gmap->iMjckc[icell]] * sim;

    // y direction
    // yp
    if (param->superbee == 1)
    {
        if ((*data)->qy[icell] > 0)
        {
            if (gmap->jj[icell] == param->ny-1)
            {sjp = (*data)->s_subs[kk][gmap->icjPkc[icell]];}
            else
            {
                sjp = tvd_superbee((*data)->s_subs[kk][icell], (*data)->s_subs[kk][gmap->icjPkc[icell]],
                    (*data)->s_subs[kk][gmap->icjPkc[gmap->icjPkc[icell]]], (*data)->qy[icell]/gmap->Ay[icell], param->dy, param->dt, param);
            }
        }
        else if ((*data)->qy[icell] < 0)
        {
            sjp = tvd_superbee((*data)->s_subs[kk][gmap->icjPkc[icell]], (*data)->s_subs[kk][icell],
                (*data)->s_subs[kk][gmap->icjMkc[icell]], (*data)->qy[icell]/gmap->Ay[icell], param->dy, param->dt, param);
        }
        else    {sjp = 0.0;}

    }
    else
    {
        if ((*data)->qy[icell] < 0)    {sjp = (*data)->s_subs[kk][icell];}
        else if ((*data)->qy[icell] > 0)   {sjp = (*data)->s_subs[kk][gmap->icjPkc[icell]];}
        else    {sjp = 0.0;}
    }
    // ym
    if (param->superbee == 1)
    {
        if ((*data)->qy[gmap->icjMkc[icell]] > 0)
        {
            sjm = tvd_superbee((*data)->s_subs[kk][gmap->icjMkc[icell]], (*data)->s_subs[kk][icell],
                (*data)->s_subs[kk][gmap->icjPkc[icell]], (*data)->qy[gmap->icjMkc[icell]]/gmap->Ay[gmap->icjMkc[icell]], param->dy, param->dt, param);
        }
        else if ((*data)->qy[gmap->icjMkc[icell]] < 0)
        {
            if (gmap->jj[icell] == 0)
            {sjm = (*data)->s_subs[kk][gmap->icjMkc[icell]];}
            else
            {
                sjm = tvd_superbee((*data)->s_subs[kk][icell], (*data)->s_subs[kk][gmap->icjMkc[icell]],
                    (*data)->s_subs[kk][gmap->icjMkc[gmap->icjMkc[icell]]], (*data)->qy[gmap->icjMkc[icell]]/gmap->Ay[gmap->icjMkc[icell]], param->dy, param->dt, param);
            }
        }
        else    {sjm = 0.0;}
    }
    else
    {
        if ((*data)->qy[gmap->icjMkc[icell]] < 0)    {sjm = (*data)->s_subs[kk][gmap->icjMkc[icell]];}
        else if ((*data)->qy[gmap->icjMkc[icell]] > 0)   {sjm = (*data)->s_subs[kk][icell];}
        else    {sjm = 0.0;}
    }
    // distp = param->dy / gmap->cosy[icell];
    // distm = param->dy / gmap->cosy[gmap->icjMkc[icell]];
    distp = param->dy;   distm = param->dy;
    // fy = ((*data)->qy[icell] * sjp - (*data)->qy[gmap->icjMkc[icell]] * sjm) * param->dx * gmap->dz3d[icell];
    // fy = (*data)->qy[icell] * sjp / distp - (*data)->qy[gmap->icjMkc[icell]] * sjm / distm;
    fy = (*data)->qy[icell] * sjp - (*data)->qy[gmap->icjMkc[icell]] * sjm;

    // z direction
    // zp
    if (param->superbee == 1)
    {
        if ((*data)->qz[icell] > 0)
        {
            if (gmap->kk[icell] == param->nz-1)
            {skp = (*data)->s_subs[kk][gmap->icjckP[icell]];}
            else
            {
                skp = tvd_superbee((*data)->s_subs[kk][icell], (*data)->s_subs[kk][gmap->icjckP[icell]],
                    (*data)->s_subs[kk][gmap->icjckP[gmap->icjckP[icell]]], (*data)->qz[icell]/gmap->Az[icell], gmap->dz3d[icell], param->dt, param);
            }
        }
        else if ((*data)->qz[icell] < 0)
        {
            skp = tvd_superbee((*data)->s_subs[kk][gmap->icjckP[icell]], (*data)->s_subs[kk][icell],
                (*data)->s_subs[kk][gmap->icjckM[icell]], (*data)->qz[icell]/gmap->Az[icell], gmap->dz3d[icell], param->dt, param);
        }
        else    {skp = 0.0;}
    }
    else
    {
        if ((*data)->qz[icell] > 0)    {skp = (*data)->s_subs[kk][gmap->icjckP[icell]];}
        else if ((*data)->qz[icell] < 0)   {skp = (*data)->s_subs[kk][icell];}
        else    {skp = 0.0;}
    }
    // zm
    if (param->superbee == 1)
    {
        if ((*data)->qz[gmap->icjckM[icell]] > 0)
        {
            skm = tvd_superbee((*data)->s_subs[kk][gmap->icjckM[icell]], (*data)->s_subs[kk][icell],
                (*data)->s_subs[kk][gmap->icjckP[icell]], (*data)->qz[gmap->icjckM[icell]]/gmap->Az[gmap->icjckM[icell]], gmap->dz3d[icell], param->dt, param);
        }
        else if ((*data)->qz[gmap->icjckM[icell]] < 0)
        {
            if (gmap->istop[icell] == 1)
            {skm = (*data)->s_subs[kk][gmap->icjckM[icell]];}
            else
            {
                skm = tvd_superbee((*data)->s_subs[kk][icell], (*data)->s_subs[kk][gmap->icjckM[icell]],
                    (*data)->s_subs[kk][gmap->icjckM[gmap->icjckM[icell]]], (*data)->qz[gmap->icjckM[icell]]/gmap->Az[gmap->icjckM[icell]], gmap->dz3d[icell], param->dt, param);
            }
        }
        else    {skm = 0.0;}
        // scalar doesn't leave subsurface domain if surface is dry
        if (param->bctype_GW[5] == 2 & gmap->istop[icell] == 1)
        {
            // scalar doesn't leave subsurface domain if surface is dry
            if (param->sim_shallowwater == 1 & (*data)->dept[gmap->top2d[icell]] <= 0.0)   {skm = 0.0;}
            else if (param->sim_shallowwater == 0)  {skm = 0.0;}
        }
    }
    else
    {
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
    }
    // moved s-s advection elsewhere, ZhiLi20220425
    if (gmap->kk[icell] == 0 | gmap->istop[icell] == 1) {
        if (param->sim_shallowwater == 1)   {skm = 0.0;}
    }
    // skm = 0.0;
    if (gmap->kk[icell] == param->nz-1) {distp = 0.5*gmap->dz3d[icell];}
    else {distp = 0.5*(gmap->dz3d[icell] + gmap->dz3d[gmap->icjckP[icell]]);}
    if (gmap->kk[icell] == 0 | gmap->istop[icell] == 1)   {distm = 0.5*gmap->dz3d[icell];}
    else {distm = 0.5*(gmap->dz3d[icell] + gmap->dz3d[gmap->icjckM[icell]]);}
    // fz = ((*data)->qz[icell] * skp - (*data)->qz[gmap->icjckM[icell]] * skm) * param->dx * param->dy;
    distm = gmap->dz3d[icell];  distp = gmap->dz3d[icell];
    // fz = (*data)->qz[icell] * skp / distp - (*data)->qz[gmap->icjckM[icell]] * skm / distm;
    fz = (*data)->qz[icell] * skp - (*data)->qz[gmap->icjckM[icell]] * skm;

    // (*data)->sm_subs[kk][icell] = (*data)->sm_subs[kk][icell] + param->dt * (fx + fy + fz);
    (*data)->sm_subs[kk][icell] = (*data)->sm_subs[kk][icell] + param->dt * (fx + fy + fz);

}

// >>>>> calculate diffusive-dispersive flux across one face
double dispersive_flux(Data **data, Map *gmap, Config *param, int icell, int kk, char* axis)
{
    double coeff, sxp, sxm, syp, sym, szp, szm, fx, fy, fz, dist;
    int iPjP, iPjM, iPkP, iPkM, iMjP, iMkP, jPkP, jPkM, jMkP;
    fx = 0.0;
    fy = 0.0;
    fz = 0.0;

    if (strcmp(axis, "x") == 0 & gmap->actv[icell] == 1)
    {
        if (gmap->actv[icell] == 1)
        {
            coeff = gmap->dz3d[icell] * param->dy * param->dx;
            if ((*data)->Kx[icell] > 0 & gmap->actv[gmap->iPjckc[icell]] == 1)
            {
                if (gmap->actv[gmap->iPjckc[icell]] == 1)
                {
                    dist = param->dx / gmap->cosx[icell];
                    fx = (*data)->Dxx[icell] * gmap->Ax[icell] * ((*data)->s_subs[kk][gmap->iPjckc[icell]] - (*data)->s_subs[kk][icell]) / dist;
                    // fx = coeff * (*data)->Dxx[icell] * ((*data)->s_subs[kk][gmap->iPjckc[icell]] - (*data)->s_subs[kk][icell]) / param->dx / param->dx;
                }
                // {fx = coeff * (*data)->Dxx[icell] * ((*data)->s_subs[kk][gmap->iPjckc[icell]] - (*data)->s_subs[kk][icell]) / param->dx / param->dx;}
            }
            if ((*data)->Ky[icell] > 0 & gmap->actv[gmap->icjPkc[icell]] == 1 & icell < param->n3ci)
            {
                iPjP = gmap->iPjckc[icell] + param->nz*param->nx;
                iPjM = gmap->iPjckc[icell] - param->nz*param->nx;
                syp = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->iPjckc[icell]] + (*data)->s_subs[kk][gmap->icjPkc[icell]] + (*data)->s_subs[kk][iPjP]);
                sym = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->iPjckc[icell]] + (*data)->s_subs[kk][gmap->icjMkc[icell]] + (*data)->s_subs[kk][iPjM]);
                // fy = coeff * (*data)->Dxy[icell] * (syp - sym) / param->dy / param->dy;
                dist = param->dy / gmap->cosy[icell];
                fy = (*data)->Dxy[icell] * gmap->Ay[icell] * (syp - sym) / dist;
            }
            if ((*data)->Kz[icell] > 0 & gmap->actv[gmap->icjckP[icell]] == 1 & icell < param->n3ci)
            {
                iPkP = gmap->iPjckc[icell] + 1;
                iPkM = gmap->iPjckc[icell] - 1;
                szp = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->iPjckc[icell]] + (*data)->s_subs[kk][gmap->icjckP[icell]] + (*data)->s_subs[kk][iPkP]);
                szm = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->iPjckc[icell]] + (*data)->s_subs[kk][gmap->icjckM[icell]] + (*data)->s_subs[kk][iPkM]);
                // fz = coeff * (*data)->Dxz[icell] * (szp - szm) / gmap->dz3d[icell] / gmap->dz3d[icell];
                dist = gmap->dz3d[icell];
                fz = (*data)->Dxz[icell] * gmap->Az[icell] * (szp - szm) / dist;
            }
        }
    }
    else if (strcmp(axis, "y") == 0)
    {
        if (gmap->actv[icell] == 1)
        {
            coeff = gmap->dz3d[icell] * param->dy * param->dx;
            if ((*data)->Kx[icell] > 0 & gmap->actv[gmap->iPjckc[icell]] == 1 & icell < param->n3ci)
            {
                iPjP = gmap->icjPkc[icell] + param->nz;
                iMjP = gmap->icjPkc[icell] - param->nz;
                sxp = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->iPjckc[icell]] + (*data)->s_subs[kk][gmap->icjPkc[icell]] + (*data)->s_subs[kk][iPjP]);
                sxm = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->iMjckc[icell]] + (*data)->s_subs[kk][gmap->icjPkc[icell]] + (*data)->s_subs[kk][iMjP]);
                // fx = coeff * (*data)->Dyx[icell] * (sxp - sxm) / param->dx / param->dx;
                dist = param->dx / gmap->cosx[icell];
                fx = (*data)->Dyx[icell] * gmap->Ay[icell] * (sxp - sxm) / dist;
            }
            if ((*data)->Ky[icell] > 0)
            {
                if (gmap->actv[gmap->icjPkc[icell]] == 1)
                {
                    // fy = coeff * (*data)->Dyy[icell] * ((*data)->s_subs[kk][gmap->icjPkc[icell]] - (*data)->s_subs[kk][icell]) / param->dy / param->dy;
                    dist = param->dy / gmap->cosy[icell];
                    fy = (*data)->Dyy[icell] * gmap->Ay[icell] * ((*data)->s_subs[kk][gmap->icjPkc[icell]] - (*data)->s_subs[kk][icell]) / dist;
                    if (gmap->jj[gmap->icjPkc[icell]] == 0 | gmap->jj[icell] == param->ny-1)  {fy = fy * 2.0;}
                }
            }
            if ((*data)->Kz[icell] > 0 & gmap->actv[gmap->icjckP[icell]] == 1 & icell < param->n3ci)
            {
                jPkP = gmap->icjPkc[icell] + 1;
                jPkM = gmap->icjPkc[icell] - 1;
                szp = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->icjPkc[icell]] + (*data)->s_subs[kk][gmap->icjckP[icell]] + (*data)->s_subs[kk][jPkP]);
                szm = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->icjPkc[icell]] + (*data)->s_subs[kk][gmap->icjckM[icell]] + (*data)->s_subs[kk][jPkM]);
                // fz = coeff * (*data)->Dyz[icell] * (szp - szm) / gmap->dz3d[icell] / gmap->dz3d[icell];
                dist = gmap->dz3d[icell];
                fz = (*data)->Dyz[icell] * gmap->Az[icell] * (szp - szm) / dist;
            }
        }

    }
    else if (strcmp(axis, "z") == 0)
    {
        if (gmap->actv[icell] == 1)
        {
            coeff = gmap->dz3d[icell] * param->dy * param->dx;
            if ((*data)->Kx[icell] > 0 & gmap->actv[gmap->iPjckc[icell]] == 1 & icell < param->n3ci)
            {
                iPkP = gmap->icjckP[icell] + param->nz;
                iMkP = gmap->icjckP[icell] - param->nz;
                sxp = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->iPjckc[icell]] + (*data)->s_subs[kk][gmap->icjckP[icell]] + (*data)->s_subs[kk][iPkP]);
                sxm = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->iMjckc[icell]] + (*data)->s_subs[kk][gmap->icjckP[icell]] + (*data)->s_subs[kk][iMkP]);
                // fx = coeff * (*data)->Dzx[icell] * (sxp - sxm) / param->dx / param->dx;
                dist = param->dx / gmap->cosx[icell];
                fx = (*data)->Dzx[icell] * gmap->Ax[icell] * (sxp - sxm) / dist;
            }
            if ((*data)->Ky[icell] > 0 & gmap->actv[gmap->icjPkc[icell]] == 1 & icell < param->n3ci)
            {
                jPkP = gmap->icjckP[icell] + param->nz*param->nx;
                jMkP = gmap->icjckP[icell] - param->nz*param->nx;
                syp = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->icjckP[icell]] + (*data)->s_subs[kk][gmap->icjPkc[icell]] + (*data)->s_subs[kk][jPkP]);
                sym = 0.25 * ((*data)->s_subs[kk][icell] + (*data)->s_subs[kk][gmap->icjckP[icell]] + (*data)->s_subs[kk][gmap->icjMkc[icell]] + (*data)->s_subs[kk][jMkP]);
                // fy = coeff * (*data)->Dzy[icell] * (syp - sym) / param->dy / param->dy;
                dist = param->dy / gmap->cosy[icell];
                fy = (*data)->Dzy[icell] * gmap->Ay[icell] * (syp - sym) / dist;
            }
            if ((*data)->Kz[icell] > 0 & gmap->actv[gmap->icjckP[icell]] == 1)
            {
                // fz = coeff * (*data)->Dzz[icell] * ((*data)->s_subs[kk][gmap->icjckP[icell]] - (*data)->s_subs[kk][icell]) / gmap->dz3d[icell] / gmap->dz3d[icell];
                fz = (*data)->Dzz[icell] * gmap->Az[icell] * ((*data)->s_subs[kk][gmap->icjckP[icell]] - (*data)->s_subs[kk][icell]) / gmap->dz3d[icell];
                if (gmap->kk[gmap->icjckP[icell]] == 0 | gmap->kk[icell] == param->nz-1)  {fz = fz * 2.0;}
            }
        }
        else
        {
            if (gmap->actv[gmap->icjckP[icell]] == 1)
            {
                if ((*data)->dept[gmap->top2d[gmap->icjckP[icell]]] > 0.0)  {
                    coeff = 2.0 * gmap->dz3d[icell] * param->dy * param->dx;
                    fz = (*data)->Dzz[icell] * gmap->Az[icell] * ((*data)->s_subs[kk][gmap->icjckP[icell]] - (*data)->s_subs[kk][icell]) / gmap->dz3d[gmap->icjckP[icell]];
                    // fz = coeff * (*data)->Dzz[icell] * ((*data)->s_subs[kk][gmap->icjckP[icell]] - (*data)->s_subs[kk][icell]) / gmap->dz3d[icell] / gmap->dz3d[icell];
                }
            }
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
                // Henry's problem, ZhiLi20210413
                // if (gmap->kk[ii] < 2)
                // {
                //     (*data)->s_subs[kk][gmap->icjPkc[ii]] = (*data)->s_subs[kk][ii];
                // }

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
            (*data)->r_rho[ii] = 1.0 + (*data)->s_subs[0][ii] * 0.000744;
            (*data)->r_visc[ii] = 1.0 / (1.0 + (*data)->s_subs[0][ii] * 0.0022);

            // (*data)->r_rho[ii] = 1.0;
            // (*data)->r_visc[ii] = 1.0;
            // if (gmap->actv[ii] == 1)
            // {
            //     (*data)->r_rho[ii] = 1.0 + (*data)->s_subs[0][ii] * 0.000744;
            //     // (*data)->r_visc[ii] = 1.0 - (*data)->s_subs[0][ii] * 0.00156;
            //     // (*data)->r_rho[ii] = 1.0 + (*data)->s_subs[0][ii] * 0.00068;
            //     (*data)->r_visc[ii] = 1.0 / (1.0 + (*data)->s_subs[0][ii] * 0.0022);
            //     // (*data)->r_visc[ii] = 1.0;
            // }
            // else if (ii > param->n3ci)
            // {
            //     (*data)->r_rho[ii] = 1.0 + (*data)->s_subs[0][ii] * 0.000744;
            //     (*data)->r_visc[ii] = 1.0 / (1.0 + (*data)->s_subs[0][ii] * 0.0022);
            // }
            // else
            // {
            //     (*data)->r_rho[ii] = 1.0;
            //     (*data)->r_visc[ii] = 1.0;
            // }

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
        (*data)->Dxx[ii] = param->difux * (*data)->wcs[ii];
        (*data)->Dxy[ii] = 0.0;
        (*data)->Dxz[ii] = 0.0;
        (*data)->Dyy[ii] = param->difuy * (*data)->wcs[ii];
        (*data)->Dyx[ii] = 0.0;
        (*data)->Dyz[ii] = 0.0;
        (*data)->Dzz[ii] = param->difuz * (*data)->wcs[ii];
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
