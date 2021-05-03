// Solves the 2D depth-integrated Navier-Stokes equation (shallow water equation)
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>

#include "configuration.h"
#include "initialize.h"
#include "map.h"
#include "mpifunctions.h"
#include "scalar.h"
#include "utility.h"

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

void solve_shallowwater(Data **data, Map *smap, Map *gmap, Config *param, int irank, int nrank);
void shallowwater_velocity(Data **data, Map *smap, Map *gmap, Config *param, int irank, int nrank);
void momentum_source(Data **data, Map *smap, Config *param);
void wind_source(Data **data, Map *smap, Config *param, int ii);
void shallowwater_rhs(Data **data, Map *smap, Config *param);
void shallowwater_mat_coeff(Data **data, Map *smap, Config *param, int irank, int nrank);
void build_shallowwater_system(Data *data, Map *smap, Config *param, QMatrix A, Vector b);
void solve_shallowwater_system(Data **data, Map *smap, QMatrix A, Vector b, Vector x, Config *param);
void enforce_surf_bc(Data **data, Map *smap, Config *param, int irank, int nrank);
void cfl_limiter(Data **data, Map *smap, Config *param);
void evaprain(Data **data, Map *smap, Config *param);
void subsurface_source(Data **data, Map *smap, Config *param);
void waterfall_location(Data **data, Map *smap, Config *param);
void update_velocity(Data **data, Map *smap, Config *param, int irank);
void waterfall_velocity(Data **data, Map *smap, Config *param);
void enforce_velo_bc(Data **data, Map *smap, Config *param, int irank, int nrank);
void interp_velocity(Data **data, Map *smap, Config *param);
void update_drag_coef(Data **data, Config *param);
void volume_by_flux(Data **data, Map *smap, Config *param);
void update_subgrid_variable(Data **data, Map *smap, Config *param, int irank, int nrank);
void subgrid_index(Data **data, Map *smap, Config *param);
void subgrid_interp_and_combine(Data **data, Map *smap, Config *param, int irank, int nrank);


// >>>>> Top level shallowwater solver
void solve_shallowwater(Data **data, Map *smap, Map *gmap, Config *param, int irank, int nrank)
{
    int ii;
    // allocate linear system
    QMatrix A;
    Q_Constr(&A, "A", param->n2ci, False, Rowws, Normal, True);
    Vector b;
    V_Constr(&b, "b", param->n2ci, Normal, True);
    Vector x;
    V_Constr(&x, "x", param->n2ci, Normal, True);
    for (ii = 0; ii < param->n2ci; ii++)    {(*data)->etan[ii] = (*data)->eta[ii];}
    enforce_surf_bc(data, smap, param, irank, nrank);
    if (param->use_mpi == 1)
    {
        mpi_exchange_surf((*data)->etan, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->eta, smap, 2, param, irank, nrank);
    }
    momentum_source(data, smap, param);
    if (param->use_mpi == 1)
    {
        mpi_exchange_surf((*data)->Ex, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Ey, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Dx, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Dy, smap, 2, param, irank, nrank);
    }
    shallowwater_rhs(data, smap, param);
    shallowwater_mat_coeff(data, smap, param, irank, nrank);
    build_shallowwater_system(*data, smap, param, A, b);
    solve_shallowwater_system(data, smap, A, b, x, param);
    enforce_surf_bc(data, smap, param, irank, nrank);
    // Update depth
    cfl_limiter(data, smap, param);
    evaprain(data, smap, param);
    update_depth(data, smap, param, irank);
    // free memory
    Q_Destr(&A);
    V_Destr(&b);
    V_Destr(&x);
}

// >>>>> Velocity update for shallowwater solver
void shallowwater_velocity(Data **data, Map *smap, Map *gmap, Config *param, int irank, int nrank)
{
    if (param->sim_groundwater == 1)
    {subsurface_source(data, smap, param);}
    // printf("-----\n");
    if (param->use_mpi == 1)
    {mpi_exchange_surf((*data)->eta, smap, 2, param, irank, nrank);}
    update_depth(data, smap, param, irank);
    if (param->use_mpi == 1)
    {
        mpi_exchange_surf((*data)->dept, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->deptx, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->depty, smap, 2, param, irank, nrank);
    }
    // Update subgrid variables
    update_subgrid_variable(data, smap, param, irank, nrank);
    volume_by_flux(data, smap, param);
    if (param->use_mpi == 1)
    {
        mpi_exchange_surf((*data)->Vs, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Vsx, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Vsy, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Asx, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Asy, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Asz, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Aszx, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Aszy, smap, 2, param, irank, nrank);
    }
    // Update bottom drag
    update_drag_coef(data, param);
    // Update velocity
    waterfall_location(data, smap, param);
    update_velocity(data, smap, param, irank);
    // waterfall_velocity(data, smap, param);
    enforce_velo_bc(data, smap, param, irank, nrank);
    // printf("Velocity NEW : velo = %f, %f\n",(*data)->uu[30],(*data)->vv[30]);
    if (param->use_mpi == 1)
    {
        mpi_exchange_surf((*data)->uu, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->vv, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Fu, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Fv, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->CDx, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->CDy, smap, 2, param, irank, nrank);
    }
    // volume_by_flux(data, smap, param);
    interp_velocity(data, smap, param);
    if (param->use_mpi == 1)
    {
        mpi_exchange_surf((*data)->uy, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->vx, smap, 2, param, irank, nrank);
    }

}

// >>>>> Momentum source term
void momentum_source(Data **data, Map *smap, Config *param)
{
    int ii;
    double advX, advY, difX, difY, facdx, facdy, velx, vely, gradp;
    for (ii = 0; ii < param->n2ci; ii++)
    {
        // advection terms
        advX = (0.5/param->dx) * (((*data)->uu[ii]+fabs((*data)->uu[ii]))*((*data)->uu[ii]-(*data)->uu[smap->iMjc[ii]]) + \
                                ((*data)->uu[ii]-fabs((*data)->uu[ii]))*((*data)->uu[smap->iPjc[ii]]-(*data)->uu[ii])) + \
               (0.5/param->dy) * (((*data)->vx[ii]+fabs((*data)->vx[ii]))*((*data)->uu[ii]-(*data)->uu[smap->icjM[ii]]) + \
                                ((*data)->vx[ii]-fabs((*data)->vx[ii]))*((*data)->uu[smap->icjP[ii]]-(*data)->uu[ii]));
        advY = (0.5/param->dx) * (((*data)->uy[ii]+fabs((*data)->uy[ii]))*((*data)->vv[ii]-(*data)->vv[smap->iMjc[ii]]) + \
                                ((*data)->uy[ii]-fabs((*data)->uy[ii]))*((*data)->vv[smap->iPjc[ii]]-(*data)->vv[ii])) + \
               (0.5/param->dy) * (((*data)->vv[ii]+fabs((*data)->vv[ii]))*((*data)->vv[ii]-(*data)->vv[smap->icjM[ii]]) + \
                                ((*data)->vv[ii]-fabs((*data)->vv[ii]))*((*data)->vv[smap->icjP[ii]]-(*data)->vv[ii]));
        if ((*data)->uu[ii] == 0 | (*data)->cflx[ii] > 0.7)   {advX = 0.0;}
        else    {if ((*data)->cflx[ii] > 0.5)    {advX = advX * (0.7 - (*data)->cflx[ii]) / (0.7 - 0.5);}}
        if ((*data)->vv[ii] == 0 | (*data)->cfly[ii] > 0.7)   {advY = 0.0;}
        else    {if ((*data)->cfly[ii] > 0.5)    {advY = advY * (0.7 - (*data)->cfly[ii]) / (0.7 - 0.5);}}
        // diffusion terms
        difX = 0.0;
        difY = 0.0;
        if ((*data)->Vsx[ii] > 0.0)
        {
            difX = (param->viscx/(*data)->Vsx[ii]/param->dx) * \
                    ((*data)->Asx[ii]*((*data)->uu[smap->iPjc[ii]] - (*data)->uu[ii]) - \
                    (*data)->Asx[ii]*((*data)->uu[ii] - (*data)->uu[smap->iMjc[ii]])) + \
                    (param->viscy/(*data)->Vsx[ii]/param->dy) * \
                    ((*data)->Asy[ii]*((*data)->uu[smap->icjP[ii]] - (*data)->uu[ii]) - \
                    (*data)->Asy[smap->icjM[ii]]*((*data)->uu[ii] - (*data)->uu[smap->icjM[ii]]));
        }
        if ((*data)->Vsy[ii] > 0.0)
        {
            difY= (param->viscx/(*data)->Vsy[ii]/param->dx) * \
                    ((*data)->Asx[ii]*((*data)->vv[smap->iPjc[ii]] - (*data)->vv[ii]) - \
                    (*data)->Asx[smap->iMjc[ii]]*((*data)->vv[ii] - (*data)->vv[smap->iMjc[ii]])) + \
                    (param->viscy/(*data)->Vsy[ii]/param->dy) * \
                    ((*data)->Asy[ii]*((*data)->vv[smap->icjP[ii]] - (*data)->vv[ii]) - \
                    (*data)->Asy[ii]*((*data)->vv[ii] - (*data)->vv[smap->icjM[ii]]));
        }
        // drag terms
        facdx = 0.0;
        facdy = 0.0;
        velx = sqrt((*data)->uu[ii]*(*data)->uu[ii] + (*data)->vx[ii]*(*data)->vx[ii]);
        vely = sqrt((*data)->uy[ii]*(*data)->uy[ii] + (*data)->vv[ii]*(*data)->vv[ii]);
        if ((*data)->Vsx[ii] > 0.0) {facdx = (*data)->Aszx[ii]/(*data)->Vsx[ii];}
        if ((*data)->Vsy[ii] > 0.0) {facdy = (*data)->Aszy[ii]/(*data)->Vsy[ii];}

        if (param->difuwave == 0)
        {
            (*data)->Dx[ii] = 1.0 / (0.5 * param->dt * (*data)->CDx[ii] * velx * facdx + 1.0);
            (*data)->Dy[ii] = 1.0 / (0.5 * param->dt * (*data)->CDy[ii] * vely * facdy + 1.0);
            // momentum source
            (*data)->Ex[ii] = (*data)->uu[ii] + param->dt * (difX - advX);
            (*data)->Ey[ii] = (*data)->vv[ii] + param->dt * (difY - advY);
            // (*data)->Ex[ii] = (*data)->uu[ii];
            // (*data)->Ey[ii] = (*data)->vv[ii];
            if (param->sim_wind == 1)   {wind_source(data, smap, param, ii);}
            (*data)->Ex[ii] = (*data)->Ex[ii] * (*data)->Dx[ii];
            (*data)->Ey[ii] = (*data)->Ey[ii] * (*data)->Dy[ii];

            //
            // (*data)->Ex[ii] = (*data)->uu[ii];
            // (*data)->Ey[ii] = (*data)->vv[ii];
        }
        else
        {
            facdx = 0.0;
            facdy = 0.0;
            if ((*data)->Aszx[ii] > 0.0) {facdx = (*data)->Vsx[ii]/(*data)->Aszx[ii];}
            if ((*data)->Aszy[ii] > 0.0) {facdy = (*data)->Vsy[ii]/(*data)->Aszy[ii];}

            gradp = 0.0;
            if ((*data)->eta[smap->iPjc[ii]] > (*data)->eta[ii] & (*data)->dept[smap->iPjc[ii]] > param->min_dept)
            {gradp += pow(((*data)->eta[smap->iPjc[ii]] - (*data)->eta[ii]) / param->dx, 2.0);}
            else if ((*data)->eta[smap->iPjc[ii]] < (*data)->eta[ii] & (*data)->dept[ii] > param->min_dept)
            {gradp += pow(((*data)->eta[smap->iPjc[ii]] - (*data)->eta[ii]) / param->dx, 2.0);}
            if ((*data)->eta[smap->icjP[ii]] > (*data)->eta[ii] & (*data)->dept[smap->icjP[ii]] > param->min_dept)
            {gradp += pow(((*data)->eta[smap->icjP[ii]] - (*data)->eta[ii]) / param->dy, 2.0);}
            else if ((*data)->eta[smap->icjP[ii]] < (*data)->eta[ii] & (*data)->dept[ii] > param->min_dept)
            {gradp += pow(((*data)->eta[smap->icjP[ii]] - (*data)->eta[ii]) / param->dy, 2.0);}
            gradp = pow(gradp, 0.5);
            // avoid gradp being too small
            if (gradp < param->min_dept / param->dx)
            {gradp = param->min_dept / param->dx;}
            // gradp = pow(pow(((*data)->eta[smap->iPjc[ii]] - (*data)->eta[ii]) / param->dx, 2.0) + \
                pow(((*data)->eta[smap->icjP[ii]] - (*data)->eta[ii]) / param->dy, 2.0), 0.5);

            velx = pow(2.0 * param->grav * gradp * facdx / (*data)->CDx[ii], 0.5);
            vely = pow(2.0 * param->grav * gradp * facdy / (*data)->CDy[ii], 0.5);

            if (velx != 0.0 & facdx != 0.0 & (*data)->CDx[ii] != 0.0)
            {(*data)->Dx[ii] = facdx / (0.5 * (*data)->CDx[ii] * velx);}
            else
            {(*data)->Dx[ii] = 1.0;}
            if (vely != 0.0 & facdy != 0.0 & (*data)->CDy[ii] != 0.0)
            {(*data)->Dy[ii] = facdy / (0.5 * (*data)->CDy[ii] * vely);}
            else
            {(*data)->Dy[ii] = 1.0;}
            // limiter on D
            // if ((*data)->Dx[ii] > param->dt)    {(*data)->Dx[ii] = param->dt;}
            // if ((*data)->Dy[ii] > param->dt)    {(*data)->Dy[ii] = param->dt;}

            // momentum source
            (*data)->Ex[ii] = 0.0;
            (*data)->Ey[ii] = 0.0;
        }
    }
}

void wind_source(Data **data, Map *smap, Config *param, int ii)
{
    double phi, omega, tau, pi = 3.1415926;
    double tauXP, tauYP;
    // phi is the wind direction from the north
    phi = (*data)->current_winddir[0] + param->north_angle;
    // omega is the wind direction in rad from positive x
    omega = phi * pi / 180.0;
    // tau is the total wind stress
    tau = param->rhoa * param->Cw * \
    ((*data)->current_windspd[0] - (*data)->uu[ii]*cos(omega) - (*data)->vv[ii]*sin(omega)) * \
    ((*data)->current_windspd[0] - (*data)->uu[ii]*cos(omega) - (*data)->vv[ii]*sin(omega));
    // apply the thin layer model when necessary
    if ((*data)->deptx[ii] < param->hD)
    {tauXP = tau * exp(param->CwT*((*data)->deptx[ii]-param->hD)/param->hD);}
    else
    {tauXP = tau;}
    if ((*data)->depty[ii] < param->hD)
    {tauYP = tau * exp(param->CwT*((*data)->depty[ii]-param->hD)/param->hD);}
    else
    {tauYP = tau;}
    // add wind drag to the source term
    if ((*data)->deptx[ii] > 0)
    {(*data)->Ex[ii] += param->dt * tauXP * cos(omega) / ((*data)->deptx[ii] * param->rhow);}
    if ((*data)->depty[ii] > 0)
    {(*data)->Ey[ii] += param->dt * tauYP * sin(omega) / ((*data)->depty[ii] * param->rhow);}
}


// >>>>> Matrix right-hand-side
void shallowwater_rhs(Data **data, Map *smap, Config *param)
{
    int ii, jj, kk;
    for (ii = 0; ii < param->n2ci; ii++)
    {
        (*data)->Srhs[ii] = (*data)->eta[ii] * (*data)->Asz[ii] - param->dt * \
            ((*data)->Asx[ii]*(*data)->Ex[ii] - (*data)->Asx[smap->iMjc[ii]]*(*data)->Ex[smap->iMjc[ii]] + \
            (*data)->Asy[ii]*(*data)->Ey[ii] - (*data)->Asy[smap->icjM[ii]]*(*data)->Ey[smap->icjM[ii]]);
    }
    // inflow as a source term
    if (param->n_inflow > 0)
    {
        for (kk = 0; kk < param->n_inflow; kk++)
        {
            for (ii = 0; ii < (*data)->inflowloc_len[kk]; ii++)
            {
                jj = (*data)->inflowloc[kk][ii];
                (*data)->Srhs[jj] += (*data)->current_inflow[kk] * param->dt / (*data)->inflowloc_len[kk];
            }
        }
    }
}

// >>>>> Matrix coefficients
void shallowwater_mat_coeff(Data **data, Map *smap, Config *param, int irank, int nrank)
{
    int ii, im, jm, kk, jj;
    double coef;
    size_t n;
    if (param->difuwave == 0)
    {coef = param->grav * param->dt * param->dt;}
    else
    {coef = param->grav * param->dt;}

    for (ii = 0; ii < param->n2ci; ii++)
    {
        im = smap->iMjc[ii];
        jm = smap->icjM[ii];
        (*data)->Sxp[ii] = 0.0;
        (*data)->Sxm[ii] = 0.0;
        (*data)->Syp[ii] = 0.0;
        (*data)->Sym[ii] = 0.0;
        if ((*data)->Vsx[ii] > 0.0)
        {(*data)->Sxp[ii] = coef * (*data)->Asx[ii] * (*data)->Asx[ii] * (*data)->Dx[ii] / (*data)->Vsx[ii];}
        if ((*data)->Vsx[im] > 0.0)
        {(*data)->Sxm[ii] = coef * (*data)->Asx[im] * (*data)->Asx[im] * (*data)->Dx[im] / (*data)->Vsx[im];}
        if ((*data)->Vsy[ii] > 0.0)
        {(*data)->Syp[ii] = coef * (*data)->Asy[ii] * (*data)->Asy[ii] * (*data)->Dy[ii] / (*data)->Vsy[ii];}
        if ((*data)->Vsy[jm] > 0.0)
        {(*data)->Sym[ii] = coef * (*data)->Asy[jm] * (*data)->Asy[jm] * (*data)->Dy[jm] / (*data)->Vsy[jm];}
        (*data)->Sct[ii] = (*data)->Asz[ii] + (*data)->Sxp[ii] + (*data)->Sxm[ii] + (*data)->Syp[ii] + (*data)->Sym[ii];

        // avoid singularity of matrix
        if ((*data)->dept[ii]==0.0)
        {
            (*data)->Sct[ii] = param->dx * param->dy;
            (*data)->Srhs[ii] = (*data)->eta[ii] * param->dx * param->dy;
            if ((*data)->uu[ii]==0.0 & (*data)->uu[im]==0.0 & (*data)->vv[ii]==0.0 & (*data)->vv[jm]==0.0)
            {
                (*data)->Sxp[ii] = 0.0;
                (*data)->Sxm[ii] = 0.0;
                (*data)->Syp[ii] = 0.0;
                (*data)->Sym[ii] = 0.0;
            }

        }
        else if ((*data)->Sct[ii] == 0.0)
        {
            (*data)->Sct[ii] = param->dx * param->dy;
            (*data)->Srhs[ii] = (*data)->eta[ii] * param->dx * param->dy;
        }
        // x-boundary
        if (smap->ii[ii] == 0)
        {
            // outer boundary
            if (irank % param->mpi_nx == 0) {(*data)->Sct[ii] -= (*data)->Sxm[ii];}
            // inner boundary (Dirichlet type)
            else    {(*data)->Srhs[ii] += (*data)->Sxm[ii] * (*data)->eta[smap->iMjc[ii]];}
        }
        else if (smap->ii[ii] == param->nx-1)
        {
            if (irank % param->mpi_nx == param->mpi_nx - 1) {(*data)->Sct[ii] -= (*data)->Sxp[ii];}
            else {(*data)->Srhs[ii] += (*data)->Sxp[ii] * (*data)->eta[smap->iPjc[ii]];}
        }
        // y-boundary
        if (smap->jj[ii] == 0)
        {
            if (irank < param->mpi_nx)  {(*data)->Sct[ii] -= (*data)->Sym[ii];}
            else    {(*data)->Srhs[ii] += (*data)->Sym[ii] * (*data)->eta[smap->icjM[ii]];}
        }
        else if (smap->jj[ii] == param->ny-1)
        {
            if (irank >= param->mpi_nx*(param->mpi_ny-1))   {(*data)->Sct[ii] -= (*data)->Syp[ii];}
            else    {(*data)->Srhs[ii] += (*data)->Syp[ii] * (*data)->eta[smap->icjP[ii]];}
        }
    }
    // enforce tidal boundary condition
    for (kk = 0; kk < param->n_tide; kk++)
    {
        if ((*data)->tideloc[kk][0] != -1)
        {
            // n = sizeof((*data)->tideloc[kk]) / sizeof((*data)->tideloc[kk][0]);
            n = (*data)->tideloc_len[kk];
            for (ii = 0; ii < n; ii++)
            {
                jj = (*data)->tideloc[kk][ii];
                (*data)->Sct[jj] = 1.0;
                (*data)->Sxp[jj] = 0.0;     (*data)->Sxm[jj] = 0.0;
                (*data)->Syp[jj] = 0.0;     (*data)->Sym[jj] = 0.0;
                (*data)->Srhs[jj] = (*data)->current_tide[kk];
            }
        }
    }
}

// >>>>> Setup the linear system of equations
void build_shallowwater_system(Data *data, Map *smap, Config *param, QMatrix A, Vector b)
{
    size_t ii, jj, kk, n;
    int im, dist;
    for (ii = 1; ii <= param->n2ci; ii++)
    {
        im = ii - 1;
        kk = 0;
        n = 5;
        if (smap->ii[im] == 0)  {n -= 1;}
        if (smap->ii[im] == param->nx-1)    {n -= 1;}
        if (smap->jj[im] == 0)  {n -= 1;}
        if (smap->jj[im] == param->ny-1)    {n -= 1;}
        Q_SetLen(&A, ii, n);
        // ym
        dist = im - smap->icjM[im];
        if (ii-dist > 0 & ii-dist <= param->n2ci)
        {Q_SetEntry(&A, ii, kk, ii-dist, -data->Sym[im]);    kk++;}
        // xm
        dist = im - smap->iMjc[im];
        if (ii-dist > 0 & ii-dist <= param->n2ci)
        {Q_SetEntry(&A, ii, kk, ii-dist, -data->Sxm[im]);    kk++;}
        // ct
        Q_SetEntry(&A, ii, kk, ii, data->Sct[im]);
        kk++;
        // xp
        dist = smap->iPjc[im] - im;
        if (ii+dist > 0 & ii+dist <= param->n2ci)
        {Q_SetEntry(&A, ii, kk, ii+dist, -data->Sxp[im]);    kk++;}
        // yp
        dist = smap->icjP[im] - im;
        if (ii+dist > 0 & ii+dist <= param->n2ci)
        {Q_SetEntry(&A, ii, kk, ii+dist, -data->Syp[im]);    kk++;}
        // rhs
        V_SetCmp(&b, ii, data->Srhs[im]);
    }

}

// >>>>> Solve the shallow water system
void solve_shallowwater_system(Data **data, Map *smap, QMatrix A, Vector b, Vector x, Config *param)
{
    size_t ii;
    V_SetAllCmp(&x, 0.0);
    SetRTCAccuracy(0.00000001);
    CGIter(&A, &x, &b, 10000000, SSORPrecond, 1);
    for (ii = 0; ii < param->n2ci; ii++)    {(*data)->eta[ii] = V_GetCmp(&x, ii+1);}

}

// >>>>> Enforce boundary condition for free surface
void enforce_surf_bc(Data **data, Map *smap, Config *param, int irank, int nrank)
{
    int ii, jj, kk;
    int n;
    // remove negative surface elevation
    for (ii = 0; ii < param->n2ci; ii++)
    {
        if ((*data)->eta[ii] < (*data)->bottom[ii])
        {(*data)->eta[ii] = (*data)->bottom[ii];}
    }
    // enforce tidal elevation
    for (kk = 0; kk < param->n_tide; kk++)
    {
        if ((*data)->tideloc[kk][0] != -1)
        {
            n = (*data)->tideloc_len[kk];
            for (ii = 0; ii < n; ii++)
            {
                jj = (*data)->tideloc[kk][ii];
                (*data)->eta[jj] = (*data)->current_tide[kk];
            }
        }
    }
    // enforce surface elevation for ghost cells
    for (ii = 0; ii < param->n2ci; ii++)
    {
        // x boundary
        if (smap->ii[ii] == 0)
        {(*data)->eta[smap->iMjc[ii]] = (*data)->eta[ii];}
        else if (smap->ii[ii] == param->nx-1)
        {(*data)->eta[smap->iPjc[ii]] = (*data)->eta[ii];}
        // y boundary
        if (smap->jj[ii] == 0)
        {(*data)->eta[smap->icjM[ii]] = (*data)->eta[ii];}
        else if (smap->jj[ii] == param->ny-1)
        {(*data)->eta[smap->icjP[ii]] = (*data)->eta[ii];}
    }
}

// >>>>> cfl limiter for stability
void cfl_limiter(Data **data, Map *smap, Config *param)
{
    int ii, wet;
    double diff;
    // restrict wetting within 1 cell
    for (ii = 0; ii < param->n2ci; ii++)
    {
        diff = (*data)->eta[ii] - (*data)->bottom[ii];

        if ((*data)->dept[ii] <= 0 & diff > 0)
        {
            wet = 0;
            // check for adjacent cells
            if ((*data)->dept[smap->iPjc[ii]] > 0.0)    {wet = 1;}
            else if ((*data)->dept[smap->iMjc[ii]] > 0.0)    {wet = 1;}
            else if ((*data)->dept[smap->icjP[ii]] > 0.0)    {wet = 1;}
            else if ((*data)->dept[smap->icjM[ii]] > 0.0)    {wet = 1;}
            // force to dry for isolated wet cells
            if (wet == 0)
            {
                (*data)->eta[ii] = (*data)->bottom[ii];
                (*data)->cfl_active[ii] = 1;
            }
        }
    }
    // remove small depth
    for (ii = 0; ii < param->n2ci; ii++)
    {
        diff = (*data)->eta[ii] - (*data)->bottom[ii];
        if (diff > 0 & diff < param->min_dept)
        {(*data)->eta[ii] = (*data)->bottom[ii];}
    }
}

// >>>>> apply rainfall and evaporation
void evaprain(Data **data, Map *smap, Config *param)
{
    int ii;
    double diff;
    // rainfall
    (*data)->rain_sum[0] += (*data)->rain[0] * param->dt;
    for (ii = 0; ii < param->n2ci; ii++)
    {
        if (smap->ii[ii] != 0 & smap->ii[ii] != param->nx-1)
        {
            if (smap->jj[ii] != 0 & smap->jj[ii] != param->ny-1)
            {if ((*data)->rain_sum[0] > param->min_dept) {(*data)->eta[ii] += (*data)->rain_sum[0];}}
        }
    }
    // evaporation
    // only apply evaporation when rainfall = 0
    if ((*data)->rain[0] == 0.0)
    {
        for (ii = 0; ii < param->n2ci; ii++)
        {(*data)->eta[ii] -= (*data)->evap[ii] * param->dt;}
        // remove negative or small depth
        for (ii = 0; ii < param->n2ci; ii++)
        {
            diff = (*data)->eta[ii] - (*data)->bottom[ii];
            if (diff < 0)
            {(*data)->eta[ii] = (*data)->bottom[ii];}
            else if (diff < param->min_dept)
            {(*data)->eta[ii] = (*data)->bottom[ii];}
        }
    }
}

// >>>>> add the exchange flux to surface domain
void subsurface_source(Data **data, Map *smap, Config *param)
{
    int ii;
    double diff, qz;
    // NOTE: seepage is calculated in the subsurface framework
    //       to convert it into the surface framework, multiple by param->wcs
    // NOTE: The above note is wrong!!! Because qz is already the Darcy flux!
    // NOTE: Which one is wrong remains undecided!!!
    for (ii = 0; ii < param->n2ci; ii++)
    {
        diff = (*data)->eta[ii] - (*data)->bottom[ii];
        // infiltration
        if ((*data)->qseepage[ii] < 0)
        {
            // (*data)->eta[ii] += (*data)->qseepage[ii] * param->dt;
            (*data)->eta[ii] += (*data)->qseepage[ii] * param->dt * param->wcs;
            (*data)->reset_seepage[ii] = 1;
        }
        // seepage
        else
        {
            if ((*data)->qseepage[ii]*param->dt*param->wcs > param->min_dept)
            {
                (*data)->eta[ii] += (*data)->qseepage[ii] * param->dt * param->wcs;
                (*data)->reset_seepage[ii] = 1;
            }
            // if ((*data)->qseepage[ii]*param->dt > param->min_dept)
            // {
            //     (*data)->eta[ii] += (*data)->qseepage[ii] * param->dt;
            //     (*data)->reset_seepage[ii] = 1;
            // }
            else
            {(*data)->reset_seepage[ii] = 0;}
        }
    }
}

// >>>>> find the waterfall locations
void waterfall_location(Data **data, Map *smap, Config *param)
{
    int ii;
    double facd;
    for (ii = 0; ii < param->n2ci; ii++)
    {
        (*data)->wtfx[ii] = 0;
        (*data)->wtfy[ii] = 0;
    }
    // wtf on xp and yp
    for (ii = 0; ii < param->n2ci; ii++)
    {
        facd = (*data)->eta[smap->iPjc[ii]] - (*data)->bottomXP[ii];
        if ((*data)->eta[ii] < (*data)->bottomXP[ii] & facd > param->wtfh)
        {(*data)->wtfx[ii] = -1;}
        facd = (*data)->eta[smap->icjP[ii]] - (*data)->bottomYP[ii];
        if ((*data)->eta[ii] < (*data)->bottomYP[ii] & facd > param->wtfh)
        {(*data)->wtfy[ii] = -1;}
    }
    // wtf on xm and ym
    for (ii = 0; ii < param->n2ci; ii++)
    {
        facd = (*data)->eta[smap->iMjc[ii]] - (*data)->bottomXP[smap->iMjc[ii]];
        if ((*data)->eta[ii] < (*data)->bottomXP[smap->iMjc[ii]] & facd > param->wtfh)
        {(*data)->wtfx[ii] = 1;}
        facd = (*data)->eta[smap->icjM[ii]] - (*data)->bottomYP[smap->icjM[ii]];
        if ((*data)->eta[ii] < (*data)->bottomYP[smap->icjM[ii]] & facd > param->wtfh)
        {(*data)->wtfy[ii] = 1;}
    }
}

// >>>>> update face velocity
void update_velocity(Data **data, Map *smap, Config *param, int irank)
{
    int ii, count;
    double effhx, effhy, coef, gradp, velx, vely, velo, epsu, epsv, facdx, facdy;
    coef = param->grav * param->dt;
    // save velocity at previous time step
    for (ii = 0; ii < param->n2ct; ii++)
    {
        (*data)->un[ii] = (*data)->uu[ii];
        (*data)->vn[ii] = (*data)->vv[ii];
    }
    // update new velocity
    for (ii = 0; ii < param->n2ci; ii++)
    {
        (*data)->uu[ii] = 0.0;
        (*data)->vv[ii] = 0.0;
        if ((*data)->Vsx[ii] > 0)   {effhx = (*data)->Asx[ii] / (*data)->Vsx[ii];}
        else {effhx = 0.0;}
        if ((*data)->Vsy[ii] > 0)   {effhy = (*data)->Asy[ii] / (*data)->Vsy[ii];}
        else {effhy = 0.0;}
        if (param->difuwave == 0)
        {
            // ignore drag inversion for velocity update -- consistent with Frehd
            // (*data)->uu[ii] = ((*data)->Ex[ii] - coef * effhx * ((*data)->eta[smap->iPjc[ii]] - (*data)->eta[ii]));
            // (*data)->vv[ii] = ((*data)->Ey[ii] - coef * effhy * ((*data)->eta[smap->icjP[ii]] - (*data)->eta[ii]));

            (*data)->uu[ii] = ((*data)->Ex[ii] - coef * effhx * ((*data)->eta[smap->iPjc[ii]] - (*data)->eta[ii])) * (*data)->Dx[ii];
            (*data)->vv[ii] = ((*data)->Ey[ii] - coef * effhy * ((*data)->eta[smap->icjP[ii]] - (*data)->eta[ii])) * (*data)->Dy[ii];
            // (*data)->uu[ii] = ((*data)->Ex[ii] - coef * effhx * ((*data)->eta[smap->iPjc[ii]] - (*data)->eta[ii]));
            // (*data)->vv[ii] = ((*data)->Ey[ii] - coef * effhy * ((*data)->eta[smap->icjP[ii]] - (*data)->eta[ii]));
        }
        else
        {
            facdx = 0.0;
            facdy = 0.0;
            // if ((*data)->Vsx[ii] > 0.0) {facdx = (*data)->Aszx[ii]/(*data)->Vsx[ii];}
            // if ((*data)->Vsy[ii] > 0.0) {facdy = (*data)->Aszy[ii]/(*data)->Vsy[ii];}
            if ((*data)->Aszx[ii] > 0.0) {facdx = (*data)->Vsx[ii]/(*data)->Aszx[ii];}
            if ((*data)->Aszy[ii] > 0.0) {facdy = (*data)->Vsy[ii]/(*data)->Aszy[ii];}
            gradp = 0.0;
            if ((*data)->eta[smap->iPjc[ii]] > (*data)->eta[ii] & (*data)->dept[smap->iPjc[ii]] > 0.0)
            {gradp += pow(((*data)->eta[smap->iPjc[ii]] - (*data)->eta[ii]) / param->dx, 2.0);}
            else if ((*data)->eta[smap->iPjc[ii]] < (*data)->eta[ii] & (*data)->dept[ii] > 0.0)
            {gradp += pow(((*data)->eta[smap->iPjc[ii]] - (*data)->eta[ii]) / param->dx, 2.0);}
            if ((*data)->eta[smap->icjP[ii]] > (*data)->eta[ii] & (*data)->dept[smap->icjP[ii]] > 0.0)
            {gradp += pow(((*data)->eta[smap->icjP[ii]] - (*data)->eta[ii]) / param->dy, 2.0);}
            else if ((*data)->eta[smap->icjP[ii]] < (*data)->eta[ii] & (*data)->dept[ii] > 0.0)
            {gradp += pow(((*data)->eta[smap->icjP[ii]] - (*data)->eta[ii]) / param->dy, 2.0);}
            gradp = pow(gradp, 0.5);

            // avoid gradp being too small
            if (gradp < param->min_dept / param->dx)
            {gradp = param->min_dept / param->dx;}

            // velx = sqrt((*data)->uu[ii]*(*data)->uu[ii] + (*data)->vx[ii]*(*data)->vx[ii]);
            // vely = sqrt((*data)->uy[ii]*(*data)->uy[ii] + (*data)->vv[ii]*(*data)->vv[ii]);

            velx = pow(2.0 * param->grav * gradp * facdx / (*data)->CDx[ii], 0.5);
            vely = pow(2.0 * param->grav * gradp * facdy / (*data)->CDy[ii], 0.5);

            if (velx != 0.0 & facdx != 0.0 & (*data)->CDx[ii] != 0.0)
            {(*data)->Dx[ii] = facdx / (0.5 * (*data)->CDx[ii] * velx);}
            else
            {(*data)->Dx[ii] = 1.0;}
            if (vely != 0.0 & facdy != 0.0 & (*data)->CDy[ii] != 0.0)
            {(*data)->Dy[ii] = facdy / (0.5 * (*data)->CDy[ii] * vely);}
            else
            {(*data)->Dy[ii] = 1.0;}
            // limiter on D
            // if ((*data)->Dx[ii] > param->dt)    {(*data)->Dx[ii] = param->dt;}
            // if ((*data)->Dy[ii] > param->dt)    {(*data)->Dy[ii] = param->dt;}

            (*data)->uu[ii] = - param->grav * effhx * ((*data)->eta[smap->iPjc[ii]] - (*data)->eta[ii]) * (*data)->Dx[ii];
            (*data)->vv[ii] = - param->grav * effhy * ((*data)->eta[smap->icjP[ii]] - (*data)->eta[ii]) * (*data)->Dy[ii];
        }
    }
    // apply various velocity limiters
    for (ii = 0; ii < param->n2ci; ii++)
    {
        // zero velocity when face area is zero
        if ((*data)->Asx[ii] < param->wtfh*param->dy) {(*data)->uu[ii] = 0.0;}
        if ((*data)->Asy[ii] < param->wtfh*param->dx) {(*data)->vv[ii] = 0.0;}
        // zero velocity out of a dry cell
        if ((*data)->dept[ii] < param->wtfh)
        {
            if ((*data)->uu[ii] > 0)    {(*data)->uu[ii] = 0.0;}
            if ((*data)->uu[smap->iMjc[ii]] < 0)    {(*data)->uu[smap->iMjc[ii]] = 0.0;}
            if ((*data)->vv[ii] > 0)    {(*data)->vv[ii] = 0.0;}
            if ((*data)->vv[smap->icjM[ii]] < 0)    {(*data)->vv[smap->icjM[ii]] = 0.0;}
        }
        // apply the cfl limiter
        if ((*data)->cfl_active[ii] == 1)
        {
            (*data)->uu[ii] = 0.0;
            (*data)->uu[smap->iMjc[ii]] = 0.0;
            (*data)->vv[ii] = 0.0;
            (*data)->vv[smap->icjM[ii]] = 0.0;
            (*data)->cfl_active[ii] = 0;
        }
    }
    // update flow rates
    for (ii = 0; ii < param->n2ci; ii++)
    {
        (*data)->Fu[ii] = (*data)->uu[ii] * (*data)->Asx[ii];
        (*data)->Fv[ii] = (*data)->vv[ii] * (*data)->Asy[ii];
        (*data)->cflx[ii] = fabs((*data)->uu[ii] * param->dt / param->dx);
        (*data)->cfly[ii] = fabs((*data)->vv[ii] * param->dt / param->dy);
        if ((*data)->cflx[ii] > 1 | (*data)->cfly[ii] > 1)
        {printf("WARNING: CFL = %f, %f for cell (%d,%d) of rank %d!\n",(*data)->cflx[ii],(*data)->cfly[ii],smap->ii[ii],smap->jj[ii],irank);}
    }

}

// >>>>> correct velocity for waterfall
void waterfall_velocity(Data **data, Map *smap, Config *param)
{
    int ii;
    double Cw = 0.7, cfl_max = 0.24;
    for (ii = 0; ii < param->n2ci; ii++)
    {
        if ((*data)->wtfx[ii] != 0)
        {
            // Weir equation
            // (*data)->uu[ii] = (*data)->wtfx[ii] * Cw * sqrt(2.0*param->grav*(*data)->Asx[ii]/param->dy);
            // Kinematic Wave
            if ((*data)->wtfx[ii] == 1)
            {
                (*data)->uu[smap->iMjc[ii]] = (*data)->wtfx[ii] * Cw * sqrt(2.0*param->grav*(*data)->Asx[smap->iMjc[ii]]/param->dy);
                // (*data)->uu[smap->iMjc[ii]] = pow(((*data)->bottom[smap->iMjc[ii]] - (*data)->bottom[ii])/param->dx, 0.5)
                //     * (*data)->Asx[smap->iMjc[ii]]/param->dy/param->manning;
                if (fabs((*data)->uu[smap->iMjc[ii]]) > cfl_max * param->dx / param->dt)
                {(*data)->uu[smap->iMjc[ii]] = (*data)->wtfx[ii] * cfl_max * param->dx / param->dt;}
            }
            else if ((*data)->wtfx[ii] == -1)
            {
                (*data)->uu[ii] = (*data)->wtfx[ii] * Cw * sqrt(2.0*param->grav*(*data)->Asx[ii]/param->dy);
                // (*data)->uu[ii] = -pow(((*data)->bottom[smap->iPjc[ii]] - (*data)->bottom[ii])/param->dx, 0.5)
                //     * (*data)->Asx[ii]/param->dy/param->manning;
                if (fabs((*data)->uu[ii]) > cfl_max * param->dx / param->dt)
                {(*data)->uu[ii] = (*data)->wtfx[ii] * cfl_max * param->dx / param->dt;}
            }
        }
        if ((*data)->wtfy[ii] != 0)
        {
            // (*data)->vv[ii] = (*data)->wtfy[ii] * Cw * sqrt(2.0*param->grav*(*data)->Asy[ii]/param->dx);
            if ((*data)->wtfy[ii] == 1)
            {
                // (*data)->vv[smap->icjM[ii]] = (*data)->wtfy[ii] * Cw * sqrt(2.0*param->grav*(*data)->Asy[smap->icjM[ii]]/param->dx);
                (*data)->vv[smap->icjM[ii]] = pow(((*data)->bottom[smap->icjM[ii]] - (*data)->bottom[ii])/param->dy, 0.5)
                    * (*data)->Asy[smap->icjM[ii]]/param->dx/param->manning;
                if (fabs((*data)->vv[smap->icjM[ii]]) > cfl_max * param->dy / param->dt)
                {(*data)->vv[smap->icjM[ii]] = (*data)->wtfy[ii] * cfl_max * param->dy / param->dt;}
            }
            else if ((*data)->wtfy[ii] == -1)
            {
                (*data)->vv[ii] = (*data)->wtfy[ii] * Cw * sqrt(2.0*param->grav*(*data)->Asy[ii]/param->dx);
                // (*data)->vv[ii] = -pow(((*data)->bottom[smap->icjP[ii]] - (*data)->bottom[ii])/param->dy, 0.5)
                    // * (*data)->Asy[ii]/param->dx/param->manning;
                if (fabs((*data)->vv[ii]) > cfl_max * param->dy / param->dt)
                {(*data)->vv[ii] = (*data)->wtfy[ii] * cfl_max * param->dy / param->dt;}
            }

        }
    }
}

// >>>>> Enforce velocity boundary conditions
void enforce_velo_bc(Data **data, Map *smap, Config *param, int irank, int nrank)
{
    int ii, kk;
    double Fw, Ftot, apply_tide;
    // update velocity at ghost cells
    for (ii = 0; ii < param->n2ci; ii++)
    {
        // x boundary
        if (smap->ii[ii] == 0)
        {
            (*data)->uu[smap->iMjc[ii]] = (*data)->uu[ii];
            (*data)->vv[smap->iMjc[ii]] = (*data)->vv[ii];
        }
        else if (smap->ii[ii] == param->nx-1)
        {
            (*data)->uu[smap->iPjc[ii]] = (*data)->uu[ii];
            (*data)->vv[smap->iPjc[ii]] = (*data)->vv[ii];
        }
        // y boundary
        if (smap->jj[ii] == 0)
        {
            (*data)->uu[smap->icjM[ii]] = (*data)->uu[ii];
            (*data)->vv[smap->icjM[ii]] = (*data)->vv[ii];
        }
        else if (smap->jj[ii] == param->ny-1)
        {
            (*data)->uu[smap->icjP[ii]] = (*data)->uu[ii];
            (*data)->vv[smap->icjP[ii]] = (*data)->vv[ii];
        }
    }
    // update velocity at tidal boundary
    for (ii = 0; ii < param->n2ci; ii++)
    {
        // check for tidal boundary
        apply_tide = 0;
        // check if on tidal boundary
        for (kk = 0; kk < param->n_tide; kk++)
        {
            if (smap->ii[ii] >= param->tide_locX[2*kk] & smap->ii[ii] <= param->tide_locX[2*kk+1])
            {if (smap->jj[ii] >= param->tide_locY[2*kk] & smap->jj[ii] <= param->tide_locY[2*kk+1])  {apply_tide = 1;}}
        }
        // update velocity
        Fw = ((*data)->eta[ii] - (*data)->etan[ii]) * (*data)->Asz[ii] / param->dt;
        if (smap->ii[ii] == 0 & irank % param->mpi_nx == 0 & apply_tide == 1)
        {
            Ftot = -(*data)->Fu[ii] - (*data)->Fv[ii] + (*data)->Fv[smap->icjM[ii]] - Fw;
            if ((*data)->Asx[ii] > 0)   {(*data)->uu[smap->iMjc[ii]] = -Ftot / (*data)->Asx[ii];}
            else    {(*data)->uu[smap->iMjc[ii]] = 0.0;}
        }
        if (smap->ii[ii] == param->nx-1 & irank%param->mpi_nx == param->mpi_nx-1 & apply_tide == 1)
        {
            Ftot = (*data)->Fu[smap->iMjc[ii]] - (*data)->Fv[ii] + (*data)->Fv[smap->icjM[ii]] - Fw;
            if ((*data)->Asx[ii] > 0)   {(*data)->uu[ii] = Ftot / (*data)->Asx[ii];}
            else    {(*data)->uu[ii] = 0.0;}
        }
        if (smap->jj[ii] == 0 & irank < param->mpi_nx & apply_tide == 1)
        {
            Ftot = -(*data)->Fv[ii] - (*data)->Fu[ii] + (*data)->Fu[smap->iMjc[ii]] - Fw;
            if ((*data)->Asy[ii] > 0)   {(*data)->vv[smap->icjM[ii]] = -Ftot / (*data)->Asy[ii];}
            else    {(*data)->vv[smap->icjM[ii]] = 0.0;}
        }
        if (smap->jj[ii] == param->ny-1 & irank > param->mpi_nx*(param->mpi_ny-1) & apply_tide == 1)
        {
            Ftot = (*data)->Fv[smap->icjM[ii]] - (*data)->Fu[ii] + (*data)->Fu[smap->iMjc[ii]] - Fw;
            if ((*data)->Asy[ii] > 0)   {(*data)->vv[ii] = Ftot / (*data)->Asy[ii];}
            else    {(*data)->vv[ii] = 0.0;}
        }
    }
}

// >>>>> interpolate velocity to get vx and uy
void interp_velocity(Data **data, Map *smap, Config *param)
{
    int ii;
    for (ii = 0; ii < param->n2ci; ii++)
    {
        // (*data)->uy[ii] = 0.25 * ((*data)->uu[ii] + (*data)->uu[smap->iMjc[ii]] + \
        //         (*data)->uu[smap->icjP[ii]] + (*data)->uu[smap->icjP[ii]+param->nx]);
        (*data)->uy[ii] = 0.25 * ((*data)->uu[ii] + (*data)->uu[smap->iMjc[ii]] + \
                (*data)->uu[smap->icjP[ii]] + (*data)->uu[smap->icjP[ii]-1]);
        (*data)->vx[ii] = 0.25 * ((*data)->vv[ii] + (*data)->vv[smap->icjM[ii]] + \
                (*data)->vv[smap->iPjc[ii]] + (*data)->vv[smap->iPjc[ii]-param->nx]);
    }
}

// >>>>> update drag coefficient
void update_drag_coef(Data **data, Config *param)
{
    double coef, effh, expo, kk = 0.41;
    int ii;
    coef = param->grav * param->manning * param->manning;
    for (ii = 0; ii < param->n2ci; ii++)
    {
        if ((*data)->Vs[ii] > 0.0)
        {
            effh = (*data)->Vs[ii] / (param->dx * param->dy);
            // apply the thin-layer drag model
            if (effh < param->hD)   {expo = 2.0 / 3.0;}
            else    {expo = 1.0 / 3.0;}
            (*data)->CDx[ii] = coef / pow(effh, expo);
            (*data)->CDy[ii] = coef / pow(effh, expo);
        }
    }
}

// >>>>> update cell volume using flux
void volume_by_flux(Data **data, Map *smap, Config *param)
{
    int ii, jj, kk;
    for (ii = 0; ii < param->n2ci; ii++)
    {
        // volume change by flux
        (*data)->Vflux[ii] = (*data)->Vsn[ii] + param->dt * \
          ((*data)->Fu[smap->iMjc[ii]] - (*data)->Fu[ii] + (*data)->Fv[smap->icjM[ii]] - (*data)->Fv[ii]);
        // subsurface source
        if (param->sim_groundwater == 1)
        {
            if ((*data)->qseepage[ii] < 0 & (*data)->Vflux[ii] > -(*data)->qseepage[ii]*param->dt*param->wcs*(*data)->Asz[ii])
            {(*data)->Vflux[ii] += (*data)->qseepage[ii] * param->dt * param->wcs * (*data)->Asz[ii];}
            else if ((*data)->qseepage[ii] > 0)
            {
                if ((*data)->dept[ii] > 0)
                {(*data)->Vflux[ii] += (*data)->qseepage[ii] * param->dt * param->wcs * (*data)->Asz[ii];}
                else if ((*data)->qseepage[ii]*param->dt*param->wcs > param->min_dept)
                {(*data)->Vflux[ii] += (*data)->qseepage[ii] * param->dt * param->wcs * (*data)->Asz[ii];}
            }
        }

        // inflow
        if (param->n_inflow > 0)
        {
            for (kk = 0; kk < param->n_inflow; kk++)
            {
                for (jj = 0; jj < (*data)->inflowloc_len[kk]; jj++)
                {
                    if (ii == (*data)->inflowloc[kk][jj] & (*data)->current_inflow[kk] > 0)
                    {(*data)->Vflux[ii] += (*data)->current_inflow[kk] * param->dt / (*data)->inflowloc_len[kk];}
                }
            }
        }
    }

}


// >>>>> update subgrid variables
void update_subgrid_variable(Data **data, Map *smap, Config *param, int irank, int nrank)
{
    int ii, jj, kk;
    if (param->use_subgrid == 1)
    {
        for (ii = 0; ii < param->n2ct; ii++)    {(*data)->Vsn[ii] = (*data)->Vs[ii];}
        subgrid_index(data, smap, param);
        subgrid_interp_and_combine(data, smap, param, irank, nrank);
    }
    else
    {
        for (ii = 0; ii < param->n2ct; ii++)
        {
            (*data)->Vsn[ii] = (*data)->Vs[ii];
            (*data)->Vs[ii] = (*data)->dept[ii] * param->dx * param->dy;
            if ((*data)->dept[ii] > 0) {(*data)->Asz[ii] = param->dx * param->dy;}
            else    {(*data)->Asz[ii] = 0.0;}
            (*data)->Asx[ii] = (*data)->deptx[ii] * param->dy;
            (*data)->Asy[ii] = (*data)->depty[ii] * param->dx;
            if (param->nx == 1) {(*data)->Asx[ii] = 0.0;}
            if (param->ny == 1) {(*data)->Asy[ii] = 0.0;}
        }
        for (ii = 0; ii < param->n2ci; ii++)
        {
            (*data)->Vsx[ii] = 0.5 * ((*data)->Vs[ii] + (*data)->Vs[smap->iPjc[ii]]);
            (*data)->Vsy[ii] = 0.5 * ((*data)->Vs[ii] + (*data)->Vs[smap->icjP[ii]]);
            (*data)->Aszx[ii] = 0.5 * ((*data)->Asz[ii] + (*data)->Asz[smap->iPjc[ii]]);
            (*data)->Aszy[ii] = 0.5 * ((*data)->Asz[ii] + (*data)->Asz[smap->icjP[ii]]);
        }
    }
    // subgrid variables along boundaries
    for (ii = 0; ii < param->nx; ii++)
    {
        (*data)->Vs[smap->jMou[ii]] = (*data)->Vs[smap->jMin[ii]];
        (*data)->Vs[smap->jPou[ii]] = (*data)->Vs[smap->jPin[ii]];
        (*data)->Asx[smap->jMou[ii]] = (*data)->Asx[smap->jMin[ii]];
        (*data)->Asx[smap->jPou[ii]] = (*data)->Asx[smap->jPin[ii]];
        (*data)->Asy[smap->jMou[ii]] = (*data)->Asy[smap->jMin[ii]];
        (*data)->Asy[smap->jPou[ii]] = (*data)->Asy[smap->jPin[ii]];
        (*data)->Asz[smap->jMou[ii]] = (*data)->Asz[smap->jMin[ii]];
        (*data)->Asz[smap->jPou[ii]] = (*data)->Asz[smap->jPin[ii]];
    }
    for (ii = 0; ii < param->ny; ii++)
    {
        (*data)->Vs[smap->iMou[ii]] = (*data)->Vs[smap->iMin[ii]];
        (*data)->Vs[smap->iPou[ii]] = (*data)->Vs[smap->iPin[ii]];
        (*data)->Asx[smap->iMou[ii]] = (*data)->Asx[smap->iMin[ii]];
        (*data)->Asx[smap->iPou[ii]] = (*data)->Asx[smap->iPin[ii]];
        (*data)->Asy[smap->iMou[ii]] = (*data)->Asy[smap->iMin[ii]];
        (*data)->Asy[smap->iPou[ii]] = (*data)->Asy[smap->iPin[ii]];
        (*data)->Asz[smap->iMou[ii]] = (*data)->Asz[smap->iMin[ii]];
        (*data)->Asz[smap->iPou[ii]] = (*data)->Asz[smap->iPin[ii]];
    }
    for (ii = 0; ii < param->nx; ii++)
    {
        (*data)->Vsy[smap->jMou[ii]] = (*data)->Vs[smap->jMin[ii]];
        (*data)->Aszy[smap->jMou[ii]] = (*data)->Asz[smap->jMin[ii]];
    }
    for (ii = 0; ii < param->ny; ii++)
    {
        (*data)->Vsx[smap->iMou[ii]] = (*data)->Vs[smap->iMin[ii]];
        (*data)->Aszx[smap->iMou[ii]] = (*data)->Asz[smap->iMin[ii]];
    }
}

// >>>>> Search for the index in pre-stored surface elevations
void subgrid_index(Data **data, Map *smap, Config *param)
{
    int ii, jj;
    for (ii = 0; ii < param->n2ci; ii++)
    {
        // if eta is increased to the next interval
        if ((*data)->eta[ii] > (*data)->layers_sub[(*data)->eta_ind[ii]+1])
        {
            for (jj = (*data)->eta_ind[ii]+1; jj < param->nlay_sub; jj++)
            {
                if ((*data)->eta[ii] < (*data)->layers_sub[jj])
                {(*data)->eta_ind[ii] = jj-1;    break;}
                (*data)->eta_ind[ii] = param->nlay_sub - 1;
            }
        }
        // if eta is decrease to the interval below
        else if ((*data)->eta[ii] < (*data)->layers_sub[(*data)->eta_ind[ii]])
        {
            for (jj = (*data)->eta_ind[ii]; jj > 0; jj--)
            {
                if ((*data)->eta[ii] > (*data)->layers_sub[jj])
                {(*data)->eta_ind[ii] = jj;     break;}
                (*data)->eta_ind[ii] = 0;
            }
        }
    }
}

// >>>>> Calculate subgrid variables for current surface elevation
void subgrid_interp_and_combine(Data **data, Map *smap, Config *param, int irank, int nrank)
{
    int ii, jj, kk;
    for (ii = 0; ii < param->n2ci; ii++)
    {
        // interpolate to update subgrid variables
        kk = (*data)->eta_ind[ii];
        if (kk == 0 | kk == param->nlay_sub-1)
        {
            // if surf out of bound, no interpolation is needed
            (*data)->Vxp[ii] = (*data)->Vxp_sub[kk][ii];
            (*data)->Vxm[ii] = (*data)->Vxm_sub[kk][ii];
            (*data)->Vyp[ii] = (*data)->Vyp_sub[kk][ii];
            (*data)->Vym[ii] = (*data)->Vym_sub[kk][ii];
            (*data)->Axp[ii] = (*data)->Axp_sub[kk][ii];
            (*data)->Axm[ii] = (*data)->Axm_sub[kk][ii];
            (*data)->Ayp[ii] = (*data)->Ayp_sub[kk][ii];
            (*data)->Aym[ii] = (*data)->Aym_sub[kk][ii];
            (*data)->Asz[ii] = (*data)->Asz_sub[kk][ii];
        }
        else
        {
            (*data)->Vxp[ii] = interp_sub((*data)->layers_sub, (*data)->Vxp_sub, (*data)->eta[ii], ii, kk);
            (*data)->Vxm[ii] = interp_sub((*data)->layers_sub, (*data)->Vxm_sub, (*data)->eta[ii], ii, kk);
            (*data)->Vyp[ii] = interp_sub((*data)->layers_sub, (*data)->Vyp_sub, (*data)->eta[ii], ii, kk);
            (*data)->Vym[ii] = interp_sub((*data)->layers_sub, (*data)->Vym_sub, (*data)->eta[ii], ii, kk);
            (*data)->Axp[ii] = interp_sub((*data)->layers_sub, (*data)->Axp_sub, (*data)->eta[ii], ii, kk);
            (*data)->Axm[ii] = interp_sub((*data)->layers_sub, (*data)->Axm_sub, (*data)->eta[ii], ii, kk);
            (*data)->Ayp[ii] = interp_sub((*data)->layers_sub, (*data)->Ayp_sub, (*data)->eta[ii], ii, kk);
            (*data)->Aym[ii] = interp_sub((*data)->layers_sub, (*data)->Aym_sub, (*data)->eta[ii], ii, kk);
            (*data)->Asz[ii] = interp_sub((*data)->layers_sub, (*data)->Asz_sub, (*data)->eta[ii], ii, kk);
            // check for inundation of each face
            if ((*data)->bottomXP[ii] > (*data)->layers_sub[kk] & (*data)->bottomXP[ii] <= (*data)->layers_sub[kk+1])
            {
                // xp
                if ((*data)->eta[ii] <= (*data)->bottomXP[ii])
                {(*data)->Axp[ii] = 0.0;}
                else
                {(*data)->Axp[ii] = (*data)->Axp_sub[kk][ii] * \
                    ((*data)->eta[ii]-(*data)->bottomXP[ii]) / ((*data)->layers_sub[kk]-(*data)->bottomXP[ii]);}
            }
            if ((*data)->bottomXP[smap->iMjc[ii]] > (*data)->layers_sub[kk] & (*data)->bottomXP[smap->iMjc[ii]] <= (*data)->layers_sub[kk+1])
            {
                // xm
                if ((*data)->eta[ii] <= (*data)->bottomXP[smap->iMjc[ii]])
                {(*data)->Axm[ii] = 0.0;}
                else
                {(*data)->Axm[ii] = (*data)->Axm_sub[kk][ii] * \
                    ((*data)->eta[ii]-(*data)->bottomXP[smap->iMjc[ii]]) / ((*data)->layers_sub[kk]-(*data)->bottomXP[smap->iMjc[ii]]);}
            }
            if ((*data)->bottomYP[ii] > (*data)->layers_sub[kk] & (*data)->bottomYP[ii] <= (*data)->layers_sub[kk+1])
            {
                // yp
                if ((*data)->eta[ii] <= (*data)->bottomYP[ii])
                {(*data)->Ayp[ii] = 0.0;}
                else
                {(*data)->Ayp[ii] = (*data)->Ayp_sub[kk][ii] * \
                    ((*data)->eta[ii]-(*data)->bottomYP[ii]) / ((*data)->layers_sub[kk]-(*data)->bottomYP[ii]);}
            }
            if ((*data)->bottomYP[smap->icjM[ii]] > (*data)->layers_sub[kk] & (*data)->bottomYP[smap->icjM[ii]] <= (*data)->layers_sub[kk+1])
            {
                // ym
                if ((*data)->eta[ii] <= (*data)->bottomYP[smap->icjM[ii]])
                {(*data)->Aym[ii] = 0.0;}
                else
                {(*data)->Aym[ii] = (*data)->Aym_sub[kk][ii] * \
                    ((*data)->eta[ii]-(*data)->bottomYP[smap->icjM[ii]]) / ((*data)->layers_sub[kk]-(*data)->bottomYP[smap->icjM[ii]]);}
            }

        }


    }

    // subgrid variables in ghost cells
    for (ii = 0; ii < param->nx; ii++)
    {
        (*data)->Vyp[smap->jMou[ii]] = (*data)->Vyp[smap->jMin[ii]];
        (*data)->Vyp[smap->jPou[ii]] = (*data)->Vyp[smap->jPin[ii]];
        (*data)->Vym[smap->jMou[ii]] = (*data)->Vym[smap->jMin[ii]];
        (*data)->Vym[smap->jPou[ii]] = (*data)->Vym[smap->jPin[ii]];
        (*data)->Ayp[smap->jMou[ii]] = (*data)->Ayp[smap->jMin[ii]];
        (*data)->Ayp[smap->jPou[ii]] = (*data)->Ayp[smap->jPin[ii]];
        (*data)->Aym[smap->jMou[ii]] = (*data)->Aym[smap->jMin[ii]];
        (*data)->Aym[smap->jPou[ii]] = (*data)->Aym[smap->jPin[ii]];
    }
    for (ii = 0; ii < param->ny; ii++)
    {
        (*data)->Vxp[smap->iMou[ii]] = (*data)->Vxp[smap->iMin[ii]];
        (*data)->Vxp[smap->iPou[ii]] = (*data)->Vxp[smap->iPin[ii]];
        (*data)->Vxm[smap->iMou[ii]] = (*data)->Vxm[smap->iMin[ii]];
        (*data)->Vxm[smap->iPou[ii]] = (*data)->Vxm[smap->iPin[ii]];
        (*data)->Axp[smap->iMou[ii]] = (*data)->Axp[smap->iMin[ii]];
        (*data)->Axp[smap->iPou[ii]] = (*data)->Axp[smap->iPin[ii]];
        (*data)->Axm[smap->iMou[ii]] = (*data)->Axm[smap->iMin[ii]];
        (*data)->Axm[smap->iPou[ii]] = (*data)->Axm[smap->iPin[ii]];
    }

    // mpi exchange
    if (param->use_mpi == 1)
    {
        mpi_exchange_surf((*data)->Vxp, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Vxm, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Vyp, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Vym, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Axp, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Axm, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Ayp, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Aym, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Asz, smap, 2, param, irank, nrank);
    }

    // combine subgrid variables
    for (ii = 0; ii < param->n2ci; ii++)
    {
        // volumes
        (*data)->Vs[ii] = (*data)->Vxp[ii] + (*data)->Vxm[ii];
        (*data)->Vsx[ii] = (*data)->Vxp[ii] + (*data)->Vxm[smap->iPjc[ii]];
        (*data)->Vsy[ii] = (*data)->Vyp[ii] + (*data)->Vym[smap->icjP[ii]];
        // areas
        (*data)->Asx[ii] = (*data)->Axp[ii];
        if ((*data)->Axm[smap->iPjc[ii]] < (*data)->Axp[ii])
        {(*data)->Asx[ii] = (*data)->Axm[smap->iPjc[ii]];}
        (*data)->Asy[ii] = (*data)->Ayp[ii];
        if ((*data)->Aym[smap->icjP[ii]] < (*data)->Ayp[ii])
        {(*data)->Asy[ii] = (*data)->Aym[smap->icjP[ii]];}
    }

    // other subgrid variables
    for (ii = 0; ii < param->n2ci; ii++)
    {
        if ((*data)->Asx[ii] > 0)
        {(*data)->Aszx[ii] = 0.5 * ((*data)->Asz[ii] + (*data)->Asz[smap->iPjc[ii]]);}
        else
        {(*data)->Aszx[ii] = 0.5 * (*data)->Asz[ii];}
        if ((*data)->Asy[ii] > 0)
        {(*data)->Aszy[ii] = 0.5 * ((*data)->Asz[ii] + (*data)->Asz[smap->icjP[ii]]);}
        else
        {(*data)->Aszy[ii] = 0.5 * (*data)->Asz[ii];}
    }
    // final check
    for (ii = 0; ii < param->n2ct; ii++)
    {
        if ((*data)->Vs[ii] <= 0.0)  {(*data)->Vs[ii] = 0.0;}
        if ((*data)->Vsx[ii] <= 0.0)  {(*data)->Vsx[ii] = 0.0;}
        if ((*data)->Vsy[ii] <= 0.0)  {(*data)->Vsy[ii] = 0.0;}
        if ((*data)->Asx[ii] <= 0.0)  {(*data)->Asx[ii] = 0.0;}
        if ((*data)->Asy[ii] <= 0.0)  {(*data)->Asy[ii] = 0.0;}
        if ((*data)->Asz[ii] <= 0.0)  {(*data)->Asz[ii] = 0.0;}
    }
}
