// Functions for groundwater solver
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

#include "linsys.h"

// #include "laspack/vector.h"
// #include "laspack/qmatrix.h"
// #include "laspack/rtc.h"

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

void solve_groundwater(Data **data, Lisy **sysm, Map *smap, Map *gmap, Config *param, int irank, int nrank);
void compute_K_face(Data **data, Map *gmap, Config *param, int irank, int nrank);
void baroclinic_face(Data **data, Map *gmap, Config *param, int irank, int nrank);
void groundwater_mat_coeff(Data **data, Map *gmap, Config *param, int irank);
void groundwater_rhs(Data **data, Map *gmap, Config *param, int irank);
void build_groundwater_system(Data *data, Map *gmap, Config *param, QMatrix A, Vector b);
void solve_groundwater_system(Data **data, Map *gmap, QMatrix A, Vector b, Vector x, Config *param);
void enforce_head_bc(Data **data, Map *gmap, Config *param);
void groundwater_flux(Data **data, Map *gmap, Config *param, int irank);
void check_room(Data **data, Map *gmap, Config *param);
void update_water_content(Data **data, Map *gmap, Config *param);
void enforce_moisture_bc(Data **data, Map *gmap, Config *param);
void reallocate_water_content(Data **data, Map *gmap, Config *param, int irank);
int check_adj_sat(Data *data, Map *gmap, Config *param, int ii);
void check_head_gradient(Data **data, Map *gmap, Config *param, int ii);
double allocate_send(Data **data, Map *gmap, Config *param, int ii, double dV);
double allocate_recv(Data **data, Map *gmap, Config *param, int ii, double dV);
void volume_by_flux_subs(Data **data, Map *gmap, Config *param);
void adaptive_time_step(Data *data, Map *gmap, Config **param, int root, int irank);

// >>>>> Top level groundwater solver <<<<<
void solve_groundwater(Data **data, Lisy **sysm, Map *smap, Map *gmap, Config *param, int irank, int nrank)
{
    int ii;
    // allocate linear system
    QMatrix A;
    Q_Constr(&A, "A", param->n3ci, False, Rowws, Normal, True);
    Vector b;
    V_Constr(&b, "b", param->n3ci, Normal, True);
    Vector x;
    V_Constr(&x, "x", param->n3ci, Normal, True);

    if ((*data)->repeat[0] == 0)
    {
        for (ii = 0; ii < param->n3ct; ii++)
        {
            (*data)->hn[ii] = (*data)->h[ii];
            (*data)->hwc[ii] = compute_hwc(*data, ii, param);
            (*data)->wcn[ii] = (*data)->wc[ii];
            (*data)->wch[ii] = compute_wch(*data, ii, param);
            (*data)->ch[ii] = compute_ch(*data, ii, param);

        }
    }
    else
    {
        for (ii = 0; ii < param->n3ct; ii++)
        {
            (*data)->h[ii] = (*data)->hn[ii];
            (*data)->hwc[ii] = compute_hwc(*data, ii, param);
            (*data)->wc[ii] = (*data)->wcn[ii];
            (*data)->wch[ii] = compute_wch(*data, ii, param);
            (*data)->ch[ii] = compute_ch(*data, ii, param);
        }
    }

    if (param->use_mpi == 1)
    {
        mpi_exchange_subsurf((*data)->hn, gmap, 2, param, irank, nrank);
        mpi_exchange_subsurf((*data)->wcn, gmap, 2, param, irank, nrank);
    }

    // >>> Predictor step
    compute_K_face(data, gmap, param, irank, nrank);
    baroclinic_face(data, gmap, param, irank, nrank);
    groundwater_mat_coeff(data, gmap, param, irank);
    groundwater_rhs(data, gmap, param, irank);
    build_groundwater_system(*data, gmap, param, A, b);
    // solve_groundwater_system(data, gmap, A, b, x, param);
    build_crs_matrix(*data, sysm, param, 1);
    cgsolve(sysm, A, (*data)->Grhs, (*data)->h, param->n3ci);
    for (ii = 0; ii < 5; ii++) {printf("(%d,%d) : h=%f\n",gmap->jj[ii],gmap->kk[ii],(*data)->h[ii]);}
    enforce_head_bc(data, gmap, param);
    if (param->use_mpi == 1)
    {mpi_exchange_subsurf((*data)->h, gmap, 2, param, irank, nrank);}

    // >>> Corrector step
    if (param->use_corrector == 1)
    {
        compute_K_face(data, gmap, param, irank, nrank);
        groundwater_flux(data, gmap, param, irank);
        check_room(data, gmap, param);
        update_water_content(data, gmap, param);
    }
    else
    {
        groundwater_flux(data, gmap, param, irank);
        for (ii = 0; ii < param->n3ci; ii++)
        // {(*data)->wc[ii] = compute_wch(*data, ii, param);}
        {(*data)->wc[ii] = (*data)->wcs[ii];}
    }
    for (ii = 0; ii < param->n3ci; ii++)
    {
        (*data)->hwc[ii] = compute_hwc(*data, ii, param);
        (*data)->wch[ii] = compute_wch(*data, ii, param);
    }
    if (param->use_mpi == 1)
    {
        mpi_exchange_subsurf((*data)->wc, gmap, 2, param, irank, nrank);
        mpi_exchange_subsurf((*data)->wch, gmap, 2, param, irank, nrank);
        mpi_exchange_subsurf((*data)->hwc, gmap, 2, param, irank, nrank);
    }

    // >>> Post-allocation step
    if (param->use_corrector == 1)
    {reallocate_water_content(data, gmap, param, irank);}

    (*data)->qbc[0] -= (*data)->qz[0] * param->dt;
    (*data)->qbc[1] -= (*data)->qz[param->nz-1] * param->dt;

    // for (ii = 0; ii < param->n3ci; ii++)
    // {
    //     if (gmap->istop[ii] == 1 & gmap->kk[ii] > 1)
    //     {
    //         printf(" (jj,kk)=(%d,%d)  :  (h,wc)=(%f,%f)  (d,hbc)=(%f,%f)  q=(%f,%f,%f,%f)  s=%f\n",
    //             gmap->jj[ii],gmap->kk[ii],(*data)->h[ii],(*data)->wc[ii],
    //             (*data)->dept[gmap->top2d[ii]],(*data)->h[gmap->icjPkc[ii]],
    //             (*data)->qy[ii],(*data)->qy[gmap->icjMkc[ii]],(*data)->qz[ii],(*data)->qz[gmap->icjckM[ii]],
    //             (*data)->s_subs[0][ii]);
    //     }
    // }
    // printf("  --------------------  \n");

    // >>> Final check of water content
    volume_by_flux_subs(data, gmap, param);
    for (ii = 0; ii < param->n3ci; ii++)
    {
        if (gmap->actv[ii] == 1)
        {
            if ((*data)->wc[ii] > (*data)->wcs[ii])
            {
                (*data)->vloss[ii] += ((*data)->wc[ii] - (*data)->wcs[ii]) * param->dx * param->dy * gmap->dz3d[ii];
                (*data)->wc[ii] = (*data)->wcs[ii];
            }
            else if ((*data)->wc[ii] < (*data)->wcr[ii])
            {
                (*data)->vloss[ii] -= ((*data)->wcr[ii] - (*data)->wc[ii]) * param->dx * param->dy * gmap->dz3d[ii];
                (*data)->wc[ii] = (*data)->wcr[ii];
            }
            // update cell volume
            (*data)->Vgn[ii] = (*data)->Vg[ii];
            (*data)->Vg[ii] = param->dx*param->dy*gmap->dz3d[ii]*(*data)->wc[ii];
        }
        else
        {(*data)->wc[ii] = 0.0; (*data)->h[ii] = -100.0;    (*data)->Vg[ii] = 0.0;}

    }
    enforce_moisture_bc(data, gmap, param);
    if (param->use_mpi == 1)
    {mpi_exchange_subsurf((*data)->wc, gmap, 2, param, irank, nrank);}
    // >>> Adaptive time stepping
    if (param->dt_adjust == 1)
    {adaptive_time_step(*data, gmap, &param, 0, irank);}
    // free memory
    Q_Destr(&A);
    V_Destr(&b);
    V_Destr(&x);

}

// >>>>> Compute hydraulic conductivity on cell faces <<<<<
void compute_K_face(Data **data, Map *gmap, Config *param, int irank, int nrank)
{
    int ii, ip, im, jp, jm;
    double Kp, Km, coef;
    // Check for lateral coupling between surface and subsurface domains
    if (param->sim_shallowwater == 1)
    {
        for (ii = 0; ii < param->n3ci; ii++)
        {
            ip = gmap->iPjckc[ii];
            im = gmap->iMjckc[ii];
            jp = gmap->icjPkc[ii];
            jm = gmap->icjMkc[ii];
            // Kx
            if (gmap->actv[ii] == 1 & gmap->actv[ip] == 0 & (*data)->dept[gmap->top2d[ip]] > 0.0)
            {if ((*data)->eta[gmap->top2d[ip]] - (*data)->offset[0] > gmap->bot3d[ii] + 0.5*gmap->dz3d[ii])  {gmap->actvXp[ii] = 1;}}
            if (gmap->actv[ii] == 1 & gmap->actv[im] == 0 & (*data)->dept[gmap->top2d[im]] > 0.0)
            {if ((*data)->eta[gmap->top2d[im]] - (*data)->offset[0] > gmap->bot3d[ii] + 0.5*gmap->dz3d[ii])  {gmap->actvXp[im] = 1;}}
            // Ky
            if (gmap->actv[ii] == 1 & gmap->actv[jp] == 0 & (*data)->dept[gmap->top2d[jp]] > 0.0)
            {if ((*data)->eta[gmap->top2d[jp]] - (*data)->offset[0] > gmap->bot3d[ii] + 0.5*gmap->dz3d[ii])  {gmap->actvYp[ii] = 1;}}
            if (gmap->actv[ii] == 1 & gmap->actv[jm] == 0 & (*data)->dept[gmap->top2d[jm]] > 0.0)
            {if ((*data)->eta[gmap->top2d[jm]] - (*data)->offset[0] > gmap->bot3d[ii] + 0.5*gmap->dz3d[ii])  {gmap->actvYp[jm] = 1;}}
        }
    }
    // conductivities for interior cells
    for (ii = 0; ii < param->n3ci; ii++)
    {
        // Kx
        Kp = compute_K(*data, (*data)->Ksx, gmap->iPjckc[ii], param);
        Km = compute_K(*data, (*data)->Ksx, ii, param);
        (*data)->Kx[ii] = 0.5 * (Kp + Km);
        if (gmap->actv[ii] == 0 | gmap->actv[gmap->iPjckc[ii]] == 0)    {(*data)->Kx[ii] = 0.0;}
        if (gmap->actvXp[ii] == 1)
        {
            if (gmap->actv[ii] == 1)    {(*data)->Kx[ii] = Km;}
            else {(*data)->Kx[ii] = Kp;}
        }
        // Ky
        Kp = compute_K(*data, (*data)->Ksy, gmap->icjPkc[ii], param);
        Km = compute_K(*data, (*data)->Ksy, ii, param);
        (*data)->Ky[ii] = 0.5 * (Kp + Km);
        if (gmap->jj[ii] == param->ny-1 & irank >= param->mpi_nx*(param->mpi_ny-1))
        {
            if (param->bctype_GW[2] == 0)   {(*data)->Ky[ii] = 0.0;}
            else if (param->bctype_GW[2] == 1)  {(*data)->Ky[ii] = Kp;}
        }
        if (gmap->actv[ii] == 0 | gmap->actv[gmap->icjPkc[ii]] == 0)    {(*data)->Ky[ii] = 0.0;}
        if (gmap->actvYp[ii] == 1)
        {
            if (gmap->actv[ii] == 1)    {(*data)->Ky[ii] = Km;}
            else {(*data)->Ky[ii] = Kp;}
        }
        // Kz
        Kp = compute_K(*data, (*data)->Ksz, gmap->icjckP[ii], param);
        Km = compute_K(*data, (*data)->Ksz, ii, param);
        if (gmap->istop[gmap->icjckP[ii]] == 1) {(*data)->Kz[ii] = Kp;}
        else if (gmap->actv[ii] == 0)   {(*data)->Kz[ii] = 0.0;}
        else if (gmap->icjckP[ii] > param->n3ci)    {(*data)->Kz[ii] = Km;}
        else    {(*data)->Kz[ii] = 0.5 * (Kp + Km);}
    }
    // conductivities on iM, jM, kM faces
    for (ii = 0; ii < param->ny*param->nz; ii++)
    {
        Kp = compute_K(*data, (*data)->Ksx, gmap->iMin[ii], param) ;
        Km = compute_K(*data, (*data)->Ksx, gmap->iMou[ii], param) ;
        (*data)->Kx[gmap->iMou[ii]] = 0.5 * (Kp + Km);
        if (param->bctype_GW[0] == 0 & (irank+1) % param->mpi_nx == 0)
        {(*data)->Kx[gmap->iPin[ii]] = 0;}
        if (param->bctype_GW[1] == 0 & irank % param->mpi_nx == 0)
        {(*data)->Kx[gmap->iMou[ii]] = 0;}
        if (gmap->actv[gmap->iMin[ii]] == 0 | gmap->actv[gmap->iMou[ii]] == 0)
        {(*data)->Kx[gmap->iMou[ii]] = 0;}

    }
    for (ii = 0; ii < param->nx*param->nz; ii++)
    {
        Kp = compute_K(*data, (*data)->Ksy, gmap->jMin[ii], param);
        Km = compute_K(*data, (*data)->Ksy, gmap->jMou[ii], param);
        (*data)->Ky[gmap->jMou[ii]] = 0.5 * (Kp + Km);
        if (param->bctype_GW[3] == 0 & irank < param->mpi_nx)
        {(*data)->Ky[gmap->jMou[ii]] = 0;}
        else if (param->bctype_GW[3] == 1 & irank < param->mpi_nx)
        {(*data)->Ky[gmap->jMou[ii]] = Km;}
        if (gmap->actv[gmap->jMin[ii]] == 0 | gmap->actv[gmap->jMou[ii]] == 0)
        {(*data)->Ky[gmap->jMou[ii]] = 0;}
    }
    for (ii = 0; ii < param->nx*param->ny; ii++)
    {
        (*data)->Kz[gmap->kMou[ii]] = 0.0;
        if (param->bctype_GW[4] == 0)   {(*data)->Kz[gmap->kPin[ii]] = 0;}
    }
    for (ii = 0; ii < param->n3ci; ii++)
    {
        if (gmap->istop[ii] == 1)
        {
            Kp = compute_K(*data, (*data)->Ksz, ii, param);
            Km = param->Ksz;
            if (param->sim_shallowwater == 1)
            {
                if ((*data)->dept[gmap->top2d[ii]] > 0)
                {(*data)->Kz[gmap->icjckM[ii]] = 0.5*(Km+Kp);}
                else if (param->bctype_GW[5] == 0)
                {(*data)->Kz[gmap->icjckM[ii]] = 0.0;}
                else
                {(*data)->Kz[gmap->icjckM[ii]] = Kp;}
            }
            else
            {
                if (param->bctype_GW[5] == 1)   {(*data)->Kz[gmap->icjckM[ii]] = Km;}
                else if (param->bctype_GW[5] == 0)  {(*data)->Kz[gmap->icjckM[ii]] = 0.0;}
                else    {(*data)->Kz[gmap->icjckM[ii]] = Kp;}
            }
        }
    }

    // Set K to zero for unsaturated side faces
    if (param->use_full3d == 0)
    {
        for (ii = 0; ii < param->n3ci; ii++)
        {
            if (gmap->actv[ii] == 1 & (*data)->wc[ii] < param->wcs)
            {
                (*data)->Kx[ii] = 0.0;
                (*data)->Kx[gmap->iMjckc[ii]] = 0.0;
                (*data)->Ky[ii] = 0.0;
                (*data)->Ky[gmap->icjMkc[ii]] = 0.0;
            }
        }
    }

    // Kuan 2019
    // for (ii = 0; ii < param->n3ci; ii++)
    // {
    //     if (gmap->istop[ii] == 1)
    //     {if (gmap->jj[ii] < 22)    {(*data)->Kz[gmap->icjckM[ii]] = 0.0;}}
    // }

    // Geng 2015
    // for (ii = 0; ii < param->n3ci; ii++)
    // {
    //     if (gmap->istop[ii] == 1)
    //     {if (gmap->jj[ii] < 20)    {(*data)->Kz[gmap->icjckM[ii]] = 0.0;}}
    //     // if (gmap->istop[ii] == 1)
    //     // {(*data)->Kz[ii] = 0.0;}
    // }

    // for (ii = 0; ii < param->n3ci; ii++)
    // {
    //     if (gmap->jj[ii] == param->ny-1)
    //     {if (gmap->kk[ii] > 15 | gmap->kk[ii] < 9)    {(*data)->Ky[ii] = 0.0;}}
    // }

    // ZhiLi20200821
    // if (strcmp(param->sim_id, " MaxP1"))
    // {
    //     for (ii = 0; ii < param->n3ci; ii++)
    //     {
    //         if (gmap->jj[ii] >= 5)  {(*data)->Ky[ii] = 0.0;}
    //         // avoid exfiltration
    //         if (gmap->ii[ii] == 0 & gmap->istop[ii] == 1)
    //         {(*data)->Kz[gmap->icjckM[ii]] = 0.0;}
    //         if (gmap->ii[ii] == 2 & gmap->istop[ii] == 1)
    //         {(*data)->Kz[gmap->icjckM[ii]] = 0.0;}
    //     }
    //     // {if (gmap->jj[ii] >= 79)  {(*data)->Ky[ii] = 0.0;}}
    // }
    //

    // Maxwell P1
    // for (ii = 0; ii < param->n3ci; ii++)
    // {
    //     if (gmap->jj[ii] >= 5)  {(*data)->Ky[ii] = 0.0;}
    //     else if (gmap->jj[ii] == 0) {(*data)->Ky[ii] = 0.0;}
    // }

    // for (ii = 0; ii < param->n3ci; ii++)
    // {
    //     if (gmap->ii[ii] <= 1 | gmap->ii[ii] >= 13)
    //     {
    //         (*data)->Kx[ii] = 0.0;
    //         (*data)->Kx[gmap->iMjckc[ii]] = 0.0;
    //         (*data)->Ky[ii] = 0.0;
    //         (*data)->Kz[ii] = 0.0;
    //         (*data)->Kz[gmap->icjckM[ii]] = 0.0;
    //     }
    //     if (gmap->jj[ii] >= 19 | gmap->jj[ii] == 0)
    //     {
    //         (*data)->Kx[ii] = 0.0;
    //         (*data)->Ky[ii] = 0.0;
    //         (*data)->Kz[ii] = 0.0;
    //     }
    // }
}

// >>>>> Compute hydraulic conductivity on cell faces <<<<<
void baroclinic_face(Data **data, Map *gmap, Config *param, int irank, int nrank)
{
    int ii;
    double Kp, Km, coef;
    // density effects
    if (param->baroclinic == 1 & param->n_scalar > 0)
    {update_rhovisc(data, gmap, param, irank);}
    // density and viscosity corrections for interior cells
    for (ii = 0; ii < param->n3ci; ii++)
    {
        (*data)->r_rhoxp[ii] = 1.0;
        (*data)->r_rhoyp[ii] = 1.0;
        (*data)->r_rhozp[ii] = 1.0;
        (*data)->r_viscxp[ii] = 1.0;
        (*data)->r_viscyp[ii] = 1.0;
        (*data)->r_visczp[ii] = 1.0;
        // x
        if (gmap->actv[ii] == 1 & gmap->actv[gmap->iPjckc[ii]] == 1)
        {
            (*data)->r_rhoxp[ii] = 0.5 * ((*data)->r_rho[ii] + (*data)->r_rho[gmap->iPjckc[ii]]);
            (*data)->r_viscxp[ii] = 0.5 * ((*data)->r_visc[ii] + (*data)->r_visc[gmap->iPjckc[ii]]);
        }
        // y
        if (gmap->jj[ii] == param->ny-1 & irank >= param->mpi_nx*(param->mpi_ny-1) & param->bctype_GW[2] == 1)
        {
            (*data)->r_rhoyp[ii] = (*data)->r_rho[gmap->icjPkc[ii]];
            (*data)->r_viscyp[ii] = (*data)->r_visc[gmap->icjPkc[ii]];
        }
        else
        {
            if (gmap->actv[ii] == 1 & gmap->actv[gmap->icjPkc[ii]] == 1)
            {
                (*data)->r_rhoyp[ii] = 0.5 * ((*data)->r_rho[ii] + (*data)->r_rho[gmap->icjPkc[ii]]);
                (*data)->r_viscyp[ii] = 0.5 * ((*data)->r_visc[ii] + (*data)->r_visc[gmap->icjPkc[ii]]);
                if (gmap->jj[ii] == param->ny-1 & irank >= param->mpi_nx*(param->mpi_ny-1) & param->bctype_GW[2] == 1)
                {
                    (*data)->r_rhoyp[ii] = (*data)->r_rho[gmap->icjPkc[ii]];
                    (*data)->r_viscyp[ii] = (*data)->r_visc[gmap->icjPkc[ii]];
                }
            }
        }
        // z
        if (gmap->actv[ii] == 1 & gmap->actv[gmap->icjckP[ii]] == 1)
        {
            (*data)->r_rhozp[ii] = 0.5 * ((*data)->r_rho[ii] + (*data)->r_rho[gmap->icjckP[ii]]);
            (*data)->r_visczp[ii] = 0.5 * ((*data)->r_visc[ii] + (*data)->r_visc[gmap->icjckP[ii]]);
        }
        // surface-subsurface interface
        if (gmap->istop[ii] == 1)
        {
            (*data)->r_rhozp[gmap->icjckM[ii]] = 0.5 * ((*data)->r_rho[ii] + (*data)->r_rho[gmap->icjckM[ii]]);
            (*data)->r_visczp[gmap->icjckM[ii]] = 0.5 * ((*data)->r_visc[ii] + (*data)->r_visc[gmap->icjckM[ii]]);
        }
    }
    // conductivities on iM, jM, kM faces
    for (ii = 0; ii < param->ny*param->nz; ii++)
    {
        if (gmap->actv[gmap->iMin[ii]] == 1)
        {
            (*data)->r_rhoxp[gmap->iMou[ii]] = 0.5 * ((*data)->r_rho[gmap->iMin[ii]] + (*data)->r_rho[gmap->iMou[ii]]);
            (*data)->r_viscxp[gmap->iMou[ii]] = 0.5 * ((*data)->r_visc[gmap->iMin[ii]] + (*data)->r_visc[gmap->iMou[ii]]);
        }
    }
    for (ii = 0; ii < param->nx*param->nz; ii++)
    {
        if (gmap->actv[gmap->jMin[ii]] == 1)
        {
            (*data)->r_rhoyp[gmap->jMou[ii]] = 0.5 * ((*data)->r_rho[gmap->jMin[ii]] + (*data)->r_rho[gmap->jMou[ii]]);
            (*data)->r_viscyp[gmap->jMou[ii]] = 0.5 * ((*data)->r_visc[gmap->jMin[ii]] + (*data)->r_visc[gmap->jMou[ii]]);
            if (irank < param->mpi_nx & param->bctype_GW[3] == 1)
            {
                (*data)->r_rhoyp[gmap->jMou[ii]] = (*data)->r_rho[gmap->jMou[ii]];
                (*data)->r_viscyp[gmap->jMou[ii]] = (*data)->r_visc[gmap->jMou[ii]];
            }
        }
    }
}

// >>>>> Compute matrix coefficients <<<<<
void groundwater_mat_coeff(Data **data, Map *gmap, Config *param, int irank)
{
    int ii, jp, jm;
    double dzf;
    for (ii = 0; ii < param->n3ci; ii++)
    {
        // coeff xp
        (*data)->Gxp[ii] = - (*data)->Kx[ii] * param->dt * (*data)->r_rhoxp[ii] * (*data)->r_viscxp[ii] / (pow(param->dx,2.0));
        if (param->sim_shallowwater == 1)
        {if (gmap->actvXp[ii] == 1)  {(*data)->Gxp[ii] = 2.0 * (*data)->Gxp[ii];}}
        // coeff xm
        (*data)->Gxm[ii] = - (*data)->Kx[gmap->iMjckc[ii]] * param->dt * (*data)->r_rhoxp[gmap->iMjckc[ii]] * (*data)->r_viscxp[gmap->iMjckc[ii]] / (pow(param->dx,2.0));
        if (param->sim_shallowwater == 1)
        {if (gmap->actvXp[gmap->iMjckc[ii]] == 1)  {(*data)->Gxm[ii] = 2.0 * (*data)->Gxm[ii];}}
        // coeff yp
        (*data)->Gyp[ii] = - (*data)->Ky[ii] * param->dt * (*data)->r_rhoyp[ii] * (*data)->r_viscyp[ii] / (pow(param->dy,2.0));
        if (gmap->jj[ii] == param->ny-1 & param->bctype_GW[2] == 1 & irank >= param->mpi_nx*(param->mpi_ny-1))     {(*data)->Gyp[ii] = 2.0 * (*data)->Gyp[ii]; }
        if (param->sim_shallowwater == 1)
        {if (gmap->actvYp[ii] == 1)  {(*data)->Gyp[ii] = 2.0 * (*data)->Gyp[ii];}}
        // coeff ym
        (*data)->Gym[ii] = - (*data)->Ky[gmap->icjMkc[ii]] * param->dt * (*data)->r_rhoyp[gmap->icjMkc[ii]] * (*data)->r_viscyp[gmap->icjMkc[ii]] / (pow(param->dy,2.0));
        if (gmap->jj[ii] == 0 & param->bctype_GW[3] == 1 & irank < param->mpi_nx)     {(*data)->Gym[ii] = 2.0 * (*data)->Gym[ii];}
        if (param->sim_shallowwater == 1)
        {if (gmap->actvYp[gmap->icjMkc[ii]] == 1)  {(*data)->Gym[ii] = 2.0 * (*data)->Gym[ii];}}
        // coeff zp
        dzf = 0.5 * (gmap->dz3d[ii] + gmap->dz3d[gmap->icjckP[ii]]);
        (*data)->Gzp[ii] = - (*data)->Kz[ii] * param->dt * (*data)->r_rhozp[ii] * (*data)->r_visczp[ii] / (gmap->dz3d[ii]*dzf);
        // coeff zm
        if (gmap->istop[ii] == 1)   {dzf = 0.5 * gmap->dz3d[ii];}
        else    {dzf = 0.5 * (gmap->dz3d[ii] + gmap->dz3d[gmap->icjckM[ii]]);}
        (*data)->Gzm[ii] = - (*data)->Kz[gmap->icjckM[ii]] * param->dt * (*data)->r_rhozp[gmap->icjckM[ii]] * (*data)->r_visczp[gmap->icjckM[ii]] / (gmap->dz3d[ii]*dzf);
        // coeff ct
        (*data)->Gct[ii] = ((*data)->ch[ii] + param->Ss*(*data)->wcn[ii]/(*data)->wcs[ii]) * (*data)->r_rho[ii];
        (*data)->Gct[ii] -= ((*data)->Gxp[ii] + (*data)->Gxm[ii] + (*data)->Gyp[ii] + (*data)->Gym[ii]);
        // ct on yp
        if (gmap->jj[ii] == param->ny-1)
        {if (param->bctype_GW[2] == 2 & irank >= param->mpi_nx*(param->mpi_ny-1))    {(*data)->Gct[ii] += (*data)->Gyp[ii];}}
        // ct on ym
        if (gmap->jj[ii] == 0)
        {if (param->bctype_GW[3] == 2 & irank < param->mpi_nx)    {(*data)->Gct[ii] += (*data)->Gym[ii];}}
        // ct on zp
        if (gmap->kk[ii] == param->nz-1 & param->bctype_GW[4] != 1)
        {(*data)->Gct[ii] = (*data)->Gct[ii];}
        else
        {(*data)->Gct[ii] -= (*data)->Gzp[ii];}
        // ct on zm
        if (gmap->istop[ii] == 1)
        {
            if (param->sim_shallowwater == 1 & (*data)->dept[gmap->top2d[ii]] > 0)
            {(*data)->Gct[ii] -= (*data)->Gzm[ii];}
            else if (param->bctype_GW[5] == 1)
            {(*data)->Gct[ii] -= (*data)->Gzm[ii];}
            // printf(" Gct = %f, Gzp = %f\n",(*data)->Gct[ii],(*data)->Gzp[ii]);
        }
        else
        {(*data)->Gct[ii] -= (*data)->Gzm[ii];}
    }
}

// >>>>> Compute right hand side <<<<<
void groundwater_rhs(Data **data, Map *gmap, Config *param, int irank)
{
    int ii;
    double dh;
    for (ii = 0; ii < param->n3ci; ii++)
    {
        // base terms
        (*data)->Grhs[ii] = ((*data)->ch[ii] + param->Ss*(*data)->wcn[ii]/(*data)->wcs[ii]) * (*data)->hn[ii] * (*data)->r_rho[ii];
        (*data)->Grhs[ii] -= param->dt * ((*data)->Kz[ii] * (*data)->r_rhozp[ii] * (*data)->r_rhozp[ii] * (*data)->r_visczp[ii] \
            - (*data)->Kz[gmap->icjckM[ii]] * (*data)->r_rhozp[gmap->icjckM[ii]] * (*data)->r_rhozp[gmap->icjckM[ii]] * (*data)->r_visczp[gmap->icjckM[ii]]) / gmap->dz3d[ii];
        // density term
        // (*data)->Grhs[ii] -= (*data)->wc[ii] * ((*data)->r_rho[ii] - (*data)->r_rhon[ii]);
        // boundary terms
        // x
        if (gmap->ii[ii] == param->nx-1)
        {(*data)->Grhs[ii] -= (*data)->Gxp[ii] * (*data)->hn[gmap->iPjckc[ii]];}
        if (gmap->ii[ii] == 0)
        {(*data)->Grhs[ii] -= (*data)->Gxm[ii] * (*data)->hn[gmap->iMjckc[ii]];}
        if (param->sim_shallowwater == 1)
        {
            // lateral surface-subsurface exchange
            if (gmap->actvXp[ii] == 1)
            {
                dh = (*data)->eta[gmap->top2d[gmap->iPjckc[ii]]] - (*data)->offset[0] - (gmap->bot3d[ii] + 0.5*gmap->dz3d[ii]);
                (*data)->Grhs[ii] -= (*data)->Gxp[ii] * dh;
            }
            if (gmap->actvXp[gmap->iMjckc[ii]] == 1)
            {
                dh = (*data)->eta[gmap->top2d[gmap->iMjckc[ii]]] - (*data)->offset[0] - (gmap->bot3d[ii] + 0.5*gmap->dz3d[ii]);
                (*data)->Grhs[ii] -= (*data)->Gxm[ii] * dh;
            }
        }
        // y
        if (gmap->jj[ii] == (param->ny-1))
        {
            // flux bc
            if (param->bctype_GW[2] == 2 & irank >= param->mpi_nx*(param->mpi_ny-1))
            {(*data)->Grhs[ii] -= (*data)->r_rhoyp[ii] * param->qyp * param->dt / param->dy;}
            else
            {(*data)->Grhs[ii] -= (*data)->Gyp[ii] * (*data)->hn[gmap->icjPkc[ii]];}
        }
        if (gmap->jj[ii] == 0)
        {

            if (param->bctype_GW[3] == 2 & irank < param->mpi_nx)
            {
                (*data)->Grhs[ii] -= (*data)->r_rhoyp[gmap->icjMkc[ii]] * param->qym * param->dt / param->dy;
            }
            else
            {(*data)->Grhs[ii] -= (*data)->Gym[ii] * (*data)->hn[gmap->icjMkc[ii]];}

        }
        if (param->sim_shallowwater == 1)
        {
            // lateral surface-subsurface exchange
            if (gmap->actvYp[ii] == 1)
            {
                dh = (*data)->eta[gmap->top2d[gmap->icjPkc[ii]]] - (*data)->offset[0] - (gmap->bot3d[ii] + 0.5*gmap->dz3d[ii]);
                (*data)->Grhs[ii] -= (*data)->Gyp[ii] * dh;
            }
            if (gmap->actvYp[gmap->icjMkc[ii]] == 1)
            {
                dh = (*data)->eta[gmap->top2d[gmap->icjMkc[ii]]] - (*data)->offset[0] - (gmap->bot3d[ii] + 0.5*gmap->dz3d[ii]);
                (*data)->Grhs[ii] -= (*data)->Gym[ii] * dh;
            }
        }
        // z
        if (gmap->kk[ii] == param->nz-1)
        {
            if (param->bctype_GW[4] == 2)
            {(*data)->Grhs[ii] += (*data)->r_rhozp[ii] * ((*data)->qbot + param->dt * (*data)->Kz[ii] * (*data)->r_rhozp[ii] * (*data)->r_visczp[ii]/ gmap->dz3d[ii]);}
            else if (param->bctype_GW[4] == 1)
            {(*data)->Grhs[ii] += (*data)->Gzp[ii] * (*data)->hn[gmap->icjckP[ii]];}
        }
        if (gmap->istop[ii] == 1)
        {
            if (param->sim_shallowwater == 1)
            {
                if ((*data)->dept[gmap->top2d[ii]] > 0)
                {(*data)->Grhs[ii] -= (*data)->Gzm[ii] * (*data)->dept[gmap->top2d[ii]];}
                else
                {
                    if (param->bctype_GW[5] == 2)
                    {
                        if ((*data)->qtop[gmap->top2d[ii]] > 0.0 & (*data)->wc[ii] < (*data)->wcs[ii])
                        {
                            (*data)->Grhs[ii] += - param->dt * (*data)->r_rhozp[gmap->icjckM[ii]] * \
                                ((*data)->qtop[gmap->top2d[ii]] + (*data)->Kz[gmap->icjckM[ii]]*(*data)->r_rhozp[gmap->icjckM[ii]]*(*data)->r_visczp[gmap->icjckM[ii]]) / gmap->dz3d[ii];
                        }
                        else if ((*data)->qtop[gmap->top2d[ii]] < 0.0 & (*data)->wc[ii] > (*data)->wcr[ii])
                        {
                            (*data)->Grhs[ii] += - param->dt * (*data)->r_rhozp[gmap->icjckM[ii]] * \
                                ((*data)->qtop[gmap->top2d[ii]] + (*data)->Kz[gmap->icjckM[ii]]*(*data)->r_rhozp[gmap->icjckM[ii]]*(*data)->r_visczp[gmap->icjckM[ii]]) / gmap->dz3d[ii];
                        }
                        else
                        {
                            (*data)->Grhs[ii] += - param->dt * (*data)->r_rhozp[gmap->icjckM[ii]] * \
                                (*data)->Kz[gmap->icjckM[ii]]*(*data)->r_rhozp[gmap->icjckM[ii]]*(*data)->r_visczp[gmap->icjckM[ii]] / gmap->dz3d[ii];
                        }
                    }
                    else if (param->bctype_GW[5] == 1)
                    {(*data)->Grhs[ii] -= (*data)->Gzm[ii] * (*data)->hn[gmap->icjckM[ii]];}
                }
            }
            else
            {
                if (param->bctype_GW[5] == 2)
                {
                    // if ((*data)->qtop > 0.0 & (*data)->wc[ii] < (*data)->wcs[ii])
                    // {(*data)->Grhs[ii] += - param->dt * ((*data)->qtop + (*data)->Kz[gmap->icjckM[ii]]) / gmap->dz3d[ii];}
                    // else if ((*data)->qtop < 0.0 & (*data)->wc[ii] > (*data)->wcr[ii])
                    // {(*data)->Grhs[ii] += - param->dt * ((*data)->qtop + (*data)->Kz[gmap->icjckM[ii]]) / gmap->dz3d[ii];}
                    // else
                    // {(*data)->Grhs[ii] += - param->dt * (*data)->Kz[gmap->icjckM[ii]] / gmap->dz3d[ii];}
                    (*data)->Grhs[ii] += - param->dt * (*data)->r_rhozp[gmap->icjckM[ii]] * \
                        ((*data)->qtop[gmap->top2d[ii]] + (*data)->Kz[gmap->icjckM[ii]] * (*data)->r_rhozp[gmap->icjckM[ii]] * (*data)->r_visczp[gmap->icjckM[ii]]) / gmap->dz3d[ii];
                    // printf("  >> Top flux = %f\n",86400.0*100.0*(*data)->qtop);
                }
                else if (param->bctype_GW[5] == 1)
                {(*data)->Grhs[ii] -= (*data)->Gzm[ii] * (*data)->htop;}
            }

            // printf(" RHS = %f, term1=%f, term2=%f, term3=%f\n",(*data)->Grhs[ii],((*data)->ch[ii] + param->Ss*(*data)->wcn[ii]/(*data)->wcs[ii]) * (*data)->hn[ii] * (*data)->r_rho[ii],
            // param->dt * ((*data)->Kz[ii] * (*data)->r_rhozp[ii] * (*data)->r_rhozp[ii] * (*data)->r_visczp[ii] \
            //     - (*data)->Kz[gmap->icjckM[ii]] * (*data)->r_rhozp[gmap->icjckM[ii]] * (*data)->r_rhozp[gmap->icjckM[ii]] * (*data)->r_visczp[gmap->icjckM[ii]]) / gmap->dz3d[ii],
            //     -(*data)->Gzm[ii] * (*data)->htop);
        }
    }
}

// >>>>> Build linear system <<<<<
void build_groundwater_system(Data *data, Map *gmap, Config *param, QMatrix A, Vector b)
{
    size_t ii, jj, kk, n, dist, irow;
    int im;
    double value, *row;
    row = malloc(8*sizeof(double));

    irow = 1;
    for (ii = 1; ii <= param->n3ci; ii++)
    {

        for (jj = 0; jj < 8; jj++)  {row[jj] = 0.0;}
        im = ii - 1;
        if (gmap->actv[im] == 0)
        {
            Q_SetLen(&A, ii, 1);
            Q_SetEntry(&A, ii, 0, ii, 1.0);
            V_SetCmp(&b, ii, data->hn[im]);
            row[3] = 1.0;   row[7] = data->hn[im];
        }
        else
        {
            kk = 0;
            n = 7;
            if (gmap->ii[im] == 0)  {n -= 1;}
            if (gmap->ii[im] == param->nx-1)    {n -= 1;}
            if (gmap->jj[im] == 0)  {n -= 1;}
            if (gmap->jj[im] == param->ny-1)    {n -= 1;}
            if (gmap->istop[im] == 1)   {n -= 1;}
            if (gmap->kk[im] == param->nz-1)    {n -= 1;}
            Q_SetLen(&A, ii, n);
            // set ym entry
            dist = im - gmap->icjMkc[im];
            if (ii-dist > 0 & ii-dist <= param->n3ci & gmap->actv[gmap->icjMkc[im]] == 1)
            {Q_SetEntry(&A, ii, kk, ii-dist, data->Gym[im]);    kk++;   row[0]=data->Gym[im];}
            // set xm entry
            dist = im - gmap->iMjckc[im];
            if (ii-dist > 0 & ii-dist <= param->n3ci & gmap->actv[gmap->iMjckc[im]] == 1)
            {Q_SetEntry(&A, ii, kk, ii-dist, data->Gxm[im]);    kk++;   row[1]=data->Gxm[im];}
            // set zm entry
            dist = im - gmap->icjckM[im];
            if (ii-dist > 0 & ii-dist <= param->n3ci & gmap->actv[gmap->icjckM[im]] == 1)
            {Q_SetEntry(&A, ii, kk, ii-dist, data->Gzm[im]);    kk++;   row[2]=data->Gzm[im];}
            // set ct entry
            Q_SetEntry(&A, ii, kk, ii, data->Gct[im]);
            row[3]=data->Gct[im];
            kk++;
            // set zp entry
            dist = gmap->icjckP[im] - im;
            if (ii+dist > 0 & ii+dist <= param->n3ci & gmap->actv[gmap->icjckP[im]] == 1)
            {Q_SetEntry(&A, ii, kk, ii+dist, data->Gzp[im]);    kk++;   row[4]=data->Gzp[im];}
            // set xp entry
            dist = gmap->iPjckc[im] - im;
            if (ii+dist > 0 & ii+dist <= param->n3ci & gmap->actv[gmap->iPjckc[im]] == 1)
            {Q_SetEntry(&A, ii, kk, ii+dist, data->Gxp[im]);    kk++;   row[5]=data->Gxp[im];}
            // set yp entry
            dist = gmap->icjPkc[im] - im;
            if (ii+dist > 0 & ii+dist <= param->n3ci & gmap->actv[gmap->icjPkc[im]] == 1)
            {Q_SetEntry(&A, ii, kk, ii+dist, data->Gyp[im]);    kk++;   row[6]=data->Gyp[im];}
            // set right hand side
            V_SetCmp(&b, ii, data->Grhs[im]);
            row[7]=data->Grhs[im];

        }

        // if (gmap->ii[im] == 0 & gmap->istop[im] == 1)
        // {
        //     printf("(%d,%d,%d) - [%f, %f, %f, (%f), %f, %f, %f][h] = [%f]\n",gmap->ii[im],gmap->jj[im],gmap->kk[im],
        //         1e5*row[0],1e5*row[1],1e5*row[2],1e5*row[3],1e5*row[4],1e5*row[5],1e5*row[6],1e5*row[7],1e5*row[8]);
        //
        // }

    }
}

// >>>>> solve linear system <<<<<
void solve_groundwater_system(Data **data, Map *gmap, QMatrix A, Vector b, Vector x, Config *param)
{
    size_t ii;
    V_SetAllCmp(&x, 0.0);
    SetRTCAccuracy(0.00000001);
    CGIter(&A, &x, &b, 10000000, SSORPrecond, 1);
    for (ii = 0; ii < param->n3ci; ii++)    {(*data)->h[ii] = V_GetCmp(&x, ii+1);}
}

// >>>>> enforce head boundary conditions
void enforce_head_bc(Data **data, Map *gmap, Config *param)
{
    int ii;
    // xp and xm boundaries
    for (ii = 0; ii < param->ny*param->nz; ii++)
    {
        (*data)->h[gmap->iPou[ii]] = (*data)->h[gmap->iPin[ii]];
        (*data)->h[gmap->iMou[ii]] = (*data)->h[gmap->iMin[ii]];
    }
    // yp and ym boundaries
    for (ii = 0; ii < param->nx*param->nz; ii++)
    {
        (*data)->h[gmap->jPou[ii]] = (*data)->h[gmap->jPin[ii]];
        (*data)->h[gmap->jMou[ii]] = (*data)->h[gmap->jMin[ii]];
    }
    // enforce hydrostatic head bc for side boundaries
    if (param->bctype_GW[2] == 1)
    {
        for (ii = 0; ii < param->nx*param->nz; ii++)
        {
            (*data)->h[gmap->jPou[ii]] = (*data)->bottom[gmap->top2d[gmap->jPin[ii]]] - gmap->bot3d[gmap->jPin[ii]] - 0.5*gmap->dz3d[gmap->jPin[ii]];
            (*data)->h[gmap->jPou[ii]] = (*data)->h[gmap->jPou[ii]] * (*data)->r_rhoyp[gmap->jPin[ii]];
            if (param->sim_shallowwater == 1)
            {(*data)->h[gmap->jPou[ii]] += ((*data)->dept[gmap->top2d[gmap->jPin[ii]]] * (*data)->r_rhoyp[gmap->jPin[ii]]);}
        }
    }
    if (param->bctype_GW[3] == 1)
    {
        for (ii = 0; ii < param->nx*param->nz; ii++)
        {
            (*data)->h[gmap->jMou[ii]] = (*data)->bottom[gmap->top2d[gmap->jMin[ii]]] - gmap->bot3d[gmap->jMin[ii]] - 0.5*gmap->dz3d[gmap->jMin[ii]];
            (*data)->h[gmap->jMou[ii]] = (*data)->h[gmap->jMou[ii]] * (*data)->r_rhoyp[gmap->jMou[ii]];
            if (param->sim_shallowwater == 1)    {(*data)->h[gmap->jMou[ii]] += (*data)->dept[gmap->top2d[gmap->jMin[ii]]];}
            // Kuan2019
            // (*data)->h[gmap->jMou[ii]] -= 0.24;
        }
    }
    // zp boundary
    if (param->bctype_GW[4] == 1)
    {for (ii = 0; ii < param->nx*param->ny; ii++)    {(*data)->h[gmap->kPou[ii]] = (*data)->hbot;}}
    else
    {for (ii = 0; ii < param->nx*param->ny; ii++)    {(*data)->h[gmap->kPou[ii]] = (*data)->h[gmap->kMin[ii]];}}
    // zm boundary
    if (param->sim_shallowwater == 1)
    {
        for (ii = 0; ii < param->n3ci; ii++)
        {if (gmap->istop[ii] == 1)   {(*data)->h[gmap->icjckM[ii]] = (*data)->dept[gmap->top2d[ii]];}}
    }
    else
    {
        if (param->bctype_GW[5] == 1)
        {
            for (ii = 0; ii < param->n3ci; ii++)
            {if (gmap->istop[ii] == 1)  {(*data)->h[gmap->icjckM[ii]] = (*data)->htop;}}
        }
        else
        {
            for (ii = 0; ii < param->n3ci; ii++)
            {if (gmap->istop[ii] == 1)  {(*data)->h[gmap->icjckM[ii]] = (*data)->h[ii];}}
        }
    }
}

// >>>>> calculate flux between cells <<<<<
void groundwater_flux(Data **data, Map *gmap, Config *param, int irank)
{
    int ii, jj;
    double dzf, vseep, dh;
    for (ii = 0; ii < param->n3ci; ii++)
    {
        dzf = 0.5 * (gmap->dz3d[ii] + gmap->dz3d[gmap->icjckP[ii]]);
        // x flux
        (*data)->qx[ii] = (*data)->Kx[ii] * (*data)->r_viscxp[ii] * ((*data)->h[gmap->iPjckc[ii]]-(*data)->h[ii]) / param->dx;
        if (gmap->actvXp[ii] == 1)
        {
            if (gmap->actv[ii] == 1)
            {
                dh = (*data)->eta[gmap->top2d[gmap->iPjckc[ii]]] - (*data)->offset[0] - (gmap->bot3d[ii] + 0.5*gmap->dz3d[ii]);
                (*data)->qx[ii] = (*data)->Kx[ii] * (*data)->r_viscxp[ii] * (dh-(*data)->h[ii]) / param->dx;
            }
            else
            {
                dh = (*data)->eta[gmap->top2d[ii]] - (*data)->offset[0] - (gmap->bot3d[gmap->iPjckc[ii]] + 0.5*gmap->dz3d[gmap->iPjckc[ii]]);
                (*data)->qx[ii] = (*data)->Kx[ii] * (*data)->r_viscxp[ii] * ((*data)->h[gmap->iPjckc[ii]]-dh) / param->dx;
            }
        }
        // y flux
        if (gmap->jj[ii] == param->ny-1 & param->bctype_GW[2] == 2 & irank >= param->mpi_nx*(param->mpi_ny-1))
        {(*data)->qy[ii] = param->qyp;}
        else if (gmap->actvYp[ii] == 1)
        {
            if (gmap->actv[ii] == 1)
            {
                dh = (*data)->eta[gmap->top2d[gmap->icjPkc[ii]]] - (*data)->offset[0] - (gmap->bot3d[ii] + 0.5*gmap->dz3d[ii]);
                (*data)->qy[ii] = (*data)->Ky[ii] * (*data)->r_viscyp[ii] * (dh-(*data)->h[ii]) / param->dy;
            }
            else
            {
                dh = (*data)->eta[gmap->top2d[ii]] - (*data)->offset[0] - (gmap->bot3d[gmap->icjPkc[ii]] + 0.5*gmap->dz3d[gmap->icjPkc[ii]]);
                (*data)->qy[ii] = (*data)->Ky[ii] * (*data)->r_viscyp[ii] * ((*data)->h[gmap->icjPkc[ii]]-dh) / param->dy;
            }
        }
        else
        {
            (*data)->qy[ii] = (*data)->Ky[ii] * (*data)->r_viscyp[ii] * ((*data)->h[gmap->icjPkc[ii]]-(*data)->h[ii]) / param->dy;
            if (gmap->jj[ii] == param->ny-1 & param->bctype_GW[2] == 1 & irank >= param->mpi_nx*(param->mpi_ny-1))  {(*data)->qy[ii] = 2.0 * (*data)->qy[ii]; }
        }
        // z flux
        if (gmap->kk[ii] == param->nz-1)
        {
            if (param->bctype_GW[4] == 2)   {(*data)->qz[ii] = (*data)->qbot;}
            else if (param->bctype_GW[4] == 3)   {(*data)->qz[ii] = - (*data)->Kz[ii] * (*data)->r_visczp[ii] * (*data)->r_rhozp[ii];}
            else
            {
                (*data)->qz[ii] = (*data)->Kz[ii] * (*data)->r_visczp[ii] * \
                    (((*data)->h[gmap->icjckP[ii]]-(*data)->h[ii]) / dzf - (*data)->r_rhozp[ii]);
            }
        }
        else
        {
            (*data)->qz[ii] = (*data)->Kz[ii] * (*data)->r_visczp[ii] * \
                (((*data)->h[gmap->icjckP[ii]]-(*data)->h[ii]) / dzf - (*data)->r_rhozp[ii]);
        }
    }
    for (ii = 0; ii < param->ny*param->nz; ii++)
    {
        (*data)->qx[gmap->iMou[ii]] = (*data)->Kx[gmap->iMou[ii]] * (*data)->r_viscxp[ii] \
                    * ((*data)->h[gmap->iMin[ii]]-(*data)->h[gmap->iMou[ii]]) / param->dx;
    }
    for (ii = 0; ii < param->nx*param->nz; ii++)
    {
        if (param->bctype_GW[3] == 2 & irank < param->mpi_nx)
        {(*data)->qy[gmap->jMou[ii]] = param->qym;}
        else
        {
            (*data)->qy[gmap->jMou[ii]] = (*data)->Ky[gmap->jMou[ii]] * (*data)->r_viscyp[ii] \
                        * ((*data)->h[gmap->jMin[ii]]-(*data)->h[gmap->jMou[ii]]) / param->dy;
            if (param->bctype_GW[3] == 1 & irank < param->mpi_nx)   {(*data)->qy[gmap->jMou[ii]] = 2.0 * (*data)->qy[gmap->jMou[ii]];}
        }
    }
    for (jj = 0; jj < param->n3ci; jj++)
    {
        if (gmap->istop[jj] == 1)
        {
            ii = gmap->icjckM[jj];
            dzf = 0.5 * gmap->dz3d[jj];
            if (param->sim_shallowwater == 1)
            {
                if ((*data)->dept[gmap->top2d[jj]] > 0)
                {
                    (*data)->qz[ii] = (*data)->Kz[ii] * (*data)->r_visczp[ii] * (((*data)->h[jj]-(*data)->h[ii]) / dzf - (*data)->r_rhozp[ii]);
                    if ((*data)->qz[ii] < 0)
                    {
                        // check if infiltration > depth available
                        // vseep = fabs((*data)->qz[ii])*param->dt;
                        // if (vseep > (*data)->dept[gmap->top2d[jj]])
                        // {(*data)->qz[ii] = -(*data)->dept[gmap->top2d[jj]]/param->dt;}
                        vseep = fabs((*data)->qz[ii])*param->dt*(*data)->wcs[jj];
                        if (vseep > (*data)->dept[gmap->top2d[jj]])
                        {(*data)->qz[ii] = -(*data)->dept[gmap->top2d[jj]]/(param->dt*(*data)->wcs[jj]);}
                        // check if infiltration > void volume in the first subsurface cell
                        // else if (vseep/(*data)->wcs[jj] > gmap->dz3d[jj]*((*data)->wcs[jj] - (*data)->wc[jj]))
                        // {(*data)->qz[ii] = -gmap->dz3d[jj]*((*data)->wcs[jj] - (*data)->wc[jj])/param->dt;}
                    }
                }
                else if (param->bctype_GW[5] == 2)
                {
                    (*data)->qz[ii] = (*data)->Kz[ii] * (*data)->r_visczp[ii] * (((*data)->h[jj]-(*data)->h[ii]) / dzf - (*data)->r_rhozp[ii]);
                    // infiltration over dry surface is not realistic for coupled surface-subsurface simulations
                    // i.e. infiltration must happen through surface-subsurface exchange, it cannot be directly
                    //      applied to the subsurface domain, ZhiLi20201116
                    if ((*data)->qz[ii] < 0.0)  {(*data)->qz[ii] = 0.0;}
                    // add evaporation
                    if ((*data)->qtop[gmap->top2d[jj]] > 0.0 & (*data)->wc[jj] > (*data)->wcr[jj])
                    {
                        if ((*data)->wc[jj] > (*data)->wcr[jj] + (*data)->qtop[gmap->top2d[jj]] * param->dt / gmap->dz3d[jj])
                        {(*data)->qz[ii] += (*data)->qtop[gmap->top2d[jj]];}
                    }
                    // add rainfall
                    else if ((*data)->qtop[gmap->top2d[jj]] < 0.0 & (*data)->wc[jj] < (*data)->wcs[jj])
                    {(*data)->qz[ii] += (*data)->qtop[gmap->top2d[jj]];}
                }
                else if (param->bctype_GW[5] == 3)
                {(*data)->qz[ii] = -(*data)->Kz[ii] * (*data)->r_visczp[ii] * (*data)->r_rhozp[ii];}
                else if (param->bctype_GW[5] == 1)
                {(*data)->qz[ii] = (*data)->Kz[ii] * (*data)->r_visczp[ii] * (((*data)->h[jj]-(*data)->h[ii]) / dzf - (*data)->r_rhozp[ii]);}
                else
                {(*data)->qz[ii] = 0.0;}
                // surface-subsurface exchange
                if ((*data)->reset_seepage[gmap->top2d[jj]] == 1)
                {
                    (*data)->qseepage[gmap->top2d[jj]] = 0.0;
                    (*data)->reset_seepage[gmap->top2d[jj]] = 0;
                }
                (*data)->qseepage[gmap->top2d[jj]] += (*data)->qz[ii];
                // if evaporation exists, evaporation does not contribute to seepage
                if (param->bctype_GW[5] == 2 & (*data)->qtop[gmap->top2d[jj]] > 0.0 & (*data)->wc[jj] > (*data)->wcr[jj])
                {(*data)->qseepage[gmap->top2d[jj]] -= (*data)->qtop[gmap->top2d[jj]];}

            }
            else
            {
                if (param->bctype_GW[5] == 1)
                {(*data)->qz[ii] = (*data)->Kz[ii] * (*data)->r_visczp[ii] * (((*data)->h[jj]-(*data)->h[ii]) / dzf - (*data)->r_rhozp[ii]);}
                else if (param->bctype_GW[5] == 2)
                {
                    if ((*data)->qtop[gmap->top2d[jj]] < 0.0 & (*data)->wc[ii] < (*data)->wcs[ii])
                    {(*data)->qz[ii] = (*data)->qtop[gmap->top2d[jj]];}
                    else if ((*data)->qtop[gmap->top2d[jj]] > 0.0 & (*data)->wc[ii] > (*data)->wcr[ii])
                    {(*data)->qz[ii] = (*data)->qtop[gmap->top2d[jj]];}
                    else
                    {(*data)->qz[ii] = 0.0;}
                }
                else if (param->bctype_GW[5] == 3)
                {(*data)->qz[ii] = -(*data)->Kz[ii] * (*data)->r_visczp[ii] * (*data)->r_rhozp[ii];}
                else
                {(*data)->qz[ii] = 0.0;}

                // Geng 2015
                // if (gmap->jj[jj] < 20)    {(*data)->qz[ii] = 0.0;}

            }
        }
    }

}

// check available room in each grid cell
void check_room(Data **data, Map *gmap, Config *param)
{
    int ii;
    for (ii = 0; ii < param->n3ci; ii++)
    {
        if (gmap->actv[ii] == 0)
        {(*data)->room[ii] = 0.0;}
        else
        {(*data)->room[ii] = ((*data)->wcs[ii]-(*data)->wc[ii]) * gmap->dz3d[ii]*param->dx*param->dy;}
    }
}

// >>>>> update water content <<<<<
void update_water_content(Data **data, Map *gmap, Config *param)
{
    int ii, jj, i_sat;
    double coeff, dqx, dqy, dqz, room_col, dz_eff;
    for (ii = 0; ii < param->n3ci; ii++)
    {
        // update water content
        coeff = (*data)->r_rho[ii] + (*data)->r_rho[ii] * param->Ss * ((*data)->h[ii]-(*data)->hn[ii]) / (*data)->wcs[ii];
        // x component
        if (param->sim_shallowwater == 1)
        {
            if (gmap->actvXp[ii] == 1 & gmap->actv[ii] == 1)
            {
                dz_eff = (*data)->eta[gmap->top2d[gmap->iPjckc[ii]]] - (*data)->offset[0] - gmap->bot3d[ii];
                dqx = param->dt * ((*data)->qx[ii]*(*data)->r_rhoxp[ii]*dz_eff/gmap->dz3d[ii]
                    - (*data)->qx[gmap->iMjckc[ii]]*(*data)->r_rhoxp[gmap->iMjckc[ii]]) / param->dx;
            }
            else if (gmap->actvXp[gmap->iMjckc[ii]] == 1 & gmap->actv[ii] == 1)
            {
                dz_eff = (*data)->eta[gmap->top2d[gmap->iMjckc[ii]]] - (*data)->offset[0] - gmap->bot3d[ii];
                dqx = param->dt * ((*data)->qx[ii]*(*data)->r_rhoxp[ii]
                    - (*data)->qx[gmap->iMjckc[ii]]*(*data)->r_rhoxp[gmap->iMjckc[ii]]*dz_eff/gmap->dz3d[ii]) / param->dx;
            }
            else
            {dqx = param->dt * ((*data)->qx[ii]*(*data)->r_rhoxp[ii] - (*data)->qx[gmap->iMjckc[ii]]*(*data)->r_rhoxp[gmap->iMjckc[ii]]) / param->dx;}
        }
        else
        {dqx = param->dt * ((*data)->qx[ii]*(*data)->r_rhoxp[ii] - (*data)->qx[gmap->iMjckc[ii]]*(*data)->r_rhoxp[gmap->iMjckc[ii]]) / param->dx;}

        // y component
        if (param->sim_shallowwater == 1)
        {
            if (gmap->actvYp[ii] == 1 & gmap->actv[ii] == 1)
            {
                dz_eff = (*data)->eta[gmap->top2d[gmap->icjPkc[ii]]] - (*data)->offset[0] - gmap->bot3d[ii];
                dqy = param->dt * ((*data)->qy[ii]*(*data)->r_rhoyp[ii]*dz_eff/gmap->dz3d[ii]
                    - (*data)->qy[gmap->icjMkc[ii]]*(*data)->r_rhoyp[gmap->icjMkc[ii]]) / param->dy;
            }
            else if (gmap->actvYp[gmap->icjMkc[ii]] == 1 & gmap->actv[ii] == 1)
            {
                dz_eff = (*data)->eta[gmap->top2d[gmap->icjMkc[ii]]] - (*data)->offset[0] - gmap->bot3d[ii];
                dqy = param->dt * ((*data)->qy[ii]*(*data)->r_rhoyp[ii]
                    - (*data)->qy[gmap->icjMkc[ii]]*(*data)->r_rhoyp[gmap->icjMkc[ii]]*dz_eff/gmap->dz3d[ii]) / param->dy;
            }
            else
            {dqy = param->dt * ((*data)->qy[ii]*(*data)->r_rhoyp[ii] - (*data)->qy[gmap->icjMkc[ii]]*(*data)->r_rhoyp[gmap->icjMkc[ii]]) / param->dy;}
        }
        else
        {dqy = param->dt * ((*data)->qy[ii]*(*data)->r_rhoyp[ii] - (*data)->qy[gmap->icjMkc[ii]]*(*data)->r_rhoyp[gmap->icjMkc[ii]]) / param->dy;}
        // z component
        dqz = param->dt * ((*data)->qz[ii]*(*data)->r_rhozp[ii] - (*data)->qz[gmap->icjckM[ii]]*(*data)->r_rhozp[gmap->icjckM[ii]]) / gmap->dz3d[ii];

        (*data)->wc[ii] = ((*data)->wcn[ii]*(*data)->r_rho[ii] + dqx + dqy + dqz) / coeff;

        // if (gmap->jj[ii] == 140 & gmap->istop[ii] == 1) {
        //     printf(" (jj,kk)=(%d,%d)  wc=%f, dept=%f, actv=(%d,%d,%d) \n",gmap->jj[ii],gmap->kk[ii],
        //         (*data)->wc[ii],(*data)->dept[gmap->top2d[ii]],gmap->actv[gmap->icjMkc[ii]],gmap->actv[ii],gmap->actv[gmap->icjPkc[ii]]);
        //     printf("                  Z - qz=%f, head=%f, dept=%f, deptYP=%f\n",
        //         (*data)->qz[gmap->icjckM[ii]],(*data)->h[ii],(*data)->dept[gmap->top2d[ii]],(*data)->dept[gmap->top2d[gmap->icjPkc[ii]]]);
        //     printf("                  Y - qy=%f,%f, (hl-h-hr)=(%f,%f,%f)\n",
        //         (*data)->qy[gmap->icjckM[ii]],(*data)->qy[ii],(*data)->h[gmap->icjMkc[ii]],(*data)->h[ii],(*data)->h[gmap->icjPkc[ii]]);
        //     printf("     -----     \n");
        // }

        // Warrick 1971
        // if (ii == 0)    {(*data)->wc[ii] = (*data)->wcs[ii];}
    }
}

// enforce moisture boundary conditions
void enforce_moisture_bc(Data **data, Map *gmap, Config *param)
{
    int ii;
    // xp and xm boundaries
    for (ii = 0; ii < param->ny*param->nz; ii++)
    {
        (*data)->wc[gmap->iPou[ii]] = (*data)->wc[gmap->iPin[ii]];
        (*data)->wc[gmap->iMou[ii]] = (*data)->wc[gmap->iMin[ii]];
    }
    // yp and ym boundaries
    for (ii = 0; ii < param->nx*param->nz; ii++)
    {
        (*data)->wc[gmap->jPou[ii]] = (*data)->wc[gmap->jPin[ii]];
        (*data)->wc[gmap->jMou[ii]] = (*data)->wc[gmap->jMin[ii]];
        if (param->bctype_GW[2] == 1)   {(*data)->wc[gmap->jPou[ii]] = (*data)->wcs[gmap->jPin[ii]];}
        if (param->bctype_GW[3] == 1)   {(*data)->wc[gmap->jMou[ii]] = (*data)->wcs[gmap->jMin[ii]];}
    }
    // zp boundary
    for (ii = 0; ii < param->nx*param->ny; ii++)
    {
        (*data)->wc[gmap->kPou[ii]] = (*data)->wc[gmap->kPin[ii]];
        (*data)->wc[gmap->kMou[ii]] = (*data)->wc[gmap->kMin[ii]];
    }
}

// >>>>> post-allocation of moisture <<<<<
void reallocate_water_content(Data **data, Map *gmap, Config *param, int irank)
{
    int ii, adj_sat, dir, alloc_type1=0, alloc_type2=0, alloc_type3=0, alloc_type4=0;
    double wcp, dV, dV_tot = 0.0;

    int flag = 0;
    for (ii = 0; ii < param->n3ci; ii++)
    {
        wcp = (*data)->wc[ii];
        if (gmap->actv[ii] == 1)
        {
            (*data)->wch[ii] = compute_wch(*data, ii, param);
            (*data)->hwc[ii] = compute_hwc(*data, ii, param);
            // over-saturated cell
            if ((*data)->wc[ii] >= (*data)->wcs[ii])
            {
                // send moisture
                if (param->post_allocate == 1)
                {
                    dV = ((*data)->wc[ii] - (*data)->wcs[ii]) * gmap->dz3d[ii] * param->dx * param->dy;
                    if (dV > 0.1 * (*data)->wcs[ii]*gmap->dz3d[ii]*param->dx*param->dy)
                    {
                        // printf("1 : ii = %d, wc = %f, dwc = %f\n",ii,(*data)->wc[ii],(*data)->wc[ii] - (*data)->wcs[ii]);
                        (*data)->repeat[0] = 1;
                    }
                    check_head_gradient(data, gmap, param, ii);
                    dV = allocate_send(data, gmap, param, ii, dV);
                    if (dV > 0)    {dV = allocate_send(data, gmap, param, ii, dV);}
                }
                dV_tot += dV;
                (*data)->wc[ii] = (*data)->wcs[ii];
                alloc_type1 += 1;
            }
            // unsaturated cell
            else
            {
                // check if a cell is adjacent to a saturated cell
                adj_sat = check_adj_sat(*data, gmap, param, ii);
                // isolated unsaturated cell
                if (adj_sat == 0)
                {
                    // if ((*data)->h[ii] < 0)
                    if ((*data)->wc[ii] < 0.9999*(*data)->wcs[ii])
                    {(*data)->h[ii] = (*data)->hwc[ii]; alloc_type2 += 1;}
                }
                // unsaturated cell adjacent to saturated cell
                else
                {
                    // receive moisture
                    if (param->post_allocate == 1)
                    {
                        if ((*data)->wch[ii] > (*data)->wc[ii])
                        {
                            dV = ((*data)->wch[ii] - (*data)->wc[ii]) * gmap->dz3d[ii] * param->dx * param->dy;
                            if (dV > 0.1 * (*data)->wcs[ii]*gmap->dz3d[ii]*param->dx*param->dy)
                            {(*data)->repeat[0] = 1;}
                            check_head_gradient(data, gmap, param, ii);
                            dV = allocate_recv(data, gmap, param, ii, dV);
                            // if (dV > 0) {dV = allocate_recv(data, gmap, param, ii, dV);}
                            alloc_type3 += 1;
                        }
                        // send moisture
                        else
                        {
                            dV = ((*data)->wc[ii] - (*data)->wch[ii]) * gmap->dz3d[ii] * param->dx * param->dy;
                            if (dV > 0.1 * (*data)->wcs[ii]*gmap->dz3d[ii]*param->dx*param->dy)
                            {
                                // printf("3 : ii = %d, wc = %f, dwc = %f\n",ii,(*data)->wc[ii],(*data)->wc[ii] - (*data)->wch[ii]);
                                (*data)->repeat[0] = 1;
                            }
                            check_head_gradient(data, gmap, param, ii);
                            dV = allocate_send(data, gmap, param, ii, dV);
                            if (dV > 0)    {dV = allocate_send(data, gmap, param, ii, dV);}
                            alloc_type4 += 1;

                        }
                    }

                    dV_tot += dV;
                    (*data)->wc[ii] = (*data)->wch[ii];
                    // if (dV < 1e-20) {
                    //     (*data)->wc[ii] = (*data)->wch[ii];
                    // }
                    // else {
                    //     (*data)->h[ii] = (*data)->hwc[ii];
                    // }
                    (*data)->room[ii] = ((*data)->wcs[ii] - (*data)->wc[ii]) * param->dx*param->dy*gmap->dz3d[ii];
                }
            }
        }
    }

    // for (ii = 0; ii < param->n3ci; ii++)
    // {
    //     if (gmap->ii[ii] == 45 & gmap->jj[ii] == 120 & gmap->istop[gmap->icjckM[ii]] == 1)
    //     {
    //         printf(" ----- ii=%d, Final wc = %f,   h = %f\n",ii,(*data)->wc[ii],(*data)->h[ii]);
    //         printf(" -----        dz=%f, dept=%f, wckM=%f, hkM=%f, seepage=%f\n",gmap->dz3d[ii],(*data)->dept[gmap->top2d[ii]], \
    //             (*data)->wc[gmap->icjckM[ii]],(*data)->h[gmap->icjckM[ii]],1e3*(*data)->qz[gmap->icjckM[gmap->icjckM[ii]]]);
    //         printf(" ========== \n\n");
    //     }
    // }
    // printf(" >> Total number of alloc type1=%d, type2=%d, type3=%d, type4=%d\n",alloc_type1,alloc_type2,alloc_type3,alloc_type4);
    // if (dV_tot/(param->dx*param->dy*gmap->dz3d[0]) > 0.005)
    // {
    //     printf(" >> Total volume loss = %f\n",dV_tot/(param->dx*param->dy*gmap->dz3d[0]));
    // }

}

// >>>>> Check if a grid cell is adjacent to a saturated cell <<<<<
int check_adj_sat(Data *data, Map *gmap, Config *param, int ii)
{
    int adj_sat = 0;
    if (gmap->actv[gmap->icjckP[ii]] == 1 & data->wc[gmap->icjckP[ii]] >= data->wcs[gmap->icjckP[ii]])
    {adj_sat = 1;}
    else if (data->Kx[ii] != 0 & data->wc[gmap->iPjckc[ii]] >= data->wcs[gmap->iPjckc[ii]])
    {adj_sat = 1;}
    else if (data->Kx[gmap->iMjckc[ii]] != 0 & data->wc[gmap->iMjckc[ii]] >= data->wcs[gmap->iMjckc[ii]])
    {adj_sat = 1;}
    else if (data->Ky[ii] != 0 & data->wc[gmap->icjPkc[ii]] >= data->wcs[gmap->icjPkc[ii]])
    {adj_sat = 1;}
    else if (data->Ky[gmap->icjMkc[ii]] != 0 & data->wc[gmap->icjMkc[ii]] >= data->wcs[gmap->icjMkc[ii]])
    {adj_sat = 1;}
    else
    {
        if (gmap->istop[ii] == 1)
        {
            if (param->sim_shallowwater == 1 & data->dept[gmap->top2d[ii]] > 0)    {adj_sat = 1;}
            else if (param->sim_shallowwater == 0 & param->bctype_GW[5] == 1 & param->htop >= 0)   {adj_sat = 1;}
        }
        else
        {
            if (gmap->actv[gmap->icjckM[ii]] == 1 & data->wc[gmap->icjckM[ii]] >= data->wcs[gmap->icjckM[ii]])
            {adj_sat = 1;}
        }
    }
    return adj_sat;
}

// >>>>> check the direction of allocation
void check_head_gradient(Data **data, Map *gmap, Config *param, int ii)
{
    int jj;
    double dh, dzp, dzm, gradx, grady, gradz, grad_tot;
    // get dh at the iP face
    dh = ((*data)->h[gmap->iPjckc[ii]]-(*data)->h[ii]) / param->dx;
    if (gmap->actv[gmap->iPjckc[ii]] == 1 & (*data)->Kx[ii] > 0)   {(*data)->dh6[0] = dh;}
    else    {(*data)->dh6[0] = 0.0;}
    // get dh at the iM face
    dh = ((*data)->h[ii]-(*data)->h[gmap->iMjckc[ii]]) / param->dx;
    if (gmap->actv[gmap->iMjckc[ii]] == 1 & (*data)->Kx[gmap->iMjckc[ii]] > 0) {(*data)->dh6[1] = dh;}
    else    {(*data)->dh6[1] = 0.0;}
    // get dh at the jP face
    dh = ((*data)->h[gmap->icjPkc[ii]]-(*data)->h[ii]) / param->dy;
    if (gmap->actv[gmap->icjPkc[ii]] == 1 & (*data)->Ky[ii] > 0)   {(*data)->dh6[2] = dh;}
    else    {(*data)->dh6[2] = 0.0;}
    // get dh at the jM face
    dh = ((*data)->h[ii]-(*data)->h[gmap->icjMkc[ii]]) / param->dy;
    if (gmap->actv[gmap->icjMkc[ii]] == 1 & (*data)->Ky[gmap->icjMkc[ii]] > 0) {(*data)->dh6[3] = dh;}
    else    {(*data)->dh6[3] = 0.0;}
    // get dh at the kP face
    if (gmap->kk[ii] == param->nz-1)
    {if (param->bctype_GW[4] == 1)   {dzp = 0.5 * gmap->dz3d[ii];}}
    else    {dzp = 0.5 * (gmap->dz3d[ii] + gmap->dz3d[gmap->icjckP[ii]]);}
    if (dzp == 0.0) {dh = 0.0;}
    else    {dh = ((*data)->h[gmap->icjckP[ii]]-(*data)->h[ii]) / dzp - 1.0;}
    if (gmap->actv[gmap->icjckP[ii]] == 1 & (*data)->Kz[ii] > 0)    {(*data)->dh6[4] = dh;}
    else    {(*data)->dh6[4] = 0.0;}
    // get dh at the kM face
    if (gmap->istop[ii] == 1)
    {
        if (param->sim_shallowwater == 1)
        {if ((*data)->dept[gmap->top2d[ii]] > 0)    {dzm = 0.5 * gmap->dz3d[ii];}}
        else
        {if (param->bctype_GW[5] == 1)   {dzm = 0.5 * gmap->dz3d[ii];}}
    }
    else    {dzm = 0.5 * (gmap->dz3d[ii] + gmap->icjckM[ii]);}
    if (dzm == 0.0) {dh = 0.0;}
    else    {dh = ((*data)->h[ii]-(*data)->h[gmap->icjckM[ii]]) / dzm - 1.0;}
    if (gmap->actv[gmap->icjckM[ii]] == 1 & (*data)->Kz[gmap->icjckM[ii]] > 0)    {(*data)->dh6[5] = dh;}
    else    {(*data)->dh6[5] = 0.0;}

    if (param->use_full3d == 0)
    {
        (*data)->dh6[0] = 0.0;
        (*data)->dh6[1] = 0.0;
        (*data)->dh6[2] = 0.0;
        (*data)->dh6[3] = 0.0;
    }

    // combine gradients in each direction
    grad_tot = 0.0;
    if ((*data)->dh6[0]*(*data)->dh6[1] >= 0)
    {gradx = 0.5*((*data)->dh6[0]+(*data)->dh6[1]);     grad_tot += fabs(gradx);}
    else    {gradx = 0.0;   grad_tot += (fabs((*data)->dh6[0]) + fabs((*data)->dh6[1]));}
    if ((*data)->dh6[2]*(*data)->dh6[3] >= 0)
    {grady = 0.5*((*data)->dh6[2]+(*data)->dh6[3]);     grad_tot += fabs(grady);}
    else    {grady = 0.0;   grad_tot += (fabs((*data)->dh6[2]) + fabs((*data)->dh6[3]));}
    if ((*data)->dh6[4]*(*data)->dh6[5] >= 0)
    {gradz = 0.5*((*data)->dh6[4]+(*data)->dh6[5]);     grad_tot += fabs(gradz);}
    else    {gradz = 0.0;   grad_tot += (fabs((*data)->dh6[4]) + fabs((*data)->dh6[5]));}
    // get moisture split ratio
    if (grad_tot > 0)
    {
        if (gradx > 0)  {(*data)->rsplit[0] = 0.0;  (*data)->rsplit[1] = gradx / grad_tot;}
        else if (gradx < 0) {(*data)->rsplit[1] = 0.0;  (*data)->rsplit[0] = -gradx / grad_tot;}
        else    {(*data)->rsplit[0] = fabs((*data)->dh6[0])/grad_tot;   (*data)->rsplit[1] = fabs((*data)->dh6[1])/grad_tot;}
        if (grady > 0)  {(*data)->rsplit[2] = 0.0;  (*data)->rsplit[3] = grady / grad_tot;}
        else if (grady < 0) {(*data)->rsplit[3] = 0.0;  (*data)->rsplit[2] = -grady / grad_tot;}
        else    {(*data)->rsplit[2] = fabs((*data)->dh6[2])/grad_tot;   (*data)->rsplit[3] = fabs((*data)->dh6[3])/grad_tot;}
        if (gradz > 0)  {(*data)->rsplit[4] = 0.0;  (*data)->rsplit[5] = gradz / grad_tot;}
        else if (gradz < 0) {(*data)->rsplit[5] = 0.0;  (*data)->rsplit[4] = -gradz / grad_tot;}
        else    {(*data)->rsplit[4] = fabs((*data)->dh6[4])/grad_tot;   (*data)->rsplit[5] = fabs((*data)->dh6[5])/grad_tot;}
    }
    else
    {for (jj = 0; jj < 6; jj++)  {(*data)->rsplit[jj] = 0.0;}}
}

// >>>>> send moisture <<<<<
double allocate_send(Data **data, Map *gmap, Config *param, int ii, double dV)
{
    int ll, rev = 0;
    double dVxp=0, dVxm=0, dVyp=0, dVym=0, dVzp=0, dVzm=0, temp, Vres=0.0;
    // send up
    if ((*data)->rsplit[5] > 0)
    {
        dVzm = dV * (*data)->rsplit[5];
        ll = ii;
        while (gmap->istop[ll] != 1)
        {
            ll -= 1;
            if ((*data)->room[ll] > 0)
            {
                if ((*data)->room[ll] > dVzm)
                {
                    (*data)->wc[ll] += dVzm / (param->dx*param->dy*gmap->dz3d[ll]);
                    (*data)->room[ll] -= dVzm;
                    (*data)->qz[ll] += dVzm / (param->dx * param->dy * param->dt);
                    dVzm = 0.0;
                    break;
                }
                else
                {
                    (*data)->wc[ll] += (*data)->room[ll] / (param->dx*param->dy*gmap->dz3d[ll]);
                    dVzm -= (*data)->room[ll];
                    (*data)->qz[ll] += (*data)->room[ll] / (param->dx * param->dy * param->dt);
                    (*data)->room[ll] = 0.0;
                }
            }
        }
        // release excess moisture from top surface
        if (dVzm > 0.0)
        {
            if (param->sim_shallowwater == 1)
            {
                // Here, if seepage is too small, we simply neglect it. ZhiLi20201211
                if ((*data)->dept[gmap->top2d[ll]] > 0)
                {
                    // Kuan2019
                    (*data)->dept[gmap->top2d[ll]] += dVzm / (param->dx*param->dy);
                    (*data)->eta[gmap->top2d[ll]] += dVzm / (param->dx*param->dy);
                    (*data)->qz[gmap->icjckM[ll]] += dVzm / (param->dx * param->dy * param->dt);
                }
                else
                {
                    if (dVzm / (param->dx*param->dy) > param->min_dept)
                    {
                        (*data)->dept[gmap->top2d[ll]] += dVzm / (param->dx*param->dy);
                        (*data)->eta[gmap->top2d[ll]] += dVzm / (param->dx*param->dy);
                        (*data)->qz[gmap->icjckM[ll]] += dVzm / (param->dx * param->dy * param->dt);
                    }
                }
                dVzm = 0.0;
            }
            else
            {
                if (param->bctype_GW[5] == 1)
                {
                    (*data)->qz[gmap->icjckM[ll]] += dVzm / (param->dx * param->dy * param->dt);
                    dVzm = 0.0;
                }
                // if top surface is no-flux or fixed q, send extra moisture back down
                else
                {
                    while (gmap->kk[ll] != param->nz-1)
                    {
                        ll += 1;
                        if ((*data)->room[ll] > 0)
                        {
                            if ((*data)->room[ll] > dVzm)
                            {
                                (*data)->wc[ll] += dVzm / (param->dx*param->dy*gmap->dz3d[ll]);
                                (*data)->qz[gmap->icjckM[ll]] -= dVzm / (param->dx * param->dy * param->dt);
                                (*data)->room[ll] -= dVzm;
                                dVzm = 0.0;
                                break;
                            }
                            else
                            {
                                (*data)->wc[ll] += (*data)->room[ll] / (param->dx*param->dy*gmap->dz3d[ll]);
                                (*data)->qz[gmap->icjckM[ll]] -= (*data)->room[ll] / (param->dx * param->dy * param->dt);
                                dVzm -= (*data)->room[ll];
                                (*data)->room[ll] = 0.0;
                            }
                        }
                    }
                }
            }
            // Print warning if too much water remains after redistribution
            if (dVzm / (param->dx*param->dy*gmap->dz3d[ii]) > 2e-2)
            {
                // printf(" WARNING : extra moisture amount %f during 'send up' for cell (%d, %d, %d)\n", \
                    dVzm / (param->dx*param->dy*gmap->dz3d[ii]), gmap->ii[ii],gmap->jj[ii],gmap->kk[ii]);
                dVzm = 0.0;
            }
        }
    }
    // send down
    if ((*data)->rsplit[4] > 0)
    {
        dVzp = dV * (*data)->rsplit[4];
        ll = ii;
        while (gmap->kk[ll] != param->nz-1)
        {
            ll += 1;
            if ((*data)->room[ll] > 0)
            {
                if ((*data)->room[ll] > dVzp)
                {
                    (*data)->wc[ll] += dVzp / (param->dx*param->dy*gmap->dz3d[ll]);
                    (*data)->qz[gmap->icjckM[ll]] -= dVzp / (param->dx * param->dy * param->dt);
                    (*data)->room[ll] -= dVzp;
                    dVzp = 0.0;
                    break;
                }
                else
                {
                    (*data)->wc[ll] += (*data)->room[ll] / (param->dx*param->dy*gmap->dz3d[ll]);
                    (*data)->qz[gmap->icjckM[ll]] -= (*data)->room[ll] / (param->dx * param->dy * param->dt);
                    dVzp -= (*data)->room[ll];
                    (*data)->room[ll] = 0.0;
                }
            }
        }
        // release excess moisture from bottom boundary
        if (dVzp > 0.0)
        {
            if (param->bctype_GW[5] == 1)   {dVzp = 0.0;}
            else
            {
                if (dVzp / (param->dx*param->dy*gmap->dz3d[ii]) > 2e-2)
                {
                    // printf(" WARNING : extra moisture amount %f during 'send down' for cell (%d, %d, %d)\n", \
                        dVzp / (param->dx*param->dy*gmap->dz3d[ii]), gmap->ii[ii],gmap->jj[ii],gmap->kk[ii]);
                    dVzp = 0.0;
                }
            }
        }
    }
    // send in x direction
    if ((*data)->rsplit[0] > 0)
    {
        dVxp = dV * (*data)->rsplit[0];
        ll = gmap->iPjckc[ii];
        if ((*data)->room[ll] > 0)
        {
            if ((*data)->room[ll] > dVxp)
            {
                (*data)->wc[ll] += dVxp / (param->dx*param->dy*gmap->dz3d[ll]);
                (*data)->qx[gmap->iMjckc[ll]] -= dVxp / (gmap->dz3d[ll] * param->dy * param->dt);
                (*data)->room[ll] -= dVxp;
                dVxp = 0.0;
            }
            else
            {
                (*data)->wc[ll] += (*data)->room[ll] / (param->dx*param->dy*gmap->dz3d[ll]);
                (*data)->qx[gmap->iMjckc[ll]] -= (*data)->room[ll] / (gmap->dz3d[ll] * param->dy * param->dt);
                dVxp -= (*data)->room[ll];
                (*data)->room[ll] = 0.0;
            }
        }
    }
    if ((*data)->rsplit[1] > 0)
    {
        dVxm = dV * (*data)->rsplit[1];
        ll = gmap->iMjckc[ii];
        if ((*data)->room[ll] > 0)
        {
            if ((*data)->room[ll] > dVxm)
            {
                (*data)->wc[ll] += dVxm / (param->dx*param->dy*gmap->dz3d[ll]);
                (*data)->qx[ll] += dVxm / (gmap->dz3d[ll] * param->dy * param->dt);
                (*data)->room[ll] -= dVxm;
                dVxm = 0.0;
            }
            else
            {
                (*data)->wc[ll] += (*data)->room[ll] / (param->dx*param->dy*gmap->dz3d[ll]);
                (*data)->qx[ll] += (*data)->room[ll] / (gmap->dz3d[ll] * param->dy * param->dt);
                dVxm -= (*data)->room[ll];
                (*data)->room[ll] = 0.0;
            }
        }
    }
    // send in y direction
    if ((*data)->rsplit[2] > 0)
    {
        dVyp = dV * (*data)->rsplit[2];
        ll = gmap->icjPkc[ii];
        if ((*data)->room[ll] > 0)
        {
            if ((*data)->room[ll] > dVyp)
            {
                (*data)->wc[ll] += dVyp / (param->dx*param->dy*gmap->dz3d[ll]);
                (*data)->qy[gmap->icjMkc[ll]] -= dVyp / (gmap->dz3d[ll] * param->dx * param->dt);
                (*data)->room[ll] -= dVyp;
                dVyp = 0.0;
            }
            else
            {
                (*data)->wc[ll] += (*data)->room[ll] / (param->dx*param->dy*gmap->dz3d[ll]);
                (*data)->qy[gmap->icjMkc[ll]] -= (*data)->room[ll] / (gmap->dz3d[ll] * param->dx * param->dt);
                dVyp -= (*data)->room[ll];
                (*data)->room[ll] = 0.0;
            }
        }
        if (dVyp > 0)
        {
            if (dVyp / (param->dx*param->dy*gmap->dz3d[ii]) > 1e-2)
            {
                // printf(" WARNING : extra moisture amount %f during 'send yp' for cell (%d, %d, %d)\n", \
                    dVyp / (param->dx*param->dy*gmap->dz3d[ii]), gmap->ii[ii],gmap->jj[ii],gmap->kk[ii]);
                dVyp = 0.0;
            }
        }
    }
    if ((*data)->rsplit[3] > 0)
    {
        dVym = dV * (*data)->rsplit[3];
        ll = gmap->icjMkc[ii];
        if ((*data)->room[ll] > 0)
        {
            if ((*data)->room[ll] > dVym)
            {
                (*data)->wc[ll] += dVym / (param->dx*param->dy*gmap->dz3d[ll]);
                (*data)->qy[ll] += dVym / (gmap->dz3d[ll] * param->dx * param->dt);
                (*data)->room[ll] -= dVym;
                dVym = 0.0;
            }
            else
            {
                (*data)->wc[ll] += (*data)->room[ll] / (param->dx*param->dy*gmap->dz3d[ll]);
                (*data)->qy[ll] += (*data)->room[ll] / (gmap->dz3d[ll] * param->dx * param->dt);
                dVym -= (*data)->room[ll];
                (*data)->room[ll] = 0.0;
            }
        }
        if (dVym > 0)
        {
            if (dVym / (param->dx*param->dy*gmap->dz3d[ii]) > 1e-2)
            {
                // printf(" WARNING : extra moisture amount %f during 'send ym' for cell (%d, %d, %d)\n", \
                    dVym / (param->dx*param->dy*gmap->dz3d[ii]), gmap->ii[ii],gmap->jj[ii],gmap->kk[ii]);
                dVym = 0.0;
            }
        }
    }
    // Check if reversed send is needed
    if (dVzp + dVzm > 0)
    {
        temp = (*data)->rsplit[5];
        (*data)->rsplit[5] = (*data)->rsplit[4];
        (*data)->rsplit[4] = temp;
        Vres += (dVzp + dVzm);
    }
    if (dVxp + dVxm > 0)
    {
        temp = (*data)->rsplit[1];
        (*data)->rsplit[1] = (*data)->rsplit[0];
        (*data)->rsplit[0] = temp;
        Vres += (dVxp + dVxm);
    }
    if (dVyp + dVym > 0)
    {
        temp = (*data)->rsplit[3];
        (*data)->rsplit[3] = (*data)->rsplit[2];
        (*data)->rsplit[2] = temp;
        Vres += (dVyp + dVym);
    }
    return Vres;
}


// >>>>> send moisture <<<<<
double allocate_recv(Data **data, Map *gmap, Config *param, int ii, double dV)
{
    int ll, dir, rev = 0;
    double dVxp=0, dVxm=0, dVyp=0, dVym=0, dVzp=0, dVzm=0, temp, Vres=0.0;
    // recv from up
    if ((*data)->rsplit[5] > 0)
    {
        dVzm = dV * (*data)->rsplit[5];
        ll = ii;
        while (gmap->istop[ll] != 1)
        {
            ll -= 1;
            if ((*data)->wc[ll] > (*data)->wcr[ll])
            {
                if (((*data)->wc[ll]-(*data)->wcr[ll]) > dVzm / (param->dx*param->dy*gmap->dz3d[ll]))
                {
                    (*data)->wc[ll] -= dVzm / (param->dx*param->dy*gmap->dz3d[ll]);
                    (*data)->qz[ll] -= dVzm / (param->dx * param->dy * param->dt);
                    (*data)->room[ll] += dVzm;
                    dVzm = 0.0;
                    break;
                }
                else
                {
                    dVzm -= ((*data)->wc[ll]-(*data)->wcr[ll]) * (param->dx*param->dy*gmap->dz3d[ll]);
                    (*data)->qz[ll] -= ((*data)->wc[ll]-(*data)->wcr[ll]) / (param->dx * param->dy * param->dt);
                    (*data)->room[ll] += ((*data)->wc[ll]-(*data)->wcr[ll]) * (param->dx*param->dy*gmap->dz3d[ll]);
                    (*data)->wc[ll] = (*data)->wcr[ll];
                }
            }
            // limite recv to 1 adjacent cell, ZhiLi20200827
            dVzm = 0.0;
            break;
        }
        // extract moisture from surface domain
        if (dVzm > 0.0)
        {
            if (param->sim_shallowwater == 1)
            {
                if ((*data)->dept[gmap->top2d[ll]] > dVzm / (param->dx*param->dy))
                {
                    (*data)->dept[gmap->top2d[ll]] -= dVzm / (param->dx*param->dy);
                    // (*data)->qz[gmap->icjckM[ll]] -= dVzm / (param->dx * param->dy * param->dt);
                    dVzm = 0.0;
                }
                else
                {(*data)->dept[gmap->top2d[ll]] = param->min_dept;}
            }
            else    {if (param->bctype_GW[5] == 1)   {dVzm = 0.0;}}
            dVzm = 0.0;
        }
    }
    // recv from bottom
    if ((*data)->rsplit[4] > 0)
    {
        dVzp = dV * (*data)->rsplit[4];
        ll = ii;
        while (gmap->kk[ll] != param->nz-1)
        {
            ll += 1;
            if ((*data)->wc[ll] > (*data)->wcr[ll])
            {
                if (((*data)->wc[ll]-(*data)->wcr[ll]) > dVzp / (param->dx*param->dy*gmap->dz3d[ll]))
                {
                    (*data)->wc[ll] -= dVzp / (param->dx*param->dy*gmap->dz3d[ll]);
                    (*data)->qz[gmap->icjckM[ll]] += dVzp / (param->dx * param->dy * param->dt);
                    (*data)->room[ll] += dVzp;
                    dVzp = 0.0;
                    break;
                }
                else
                {
                    dVzp -= ((*data)->wc[ll]-(*data)->wcr[ll]) * (param->dx*param->dy*gmap->dz3d[ll]);
                    (*data)->qz[gmap->icjckM[ll]] += ((*data)->wc[ll]-(*data)->wcr[ll]) / (param->dx * param->dy * param->dt);
                    (*data)->room[ll] += ((*data)->wc[ll]-(*data)->wcr[ll]) * (param->dx*param->dy*gmap->dz3d[ll]);
                    (*data)->wc[ll] = (*data)->wcr[ll];
                    // Limit recv within 1 neighbor cell to avoid instability! ZhiLi20200621
                    dVzp = 0.0;
                }
            }
            // limite recv to 1 adjacent cell, ZhiLi20200827
            dVzp = 0.0;
            break;
        }
        if (dVzp > 0.0) {if (param->bctype_GW[4] == 1)   {dVzp = 0.0;}}
    }
    // // send in x direction
    // if ((*data)->rsplit[0] > 0)
    // {
    //     dVxp = dV * (*data)->rsplit[0];
    //     ll = gmap->iPjckc[ii];
    //     if ((*data)->wc[ll] > (*data)->wcr[ll])
    //     {
    //         if (((*data)->wc[ll]-(*data)->wcr[ll]) > dVxp / (param->dx*param->dy*gmap->dz3d[ll]))
    //         {
    //             (*data)->wc[ll] -= dVxp / (param->dx*param->dy*gmap->dz3d[ll]);
    //             (*data)->qx[gmap->iMjckc[ll]] += dVxp / (gmap->dz3d[ll] * param->dy * param->dt);
    //             (*data)->room[ll] += dVxp;
    //             dVxp = 0.0;
    //         }
    //         else
    //         {
    //             dVxp -= ((*data)->wc[ll]-(*data)->wcr[ll]) * (param->dx*param->dy*gmap->dz3d[ll]);
    //             (*data)->qx[gmap->iMjckc[ll]] += ((*data)->wc[ll]-(*data)->wcr[ll]) / (gmap->dz3d[ll] * param->dy * param->dt);
    //             (*data)->room[ll] += ((*data)->wc[ll]-(*data)->wcr[ll]) * (param->dx*param->dy*gmap->dz3d[ll]);
    //             (*data)->wc[ll] = (*data)->wcr[ll];
    //         }
    //     }
    // }
    // if ((*data)->rsplit[1] > 0)
    // {
    //     dVxm = dV * (*data)->rsplit[1];
    //     ll = gmap->iMjckc[ii];
    //     if ((*data)->wc[ll] > (*data)->wcr[ll])
    //     {
    //         if (((*data)->wc[ll]-(*data)->wcr[ll]) > dVxm / (param->dx*param->dy*gmap->dz3d[ll]))
    //         {
    //             (*data)->wc[ll] -= dVxm / (param->dx*param->dy*gmap->dz3d[ll]);
    //             (*data)->qx[ll] -= dVxm / (gmap->dz3d[ll] * param->dy * param->dt);
    //             (*data)->room[ll] += dVxm;
    //             dVxm = 0.0;
    //         }
    //         else
    //         {
    //             dVxm -= ((*data)->wc[ll]-(*data)->wcr[ll]) * (param->dx*param->dy*gmap->dz3d[ll]);
    //             (*data)->qx[ll] -= ((*data)->wc[ll]-(*data)->wcr[ll]) / (gmap->dz3d[ll] * param->dy * param->dt);
    //             (*data)->room[ll] += ((*data)->wc[ll]-(*data)->wcr[ll]) * (param->dx*param->dy*gmap->dz3d[ll]);
    //             (*data)->wc[ll] = (*data)->wcr[ll];
    //         }
    //     }
    // }
    // send in y direction
    if ((*data)->rsplit[2] > 0)
    {
        dVyp = dV * (*data)->rsplit[2];
        ll = gmap->icjPkc[ii];
        if ((*data)->wc[ll] > (*data)->wcr[ll])
        {
            if (((*data)->wc[ll]-(*data)->wcr[ll]) > dVyp / (param->dx*param->dy*gmap->dz3d[ll]))
            {
                (*data)->wc[ll] -= dVyp / (param->dx*param->dy*gmap->dz3d[ll]);
                (*data)->qy[gmap->icjMkc[ll]] += dVyp / (gmap->dz3d[ll] * param->dx * param->dt);
                (*data)->room[ll] += dVyp;
                dVyp = 0.0;
            }
            else
            {
                dVyp -= ((*data)->wc[ll]-(*data)->wcr[ll]) * (param->dx*param->dy*gmap->dz3d[ll]);
                (*data)->qy[gmap->icjMkc[ll]] += ((*data)->wc[ll]-(*data)->wcr[ll]) / (gmap->dz3d[ll] * param->dx * param->dt);
                (*data)->room[ll] += ((*data)->wc[ll]-(*data)->wcr[ll]) * (param->dx*param->dy*gmap->dz3d[ll]);
                (*data)->wc[ll] = (*data)->wcr[ll];
            }
        }
        // limite recv to 1 adjacent cell, ZhiLi20200827
        dVyp = 0.0;
    }
    if ((*data)->rsplit[3] > 0)
    {
        dVym = dV * (*data)->rsplit[3];
        ll = gmap->icjMkc[ii];
        if ((*data)->wc[ll] > (*data)->wcr[ll])
        {
            if (((*data)->wc[ll]-(*data)->wcr[ll]) > dVym / (param->dx*param->dy*gmap->dz3d[ll]))
            {
                (*data)->wc[ll] -= dVym / (param->dx*param->dy*gmap->dz3d[ll]);
                (*data)->qy[ll] -= dVym / (gmap->dz3d[ll] * param->dx * param->dt);
                (*data)->room[ll] += dVym;
                dVym = 0.0;
            }
            else
            {
                dVym -= ((*data)->wc[ll]-(*data)->wcr[ll]) * (param->dx*param->dy*gmap->dz3d[ll]);
                (*data)->qy[ll] -= ((*data)->wc[ll]-(*data)->wcr[ll]) / (gmap->dz3d[ll] * param->dx * param->dt);
                (*data)->room[ll] += ((*data)->wc[ll]-(*data)->wcr[ll]) * (param->dx*param->dy*gmap->dz3d[ll]);
                (*data)->wc[ll] = (*data)->wcr[ll];
            }
        }
        // limite recv to 1 adjacent cell, ZhiLi20200827
        dVym = 0.0;
    }
    // Check if reversed send is needed
    if (dVzp + dVzm > 0)
    {
        temp = (*data)->rsplit[5];
        (*data)->rsplit[5] = (*data)->rsplit[4];
        (*data)->rsplit[4] = temp;
        Vres += (dVzp + dVzm);
    }
    if (dVxp + dVxm > 0)
    {
        temp = (*data)->rsplit[1];
        (*data)->rsplit[1] = (*data)->rsplit[0];
        (*data)->rsplit[0] = temp;
        Vres += (dVxp + dVxm);
    }
    if (dVyp + dVym > 0)
    {
        temp = (*data)->rsplit[3];
        (*data)->rsplit[3] = (*data)->rsplit[2];
        (*data)->rsplit[2] = temp;
        Vres += (dVyp + dVym);
    }
    return Vres;
}

// >>>>> update cell volume using flux
void volume_by_flux_subs(Data **data, Map *gmap, Config *param)
{
    int ii, jj, kk;
    for (ii = 0; ii < param->n3ci; ii++)
    {
      (*data)->Vgflux[ii] = (*data)->Vg[ii] + param->dt *
        ((-(*data)->qx[gmap->iMjckc[ii]] + (*data)->qx[ii]) * param->dy * gmap->dz3d[ii] +
        (-(*data)->qy[gmap->icjMkc[ii]] + (*data)->qy[ii]) * param->dx * gmap->dz3d[ii] +
        (-(*data)->qz[gmap->icjckM[ii]] + (*data)->qz[ii]) * param->dy * param->dx);
    }
}


// >>>>> adjust dt for groundwater solver <<<<<
void adaptive_time_step(Data *data, Map *gmap, Config **param, int root, int irank)
{
    double qin, qou, r_red, r_inc, dq, dq_max, dKdwc, dt_Co, dt_Comin, wcmax, hmax, rmax;
    int ii, imax = 0, jmax = 0, kmax = 0, unsat = 0, sat = 0, itop=0;
    r_red = 0.9;
    r_inc = 1.1;
    dq_max = 0.0;
    dt_Comin = 1e8;
    // check if both sat and unsat zones exist
    for (ii = 0; ii < (*param)->n3ci; ii++)
    {if (data->wc[ii] >= (*param)->wcs) {sat = 1;  break;}}
    for (ii = 0; ii < (*param)->n3ci; ii++)
    {if (data->wc[ii] < (*param)->wcs) {unsat = 1;  break;}}
    // adjust dt
    for (ii = 0; ii < (*param)->n3ci; ii++)
    {
        if (gmap->actv[ii] == 1 & gmap->istop[ii] == 0)
        {
            qin = data->qx[ii] + data->qy[ii] + data->qz[ii];
            qou = data->qx[gmap->iMjckc[ii]] + data->qy[gmap->icjMkc[ii]] + data->qz[gmap->icjckM[ii]];
            dq = fabs(qin - qou) * (*param)->dt / gmap->dz3d[ii];
            if (dq > dq_max)    {
                dq_max = dq; imax=gmap->ii[ii]; jmax=gmap->jj[ii]; kmax=gmap->kk[ii]; itop=gmap->istop[ii];
                wcmax = data->wc[ii];    hmax = data->h[ii];  rmax = data->qtop[gmap->top2d[ii]];
            }
            if (data->wc[ii] < (*param)->wcs)
            {
                dKdwc = compute_dKdwc(data, data->Ksz, ii, *param);
                dt_Co = (*param)->Co_max * gmap->dz3d[ii] / dKdwc;
                if (dt_Co < dt_Comin)   {dt_Comin = dt_Co;}
            }
        }
    }
    if (dq_max > 0.02)
    {
        (*param)->dt = (*param)->dt * r_red;
    }
    else if (dq_max < 0.01) {(*param)->dt = (*param)->dt * r_inc;}
    if ((*param)->dt > (*param)->dt_max)  {(*param)->dt = (*param)->dt_max;}
    if ((*param)->dt < (*param)->dt_min)
    {
        (*param)->dt = (*param)->dt_min;
        printf(" >>> Min dt is reached at (%d,%d,%d): dt=%f, dqmax=%f, itop=%d, wc=%f, h=%f, rain=%f\n",
            imax,jmax,kmax,(*param)->dt, dq_max,itop, wcmax, hmax, rmax);
    }
    // if (sat == 1 & unsat == 0) {(*param)->dt = (*param)->dt_max;}
    // apply the Courant number criteria
    // if ((*param)->dt > dt_Comin & sat == 1 & unsat == 1)   {(*param)->dt = dt_Comin;}
    if ((*param)->dt > dt_Comin)   {(*param)->dt = dt_Comin;}
    if ((*param)->dt > (*param)->dt_max)  {(*param)->dt = (*param)->dt_max;}
    if ((*param)->dt < (*param)->dt_min)  {
        (*param)->dt = (*param)->dt_min;
    }
    // printf("  >> dt by CO = %f\n",dt_Comin);
    // unify dt for all processes
    if ((*param)->use_mpi == 1)
    {
        mpi_gather_double(&(*param)->dt_root[0], &(*param)->dt, 1, root);
        (*param)->dt = getMin((*param)->dt_root, (*param)->mpi_nx*(*param)->mpi_ny);
        mpi_bcast_double(&(*param)->dt, 1, root);
    }
}
