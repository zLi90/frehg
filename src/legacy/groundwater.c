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
#include "subroutines.h"
#include "utility.h"

// #include "../laspack/vector.h"
// #include "../laspack/qmatrix.h"
// #include "../laspack/rtc.h"

#include "../laspack/errhandl.h"
#include "../laspack/vector.h"
#include "../laspack/matrix.h"
#include "../laspack/qmatrix.h"
#include "../laspack/operats.h"
#include "../laspack/factor.h"
#include "../laspack/precond.h"
#include "../laspack/eigenval.h"
#include "../laspack/rtc.h"
#include "../laspack/itersolv.h"
#include "../laspack/mlsolv.h"

void solve_groundwater(Data **data, Map *smap, Map *gmap, Config *param, int irank, int nrank);
void compute_K_face(Data **data, Map *gmap, Config *param, int irank, int nrank);
void baroclinic_face(Data **data, Map *gmap, Config *param, int irank, int nrank);
void jacobian_mat_coeff(Data **data, Map *gmap, Config *param, int irank);
void newton_iter(Data **data, Map *gmap, QMatrix A, Vector b, Vector x, Config *param, int irank);
void groundwater_mat_coeff(Data **data, Map *gmap, Config *param, int irank);
void groundwater_rhs(Data **data, Map *gmap, Config *param, int irank);
void build_groundwater_system(Data *data, Map *gmap, Config *param, QMatrix A, Vector b);
void solve_groundwater_system(Data **data, Map *gmap, QMatrix A, Vector b, Vector x, Config *param);
void enforce_head_bc(Data **data, Map *gmap, Config *param);
void groundwater_flux(Data **data, Map *gmap, Config *param, int irank);
void pseudo_seepage(Data **data, Map *gmap, Config *param, double dt);
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
void solve_groundwater(Data **data, Map *smap, Map *gmap, Config *param, int irank, int nrank)
{
    int ii;
    // allocate linear system
    QMatrix A;
    Q_Constr(&A, "A", gmap->nactv, False, Rowws, Normal, True);
    Vector b;
    V_Constr(&b, "b", gmap->nactv, Normal, True);
    Vector x;
    V_Constr(&x, "x", gmap->nactv, Normal, True);

    for (ii = 0; ii < param->n3ct; ii++)
    {
        (*data)->hn[ii] = (*data)->h[ii];
        (*data)->hwc[ii] = compute_hwc(*data, ii, param);
        (*data)->wcn[ii] = (*data)->wc[ii];
        (*data)->wch[ii] = compute_wch(*data, (*data)->h[ii], ii, param);
        (*data)->ch[ii] = compute_ch(*data, ii, param);
        // if (ii < 0)    {printf("wc = %f\n",(*data)->wc[ii]);}
    }

    if (param->use_mpi == 1)
    {
        mpi_exchange_subsurf((*data)->hn, gmap, 2, param, irank, nrank);
        mpi_exchange_subsurf((*data)->wcn, gmap, 2, param, irank, nrank);
    }

    // >>> Predictor step
    compute_K_face(data, gmap, param, irank, nrank);
    baroclinic_face(data, gmap, param, irank, nrank);

    // groundwater_mat_coeff(data, gmap, param, irank);
    // groundwater_rhs(data, gmap, param, irank);
    if (param->iter_solve == 1) {
        newton_iter(data, gmap, A, b, x, param, irank);
        if (param->converge == 0)  {
            mpi_print("   >> Switching to PCA ... !", irank);
            param->iter_solve = 0;
            param->dtg = param->dt_min;
            // param->dt = param->dtg;
            groundwater_mat_coeff(data, gmap, param, irank);
            groundwater_rhs(data, gmap, param, irank);
            build_groundwater_system(*data, gmap, param, A, b);
            solve_groundwater_system(data, gmap, A, b, x, param);
        }
    }
    else {
        // if (param->dtg > 10.0*param->dt_min & param->dtg == param->dtn & param->dtg <= param->dt_max)   {
        //     mpi_print("   >> Switching to NEWTON ... !", irank);
        //     param->iter_solve = 1;
        //     newton_iter(data, gmap, A, b, x, param, irank);
        //     if (param->converge == 0)  {
        //         mpi_print("   >> Switching to PCA ... !", irank);
        //         param->iter_solve = 0;
        //         param->dtg = param->dt_min;
        //         // param->dt = param->dtg;
        //         groundwater_mat_coeff(data, gmap, param, irank);
        //         groundwater_rhs(data, gmap, param, irank);
        //         build_groundwater_system(*data, gmap, param, A, b);
        //         solve_groundwater_system(data, gmap, A, b, x, param);
        //     }
        // }
        // else {
        //     groundwater_mat_coeff(data, gmap, param, irank);
        //     groundwater_rhs(data, gmap, param, irank);
        //     build_groundwater_system(*data, gmap, param, A, b);
        //     solve_groundwater_system(data, gmap, A, b, x, param);
        // }
        groundwater_mat_coeff(data, gmap, param, irank);
        groundwater_rhs(data, gmap, param, irank);
        build_groundwater_system(*data, gmap, param, A, b);
        solve_groundwater_system(data, gmap, A, b, x, param);
    }
    enforce_head_bc(data, gmap, param);
    if (param->use_mpi == 1)
    {mpi_exchange_subsurf((*data)->h, gmap, 2, param, irank, nrank);}

    // >>> Corrector step
    compute_K_face(data, gmap, param, irank, nrank);
    groundwater_flux(data, gmap, param, irank);
    if (param->use_corrector == 1 & param->iter_solve == 0)
    {
        check_room(data, gmap, param);
        update_water_content(data, gmap, param);
    }
    else
    {
        for (ii = 0; ii < param->n3ci; ii++)
        {(*data)->wc[ii] = compute_wch(*data, (*data)->h[ii], ii, param);}
    }
    for (ii = 0; ii < param->n3ci; ii++)
    {
        (*data)->hwc[ii] = compute_hwc(*data, ii, param);
        (*data)->wch[ii] = compute_wch(*data, (*data)->h[ii], ii, param);
    }
    if (param->use_mpi == 1)
    {
        mpi_exchange_subsurf((*data)->wc, gmap, 2, param, irank, nrank);
        mpi_exchange_subsurf((*data)->wch, gmap, 2, param, irank, nrank);
        mpi_exchange_subsurf((*data)->hwc, gmap, 2, param, irank, nrank);
    }

    // >>> Post-allocation step
    if (param->use_corrector == 1 & param->iter_solve == 0)
    {reallocate_water_content(data, gmap, param, irank);}

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
    if (param->sync_coupling == 1 | param->sim_shallowwater == 0)  {param->dt = param->dtg;}
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
    double hp, hm, ksp, ksm;
    // conductivities for interior cells
    for (ii = 0; ii < param->n3ci; ii++)
    {
        // Kx
        hp = (*data)->h[gmap->iPjckc[ii]];  ksp = (*data)->Ksx[gmap->iPjckc[ii]];
        hm = (*data)->h[ii];  ksm = (*data)->Ksx[ii];
        Kp = compute_K(*data, hp, ksp, gmap->iPjckc[ii], param);
        Km = compute_K(*data, hm, ksm, ii, param);
        (*data)->Kx[ii] = 0.5 * (Kp + Km);
        // Ky
        hp = (*data)->h[gmap->icjPkc[ii]];  ksp = (*data)->Ksy[gmap->icjPkc[ii]];
        hm = (*data)->h[ii];  ksm = (*data)->Ksy[ii];
        Kp = compute_K(*data, hp, ksp, gmap->icjPkc[ii], param);
        Km = compute_K(*data, hm, ksm, ii, param);
        (*data)->Ky[ii] = 0.5 * (Kp + Km);
        if (gmap->jj[ii] == param->ny-1 & irank >= param->mpi_nx*(param->mpi_ny-1))
        {
            if (param->bctype_GW[2] == 0)   {(*data)->Ky[ii] = 0.0;}
            else if (param->bctype_GW[2] == 1)  {(*data)->Ky[ii] = Km;}
        }
        // Kz
        hp = (*data)->h[gmap->icjckP[ii]];  ksp = (*data)->Ksz[gmap->icjckP[ii]];
        hm = (*data)->h[ii];  ksm = (*data)->Ksz[ii];
        Kp = compute_K(*data, hp, ksp, gmap->icjckP[ii], param);
        Km = compute_K(*data, hm, ksm, ii, param);
        if (gmap->istop[gmap->icjckP[ii]] == 1) {(*data)->Kz[ii] = Kp;}
        else if (gmap->actv[ii] == 0)   {(*data)->Kz[ii] = 0.0;}
        else if (gmap->icjckP[ii] > param->n3ci)    {(*data)->Kz[ii] = Km;}
        else    {(*data)->Kz[ii] = 0.5 * (Kp + Km);}
    }
    // conductivities on iM, jM, kM faces
    for (ii = 0; ii < param->ny*param->nz; ii++)
    {
        hp = (*data)->h[gmap->iMin[ii]];  ksp = (*data)->Ksx[gmap->iMin[ii]];
        hm = (*data)->h[gmap->iMou[ii]];  ksm = (*data)->Ksx[gmap->iMou[ii]];
        Kp = compute_K(*data, hp, ksp, gmap->iMin[ii], param);
        // Km = compute_K(*data, hm, ksm, gmap->iMou[ii], param);
        Km = Kp;
        (*data)->Kx[gmap->iMou[ii]] = 0.5 * (Kp + Km);
        if (param->bctype_GW[0] == 0 & (irank+1) % param->mpi_nx == 0)
        {(*data)->Kx[gmap->iPin[ii]] = 0;}
        if (param->bctype_GW[1] == 0 & irank % param->mpi_nx == 0)
        {(*data)->Kx[gmap->iMou[ii]] = 0;}
    }
    for (ii = 0; ii < param->nx*param->nz; ii++)
    {
        hp = (*data)->h[gmap->jMin[ii]];  ksp = (*data)->Ksy[gmap->jMin[ii]];
        hm = (*data)->h[gmap->jMou[ii]];  ksm = (*data)->Ksy[gmap->jMou[ii]];
        Kp = compute_K(*data, hp, ksp, gmap->jMin[ii], param);
        // Km = compute_K(*data, hm, ksm, gmap->jMou[ii], param);
        Km = Kp;
        (*data)->Ky[gmap->jMou[ii]] = 0.5 * (Kp + Km);
        if (param->bctype_GW[3] == 0 & irank < param->mpi_nx)
        {(*data)->Ky[gmap->jMou[ii]] = 0;}
        else if (param->bctype_GW[3] == 1 & irank < param->mpi_nx)
        {(*data)->Ky[gmap->jMou[ii]] = Km;}
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
            hp = (*data)->h[ii];  ksp = (*data)->Ksz[ii];
            Kp = compute_K(*data, hp, ksp, ii, param);
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

    // seal inactive faces
    for (ii = 0; ii < param->n3ci; ii++)    {
        if (gmap->actv[ii] == 0)    {
            (*data)->Kx[ii] = 0.0;  (*data)->Kx[gmap->iMjckc[ii]] = 0.0;
            (*data)->Ky[ii] = 0.0;  (*data)->Ky[gmap->icjMkc[ii]] = 0.0;
            (*data)->Kz[ii] = 0.0;  (*data)->Kz[gmap->icjckM[ii]] = 0.0;
        }
    }

    // MAXP1
    // for (ii = 0; ii < param->n3ci; ii++)
    // {
    //     if (gmap->jj[ii] > 4)  {(*data)->Ky[ii] = 0.0; (*data)->Kz[ii] = 0.0; (*data)->Ksy[ii] = 0.0; (*data)->Ksz[ii] = 0.0;}
    //     // else if (gmap->jj[ii] == 0) {(*data)->Ky[ii] = 0.0;}
    // }

    // vcatch
    // for (ii = 0; ii < param->n3ci; ii++)
    // {
    //     if (gmap->jj[ii] >= 19)  {
    //         (*data)->Kx[ii] = 0.0; (*data)->Ksx[ii] = 0.0;
    //         (*data)->Ky[ii] = 0.0; (*data)->Kz[gmap->icjckM[ii]] = 0.0;
    //         (*data)->Ksy[ii] = 0.0; (*data)->Ksz[gmap->icjckM[ii]] = 0.0;
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


// >>>>> Build the Jacobian coefficients
void jacobian_mat_coeff(Data **data, Map *gmap, Config *param, int irank)
{
    int ii, jp, jm;
    double dzf;
    for (ii = 0; ii < param->n3ci; ii++)
    {
        if (gmap->actv[ii] == 1)    {
            (*data)->resi[ii] = compute_residual(data, gmap, ii, "none", param, irank);
            (*data)->Grhs[ii] = -(*data)->resi[ii];
            (*data)->Gxp[ii] = compute_jacobian_fd(data, gmap, ii, "xp", param, irank);
            (*data)->Gxm[ii] = compute_jacobian_fd(data, gmap, ii, "xm", param, irank);
            (*data)->Gyp[ii] = compute_jacobian_fd(data, gmap, ii, "yp", param, irank);
            (*data)->Gym[ii] = compute_jacobian_fd(data, gmap, ii, "ym", param, irank);
            (*data)->Gzp[ii] = compute_jacobian_fd(data, gmap, ii, "zp", param, irank);
            (*data)->Gzm[ii] = compute_jacobian_fd(data, gmap, ii, "zm", param, irank);
            (*data)->Gct[ii] = compute_jacobian_fd(data, gmap, ii, "ct", param, irank);
        }
        else {
            (*data)->Gxp[ii] = 0.0; (*data)->Gxm[ii] = 0.0;
            (*data)->Gyp[ii] = 0.0; (*data)->Gym[ii] = 0.0;
            (*data)->Gzp[ii] = 0.0; (*data)->Gzm[ii] = 0.0;
            (*data)->Grhs[ii] = 0.0;
            (*data)->Gct[ii] = 1.0;
        }
    }
}

// >>>>> The iterative Newton loop
void newton_iter(Data **data, Map *gmap, QMatrix A, Vector b, Vector x, Config *param, int irank) {
    int ii, iter, iter_max = 15;
    double eps, epsn, eps_min = 1e-8;
    for (ii = 0; ii < param->n3ci; ii++)    {
        (*data)->hnm[ii] = (*data)->hn[ii];
        (*data)->hn[ii] = (*data)->h[ii];
        (*data)->hp[ii] = (*data)->h[ii];
        (*data)->wcn[ii] = (*data)->wc[ii];
    }

    iter = 1;
    eps = 1.0;
    epsn = eps;
    param->converge = 1;
    param->n_iter = 0;
    while (iter < iter_max & sqrt(eps) > 1e-5 + 1e-5*sqrt(eps_min)) {
        // build Jacobian Matrix
        jacobian_mat_coeff(data, gmap, param, irank);
        build_groundwater_system(*data, gmap, param, A, b);
        solve_groundwater_system(data, gmap, A, b, x, param);
        if (iter == 1)  {
            for (ii = 0; ii < param->n3ci; ii++)    {
                if (gmap->actv[ii] == 1)    {(*data)->h[ii] += (*data)->h_incr[ii];}
                (*data)->hp[ii] = (*data)->h[ii];
            }
        }
        else {
            eps = 0.0;
            eps_min = 0.0;
            for (ii = 0; ii < param->n3ci; ii++)    {
                if (gmap->actv[ii] == 1)    {
                    (*data)->h[ii] += (*data)->h_incr[ii];
                    eps += pow(((*data)->h[ii] - (*data)->hp[ii]), 2.0);
                    eps_min += pow((*data)->h[ii], 2.0);
                }
                (*data)->hp[ii] = (*data)->h[ii];
            }
            if (fabs(eps-epsn) < 1e-6)    {eps = 0.0;}
        }
        printf("     >>> Newton iter : %d    EPS : %f<<< \n",iter,eps);
        param->n_iter = iter;
        if (iter == iter_max-1 | isnan(eps) | eps > 1e4)  {
            printf("     >>> No converge of Newton \n");
            param->converge = 0;
            for (ii = 0; ii < param->n3ci; ii++)    {
                (*data)->h[ii] = (*data)->hn[ii];
            }
            break;
        }
        if (iter > 1)   {epsn = eps;}
        iter += 1;
    }

}


// >>>>> Compute matrix coefficients <<<<<
void groundwater_mat_coeff(Data **data, Map *gmap, Config *param, int irank)
{
    int ii, jp, jm;
    double dzf;
    for (ii = 0; ii < param->n3ci; ii++)
    {
        if (gmap->actv[ii] == 1)    {
            // coeff xp
            (*data)->Gxp[ii] = - (*data)->Kx[ii] * param->dtg * (*data)->r_rhoxp[ii] * (*data)->r_viscxp[ii] * gmap->cosx[ii] / (pow(param->dx,2.0));
            (*data)->Gxp[ii] = (*data)->Gxp[ii] * gmap->Ax[ii] * param->dx;
            // coeff xm
            (*data)->Gxm[ii] = - (*data)->Kx[gmap->iMjckc[ii]] * param->dtg * (*data)->r_rhoxp[gmap->iMjckc[ii]] * (*data)->r_viscxp[gmap->iMjckc[ii]] * gmap->cosx[gmap->iMjckc[ii]] / (pow(param->dx,2.0));
            (*data)->Gxm[ii] = (*data)->Gxm[ii] * gmap->Ax[gmap->iMjckc[ii]] * param->dx;
            // coeff yp
            (*data)->Gyp[ii] = - (*data)->Ky[ii] * param->dtg * (*data)->r_rhoyp[ii] * (*data)->r_viscyp[ii] * gmap->cosy[ii] / (pow(param->dy,2.0));
            (*data)->Gyp[ii] = (*data)->Gyp[ii] * gmap->Ay[ii] * param->dy;
            if (gmap->jj[ii] == param->ny-1 & param->bctype_GW[2] == 1 & irank >= param->mpi_nx*(param->mpi_ny-1))     {(*data)->Gyp[ii] = 2.0 * (*data)->Gyp[ii]; }
            // coeff ym
            (*data)->Gym[ii] = - (*data)->Ky[gmap->icjMkc[ii]] * param->dtg * (*data)->r_rhoyp[gmap->icjMkc[ii]] * (*data)->r_viscyp[gmap->icjMkc[ii]] * gmap->cosy[gmap->icjMkc[ii]] / (pow(param->dy,2.0));
            (*data)->Gym[ii] = (*data)->Gym[ii] * gmap->Ay[gmap->icjMkc[ii]] * param->dy;
            if (gmap->jj[ii] == 0 & param->bctype_GW[3] == 1 & irank < param->mpi_nx)     {(*data)->Gym[ii] = 2.0 * (*data)->Gym[ii];}
            // coeff zp
            dzf = 0.5 * (gmap->dz3d[ii] + gmap->dz3d[gmap->icjckP[ii]]);
            (*data)->Gzp[ii] = - (*data)->Kz[ii] * param->dtg * (*data)->r_rhozp[ii] * (*data)->r_visczp[ii] / (gmap->dz3d[ii]*dzf);
            (*data)->Gzp[ii] = (*data)->Gzp[ii] * gmap->V[ii];
            // coeff zm
            if (gmap->istop[ii] == 1)   {dzf = 0.5 * gmap->dz3d[ii];}
            else    {dzf = 0.5 * (gmap->dz3d[ii] + gmap->dz3d[gmap->icjckM[ii]]);}
            (*data)->Gzm[ii] = - (*data)->Kz[gmap->icjckM[ii]] * param->dtg * (*data)->r_rhozp[gmap->icjckM[ii]] * (*data)->r_visczp[gmap->icjckM[ii]] / (gmap->dz3d[ii]*dzf);
            (*data)->Gzm[ii] = (*data)->Gzm[ii] * gmap->V[ii];

            // coeff ct
            (*data)->Gct[ii] = ((*data)->ch[ii] + param->Ss*(*data)->wcn[ii]/(*data)->wcs[ii]) * (*data)->r_rho[ii];
            (*data)->Gct[ii] = (*data)->Gct[ii] * gmap->V[ii];
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
        else {
            (*data)->Gxp[ii] = 0.0; (*data)->Gxm[ii] = 0.0;
            (*data)->Gyp[ii] = 0.0; (*data)->Gym[ii] = 0.0;
            (*data)->Gzp[ii] = 0.0; (*data)->Gzm[ii] = 0.0;
            (*data)->Gct[ii] = 1.0;
        }

    }
}

// >>>>> Compute right hand side <<<<<
void groundwater_rhs(Data **data, Map *gmap, Config *param, int irank)
{
    int ii, sign;
    double dh;
    for (ii = 0; ii < param->n3ci; ii++)
    {
        if (gmap->actv[ii] == 1)    {
            // base terms
            (*data)->Grhs[ii] = ((*data)->ch[ii] + param->Ss*(*data)->wcn[ii]/(*data)->wcs[ii]) * (*data)->hn[ii] * (*data)->r_rho[ii] * gmap->V[ii];

            (*data)->Grhs[ii] -= param->dtg * gmap->V[ii] *
                (*data)->Kz[ii] * (*data)->r_rhozp[ii] * (*data)->r_rhozp[ii] * (*data)->r_visczp[ii] / gmap->dz3d[ii];

            (*data)->Grhs[ii] += param->dtg * gmap->V[ii] *
                (*data)->Kz[gmap->icjckM[ii]] * (*data)->r_rhozp[gmap->icjckM[ii]] * (*data)->r_rhozp[gmap->icjckM[ii]] * (*data)->r_visczp[gmap->icjckM[ii]] / gmap->dz3d[ii];
            // terrain terms
            if (param->follow_terrain == 1) {

                if (gmap->bot3d[gmap->iPjckc[ii]] > gmap->bot3d[ii])    {sign = 1;}   else    {sign = -1;}
                (*data)->Grhs[ii] += sign * param->dtg * gmap->Ax[ii] * param->dx * (*data)->Kx[ii] * gmap->sinx[ii] / param->dx;
                if (gmap->bot3d[ii] > gmap->bot3d[gmap->iMjckc[ii]])    {sign = -1;}   else    {sign = 1;}
                (*data)->Grhs[ii] += sign * param->dtg * gmap->Ax[gmap->iMjckc[ii]] * param->dx * (*data)->Kx[gmap->iMjckc[ii]] * gmap->sinx[gmap->iMjckc[ii]] / param->dx;
                if (gmap->bot3d[gmap->icjPkc[ii]] > gmap->bot3d[ii])    {sign = 1;}   else    {sign = -1;}
                (*data)->Grhs[ii] += sign * param->dtg * gmap->Ay[ii] * param->dy *(*data)->Ky[ii] * gmap->siny[ii] / param->dy;
                if (gmap->bot3d[ii] > gmap->bot3d[gmap->icjMkc[ii]])    {sign = -1;}   else    {sign = 1;}
                (*data)->Grhs[ii] += sign * param->dtg * gmap->Ay[gmap->icjMkc[ii]] * param->dy * (*data)->Ky[gmap->icjMkc[ii]] * gmap->siny[gmap->icjMkc[ii]] / param->dy;
            }
            // density term
            // (*data)->Grhs[ii] -= (*data)->wc[ii] * ((*data)->r_rho[ii] - (*data)->r_rhon[ii]);
            // boundary terms
            // x
            if (gmap->ii[ii] == param->nx-1)
            {(*data)->Grhs[ii] -= (*data)->Gxp[ii] * (*data)->hn[gmap->iPjckc[ii]];}
            if (gmap->ii[ii] == 0)
            {(*data)->Grhs[ii] -= (*data)->Gxm[ii] * (*data)->hn[gmap->iMjckc[ii]];}
            // y
            if (gmap->jj[ii] == (param->ny-1))
            {
                // flux bc
                if (param->bctype_GW[2] == 2 & irank >= param->mpi_nx*(param->mpi_ny-1))
                {(*data)->Grhs[ii] -= (*data)->r_rhoyp[ii] * gmap->Ay[ii] * param->dy * param->qyp * param->dtg / param->dy;}
                else
                {(*data)->Grhs[ii] -= (*data)->Gyp[ii] * (*data)->hn[gmap->icjPkc[ii]];}
            }
            if (gmap->jj[ii] == 0)
            {

                if (param->bctype_GW[3] == 2 & irank < param->mpi_nx)
                {
                    (*data)->Grhs[ii] -= (*data)->r_rhoyp[gmap->icjMkc[ii]] * gmap->Ay[gmap->icjMkc[ii]] * param->dy * param->qym * param->dtg / param->dy;
                }
                else
                {(*data)->Grhs[ii] -= (*data)->Gym[ii] * (*data)->hn[gmap->icjMkc[ii]];}

            }
            // z
            if (gmap->kk[ii] == param->nz-1)
            {
                if (param->bctype_GW[4] == 2)
                {(*data)->Grhs[ii] += (*data)->r_rhozp[ii] * gmap->V[ii] * ((*data)->qbot + param->dtg * (*data)->Kz[ii] * (*data)->r_rhozp[ii] * (*data)->r_visczp[ii]/ gmap->dz3d[ii]);}
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
                                (*data)->Grhs[ii] += - param->dtg * (*data)->r_rhozp[gmap->icjckM[ii]] * gmap->V[ii] * \
                                    ((*data)->qtop[gmap->top2d[ii]] + (*data)->Kz[gmap->icjckM[ii]]*(*data)->r_rhozp[gmap->icjckM[ii]]*(*data)->r_visczp[gmap->icjckM[ii]]) / gmap->dz3d[ii];
                            }
                            else if ((*data)->qtop[gmap->top2d[ii]] < 0.0 & (*data)->wc[ii] > (*data)->wcr[ii])
                            {
                                (*data)->Grhs[ii] += - param->dtg * (*data)->r_rhozp[gmap->icjckM[ii]] * gmap->V[ii] * \
                                    ((*data)->qtop[gmap->top2d[ii]] + (*data)->Kz[gmap->icjckM[ii]]*(*data)->r_rhozp[gmap->icjckM[ii]]*(*data)->r_visczp[gmap->icjckM[ii]]) / gmap->dz3d[ii];
                            }
                            else
                            {
                                (*data)->Grhs[ii] += - param->dtg * (*data)->r_rhozp[gmap->icjckM[ii]] * gmap->V[ii] * \
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
                        (*data)->Grhs[ii] += - param->dtg * (*data)->r_rhozp[gmap->icjckM[ii]] * gmap->V[ii] * \
                            ((*data)->qtop[gmap->top2d[ii]] + (*data)->Kz[gmap->icjckM[ii]] * (*data)->r_rhozp[gmap->icjckM[ii]] * (*data)->r_visczp[gmap->icjckM[ii]]) / gmap->dz3d[ii];
                    }
                    else if (param->bctype_GW[5] == 1)
                    {(*data)->Grhs[ii] -= (*data)->Gzm[ii] * (*data)->htop;}
                }
            }
        }
        else    {
            (*data)->Grhs[ii] = (*data)->hn[ii];
        }

    }
}


// >>>>> Build linear system <<<<<
void build_groundwater_system(Data *data, Map *gmap, Config *param, QMatrix A, Vector b)
{
    size_t ii, jj, kk, n, irow;
    int im, dist;
    double value, *row;
    row = malloc(8*sizeof(double));
    for (irow = 1; irow <= gmap->nactv; irow++) {
        ii = gmap->mat2dom[irow-1];
        kk = 0;
        n = 7;
        if (gmap->actv[gmap->iMjckc[ii]] == 0 | gmap->ii[ii] == 0)  {n -= 1;}
        if (gmap->actv[gmap->iPjckc[ii]] == 0 | gmap->ii[ii] == param->nx-1)    {n -= 1;}
        if (gmap->actv[gmap->icjMkc[ii]] == 0 | gmap->jj[ii] == 0)  {n -= 1;}
        if (gmap->actv[gmap->icjPkc[ii]] == 0 | gmap->jj[ii] == param->ny-1)    {n -= 1;}
        if (gmap->actv[gmap->icjckM[ii]] == 0 | gmap->istop[ii] == 1)   {n -= 1;}
        if (gmap->actv[gmap->icjckP[ii]] == 0 | gmap->kk[ii] == param->nz-1)    {n -= 1;}
        Q_SetLen(&A, irow, n);
        // set ym entry
        dist = gmap->distym[irow-1];
        if (gmap->jj[ii] > 0 & dist > 0)    {Q_SetEntry(&A, irow, kk, irow-dist, data->Gym[ii]);    kk++;}
        // // set xm entry
        dist = gmap->distxm[irow-1];
        if (gmap->ii[ii] > 0 & dist > 0)    {Q_SetEntry(&A, irow, kk, irow-dist, data->Gxm[ii]);    kk++;}
        // // set zm entry
        dist = gmap->distzm[irow-1];
        if (gmap->kk[ii] > 0 & gmap->istop[ii] == 0)    {Q_SetEntry(&A, irow, kk, irow-dist, data->Gzm[ii]);    kk++;}
        // ct
        Q_SetEntry(&A, irow, kk, irow, data->Gct[ii]);  kk++;
        // set zp entry
        dist = gmap->distzp[irow-1];
        if (gmap->kk[ii] < param->nz-1)    {Q_SetEntry(&A, irow, kk, irow+dist, data->Gzp[ii]);    kk++;}
        // // set xp entry
        dist = gmap->distxp[irow-1];
        if (gmap->ii[ii] < param->nx-1 & dist > 0)    {Q_SetEntry(&A, irow, kk, irow+dist, data->Gxp[ii]);    kk++;}
        // // set yp entry
        dist = gmap->distyp[irow-1];
        if (gmap->jj[ii] < param->ny-1 & dist > 0)    {Q_SetEntry(&A, irow, kk, irow+dist, data->Gyp[ii]);    kk++;}
        // set rhs
        V_SetCmp(&b, irow, data->Grhs[ii]);
    }
}

// >>>>> solve linear system <<<<<
void solve_groundwater_system(Data **data, Map *gmap, QMatrix A, Vector b, Vector x, Config *param)
{
    size_t ii, irow;
    V_SetAllCmp(&x, 0.0);
    SetRTCAccuracy(0.00000001);
    if (param->iter_solve == 0) {
        CGIter(&A, &x, &b, 1000, SSORPrecond, 1);
        for (irow = 0; irow < gmap->nactv; irow++) {
            ii = gmap->mat2dom[irow];
            (*data)->h[ii] = V_GetCmp(&x, irow+1);
        }
    }
    else {
        GMRESIter(&A, &x, &b, 1000, ILUPrecond, 1);
        for (irow = 0; irow < gmap->nactv; irow++) {
            ii = gmap->mat2dom[irow];
            (*data)->h_incr[ii] = V_GetCmp(&x, irow+1);
        }
    }
}

// >>>>> enforce head boundary conditions
void enforce_head_bc(Data **data, Map *gmap, Config *param)
{
    int ii;
    // for (ii = 0; ii < param->n3ci; ii++)    {
    //     if (gmap->actv[ii] == 1 & (*data)->h[ii] < -1500.0) {(*data)->h[ii] = -1500.0;}
    // }
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
        if (gmap->actv[ii] == 1)    {
            // x flux
            (*data)->qx[ii] = darcy_flux(*data, gmap, ii, 0.0, 0.0, "x", "none", param, irank);
            // y flux
            (*data)->qy[ii] = darcy_flux(*data, gmap, ii, 0.0, 0.0, "y", "none", param, irank);
            // z flux
            (*data)->qz[ii] = darcy_flux(*data, gmap, ii, 0.0, 0.0, "z", "none", param, irank);
        }
    }
    for (ii = 0; ii < param->ny*param->nz; ii++)
    {(*data)->qx[gmap->iMou[ii]] = darcy_flux(*data, gmap, ii, 0.0, 0.0, "x", "back", param, irank);}
    for (ii = 0; ii < param->nx*param->nz; ii++)
    {(*data)->qy[gmap->jMou[ii]] = darcy_flux(*data, gmap, ii, 0.0, 0.0, "y", "back", param, irank);}
    for (jj = 0; jj < param->n3ci; jj++)
    {
        if (gmap->actv[jj] == 1 & gmap->istop[jj] == 1)
        {
            ii = gmap->icjckM[jj];
            (*data)->qz[ii] = darcy_flux(*data, gmap, jj, 0.0, 0.0, "z", "back", param, irank);
            // dzf = 0.5 * gmap->dz3d[jj];
            if (param->sim_shallowwater == 1)
            {
                if ((*data)->qz[ii] < 0.0)  {
                    if ((*data)->wcs[jj] - (*data)->wc[jj] < fabs((*data)->qz[ii])*param->dtg/gmap->V[ii]) {
                        if ((*data)->dept[gmap->top2d[jj]] <= param->min_dept)  {
                            (*data)->eta[gmap->top2d[jj]] += fabs((*data)->qz[ii])*param->dtg/gmap->Az[ii];
                            (*data)->dept[gmap->top2d[jj]] += fabs((*data)->qz[ii])*param->dtg/gmap->Az[ii];
                        }
                    }
                }
                (*data)->qseepage_old[gmap->top2d[jj]] = (*data)->qseepage[gmap->top2d[jj]];
                if ((*data)->reset_seepage[gmap->top2d[jj]] == 1)   {
                    (*data)->qseepage[gmap->top2d[jj]] = 0.0;
                    (*data)->reset_seepage[gmap->top2d[jj]] = 0;
                }
                (*data)->qseepage[gmap->top2d[jj]] += (*data)->qz[ii]/gmap->Az[ii];
                // if evaporation exists, evaporation does not contribute to seepage
                if (param->bctype_GW[5] == 2 & (*data)->qtop[gmap->top2d[jj]] > 0.0
                    & (*data)->wc[jj] > (*data)->wcr[jj] & (*data)->dept[gmap->top2d[jj]] <= param->min_dept)
                {
                    (*data)->qseepage[gmap->top2d[jj]] -= (*data)->qtop[gmap->top2d[jj]];
                }
                (*data)->qss[gmap->top2d[jj]] = (*data)->qseepage[gmap->top2d[jj]];

            }
        }
    }

}

// approximate seepage when dts != dtg
void pseudo_seepage(Data **data, Map *gmap, Config *param, double dt)  {
    int ii, jj;
    double q_max = 2e-4, delt = param->dtg, delt_new;
    for (ii = 0; ii < param->n3ci; ii++)    {
        if (gmap->istop[ii] == 1)   {
            jj = gmap->top2d[ii];
            // (*data)->qss[jj] = (*data)->qseepage_old[jj] + dt * ((*data)->qseepage[jj] - (*data)->qseepage_old[jj]);

            if ((*data)->dept[jj] > 0.0 & (*data)->eta[jj] != (*data)->etan[jj])    {
                delt_new = fabs(q_max * gmap->dz3d[ii] * param->dt / (*data)->Kz[gmap->icjckM[ii]] / ((*data)->eta[jj]-(*data)->etan[jj]));
                if (delt_new < delt)  {delt = delt_new;}
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
    if (param->iter_solve == 0) {
        for (ii = 0; ii < param->n3ci; ii++)
        {
            if (gmap->actv[ii] == 1)    {
                // update water content
                coeff = gmap->V[ii] * ((*data)->r_rho[ii] + (*data)->r_rho[ii] * param->Ss * ((*data)->h[ii]-(*data)->hn[ii]) / (*data)->wcs[ii]);
                // x component
                dqx = param->dtg * ((*data)->qx[ii]*(*data)->r_rhoxp[ii] - (*data)->qx[gmap->iMjckc[ii]]*(*data)->r_rhoxp[gmap->iMjckc[ii]]);

                // y component
                dqy = param->dtg * ((*data)->qy[ii]*(*data)->r_rhoyp[ii] - (*data)->qy[gmap->icjMkc[ii]]*(*data)->r_rhoyp[gmap->icjMkc[ii]]);
                // z component
                dqz = param->dtg * ((*data)->qz[ii]*(*data)->r_rhozp[ii] - (*data)->qz[gmap->icjckM[ii]]*(*data)->r_rhozp[gmap->icjckM[ii]]);

                (*data)->wc[ii] = ((*data)->wcn[ii]*(*data)->r_rho[ii]*gmap->V[ii] + dqx + dqy + dqz) / coeff;
            }
            // Warrick 1971
            // if (ii == 0)    {(*data)->wc[ii] = (*data)->wcs[ii];}
        }
    }
    else {
        for (ii = 0; ii < param->n3ci; ii++)    {(*data)->wc[ii] = compute_wch(*data, (*data)->h[ii], ii, param);}
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
            (*data)->wch[ii] = compute_wch(*data, (*data)->h[ii], ii, param);
            (*data)->hwc[ii] = compute_hwc(*data, ii, param);
            // over-saturated cell
            if ((*data)->wc[ii] >= (*data)->wcs[ii])
            {
                // send moisture
                if (param->post_allocate == 1)
                {
                    dV = ((*data)->wc[ii] - (*data)->wcs[ii]) * gmap->dz3d[ii] * param->dx * param->dy;
                    if (dV > 0.1 * (*data)->wcs[ii]*gmap->dz3d[ii]*param->dx*param->dy)
                    {(*data)->repeat[0] = 1;}
                    check_head_gradient(data, gmap, param, ii);
                    dV = allocate_send(data, gmap, param, ii, dV);
                    // if (dV > 0)    {dV = allocate_send(data, gmap, param, ii, dV);}
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
                            // dV = allocate_recv(data, gmap, param, ii, dV);
                            // if (dV > 0) {dV = allocate_recv(data, gmap, param, ii, dV);}
                            alloc_type3 += 1;
                        }
                        // send moisture
                        else
                        {
                            dV = ((*data)->wc[ii] - (*data)->wch[ii]) * gmap->dz3d[ii] * param->dx * param->dy;
                            if (dV > 0.1 * (*data)->wcs[ii]*gmap->dz3d[ii]*param->dx*param->dy)
                            {
                                (*data)->repeat[0] = 1;
                            }
                            check_head_gradient(data, gmap, param, ii);
                            // dV = allocate_send(data, gmap, param, ii, dV);
                            // if (dV > 0)    {dV = allocate_send(data, gmap, param, ii, dV);}
                            alloc_type4 += 1;

                        }
                    }

                    dV_tot += dV;
                    (*data)->wc[ii] = (*data)->wch[ii];
                    (*data)->room[ii] = ((*data)->wcs[ii] - (*data)->wc[ii]) * param->dx*param->dy*gmap->dz3d[ii];
                }
            }
        }
    }

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
                    (*data)->qz[ll] += dVzm * gmap->Az[ll] / (param->dx * param->dy * param->dtg);
                    dVzm = 0.0;
                    break;
                }
                else
                {
                    (*data)->wc[ll] += (*data)->room[ll] / (param->dx*param->dy*gmap->dz3d[ll]);
                    dVzm -= (*data)->room[ll];
                    (*data)->qz[ll] += (*data)->room[ll] * gmap->Az[ll] / (param->dx * param->dy * param->dtg);
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
                    (*data)->qz[gmap->icjckM[ll]] += dVzm * gmap->Az[gmap->icjckM[ll]] / (param->dx * param->dy * param->dtg);
                }
                else
                {
                    if (dVzm / (param->dx*param->dy) > param->min_dept)
                    {
                        (*data)->dept[gmap->top2d[ll]] += dVzm / (param->dx*param->dy);
                        (*data)->eta[gmap->top2d[ll]] += dVzm / (param->dx*param->dy);
                        (*data)->qz[gmap->icjckM[ll]] += dVzm * gmap->Az[gmap->icjckM[ll]] / (param->dx * param->dy * param->dtg);
                    }
                }
                dVzm = 0.0;
            }
            else
            {
                if (param->bctype_GW[5] == 1)
                {
                    (*data)->qz[gmap->icjckM[ll]] += dVzm * gmap->Az[ll] / (param->dx * param->dy * param->dtg);
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
                                (*data)->qz[gmap->icjckM[ll]] -= dVzm * gmap->Az[gmap->icjckM[ll]] / (param->dx * param->dy * param->dtg);
                                (*data)->room[ll] -= dVzm;
                                dVzm = 0.0;
                                break;
                            }
                            else
                            {
                                (*data)->wc[ll] += (*data)->room[ll] / (param->dx*param->dy*gmap->dz3d[ll]);
                                (*data)->qz[gmap->icjckM[ll]] -= (*data)->room[ll] * gmap->Az[gmap->icjckM[ll]] / (param->dx * param->dy * param->dtg);
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
                    (*data)->qz[gmap->icjckM[ll]] -= dVzp * gmap->Az[gmap->icjckM[ll]] / (param->dx * param->dy * param->dtg);
                    (*data)->room[ll] -= dVzp;
                    dVzp = 0.0;
                    break;
                }
                else
                {
                    (*data)->wc[ll] += (*data)->room[ll] / (param->dx*param->dy*gmap->dz3d[ll]);
                    (*data)->qz[gmap->icjckM[ll]] -= (*data)->room[ll] * gmap->Az[gmap->icjckM[ll]] / (param->dx * param->dy * param->dtg);
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
                (*data)->qx[gmap->iMjckc[ll]] -= dVxp * gmap->Ax[gmap->iMjckc[ll]] / (gmap->dz3d[ll] * param->dy * param->dtg);
                (*data)->room[ll] -= dVxp;
                dVxp = 0.0;
            }
            else
            {
                (*data)->wc[ll] += (*data)->room[ll] / (param->dx*param->dy*gmap->dz3d[ll]);
                (*data)->qx[gmap->iMjckc[ll]] -= (*data)->room[ll] * gmap->Ax[gmap->iMjckc[ll]] / (gmap->dz3d[ll] * param->dy * param->dtg);
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
                (*data)->qx[ll] += dVxm * gmap->Ax[ll] / (gmap->dz3d[ll] * param->dy * param->dtg);
                (*data)->room[ll] -= dVxm;
                dVxm = 0.0;
            }
            else
            {
                (*data)->wc[ll] += (*data)->room[ll] / (param->dx*param->dy*gmap->dz3d[ll]);
                (*data)->qx[ll] += (*data)->room[ll] * gmap->Ax[ll] / (gmap->dz3d[ll] * param->dy * param->dtg);
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
                (*data)->qy[gmap->icjMkc[ll]] -= dVyp * gmap->Ay[gmap->icjMkc[ll]] / (gmap->dz3d[ll] * param->dx * param->dtg);
                (*data)->room[ll] -= dVyp;
                dVyp = 0.0;
            }
            else
            {
                (*data)->wc[ll] += (*data)->room[ll] / (param->dx*param->dy*gmap->dz3d[ll]);
                (*data)->qy[gmap->icjMkc[ll]] -= (*data)->room[ll] * gmap->Ay[gmap->icjMkc[ll]] / (gmap->dz3d[ll] * param->dx * param->dtg);
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
                (*data)->qy[ll] += dVym * gmap->Ay[ll] / (gmap->dz3d[ll] * param->dx * param->dtg);
                (*data)->room[ll] -= dVym;
                dVym = 0.0;
            }
            else
            {
                (*data)->wc[ll] += (*data)->room[ll] / (param->dx*param->dy*gmap->dz3d[ll]);
                (*data)->qy[ll] += (*data)->room[ll] * gmap->Ay[ll] / (gmap->dz3d[ll] * param->dx * param->dtg);
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
                    (*data)->qz[ll] -= dVzm / (param->dx * param->dy * param->dtg);
                    (*data)->room[ll] += dVzm;
                    dVzm = 0.0;
                    break;
                }
                else
                {
                    dVzm -= ((*data)->wc[ll]-(*data)->wcr[ll]) * (param->dx*param->dy*gmap->dz3d[ll]);
                    (*data)->qz[ll] -= ((*data)->wc[ll]-(*data)->wcr[ll]) / (param->dx * param->dy * param->dtg);
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
                    // (*data)->qz[gmap->icjckM[ll]] -= dVzm / (param->dx * param->dy * param->dtg);
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
                    (*data)->qz[gmap->icjckM[ll]] += dVzp / (param->dx * param->dy * param->dtg);
                    (*data)->room[ll] += dVzp;
                    dVzp = 0.0;
                    break;
                }
                else
                {
                    dVzp -= ((*data)->wc[ll]-(*data)->wcr[ll]) * (param->dx*param->dy*gmap->dz3d[ll]);
                    (*data)->qz[gmap->icjckM[ll]] += ((*data)->wc[ll]-(*data)->wcr[ll]) / (param->dx * param->dy * param->dtg);
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
    //             (*data)->qx[gmap->iMjckc[ll]] += dVxp / (gmap->dz3d[ll] * param->dy * param->dtg);
    //             (*data)->room[ll] += dVxp;
    //             dVxp = 0.0;
    //         }
    //         else
    //         {
    //             dVxp -= ((*data)->wc[ll]-(*data)->wcr[ll]) * (param->dx*param->dy*gmap->dz3d[ll]);
    //             (*data)->qx[gmap->iMjckc[ll]] += ((*data)->wc[ll]-(*data)->wcr[ll]) / (gmap->dz3d[ll] * param->dy * param->dtg);
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
    //             (*data)->qx[ll] -= dVxm / (gmap->dz3d[ll] * param->dy * param->dtg);
    //             (*data)->room[ll] += dVxm;
    //             dVxm = 0.0;
    //         }
    //         else
    //         {
    //             dVxm -= ((*data)->wc[ll]-(*data)->wcr[ll]) * (param->dx*param->dy*gmap->dz3d[ll]);
    //             (*data)->qx[ll] -= ((*data)->wc[ll]-(*data)->wcr[ll]) / (gmap->dz3d[ll] * param->dy * param->dtg);
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
                (*data)->qy[gmap->icjMkc[ll]] += dVyp / (gmap->dz3d[ll] * param->dx * param->dtg);
                (*data)->room[ll] += dVyp;
                dVyp = 0.0;
            }
            else
            {
                dVyp -= ((*data)->wc[ll]-(*data)->wcr[ll]) * (param->dx*param->dy*gmap->dz3d[ll]);
                (*data)->qy[gmap->icjMkc[ll]] += ((*data)->wc[ll]-(*data)->wcr[ll]) / (gmap->dz3d[ll] * param->dx * param->dtg);
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
                (*data)->qy[ll] -= dVym / (gmap->dz3d[ll] * param->dx * param->dtg);
                (*data)->room[ll] += dVym;
                dVym = 0.0;
            }
            else
            {
                dVym -= ((*data)->wc[ll]-(*data)->wcr[ll]) * (param->dx*param->dy*gmap->dz3d[ll]);
                (*data)->qy[ll] -= ((*data)->wc[ll]-(*data)->wcr[ll]) / (gmap->dz3d[ll] * param->dx * param->dtg);
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
      (*data)->Vgflux[ii] = (*data)->Vg[ii] + param->dtg *
        ((-(*data)->qx[gmap->iMjckc[ii]] + (*data)->qx[ii])  +
        (-(*data)->qy[gmap->icjMkc[ii]] + (*data)->qy[ii])  +
        (-(*data)->qz[gmap->icjckM[ii]] + (*data)->qz[ii]) );

    }
}




// >>>>> adjust dt for groundwater solver <<<<<
void adaptive_time_step(Data *data, Map *gmap, Config **param, int root, int irank)
{
    double qin, qou, r_red, r_inc, dq, dq_max, dKdwc, dt_Co, dt_Comin, sqerr;
    double err, err_max=0.0, err0=1e-3, errm=1e-9, s=0.9, rmin=0.1, rmax=4.0;
    double dt_eta, dt_sync=(*param)->dt_max, q_target=5e-4;
    int ii, jj, unsat = 0, sat = 0, itop=0;
    r_red = 0.75;
    r_inc = 1.25;
    dq_max = 0.0;
    dt_Comin = 1e8;
    // check if both sat and unsat zones exist
    for (ii = 0; ii < (*param)->n3ci; ii++)
    {if (data->wc[ii] >= (*param)->wcs) {sat = 1;  break;}}
    for (ii = 0; ii < (*param)->n3ci; ii++)
    {if (data->wc[ii] < (*param)->wcs) {unsat = 1;  break;}}
    // adjust dt
    (*param)->dtn = (*param)->dtg;
    for (ii = 0; ii < (*param)->n3ci; ii++)
    {
        err = 0.5*(*param)->dtg*fabs((data->h[ii]-data->hn[ii])/(*param)->dtg - (data->hn[ii]-data->hnm[ii])/(*param)->dtn);
        if (err > err_max)    {err_max = err;}
        if (gmap->actv[ii] == 1 & gmap->istop[ii] == 0)
        {
            qin = data->qx[ii]/gmap->Ax[ii] + data->qy[ii]/gmap->Ay[ii] + data->qz[ii]/gmap->Az[ii];
            qou = data->qx[gmap->iMjckc[ii]]/gmap->Ax[ii] + data->qy[gmap->icjMkc[ii]]/gmap->Ay[ii] + data->qz[gmap->icjckM[ii]]/gmap->Az[ii];
            dq = fabs(qin - qou) * (*param)->dtg / gmap->dz3d[ii];
            if (dq > dq_max)    {
                dq_max = dq;
            }
            if (data->wc[ii] < (*param)->wcs)
            {
                dKdwc = compute_dKdwc(data, data->Ksz, ii, *param);
                dt_Co = (*param)->Co_max * gmap->dz3d[ii] / dKdwc;
                if (dt_Co < dt_Comin)   {dt_Comin = dt_Co;}
            }
        }
    }
    // adjust dt due to async coupling
    if ((*param)->sync_coupling == 0)  {
        for (ii = 0; ii < (*param)->n3ci; ii++) {
            if (gmap->istop[ii] == 1)   {
                jj = gmap->top2d[ii];
                if (data->dept[jj] > 0.0 & data->eta[jj] != data->etan[jj])    {
                    dt_eta = fabs(q_target * gmap->dz3d[ii] * (*param)->dt / data->Kz[gmap->icjckM[ii]] / (data->eta[jj]-data->etan[jj]));
                    if (dt_eta < dt_sync)  {dt_sync = dt_eta;}
                }
            }
        }
    }
    // adjust dt
    if ((*param)->iter_solve == 1)  {
        if ((*param)->n_iter < 4)   {
            (*param)->dtg = (*param)->dtg * 2.0;
        }
        else if ((*param)->n_iter > 8)  {
            (*param)->dtg = (*param)->dtg * 0.5;
        }
    }
    else {
        // (*param)->dtg = (*param)->dtg * r_inc;
        if (dq_max > 0.02)  {(*param)->dtg = (*param)->dtg * r_red;}
        else if (dq_max < 0.01) {(*param)->dtg = (*param)->dtg * r_inc;}
        if ((*param)->dtg > dt_Comin)   {(*param)->dtg = dt_Comin;}
    }
    // if ((*param)->sync_coupling == 0)  {
    //     if ((*param)->dtg > dt_sync)    {(*param)->dtg = dt_sync;}
    // }
    if ((*param)->dtg > (*param)->dt_max)  {(*param)->dtg = (*param)->dt_max;}
    if ((*param)->dtg < (*param)->dt_min)  {(*param)->dtg = (*param)->dt_min;}
    // unify dt for all processes
    if ((*param)->use_mpi == 1)
    {
        mpi_gather_double(&(*param)->dt_root[0], &(*param)->dtg, 1, root);
        (*param)->dtg = getMin((*param)->dt_root, (*param)->mpi_nx*(*param)->mpi_ny);
        mpi_bcast_double(&(*param)->dtg, 1, root);
    }
}
