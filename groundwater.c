
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>

// -----------------------------------------------------------------------------
// Groundwater module that solves saturated groundwater flow equation following
// MODFLOW.
// - ZhiLi 2019-01-09 -
// -----------------------------------------------------------------------------

#include "configuration.h"
#include "bathymetry.h"
#include "initialize.h"
#include "map.h"
#include "mpifunctions.h"
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

void groundwaterExchange(Data **data, Ground **ground, Maps *map, Gmaps *gmap, Config *setting, int irank, int nrank);
void computeConductance(Ground **ground, Data *data, Gmaps *gmap, Config *setting, int irank, int nrank);
void groundMatrixCoeff(Ground **ground, Data **data, Gmaps *gmap, Config *setting);
void setupGroundMatrix(Ground *ground, Gmaps *gmap, Config *setting, QMatrix A);
void solveGroundMatrix(Ground *ground, Gmaps *gmap, Config *setting, QMatrix A, Vector x, Vector z);
void getHead(Ground **ground, Gmaps *gmap, Config *setting, Vector x);
void enforceSidewallBC(Ground **ground, Gmaps *gmap, Config *setting);
void computeSeepage(Data **data, Ground **ground, Maps *map, Gmaps *gmap, Config *setting);
void computeFlowRate(Ground **ground, Data *data, Gmaps *gmap, Config *setting);
void updateWaterContent(Ground **ground, Gmaps *gmap, Config *setting);
void groundScalarAdvection(Ground **ground, Data *data, Gmaps *gmap, Config *setting);
void groundScalarDiffusion(Ground **ground, Data *data, Gmaps *gmap, Config *setting);
void groundScalarTransport(Ground **ground, Data **data, Gmaps *gmap, Maps *map, Config *setting, int irank, int nrank);
double updateSaturation(double h, Config *setting);
double updateDSDH(double h, Config *setting);
double hydraulicCond(double h, double Ks, Config *setting);
double faceCond(double K1, double K2);

// =========== Top-level groundwater exchange ==========
void groundwaterExchange(Data **data, Ground **ground, Maps *map, Gmaps *gmap, Config *setting, int irank, int nrank)
{

    QMatrix AA;
    Q_Constr(&AA, "A", setting->N3ci, False, Rowws, Normal, True);
    Vector zz;
    V_Constr(&zz, "z", setting->N3ci, Normal, True);
    Vector xx;
    V_Constr(&xx, "x", setting->N3ci, Normal, True);
    computeConductance(ground, *data, gmap, setting, irank, nrank);
    groundMatrixCoeff(ground, data, gmap, setting);
    setupGroundMatrix(*ground, gmap, setting, AA);
    solveGroundMatrix(*ground, gmap, setting, AA, xx, zz);
    getHead(ground, gmap, setting, xx);
    enforceSidewallBC(ground, gmap, setting);
    Q_Destr(&AA);
    V_Destr(&xx);
    V_Destr(&zz);

    if (setting->useMPI == 1)
    {mpiexchangeGround((*ground)->h, gmap, setting, irank, nrank);}
    computeSeepage(data, ground, map, gmap, setting);
    computeFlowRate(ground, *data, gmap, setting);
    updateWaterContent(ground, gmap, setting);
}

// ===== Function --- Calculate conductivity =====
double hydraulicCond(double h, double Ks, Config *setting)
{
  double n, m, ah, K, S, Kr;
  n = setting->a1;
  m = 1.0 - 1.0 / n;
  ah = fabs(setting->a2 * h);
  S = pow(1 + pow(ah,n), -m);
  Kr = pow(S,0.5) * pow(1-pow(1-pow(S,1.0/m),m), 2.0);
  if (Kr > 1.0)  {Kr = 1.0;}
  if (Kr < 0.01)   {Kr = 0.0;}
  K = Ks * Kr;
  if (Ks == 0.0)    {K = 0.0;}
  return K;
}

// ===== Average to get conductivity on face =====
double faceCond(double K1, double K2)
{return 0.5 * (K1 + K2);}

// FUNCTION --- Calculate saturation
double updateSaturation(double h, Config *setting)
{
    double n, m, ah, wc, wcr, wcs, S;
    n = setting->a1;
    m = 1.0 - 1.0 / n;
    ah = fabs(setting->a2 * h);
    wcr = setting->Sres * setting->porosity;
    wcs = setting->porosity;
    S = pow(1 + pow(ah,n), -m);
    wc = wcr + (wcs - wcr) * S;
    if (wc > wcs)   {wc = wcs;}
    if (wc < wcr)   {wc = wcr;}
    return wc;
}

// FUNCTION --- Calculate dS/dh
double updateDSDH(double h, Config *setting)
{
    double n, m, ah, dSdh, wcr, wcs;
    n = setting->a1;
    m = 1.0 - 1.0 / n;
    ah = fabs(setting->a2 * h);
    wcr = setting->Sres * setting->porosity;
    wcs = setting->porosity;
    dSdh = (setting->a2 * n * m * (wcs-wcr) * pow(ah,n-1)) / pow((1+pow(ah,n)),m+1);
    return dSdh;
}

// =============== Compute conductance ===============
void computeConductance(Ground **ground, Data *data, Gmaps *gmap, Config *setting, int irank, int nrank)
{
    int ii;
    double coef, dzavg, ah, n, Kc, Kp, Km;
    coef = setting->dtg * setting->porosity;
    // store layer information
    for (ii = 0; ii < setting->N2ci; ii++)
    {(*ground)->nlay[ii] = gmap->nlay[ii];}
    // specific storage
    (*ground)->SS = 1000.0 * setting->g * (setting->compS + setting->porosity * setting->compW);
    // loop over all cells to compute conductance
    for (ii = 0; ii < setting->N3ci; ii++)
    {
        (*ground)->hOld[ii] = (*ground)->h[ii];
        // ===== Cx =====
        Kc = hydraulicCond((*ground)->h[ii], setting->Kxx, setting);
        if (gmap->iPjckc[ii] <= setting->N3ci)
        {Kp = hydraulicCond((*ground)->h[gmap->iPjckc[ii]], setting->Kxx, setting);}
        else
        {Kp = Kc;}
        (*ground)->Cx[ii] = setting->dtg * faceCond(Kc,Kp) / (setting->dx * setting->dx);
        if (gmap->iMjckc[ii] >= setting->N3ci)
        {(*ground)->Cx[gmap->iMjckc[ii]] = setting->dtg * Kc / (setting->dx * setting->dx);}
        // ===== Cy =====
        Kc = hydraulicCond((*ground)->h[ii], setting->Kyy, setting);
        if (gmap->icjPkc[ii] <= setting->N3ci)
        {Kp = hydraulicCond((*ground)->h[gmap->icjPkc[ii]], setting->Kyy, setting);}
        else
        {Kp = Kc;}
        (*ground)->Cy[ii] = setting->dtg * faceCond(Kc,Kp) / (setting->dy * setting->dy);
        if (gmap->icjMkc[ii] >= setting->N3ci)
        {(*ground)->Cy[gmap->icjMkc[ii]] = setting->dtg * Kc / (setting->dy * setting->dy);}
        // ===== Cz ======
        // NOTE: Unlike Cx and Cy, Cz[ii] is its kM face (upward face)
        if (gmap->actv[ii] == 1)
        {
            // if (gmap->ii[ii] == 1 & gmap->jj[ii] == 0)
            // {
            //     printf("ii,h = %d,%f\n",ii,(*ground)->h[ii]);
            // }
            Kc = hydraulicCond((*ground)->h[ii], setting->Kzz, setting);
            if (gmap->istop[ii] == 1)
            {
                if (data->depth[ii] > 0.0)
                {(*ground)->Cz[ii] = setting->dtg * Kc / (0.5 * gmap->dz3d[ii] * gmap->dz3d[ii]);}
                else
                {(*ground)->Cz[ii] = 0.0;}
            }
            else
            {
                Km = hydraulicCond((*ground)->h[gmap->icjckM[ii]], setting->Kzz, setting);
                (*ground)->Cz[ii] = setting->dtg * faceCond(Kc,Km) / (gmap->dz3d[ii] * 0.5*(gmap->dz3d[ii]+gmap->dz3d[gmap->icjckM[ii]]));
            }
        }
        else
        {(*ground)->Cz[ii] = 0.0;}
        // ===== cell volume =====
        (*ground)->V[ii] = setting->dx * setting->dy * (*ground)->wc[ii] * gmap->dz3d[ii];

    }
    // zero conductance for non-active faces
    for (ii = 0; ii < setting->N3ci; ii++)
    {
        if (gmap->actv[ii] == 0)
        {
            (*ground)->Cx[ii] = 0.0;
            (*ground)->Cx[gmap->iMjckc[ii]] = 0.0;
            (*ground)->Cy[ii] = 0.0;
            (*ground)->Cy[gmap->icjMkc[ii]] = 0.0;
            (*ground)->Cz[ii] = 0.0;
            (*ground)->Cz[gmap->icjckP[ii]] = 0.0;
        }
        else
        {
            if (gmap->iPjckc[ii] >= setting->N3ci)
            {(*ground)->Cx[ii] = 0.0;}
            if (gmap->iMjckc[ii] >= setting->N3ci)
            {(*ground)->Cx[gmap->iMjckc[ii]] = 0.0;}
            if (gmap->icjPkc[ii] >= setting->N3ci)
            {(*ground)->Cy[ii] = 0.0;}
            if (gmap->icjMkc[ii] >= setting->N3ci)
            {(*ground)->Cy[gmap->icjMkc[ii]] = 0.0;}
            if (gmap->icjckP[ii] >= setting->N3ci)
            {(*ground)->Cz[gmap->icjckP[ii]] = 0.0;}
        }
    }

}

// =============== Compute the coefficients of the matrix A ===============
void groundMatrixCoeff(Ground **ground, Data **data, Gmaps *gmap, Config *setting)
{
    int ii,jj,kk;
    double depth, dwdh, Ve, wc, Kp, Km, Kc, Kpf, Kmf, allV, Vvoid;
    for (ii = 0; ii < setting->N3ci; ii++)
    {
        if (gmap->actv[ii] == 1)
        {
            (*ground)->GnCt[ii] = (*ground)->SS;
            if (setting->useUnSat == 1)
            {
                wc = updateSaturation((*ground)->h[ii], setting);
                dwdh = updateDSDH((*ground)->h[ii], setting);
                (*ground)->GnCt[ii] = (*ground)->GnCt[ii] * wc / setting->porosity + dwdh;
            }
            // XP
            if ((*ground)->Cx[ii] > 0.0)
            {
                (*ground)->GnXP[ii] = (*ground)->Cx[ii];
                (*ground)->GnCt[ii] += (*ground)->GnXP[ii];
            }
            else
            {(*ground)->GnXP[ii] = 0.0;}
            // XM
            if ((*ground)->Cx[gmap->iMjckc[ii]] > 0.0)
            {
                (*ground)->GnXM[ii] = (*ground)->Cx[gmap->iMjckc[ii]];
                (*ground)->GnCt[ii] += (*ground)->GnXM[ii];
            }
            else
            {(*ground)->GnXM[ii] = 0.0;}
            // YP
            if ((*ground)->Cy[ii] > 0.0)
            {
                (*ground)->GnYP[ii] = (*ground)->Cy[ii];
                (*ground)->GnCt[ii] += (*ground)->GnYP[ii];
            }
            else
            {(*ground)->GnYP[ii] = 0.0;}
            // YM
            if ((*ground)->Cy[gmap->icjMkc[ii]] > 0.0)
            {
                (*ground)->GnYM[ii] = (*ground)->Cy[gmap->icjMkc[ii]];
                (*ground)->GnCt[ii] += (*ground)->GnYM[ii];
            }
            else
            {(*ground)->GnYM[ii] = 0.0;}
            // ZP
            if (gmap->icjckP[ii] >= 0 & gmap->icjckP[ii] < setting->N3ci)
            {
                (*ground)->GnZP[ii] = (*ground)->Cz[gmap->icjckP[ii]];
                (*ground)->GnCt[ii] += (*ground)->GnZP[ii];
            }
            else
            {(*ground)->GnZP[ii] = 0.0;}
            // ZM
            if (gmap->actv[gmap->iMjckc[ii]] == 1)
            {
                (*ground)->GnZM[ii] = (*ground)->Cz[ii];
                (*ground)->GnCt[ii] += (*ground)->GnZM[ii];
            }
            else
            {(*ground)->GnZM[ii] = 0.0;}
        }
        else
        {
            (*ground)->GnCt[ii] = (*ground)->SS;
            (*ground)->GnXP[ii] = 0.0;
            (*ground)->GnXM[ii] = 0.0;
            (*ground)->GnYP[ii] = 0.0;
            (*ground)->GnYM[ii] = 0.0;
            (*ground)->GnZP[ii] = 0.0;
            (*ground)->GnZM[ii] = 0.0;
        }
        // Source term B
        (*ground)->B[ii] = (*ground)->h[ii] * (*ground)->SS;
        // Unsaturated zone
        if (setting->useUnSat == 1)
        {
            wc = updateSaturation((*ground)->h[ii], setting);
            dwdh = updateDSDH((*ground)->h[ii], setting);
            (*ground)->B[ii] = (*ground)->B[ii] * wc / setting->porosity + dwdh * (*ground)->h[ii];
            // evaporation
            if (setting->useEvap == 1 & gmap->istop[ii] == 1 & gmap->actv[ii] == 1)
            {
                Ve = setting->qe * setting->dtg * setting->dx * setting->dy * setting->porosity;
                kk = ii;
                while (Ve > 0.0)
                {
                    Vvoid = gmap->dz3d[kk] * setting->dx * setting->dy * setting->porosity;
                    allV = (*ground)->V[kk] - Vvoid * setting->Sres;
                    if (Ve <= allV)
                    {
                        (*ground)->B[kk] -= setting->qe * setting->dtg / gmap->dz3d[kk];
                        Ve = 0.0;
                    }
                    else
                    {
                        (*ground)->B[kk] -= allV / (setting->dx * setting->dy * setting->porosity);
                        Ve -= allV;
                        kk = gmap->icjckP[kk];
                    }
                }
            }
        }

        // Gravity terms and surface/subsurface interface BC
        Kc = hydraulicCond((*ground)->h[ii], setting->Kzz, setting);
        if (gmap->icjckP[ii] < setting->N3ci & gmap->icjckP[ii] >= 0)
        {
            Kp = hydraulicCond((*ground)->h[gmap->icjckP[ii]], setting->Kzz, setting);
            Kpf = faceCond(Kc,Kp);
        }
        else
        {Kpf = 0.0;}
        if (gmap->istop[ii] != 1)
        {
            Km = hydraulicCond((*ground)->h[gmap->icjckM[ii]], setting->Kzz, setting);
            Kmf = faceCond(Kc,Km);
        }
        else
        {Kmf = Kc;}
        (*ground)->B[ii] += setting->dtg * (Kmf - Kpf) / gmap->dz3d[ii];


        // Add surface flow BC to B
        if (gmap->actv[ii] == 1 & gmap->istop[ii] == 1)
        {
            if ((*data)->depth[gmap->top2D[ii]] > 0)
            {(*ground)->B[ii] += (*ground)->Cz[ii] * (*data)->depth[gmap->top2D[ii]];}
            else
            {(*ground)->B[ii] -= setting->dtg * Kmf / gmap->dz3d[ii];}
        }

        if (gmap->actv[ii] == 0)
        {(*ground)->B[ii] = (*ground)->SS;}

    }
}

// =============== Generate Matrix A and z ===============
void setupGroundMatrix(Ground *ground, Gmaps *gmap, Config *setting, QMatrix A)
{
    // Assume no flow BC for XP, XM, YP, YM and ZP boundary
    // Assume Direchlet BC for ZM boundary
    // Q_SetLen(matrix, row, num_elem)
    // Q_SetEntry(matrix, row, index, col, value)
    size_t ii, kk, N, dist, row;
    int jj;
	row = 0;
    for (ii = 1; ii <= setting->N3ci; ii++)
    {
		//if (gmap->actv[ii] == 1)
		//{
		    jj = ii - 1;
		    N = 1;
		    kk = 0;
		    // compute number of elements in this row
		    if (ground->GnZP[jj] != 0) {N += 1;}
		    if (ground->GnZM[jj] != 0 & gmap->istop[jj] != 1) {N += 1;}
		    if (ground->GnXP[jj] != 0) {N += 1;}
		    if (ground->GnXM[jj] != 0) {N += 1;}
		    if (ground->GnYP[jj] != 0) {N += 1;}
		    if (ground->GnYM[jj] != 0) {N += 1;}
		    Q_SetLen(&A, ii, N);
		    // set YM entry
		    if (ground->GnYM[jj] != 0)
		    {
		        dist = jj - gmap->icjMkc[jj];
		        Q_SetEntry(&A, ii, kk, ii-dist, -ground->GnYM[jj]);
		        kk++;
		    }
		    // set XM entry
		    if (ground->GnXM[jj] != 0)
		    {
		        dist = jj - gmap->iMjckc[jj];
		        Q_SetEntry(&A, ii, kk, ii-dist, -ground->GnXM[jj]);
		        kk++;
		    }
		    // set ZM entry
		    if (ground->GnZM[jj] != 0 & gmap->istop[jj] != 1)
		    {
		        Q_SetEntry(&A, ii, kk, ii-1, -ground->GnZM[jj]);
		        kk++;
		    }
		    // set Ct entry
		    Q_SetEntry(&A, ii, kk, ii, ground->GnCt[jj]);
		    kk++;
		    // set ZP entry
		    if (ground->GnZP[jj] != 0)
		    {
		        Q_SetEntry(&A, ii, kk, ii+1, -ground->GnZP[jj]);
		        kk++;
		    }
		    // set XP entry
		    if (ground->GnXP[jj] != 0)
		    {
		        dist = gmap->iPjckc[jj] - jj;
		        if (ii+dist < setting->N3ci)
		        {Q_SetEntry(&A, ii, kk, ii+dist, -ground->GnXP[jj]);}
		        kk++;
		    }
		    // set YP entry
		    if (ground->GnYP[jj] != 0)
		    {
		        dist = gmap->icjPkc[jj] - jj;
		        if (ii+dist < setting->N3ci)
		        {Q_SetEntry(&A, ii, kk, ii+dist, -ground->GnYP[jj]);}
		        kk++;
		    }
		//}
    }
}

// =============== Solve the linear system Ax = z ===============
void solveGroundMatrix(Ground *ground, Gmaps *gmap, Config *setting, QMatrix A, Vector x, Vector z)
{

    size_t ii;
    // create the z Vector
    for (ii = 1; ii <= setting->N3ci; ii++)
    {
		//if (gmap->actv[ii] == 1)
		  V_SetCmp(&z, ii, ground->B[ii-1]);
	   }

     // for (ii = 0; ii < setting->N3ci; ii++)
     // {
     //     if (gmap->ii[ii] == 1){
     //     printf("jj,kk,actv,Gc,xp,xm,yp,ym,zp,zm,B = %d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f\n",gmap->jj[ii], \
     //            gmap->kk[ii],gmap->actv[ii],ground->GnCt[ii], \
     //            ground->GnXP[ii],ground->GnXM[ii],ground->GnYP[ii],ground->GnYM[ii], \
     //            ground->GnZP[ii],ground->GnZM[ii],ground->B[ii]);}
     //     if (gmap->kk[ii] == 5) {printf("-------------------- \n");}
     // }
    // initialize x
    V_SetAllCmp(&x, 0.0);
    // set stopping criteria
    SetRTCAccuracy(setting->eps);
    // solve the linear system
    //JacobiIter(&A, &x, &z, setting->maxIter, NULL, 1.1);
    CGIter(&A, &x, &z, setting->maxIter, SSORPrecond, 1);
}

// =============== Assign x back to head ===============
void getHead(Ground **ground, Gmaps *gmap, Config *setting, Vector x)
{
    size_t ii;
    for (ii = 0; ii < setting->N3ci; ii++)
    {(*ground)->h[ii] = V_GetCmp(&x, ii+1);}
    // for (ii = 0; ii < setting->N3ci; ii++)
    // {
    //     if (gmap->ii[ii] == 1){
    //     printf("ii,jj,kk,actv,h = %zu,%d,%d,%d,%f\n",ii,gmap->jj[ii], \
    //            gmap->kk[ii],gmap->actv[ii],(*ground)->h[ii]);}
    //     if (gmap->kk[ii] == 5) {printf("-------------------- \n");}
    // }

}

// =============== Neumann BC for h along sidewall ===============
void enforceSidewallBC(Ground **ground, Gmaps *gmap, Config *setting)
{
	int ii;
	for (ii = 0; ii < setting->N3ci; ii++)
	{
		if (gmap->ii[ii] == 0)	{(*ground)->h[gmap->iMjckc[ii]] = (*ground)->h[ii];}
		if (gmap->ii[ii] == setting->nx-1)	{(*ground)->h[gmap->iPjckc[ii]] = (*ground)->h[ii];}
		if (gmap->jj[ii] == 0)	{(*ground)->h[gmap->icjMkc[ii]] = (*ground)->h[ii];}
		if (gmap->jj[ii] == setting->ny-1)	{(*ground)->h[gmap->icjPkc[ii]] = (*ground)->h[ii];}
	}
}

// ========== Compute flow rate from subsurface to surface ==========
void computeSeepage(Data **data, Ground **ground, Maps *map, Gmaps *gmap, Config *setting)
{
    int ii, kk;
    float d;
    for (ii = 0; ii < setting->N3ci; ii++)
    {
        // find 2D index of top layer 3D cell
        if (gmap->istop[ii] == 1)
        {
            kk = gmap->top2D[ii];
            // update head boundary condition
            (*ground)->h[gmap->icjckM[ii]] = (*data)->depth[kk];
            // compute flux (positive = upward)
            (*data)->Qseep[kk] = (*ground)->Cz[ii] * (*ground)->V[ii] * ((*ground)->h[ii] - (*data)->depth[kk]) / setting->dtg;
        }
    }
}

// ========== Compute flow rate for subsurface faces ==========
void computeFlowRate(Ground **ground, Data *data, Gmaps *gmap, Config *setting)
{
    int ii;
    for (ii = 0; ii < setting->N3ci; ii++)
    {
        // assume no lateral exchange
        // Quu
        if (gmap->actv[ii] == 1 & gmap->actv[gmap->iPjckc[ii]] == 1 & \
                gmap->iPjckc[ii] < setting->N3ci)
        {(*ground)->Quu[ii] = (*ground)->Cx[ii] * ((*ground)->h[ii] - (*ground)->h[gmap->iPjckc[ii]]);}
        else
        {(*ground)->Quu[ii] = 0.0;}
        // Qvv
        if (gmap->actv[ii] == 1 & gmap->actv[gmap->icjPkc[ii]] == 1 & \
                gmap->icjPkc[ii] < setting->N3ci)
        {(*ground)->Qvv[ii] = (*ground)->Cy[ii] * ((*ground)->h[ii] - (*ground)->h[gmap->icjPkc[ii]]);}
        else
        {(*ground)->Qvv[ii] = 0.0;}
        // Qww
        if (gmap->actv[ii] == 1)
        {
            if (gmap->istop[ii] == 1)
            {
                if (data->depth[gmap->top2D[ii]] == 0)
                {(*ground)->Qww[ii] = 0.0;}
                else
                {(*ground)->Qww[ii] = setting->dtg * data->Qseep[gmap->top2D[ii]];}
            }
            else
            {(*ground)->Qww[ii] = (*ground)->Cz[ii] * ((*ground)->h[ii] - (*ground)->h[gmap->icjckM[ii]]);}
            if (gmap->icjckP[ii] >= setting->N3ci)
            {(*ground)->Qww[gmap->icjckP[ii]] = 0.0;}
        }
        else
        {(*ground)->Qww[ii] = 0.0;}
    }
}

// ========== Update water content ==========
void updateWaterContent(Ground **ground, Gmaps *gmap, Config *setting)
{
    int ii;
    double a, m, n, wcs, wcr;
    n = setting->a1;
    a = setting->a2;
    m = 1.0 - 1.0 / n;
    wcr = setting->Sres * setting->porosity;
    wcs = setting->porosity;
    for (ii = 0; ii < setting->N3ci; ii++)
    {
        if (gmap->actv[ii] == 1)
        {
            (*ground)->wc[ii] -= setting->dtg * ((*ground)->Quu[gmap->iMjckc[ii]] - (*ground)->Quu[ii]);
            (*ground)->wc[ii] -= setting->dtg * ((*ground)->Qvv[gmap->icjMkc[ii]] - (*ground)->Qvv[ii]);
            (*ground)->wc[ii] -= setting->dtg * ((*ground)->Qww[ii] - (*ground)->Qww[gmap->icjckP[ii]]);
            (*ground)->wc[ii] = (*ground)->wc[ii] / (1.0 + (*ground)->SS * ((*ground)->h[ii]-(*ground)->hOld[ii]) / setting->porosity);
            if ((*ground)->wc[ii] >= wcs)
            {(*ground)->wc[ii] = wcs;}
            else
            {
                if ((*ground)->wc[ii] < wcr+0.01)
                {(*ground)->wc[ii] = wcr+0.01;}
                (*ground)->h[ii] = -(1.0/a) * pow((pow((wcs-wcr)/((*ground)->wc[ii]-wcr),1.0/m) - 1.0),1.0/n);
            }
        }
    }
}

// ========== Subsurface scalar advection ==========
void groundScalarAdvection(Ground **ground, Data *data, Gmaps *gmap, Config *setting)
{
	int ii;
	for (ii = 0; ii < setting->N3ci; ii++)
	{
		double SiP, SiM, SjP, SjM, SkP, SkM;
		// scalar on iP face
		if ((*ground)->Quu[ii] != 0)
		{
			if ((*ground)->Quu[ii] > 0) {SiP = (*ground)->S[ii];}
			else {SiP = (*ground)->S[gmap->iPjckc[ii]];}
		}
		else
		{SiP = 0;}
		// scalar on iM face
		if ((*ground)->Quu[gmap->iMjckc[ii]] != 0)
		{
			if ((*ground)->Quu[gmap->iMjckc[ii]] > 0) {SiM = (*ground)->S[gmap->iMjckc[ii]];}
			else {SiM = (*ground)->S[ii];}
		}
		else
		{SiM = 0;}
        // scalar on jP face
		if ((*ground)->Qvv[ii] != 0)
		{
			if ((*ground)->Qvv[ii] > 0) {SjP = (*ground)->S[ii];}
			else {SjP = (*ground)->S[gmap->icjPkc[ii]];}
		}
		else
		{SiM = 0;}
		// scalar on jM face
		if ((*ground)->Qvv[gmap->icjMkc[ii]] != 0)
		{
			if ((*ground)->Qvv[gmap->icjMkc[ii]] > 0) {SjM = (*ground)->S[gmap->icjMkc[ii]];}
			else {SjM = (*ground)->S[ii];}
		}
		else
		{SjM = 0;}
		// scalar on kM face
		if ((*ground)->Qww[ii] != 0)
		{
			if ((*ground)->Qww[ii] > 0) {SkM = (*ground)->S[ii];}
			else
      {
        if (gmap->istop[ii] == 1)
        {SkM = data->S[gmap->top2D[ii]];}
        else
        {SkM = (*ground)->S[gmap->icjckM[ii]];}
      }
		}
		else
		{SkM = 0;}
		// scalar on kP face
		if ((*ground)->Qww[gmap->icjckP[ii]] != 0)
		{
			if ((*ground)->Qww[gmap->icjckP[ii]] > 0) {SkP = (*ground)->S[gmap->icjckP[ii]];}
			else {SkP = (*ground)->S[ii];}
		}
		else
		{SkP = 0;}
		// Advective transport
		(*ground)->Sm[ii] = (*ground)->Sm[ii] + setting->dtg * \
			(-(*ground)->Quu[ii] * SiP + (*ground)->Quu[gmap->iMjckc[ii]] * SiM - \
			(*ground)->Qvv[ii] * SjP + (*ground)->Qvv[gmap->icjMkc[ii]] * SjM - \
			(*ground)->Qww[ii] * SkM + (*ground)->Qww[gmap->icjckP[ii]] * SkP);
	}

}

// ========== Subsurface scalar diffusion ==========
void groundScalarDiffusion(Ground **ground, Data *data, Gmaps *gmap, Config *setting)
{
	int ii;
	double AiP, AiM, AjP, AjM, AkP, AkM;
	double JiP, JiM, JjP, JjM, JkP, JkM;
	double dz;
	for (ii = 0; ii < setting->N3ci; ii++)
	{
		// compute face area
		AiP = setting->dy * dz * setting->porosity;
		AiM = setting->dy * dz * setting->porosity;
		AjP = setting->dx * dz * setting->porosity;
		AjM = setting->dx * dz * setting->porosity;
		AkP = setting->dy * setting->dx * setting->porosity;
		AkM = setting->dy * setting->dx * setting->porosity;
		// iP face flux
		if ((*ground)->Cx[ii] > 0 & gmap->iPjckc[ii] < setting->N3ci)
		{JiP = AiP * ((*ground)->S[gmap->iPjckc[ii]] - (*ground)->S[ii]) / setting->dx;}
		else
		{JiP = 0;}
		// iM face flux
		if ((*ground)->Cx[gmap->iMjckc[ii]] > 0 & gmap->iMjckc[ii] < setting->N3ci)
		{JiM = AiM * ((*ground)->S[ii] - (*ground)->S[gmap->iMjckc[ii]]) / setting->dx;}
		else
		{JiM = 0;}
		// jP face flux
		if ((*ground)->Cy[ii] > 0 & gmap->icjPkc[ii] < setting->N3ci)
		{JjP = AjP * ((*ground)->S[gmap->icjPkc[ii]] - (*ground)->S[ii]) / setting->dy;}
		else
		{JjP = 0;}
		// jM face flux
		if ((*ground)->Cy[gmap->icjMkc[ii]] > 0 & gmap->icjMkc[ii] < setting->N3ci)
		{JjM = AjM * ((*ground)->S[ii] - (*ground)->S[gmap->icjMkc[ii]]) / setting->dy;}
		else
		{JjM = 0;}
		// kP face flux
		if ((*ground)->Cz[gmap->icjckP[ii]] > 0 & gmap->icjckP[ii] >= 0)
		{
      dz = 0.5 * (gmap->dz3d[ii] + gmap->dz3d[gmap->icjckP[ii]]);
      JkP = AkP * ((*ground)->S[ii] - (*ground)->S[gmap->icjckP[ii]]) / dz;
    }
		else
		{JkP = 0;}
		// kM face flux
		if ((*ground)->Cz[ii] > 0)
		{
      dz = 0.5 * (gmap->dz3d[ii] + gmap->dz3d[gmap->icjckM[ii]]);
      if (gmap->istop[ii] == 1)
      {JkM = 2.0 * AkM * (data->S[gmap->top2D[ii]] - (*ground)->S[ii]) / gmap->dz3d[ii];}
      else
      {JkM = AkM * ((*ground)->S[gmap->icjckM[ii]] - (*ground)->S[ii]) / dz;}
    }
		else
		{JkM = 0;}
		// diffusive transport
		(*ground)->Sm[ii] = (*ground)->Sm[ii] + setting->dtg * \
			(setting->Kmx * (JiP - JiM) + setting->Kmy * (JjP - JjM) + setting->Kmx * (JkM - JkP));
	}
}

// ========== Subsurface scalar transport ==========
void groundScalarTransport(Ground **ground, Data **data, Gmaps *gmap, Maps *map, Config *setting, int irank, int nrank)
{
	int ii, jj, *N, flag;
	double *Smax, *Smin, *Sarr;
	N = malloc(setting->N3ci * sizeof(int));
	Smax = malloc(setting->N3ci * sizeof(double));
	Smin = malloc(setting->N3ci * sizeof(double));
  // estimate scalar mass
  for (ii = 0; ii < setting->N3ci; ii++)
  {(*ground)->Sm[ii] = (*ground)->S[ii] * (*ground)->V[ii];}
	// scalar limiter
	for (ii = 0; ii < setting->N3ci; ii++)
	{
		int iP=0, iM=0, jP=0, jM=0, kP=0, kM=0;
		N[ii] = 0;
		if ((*ground)->Cx[ii] > 0 & gmap->iPjckc[ii] < setting->N3ci) {N[ii] += 1; iP = 1;}
		if ((*ground)->Cx[gmap->iMjckc[ii]] > 0 & gmap->iMjckc[ii] < setting->N3ci) {N[ii] += 1; iM = 1;}
		if ((*ground)->Cy[ii] > 0 & gmap->icjPkc[ii] < setting->N3ci) {N[ii] += 1; jP = 1;}
		if ((*ground)->Cy[gmap->icjMkc[ii]] > 0 & gmap->icjMkc[ii] < setting->N3ci) {N[ii] += 1; jM = 1;}
		if ((*ground)->Cz[ii] > 0) {N[ii] += 1; kM = 1;}
		if (gmap->icjckP[ii] >= 0 & gmap->icjckP[ii] < setting->N3ci) {N[ii] += 1; kP = 1;}
		if (N[ii] > 0)
		{
			int jj = 0;
			Sarr = malloc(N[ii] * sizeof(double));
			if (iP == 1) {Sarr[jj] = (*ground)->S[gmap->iPjckc[ii]];	jj++;}
			if (iM == 1) {Sarr[jj] = (*ground)->S[gmap->iMjckc[ii]];	jj++;}
			if (jP == 1) {Sarr[jj] = (*ground)->S[gmap->icjPkc[ii]];	jj++;}
			if (jM == 1) {Sarr[jj] = (*ground)->S[gmap->icjMkc[ii]];	jj++;}
			if (kP == 1) {Sarr[jj] = (*ground)->S[gmap->icjckP[ii]];	jj++;}
			if (kM == 1) {Sarr[jj] = (*ground)->S[gmap->icjckM[ii]];	jj++;}
			Smax[ii] = getMax(Sarr, N[ii]);
			Smin[ii] = getMin(Sarr, N[ii]);
			iP = 0;
			iM = 0;
			jP = 0;
			jM = 0;
			kP = 0;
			kM = 0;
			free(Sarr);
		}
	}
	// advection
	groundScalarAdvection(ground, *data, gmap, setting);
	// diffusion
	groundScalarDiffusion(ground, *data, gmap, setting);
	// update scalar
    flag = 0;
	for (ii = 0; ii < setting->N3ci; ii++)
	{
		if ((*ground)->V[ii] > 0)
		{(*ground)->S[ii] = (*ground)->Sm[ii] / (*ground)->V[ii];}
		else
		{(*ground)->S[ii] = 0;}

		if (N[ii] > 0)
		{
      // CFL checking
      if (gmap->istop[ii] == 1 & (*ground)->S[ii] > Smax[ii])
      {flag = 1;}
			if ((*ground)->S[ii] > Smax[ii]) {(*ground)->S[ii] = Smax[ii];}
			if ((*ground)->S[ii] < Smin[ii]) {(*ground)->S[ii] = Smin[ii];}
		}

	}
    if (flag == 1)
    {
        printf("WARNING: Groundwater time step might be too big for scalar transport!\n");
    }
	// set top boundary condition
	for (ii = 0; ii < setting->N3ci; ii++)
	{
		if (gmap->istop[ii] == 1)
		{
			jj = gmap->top2D[ii];
			(*ground)->S[gmap->icjckM[ii]] = (*data)->S[jj];
			(*data)->Ssub[jj] = (*ground)->S[ii];
			if ((*ground)->Cz[ii] > 0)
			{(*data)->dz[jj] = gmap->htop[jj];}
			else
			{(*data)->dz[jj] = 0;}
		}
	}
	// set side boundary condition
	for (ii = 0; ii < setting->N3ci; ii++)
	{
		if (gmap->ii[ii] == 0)	{(*ground)->S[gmap->iMjckc[ii]] = (*ground)->S[ii];}
		if (gmap->ii[ii] == setting->nx-1)	{(*ground)->S[gmap->iPjckc[ii]] = (*ground)->S[ii];}
		if (gmap->jj[ii] == 0)	{(*ground)->S[gmap->icjMkc[ii]] = (*ground)->S[ii];}
		if (gmap->jj[ii] == setting->ny-1)	{(*ground)->S[gmap->icjPkc[ii]] = (*ground)->S[ii];}
        if (gmap->ii[ii] == 0)	{(*ground)->Sm[gmap->iMjckc[ii]] = (*ground)->Sm[ii];}
		if (gmap->ii[ii] == setting->nx-1)	{(*ground)->Sm[gmap->iPjckc[ii]] = (*ground)->Sm[ii];}
		if (gmap->jj[ii] == 0)	{(*ground)->Sm[gmap->icjMkc[ii]] = (*ground)->Sm[ii];}
		if (gmap->jj[ii] == setting->ny-1)	{(*ground)->Sm[gmap->icjPkc[ii]] = (*ground)->Sm[ii];}
	}

	free(N);
	free(Smax);
	free(Smin);

	// exchange between threads
	if (setting->useMPI == 1)
	{
		mpiexchangeGround((*ground)->S, gmap, setting, irank, nrank);
		mpiexchangeGround((*ground)->Sm, gmap, setting, irank, nrank);
	}
}
