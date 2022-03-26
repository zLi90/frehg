// Sub-functions called in the solver
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<mpi.h>

#include"configuration.h"
#include"initialize.h"
#include"map.h"
#include"mpifunctions.h"

double darcy_flux(Data *data, Map *gmap, int ic, double dhc, double dhs, char *axis, char *dir, Config *param, int irank);
double compute_residual(Data **data, Map *gmap, int ii, char *jaco, Config *param, int irank);
double compute_jacobian_fd(Data **data, Map *gmap, int ii, char* axis, Config *param, int irank);
double compute_wch(Data *data, double h, int ii, Config *param);
double compute_hwc(Data *data, int ii, Config *param);
double compute_ch(Data *data, int ii, Config *param);
double compute_K(Data *data, double h, double Ks, int ii, Config *param);
// double compute_K(Data *data, double *Ksat, int ii, Config *param);
double compute_dKdwc(Data *data, double *Ksat, int ii, Config *param);
double compute_dKdh(Data *data, double *Ksat, int ii, Config *param);
double compute_dwcdh(Data *data, int ii, Config *param);
double tvd_superbee(double sp, double sc, double sm, double u, double delta, double dt, Config *param);

// >>>>> Compute Darcy flux between faces <<<<<
double darcy_flux(Data *data, Map *gmap, int ic, double dhc, double dhs, char *axis, char *dir, Config *param, int irank)
{
    int is, sign;
    double q, hc, hs, k, kc, ks, kface, visc, rho, delta, cos, sin, vseep;
    double Ks;
    if (strcmp(axis, "x") == 0)    {
        if (strcmp(dir, "back") == 0 & gmap->ii[ic] == 0)  {
            is = ic;    ic = gmap->iMjckc[is];  delta = 0.5*param->dx;
            sign = 0;
        }
        else {
            is = gmap->iPjckc[ic];  delta = param->dx;
            if (gmap->bot3d[is] > gmap->bot3d[ic])    {sign = 1;}   else    {sign = -1;}
            if (gmap->ii[ic] == param->nx-1)    {delta = 0.5*param->dx;}
        }
        visc = data->r_viscxp[ic];   rho = 0.0;
        sin = gmap->sinx[ic];   cos = gmap->cosx[ic];
        Ks = data->Ksx[ic];
    }
    else if (strcmp(axis, "y") == 0)    {
        if (strcmp(dir, "back") == 0 & gmap->jj[ic] == 0)  {
            is = ic;    ic = gmap->icjMkc[is];
            if (irank < param->mpi_nx)    {delta = 0.5*param->dy; sign = 0;}
            else {
                delta = param->dy;
                if (gmap->bot3d[is] > gmap->bot3d[ic])    {sign = 1;}
                else    {sign = -1;}
            }
        }
        else {
            is = gmap->icjPkc[ic];  delta = param->dy;
            if (gmap->bot3d[is] > gmap->bot3d[ic])    {sign = 1;}   else    {sign = -1;}
            if (gmap->jj[ic] == param->ny-1)    {
                if (irank >= param->mpi_nx*(param->mpi_ny-1))   {delta = 0.5*param->dy;}
                else {delta = param->dy;}
            }
        }
        visc = data->r_viscyp[ic];   rho = 0.0;
        sin = gmap->siny[ic];   cos = gmap->cosy[ic];
        Ks = data->Ksy[ic];
    }
    else if (strcmp(axis, "z") == 0)    {
        if (strcmp(dir, "back") == 0 & gmap->kk[ic] == 0)  {
            is = ic;    ic = gmap->icjckM[is];  delta = 0.5*gmap->dz3d[is];
            Ks = data->Ksz[is];
            if (param->bctype_GW[5] == 0)   {Ks = 0.0;}
        }
        else {
            is = gmap->icjckP[ic];  delta = 0.5*(gmap->dz3d[ic]+gmap->dz3d[is]);
            Ks = data->Ksz[is];
            if (gmap->kk[ic] == param->nz-1)    {
                delta = 0.5*gmap->dz3d[ic]; if (param->bctype_GW[4] == 0)   {Ks = 0.0;}
            }
        }
        visc = data->r_visczp[ic];   rho = data->r_rhozp[ic];
        sin = 0.0;   cos = 1.0; sign = 0;
    }
    hc = data->h[ic]+dhc;   hs = data->h[is]+dhs;

    if (is < param->n3ci)  {
        ks = compute_K(data, hs, Ks, is, param);
        if (ic < param->n3ci)   {kc = compute_K(data, hc, Ks, ic, param);}
        else {kc = ks;}
    }
    else {
        kc = compute_K(data, hc, Ks, ic, param);
        ks = kc;
    }
    kface = 0.5 * (kc + ks);
    // saturated top
    if (strcmp(axis, "z") == 0 & strcmp(dir, "back") == 0)  {kface = data->Ksz[is];}

    q = kface * visc * cos * ((hs - hc)/delta - rho) + sign * kface * sin;
    // enforce boundary condition
    if (strcmp(axis, "y") == 0) {
        if (strcmp(dir, "back") == 0)   {
            if (gmap->jj[ic] == 0 & irank < param->mpi_nx) {
                if (param->bctype_GW[3] == 2)   {q = param->qym;}
                else if (param->bctype_GW[3] == 0)  {q = 0.0;}
            }
        }
        else  {
            if (gmap->jj[ic] == param->ny-1 & irank >= param->mpi_nx*(param->mpi_ny-1)) {
                if (param->bctype_GW[2] == 2)   {q = param->qyp;}
                else if (param->bctype_GW[2] == 0)  {q = 0.0;}
            }
        }
    }
    if (strcmp(axis, "z") == 0) {
        if (strcmp(dir, "back") == 0 & gmap->istop[is] == 1)   {
            // surface-subsurface exchange
            if (param->sim_shallowwater == 1)   {
                if (data->dept[gmap->top2d[is]] > 0.0)  {
                    // flux cannot exceed surface water available
                    vseep = fabs(q)*param->dtg*data->wcs[is];
                    if (q < 0.0 & vseep > data->dept[gmap->top2d[is]])    {
                        q = -data->dept[gmap->top2d[is]]/(param->dtg*data->wcs[is]);
                    }
                }
                else if (param->bctype_GW[5] == 2)   {
                    // infiltration from dry land is prohibited
                    if (q <= 0.0)  {q = 0.0;}
                    // evaporation and rainfall
                    if (data->qtop[gmap->top2d[is]] > 0.0 & data->wc[is] > data->wcr[is])
                    {q += data->qtop[gmap->top2d[is]];}
                    else if (data->qtop[gmap->top2d[is]] < 0.0)
                    {q += data->qtop[gmap->top2d[is]];}
                }
                else if (param->bctype_GW[5] == 0)   {q = 0.0;}
            }
            else {
                if (param->bctype_GW[5] == 2)   {q = data->qtop[gmap->top2d[is]];}
                else if (param->bctype_GW[5] == 0)  {q = 0.0;}
            }
        }
        else if (gmap->kk[ic] == param->nz-1)   {
            if (param->bctype_GW[4] == 0)   {q = 0.0;}
            else if (param->bctype_GW[4] == 2)  {q = data->qbot;}
            else if (param->bctype_GW[5] == 3)  {q = -kface * visc * rho;}
        }
    }
    return q;
}

// >>>>> Residual for the Richards equation <<<<<
double compute_residual(Data **data, Map *gmap, int ii, char *jaco, Config *param, int irank)
{
    double wcp, resi, dhc, dhm, dhp, dh, dh_min=1e-8;
    // get step size
    if (abs((*data)->h[ii]) < 1.0)  {
        if ((*data)->h[ii] > 0) {dh = dh_min;}
        else {dh = -dh_min;}
    }
    else {dh = (*data)->h[ii] * dh_min;}
    dhc = 0.0;
    if (strcmp(jaco, "ct") == 0)    {dhc = dh;}
    wcp = compute_wch(*data, (*data)->h[ii]+dhc, ii, param);
    resi = param->Ss*wcp*((*data)->h[ii]+dhc-(*data)->hn[ii])/(*data)->wcs[ii]/param->dtg;
    resi += (wcp-(*data)->wcn[ii])/param->dtg;
    // >>> flux in x
    dhp = 0.0;  dhm = 0.0; dhc = 0.0;
    if (strcmp(jaco, "ct") == 0)    {dhc = dh;}
    else if (strcmp(jaco, "xp") == 0)   {dhp = dh;}
    else if (strcmp(jaco, "xm") == 0)   {dhm = dh;}
    resi -= darcy_flux(*data, gmap, ii, dhc, dhp, "x", "none", param, irank)/param->dx;
    if (gmap->ii[ii] == 0)  {resi += darcy_flux(*data, gmap, ii, dhm, dhc, "x", "back", param, irank)/param->dx;}
    else {resi += darcy_flux(*data, gmap, gmap->iMjckc[ii], dhm, dhc, "x", "none", param, irank)/param->dx;}
    // >>> flux in y
    dhp = 0.0;  dhm = 0.0; dhc = 0.0;
    if (strcmp(jaco, "ct") == 0)    {dhc = dh;}
    else if (strcmp(jaco, "yp") == 0)   {dhp = dh;}
    else if (strcmp(jaco, "ym") == 0)   {dhm = dh;}
    resi -= darcy_flux(*data, gmap, ii, dhc, dhp, "y", "none", param, irank)/param->dy;
    if (gmap->jj[ii] == 0)  {resi += darcy_flux(*data, gmap, ii, dhm, dhc, "y", "back", param, irank)/param->dy;}
    else {resi += darcy_flux(*data, gmap, gmap->icjMkc[ii], dhm, dhc, "y", "none", param, irank)/param->dy;}
    // >>> flux in z
    dhp = 0.0;  dhm = 0.0; dhc = 0.0;
    if (strcmp(jaco, "ct") == 0)    {dhc = dh;}
    else if (strcmp(jaco, "zp") == 0)   {dhp = dh;}
    else if (strcmp(jaco, "zm") == 0)   {dhm = dh;}
    resi -= darcy_flux(*data, gmap, ii, dhc, dhp, "z", "none", param, irank) / gmap->dz3d[ii];
    if (gmap->kk[ii] == 0)  {resi += darcy_flux(*data, gmap, ii, dhm, dhc, "z", "back", param, irank)/ gmap->dz3d[ii];}
    else {resi += darcy_flux(*data, gmap, gmap->icjckM[ii], dhm, dhc, "z", "none", param, irank)/ gmap->dz3d[ii];}

    return resi;
}

double compute_jacobian_fd(Data **data, Map *gmap, int ii, char* axis, Config *param, int irank)
{
    double jaco, dh, dh_min=1e-8;
    // get step size
    if (abs((*data)->h[ii]) < 1.0)  {
        if ((*data)->h[ii] > 0) {dh = dh_min;}
        else {dh = -dh_min;}
    }
    else {dh = (*data)->h[ii] * dh_min;}
    if (strcmp(axis, "ct") == 0)    {
        jaco = (compute_residual(data, gmap, ii, "ct", param, irank) - compute_residual(data, gmap, ii, "none", param, irank)) / dh;
    }
    else if (strcmp(axis, "xp") == 0)  {
        jaco = (compute_residual(data, gmap, ii, "xp", param, irank) - compute_residual(data, gmap, ii, "none", param, irank)) / dh;
    }
    else if (strcmp(axis, "xm") == 0)  {
        jaco = (compute_residual(data, gmap, ii, "xm", param, irank) - compute_residual(data, gmap, ii, "none", param, irank)) / dh;
    }
    else if (strcmp(axis, "yp") == 0)  {
        jaco = (compute_residual(data, gmap, ii, "yp", param, irank) - compute_residual(data, gmap, ii, "none", param, irank)) / dh;
    }
    else if (strcmp(axis, "ym") == 0)  {
        jaco = (compute_residual(data, gmap, ii, "ym", param, irank) - compute_residual(data, gmap, ii, "none", param, irank)) / dh;
    }
    else if (strcmp(axis, "zp") == 0)  {
        jaco = (compute_residual(data, gmap, ii, "zp", param, irank) - compute_residual(data, gmap, ii, "none", param, irank)) / dh;
    }
    else if (strcmp(axis, "zm") == 0)  {
        jaco = (compute_residual(data, gmap, ii, "zm", param, irank) - compute_residual(data, gmap, ii, "none", param, irank)) / dh;
    }
    return jaco;
}



// >>>>> Compute water content from h using water retention curve <<<<<
double compute_wch(Data *data, double h, int ii, Config *param)
{
    double wc, m, s, wcs, wcm;
    // h = data->h[ii];
    if (param->use_vg == 1) {
        m = 1.0 - 1.0/data->vgn[ii];
        if (param->use_mvg == 1)
        {wcm = data->wcr[ii] + (data->wcs[ii]-data->wcr[ii])*pow((1.0 + pow(fabs(param->aev)*data->vga[ii],data->vgn[ii])), m);}
        else
        {wcm = data->wcs[ii];}
        s = pow(1.0 + pow(fabs(data->vga[ii]*h), data->vgn[ii]), -m);

        if (h > param->aev)  {wc = data->wcs[ii];}
        else    {wc = data->wcr[ii] + (wcm-data->wcr[ii]) * s;}
    }
    else
    {wc = data->wcr[ii] + (data->wcs[ii]-data->wcr[ii]) * exp(0.1634*h);}

    if (wc > data->wcs[ii])    {wc = data->wcs[ii];}
    else if (wc < data->wcr[ii])   {wc = data->wcr[ii];}

    // Warrick 1971
    // if (h < 0.0)
    // {
    //     if (h <= -0.29484)  {wc = 0.6829 - 0.09524*log(fabs(h)*100.0);}
    //     else    {wc = 0.4531 - 0.02732*log(fabs(h)*100.0);}
    //     if (h > param->aev)  {wc = data->wcs[ii];}
    //     if (wc > data->wcs[ii])    {wc = data->wcs[ii];}
    //     else if (wc < data->wcr[ii])   {wc = data->wcr[ii];}
    // }

    return wc;
}

// >>>>> Compute h from water content using water retention curve <<<<<
double compute_hwc(Data *data, int ii, Config *param)
{
    double wc, wcs, h, m, wcm, eps = 1e-7;
    wc = data->wc[ii];

    if (param->use_vg == 1)  {
        m = 1.0 - 1.0/data->vgn[ii];
        if (param->use_mvg == 1)
        {wcm = data->wcr[ii] + (data->wcs[ii]-data->wcr[ii])*pow((1.0 + pow(fabs(param->aev)*data->vga[ii],data->vgn[ii])), m);}
        else
        {wcm = data->wcs[ii];}

        if (wc - data->wcr[ii] < eps)  {wc = data->wcr[ii] + eps;}
        if (wc < data->wcs[ii])
        {
            h = -(1.0/data->vga[ii]) *
                pow(pow((wcm - data->wcr[ii])/(wc - data->wcr[ii]),(1.0/m)) - 1.0,(1.0/data->vgn[ii]));
        }
        else    {h = 0.0;}
    }
    else
    {h = log((wc - data->wcr[ii])/(data->wcs[ii] - data->wcr[ii])) / 0.1634;}

    // Warrick 1971
    // if (wc < 0.38-eps)
    // {
    //     if (wc < 0.36/0.38)    {h = -exp((2.95 - wc/0.38)/0.25) / 9810.0;}
    //     else    {h = -exp((1.52 - wc/0.38)/0.07) / 9810.0;}
    // }

    return h;
}

// >>>>> Compute specific capacity <<<<<
double compute_ch(Data *data, int ii, Config *param)
{
    double h, wcs, c, m, wcm, deno, nume;
    h = data->h[ii];
    m = 1.0 - 1.0/data->vgn[ii];
    if (param->use_mvg == 1)
    {wcm = data->wcr[ii] + (data->wcs[ii]-data->wcr[ii])*pow((1.0 + pow(fabs(param->aev)*data->vga[ii],data->vgn[ii])), m);}
    else
    {wcm = data->wcs[ii];}

    nume = data->vga[ii]*data->vgn[ii]*m * (wcm-data->wcr[ii]) *
        pow(fabs(data->vga[ii]*h),data->vgn[ii]-1);
    deno = pow((1.0 + pow(fabs(data->vga[ii]*h),data->vgn[ii])), m+1);
    c = nume / deno;

    // Warrick 1971
    // if (h < 0.0)
    // {
    //     if (h <= -0.29484)  {c = 0.09524 / (fabs(h));}
    //     else    {c = 0.02732 / (fabs(h));}
    // }

    if (param->use_mvg == 1)
    {if (h > param->aev)  {c = 0.0;}}
    else
    {if (h > 0.0)  {c = 0.0;}}
    return c;
}
// >>>>> Compute hydraulic conductivity <<<<<
double compute_K(Data *data, double h, double Ks, int ii, Config *param)
{
    double m, wcs, Keff, s, wcm, nume, deno, ratio;
    if (param->use_vg == 1) {
        m = 1.0 - 1.0/data->vgn[ii];
        s = pow(1.0 + pow(fabs(data->vga[ii]*h), data->vgn[ii]), -m);
        if (param->use_mvg == 1)
        {wcm = data->wcr[ii] + (data->wcs[ii]-data->wcr[ii])*pow((1.0 + pow(fabs(param->aev)*data->vga[ii],data->vgn[ii])), m);}
        else
        {wcm = data->wcs[ii];}
        nume = 1.0-pow(1.0-pow(s*(data->wcs[ii]-data->wcr[ii])/(wcm-data->wcr[ii]),1.0/m),m);
        deno = 1.0-pow(1.0-pow((data->wcs[ii]-data->wcr[ii])/(wcm-data->wcr[ii]),1.0/m),m);

        if (param->use_mvg == 1)
        {
            if (deno == 0.0)    {Keff = Ks;}
            else    {Keff = Ks * pow(s,0.5) * pow(nume/deno, 2.0);}
        }
        else
        {Keff = Ks * pow(s,0.5) * pow(1-pow(1-pow(s,1.0/m),m), 2.0);}
    }
    else {Keff = Ks * exp(0.1634 * h);}

    if (Keff > Ks)  {Keff = Ks;}
    if (param->use_mvg == 1)    {if (h > param->aev)  {Keff = Ks;}}
    else    {if (h > 0.0)  {Keff = Ks;}}
    return Keff;
}

// // >>>>> Compute hydraulic conductivity <<<<<
// double compute_K(Data *data, double *Ksat, int ii, Config *param)
// {
//     double m, wcs, Keff, s, h, Ks, wcm, nume, deno, ratio;
//     Ks = Ksat[ii];
//     h = data->h[ii];
//
//     if (param->use_vg == 1) {
//         m = 1.0 - 1.0/data->vgn[ii];
//         s = pow(1.0 + pow(fabs(data->vga[ii]*h), data->vgn[ii]), -m);
//         if (param->use_mvg == 1)
//         {wcm = data->wcr[ii] + (data->wcs[ii]-data->wcr[ii])*pow((1.0 + pow(fabs(param->aev)*data->vga[ii],data->vgn[ii])), m);}
//         else
//         {wcm = data->wcs[ii];}
//         nume = 1.0-pow(1.0-pow(s*(data->wcs[ii]-data->wcr[ii])/(wcm-data->wcr[ii]),1.0/m),m);
//         deno = 1.0-pow(1.0-pow((data->wcs[ii]-data->wcr[ii])/(wcm-data->wcr[ii]),1.0/m),m);
//
//         if (param->use_mvg == 1)
//         {
//             if (deno == 0.0)    {Keff = Ks;}
//             else    {Keff = Ks * pow(s,0.5) * pow(nume/deno, 2.0);}
//         }
//         else
//         {Keff = Ks * pow(s,0.5) * pow(1-pow(1-pow(s,1.0/m),m), 2.0);}
//     }
//     else {Keff = Ks * exp(0.1634 * h);}
//
//     // Warrick 1971
//     // if (h < 0.0)
//     // {
//     //     if (h < -0.29484)   {Keff = 19.34 * 1e4 * pow(fabs(h)*100.0, -3.4095) / (1e2 * 86400.0);}
//     //     else    {Keff = 516.8 * pow(fabs(h)*100.0, -0.97814) / (1e2 * 86400.0);}
//     // }
//
//     if (Keff > Ks)  {Keff = Ks;}
//     if (param->use_mvg == 1)
//     {if (h > param->aev)  {Keff = Ks;}}
//     else
//     {if (h > 0.0)  {Keff = Ks;}}
//
//     return Keff;
// }

// >>>>> Compute dKdwc for adaptive time stepping <<<<<
double compute_dKdwc(Data *data, double *Ksat, int ii, Config *param)
{
    double term1, term2, term0, m, s, dKdwc, wc, Ks, wcm, c1, c2;
    Ks = Ksat[ii];
    wc = data->wc[ii];
    if (param->use_mvg == 1)
    {wcm = data->wcr[ii] + (data->wcs[ii]-data->wcr[ii])*pow((1.0 + pow(fabs(param->aev)*data->vga[ii],data->vgn[ii])), m);}
    else
    {wcm = data->wcs[ii];}

    // use the lambda limiter
    if (wc > 0.9999 * data->wcs[ii] & wc < data->wcs[ii])
    {wc = 0.9999 * data->wcs[ii];}
    m = 1.0 - 1.0/data->vgn[ii];
    s = (wc - data->wcr[ii]) / (data->wcs[ii] - data->wcr[ii]);
    if (param->use_mvg == 0)
    {
        // calculate dKdwc
        term0 = pow(1.0 - pow(s,1.0/m), m);
        term1 = 0.5 * Ks * pow(s,-0.5) * (1.0 - term0) * (1.0 - term0);
        term2 = 2.0 * Ks * pow(s,(2.0-m)/(2.0*m)) * (1.0 - term0) * pow(1.0 - pow(s,1.0/m), m-1);
    }
    else
    {
        c2 = (data->wcs[ii] - data->wcr[ii])/(wcm - data->wcr[ii]);
        c1 = 1.0 / pow(1.0 - pow((1.0 - pow(c2,1.0/m)),m), 2.0);
        term0 = pow(1.0 - pow(c2*s,1.0/m), m);
        term1 = 0.5 * c1 * Ks * pow(s,-0.5) * (1.0 - term0) * (1.0 - term0);
        term2 = 2.0 * c1 * Ks * c2 * pow(s,0.5) * pow(c2*s,1.0/m-1.0) * (1.0 - term0) * pow(1.0 - pow(c2*s,1.0/m), m-1);
    }
    dKdwc = (term1 + term2) / (data->wcs[ii] - data->wcr[ii]);
    return dKdwc;
}

// >>>>> Compute dKdh for the Jacobian <<<<<
double compute_dKdh(Data *data, double *Ksat, int ii, Config *param)
{
    double dkdh, k1, k2, h, m, s, wcm, deno, nume, dh=1e-3;
    h = data->h[ii]+dh;
    k1 = compute_K(data, h, Ksat[ii], ii, param);
    if (param->use_vg == 1) {
        m = 1.0 - 1.0/data->vgn[ii];
        s = pow(1.0 + pow(fabs(data->vga[ii]*h), data->vgn[ii]), -m);

        if (param->use_mvg == 1)
        {wcm = data->wcr[ii] + (data->wcs[ii]-data->wcr[ii])*pow((1.0 + pow(fabs(param->aev)*data->vga[ii],data->vgn[ii])), m);}
        else
        {wcm = data->wcs[ii];}
        nume = 1.0-pow(1.0-pow(s*(data->wcs[ii]-data->wcr[ii])/(wcm-data->wcr[ii]),1.0/m),m);
        deno = 1.0-pow(1.0-pow((data->wcs[ii]-data->wcr[ii])/(wcm-data->wcr[ii]),1.0/m),m);

        if (param->use_mvg == 1)
        {
            if (deno == 0.0)    {k2 = Ksat[ii];}
            else    {k2 = Ksat[ii] * pow(s,0.5) * pow(nume/deno, 2.0);}
        }
        else
        {k2 = Ksat[ii] * pow(s,0.5) * pow(1-pow(1-pow(s,1.0/m),m), 2.0);}
    }
    else {k2 = Ksat[ii] * exp(0.1634 * h);}
    if (k2 > Ksat[ii] || h >= 0.0)  {k2 = Ksat[ii];}
    if (data->h[ii] > -dh)    {dkdh = 0.0;}
    else {dkdh = 0.5 * (k2 - k1) / dh;}
    return dkdh;
}

// >>>>> Compute dwcdh for the Jacobian <<<<<
double compute_dwcdh(Data *data, int ii, Config *param)
{
    double dwcdh, wc1, wc2, h, m, s, wcm, dh=1e-3;
    h = data->h[ii]+dh;
    wc1 = compute_wch(data, h, ii, param);

    if (param->use_vg == 1) {
        m = 1.0 - 1.0/data->vgn[ii];
        if (param->use_mvg == 1)
        {wcm = data->wcr[ii] + (data->wcs[ii]-data->wcr[ii])*pow((1.0 + pow(fabs(param->aev)*data->vga[ii],data->vgn[ii])), m);}
        else
        {wcm = data->wcs[ii];}
        s = pow(1.0 + pow(fabs(data->vga[ii]*h), data->vgn[ii]), -m);
        if (h > param->aev)  {wc2 = data->wcs[ii];}
        else    {wc2 = data->wcr[ii] + (wcm-data->wcr[ii]) * s;}
        wc2 = data->wcr[ii] + (data->wcs[ii]-data->wcr[ii]) * s;
    }
    else {wc2 = data->wcr[ii] + (data->wcs[ii]-data->wcr[ii]) * exp(0.1634*h);}
    if (data->h[ii] > -dh)    {dwcdh = 0.0;}
    else {dwcdh = (wc2 - wc1) / dh;}
    return dwcdh;
}

// >>>>> 2nd-order TVD scheme with superbee limiter
double tvd_superbee(double sp, double sc, double sm, double u, double delta, double dt, Config *param)
{
    double phi, r, coef, r1, r2;
    coef = fabs(u) * dt / delta;
    r1 = 1.0;
    r2 = 2.0;
    phi = 0.0;
    if (sp != sc)
    {
        r = (sc - sm) / (sp - sc);
        if (2.0*r < 1.0)  {r1 = 2.0*r;}
        if (r < 2.0)    {r2 = r;}
        if (r1 > 0.0 | r2 > 0.0)
        {
            if (r1 > r2)    {phi = r1;}
            else    {phi = r2;}
        }
    }
    return sc + 0.5*phi*(1.0-coef)*(sp - sc);
}
