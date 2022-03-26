#include"configuration.h"
#include"initialize.h"
#include"map.h"
#include"mpifunctions.h"

#ifndef SUBROUTINES_H
#define SUBROUTINES_H

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


#endif
