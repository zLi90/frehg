// Header file for initialize.c
#include"map.h"
#include"configuration.h"

#ifndef INITIALIZE_H
#define INITIALIZE_H



typedef struct Data
{
    // bathymetry related fields
    double *bottom, *bottom_root, *offset, *bottomXP, *bottomYP;
    // surface domain
    int *reset_seepage;
    double *uu, *uy, *vv, *vx, *eta, *etan, *dept, *deptx, *depty, *qseepage;
    double *uu_root, *vv_root, *eta_root, *dept_root, *seep_root;
    double *uu_out, *vv_out, *eta_out, *dept_out, *seep_out;
    double *Fu, *Fv, *Ex, *Ey, *Dx, *Dy, *CDx, *CDy, *wtfx, *wtfy, *cflx, *cfly, *cfl_active;
    double *Vs, *Vsx, *Vsy, *Asx, *Asy, *Asz, *Aszx, *Aszy;
    // subsurface domain
    double *h, *hn, *hp, *hwc, *wc, *wcn, *wcp, *wch, *h_root, *wc_root, *dh6, *rsplit;
    double *vloss, *vloss_root, *room, qtop, qbot, hbot, htop;
    double *Kx, *Ky, *Kz, *qx, *qy, *qz, *qx_root, *qy_root, *qz_root, *Vg, *ch;
    double *h_out, *wc_out, *qx_out, *qy_out, *qz_out;
    double *wcs, *wcr, *vga, *vgn, *Ksz, *Ksx, *Ksy;
    double *t_out, *qbc;
    int *repeat;
    // linear system
    double *Sct, *Srhs, *Sxp, *Syp, *Sxm, *Sym;
    double *Gct, *Grhs, *Gxp, *Gxm, *Gyp, *Gym, *Gzp, *Gzm;
    // boundary conditions
    double *tide, *tide1, *t_tide1, *tide2, *t_tide2;
    double *rain, *evap, *q_rain, *t_rain, *q_evap, *t_evap, *rain_sum;
    double *inflow, *t_inflow;
}Data;

#endif

void init(Data **data, Map **smap, Map **gmap, Config **param, int irank, int nrank);
void init_domain(Config **param);
void init_Data(Data **data, Config *param);
void ic_surface(Data **data, Map *smap, Config *param, int irank, int nrank);
void bc_surface(Data **data, Map *smap, Config *param, int irank);
void update_depth(Data **data, Map *smap, Config *param, int irank);
void read_bathymetry(Data **data, Config *param, int irank, int nrank);
void boundary_bath(Data **data, Map *smap, Config *param, int irank, int nrank);
void ic_subsurface(Data **data, Map *gmap, Config *param);
