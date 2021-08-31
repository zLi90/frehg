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
    double *uu, *un, *uy, *vv, *vn, *vx, *eta, *etan, *dept, *deptx, *depty, *qseepage;
    double *uu_root, *vv_root, *un_root, *vn_root, *eta_root, *dept_root, *seep_root;
    double *uu_out, *vv_out, *un_out, *vn_out, *eta_out, *dept_out, *seep_out;
    double *Fu, *Fv, *Ex, *Ey, *Dx, *Dy, *CDx, *CDy, *wtfx, *wtfy, *cflx, *cfly, *cfl_active;
    double *Vs, *Vsn, *Vflux, *Vsx, *Vsy, *Asx, *Asy, *Asz, *Aszx, *Aszy;
    // subgrid settings
    int *eta_ind;
    double *layers_sub;
    double *edges_root, *Vxp_sub_root, *Vxm_sub_root, *Vyp_sub_root, *Vym_sub_root;
    double *Axp_sub_root, *Axm_sub_root, *Ayp_sub_root, *Aym_sub_root, *Asz_sub_root;
    double **Vxp_sub, **Vxm_sub, **Vyp_sub, **Vym_sub;
    double **Axp_sub, **Axm_sub, **Ayp_sub, **Aym_sub, **Asz_sub;
    double *Vxp, *Vxm, *Vyp, *Vym, *Axp, *Axm, *Ayp, *Aym;
    // subsurface domain
    double *h, *hn, *hp, *hwc, *wc, *wcn, *wcp, *wch, *h_root, *wc_root, *dh6, *rsplit;
    double *vloss, *vloss_root, *room, *qtop, qbot, hbot, htop;
    double *Kx, *Ky, *Kz, *qx, *qy, *qz, *qx_root, *qy_root, *qz_root, *Vg, *Vgn, *Vgflux, *ch;
    double *h_out, *wc_out, *qx_out, *qy_out, *qz_out;
    double *wcs, *wcr, *vga, *vgn, *Ksz, *Ksx, *Ksy;
    double *t_out, *qbc;
    double *r_rho, *r_rhon, *r_visc, *r_rhoxp, *r_rhoyp, *r_rhozp, *r_viscxp, *r_viscyp, *r_visczp;
    int *repeat;
    // linear system
    double *Sct, *Srhs, *Sxp, *Syp, *Sxm, *Sym;
    double *Gct, *Grhs, *Gxp, *Gxm, *Gyp, *Gym, *Gzp, *Gzm;
    // boundary conditions
    int **tideloc, *tideloc_len, **inflowloc, *inflowloc_len;
    double **tide, **t_tide, *current_tide;
    double *rain, *evap, *q_rain, *t_rain, *q_evap, *t_evap, *rain_sum;
    double *rain_data, *evap_data, *current_rain, *current_evap;
    double **inflow, **t_inflow, *current_inflow;
    double *wind_dir, *wind_spd, *current_windspd, *current_winddir, *t_wind;
    // scalar
    double **s_surf, **sm_surf, **s_surf_root, **s_surf_out, **s_surfkP;
    double **s_subs, **sm_subs, **s_subs_root, **s_subs_out;
    double ***s_tide, ***t_s_tide, **current_s_tide;
    double ***s_inflow, ***t_s_inflow, **current_s_inflow;
    double **sseepage;
    double *Dxx, *Dxy, *Dxz, *Dyy, *Dyx, *Dyz, *Dzz, *Dzx, *Dzy;
}Data;

#endif

void init(Data **data, Map **smap, Map **gmap, Config **param, int irank, int nrank);
void init_domain(Config **param);
void init_Data(Data **data, Config *param);
void ic_surface(Data **data, Map *smap, Map *gmap, Config *param, int irank, int nrank);
void bc_surface(Data **data, Map *smap, Config *param, int irank);
void get_BC_location(int **loc, int *loc_len, Config *param, int irank, int n_bc, int *locX, int *locY);
void update_depth(Data **data, Map *smap, Config *param, int irank);
void read_bathymetry(Data **data, Config *param, int irank, int nrank);
void boundary_bath(Data **data, Map *smap, Config *param, int irank, int nrank);
void ic_subsurface(Data **data, Map *gmap, Config *param, int irank, int nrank);
void restart_subsurface(double *ic_array, char *fname, Config *param, int irank);
void init_subgrid(Data **data, Map *smap, Config *param, int irank, int nrank);
