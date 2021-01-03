// Head file for configuration.c

#ifndef CONFIGURATION_H
#define CONFIGURATION_H

typedef struct Config
{
    // Directory
    char finput[100], foutput[100], sim_id[6];
    // Domain Geometry
    int NX, NY, NZ, nx, ny, nz, use_mpi, mpi_nx, mpi_ny;
    int n2ci, n2ct, N2CI, n3ci, n3ct, N3CI;
    double dx, dy, dz, botZ, dz_incre;
    // Time Control
    int NT;
    double dt, Tend, dt_out, *dt_root;
    // Bathymetry
    int bath_file;
    // parameters
    double grav, viscx, viscy, rhoa, rhow;
    double min_dept, wtfh, hD, manning;
    // Surface water
    int sim_shallowwater, eta_file, uv_file, difuwave;
    double init_eta, *init_tide, *init_inflow, q_evap, q_rain;
    int *bctype_SW, *inflow_locX, *inflow_locY;
    int n_tide, n_inflow, *tide_locX, *tide_locY;
    int *tide_dat_len, *inflow_dat_len;
    int *tide_file, *inflow_file, evap_file, rain_file, evap_model;
    int sim_wind, wind_file;
    double init_winddir, init_windspd, Cw, CwT, north_angle;
    // subgrid model
    int use_subgrid;
    // Groundwater
    int sim_groundwater, dt_adjust, use_corrector, post_allocate, use_mvg, use_full3d;
    double init_h, init_wc, init_wt_abs, init_wt_rel, qtop, qbot, htop, hbot, aev;
    double dt_max, dt_min, Co_max, Ksx, Ksy, Ksz, Ss, wcr, wcs, soil_a, soil_n;
    int *bctype_GW, h_file, wc_file;
    // Scalar
    int n_scalar, *scalar_surf_file, *scalar_tide_datlen, *scalar_tide_file, *scalar_inflow_datlen, *scalar_inflow_file;
    int *scalar_subs_file, baroclinic;
    double difux, difuy, difuz, disp_lat, disp_lon;
    double *init_s_surf, *init_s_subs;
    double *s_tide, *s_inflow;

}Config;

#endif


void read_input(Config **param);
