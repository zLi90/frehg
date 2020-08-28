// Head file for configuration.c

#ifndef CONFIGURATION_H
#define CONFIGURATION_H

typedef struct Config
{
    // Directory
    char finput[50], foutput[50], sim_id[6];
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
    double grav, viscx, viscy;
    double min_dept, wtfh, hD, manning;
    // Surface water
    int sim_shallowwater, eta_file, uv_file;
    double init_eta, *init_tide, q_evap, q_rain;
    int *bctype_SW, *inflow_loc;
    int n_tide, n_tide1, n_tide2, n_inflow, *tide_locX, *tide_locY, tide_file, evap_file, rain_file;
    // subgrid model
    int use_subgrid;
    // Groundwater
    int sim_groundwater, dt_adjust, use_corrector, post_allocate, use_mvg, use_full3d;
    double init_h, init_wc, init_wt_abs, init_wt_rel, qtop, qbot, htop, hbot, aev;
    double dt_max, dt_min, Co_max, Ksx, Ksy, Ksz, Ss, wcr, wcs, soil_a, soil_n;
    int *bctype_GW;

}Config;

#endif


void read_input(Config **param);
