// read configuration
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>


#include"utility.h"


void read_input(Config **param);


// >>>>> Read all user settings <<<<<
void read_input(Config **param)
{
    *param = malloc(sizeof(Config));
    // Directory
    strcpy((*param)->finput, read_one_input("finput", "input"));
    strcpy((*param)->foutput, read_one_input("foutput", "input"));
    strcpy((*param)->sim_id, read_one_input("sim_id", "input"));


    // Domain geometry
    (*param)->NX = (int) read_one_input_double("NX", "input");
    (*param)->NY = (int) read_one_input_double("NY", "input");
    (*param)->botZ = read_one_input_double("botZ", "input");
    (*param)->dx = read_one_input_double("dx", "input");
    (*param)->dy = read_one_input_double("dy", "input");
    (*param)->dz = read_one_input_double("dz", "input");
    (*param)->dz_incre = read_one_input_double("dz_incre", "input");
    (*param)->use_mpi = (int) read_one_input_double("use_mpi", "input");
    (*param)->mpi_nx = (int) read_one_input_double("mpi_nx", "input");
    (*param)->mpi_ny = (int) read_one_input_double("mpi_ny", "input");

    // Time control
    (*param)->dt = read_one_input_double("dt", "input");
    (*param)->Tend = read_one_input_double("Tend", "input");
    (*param)->NT = (int) read_one_input_double("NT", "input");
    (*param)->dt_out = read_one_input_double("dt_out", "input");
    if ((*param)->use_mpi == 1)
    {(*param)->dt_root = malloc((*param)->mpi_nx*(*param)->mpi_ny*sizeof(double));}

    // Bathymetry
    (*param)->bath_file = (int) read_one_input_double("bath_file", "input");

    // Parameters
    (*param)->min_dept = read_one_input_double("min_dept", "input");
    (*param)->wtfh = read_one_input_double("wtfh", "input");
    (*param)->hD = read_one_input_double("hD", "input");
    (*param)->manning = read_one_input_double("manning", "input");
    (*param)->grav = read_one_input_double("grav", "input");
    (*param)->viscx = read_one_input_double("viscx", "input");
    (*param)->viscy = read_one_input_double("viscy", "input");
    (*param)->rhoa = read_one_input_double("rhoa", "input");
    (*param)->rhow = read_one_input_double("rhow", "input");

    // Wind parameters
    (*param)->sim_wind = (int) read_one_input_double("sim_wind", "input");
    (*param)->wind_file = (int) read_one_input_double("wind_file", "input");
    (*param)->init_windspd = read_one_input_double("init_windspd", "input");
    (*param)->init_winddir = read_one_input_double("init_winddir", "input");
    (*param)->Cw = read_one_input_double("Cw", "input");
    (*param)->CwT = read_one_input_double("CwT", "input");
    (*param)->north_angle = read_one_input_double("north_angle", "input");

    // Shallow water solver
    (*param)->sim_shallowwater = (int) read_one_input_double("sim_shallowwater", "input");
    (*param)->difuwave = (int) read_one_input_double("difuwave", "input");
    (*param)->bctype_SW = read_one_input_array("bctype_SW", "input", 4);

    (*param)->init_eta = read_one_input_double("init_eta", "input");
    (*param)->eta_file = (int) read_one_input_double("eta_file", "input");
    (*param)->uv_file = (int) read_one_input_double("uv_file", "input");

    (*param)->n_tide = (int) read_one_input_double("n_tide", "input");
    (*param)->tide_file = read_one_input_array("tide_file", "input", (*param)->n_tide);
    (*param)->init_tide = read_one_input_array_double("init_tide", "input", (*param)->n_tide);
    (*param)->tide_locX = read_one_input_array("tide_locX", "input", 2*(*param)->n_tide);
    (*param)->tide_locY = read_one_input_array("tide_locY", "input", 2*(*param)->n_tide);
    (*param)->tide_dat_len = read_one_input_array("tide_dat_len", "input", (*param)->n_tide);

    (*param)->evap_file = (int) read_one_input_double("evap_file", "input");
    (*param)->evap_model = (int) read_one_input_double("evap_model", "input");
    (*param)->q_evap = read_one_input_double("q_evap", "input");
    (*param)->rain_file = (int) read_one_input_double("rain_file", "input");
    (*param)->q_rain = read_one_input_double("q_rain", "input");

    (*param)->n_inflow = (int) read_one_input_double("n_inflow", "input");
    (*param)->inflow_locX = read_one_input_array("inflow_locX", "input", 2*(*param)->n_inflow);
    (*param)->inflow_locY = read_one_input_array("inflow_locY", "input", 2*(*param)->n_inflow);
    (*param)->inflow_file = read_one_input_array("inflow_file", "input", (*param)->n_inflow);
    (*param)->inflow_dat_len = read_one_input_array("inflow_dat_len", "input", (*param)->n_inflow);
    (*param)->init_inflow = read_one_input_array_double("init_inflow", "input", (*param)->n_inflow);

    // subgrid model
    (*param)->use_subgrid = (int) read_one_input_double("use_subgrid", "input");
    (*param)->r_sub = read_one_input_double("r_sub", "input");
    (*param)->eta_sub_min = read_one_input_double("eta_sub_min", "input");
    (*param)->eta_sub_max = read_one_input_double("eta_sub_max", "input");
    (*param)->deta_sub = read_one_input_double("deta_sub", "input");

    // Groundwater solver
    (*param)->sim_groundwater = (int) read_one_input_double("sim_groundwater", "input");
    (*param)->use_full3d = (int) read_one_input_double("use_full3d", "input");
    (*param)->dt_adjust = (int) read_one_input_double("dt_adjust", "input");
    (*param)->use_corrector = (int) read_one_input_double("use_corrector", "input");
    (*param)->post_allocate = (int) read_one_input_double("post_allocate", "input");
    (*param)->use_mvg = (int) read_one_input_double("use_mvg", "input");
    (*param)->aev = read_one_input_double("aev", "input");
    (*param)->dt_max = read_one_input_double("dt_max", "input");
    (*param)->dt_min = read_one_input_double("dt_min", "input");
    (*param)->Co_max = read_one_input_double("Co_max", "input");
    (*param)->Ksx = read_one_input_double("Ksx", "input");
    (*param)->Ksy = read_one_input_double("Ksy", "input");
    (*param)->Ksz = read_one_input_double("Ksz", "input");
    (*param)->Ss = read_one_input_double("Ss", "input");
    (*param)->wcs = read_one_input_double("wcs", "input");
    (*param)->wcr = read_one_input_double("wcr", "input");
    (*param)->soil_a = read_one_input_double("soil_a", "input");
    (*param)->soil_n = read_one_input_double("soil_n", "input");
    // groundwater initial condition
    (*param)->init_wc = read_one_input_double("init_wc", "input");
    (*param)->init_h = read_one_input_double("init_h", "input");
    (*param)->init_wt_abs = read_one_input_double("init_wt_abs", "input");
    (*param)->init_wt_rel = read_one_input_double("init_wt_rel", "input");
    (*param)->h_file = (int) read_one_input_double("h_file", "input");
    (*param)->wc_file = (int) read_one_input_double("wc_file", "input");
    (*param)->qtop = read_one_input_double("qtop", "input");
    (*param)->qbot = read_one_input_double("qbot", "input");
    (*param)->qyp = read_one_input_double("qyp", "input");
    (*param)->qym = read_one_input_double("qym", "input");
    (*param)->htop = read_one_input_double("htop", "input");
    (*param)->hbot = read_one_input_double("hbot", "input");
    // groundwater boundary condition
    (*param)->bctype_GW = read_one_input_array("bctype_GW", "input", 6);

    // Scalar transport
    (*param)->n_scalar = (int) read_one_input_double("n_scalar", "input");
    (*param)->baroclinic = (int) read_one_input_double("baroclinic", "input");
    (*param)->superbee = (int) read_one_input_double("superbee", "input");

    (*param)->scalar_surf_file = read_one_input_array("scalar_surf_file", "input", (*param)->n_scalar);
    (*param)->scalar_tide_file = read_one_input_array("scalar_tide_file", "input", (*param)->n_scalar*(*param)->n_tide);
    (*param)->scalar_tide_datlen = read_one_input_array("scalar_tide_datlen", "input", (*param)->n_scalar*(*param)->n_tide);
    (*param)->scalar_inflow_file = read_one_input_array("scalar_inflow_file", "input", (*param)->n_scalar*(*param)->n_inflow);
    (*param)->scalar_inflow_datlen = read_one_input_array("scalar_inflow_datlen", "input", (*param)->n_scalar*(*param)->n_inflow);

    (*param)->scalar_subs_file = read_one_input_array("scalar_subs_file", "input", (*param)->n_scalar);

    (*param)->init_s_surf = read_one_input_array_double("init_s_surf", "input", (*param)->n_scalar);
    (*param)->init_s_subs = read_one_input_array_double("init_s_subs", "input", (*param)->n_scalar);
    (*param)->s_tide = read_one_input_array_double("s_tide", "input", (*param)->n_scalar*(*param)->n_tide);
    (*param)->s_inflow = read_one_input_array_double("s_inflow", "input", (*param)->n_scalar*(*param)->n_tide);
    (*param)->s_yp = read_one_input_array_double("s_yp", "input", (*param)->n_scalar);
    (*param)->s_ym = read_one_input_array_double("s_ym", "input", (*param)->n_scalar);
    (*param)->difux = read_one_input_double("difux", "input");
    (*param)->difuy = read_one_input_double("difuy", "input");
    (*param)->difuz = read_one_input_double("difuz", "input");
    (*param)->disp_lon = read_one_input_double("disp_lon", "input");
    (*param)->disp_lat = read_one_input_double("disp_lat", "input");

}
