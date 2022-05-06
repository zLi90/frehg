// The time stepping loop
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

#include "configuration.h"
#include "groundwater.h"
#include "initialize.h"
#include "map.h"
#include "mpifunctions.h"
#include "shallowwater.h"
#include "scalar.h"
#include "subroutines.h"
#include "utility.h"


void solve(Data **data, Map *smap, Map *gmap, Config *param, int irank, int nrank);
void get_current_bc(Data **data, Config *param, double t_current);
void get_evaprain(Data **data, Map *gmap, Config *param, double t_current);
double get_mass(Data **data, Map *gmap, Config *param);
void print_end_info(Data **data, Map *smap, Map *gmap, Config *param, int irank);

void solve(Data **data, Map *smap, Map *gmap, Config *param, int irank, int nrank)
{
    int t_save, tday, ii, jj, kk, tt = 1, sub_advance = 0;
    float t0, t1, tstep, dt_max, last_save = 0.0, t_current = 0.0, t_subsurface = 0.0;
    double max_CFLx, max_CFLy, max_CFL, *max_CFL_root, tot_mass;
    double dt_pseudo;
    // save initial condition
    dt_max = param->dt;
    mpi_print(" >>> Writing initial conditions into output files !", irank);
    write_output(data, gmap, param, 0, 0, irank);
    // begin time stepping
    mpi_print(" >>> Beginning Time loop !", irank);
    if (param->sim_groundwater == 0 | param->sync_coupling == 1)    {param->dtg = param->dt;}
    else {param->dtg = param->dt_min;}

    while (t_current < param->Tend)
    {
        if (irank == 0) {t0 = clock();}
        t_current += param->dt;
        if (param->sim_groundwater == 0 | param->sync_coupling == 1)    {t_subsurface = t_current;}
        else {if (sub_advance == 1) {t_subsurface += param->dtg;}}
        // else {if (t_current == param->dt | sub_advance == 1) {t_subsurface += param->dtg;}}
        // get boundary condition
        get_current_bc(data, param, t_current);
        get_evaprain(data, gmap, param, t_current);
        // execute solvers
        if (param->sim_shallowwater == 1)
        {solve_shallowwater(data, smap, gmap, param, irank, nrank);}

        if (param->sim_groundwater == 1)
        {
            for (ii = 0; ii < param->n3ci; ii++)    {
                if (gmap->istop[ii] == 1)   {
                    jj = gmap->icjckM[ii];
                    (*data)->rss[gmap->top2d[ii]] = 0.0;
                    if ((*data)->eta[gmap->top2d[ii]] - (*data)->bottom[gmap->top2d[ii]] > 0.0)   {
                        (*data)->rss[gmap->top2d[ii]] = fabs((param->dtg*(*data)->Kz[jj]*
                            ((*data)->eta[gmap->top2d[ii]] - (*data)->etan[gmap->top2d[ii]])/(param->dt*gmap->dz3d[ii])));
                    }
                }
            }


            if (param->sync_coupling == 1)
            {solve_groundwater(data, smap, gmap, param, irank, nrank);}
            else {
                sub_advance = 0;
                while (t_subsurface + param->dtg <= t_current)   {
                    t_subsurface += param->dtg;
                    solve_groundwater(data, smap, gmap, param, irank, nrank);
                    if (irank == 0) {printf("   >>> Subsurface executed with dt = %f \n",param->dtg);}
                }
                // if (t_subsurface <= t_current - param->dtg) {
                //     // solve_groundwater(data, smap, gmap, param, irank, nrank);
                //     // if (irank == 0) {printf("   >>> Subsurface executed with dt = %f \n",param->dtg);}
                //     while (t_subsurface + param->dtg <= t_current)   {
                //         t_subsurface += param->dtg;
                //         solve_groundwater(data, smap, gmap, param, irank, nrank);
                //         if (irank == 0) {printf("   >>> Subsurface executed with dt = %f \n",param->dtg);}
                //     }
                //     sub_advance = 1;
                // }
                // else {
                //     if (param->sim_shallowwater == 1)   {
                //         dt_pseudo = (t_current - t_subsurface + param->dtg) / param->dtg;
                //         // dt_pseudo = t_current - t_subsurface + param->dtg;
                //         // pseudo_seepage(data, gmap, param, dt_pseudo);
                //     }
                // }
            }
            // param->dtg = param->dt_max;
        }
        if (param->sim_shallowwater == 1)
        {shallowwater_velocity(data, smap, gmap, param, irank, nrank);}
        // check CFL number
        max_CFLx = getMax((*data)->cflx, param->n2ci);
        max_CFLy = getMax((*data)->cfly, param->n2ci);
        if (max_CFLx > max_CFLy)    {max_CFL = max_CFLx;}
        else    {max_CFL = max_CFLy;}
        if (param->use_mpi == 1)
        {
            max_CFL_root = malloc(param->mpi_ny*param->mpi_nx*sizeof(double));
            mpi_gather_double(max_CFL_root, &max_CFL, 1, 0);
            if (irank == 0) {max_CFL = getMax(max_CFL_root, param->mpi_ny*param->mpi_nx);   free(max_CFL_root);}
            mpi_bcast_double(&max_CFL, 1, 0);
        }

        // scalar transport
        if (param->n_scalar > 0 & param->sim_shallowwater == 1)
        {for (kk = 0; kk < param->n_scalar; kk++)    {scalar_shallowwater(data, smap, param, irank, nrank, kk);}}
        if (param->n_scalar > 0 & param->sim_groundwater == 1)
        {for (kk = 0; kk < param->n_scalar; kk++)    {scalar_groundwater(data, gmap, param, irank, nrank, kk);}}
        // reset rainfall
        if ((*data)->rain_sum[0] > param->min_dept)     {(*data)->rain_sum[0] = 0.0;}

        // model output
        if (fabs(t_current - last_save - param->dt_out) <= param->dt)
        {
            t_save = round(t_current / param->dt_out) * param->dt_out;
            last_save = t_save;
            write_output(data, gmap, param, t_save, 0, irank);
        }
        // output total mass
        tot_mass = get_mass(data, gmap, param);
        write_monitor_out(tot_mass, "mass", 0, param);
        // output monitored variables
        if (param->n_monitor > 0)
        {
            for (ii = 0; ii < param->n_monitor; ii++)
            {
                if (irank == (*data)->monitor_rank[ii])
                {
                    if (param->sim_shallowwater == 1)   {
                        write_monitor_out((*data)->eta[(*data)->monitor[ii]] - (*data)->offset[0], "surf", ii, param);
                        write_monitor_out((*data)->dept[(*data)->monitor[ii]], "depth", ii, param);
                        write_monitor_out((*data)->vv[(*data)->monitor[ii]], "vv", ii, param);
                        if (param->n_scalar > 0)
                        {
                            for (kk = 0; kk < param->n_scalar; kk++)
                            {if (kk == 0) {write_monitor_out((*data)->s_surf[kk][(*data)->monitor[ii]], "scalar1", ii, param);}}
                        }
                    }
                    if (param->sim_groundwater == 1)    {
                        write_monitor_out(1e7*(*data)->qss[(*data)->monitor[ii]], "seepage", ii, param);
                        write_monitor_out(1e7*(*data)->rss[(*data)->monitor[ii]], "r_async", ii, param);
                    }
                }
            }
        }
        // if (max_CFL > 1.0)  {param->dt = param->dt / 2.0;}
        // report time
        if (irank == 0)
        {
          append_to_file("timestep", t_current, param);
          t1 = clock();
          tstep = (t1 - t0)/(float)CLOCKS_PER_SEC;
          printf(" >>>>> Step %d (%f of %.1f sec) completed, new dt = %f, cost = %.4f sec... \n",tt,t_current,param->Tend,param->dt,tstep);
        }
        tt += 1;
    }
    print_end_info(data, smap, gmap, param, irank);
}

// >>>>> Get BC at the current time step
void get_current_bc(Data **data, Config *param, double t_current)
{
    int kk, ss, glob_ind, ind;
    // wind
    if (param->sim_wind == 1 & param->wind_file == 1)
    {
        if (param->wind_file == 1)
        {
            (*data)->current_windspd[0] = interp_bc((*data)->t_wind,(*data)->wind_spd,t_current,param->wind_dat_len);
            (*data)->current_winddir[0] = (*data)->wind_dir[0];
            if (t_current > 0.0)
            {
                ind = 0;
                while (t_current > (*data)->t_wind[ind])
                {
                    ind += 1;
                    if (ind >= param->wind_dat_len) {break;}
                }
                (*data)->current_winddir[0] = (*data)->wind_dir[ind-1];
            }
        }
        else
        {
            (*data)->current_windspd[0] = (*data)->wind_spd[0];
            (*data)->current_winddir[0] = (*data)->wind_dir[0];
        }
    }
    // tide
    if (param->n_tide > 0)
    {
        for (kk = 0; kk < param->n_tide; kk++)
        {
            if (param->tide_file[kk] == 0)
            {(*data)->current_tide[kk] = param->init_tide[kk] + (*data)->offset[0];}
            else
            {(*data)->current_tide[kk] = interp_bc((*data)->t_tide[kk],(*data)->tide[kk],t_current,param->tide_dat_len[kk])+(*data)->offset[0];}
        }
    }
    // inflow
    if (param->n_inflow > 0)
    {
        for (kk = 0; kk < param->n_inflow; kk++)
        {
            if (param->inflow_file[kk] == 0)
            {(*data)->current_inflow[kk] = param->init_inflow[kk];}
            else
            {(*data)->current_inflow[kk] = interp_bc((*data)->t_inflow[kk],(*data)->inflow[kk],t_current,param->inflow_dat_len[kk]);}

        }
    }
    // scalar
    if (param->n_scalar > 0 & param->n_tide > 0)
    {
        for (ss = 0; ss < param->n_scalar; ss++)
        {
            for (kk = 0; kk < param->n_tide; kk++)
            {
                glob_ind = ss * param->n_tide + kk;
                if (param->scalar_tide_file[glob_ind] == 0)
                {(*data)->current_s_tide[ss][kk] = param->s_tide[glob_ind];}
                else
                {
                    (*data)->current_s_tide[ss][kk] =
                        interp_bc((*data)->t_s_tide[ss][kk],(*data)->s_tide[ss][kk],t_current,param->scalar_tide_datlen[glob_ind]);
                }
            }
        }
    }
    if (param->n_scalar > 0 & param->n_inflow > 0)
    {
        for (ss = 0; ss < param->n_scalar; ss++)
        {
            for (kk = 0; kk < param->n_inflow; kk++)
            {
                glob_ind = ss * param->n_inflow + kk;
                if (param->scalar_inflow_file[glob_ind] == 0)
                {(*data)->current_s_inflow[ss][kk] = param->s_inflow[glob_ind];}
                else
                {
                    (*data)->current_s_inflow[ss][kk] =
                        interp_bc((*data)->t_s_inflow[ss][kk],(*data)->s_inflow[ss][kk],t_current,param->scalar_inflow_datlen[glob_ind]);
                }
            }
        }
    }
}

// >>>>> Update rainfall and evaporation flux
void get_evaprain(Data **data, Map *gmap, Config *param, double t_current)
{
    int ii, jj;
    double temp, velo, pres, humi, t_residual;
    double rhoa, rhow, esat, qsat, resi, alpha, qsuf, qsub;
    // rainfall
    // For now, only consider uniform rainfall, ZhiLi20201117
    if (param->rain_file == 0)
    {(*data)->rain[0] = param->q_rain;}
    else
    {
        (*data)->current_rain[0] = interp_bc((*data)->t_rain,(*data)->rain_data,t_current,param->rain_dat_len);
        (*data)->rain[0] = (*data)->current_rain[0];

    }
    // evaporation
    // read evap data
    if (param->evap_file == 1)
    {(*data)->current_evap[0] = interp_bc((*data)->t_evap,(*data)->evap_data,t_current,param->evap_dat_len);}
    else
    {(*data)->current_evap[0] = param->q_evap;}

    // remove evaporation at night
    // t_residual = (int)floor(t_current) % 86400;
    // if (t_residual < 6.0 * 3600 | t_residual > 18.0 * 3600)
    // {(*data)->current_evap[0] = 0.0;}

    // (*data)->current_evap[0] = (*data)->current_evap[0] / 2.0;

    // use the aerodynamic evap model
    if (param->sim_groundwater == 1)
    {
        for (ii = 0; ii < param->n2ci; ii++)
        {(*data)->qtop[ii] = param->qtop;}

        if (param->evap_model == 0)
        {
            for (ii = 0; ii < param->n2ci; ii++)
            {(*data)->evap[ii] = (*data)->current_evap[0];}
        }
        else if (param->evap_model == 1)
        {
            if (param->sim_wind == 0)   {velo = 1.0;}
            else    {velo = (*data)->wind_spd[0];}   //if (velo>5.0)   {velo=5.0;}
            // else    {velo = (*data)->wind_spd[0];}
            temp = 20.0;
            pres = 101.325;
            humi = 0.0029;
            rhoa = param->rhoa;
            rhow = param->rhow;
            esat = 0.6108 * exp(17.27 * temp / (temp + 237.3));
            qsat = 0.622 * esat / (pres - 0.376 * esat);
            resi = 94.909 * pow(velo, -0.9036);
            for (ii = 0; ii < param->n3ci; ii++)
            {
                if (gmap->istop[ii] == 1)
                {
                    jj = gmap->top2d[ii];
                    if ((*data)->dept[jj] > 0.0 & param->sim_shallowwater == 1)
                    {(*data)->evap[jj] = (*data)->current_evap[0];}
                    else
                    {
                        alpha = 1.8 * ((*data)->wc[ii] - param->wcr) / ((*data)->wc[ii] - param->wcr + 0.3);
                        if (alpha > 1.0)    {alpha = 1.0;}
                        qsuf = alpha * qsat;
                        (*data)->evap[jj] = rhoa * (qsuf - humi) / (rhow * resi);
                        // for ex4 - Toy delta
                        // (*data)->evap[jj] = (*data)->evap[jj] * 50.0;
                        if ((*data)->evap[jj] < 0.0)    {(*data)->evap[jj] = 0.0;}
                    }
                }
            }
        }

        if (param->evap_file == 1 | param->evap_model == 1)
        {
            for (ii = 0; ii < param->n2ci; ii++)    {(*data)->qtop[ii] += (*data)->evap[ii] / (*data)->wcs[ii];}
        }
        if ((*data)->rain[0] != 0.0)
        {
            // A temporary limiter on max rainfall rate, ZhiLi20211209
            // if ((*data)->rain[0] > 1e-5)    {(*data)->rain[0] = 1e-5;}
            for (ii = 0; ii < param->n2ci; ii++)    {
                if (gmap->jj[ii] < param->ny-1) {
                    (*data)->qtop[ii] -= (*data)->rain[0];
                }
            }
        }
    }
    else
    {
        for (ii = 0; ii < param->n2ci; ii++)
        {(*data)->evap[ii] = (*data)->current_evap[0];}
    }

    // for Geng2015
    // for (ii = 0; ii < param->n3ci; ii++)
    // {
    //     if (gmap->istop[ii] == 1)
    //     {
    //         if (gmap->top2d[ii] < 20)    {(*data)->qtop[gmap->top2d[ii]] = 0.0;}
    //     }
    // }
}

// >>>>> Calculate total mass in the domain <<<<<
double get_mass(Data **data, Map *gmap, Config *param) {
    int ii;
    double mass = 0.0;
    // actually we only calculate volume
    if (param->sim_shallowwater == 1)   {
        for (ii = 0; ii < param->n2ci; ii++)    {
            mass += (*data)->Vs[ii];
        }
    }
    if (param->sim_groundwater == 1)    {
        for (ii = 0; ii < param->n3ci; ii++)    {
            mass += (*data)->wc[ii] * param->dx * param->dy * gmap->dz3d[ii];
        }
    }
    return mass;
}

// >>>>> Print information upon completion <<<<<
void print_end_info(Data **data, Map *smap, Map *gmap, Config *param, int irank)
{
    int ii, root=0;
    double vloss_tot;

    mpi_print("\n >>>>> Simulation completed! <<<<<",root);

    if (irank == root)  {printf(" >> Domain dimension = (%d, %d)\n",param->nx,param->ny);}

    if (param->sim_groundwater == 1)
    {
        if (irank == root)  {printf(" >> Total number of vertical layers = %d\n",param->nz);}
        // mass loss
        if (param->use_mpi == 1)    {mpi_gather_double((*data)->vloss_root, (*data)->vloss, param->n3ci, root);}
        else    {for (ii = 0; ii < param->N3CI; ii++)    {(*data)->vloss_root[ii] = (*data)->vloss[ii];}}
        if (irank == root)
        {
            for (ii = 0; ii < param->N3CI; ii++)    {vloss_tot += (*data)->vloss_root[ii];}
            printf(" >> Total volume loss = %f m^3\n",vloss_tot);
        }

    }

}
