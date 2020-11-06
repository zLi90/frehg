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
#include "utility.h"


void solve(Data **data, Map *smap, Map *gmap, Config *param, int irank, int nrank);
void get_current_bc(Data **data, Config *param, double t_current);
void get_evaprain(Data **data, Config *param);
void print_end_info(Data **data, Map *smap, Map *gmap, Config *param, int irank);

void solve(Data **data, Map *smap, Map *gmap, Config *param, int irank, int nrank)
{
    int t_save, tday, ii, kk, tt = 1;
    float t0, t1, tstep, last_save = 0.0, t_current = 0.0;
    // save initial condition
    write_output(data, gmap, param, 0, 0, irank);
    // begin time stepping
    mpi_print(" >>> Beginning Time loop !", irank);
    while (t_current < param->Tend)
    {
        if (irank == 0) {t0 = clock();}
        (*data)->repeat[0] = 0;
        t_current += param->dt;
        // get boundary condition
        // get_tide(data, param, t_current);
        get_current_bc(data, param, t_current);
        get_evaprain(data, param);
        // execute solvers
        if (param->sim_shallowwater == 1)
        {solve_shallowwater(data, smap, gmap, param, irank, nrank);}

        if (param->sim_groundwater == 1)
        {
            solve_groundwater(data, smap, gmap, param, irank, nrank);
            if ((*data)->repeat[0] == 1)
            {
                if (param->dt_adjust == 1 & param->dt > param->dt_min)
                {
                    printf("   >>> Repeat the previous step with dt from %f --> %f!\n",param->dt,param->dt*0.5);
                    param->dt = 0.5 * param->dt;
                    t_current -= 0.5 * param->dt;
                    solve_groundwater(data, smap, gmap, param, irank, nrank);
                }
                else
                {mpi_print("  >>> CFL limiter violated for groundwater solver! Should reduce dt!\n",irank);}
            }
        }

        if (param->sim_shallowwater == 1)
        {shallowwater_velocity(data, smap, gmap, param, irank, nrank);}
        // printf("tidal elevation = %f,   surf = %f\n",(*data)->current_tide[0],(*data)->eta[8686]);

        // scalar transport
        if (param->n_scalar > 0 & param->sim_shallowwater == 1)
        {for (kk = 0; kk < param->n_scalar; kk++)    {scalar_shallowwater(data, smap, param, irank, nrank, kk);}}
        if (param->n_scalar > 0 & param->sim_groundwater == 1)
        {for (kk = 0; kk < param->n_scalar; kk++)    {scalar_groundwater(data, gmap, param, irank, nrank, kk);}}

        // model output
        if (fabs(t_current - last_save - param->dt_out) < 0.5*param->dt)
        {
            t_save = round(t_current / param->dt_out) * param->dt_out;
            last_save = t_current;
            write_output(data, gmap, param, t_save, 0, irank);
        }
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
    // printf("  >> Qin = %f, Qout = %f\n",(*data)->qbc[0],(*data)->qbc[1]);
    print_end_info(data, smap, gmap, param, irank);
}

// >>>>> Get BC at the current time step
void get_current_bc(Data **data, Config *param, double t_current)
{
    int kk, ss, glob_ind;
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
                if (param->scalar_file[glob_ind] == 0)
                {(*data)->current_s_tide[ss][kk] = param->s_tide[glob_ind];}
                else
                {
                    (*data)->current_s_tide[ss][kk] =
                        interp_bc((*data)->t_s_tide[ss][kk],(*data)->s_tide[ss][kk],t_current,param->scalar_dat_len[glob_ind]);
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
                if (param->scalar_file[glob_ind] == 0)
                {(*data)->current_s_inflow[ss][kk] = param->s_inflow[glob_ind];}
                else
                {
                    (*data)->current_s_inflow[ss][kk] =
                        interp_bc((*data)->t_s_inflow[ss][kk],(*data)->s_inflow[ss][kk],t_current,param->scalar_dat_len[glob_ind]);
                }
            }
        }
    }
}

// >>>>> Update rainfall and evaporation flux
void get_evaprain(Data **data, Config *param)
{
    if (param->rain_file == 0)
    {(*data)->rain[0] = param->q_rain;}
    else
    {
        // interpolate rainfall data
    }
    if (param->evap_file == 0)
    {(*data)->evap[0] = param->q_evap;}
    else
    {
        // interpolate evaporation data
    }
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
