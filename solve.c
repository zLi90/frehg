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
#include "utility.h"


void solve(Data **data, Map *smap, Map *gmap, Config *param, int irank, int nrank);
void get_tide(Data **data, Config *param);
void get_evaprain(Data **data, Config *param);
void print_end_info(Data **data, Map *smap, Map *gmap, Config *param, int irank);

void solve(Data **data, Map *smap, Map *gmap, Config *param, int irank, int nrank)
{
    int t_save, tday, ii, tt = 1;
    float t0, t1, tstep, last_save = 0.0, t_current = 0.0;
    double q1, q2, tq1, tq2, ms = 1.157e-7, spd = 86400.0;
    double qflux[36] = {0.0, 0.0, -1.0, -0.6, -0.8, -0.5, -0.1, -0.3, -0.6, -1.2, -1.2, 0.1, 0.05, 0.1,
            0.15, 0.05, 0.1, 0.1, -0.5, 0.05, 0.07, 0.05, 0.04, 0.08, 0.1, 0.05, 0.03, -0.6, -1.5, -1.5,
            -0.7, 0.2, 0.15, 0.05, 0.05, 0.0};
    // save initial condition
    write_output(data, gmap, param, 0, 0, irank);
    // begin time stepping
    mpi_print(" >>> Beginning Time loop !", irank);
    while (t_current < param->Tend)
    {
        tday = floor(t_current/spd);
        tq1 = tday*spd;
        tq2 = (tday+1)*spd;
        q1 = qflux[tday]*ms;
        q2 = qflux[tday+1]*ms;

        if (irank == 0) {t0 = clock();}
        (*data)->repeat[0] = 0;
        t_current += param->dt;
        // get boundary condition
        get_tide(data, param);
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

// >>>>> Update tidal elevation
void get_tide(Data **data, Config *param)
{
    int kk;
    if (param->tide_file == 0)
    {
        for (kk = 0; kk < param->n_tide; kk++)
        {(*data)->tide[kk] = param->init_tide[kk] + (*data)->offset[0];}
    }
    else
    {
        // interpolate tide for current time level
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
