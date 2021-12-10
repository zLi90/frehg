// Initialize model runs
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<mpi.h>

#include"configuration.h"
#include"initialize.h"
#include"groundwater.h"
#include"mpifunctions.h"
#include"map.h"
#include"shallowwater.h"
#include"scalar.h"
#include"solve.h"
#include"utility.h"


void init(Data **data, Map **smap, Map **gmap, Config **param, int irank, int nrank);
void init_domain(Config **param);
void init_Data(Data **data, Config *param);
void read_bathymetry(Data **data, Config *param, int irank, int nrank);
void boundary_bath(Data **data, Map *smap, Config *param, int irank, int nrank);
void ic_surface(Data **data, Map *smap, Map *gmap, Config *param, int irank, int nrank);
void bc_surface(Data **data, Map *smap, Config *param, int irank);
void get_BC_location(int **loc, int *loc_len, Config *param, int irank, int n_bc, int *locX, int *locY);
void update_depth(Data **data, Map *smap, Config *param, int irank);
void ic_subsurface(Data **data, Map *gmap, Config *param, int irank, int nrank);
void restart_subsurface(double *ic_array, char *fname, Config *param, int irank);
void init_subgrid(Data **data, Map *smap, Config *param, int irank, int nrank);

// >>>>> Initialize FREHG <<<<<
void init(Data **data, Map **smap, Map **gmap, Config **param, int irank, int nrank)
{
    int ii;
    // domain partition
    init_domain(param);
    // generate bathymetry
    read_bathymetry(data, *param, irank, irank);
    // build maps
    build_surf_map(smap, *param);
    build_subsurf_map(gmap, *smap, (*data)->bottom, (*data)->offset, *param, irank);
    mpi_print(" >>> Connection maps built !", irank);
    // boundary bathymetry
    boundary_bath(data, *smap, *param, irank, nrank);
    // initialize data array
    init_Data(data, *param);
    mpi_print(" >>> Data array initialized !", irank);
    // boundary condition for shallow water solver
    bc_surface(data, *smap, *param, irank);
    mpi_print(" >>> Boundary condition initialized !", irank);
    get_current_bc(data, *param, 0.0);
    mpi_print(" >>> Boundary data read !", irank);
    // initial condition for shallow water solver
    ic_surface(data, *smap, *gmap, *param, irank, nrank);
    mpi_print(" >>> Initial conditions constructed for surface domain !", irank);
    update_drag_coef(data, *param);
    // subgrid model
    // if ((*param)->use_subgrid)
    // {init_subgrid(data, *smap, *param, irank, nrank);}
    // initial condition for groundwater solver
    if ((*param)->sim_groundwater == 1)
    {
        ic_subsurface(data, *gmap, *param, irank, nrank);
        mpi_print(" >>> Initial conditions constructed for subsurface domain !", irank);
        // boundary condition for groundwater solver
        if ((*param)->n_scalar > 0)
        {
            for (ii = 0; ii < (*param)->n_scalar; ii++)    {enforce_scalar_bc(data, *gmap, *param, ii, irank);}
            if ((*param)->baroclinic == 1)
            {update_rhovisc(data, *gmap, *param, irank);}
        }
        enforce_head_bc(data, *gmap, *param);
        mpi_print(" >>> Initial conditions applied !", irank);
    }

    mpi_print(" >>> Initialization completed !", irank);

}

// >>>>> Initialize domain partition and bathymetry
void init_domain(Config **param)
{
    int ii;
    // total number of grid cells
    (*param)->nx = (*param)->NX / (*param)->mpi_nx;
    (*param)->ny = (*param)->NY / (*param)->mpi_ny;
    (*param)->n2ci = (*param)->nx * (*param)->ny;
    (*param)->n2ct = ((*param)->nx + 2) * ((*param)->ny + 2);
    (*param)->N2CI = (*param)->NX * (*param)->NY;
}

// >>>>> Read bathymetry <<<<<
void read_bathymetry(Data **data, Config *param, int irank, int nrank)
{
    int ii, jj, col, row, xrank, yrank;
    char fullname[100];
    double z_min;
    *data = malloc(sizeof(Data));
    (*data)->bottom = malloc(param->n2ct*sizeof(double));
    (*data)->bottomXP = malloc(param->n2ct*sizeof(double));
    (*data)->bottomYP = malloc(param->n2ct*sizeof(double));
    (*data)->bottom_root = malloc(param->N2CI*sizeof(double));
    (*data)->offset = malloc(1*sizeof(double));
    for (ii = 0; ii < param->n2ct; ii++)    {(*data)->bottom[ii] = 0.0;}
    // assign internal bathymetry
    if (param->bath_file == 1)
    {
        strcpy(fullname, param->finput);
        if (param->use_subgrid == 0)    {strcat(fullname, "bath");}
        else    {strcat(fullname, "bath_sub");}
        load_data((*data)->bottom_root, fullname, param->N2CI);
        // calculate bathymetry offset
        z_min = getMin((*data)->bottom_root, param->N2CI);
        if (z_min >= 0) {(*data)->offset[0] = 0;}
        else    {(*data)->offset[0] = -z_min;}
        // bathymetry for each rank
        root_to_rank((*data)->bottom_root, (*data)->bottom, param, irank, nrank, 0, (*data)->offset[0]);

        // for (ii = 0; ii < param->n2ci; ii++)
        // {
        //     xrank = irank % param->mpi_nx;
        //     yrank = floor(irank/param->mpi_nx);
        //     // jj = yrank * param->mpi_nx + xrank;
        //     // (*data)->bottom[ii] = (*data)->bottom_root[jj*param->n2ci+ii] + (*data)->offset[0];
        //     col = floor(ii / param->nx);
        //     row = ii % param->nx;
        //     jj = yrank*param->mpi_nx*param->n2ci + col*param->NX + xrank*param->nx + row;
        //     (*data)->bottom[ii] = (*data)->bottom_root[jj] + (*data)->offset[0];
        // }
    }
    else
    {
        (*data)->offset[0] = 0.0;
        for (ii = 0; ii < param->n2ci; ii++)    {(*data)->bottom[ii] = 0.0;}
    }
}

// >>>>> Set boundary bathymetry <<<<<
void boundary_bath(Data **data, Map *smap, Config *param, int irank, int nrank)
{
    int ii;
    int jj, col, row, xrank, yrank;
    char fullname[100];
    for (ii = 0; ii < param->nx; ii++)
    {
        (*data)->bottom[smap->jMou[ii]] = (*data)->bottom[smap->jMin[ii]];
        (*data)->bottom[smap->jPou[ii]] = (*data)->bottom[smap->jPin[ii]];
    }
    for (ii = 0; ii < param->ny; ii++)
    {
        (*data)->bottom[smap->iMou[ii]] = (*data)->bottom[smap->iMin[ii]];
        (*data)->bottom[smap->iPou[ii]] = (*data)->bottom[smap->iPin[ii]];
    }
    if (param->use_mpi == 1)    {mpi_exchange_surf((*data)->bottom, smap, 2, param, irank, nrank);}
    // bathymetry at cell edges
    if (param->use_subgrid == 0)
    {
        for (ii = 0; ii < param->n2ci; ii++)
        {
            (*data)->bottomXP[ii] = (*data)->bottom[ii];
            if ((*data)->bottom[smap->iPjc[ii]] > (*data)->bottom[ii])
            {(*data)->bottomXP[ii] = (*data)->bottom[smap->iPjc[ii]];}
            (*data)->bottomYP[ii] = (*data)->bottom[ii];
            if ((*data)->bottom[smap->icjP[ii]] > (*data)->bottom[ii])
            {(*data)->bottomYP[ii] = (*data)->bottom[smap->icjP[ii]];}
        }
    }
    else
    {
        (*data)->edges_root = malloc(2*param->N2CI*sizeof(double));
        strcpy(fullname, param->finput);
        strcat(fullname, "bath_edges");
        load_data((*data)->edges_root, fullname, 2*param->N2CI);
        // edges for each rank
        root_to_rank((*data)->edges_root, (*data)->bottomXP, param, irank, nrank, 0, (*data)->offset[0]);
        root_to_rank((*data)->edges_root, (*data)->bottomYP, param, irank, nrank, param->N2CI, (*data)->offset[0]);
    }
    // edge bathymetry at boundaries
    for (ii = 0; ii < param->nx; ii++)
    {
        (*data)->bottomYP[smap->jPou[ii]] = (*data)->bottomYP[smap->jPin[ii]];
        (*data)->bottomYP[smap->jMou[ii]] = (*data)->bottom[smap->jMin[ii]];
        if ((*data)->bottom[smap->jMou[ii]] > (*data)->bottom[smap->jMin[ii]])
        {(*data)->bottomYP[smap->jMou[ii]] = (*data)->bottom[smap->jMou[ii]];}
    }
    for (ii = 0; ii < param->ny; ii++)
    {
        (*data)->bottomXP[smap->iPou[ii]] = (*data)->bottomXP[smap->iPin[ii]];
        (*data)->bottomXP[smap->iMou[ii]] = (*data)->bottom[smap->iMin[ii]];
        if ((*data)->bottom[smap->iMou[ii]] > (*data)->bottom[smap->iMin[ii]])
        {(*data)->bottomXP[smap->iMou[ii]] = (*data)->bottom[smap->iMou[ii]];}
    }
}

// >>>>> Initialize data array <<<<<
void init_Data(Data **data, Config *param)
{
    int ii, n2_root, n3_root;
    n2_root = param->n2ci*param->mpi_nx*param->mpi_ny;
    n3_root = param->n3ci*param->mpi_nx*param->mpi_ny;
    // surface fields
    (*data)->uu = malloc(param->n2ct*sizeof(double));
    (*data)->un = malloc(param->n2ct*sizeof(double));
    (*data)->uy = malloc(param->n2ct*sizeof(double));
    (*data)->vv = malloc(param->n2ct*sizeof(double));
    (*data)->vn = malloc(param->n2ct*sizeof(double));
    (*data)->vx = malloc(param->n2ct*sizeof(double));
    (*data)->eta = malloc(param->n2ct*sizeof(double));
    (*data)->etan = malloc(param->n2ct*sizeof(double));
    (*data)->dept = malloc(param->n2ct*sizeof(double));
    (*data)->deptx = malloc(param->n2ct*sizeof(double));
    (*data)->depty = malloc(param->n2ct*sizeof(double));
    (*data)->Fu = malloc(param->n2ct*sizeof(double));
    (*data)->Fv = malloc(param->n2ct*sizeof(double));
    (*data)->Ex = malloc(param->n2ct*sizeof(double));
    (*data)->Ey = malloc(param->n2ct*sizeof(double));
    (*data)->Dx = malloc(param->n2ct*sizeof(double));
    (*data)->Dy = malloc(param->n2ct*sizeof(double));
    (*data)->CDx = malloc(param->n2ct*sizeof(double));
    (*data)->CDy = malloc(param->n2ct*sizeof(double));
    (*data)->Vs = malloc(param->n2ct*sizeof(double));
    (*data)->Vsn = malloc(param->n2ct*sizeof(double));
    (*data)->Vflux = malloc(param->n2ci*sizeof(double));
    (*data)->Vsx = malloc(param->n2ct*sizeof(double));
    (*data)->Vsy = malloc(param->n2ct*sizeof(double));
    (*data)->Asz = malloc(param->n2ct*sizeof(double));
    (*data)->Aszx = malloc(param->n2ct*sizeof(double));
    (*data)->Aszy = malloc(param->n2ct*sizeof(double));
    (*data)->Asx = malloc(param->n2ct*sizeof(double));
    (*data)->Asy = malloc(param->n2ct*sizeof(double));
    (*data)->wtfx = malloc(param->n2ct*sizeof(double));
    (*data)->wtfy = malloc(param->n2ct*sizeof(double));

    // subgrid variables
    if (param->use_subgrid == 1)
    {
        (*data)->Vxp = malloc(param->n2ct*sizeof(double));
        (*data)->Vxm = malloc(param->n2ct*sizeof(double));
        (*data)->Vyp = malloc(param->n2ct*sizeof(double));
        (*data)->Vym = malloc(param->n2ct*sizeof(double));
        (*data)->Axp = malloc(param->n2ct*sizeof(double));
        (*data)->Axm = malloc(param->n2ct*sizeof(double));
        (*data)->Ayp = malloc(param->n2ct*sizeof(double));
        (*data)->Aym = malloc(param->n2ct*sizeof(double));

        // calculate pre-stored surface elevations
        param->nlay_sub = ceil((param->eta_sub_max - param->eta_sub_min) / param->deta_sub);
        (*data)->layers_sub = malloc(param->nlay_sub*sizeof(double));
        (*data)->layers_sub[0] = param->eta_sub_min + (*data)->offset[0];
        for (ii = 1; ii < param->nlay_sub; ii++)
        {(*data)->layers_sub[ii] = (*data)->layers_sub[ii-1] + param->deta_sub;}

        (*data)->Vxp_sub_root = malloc(n2_root*param->nlay_sub*sizeof(double));
        (*data)->Vxm_sub_root = malloc(n2_root*param->nlay_sub*sizeof(double));
        (*data)->Vyp_sub_root = malloc(n2_root*param->nlay_sub*sizeof(double));
        (*data)->Vym_sub_root = malloc(n2_root*param->nlay_sub*sizeof(double));
        (*data)->Axp_sub_root = malloc(n2_root*param->nlay_sub*sizeof(double));
        (*data)->Axm_sub_root = malloc(n2_root*param->nlay_sub*sizeof(double));
        (*data)->Ayp_sub_root = malloc(n2_root*param->nlay_sub*sizeof(double));
        (*data)->Aym_sub_root = malloc(n2_root*param->nlay_sub*sizeof(double));
        (*data)->Asz_sub_root = malloc(n2_root*param->nlay_sub*sizeof(double));
        (*data)->Vxp_sub = malloc(param->nlay_sub*sizeof(double));
        (*data)->Vyp_sub = malloc(param->nlay_sub*sizeof(double));
        (*data)->Vxm_sub = malloc(param->nlay_sub*sizeof(double));
        (*data)->Vym_sub = malloc(param->nlay_sub*sizeof(double));
        (*data)->Axp_sub = malloc(param->nlay_sub*sizeof(double));
        (*data)->Ayp_sub = malloc(param->nlay_sub*sizeof(double));
        (*data)->Axm_sub = malloc(param->nlay_sub*sizeof(double));
        (*data)->Aym_sub = malloc(param->nlay_sub*sizeof(double));
        (*data)->Asz_sub = malloc(param->nlay_sub*sizeof(double));
        for (ii = 0; ii < param->nlay_sub; ii++)
        {
            (*data)->Vxp_sub[ii] = malloc(param->n2ci*sizeof(double));
            (*data)->Vxm_sub[ii] = malloc(param->n2ci*sizeof(double));
            (*data)->Vyp_sub[ii] = malloc(param->n2ci*sizeof(double));
            (*data)->Vym_sub[ii] = malloc(param->n2ci*sizeof(double));
            (*data)->Axp_sub[ii] = malloc(param->n2ci*sizeof(double));
            (*data)->Axm_sub[ii] = malloc(param->n2ci*sizeof(double));
            (*data)->Ayp_sub[ii] = malloc(param->n2ci*sizeof(double));
            (*data)->Aym_sub[ii] = malloc(param->n2ci*sizeof(double));
            (*data)->Asz_sub[ii] = malloc(param->n2ci*sizeof(double));
        }
        // index of layer for each surface cell
        (*data)->eta_ind = malloc(param->n2ci*sizeof(double));
    }

    (*data)->uu_root = malloc(n2_root*sizeof(double));
    (*data)->vv_root = malloc(n2_root*sizeof(double));
    (*data)->un_root = malloc(n2_root*sizeof(double));
    (*data)->vn_root = malloc(n2_root*sizeof(double));
    (*data)->eta_root = malloc(n2_root*sizeof(double));
    (*data)->dept_root = malloc(n2_root*sizeof(double));
    (*data)->seep_root = malloc(n2_root*sizeof(double));
    (*data)->uu_out = malloc(n2_root*sizeof(double));
    (*data)->vv_out = malloc(n2_root*sizeof(double));
    (*data)->un_out = malloc(n2_root*sizeof(double));
    (*data)->vn_out = malloc(n2_root*sizeof(double));
    (*data)->eta_out = malloc(n2_root*sizeof(double));
    (*data)->dept_out = malloc(n2_root*sizeof(double));
    (*data)->seep_out = malloc(n2_root*sizeof(double));

    (*data)->reset_seepage = malloc(param->n2ci*sizeof(int));
    (*data)->qseepage = malloc(param->n2ci*sizeof(double));
    (*data)->cflx = malloc(param->n2ci*sizeof(double));
    (*data)->cfly = malloc(param->n2ci*sizeof(double));
    (*data)->cfl_active = malloc(param->n2ci*sizeof(double));

    (*data)->rain = malloc(1*sizeof(double));
    (*data)->rain_sum = malloc(1*sizeof(double));
    (*data)->rain_sum[0] = 0.0;
    (*data)->evap = malloc(param->n2ci*sizeof(double));

    // subsurface fields
    (*data)->repeat = malloc(1*sizeof(int));
    if (param->sim_groundwater == 1)
    {
        (*data)->h = malloc(param->n3ct*sizeof(double));
        (*data)->hp = malloc(param->n3ct*sizeof(double));
        (*data)->hn = malloc(param->n3ct*sizeof(double));
        (*data)->hwc = malloc(param->n3ct*sizeof(double));
        (*data)->wc = malloc(param->n3ct*sizeof(double));
        (*data)->wcn = malloc(param->n3ct*sizeof(double));
        (*data)->wcp = malloc(param->n3ct*sizeof(double));
        (*data)->wch = malloc(param->n3ct*sizeof(double));
        (*data)->ch = malloc(param->n3ct*sizeof(double));
        (*data)->wcs = malloc(param->n3ct*sizeof(double));
        (*data)->wcr = malloc(param->n3ct*sizeof(double));
        (*data)->vga = malloc(param->n3ct*sizeof(double));
        (*data)->vgn = malloc(param->n3ct*sizeof(double));
        (*data)->Ksx = malloc(param->n3ct*sizeof(double));
        (*data)->Ksy = malloc(param->n3ct*sizeof(double));
        (*data)->Ksz = malloc(param->n3ct*sizeof(double));
        (*data)->Kx = malloc(param->n3ct*sizeof(double));
        (*data)->Ky = malloc(param->n3ct*sizeof(double));
        (*data)->Kz = malloc(param->n3ct*sizeof(double));
        (*data)->r_rho = malloc(param->n3ct*sizeof(double));
        (*data)->r_rhoxp = malloc(param->n3ct*sizeof(double));
        (*data)->r_rhoyp = malloc(param->n3ct*sizeof(double));
        (*data)->r_rhozp = malloc(param->n3ct*sizeof(double));
        (*data)->r_rhon = malloc(param->n3ct*sizeof(double));
        (*data)->r_visc = malloc(param->n3ct*sizeof(double));
        (*data)->r_viscxp = malloc(param->n3ct*sizeof(double));
        (*data)->r_viscyp = malloc(param->n3ct*sizeof(double));
        (*data)->r_visczp = malloc(param->n3ct*sizeof(double));
        (*data)->qx = malloc(param->n3ct*sizeof(double));
        (*data)->qy = malloc(param->n3ct*sizeof(double));
        (*data)->qz = malloc(param->n3ct*sizeof(double));
        (*data)->Vg = malloc(param->n3ct*sizeof(double));
        (*data)->Vgn = malloc(param->n3ct*sizeof(double));
        (*data)->Vgflux = malloc(param->n3ci*sizeof(double));
        (*data)->room = malloc(param->n3ct*sizeof(double));
        (*data)->vloss = malloc(param->n3ci*sizeof(double));
        (*data)->h_root = malloc(n3_root*sizeof(double));
        (*data)->wc_root = malloc(n3_root*sizeof(double));
        (*data)->vloss_root = malloc(n3_root*sizeof(double));
        (*data)->qx_root = malloc(n3_root*sizeof(double));
        (*data)->qy_root = malloc(n3_root*sizeof(double));
        (*data)->qz_root = malloc(n3_root*sizeof(double));
        (*data)->h_out = malloc(n3_root*sizeof(double));
        (*data)->wc_out = malloc(n3_root*sizeof(double));
        (*data)->qx_out = malloc(n3_root*sizeof(double));
        (*data)->qy_out = malloc(n3_root*sizeof(double));
        (*data)->qz_out = malloc(n3_root*sizeof(double));
        (*data)->dh6 = malloc(6*sizeof(double));
        (*data)->rsplit = malloc(6*sizeof(double));
        (*data)->qbc = malloc(2*sizeof(double));
        (*data)->qbc[0] = 0.0;
        (*data)->qbc[1] = 0.0;
        (*data)->qtop = malloc(param->n2ci*sizeof(double));
    }


    // scalar
    if (param->n_scalar > 0)
    {
        (*data)->s_surf = malloc(param->n_scalar*sizeof(double *));
        (*data)->s_subs = malloc(param->n_scalar*sizeof(double *));
        (*data)->sm_surf = malloc(param->n_scalar*sizeof(double *));
        (*data)->sm_subs = malloc(param->n_scalar*sizeof(double *));
        (*data)->s_surf_root = malloc(param->n_scalar*sizeof(double *));
        (*data)->s_subs_root = malloc(param->n_scalar*sizeof(double *));
        (*data)->s_surf_out = malloc(param->n_scalar*sizeof(double *));
        (*data)->s_subs_out = malloc(param->n_scalar*sizeof(double *));
        (*data)->sseepage = malloc(param->n_scalar*sizeof(double *));
        (*data)->s_surfkP = malloc(param->n_scalar*sizeof(double *));
        for (ii = 0; ii < param->n_scalar; ii++)
        {
            (*data)->s_surf[ii] = malloc(param->n2ct*sizeof(double));
            (*data)->s_subs[ii] = malloc(param->n3ct*sizeof(double));
            (*data)->sm_surf[ii] = malloc(param->n2ct*sizeof(double));
            (*data)->sm_subs[ii] = malloc(param->n3ct*sizeof(double));
            (*data)->s_surf_root[ii] = malloc(n2_root*sizeof(double));
            (*data)->s_subs_root[ii] = malloc(n3_root*sizeof(double));
            (*data)->s_surf_out[ii] = malloc(n2_root*sizeof(double));
            (*data)->s_subs_out[ii] = malloc(n3_root*sizeof(double));
            (*data)->sseepage[ii] = malloc(param->n2ci*sizeof(double));
            (*data)->s_surfkP[ii] = malloc(param->n2ci*sizeof(double));
        }
        (*data)->Dxx = malloc(param->n3ct*sizeof(double));
        (*data)->Dxy = malloc(param->n3ct*sizeof(double));
        (*data)->Dxz = malloc(param->n3ct*sizeof(double));
        (*data)->Dyy = malloc(param->n3ct*sizeof(double));
        (*data)->Dyx = malloc(param->n3ct*sizeof(double));
        (*data)->Dyz = malloc(param->n3ct*sizeof(double));
        (*data)->Dzz = malloc(param->n3ct*sizeof(double));
        (*data)->Dzx = malloc(param->n3ct*sizeof(double));
        (*data)->Dzy = malloc(param->n3ct*sizeof(double));
    }

    // linear system solver
    (*data)->Sct = malloc(param->n2ci*sizeof(double));
    (*data)->Sxp = malloc(param->n2ci*sizeof(double));
    (*data)->Sxm = malloc(param->n2ci*sizeof(double));
    (*data)->Syp = malloc(param->n2ci*sizeof(double));
    (*data)->Sym = malloc(param->n2ci*sizeof(double));
    (*data)->Srhs = malloc(param->n2ci*sizeof(double));

    (*data)->Gct = malloc(param->n3ci*sizeof(double));
    (*data)->Gxp = malloc(param->n3ci*sizeof(double));
    (*data)->Gxm = malloc(param->n3ci*sizeof(double));
    (*data)->Gyp = malloc(param->n3ci*sizeof(double));
    (*data)->Gym = malloc(param->n3ci*sizeof(double));
    (*data)->Gzp = malloc(param->n3ci*sizeof(double));
    (*data)->Gzm = malloc(param->n3ci*sizeof(double));
    (*data)->Grhs = malloc(param->n3ci*sizeof(double));
}

// >>>>> Initial condition for shallow water solver <<<<<
void ic_surface(Data **data, Map *smap, Map *gmap, Config *param, int irank, int nrank)
{
    int ii, jj, kk, col, row, xrank, yrank;
    char fid[2];
    // get monitoring locations
    if (param->n_monitor > 0)
    {
        (*data)->monitor_rank = malloc(param->n_monitor * sizeof(int));
        (*data)->monitor = malloc(param->n_monitor * sizeof(int));
        for (ii = 0; ii < param->n_monitor; ii++)
        {
            row = floor(param->monitor_locX[ii] / param->nx);
            col = floor(param->monitor_locY[ii] / param->ny);
            (*data)->monitor_rank[ii] = col * param->mpi_nx + row;
            (*data)->monitor[ii] = (param->monitor_locY[ii] - col*param->ny) * param->nx + (param->monitor_locX[ii] - row*param->nx);
        }
    }
    // get rainfall / evaporation rate
    get_evaprain(data, gmap, param, 0.0);
    // initial surface elevation
    if (param->eta_file == 0)
    {
        for (ii = 0; ii < param->N2CI; ii++)
        {
            (*data)->eta_root[ii] = param->init_eta + (*data)->offset[0];
            // if (smap->jj[ii] < 85)
            // {
            //     (*data)->eta_root[ii] = 0.4 + (*data)->offset[0];
            // }
        }
    }
    else
    {
        char fullname[50];
        strcpy(fullname, param->finput);
        strcat(fullname, "surf_ic");
        load_data((*data)->eta_root, fullname, param->N2CI);
        for (ii = 0; ii < param->N2CI; ii++)
        {(*data)->eta_root[ii] = (*data)->eta_root[ii] + (*data)->offset[0];}
    }
    // initial velocity
    if (param->uv_file == 0)
    {
        for (ii = 0; ii < param->N2CI; ii++)
        {(*data)->uu_root[ii] = 0.0;    (*data)->vv_root[ii] = 0.0;}
    }
    else
    {
        char fullname1[50];
        strcpy(fullname1, param->finput);
        strcat(fullname1, "uu_ic");
        load_data((*data)->uu_root, fullname1, param->N2CI);
        char fullname2[50];
        strcpy(fullname2, param->finput);
        strcat(fullname2, "vv_ic");
        load_data((*data)->vv_root, fullname2, param->N2CI);
    }
    // initial scalar
    if (param->n_scalar > 0)
    {
        for (kk = 0; kk < param->n_scalar; kk++)
        {
            if (param->scalar_surf_file[kk] == 0)
            {
                for (ii = 0; ii < param->N2CI; ii++)
                {(*data)->s_surf_root[kk][ii] = param->init_s_surf[kk];}
            }
            else
            {
                char fullname[50];
                strcpy(fullname, param->finput);
                strcat(fullname, "scalar_surf_ic");
                sprintf(fid, "%d", kk+1);
                strcat(fullname, fid);
                load_data((*data)->s_surf_root[kk], fullname, param->N2CI);
                for (ii = 0; ii < param->N2CI; ii++)
                {if (isnan((*data)->s_surf_root[kk][ii]))    {(*data)->s_surf_root[kk][ii] = 0.0;}}
            }
        }
    }
    // assign value to irank
    for (ii = 0; ii < param->n2ci; ii++)
    {
        xrank = irank % param->mpi_nx;
        yrank = floor(irank/param->mpi_nx);
        // jj = yrank * param->mpi_nx + xrank;
        col = floor(ii / param->nx);
        row = ii % param->nx;
        jj = yrank*param->mpi_nx*param->n2ci + col*param->NX + xrank*param->nx + row;

        (*data)->eta[ii] = (*data)->eta_root[jj];
        (*data)->uu[ii] = (*data)->uu_root[jj];
        (*data)->vv[ii] = (*data)->vv_root[jj];
        if ((*data)->eta[ii] < (*data)->bottom[ii])
        {(*data)->eta[ii] = (*data)->bottom[ii];}
        if (param->n_scalar > 0)
        {
            for (kk = 0; kk < param->n_scalar; kk++)
            {
                (*data)->s_surf[kk][ii] = (*data)->s_surf_root[kk][jj];
                // if ((*data)->eta[ii] <= (*data)->bottom[ii]+(*data)->offset[0])   {(*data)->s_surf[kk][ii] = 0.0;}
                if ((*data)->eta[ii] <= (*data)->bottom[ii])   {(*data)->s_surf[kk][ii] = 0.0;}
            }
        }
    }

    // enforce boundary conditions
    if (param->use_mpi == 1 & param->sim_shallowwater == 1)
    {
        mpi_exchange_surf((*data)->eta, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->uu, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->vv, smap, 2, param, irank, nrank);
        if (param->n_scalar > 0)
        {
            for (kk = 0; kk < param->n_scalar; kk++)
            {mpi_exchange_surf((*data)->s_surf[kk], smap, 2, param, irank, nrank);}
        }
    }
    enforce_surf_bc(data, smap, param, irank, nrank);
    // calculate depth
    update_depth(data, smap, param, irank);
    if (param->use_mpi == 1 & param->sim_shallowwater == 1)
    {
        mpi_exchange_surf((*data)->dept, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->deptx, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->depty, smap, 2, param, irank, nrank);
    }
    // calculate subgrid areas
    if (param->use_subgrid == 1)
    {
        init_subgrid(data, smap, param, irank, nrank);
    }
    else
    {
        for (ii = 0; ii < param->n2ct; ii++)
        {
            (*data)->Vs[ii] = (*data)->dept[ii] * param->dx * param->dy;
            (*data)->Vsn[ii] = (*data)->dept[ii] * param->dx * param->dy;
            if ((*data)->dept[ii] > 0) {(*data)->Asz[ii] = param->dx * param->dy;}
            else    {(*data)->Asz[ii] = 0.0;}
            (*data)->Asx[ii] = (*data)->deptx[ii] * param->dy;
            (*data)->Asy[ii] = (*data)->depty[ii] * param->dx;
        }
        for (ii = 0; ii < param->n2ci; ii++)
        {
            (*data)->Vsx[ii] = 0.5 * ((*data)->Vs[ii] + (*data)->Vs[smap->iPjc[ii]]);
            (*data)->Vsy[ii] = 0.5 * ((*data)->Vs[ii] + (*data)->Vs[smap->icjP[ii]]);
            (*data)->Aszx[ii] = 0.5 * ((*data)->Asz[ii] + (*data)->Asz[smap->iPjc[ii]]);
            (*data)->Aszy[ii] = 0.5 * ((*data)->Asz[ii] + (*data)->Asz[smap->icjP[ii]]);
        }
        for (ii = 0; ii < param->nx; ii++)
        {
            (*data)->Vsy[smap->jMou[ii]] = (*data)->Vs[smap->jMin[ii]];
            (*data)->Aszy[smap->jMou[ii]] = (*data)->Asz[smap->jMin[ii]];
        }
        for (ii = 0; ii < param->ny; ii++)
        {
            (*data)->Vsx[smap->iMou[ii]] = (*data)->Vs[smap->iMin[ii]];
            (*data)->Aszx[smap->iMou[ii]] = (*data)->Asz[smap->iMin[ii]];
        }
    }

    if (param->use_mpi == 1 & param->sim_shallowwater == 1)
    {
        mpi_exchange_surf((*data)->Vs, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Vsx, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Vsy, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Asx, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Asy, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Asz, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Aszx, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->Aszy, smap, 2, param, irank, nrank);
    }
    // scalar mass
    if (param->n_scalar > 0)
    {
        for (kk = 0; kk < param->n_scalar; kk++)
        {
            for (ii = 0; ii < param->n2ct; ii++)
            {(*data)->sm_surf[kk][ii] = (*data)->s_surf[kk][ii] * (*data)->Vs[ii];}
            for (ii = 0; ii < param->n2ci; ii++)
            {(*data)->sseepage[kk][ii] = 0.0;}
        }
    }
    // initialize other relevant variables
    for (ii = 0; ii < param->n2ci; ii++)
    {
        (*data)->qseepage[ii] = 0.0;
        (*data)->reset_seepage[ii] = 1;
    }
    for (ii = 0; ii < param->n2ct; ii++)
    {
        (*data)->CDx[ii] = 0.0;  (*data)->CDy[ii] = 0.0;
        (*data)->Dx[ii] = 0.0;  (*data)->Dy[ii] = 0.0;
        (*data)->Ex[ii] = 0.0;  (*data)->Ey[ii] = 0.0;
        (*data)->Fv[ii] = 0.0;  (*data)->Fv[ii] = 0.0;
        (*data)->uy[ii] = 0.0;  (*data)->vx[ii] = 0.0;
        (*data)->un[ii] = (*data)->uu[ii];
        (*data)->vn[ii] = (*data)->vv[ii];
        (*data)->wtfx[ii] = 0;  (*data)->wtfy[ii] = 0;
    }
}

// >>>>> Boundary condition for shallowwater solver <<<<<
void bc_surface(Data **data, Map *smap, Config *param, int irank)
{
    int kk, ii, ss, glob_ind, locX1, locX2, locY1, locY2, rank0, rank1, locx1, rankx, ranky, ncell, locy1, locy2;
    char fid[2];
    // get location of the tidal cells
    (*data)->tideloc = malloc(param->n_tide*sizeof(int *));
    (*data)->tideloc_len = malloc(param->n_tide*sizeof(int));
    get_BC_location((*data)->tideloc, (*data)->tideloc_len, param, irank, param->n_tide, param->tide_locX, param->tide_locY);
    // load tidal boundary condition
    (*data)->tide = malloc(param->n_tide*sizeof(double *));
    (*data)->t_tide = malloc(param->n_tide*sizeof(double *));
    (*data)->current_tide = malloc(param->n_tide*sizeof(double));
    for (kk = 0; kk < param->n_tide; kk++)
    {
        if (param->tide_file[kk] == 1)
        {
            char fullname[50];
            strcpy(fullname, param->finput);
            strcat(fullname, "tide");
            sprintf(fid, "%d", kk+1);
            strcat(fullname, fid);
            if (exist(fullname))
            {
                (*data)->tide[kk] = malloc(param->tide_dat_len[kk]*sizeof(double));
                (*data)->t_tide[kk] = malloc(param->tide_dat_len[kk]*sizeof(double));
                load_bc((*data)->tide[kk], (*data)->t_tide[kk], fullname, param->tide_dat_len[kk]);
            }
        }
        else
        {
            (*data)->tide[kk] = malloc(1*sizeof(double));
            (*data)->t_tide[kk] = malloc(1*sizeof(double));
            (*data)->tide[kk][0] = param->init_tide[kk];
            (*data)->t_tide[kk][0] = 0.0;
        }
    }
    // get inflow locations
    (*data)->inflowloc = malloc(param->n_inflow*sizeof(int *));
    (*data)->inflowloc_len = malloc(param->n_inflow*sizeof(int));
    get_BC_location((*data)->inflowloc, (*data)->inflowloc_len, param, irank, param->n_inflow, param->inflow_locX, param->inflow_locY);
    // load inflow boundary condition
    (*data)->inflow = malloc(param->n_inflow*sizeof(double *));
    (*data)->t_inflow = malloc(param->n_inflow*sizeof(double *));
    (*data)->current_inflow = malloc(param->n_inflow*sizeof(double));
    for (kk = 0; kk < param->n_inflow; kk++)
    {
        if (param->inflow_file[kk] == 1)
        {
            char fullname[150];
            strcpy(fullname, param->finput);
            strcat(fullname, "inflow");
            sprintf(fid, "%d", kk+1);
            strcat(fullname, fid);
            if (exist(fullname))
            {
                (*data)->inflow[kk] = malloc(param->inflow_dat_len[kk]*sizeof(double));
                (*data)->t_inflow[kk] = malloc(param->inflow_dat_len[kk]*sizeof(double));
                load_bc((*data)->inflow[kk], (*data)->t_inflow[kk], fullname, param->inflow_dat_len[kk]);
            }
        }
        else
        {
            (*data)->inflow[kk] = malloc(1*sizeof(double));
            (*data)->t_inflow[kk] = malloc(1*sizeof(double));
            (*data)->inflow[kk][0] = param->init_inflow[kk];
            (*data)->t_inflow[kk][0] = 0.0;
        }
    }

    // load wind boundary condition (As of 20210502, wind must be applied uniformly to the entire domain)
    (*data)->current_windspd = malloc(1*sizeof(double));
    (*data)->current_winddir = malloc(1*sizeof(double));
    if (param->sim_wind == 1)
    {
        if (param->wind_file == 1)
        {
            char fullname[150];
            strcpy(fullname, param->finput);
            strcat(fullname, "windspd");
            if (exist(fullname))
            {
                (*data)->wind_spd = malloc(param->wind_dat_len*sizeof(double));
                (*data)->t_wind = malloc(param->wind_dat_len*sizeof(double));
                load_bc((*data)->wind_spd, (*data)->t_wind, fullname, param->wind_dat_len);
            }
            char dirname[150];
            strcpy(dirname, param->finput);
            strcat(dirname, "winddir");
            if (exist(dirname))
            {
                (*data)->wind_dir = malloc(param->wind_dat_len*sizeof(double));
                load_bc((*data)->wind_dir, (*data)->t_wind, dirname, param->wind_dat_len);
            }
        }
        else
        {
            (*data)->wind_spd = malloc(1*sizeof(double));
            (*data)->wind_dir = malloc(1*sizeof(double));
            (*data)->wind_dir[0] = param->init_winddir;
            (*data)->wind_spd[0] = param->init_windspd;
        }
    }

    // load evaporation
    (*data)->current_evap = malloc(1*sizeof(double));
    if (param->evap_file == 1)
    {
        char fullname[150];
        strcpy(fullname, param->finput);
        strcat(fullname, "evap");
        if (exist(fullname))
        {
            (*data)->evap_data = malloc(param->evap_dat_len*sizeof(double));
            (*data)->t_evap = malloc(param->evap_dat_len*sizeof(double));
            load_bc((*data)->evap_data, (*data)->t_evap, fullname, param->evap_dat_len);
        }
    }

    // load rainfall
    (*data)->current_rain = malloc(1*sizeof(double));
    if (param->rain_file == 1)
    {
        char fullname[150];
        strcpy(fullname, param->finput);
        strcat(fullname, "rain");
        if (exist(fullname))
        {
            (*data)->rain_data = malloc(param->rain_dat_len*sizeof(double));
            (*data)->t_rain = malloc(param->rain_dat_len*sizeof(double));
            load_bc((*data)->rain_data, (*data)->t_rain, fullname, param->rain_dat_len);
        }
    }

    // load scalar boundary condition for tide
    (*data)->s_tide = malloc(param->n_scalar*sizeof(double **));
    (*data)->t_s_tide = malloc(param->n_scalar*sizeof(double **));
    (*data)->current_s_tide = malloc(param->n_scalar*sizeof(double *));
    for (ss = 0; ss < param->n_scalar; ss++)
    {
        (*data)->s_tide[ss] = malloc(param->n_tide*sizeof(double *));
        (*data)->t_s_tide[ss] = malloc(param->n_tide*sizeof(double *));
        (*data)->current_s_tide[ss] = malloc(param->n_tide*sizeof(double));
        for (kk = 0; kk < param->n_tide; kk++)
        {
            glob_ind = ss * param->n_tide + kk;
            if (param->scalar_tide_file[glob_ind] == 1)
            {
                char fullname[150];
                strcpy(fullname, param->finput);
                strcat(fullname, "scalar");
                sprintf(fid, "%d", ss+1);
                strcat(fullname, fid);
                strcat(fullname, "_tide");
                sprintf(fid, "%d", kk+1);
                strcat(fullname, fid);
                if (exist(fullname))
                {
                    (*data)->s_tide[ss][kk] = malloc(param->scalar_tide_datlen[glob_ind]*sizeof(double));
                    (*data)->t_s_tide[ss][kk] = malloc(param->scalar_tide_datlen[glob_ind]*sizeof(double));
                    load_bc((*data)->s_tide[ss][kk], (*data)->t_s_tide[ss][kk], fullname, param->scalar_tide_datlen[glob_ind]);
                }
            }
            else
            {
                (*data)->s_tide[ss][kk] = malloc(1*sizeof(double));
                (*data)->t_s_tide[ss][kk] = malloc(1*sizeof(double));
                (*data)->s_tide[ss][kk][0] = param->s_tide[glob_ind];
                (*data)->t_s_tide[ss][kk][0] = 0.0;
            }
        }
    }
    //scalar BC for inflow (now only support constant scalar for inflow )
    (*data)->s_inflow = malloc(param->n_scalar*sizeof(double **));
    (*data)->t_s_inflow = malloc(param->n_scalar*sizeof(double **));
    (*data)->current_s_inflow = malloc(param->n_scalar*sizeof(double *));
    for (ss = 0; ss < param->n_scalar; ss++)
    {
        (*data)->s_inflow[ss] = malloc(param->n_inflow*sizeof(double *));
        (*data)->t_s_inflow[ss] = malloc(param->n_inflow*sizeof(double *));
        (*data)->current_s_inflow[ss] = malloc(param->n_inflow*sizeof(double));
        for (kk = 0; kk < param->n_inflow; kk++)
        {
            glob_ind = ss * param->n_inflow + kk;
            if (param->scalar_inflow_file[glob_ind] == 1)
            {
                mpi_print("WARNING: For now, inflow scalar concentration must be a constant!", irank);
            }
            else
            {
                (*data)->s_inflow[ss][kk] = malloc(1*sizeof(double));
                (*data)->t_s_inflow[ss][kk] = malloc(1*sizeof(double));
                (*data)->s_inflow[ss][kk][0] = param->s_inflow[glob_ind];
                (*data)->t_s_inflow[ss][kk][0] = 0.0;
            }
        }
    }
}

// >>>>> Get the cell index for applying the boundary condition
void get_BC_location(int **loc, int *loc_len, Config *param, int irank, int n_bc, int *locX, int *locY)
{
    int ii, jj, kk, ll, ind, n_glob, xrank, yrank, col, row;
    int locX1, locX2, locY1, locY2;
    int *loc_glob;
    // loc = malloc(n_bc*sizeof(int *));
    // loc_len = malloc(n_bc*sizeof(int));
    for (kk = 0; kk < n_bc; kk++)
    {
        loc_len[kk] = 0;
        // get the global (i,j) index of the BC region
        locX1 = locX[2*kk];     locX2 = locX[2*kk+1];
        locY1 = locY[2*kk];     locY2 = locY[2*kk+1];
        // get the global 1D index
        n_glob = (locX2 - locX1 + 1) * (locY2 - locY1 + 1);
        loc_glob = malloc(n_glob*sizeof(int));
        ll = 0;
        for (jj = locY1; jj <= locY2; jj++)
        {
            for (ii = locX1; ii <= locX2; ii++)
            {loc_glob[ll] = jj * param->NX + ii;   ll+=1;}
        }
        // assign global index to local ranks
        xrank = irank % param->mpi_nx;
        yrank = floor(irank/param->mpi_nx);
        for (ii = 0; ii < param->n2ci; ii++)
        {
            col = floor(ii / param->nx);
            row = ii % param->nx;
            // get global index of the local cell
            ll = yrank*param->mpi_nx*param->n2ci + col*param->NX + xrank*param->nx + row;
            // check if the local cell is in the global BC region
            // count the number of BC cells in the local rank
            for (jj = 0; jj < n_glob; jj++)
            {if (ll == loc_glob[jj]) {loc_len[kk] += 1;}}
        }
        // assign the BC index
        if (loc_len[kk] == 0)
        {
            loc[kk] = malloc(1*sizeof(int));
            loc[kk][0] = -1;
        }
        else
        {
            loc[kk] = malloc(loc_len[kk]*sizeof(int));
            ind = 0;
            for (ii = 0; ii < param->n2ci; ii++)
            {
                col = floor(ii / param->nx);
                row = ii % param->nx;
                ll = yrank*param->mpi_nx*param->n2ci + col*param->NX + xrank*param->nx + row;
                for (jj = 0; jj < n_glob; jj++)
                {if (ll == loc_glob[jj]) {loc[kk][ind] = ii; ind += 1;}}
            }
        }
        free(loc_glob);
    }
}


// >>>>> Calculate surface depth from bathymetry and surface elevation <<<<<
void update_depth(Data **data, Map *smap, Config *param, int irank)
{
    int ii;
    double eta_hi, bot_hi, diff;
    // remove small depth
    for (ii = 0; ii < param->n2ci; ii++)
    {
        diff = (*data)->eta[ii] - (*data)->bottom[ii];
        if (diff > 0 & diff < param->min_dept)
        {(*data)->eta[ii] = (*data)->bottom[ii];}
    }
    // center depth
    for (ii = 0; ii < param->n2ct; ii++)
    {
        (*data)->dept[ii] = (*data)->eta[ii] - (*data)->bottom[ii];
        if ((*data)->dept[ii] <= param->min_dept)   {(*data)->dept[ii] = 0.0;}
        (*data)->deptx[ii] = 0.0;
        (*data)->depty[ii] = 0.0;
    }
    // internal face depth
    for (ii = 0; ii < param->n2ci; ii++)
    {
        // x
        eta_hi = (*data)->eta[ii];
        bot_hi = (*data)->bottom[ii];
        if ((*data)->eta[smap->iPjc[ii]] > eta_hi)  {eta_hi = (*data)->eta[smap->iPjc[ii]];}
        if ((*data)->bottom[smap->iPjc[ii]] > bot_hi)   {bot_hi = (*data)->bottom[smap->iPjc[ii]];}
        (*data)->deptx[ii] = eta_hi - bot_hi;

        // y
        eta_hi = (*data)->eta[ii];
        bot_hi = (*data)->bottom[ii];
        if ((*data)->eta[smap->icjP[ii]] > eta_hi)  {eta_hi = (*data)->eta[smap->icjP[ii]];}
        if ((*data)->bottom[smap->icjP[ii]] > bot_hi)   {bot_hi = (*data)->bottom[smap->icjP[ii]];}
        (*data)->depty[ii] = eta_hi - bot_hi;
    }
    // boundary face depth
    for (ii = 0; ii < param->nx; ii++)
    {
        // ym
        eta_hi = (*data)->eta[smap->jMin[ii]];
        bot_hi = (*data)->bottom[smap->jMin[ii]];
        if ((*data)->eta[smap->jMou[ii]] > eta_hi)  {eta_hi = (*data)->eta[smap->jMou[ii]];}
        if ((*data)->bottom[smap->jMou[ii]] > bot_hi)   {bot_hi = (*data)->bottom[smap->jMou[ii]];}
        (*data)->depty[smap->jMou[ii]] = eta_hi - bot_hi;
        // yp
        eta_hi = (*data)->eta[smap->jPin[ii]];
        bot_hi = (*data)->bottom[smap->jPin[ii]];
        if ((*data)->eta[smap->jPou[ii]] > eta_hi)  {eta_hi = (*data)->eta[smap->jPou[ii]];}
        if ((*data)->bottom[smap->jPou[ii]] > bot_hi)   {bot_hi = (*data)->bottom[smap->jPou[ii]];}
        (*data)->depty[smap->jPou[ii]] = eta_hi - bot_hi;
    }
    for (ii = 0; ii < param->ny; ii++)
    {
        // xm
        eta_hi = (*data)->eta[smap->iMin[ii]];
        bot_hi = (*data)->bottom[smap->iMin[ii]];
        if ((*data)->eta[smap->iMou[ii]] > eta_hi)  {eta_hi = (*data)->eta[smap->iMou[ii]];}
        if ((*data)->bottom[smap->iMou[ii]] > bot_hi)   {bot_hi = (*data)->bottom[smap->iMou[ii]];}
        (*data)->deptx[smap->iMou[ii]] = eta_hi - bot_hi;
        // xp
        eta_hi = (*data)->eta[smap->iPin[ii]];
        bot_hi = (*data)->bottom[smap->iPin[ii]];
        if ((*data)->eta[smap->iPou[ii]] > eta_hi)  {eta_hi = (*data)->eta[smap->iPou[ii]];}
        if ((*data)->bottom[smap->iPou[ii]] > bot_hi)   {bot_hi = (*data)->bottom[smap->iPou[ii]];}
        (*data)->deptx[smap->iPou[ii]] = eta_hi - bot_hi;
    }
    // zero depth at outer boundaries
    // To Be implemented, ZhiLi20200622

    // zero negative depth
    for (ii = 0; ii < param->n2ct; ii++)
    {
        if ((*data)->deptx[ii] < 0)    {(*data)->deptx[ii] = 0.0;}
        if ((*data)->depty[ii] < 0)    {(*data)->depty[ii] = 0.0;}
    }

}


// >>>>> Initial condition for groundwater solver <<<<<
void ic_subsurface(Data **data, Map *gmap, Config *param, int irank, int nrank)
{
    int ii, kk;
    double h_incre, zwt;
    // subsurface boundary conditions
    (*data)->qbot = param->qbot;
    (*data)->htop = param->htop;
    (*data)->hbot = param->hbot;

    for (ii = 0; ii < param->n2ci; ii++)
    {
        (*data)->qtop[ii] = ((*data)->evap[ii] - (*data)->rain[0]) / param->wcs;
    }

    // initialize soil properties
    for (ii = 0; ii < param->n3ct; ii++)
    {
        // homogeneous soil properties for now, Zhi Li 20200724
        (*data)->vga[ii] = param->soil_a;
        (*data)->vgn[ii] = param->soil_n;
        (*data)->wcr[ii] = param->wcr;
        (*data)->wcs[ii] = param->wcs;
        (*data)->Ksx[ii] = param->Ksx;
        (*data)->Ksy[ii] = param->Ksy;
        (*data)->Ksz[ii] = param->Ksz;
        // Maina heterogeneous problem
        // if (gmap->kk[ii] > 600 & gmap->kk[ii] < 1200)
        // if (gmap->kk[ii] > 120 & gmap->kk[ii] < 240)
        // {
        //     (*data)->vga[ii] = 1.04;
        //     (*data)->vgn[ii] = 1.395;
        //     (*data)->wcr[ii] = 0.106;
        //     (*data)->wcs[ii] = 0.469;
        //     (*data)->Ksz[ii] = 0.00000151;
        // }
    }
    // if init_wc within [wcr, wcs], initialize domain with init_wc
    if (param->init_wc >= param->wcr & param->init_wc <= param->wcs)
    {
        for (ii = 0; ii < param->n3ct; ii++)
        {(*data)->wc[ii] = param->init_wc;}
        // if fully saturated, use hydrostatic h
        if (param->init_wc == param->wcs)
        {
            for (ii = 0; ii < param->n3ct; ii++)
            {
                if (gmap->actv[ii] == 1)
                {
                    h_incre = (*data)->bottom[gmap->top2d[ii]] - gmap->bot3d[ii] - 0.5*gmap->dz3d[ii];
                    (*data)->h[ii] = (*data)->dept[gmap->top2d[ii]] + h_incre;
                }
                else
                {(*data)->h[ii] = 0.0;}
            }
        }
        // if unsaturated, get h from retention curve
        else
        {
            for (ii = 0; ii < param->n3ct; ii++)
            {
                (*data)->h[ii] = compute_hwc(*data, ii, param);
                // h_incre = (*data)->bottom[gmap->top2d[ii]] - gmap->bot3d[ii] - 0.5*gmap->dz3d[ii];
                // (*data)->h[ii] = compute_hwc(*data, ii, param) + (*data)->dept[gmap->top2d[ii]] + h_incre;
            }
        }
    }
    // else if init_h < 0, initialize domain with init_h
    else if (param->init_h <= 0)
    {
        for (ii = 0; ii < param->n3ct; ii++)
        {
            (*data)->h[ii] = param->init_h;
            (*data)->wc[ii] = compute_wch(*data, ii, param);
        }
    }
    // else, initialize with prescribed water table
    else
    {
        // if init_wt_rel > 0, wt is relative to surface
        if (param->init_wt_rel > 0)
        {
            for (ii = 0; ii < param->n3ct; ii++)
            {
                zwt = (*data)->bottom[gmap->top2d[ii]] - param->init_wt_rel;
                if (gmap->bot3d[ii] < zwt)
                {
                    (*data)->wc[ii] = param->wcs;
                    (*data)->h[ii] = zwt - gmap->bot3d[ii] - 0.5*gmap->dz3d[ii];
                }
                else
                {
                    // (*data)->h[ii] = zwt - gmap->bot3d[ii] - 0.5*gmap->dz3d[ii];
                    // (*data)->wc[ii] = compute_wch(*data, ii, param);
                    // (*data)->wc[ii] = param->wcr +
                        // (param->wcs-param->wcr)*((*data)->bottom[gmap->top2d[ii]]-gmap->bot3d[ii])/zwt + 0.01;

                    (*data)->wc[ii] = param->wcr + 0.02;
                    // (*data)->wc[ii] = 0.38;
                    (*data)->h[ii] = compute_hwc(*data, ii, param);
                }
            }
        }
        // else, wt is at fixed elevation
        else
        {
            for (ii = 0; ii < param->n3ct; ii++)
            {
                if (gmap->bot3d[ii] < param->init_wt_abs)
                {
                    (*data)->wc[ii] = param->wcs;
                    (*data)->h[ii] = param->init_wt_abs - gmap->bot3d[ii] - 0.5*gmap->dz3d[ii];
                }
                else
                {
                    (*data)->h[ii] = param->init_wt_abs - gmap->bot3d[ii] - 0.5*gmap->dz3d[ii];
                    (*data)->wc[ii] = compute_wch(*data, ii, param);
                }
                // if (gmap->jj[ii] == 10)
                // {
                //     printf(" !!!!! kk=%d, h->wc : %f->%f, bot=%f, dept=%f \n",gmap->kk[ii],(*data)->h[ii],(*data)->wc[ii],gmap->bot3d[ii],(*data)->dept[gmap->top2d[ii]]);
                // }
            }
        }
    }

    // read initial condidtions from file
    if (param->h_file == 1)
    {restart_subsurface((*data)->h, "head_ic", param, irank);}
    if (param->wc_file == 1)
    {restart_subsurface((*data)->wc, "moisture_ic", param, irank);}

    // fixed head boundary condition
    for (ii = 0; ii < param->n3ci; ii++)
    {
        if (gmap->istop[ii] == 1 & param->bctype_GW[5] == 1)
        {(*data)->h[gmap->icjckM[ii]] = (*data)->htop;}
        if (gmap->kk[ii] == param->nz-1 & param->bctype_GW[4] == 1)
        {(*data)->h[gmap->icjckP[ii]] = (*data)->hbot;}
        (*data)->vloss[ii] = 0.0;
    }
    // initialize other relavant fields
    for (ii = 0; ii < param->n3ct; ii++)
    {
        (*data)->hn[ii] = (*data)->h[ii];
        (*data)->hp[ii] = (*data)->h[ii];
        (*data)->hwc[ii] = compute_hwc(*data, ii, param);
        (*data)->wcn[ii] = (*data)->wc[ii];
        (*data)->wcp[ii] = (*data)->wc[ii];
        (*data)->wch[ii] = compute_wch(*data, ii, param);
        (*data)->Kx[ii] = compute_K(*data, (*data)->Ksx, ii, param);
        (*data)->Ky[ii] = compute_K(*data, (*data)->Ksy, ii, param);
        (*data)->Kz[ii] = compute_K(*data, (*data)->Ksz, ii, param);
        (*data)->qx[ii] = 0.0;
        (*data)->qy[ii] = 0.0;
        (*data)->qz[ii] = 0.0;
        (*data)->r_rho[ii] = 1.0;
        (*data)->r_rhon[ii] = 1.0;
        (*data)->r_rhoxp[ii] = 1.0;
        (*data)->r_rhoyp[ii] = 1.0;
        (*data)->r_rhozp[ii] = 1.0;
        (*data)->r_visc[ii] = 1.0;
        (*data)->r_viscxp[ii] = 1.0;
        (*data)->r_viscyp[ii] = 1.0;
        (*data)->r_visczp[ii] = 1.0;
        (*data)->Vg[ii] = param->dx*param->dy*gmap->dz3d[ii]*(*data)->wc[ii];
        (*data)->Vgn[ii] = param->dx*param->dy*gmap->dz3d[ii]*(*data)->wc[ii];
    }

    // initial scalar
    if (param->n_scalar > 0 & param->sim_groundwater == 1)
    {
        for (kk = 0; kk < param->n_scalar; kk++)
        {
            if (param->scalar_subs_file[kk] == 1)
            {
                char fscaname[50], fid[2];
                strcpy(fscaname, "scalar_subs_ic");
                sprintf(fid, "%d", kk+1);
                strcat(fscaname, fid);
                restart_subsurface((*data)->s_subs[kk], fscaname, param, irank);
            }
            else
            {

                for (ii = 0; ii < param->n3ct; ii++)
                {
                    (*data)->s_subs[kk][ii] = param->init_s_subs[kk];
                    // if (gmap->jj[ii] == 2)   {(*data)->s_subs[kk][ii] = 25.0;}
                }
                // ex5_baroclinic1d, ZhiLi20210411
                // (*data)->s_subs[kk][0] = 25.0;
            }
            for (ii = 0; ii < param->n3ct; ii++)
            {(*data)->sm_subs[kk][ii] = (*data)->s_subs[kk][ii] * (*data)->Vg[ii];}
            for (ii = 0; ii < param->n3ci; ii++)
            {
                if (gmap->istop[ii] == 1 & param->sim_shallowwater == 1)
                {
                    (*data)->s_subs[kk][gmap->icjckM[ii]] = (*data)->s_surf[kk][gmap->top2d[ii]];
                    (*data)->sm_subs[kk][gmap->icjckM[ii]] = (*data)->s_subs[kk][gmap->icjckM[ii]] * (*data)->Vg[ii];
                    (*data)->s_surfkP[kk][gmap->top2d[ii]] = (*data)->s_subs[kk][ii];
                }
            }
            // density coefficients
            if (param->baroclinic == 1 & param->n_scalar == 1 & kk == 0)
            {
                for (ii = 0; ii < param->n3ct; ii++)
                {
                    // (*data)->r_rho[ii] = 1.0 + (*data)->s_subs[0][ii] * 0.00065;
                    // (*data)->r_visc[ii] = 1.0 - (*data)->s_subs[0][ii] * 0.0015;
                    (*data)->r_rho[ii] = 1.0 + (*data)->s_subs[0][ii] * 0.00078;
                    (*data)->r_visc[ii] = 1.0 / (1.0 + (*data)->s_subs[0][ii] * 0.0022);
                    (*data)->r_rhon[ii] = (*data)->r_rho[ii];
                }
            }
        }
    }
    if (param->use_mpi == 1 & param->sim_groundwater == 1)
    {
        mpi_exchange_subsurf((*data)->wc, gmap, 2, param, irank, nrank);
        mpi_exchange_subsurf((*data)->h, gmap, 2, param, irank, nrank);
        if (param->n_scalar > 0)
        {
            for (kk = 0; kk < param->n_scalar; kk++)
            {
                mpi_exchange_subsurf((*data)->s_subs[kk], gmap, 2, param, irank, nrank);
                mpi_exchange_subsurf((*data)->sm_subs[kk], gmap, 2, param, irank, nrank);
            }
        }
    }
}

// >>>>> Read subsurface initial condition from file
void restart_subsurface(double *ic_array, char *fname, Config *param, int irank)
{
    int ii, jj, kk, xrank, yrank, col, row, count;
    char fullname[150];
    double *ic_root = malloc(param->N3CI*sizeof(double));
    strcpy(fullname, param->finput);
    strcat(fullname, fname);
    load_data(ic_root, fullname, param->N3CI);

    count = 0;
    for (ii = 0; ii < param->n2ci; ii++)
    {
        xrank = irank % param->mpi_nx;
        yrank = floor(irank/param->mpi_nx);
        col = floor(ii / param->nx);
        row = ii % param->nx;
        jj = yrank*param->mpi_nx*param->n2ci + col*param->NX + xrank*param->nx + row;
        for (kk = 0; kk < param->nz; kk++)
        {
            ic_array[count] = ic_root[jj*param->nz + kk];
            count += 1;
        }
    }
}


// >>>>> Initialize subgrid variables
void init_subgrid(Data **data, Map *smap, Config *param, int irank, int nrank)
{
    int ii, jj, kk, offind;
    char fullname[100];
    // load subgrid variables and assign to each rank
    strcpy(fullname, param->finput);
    strcat(fullname, "Vxp");
    load_data((*data)->Vxp_sub_root, fullname, param->N2CI*param->nlay_sub);
    strcpy(fullname, param->finput);
    strcat(fullname, "Vxm");
    load_data((*data)->Vxm_sub_root, fullname, param->N2CI*param->nlay_sub);
    strcpy(fullname, param->finput);
    strcat(fullname, "Vyp");
    load_data((*data)->Vyp_sub_root, fullname, param->N2CI*param->nlay_sub);
    strcpy(fullname, param->finput);
    strcat(fullname, "Vym");
    load_data((*data)->Vym_sub_root, fullname, param->N2CI*param->nlay_sub);
    strcpy(fullname, param->finput);
    strcat(fullname, "Axp");
    load_data((*data)->Axp_sub_root, fullname, param->N2CI*param->nlay_sub);
    strcpy(fullname, param->finput);
    strcat(fullname, "Axm");
    load_data((*data)->Axm_sub_root, fullname, param->N2CI*param->nlay_sub);
    strcpy(fullname, param->finput);
    strcat(fullname, "Ayp");
    load_data((*data)->Ayp_sub_root, fullname, param->N2CI*param->nlay_sub);
    strcpy(fullname, param->finput);
    strcat(fullname, "Aym");
    load_data((*data)->Aym_sub_root, fullname, param->N2CI*param->nlay_sub);
    strcpy(fullname, param->finput);
    strcat(fullname, "Asz");
    load_data((*data)->Asz_sub_root, fullname, param->N2CI*param->nlay_sub);
    for (ii = 0; ii < param->nlay_sub; ii++)
    {
        offind = ii * param->N2CI;
        root_to_rank((*data)->Vxp_sub_root, (*data)->Vxp_sub[ii], param, irank, nrank, offind, 0.0);
        root_to_rank((*data)->Vxm_sub_root, (*data)->Vxm_sub[ii], param, irank, nrank, offind, 0.0);
        root_to_rank((*data)->Vyp_sub_root, (*data)->Vyp_sub[ii], param, irank, nrank, offind, 0.0);
        root_to_rank((*data)->Vym_sub_root, (*data)->Vym_sub[ii], param, irank, nrank, offind, 0.0);
        root_to_rank((*data)->Axp_sub_root, (*data)->Axp_sub[ii], param, irank, nrank, offind, 0.0);
        root_to_rank((*data)->Axm_sub_root, (*data)->Axm_sub[ii], param, irank, nrank, offind, 0.0);
        root_to_rank((*data)->Ayp_sub_root, (*data)->Ayp_sub[ii], param, irank, nrank, offind, 0.0);
        root_to_rank((*data)->Aym_sub_root, (*data)->Aym_sub[ii], param, irank, nrank, offind, 0.0);
        root_to_rank((*data)->Asz_sub_root, (*data)->Asz_sub[ii], param, irank, nrank, offind, 0.0);
    }
    // initialize index of eta
    for (ii = 0; ii < param->n2ci; ii++)
    {
        if ((*data)->eta[ii] < (*data)->layers_sub[0])
        {(*data)->eta_ind[ii] = 0;}
        else if ((*data)->eta[ii] >= (*data)->layers_sub[param->nlay_sub-1])
        {(*data)->eta_ind[ii] = param->nlay_sub-1;}
        else
        {
            for (jj = 0; jj < param->nlay_sub-1; jj++)
            {
                if ((*data)->eta[ii] >= (*data)->layers_sub[jj] & (*data)->eta[ii] < (*data)->layers_sub[jj+1])
                {
                    (*data)->eta_ind[ii] = jj;
                    break;
                }
            }
        }
        // error checking
        // if ((*data)->dept[ii] > 0)
        // {
        //     if ((*data)->eta_ind[ii] == 0 || (*data)->eta_ind[ii] == param->nlay_sub-1)
        //     {
        //         printf("WARNING: Surface elevation exceeds range of look-up table!\n");
        //         printf("        ---> at rank %d, cell %d (1D), eta=%f\n", irank, ii, (*data)->eta[ii]);
        //     }
        // }
    }
    // extract and combine subgrid variables for a given surface elevation
    subgrid_interp_and_combine(data, smap, param, irank, nrank);
    // subgrid variables along boundaries
    for (ii = 0; ii < param->nx; ii++)
    {
        (*data)->Vs[smap->jMou[ii]] = (*data)->Vs[smap->jMin[ii]];
        (*data)->Vs[smap->jPou[ii]] = (*data)->Vs[smap->jPin[ii]];
        (*data)->Asx[smap->jMou[ii]] = (*data)->Asx[smap->jMin[ii]];
        (*data)->Asx[smap->jPou[ii]] = (*data)->Asx[smap->jPin[ii]];
        (*data)->Asy[smap->jMou[ii]] = (*data)->Asy[smap->jMin[ii]];
        (*data)->Asy[smap->jPou[ii]] = (*data)->Asy[smap->jPin[ii]];
        (*data)->Asz[smap->jMou[ii]] = (*data)->Asz[smap->jMin[ii]];
        (*data)->Asz[smap->jPou[ii]] = (*data)->Asz[smap->jPin[ii]];
    }
    for (ii = 0; ii < param->ny; ii++)
    {
        (*data)->Vs[smap->iMou[ii]] = (*data)->Vs[smap->iMin[ii]];
        (*data)->Vs[smap->iPou[ii]] = (*data)->Vs[smap->iPin[ii]];
        (*data)->Asx[smap->iMou[ii]] = (*data)->Asx[smap->iMin[ii]];
        (*data)->Asx[smap->iPou[ii]] = (*data)->Asx[smap->iPin[ii]];
        (*data)->Asy[smap->iMou[ii]] = (*data)->Asy[smap->iMin[ii]];
        (*data)->Asy[smap->iPou[ii]] = (*data)->Asy[smap->iPin[ii]];
        (*data)->Asz[smap->iMou[ii]] = (*data)->Asz[smap->iMin[ii]];
        (*data)->Asz[smap->iPou[ii]] = (*data)->Asz[smap->iPin[ii]];
    }
    for (ii = 0; ii < param->nx; ii++)
    {
        (*data)->Vsy[smap->jMou[ii]] = (*data)->Vs[smap->jMin[ii]];
        (*data)->Aszy[smap->jMou[ii]] = (*data)->Asz[smap->jMin[ii]];
    }
    for (ii = 0; ii < param->ny; ii++)
    {
        (*data)->Vsx[smap->iMou[ii]] = (*data)->Vs[smap->iMin[ii]];
        (*data)->Aszx[smap->iMou[ii]] = (*data)->Asz[smap->iMin[ii]];
    }
    for (ii = 0; ii < param->n2ct; ii++)    {(*data)->Vsn[ii] = (*data)->Vs[ii];}

}
