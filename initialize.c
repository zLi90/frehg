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
#include"solve.h"
#include"utility.h"


void init(Data **data, Map **smap, Map **gmap, Config **param, int irank, int nrank);
void init_domain(Config **param);
void init_Data(Data **data, Config *param);
void read_bathymetry(Data **data, Config *param, int irank, int nrank);
void boundary_bath(Data **data, Map *smap, Config *param, int irank, int nrank);
void ic_surface(Data **data, Map *smap, Config *param, int irank, int nrank);
void bc_surface(Data **data, Map *smap, Config *param, int irank);
void update_depth(Data **data, Map *smap, Config *param, int irank);
void ic_subsurface(Data **data, Map *gmap, Config *param);

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
    // initial condition for shallow water solver
    ic_surface(data, *smap, *param, irank, nrank);
    update_drag_coef(data, *param);
    // initial condition for groundwater solver
    ic_subsurface(data, *gmap, *param);
    // boundary condition for shallow water solver
    enforce_head_bc(data, *gmap, *param);

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
    char fullname[20];
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
        if (param->use_subgrid == 0)
        {
            // load the bathymetry file
            strcpy(fullname, "bath");
            load_data((*data)->bottom_root, fullname, param->N2CI);
            // calculate bathymetry offset
            z_min = getMin((*data)->bottom_root, param->N2CI);
            if (z_min >= 0) {(*data)->offset[0] = 0;}
            else    {(*data)->offset[0] = -z_min;}
            // bathymetry for each rank
            for (ii = 0; ii < param->n2ci; ii++)
            {
                xrank = irank % param->mpi_nx;
                yrank = floor(irank/param->mpi_nx);
                // jj = yrank * param->mpi_nx + xrank;
                // (*data)->bottom[ii] = (*data)->bottom_root[jj*param->n2ci+ii] + (*data)->offset[0];
                col = floor(ii / param->nx);
                row = ii % param->nx;
                jj = yrank*param->mpi_nx*param->n2ci + col*param->NX + xrank*param->nx + row;
                (*data)->bottom[ii] = (*data)->bottom_root[jj] + (*data)->offset[0];
            }
        }
        else
        {
            mpi_print(" >>> ERROR: Load subgrid bathymetry is disabled!", irank);
        }
    }
    else
    {(*data)->bottom[ii] = 0.0;}
}

// >>>>> Set boundary bathymetry <<<<<
void boundary_bath(Data **data, Map *smap, Config *param, int irank, int nrank)
{
    int ii;
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
    for (ii = 0; ii < param->n2ci; ii++)
    {
        (*data)->bottomXP[ii] = (*data)->bottom[ii];
        if ((*data)->bottom[smap->iPjc[ii]] > (*data)->bottom[ii])
        {(*data)->bottomXP[ii] = (*data)->bottom[smap->iPjc[ii]];}
        (*data)->bottomYP[ii] = (*data)->bottom[ii];
        if ((*data)->bottom[smap->icjP[ii]] > (*data)->bottom[ii])
        {(*data)->bottomYP[ii] = (*data)->bottom[smap->icjP[ii]];}
    }
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
    (*data)->uy = malloc(param->n2ct*sizeof(double));
    (*data)->vv = malloc(param->n2ct*sizeof(double));
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
    (*data)->Vsx = malloc(param->n2ct*sizeof(double));
    (*data)->Vsy = malloc(param->n2ct*sizeof(double));
    (*data)->Asz = malloc(param->n2ct*sizeof(double));
    (*data)->Aszx = malloc(param->n2ct*sizeof(double));
    (*data)->Aszy = malloc(param->n2ct*sizeof(double));
    (*data)->Asx = malloc(param->n2ct*sizeof(double));
    (*data)->Asy = malloc(param->n2ct*sizeof(double));
    (*data)->wtfx = malloc(param->n2ct*sizeof(double));
    (*data)->wtfy = malloc(param->n2ct*sizeof(double));

    (*data)->uu_root = malloc(n2_root*sizeof(double));
    (*data)->vv_root = malloc(n2_root*sizeof(double));
    (*data)->eta_root = malloc(n2_root*sizeof(double));
    (*data)->dept_root = malloc(n2_root*sizeof(double));
    (*data)->seep_root = malloc(n2_root*sizeof(double));
    (*data)->uu_out = malloc(n2_root*sizeof(double));
    (*data)->vv_out = malloc(n2_root*sizeof(double));
    (*data)->eta_out = malloc(n2_root*sizeof(double));
    (*data)->dept_out = malloc(n2_root*sizeof(double));
    (*data)->seep_out = malloc(n2_root*sizeof(double));

    (*data)->reset_seepage = malloc(param->n2ci*sizeof(int));
    (*data)->qseepage = malloc(param->n2ci*sizeof(double));
    (*data)->cflx = malloc(param->n2ci*sizeof(double));
    (*data)->cfly = malloc(param->n2ci*sizeof(double));
    (*data)->cfl_active = malloc(param->n2ci*sizeof(double));

    // subsurface fields
    (*data)->repeat = malloc(1*sizeof(int));
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
    (*data)->qx = malloc(param->n3ct*sizeof(double));
    (*data)->qy = malloc(param->n3ct*sizeof(double));
    (*data)->qz = malloc(param->n3ct*sizeof(double));
    (*data)->Vg = malloc(param->n3ct*sizeof(double));
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
void ic_surface(Data **data, Map *smap, Config *param, int irank, int nrank)
{
    int ii, jj, col, row, xrank, yrank;
    // initial boundary condition
    (*data)->tide = malloc(param->n_tide*sizeof(double));
    get_tide(data, param);
    (*data)->rain = malloc(1*sizeof(double));
    (*data)->rain_sum = malloc(1*sizeof(double));
    (*data)->evap = malloc(1*sizeof(double));
    (*data)->rain_sum[0] = 0.0;
    get_evaprain(data, param);
    // initial surface elevation
    if (param->eta_file == 0)
    {
        for (ii = 0; ii < param->N2CI; ii++)
        {(*data)->eta_root[ii] = param->init_eta + (*data)->offset[0];}
    }
    else
    {
        char fullname[20];
        strcpy(fullname, "surf_ic");
        load_data((*data)->eta_root, fullname, param->N2CI);
    }
    // initial velocity
    if (param->uv_file == 0)
    {
        for (ii = 0; ii < param->N2CI; ii++)
        {(*data)->uu_root[ii] = 0.0;    (*data)->vv_root[ii] = 0.0;}
    }
    else
    {
        char fullname1[20];
        strcpy(fullname1, "uu_ic");
        load_data((*data)->uu_root, fullname1, param->N2CI);
        char fullname2[20];
        strcpy(fullname2, "vv_ic");
        load_data((*data)->vv_root, fullname2, param->N2CI);
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
    }
    // enforce boundary conditions
    if (param->use_mpi == 1)
    {
        mpi_exchange_surf((*data)->eta, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->uu, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->vv, smap, 2, param, irank, nrank);
    }
    enforce_surf_bc(data, smap, param, irank, nrank);
    // calculate depth
    update_depth(data, smap, param, irank);
    if (param->use_mpi == 1)
    {
        mpi_exchange_surf((*data)->dept, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->deptx, smap, 2, param, irank, nrank);
        mpi_exchange_surf((*data)->depty, smap, 2, param, irank, nrank);
    }
    // calculate subgrid areas
    if (param->use_subgrid == 1)
    {
        mpi_print("WARNING: subgrid functions have not been implemented!",0);
    }
    else
    {
        for (ii = 0; ii < param->n2ct; ii++)
        {
            (*data)->Vs[ii] = (*data)->dept[ii] * param->dx * param->dy;
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
    if (param->use_mpi == 1)
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
    // initialize other relevant variables
    for (ii = 0; ii < param->n2ci; ii++)
    {
        (*data)->qseepage[ii] = 0.0;
        (*data)->reset_seepage[ii] = 1;
    }
}

// >>>>> Boundary condition for shallowwater solver <<<<<
void bc_surface(Data **data, Map *smap, Config *param, int irank)
{
    if (exist("tide1"))
    {load_bc((*data)->tide1, (*data)->t_tide1, "tide1", param->n_tide1);}
    if (exist("tide2"))
    {load_bc((*data)->tide2, (*data)->t_tide2, "tide2", param->n_tide2);}

    if (param->inflow_loc[0] >= 0 & param->inflow_loc[1] >= 0 & param->inflow_loc[2] >= 0)
    {load_bc((*data)->inflow, (*data)->t_inflow, "inflow", param->n_inflow);}

}


// >>>>> Calculate surface depth from bathymetry and surface elevation <<<<<
void update_depth(Data **data, Map *smap, Config *param, int irank)
{
    int ii;
    double eta_hi, bot_hi;
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
void ic_subsurface(Data **data, Map *gmap, Config *param)
{
    int ii;
    double h_incre, zwt;
    (*data)->qtop = param->qtop;
    (*data)->qbot = param->qbot;
    (*data)->htop = param->htop;
    (*data)->hbot = param->hbot;
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
                    h_incre = (*data)->bottom[gmap->top2d[ii]] - gmap->bot3d[ii] + 0.5*gmap->dz3d[ii];
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
            {(*data)->h[ii] = compute_hwc(*data, ii, param);}
        }
    }
    // else if init_h < 0, initialize domain with init_h
    else if (param->init_h < 0)
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
                    (*data)->wc[ii] = param->wcr + 0.03;
                    // (*data)->wc[ii] = 0.36;
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
            }
        }
    }
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
        (*data)->Vg[ii] = param->dx*param->dy*gmap->dz3d[ii]*(*data)->wc[ii];
    }

}
