// Utility Functions
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include"configuration.h"
#include"initialize.h"
#include"map.h"
#include"mpifunctions.h"

double compute_wch(Data *data, int ii, Config *param);
double compute_hwc(Data *data, int ii, Config *param);
double compute_ch(Data *data, int ii, Config *param);
double compute_K(Data *data, double *Ksat, int ii, Config *param);
double compute_dKdwc(Data *data, double *Ksat, int ii, Config *param);
char* read_one_input(char field[], char fname[]);
double read_one_input_double(char field[], char fname[]);
int * read_one_input_array(char field[], char fname[], int n);
double * read_one_input_array_double(char field[], char fname[], int n);
void load_data(double *arr, char filename[], int n);
void load_bc(double *value, double *tVec, char filename[], int n);
int get_index_sort(double *arr, double targ, double itval, int n);
void reorder_surf(double *out, double *root, Config *param);
void reorder_subsurf(double *out, double *root, Map *gmap, Config *param);
void write_output(Data **data, Map *gmap, Config *param, int tt, int root, int irank);
void write_one_file(double *ally, char *filename, Config *param, int tt, int n);
void append_to_file(char *filename, double val, Config *param);
double interp_bc(double *tVec, double *value, double t_current, int n_dat);
void mpi_print(char pstr[], int irank);
double getMin(double *arr, int n);
double getMax(double *arr, int n);
int exist(char *fname);

// >>>>> Compute water content from h using water retention curve <<<<<
double compute_wch(Data *data, int ii, Config *param)
{
    double wc, m, s, h, wcs, wcm;
    h = data->h[ii];
    m = 1.0 - 1.0/data->vgn[ii];
    if (param->use_mvg == 1)
    {wcm = data->wcr[ii] + (data->wcs[ii]-data->wcr[ii])*pow((1.0 + pow(fabs(param->aev)*data->vga[ii],data->vgn[ii])), m);}
    else
    {wcm = data->wcs[ii];}
    s = pow(1.0 + pow(fabs(data->vga[ii]*h), data->vgn[ii]), -m);

    if (h > param->aev)  {wc = data->wcs[ii];}
    else    {wc = data->wcr[ii] + (wcm-data->wcr[ii]) * s;}
    if (wc > data->wcs[ii])    {wc = data->wcs[ii];}
    else if (wc < data->wcr[ii])   {wc = data->wcr[ii];}

    return wc;
}

// >>>>> Compute h from water content using water retention curve <<<<<
double compute_hwc(Data *data, int ii, Config *param)
{
    double wc, wcs, h, m, wcm, eps = 1e-7;
    wc = data->wc[ii];
    m = 1.0 - 1.0/data->vgn[ii];

    if (param->use_mvg == 1)
    {wcm = data->wcr[ii] + (data->wcs[ii]-data->wcr[ii])*pow((1.0 + pow(fabs(param->aev)*data->vga[ii],data->vgn[ii])), m);}
    else
    {wcm = data->wcs[ii];}

    if (wc - data->wcr[ii] < eps)  {wc = data->wcr[ii] + eps;}
    if (wc < data->wcs[ii])
    {
        h = -(1.0/data->vga[ii]) *
            pow(pow((wcm - data->wcr[ii])/(wc - data->wcr[ii]),(1.0/m)) - 1.0,(1.0/data->vgn[ii]));
    }
    else    {h = 0.0;}

    if (isnan(h))
    {
        printf("h, wcm, wc = %f, %f, %f\n",h,wcm,wc);
    }
    return h;
}

// >>>>> Compute specific capacity <<<<<
double compute_ch(Data *data, int ii, Config *param)
{
    double h, wcs, c, m, wcm, deno, nume;
    h = data->h[ii];
    m = 1.0 - 1.0/data->vgn[ii];
    if (param->use_mvg == 1)
    {wcm = data->wcr[ii] + (data->wcs[ii]-data->wcr[ii])*pow((1.0 + pow(fabs(param->aev)*data->vga[ii],data->vgn[ii])), m);}
    else
    {wcm = data->wcs[ii];}

    nume = data->vga[ii]*data->vgn[ii]*m * (wcm-data->wcr[ii]) *
        pow(fabs(data->vga[ii]*h),data->vgn[ii]-1);
    deno = pow((1.0 + pow(fabs(data->vga[ii]*h),data->vgn[ii])), m+1);
    c = nume / deno;
    if (param->use_mvg == 1)
    {if (h > param->aev)  {c = 0.0;}}
    else
    {if (h > 0.0)  {c = 0.0;}}
    return c;
}

// >>>>> Compute hydraulic conductivity <<<<<
double compute_K(Data *data, double *Ksat, int ii, Config *param)
{
    double m, wcs, Keff, s, h, Ks, wcm, nume, deno, ratio;
    Ks = Ksat[ii];
    h = data->h[ii];
    m = 1.0 - 1.0/data->vgn[ii];
    s = pow(1.0 + pow(fabs(data->vga[ii]*h), data->vgn[ii]), -m);
    if (param->use_mvg == 1)
    {wcm = data->wcr[ii] + (data->wcs[ii]-data->wcr[ii])*pow((1.0 + pow(fabs(param->aev)*data->vga[ii],data->vgn[ii])), m);}
    else
    {wcm = data->wcs[ii];}
    nume = 1.0-pow(1.0-pow(s*(data->wcs[ii]-data->wcr[ii])/(wcm-data->wcr[ii]),1.0/m),m);
    deno = 1.0-pow(1.0-pow((data->wcs[ii]-data->wcr[ii])/(wcm-data->wcr[ii]),1.0/m),m);

    if (param->use_mvg == 1)
    {
        if (deno == 0.0)    {Keff = Ks;}
        else    {Keff = Ks * pow(s,0.5) * pow(nume/deno, 2.0);}
    }
    else
    {Keff = Ks * pow(s,0.5) * pow(1-pow(1-pow(s,1.0/m),m), 2.0);}
    if (Keff > Ks)  {Keff = Ks;}
    if (param->use_mvg == 1)
    {if (h > param->aev)  {Keff = Ks;}}
    else
    {if (h > 0.0)  {Keff = Ks;}}

    if (isnan(Keff))
    {
        printf("ii=%d, nume=%f, deno=%f, h=%f, s=%f, (wcs,wcr,wcm)=(%f,%f,%f)\n",
            ii,nume,deno,h,s,data->wcs[ii],data->wcr[ii],wcm);
    }

    return Keff;
}

// >>>>> Compute dKdwc for adaptive time stepping <<<<<
double compute_dKdwc(Data *data, double *Ksat, int ii, Config *param)
{
    double term1, term2, term0, m, s, dKdwc, wc, Ks, wcm, c1, c2;
    Ks = Ksat[ii];
    wc = data->wc[ii];
    if (param->use_mvg == 1)
    {wcm = data->wcr[ii] + (data->wcs[ii]-data->wcr[ii])*pow((1.0 + pow(fabs(param->aev)*data->vga[ii],data->vgn[ii])), m);}
    else
    {wcm = data->wcs[ii];}

    // use the lambda limiter
    if (wc > 0.9999 * data->wcs[ii] & wc < data->wcs[ii])
    {wc = 0.9999 * data->wcs[ii];}
    m = 1.0 - 1.0/data->vgn[ii];
    s = (wc - data->wcr[ii]) / (data->wcs[ii] - data->wcr[ii]);
    if (param->use_mvg == 0)
    {
        // calculate dKdwc
        term0 = pow(1.0 - pow(s,1.0/m), m);
        term1 = 0.5 * Ks * pow(s,-0.5) * (1.0 - term0) * (1.0 - term0);
        term2 = 2.0 * Ks * pow(s,(2.0-m)/(2.0*m)) * (1.0 - term0) * pow(1.0 - pow(s,1.0/m), m-1);
    }
    else
    {
        c2 = (data->wcs[ii] - data->wcr[ii])/(wcm - data->wcr[ii]);
        c1 = 1.0 / pow(1.0 - pow((1.0 - pow(c2,1.0/m)),m), 2.0);
        term0 = pow(1.0 - pow(c2*s,1.0/m), m);
        term1 = 0.5 * c1 * Ks * pow(s,-0.5) * (1.0 - term0) * (1.0 - term0);
        term2 = 2.0 * c1 * Ks * c2 * pow(s,0.5) * pow(c2*s,1.0/m-1.0) * (1.0 - term0) * pow(1.0 - pow(c2*s,1.0/m), m-1);
    }
    dKdwc = (term1 + term2) / (data->wcs[ii] - data->wcr[ii]);
    return dKdwc;
}

// >>>>> Read one user input field <<<<<
char* read_one_input(char field[], char fname[])
{
    int n = 100, found = 0;
    FILE *fid;
    char *elem, *ptr, arr[n], *out;
    size_t len;
    fid = fopen(fname, "r");
    while (fgets(arr, n, fid) != NULL)
    {
        char *elem = strtok(arr, " ");
        if (strcmp(elem, field) == 0)
        {
            found = 1;
            elem = strtok(NULL, " ");
            out = strtok(NULL, " ");
            break;
        }
    }
    if (found == 0)
    {printf("WARNING: Field %s not found in the input file!\n", field); out = "0";}
    else
    {len = strlen(out);     out[len-1] = 0;}
    return out;
}

// >>>>> Convert user input field to double <<<<<
double read_one_input_double(char field[], char fname[])
{
    char *ptr, out_str[50];
    strcpy(out_str, read_one_input(field, fname));
    return strtod(out_str, &ptr);
}

// >>>>> Read one array from user input <<<<<
int * read_one_input_array(char field[], char fname[], int n)
{
    char *ptr, *elem, out_str[50];
    int ii = 0, *out_arr;
    out_arr = malloc(n*sizeof(int));
    strcpy(out_str, read_one_input(field, fname));
    elem = strtok(out_str, ",");
    while (elem != NULL)
    {
        out_arr[ii] = strtod(elem, &ptr);
        ii += 1;
        elem = strtok(NULL, ",");
    }
    return out_arr;
}

double * read_one_input_array_double(char field[], char fname[], int n)
{
    char *ptr, *elem, out_str[50];
    int ii = 0;
    double *out_arr;
    out_arr = malloc(n*sizeof(double));
    strcpy(out_str, read_one_input(field, fname));
    elem = strtok(out_str, ",");
    while (elem != NULL)
    {
        out_arr[ii] = strtod(elem, &ptr);
        ii += 1;
        elem = strtok(NULL, ",");
    }
    return out_arr;
}

// >>>>> Read input data from data file <<<<<
void load_data(double *arr, char filename[], int n)
{
    int ii;
    FILE *fid;
    fid = fopen(filename, "r");
    if (fid == NULL)
    {printf("WARNING: Unable to open the data file: %s! \n",filename);}
    for (ii = 0; ii < n; ii++)  {fscanf(fid, "%lf,", &arr[ii]);}
    fclose(fid);
}

// >>>>> Read boundary condtion as time series data  <<<<<
void load_bc(double *value, double *tVec, char filename[], int n)
{
    int ii;
    double *bc_array = malloc(2*n*sizeof(double));
    load_data(bc_array, filename, 2*n);
    for (ii = 0; ii < n; ii++)
    {
        tVec[ii] = bc_array[2*ii];
        value[ii] = bc_array[2*ii+1];
    }
    free(bc_array);
}

// >>>>> Get index of an element in an sorted array
int get_index_sort(double *arr, double targ, double itval, int n)
{
    int count = 0, ind = 0;
    while (count < n-1)
    {
        if (arr[count] <= targ & arr[count+1] > targ)
        {
            ind = count;
            break;
        }
        count += 1;
        if (count == n-1)
        {printf("WARNING: No time index was found for BC...\n");}
    }
    return ind;
}

// >>>>> Reorder after mpi_gather <<<<<
void reorder_surf(double *out, double *root, Config *param)
{
    int ii, jj, xrank, yrank, col, istart;
    ii = 0;
    while (ii < param->N2CI)
    {
        xrank = floor((ii%param->NX) / param->nx);
        yrank = floor(ii / (param->NX*param->ny));
        col = floor(ii/param->NX) - yrank*param->ny;
        istart = xrank*param->n2ci + col*param->nx + yrank*param->n2ci*param->mpi_nx;
        for (jj = 0; jj < param->nx; jj++)
        {out[ii+jj] = root[istart+jj];}
        ii += param->nx;
    }
}

// >>>>> Reorder for subsurface domain <<<<<
void reorder_subsurf(double *out, double *root, Map *gmap, Config *param)
{
    int ii, jj, kk, ii2d, xrank, yrank, col, istart;
    ii = 0;
    while (ii < param->N3CI)
    {
        // ii2d = gmap->top2d[ii];
        ii2d = floor(ii / param->nz);
        xrank = floor((ii2d%param->NX) / param->nx);
        yrank = floor(ii2d / (param->NX*param->ny));
        col = floor(ii2d/param->NX) - yrank*param->ny;
        istart = xrank*param->n2ci + col*param->nx + yrank*param->n2ci*param->mpi_nx;
        // printf("ii, ii2d, istart, nx, nz, N3Ci = %d, %d, %d, %d, %d, %d\n",ii,ii2d,istart,param->nx, param->nz, param->N3CI);


        for (jj = 0; jj < param->nx; jj++)
        {
            for (kk = 0; kk < param->nz; kk++)
            {out[ii+jj*param->nz+kk] = root[istart*param->nz+jj*param->nz+kk];}
        }
        ii += param->nx*param->nz;
    }
}

// >>>>> Output model results <<<<<
void write_output(Data **data, Map *gmap, Config *param, int tt, int root, int irank)
{
    int ii, kk;
    char fid[2];
    // Combine all ranks at root
    if (param->use_mpi == 1)
    {
        if (param->sim_shallowwater == 1)
        {
            mpi_gather_double((*data)->eta_root, (*data)->eta, param->n2ci, root);
            mpi_gather_double((*data)->dept_root, (*data)->dept, param->n2ci, root);
            mpi_gather_double((*data)->uu_root, (*data)->uu, param->n2ci, root);
            mpi_gather_double((*data)->vv_root, (*data)->vv, param->n2ci, root);
            if (param->sim_groundwater == 1)
            {mpi_gather_double((*data)->seep_root, (*data)->qseepage, param->n2ci, root);}
            if (param->n_scalar > 0)
            {
                for (kk = 0; kk < param->n_scalar; kk++)
                {mpi_gather_double((*data)->s_surf_root[kk], (*data)->s_surf[kk], param->n2ci, root);}
            }
            if (irank == root)
            {
                for (ii = 0; ii < param->N2CI; ii++)
                {(*data)->eta_root[ii] -= (*data)->offset[0];}
                reorder_surf((*data)->eta_out, (*data)->eta_root, param);
                reorder_surf((*data)->dept_out, (*data)->dept_root, param);
                reorder_surf((*data)->uu_out, (*data)->uu_root, param);
                reorder_surf((*data)->vv_out, (*data)->vv_root, param);
                if (param->n_scalar > 0)
                {
                    for (kk = 0; kk < param->n_scalar; kk++)
                    {reorder_surf((*data)->s_surf_out[kk], (*data)->s_surf_root[kk], param);}
                }
                if (param->sim_groundwater == 1)
                {
                    for (ii = 0; ii < param->N2CI; ii++)
                    {(*data)->seep_root[ii] = (*data)->seep_root[ii]*8.64e7;}
                    reorder_surf((*data)->seep_out, (*data)->seep_root, param);
                }
            }
        }
        if (param->sim_groundwater == 1)
        {
            mpi_gather_double((*data)->h_root, (*data)->h, param->n3ci, root);
            mpi_gather_double((*data)->wc_root, (*data)->wc, param->n3ci, root);
            mpi_gather_double((*data)->qx_root, (*data)->qx, param->n3ci, root);
            mpi_gather_double((*data)->qy_root, (*data)->qy, param->n3ci, root);
            mpi_gather_double((*data)->qz_root, (*data)->qz, param->n3ci, root);
            if (param->n_scalar > 0)
            {
                for (kk = 0; kk < param->n_scalar; kk++)
                {mpi_gather_double((*data)->s_subs_root[kk], (*data)->s_subs[kk], param->n3ci, root);}
            }
            if (irank == root)
            {
                // flow rate converted to [mm/d]
                for (ii = 0; ii < param->N3CI; ii++)
                {
                    (*data)->qx_root[ii] = (*data)->qx_root[ii]*8.64e7;
                    (*data)->qy_root[ii] = (*data)->qy_root[ii]*8.64e7;
                    (*data)->qz_root[ii] = (*data)->qz_root[ii]*8.64e7;
                }
                reorder_subsurf((*data)->h_out, (*data)->h_root, gmap, param);
                reorder_subsurf((*data)->wc_out, (*data)->wc_root, gmap, param);
                reorder_subsurf((*data)->qx_out, (*data)->qx_root, gmap, param);
                reorder_subsurf((*data)->qy_out, (*data)->qy_root, gmap, param);
                reorder_subsurf((*data)->qz_out, (*data)->qz_root, gmap, param);
                if (param->n_scalar > 0)
                {
                    for (kk = 0; kk < param->n_scalar; kk++)
                    {reorder_subsurf((*data)->s_subs_out[kk], (*data)->s_subs_root[kk], gmap, param);}
                }
            }
        }

    }
    else
    {
        if (param->sim_shallowwater == 1)
        {
            for (ii = 0; ii < param->n2ci; ii++)
            {
                (*data)->eta_out[ii] = (*data)->eta[ii] - (*data)->offset[0];;
                (*data)->dept_out[ii] = (*data)->dept[ii];
                (*data)->uu_out[ii] = (*data)->uu[ii];
                (*data)->vv_out[ii] = (*data)->vv[ii];
                if (param->sim_groundwater == 1)
                {(*data)->seep_out[ii] = (*data)->qseepage[ii]*8.64e7;}
                if (param->n_scalar > 0)
                {
                    for (kk = 0; kk < param->n_scalar; kk++)
                    {(*data)->s_surf_out[kk][ii] = (*data)->s_surf[kk][ii];}
                }
            }
        }
        if (param->sim_groundwater == 1)
        {
            for (ii = 0; ii < param->n3ci; ii++)
            {
                (*data)->h_out[ii] = (*data)->h[ii];
                (*data)->wc_out[ii] = (*data)->wc[ii];
                (*data)->qx_out[ii] = (*data)->qx[ii]*8.64e7;
                (*data)->qy_out[ii] = (*data)->qy[ii]*8.64e7;
                (*data)->qz_out[ii] = (*data)->qz[ii]*8.64e7;
                if (param->n_scalar > 0)
                {
                    for (kk = 0; kk < param->n_scalar; kk++)
                    {(*data)->s_subs_out[kk][ii] = (*data)->s_subs[kk][ii];}
                }
            }
        }
    }
    // Write the output files
    if (param->sim_shallowwater == 1 & irank == root)
    {
        write_one_file((*data)->eta_out, "surf", param, tt, param->N2CI);
        write_one_file((*data)->dept_out, "depth", param, tt, param->N2CI);
        write_one_file((*data)->uu_out, "uu", param, tt, param->N2CI);
        write_one_file((*data)->vv_out, "vv", param, tt, param->N2CI);

        if (param->sim_groundwater == 1)
        {write_one_file((*data)->seep_out, "seepage", param, tt, param->N2CI);}
        if (param->n_scalar > 0)
        {
            for (kk = 0; kk < param->n_scalar; kk++)
            {
                char fullname[50];
                strcpy(fullname, "scalar_surf");
                sprintf(fid, "%d", kk+1);
                strcat(fullname, fid);
                strcat(fullname, "_");
                write_one_file((*data)->s_surf_out[kk], fullname, param, tt, param->N2CI);
            }
        }

    }
    if (param->sim_groundwater == 1 & irank == root)
    {
        write_one_file((*data)->h_out, "head", param, tt, param->N3CI);
        write_one_file((*data)->wc_out, "moisture", param, tt, param->N3CI);
        write_one_file((*data)->qx_out, "qx", param, tt, param->N3CI);
        write_one_file((*data)->qy_out, "qy", param, tt, param->N3CI);
        write_one_file((*data)->qz_out, "qz", param, tt, param->N3CI);
        if (param->n_scalar > 0)
        {
            for (kk = 0; kk < param->n_scalar; kk++)
            {
                char fullname[50];
                strcpy(fullname, "scalar_subs");
                sprintf(fid, "%d", kk+1);
                strcat(fullname, fid);
                strcat(fullname, "_");
                write_one_file((*data)->s_subs_out[kk], fullname, param, tt, param->N3CI);
            }
        }
    }
}

// >>>>> Write one output file <<<<<
void write_one_file(double *ally, char *filename, Config *param, int tt, int n)
{
    int ii;
    FILE *fp;
    // convert time step into string
    char tstr[10];
    sprintf(tstr, "%d", tt);
    // create filename for saving
    char *fullname = malloc(100);
    strcpy(fullname, param->foutput);
    strcat(fullname, filename);
    strcat(fullname, "_");
    strcat(fullname, tstr);
    // open and write to file
    fp = fopen(fullname, "w");
    for (ii = 0; ii < n; ii++)    {fprintf(fp, "%6.6f \n", ally[ii]);}
    fclose(fp);
    free(fullname);
}

// >>>>> Append to one output file
void append_to_file(char *filename, double val, Config *param)
{
    FILE *fp;
    char *fullname = malloc(100);
    strcpy(fullname, param->foutput);
    strcat(fullname, filename);
    if (exist(fullname))    {fp = fopen(fullname, "a");}
    else    {fp = fopen(fullname, "w");}
    fprintf(fp, "%6.6f \n", val);
    fclose(fp);
    free(fullname);
}

// >>>>> Interpolate to get boundary condition at current time step
double interp_bc(double *tVec, double *value, double t_current, int n_dat)
{
    int ind;
    size_t n;
    double out_val;
    if (t_current == 0.0)
    {out_val = value[0];}
    else
    {
        ind = 1;
        if (t_current > tVec[ind])
        {
            ind += 1;
            while (t_current > tVec[ind])
            {
                ind += 1;
                if (ind >= n_dat) {break;}
            }
        }
        out_val = value[ind-1] + (value[ind]-value[ind-1]) *
            (t_current-tVec[ind-1]) / (tVec[ind]-tVec[ind-1]);
    }
    return out_val;
}

// >>>>> Print at the root rank <<<<<
void mpi_print(char pstr[], int irank)
{if (irank == 0)     {printf("%s \n",pstr);}}

// >>>>> Get minimum of an array <<<<<
double getMin(double *arr, int n)
{
    int ii;
    double vmin = arr[0];
    for (ii = 1; ii < n; ii++)
    {if (arr[ii] < vmin) {vmin = arr[ii];}}
    return vmin;
}

// >>>>> Get maximum of an array <<<<<
double getMax(double *arr, int n)
{
    int ii;
    double vmax = arr[0];
    for (ii = 1; ii < n; ii++)
    {if (arr[ii] > vmax) {vmax = arr[ii];}}
    return vmax;
}

// >>>>> Check if a file exists <<<<<
int exist(char *fname)
{
    FILE *fid;
    if ((fid = fopen(fname, "r")))
    {fclose(fid);   return 1;}
    return 0;
}
