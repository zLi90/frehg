// Functions used to build and solve the linear system
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<omp.h>

#include "configuration.h"
#include "initialize.h"
#include "linsys.h"
#include "map.h"
#include "utility.h"

#include "laspack/errhandl.h"
#include "laspack/vector.h"
#include "laspack/matrix.h"
#include "laspack/qmatrix.h"
#include "laspack/operats.h"
#include "laspack/factor.h"
#include "laspack/precond.h"
#include "laspack/eigenval.h"
#include "laspack/rtc.h"
#include "laspack/itersolv.h"
#include "laspack/mlsolv.h"

void cgsolve(Lisy **sysm, QMatrix A, double *rhs, double *out, int nrow);
void decompose(Lisy **sysm);
void init_crs_matrix(Data *data, Lisy **sysm, Config *param, int nrow, int domain);
void build_crs_matrix(Data *data, Lisy **sysm, Config *param, int domain);


void cgsolve(Lisy **sysm, QMatrix A, double *rhs, double *out, int nrow)
{
    double rho, rho0, alpha, beta, deno, eps = 1.0, eps0 = 1e-10;
    int ii, id, nd, iter = 1, iter_max = 10000000;

    double *val = (*sysm)->val;
    int *col = (*sysm)->col;
    int *row_ptr = (*sysm)->row_ptr;

    // create preconditioner
    Vector prec;
    V_Constr(&prec, "prec", nrow, Normal, True);
    Vector rvec;
    V_Constr(&rvec, "rvec", nrow, Normal, True);

    // initial guess
    srand(time(NULL));
    for (ii = 0; ii < nrow; ii++)   {
        (*sysm)->x[ii] = rand() % 10;
        (*sysm)->p[ii] = 0.0;
    }
    // get initial residual
    for (ii = 0; ii < nrow; ii++)   {
        (*sysm)->resi[ii] = 0.0;
        for (int jj = row_ptr[ii]; jj < row_ptr[ii+1]; jj++)    {
            (*sysm)->resi[ii] += val[jj] * (*sysm)->x[col[jj]];
        }
        (*sysm)->resi[ii] = rhs[ii] - (*sysm)->resi[ii];
    }

    // iterate to update solution
    while (eps > eps0 & iter < iter_max)    {
        // preconditioner
        // #pragma omp parallel for
        for (ii = 0; ii < nrow; ii++)   {V_SetCmp(&rvec, ii+1, (*sysm)->resi[ii]);}
        SSORPrecond(&A, &prec, &rvec, 1.0);
        // #pragma omp parallel for
        for (ii = 0; ii < nrow; ii++)    {(*sysm)->z[ii] = V_GetCmp(&prec, ii+1);}

        rho = 0.0;
        deno = 0.0;
        // rho = r * z
        #pragma omp parallel for reduction (+:rho)
        for (ii = 0; ii < nrow; ii++)   {
            rho += (*sysm)->resi[ii] * (*sysm)->z[ii];
        }
        // p = beta*p + z
        if (iter == 1)  {beta = 0.0;}
        else    {beta = rho / rho0;}
        rho0 = rho;
        #pragma omp parallel for
        for (ii = 0; ii < nrow; ii++)   {(*sysm)->p[ii] = beta * (*sysm)->p[ii] + (*sysm)->z[ii];}
        // q = Ap
        #pragma omp parallel for reduction (+:deno)
        for (ii = 0; ii < nrow; ii++)   {
            (*sysm)->q[ii] = 0.0;
            for (int jj = row_ptr[ii]; jj < row_ptr[ii+1]; jj++)    {
                (*sysm)->q[ii] += val[jj] * (*sysm)->p[col[jj]];
            }
            deno += (*sysm)->p[ii] * (*sysm)->q[ii];
        }
        // alpha = rho / pq
        alpha = rho / deno;
        #pragma omp parallel
        {
            // x = x + alpha * p
            #pragma omp for nowait
            for (ii = 0; ii < nrow; ii++)   {(*sysm)->x[ii] += alpha * (*sysm)->p[ii];}
            // r = r - alpha * q
            #pragma omp for nowait
            for (ii = 0; ii < nrow; ii++)   {(*sysm)->resi[ii] -= alpha * (*sysm)->q[ii];}
        }
        eps = getMax((*sysm)->resi, nrow);
        // printf("   --pcgsolve--  iter %d is completed! eps=%f\n",iter,eps);
        iter += 1;
    }
    #pragma omp parallel for
    for (ii = 0; ii < nrow; ii++)   {out[ii] = (*sysm)->x[ii];}
    printf("   --pcgsolve--  iter %d is completed!\n",iter);
    V_Destr(&rvec);
    V_Destr(&prec);
}


void decompose(Lisy **sysm)
{
    int ii, jj;
    int nrow = (*sysm)->nrow[0];
    int nnz = (*sysm)->nnz[0];
    int *col = (*sysm)->col;
    for (ii = 0; ii < nnz; ii++)    {
        (*sysm)->diag[ii] = 0.0;
        (*sysm)->ltri[ii] = 0.0;
        (*sysm)->utri[ii] = 0.0;
    }
    for (ii = 0; ii < nrow; ii++)   {
        for (jj = (*sysm)->row_ptr[ii]; jj < (*sysm)->row_ptr[ii+1]; jj++)    {
            if (col[jj] == ii)  {(*sysm)->diag[jj] = (*sysm)->val[jj];}
            else if (col[jj] > ii)  {(*sysm)->utri[jj] = (*sysm)->val[jj];}
            else    {(*sysm)->ltri[jj] = (*sysm)->val[jj];}
        }
    }
}


void build_crs_matrix(Data *data, Lisy **sysm, Config *param, int domain)
{
    int irow, kk, col_ind, nrow = (*sysm)->nrow[0];
    kk = 0;
    if ((*sysm)->nbands[0] == 7)   {
        for (irow = 0; irow < (*sysm)->nrow[0]; irow++)    {
            col_ind = irow - param->nx*param->nz;
            if (col_ind >= 0)    {(*sysm)->val[kk] = data->Gym[irow];    kk+=1;}
            col_ind = irow - param->nz;
            if (col_ind >= 0)    {(*sysm)->val[kk] = data->Gxm[irow];    kk+=1;}
            col_ind = irow - 1;
            if (col_ind >= 0)    {(*sysm)->val[kk] = data->Gzm[irow];    kk+=1;}
            (*sysm)->val[kk] = data->Gct[irow];    kk+=1;
            col_ind = irow + 1;
            if (col_ind < nrow)    {(*sysm)->val[kk] = data->Gzp[irow];    kk+=1;}
            col_ind = irow + param->nz;
            if (col_ind < nrow)    {(*sysm)->val[kk] = data->Gxp[irow];    kk+=1;}
            col_ind = irow + param->nx*param->nz;
            if (col_ind < nrow)    {(*sysm)->val[kk] = data->Gyp[irow];    kk+=1;}
        }
    }
    else if ((*sysm)->nbands[0] == 5)   {
        if (domain == 0)    {
            for (irow = 0; irow < (*sysm)->nrow[0]; irow++)    {
                col_ind = irow - param->nx;
                if (col_ind >= 0)    {(*sysm)->val[kk] = data->Sym[irow];    kk+=1;}
                col_ind = irow - 1;
                if (col_ind >= 0)    {(*sysm)->val[kk] = data->Sxm[irow];    kk+=1;}
                (*sysm)->val[kk] = data->Sct[irow];    kk+=1;
                col_ind = irow + 1;
                if (col_ind < nrow)    {(*sysm)->val[kk] = data->Sxp[irow];    kk+=1;}
                col_ind = irow + param->nx;
                if (col_ind < nrow)    {(*sysm)->val[kk] = data->Syp[irow];    kk+=1;}
            }
        }
        else {
            for (irow = 0; irow < (*sysm)->nrow[0]; irow++)    {
                col_ind = irow - param->nz;
                if (col_ind >= 0)    {(*sysm)->val[kk] = data->Gym[irow];    kk+=1;}
                col_ind = irow - 1;
                if (col_ind >= 0)    {(*sysm)->val[kk] = data->Gzm[irow];    kk+=1;}
                (*sysm)->val[kk] = data->Gct[irow];    kk+=1;
                col_ind = irow + 1;
                if (col_ind < nrow)    {(*sysm)->val[kk] = data->Gzp[irow];    kk+=1;}
                col_ind = irow + param->nz;
                if (col_ind < nrow)    {(*sysm)->val[kk] = data->Gyp[irow];    kk+=1;}
            }
        }
    }
    else if ((*sysm)->nbands[0] == 3)   {
        if (domain == 0)    {
            for (irow = 0; irow < (*sysm)->nrow[0]; irow++)    {
                col_ind = irow - 1;
                if (col_ind >= 0)    {(*sysm)->val[kk] = data->Sym[irow];    kk+=1;}
                (*sysm)->val[kk] = data->Sct[irow];    kk+=1;
                col_ind = irow + 1;
                if (col_ind < nrow)    {(*sysm)->val[kk] = data->Syp[irow];    kk+=1;}
            }
        }
        else {
            for (irow = 0; irow < (*sysm)->nrow[0]; irow++)    {
                col_ind = irow - 1;
                if (col_ind >= 0)    {(*sysm)->val[kk] = data->Gzm[irow];    kk+=1;}
                (*sysm)->val[kk] = data->Gct[irow];    kk+=1;
                col_ind = irow + 1;
                if (col_ind < nrow)    {(*sysm)->val[kk] = data->Gzp[irow];    kk+=1;}
            }
        }
    }

}


void init_crs_matrix(Data *data, Lisy **sysm, Config *param, int nrow, int domain)
{
    int nnz, ii, jj, kk, nbands, col_ind, irow, ncol, ncol_old;
    *sysm = malloc(sizeof(Lisy));
    // domain = 0 means shallowwater, domain = 1 means groundwater
    if (domain == 0)    {
        // if 1D, nx must be 1
        if (param->nx == 1) {nnz = nrow + 2.0*(nrow-1); nbands = 3;}
        else    {nnz = nrow + 2.0*(nrow-1) + 2.0*(nrow-param->nx);  nbands = 5;}
    }
    else if (domain == 1)   {
        if (param->nx == 1)    {
            // if 1D for subsurface, must be nz > 1
            if (param->ny == 1) {nnz = nrow + 2.0*(nrow-1); nbands = 3;}
            // if 2D, must be ny-nz
            else    {nnz = nrow + 2.0*(nrow-1) + 2.0*(nrow-param->nz);  nbands = 5;}
        }
        // if nx != 1, must be 3D
        else    {nnz = nrow + 2.0*(nrow-1) + 2.0*(nrow-param->nz) + 2.0*(nrow-param->nz*param->nx); nbands = 7;}
    }
    (*sysm)->val = malloc(nnz*sizeof(double));
    (*sysm)->prec = malloc(nnz*sizeof(double));
    (*sysm)->diag = malloc(nnz*sizeof(double));
    (*sysm)->ltdi = malloc(nnz*sizeof(double));
    (*sysm)->utdi = malloc(nnz*sizeof(double));
    (*sysm)->ltri = malloc(nnz*sizeof(double));
    (*sysm)->utri = malloc(nnz*sizeof(double));
    (*sysm)->col = malloc(nnz*sizeof(int));
    (*sysm)->row_ptr = malloc((nrow+1)*sizeof(int));
    (*sysm)->rhs = malloc(nrow*sizeof(double));
    (*sysm)->x = malloc(nrow*sizeof(double));
    (*sysm)->z = malloc(nrow*sizeof(double));
    (*sysm)->p = malloc(nrow*sizeof(double));
    (*sysm)->q = malloc(nrow*sizeof(double));
    (*sysm)->resi = malloc(nrow*sizeof(double));
    (*sysm)->nrow = malloc(1*sizeof(int));
    (*sysm)->nnz = malloc(1*sizeof(int));
    (*sysm)->nbands = malloc(1*sizeof(int));
    (*sysm)->nrow[0] = nrow;
    (*sysm)->nnz[0] = nnz;
    (*sysm)->nbands[0] = nbands;
    for (ii = 0; ii < nnz; ii++)    {(*sysm)->val[ii] = 0.0;}
    for (ii = 0; ii < nrow; ii++)    {(*sysm)->x[ii] = 0.0;    (*sysm)->rhs[ii] = 0.0;}
    // set up row and col index
    kk = 0;
    if (nbands == 7)    {
        ncol_old = 0;
        for (irow = 0; irow < nrow; irow++) {
            ncol = 0;
            col_ind = irow - param->nx*param->nz;
            if (col_ind >= 0)    {(*sysm)->col[kk] = col_ind;    kk+=1;  ncol+=1;}
            col_ind = irow - param->nz;
            if (col_ind >= 0)    {(*sysm)->col[kk] = col_ind;    kk+=1;  ncol+=1;}
            col_ind = irow - 1;
            if (col_ind >= 0)    {(*sysm)->col[kk] = col_ind;    kk+=1;  ncol+=1;}
            (*sysm)->col[kk] = irow;    kk+=1;  ncol+=1;
            col_ind = irow + 1;
            if (col_ind < nrow)    {(*sysm)->col[kk] = col_ind;    kk+=1;  ncol+=1;}
            col_ind = irow + param->nz;
            if (col_ind < nrow)    {(*sysm)->col[kk] = col_ind;    kk+=1;  ncol+=1;}
            col_ind = irow + param->nx*param->nz;
            if (col_ind < nrow)    {(*sysm)->col[kk] = col_ind;    kk+=1;  ncol+=1;}
            // calculate row_ptr
            if (irow == 0)    {(*sysm)->row_ptr[irow] = 0;}
            else    {(*sysm)->row_ptr[irow] = (*sysm)->row_ptr[irow-1] + ncol_old;}
            if (irow == nrow-1)    {(*sysm)->row_ptr[irow+1] = (*sysm)->row_ptr[irow]+ncol;}
            ncol_old = ncol;
        }
    }
    else if (nbands == 5)   {
        ncol_old = 0;
        for (irow = 0; irow < nrow; irow++) {
            ncol = 0;
            if (domain == 0)    {col_ind = irow - param->nx;}
            else    {col_ind = irow - param->nz;}
            if (col_ind >= 0)    {(*sysm)->col[kk] = col_ind;    kk+=1;  ncol+=1;}
            col_ind = irow - 1;
            if (col_ind >= 0)    {(*sysm)->col[kk] = col_ind;    kk+=1;  ncol+=1;}
            (*sysm)->col[kk] = irow;    kk+=1;  ncol+=1;
            col_ind = irow + 1;
            if (col_ind < nrow)    {(*sysm)->col[kk] = col_ind;    kk+=1;  ncol+=1;}
            if (domain == 0)    {col_ind = irow + param->nx;}
            else    {col_ind = irow + param->nz;}
            if (col_ind < nrow)    {(*sysm)->col[kk] = col_ind;    kk+=1;  ncol+=1;}
            // calculate row_ptr
            if (irow == 0)    {(*sysm)->row_ptr[irow] = 0;}
            else    {(*sysm)->row_ptr[irow] = (*sysm)->row_ptr[irow-1] + ncol_old;}
            if (irow == nrow-1)    {(*sysm)->row_ptr[irow+1] = (*sysm)->row_ptr[irow]+ncol;}
            ncol_old = ncol;
        }
    }
    else if (nbands == 3)   {
        ncol_old = 0;
        for (irow = 0; irow < nrow; irow++) {
            ncol = 0;
            col_ind = irow - 1;
            if (col_ind >= 0)    {(*sysm)->col[kk] = col_ind;    kk+=1;  ncol+=1;}
            (*sysm)->col[kk] = irow;    kk+=1;  ncol+=1;
            col_ind = irow + 1;
            if (col_ind < nrow)    {(*sysm)->col[kk] = col_ind;    kk+=1;  ncol+=1;}
            // calculate row_ptr
            if (irow == 0)    {(*sysm)->row_ptr[irow] = 0;}
            else    {(*sysm)->row_ptr[irow] = (*sysm)->row_ptr[irow-1] + ncol_old;}
            if (irow == nrow-1)    {(*sysm)->row_ptr[irow+1] = (*sysm)->row_ptr[irow]+ncol;}
            ncol_old = ncol;
        }
    }
}
