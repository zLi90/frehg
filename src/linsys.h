// Header file for linsys.c
#include"configuration.h"

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


#ifndef LINSYS_H
#define LINSYS_H



typedef struct Lisy
{
    double *val, *prec, *diag, *ltdi, *utdi, *ltri, *utri, *rhs, *x;
    int *col, *row_ptr, *nrow, *nnz, *nbands;
    double *p, *p0, *q, *z, *resi;
}Lisy;

#endif

void cgsolve(Lisy **sysm, QMatrix A, double *rhs, double *out, int nrow);
void decompose(Lisy **sysm, double omega);
void init_crs_matrix(Data *data, Lisy **sysm, Config *param, int nrow, int domain);
void build_crs_matrix(Data *data, Lisy **sysm, Config *param, int domain);
