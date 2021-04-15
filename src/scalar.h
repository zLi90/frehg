// Header file for scalar.c
#include"initialize.h"
#include"map.h"
#include"configuration.h"

#ifndef SCALAR_H
#define SCALAR_H

#endif

void scalar_shallowwater(Data **data, Map *smap, Config *param, int irank, int nrank, int kk);
void scalar_groundwater(Data **data, Map *gmap, Config *param, int irank, int nrank, int kk);
void enforce_scalar_bc(Data **data, Map *gmap, Config *param, int kk, int irank);
void update_rhovisc(Data **data, Map *gmap, Config *param, int irank);
void dispersion_tensor(Data **data, Map *gmap, Config *param);
