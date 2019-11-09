
#include "configuration.h"

#ifndef BATHYMETRY_H
#define BATHYMETRY_H

typedef struct Bath
{
  double *allZ, *allEdgeX, *allEdgeY;
  double *bottomZ, *edgeX, *edgeY;
  double *bottomXP;
  double *bottomYP;
  double *offset;
}Bath;

#endif

void ReadBathymetry(Bath **bath, Config *setting, int irank);
void InitBathymetry(Bath *bath, Config *setting, int irank);
