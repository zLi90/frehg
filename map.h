#include "bathymetry.h"

#ifndef MAP_H
#define MAP_H

typedef struct Maps
{
  int *cntr, *trps, *sprt, *iPjc, *iMjc, *icjP, *icjM, *ii2d, *jj2d;
  int *iPjP, *iPjM, *iMjP, *iMjM;
  int *iPbd, *iPgt, *iMbd, *iMgt, *jPbd, *jPgt, *jMbd, *jMgt;
  int *iPPjc, *iMMjc, *icjPP, *icjMM;
}Maps;

typedef struct Gmaps
{
    int *cntr, *ii, *jj, *kk, maxLay, ngtiP, ngtiM, ngtjP, ngtjM, *actv;
    int *nlay, *istop, *top2D;
    int *iPjckc, *iMjckc, *icjPkc, *icjMkc, *icjckP, *icjckM;
    int *iPbd, *iPgt, *iMbd, *iMgt, *jPbd, *jPgt, *jMbd, *jMgt;
    int *kPbd, *kPgt, *kMbd, *kMgt;
    double *htop, *dz3d, *bot3d;
}Gmaps;

#endif

void createMaps(Maps **map, Config *setting);
void createGmaps(Gmaps **gmap, Bath *bath, Maps *map, Config *setting);
