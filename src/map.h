// Header file for map.c


#ifndef MAP_H
#define MAP_H


typedef struct Map
{
    // surface maps
    int *cntr, *iPjc, *iMjc, *icjP, *icjM, *ii, *jj;
    int *iPjP, *iPjM, *iMjP, *iMjM;
    int *iPin, *iPou, *iMin, *iMou, *jPin, *jPou, *jMin, *jMou;
    // subsurface maps
    int *iPjckc, *iMjckc, *icjPkc, *icjMkc, *icjckP, *icjckM, *kk;
    int *actv, *actvXp, *actvYp, *istop, *top2d, *kPin, *kPou, *kMin, *kMou;
    int nactv;
    double *bot1d, *bot3d, *dz3d;
    double *zcntr, *zcntr_root, *zcntr_out;
}Map;

#endif


void build_surf_map(Map **map, Config *param);
void build_subsurf_map(Map **map, Map *smap, double *bath, double *offset, Config *param, int irank);
