// Header file for solve.c
#include "configuration.h"
#include "initialize.h"
#include "map.h"

void solve(Data **data, Map *smap, Map *gmap, Config *param, int irank, int nrank);
void get_current_bc(Data **data, Config *param, double t_current);
void get_evaprain(Data **data, Map *gmap, Config *param, double t_current);
void print_end_info(Data **data, Map *smap, Map *gmap, Config *param, int irank);
