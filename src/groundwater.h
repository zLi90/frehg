// Header file for groundwater.c
#include "configuration.h"
#include "initialize.h"
#include "map.h"

// #include "../laspack/vector.h"
// #include "../laspack/qmatrix.h"
// #include "../laspack/rtc.h"

#include "../laspack/errhandl.h"
#include "../laspack/vector.h"
#include "../laspack/matrix.h"
#include "../laspack/qmatrix.h"
#include "../laspack/operats.h"
#include "../laspack/factor.h"
#include "../laspack/precond.h"
#include "../laspack/eigenval.h"
#include "../laspack/rtc.h"
#include "../laspack/itersolv.h"
#include "../laspack/mlsolv.h"

void solve_groundwater(Data **data, Map *smap, Map *gmap, Config *param, int irank, int nrank);
void compute_K_face(Data **data, Map *gmap, Config *param, int irank, int nrank);
void baroclinic_face(Data **data, Map *gmap, Config *param, int irank, int nrank);
void groundwater_mat_coeff(Data **data, Map *gmap, Config *param, int irank);
void groundwater_rhs(Data **data, Map *gmap, Config *param, int irank);
void build_groundwater_system(Data *data, Map *gmap, Config *param, QMatrix A, Vector b);
void solve_groundwater_system(Data **data, Map *gmap, QMatrix A, Vector b, Vector x, Config *param);
void enforce_head_bc(Data **data, Map *gmap, Config *param);
void groundwater_flux(Data **data, Map *gmap, Config *param, int irank);
void check_room(Data **data, Map *gmap, Config *param);
void update_water_content(Data **data, Map *gmap, Config *param);
void enforce_moisture_bc(Data **data, Map *gmap, Config *param);
void reallocate_water_content(Data **data, Map *gmap, Config *param, int irank);
int check_adj_sat(Data *data, Map *gmap, Config *param, int ii);
void check_head_gradient(Data **data, Map *gmap, Config *param, int ii);
double allocate_send(Data **data, Map *gmap, Config *param, int ii, double dV);
double allocate_recv(Data **data, Map *gmap, Config *param, int ii, double dV);
void volume_by_flux_subs(Data **data, Map *gmap, Config *param);
void adaptive_time_step(Data *data, Map *gmap, Config **param, int root, int irank);
