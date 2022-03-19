// Header file for shallowwater.c
#include "configuration.h"
#include "initialize.h"
#include "map.h"


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

void solve_shallowwater(Data **data, Map *smap, Map *gmap, Config *param, int irank, int nrank);
void shallowwater_velocity(Data **data, Map *smap, Map *gmap, Config *param, int irank, int nrank);
void momentum_source(Data **data, Map *smap, Config *param);
void wind_source(Data **data, Map *smap, Config *param, int ii);
void shallowwater_rhs(Data **data, Map *smap, Config *param);
void shallowwater_mat_coeff(Data **data, Map *smap, Config *param, int irank, int nrank);
void build_shallowwater_system(Data *data, Map *smap, Config *param, QMatrix A, Vector b);
void solve_shallowwater_system(Data **data, Map *smap, QMatrix A, Vector b, Vector x, Config *param);
void enforce_surf_bc(Data **data, Map *smap, Config *param, int irank, int nrank);
void cfl_limiter(Data **data, Map *smap, Config *param);
void evaprain(Data **data, Map *smap, Config *param);
void subsurface_source(Data **data, Map *smap, Config *param);
void waterfall_location(Data **data, Map *smap, Config *param);
void update_velocity(Data **data, Map *smap, Config *param, int irank);
void waterfall_velocity(Data **data, Map *smap, Config *param);
void enforce_velo_bc(Data **data, Map *smap, Config *param, int irank, int nrank);
void interp_velocity(Data **data, Map *smap, Config *param);
void update_drag_coef(Data **data, Config *param);
void update_subgrid_variable(Data **data, Map *smap, Config *param, int irank, int nrank);
void subgrid_index(Data **data, Map *smap, Config *param);
void subgrid_interp_and_combine(Data **data, Map *smap, Config *param, int irank, int nrank);
void volume_by_flux(Data **data, Map *smap, Config *param);
