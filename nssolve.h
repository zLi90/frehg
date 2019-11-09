

#include "configuration.h"
#include "bathymetry.h"
#include "map.h"
#include "subgrid.h"

#include "laspack/vector.h"
#include "laspack/qmatrix.h"

void SolveAll(Data **data, Sub **sub, Ground **ground, Maps *map, Gmaps *gmap, BC *bc, Bath *bath, Config *setting, int irank, int nrank);
void oneCompleteStep(Data **data, Sub **sub, Maps *map, BC *bc, Bath *bath, Config *setting, int tt, \
  QMatrix A, Vector x, Vector z, int itank, int nrank);
void solveSourceTerm(Data **data, Sub *sub, Maps *map, BC *bc, Bath *bath, Config *setting, int tt, \
  int irank, int nrank);
void solveFreeSurface(Data **data, Sub *sub, Maps *map, BC *bc, Bath *bath, Config *setting, int tt, \
  QMatrix A, Vector x, Vector z, int irank, int nrank);
void updateData(Data **data, Sub **sub, Maps *map, BC *bc, Bath *bath, Config *setting, int tt, \
  int irank, int nrank);
