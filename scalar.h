#include "configuration.h"
#include "bathymetry.h"
#include "map.h"
#include "subgrid.h"

void scalarTransport(Data **data, Maps *map, Bath *bath, BC *bc, Sub *sub, Config *setting, int tt, int irank, int nrank);
void scalarMass(Data **data, Config *setting);
void scalarAdvection(Data **data, Maps *map, Config *setting);
void scalarDiffusion(Data **data, Maps *map, Config *setting);
void updateScalar(Data **data, BC *bc, Bath *bath, Maps *map, Config *setting, int tt, int irank, int nrank);
void enforceScalarBC(Data **data, BC *bc, Maps *map, Config *setting, int tt, int irank, int nrank);
void scalarSubStepping(Data **data, Bath *bath, BC *bc, Maps *map, Config *setting, int tt, int irank, int nrank);
