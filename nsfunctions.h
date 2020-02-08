

#include "configuration.h"
#include "bathymetry.h"
#include "map.h"

#include "laspack/vector.h"
#include "laspack/qmatrix.h"

void updateAllFaceElev(Bath *bath, Maps *map, Config *setting);
void updateAllFaceDepth(Data **data, Bath *bath, Maps *map, Config *setting, int irank, int nrank);
void updateDepth(Data **data, Maps *map, Bath *bath, Config *setting);
void updateVelocity(Data **data, Maps *map, Config *setting);
void enforceVelocityBC(Data **data, Maps *map, Config *setting);
void enforceFreeSurfaceBC(Data **data, Maps *map, Config *setting);
void velocityInterp(Data **data, Config *setting);
void computeFaceFlowRate(Data **data, Maps *map, Config *setting);
void computeVolume(Data **data, Config *setting);
void volumeByFlux(Data **data, Maps *map, BC *bc, Config *setting, int tt, int irank);
void enforceTidalBC(Data **data, Bath *bath, BC *bc, Config *setting, int tt, int irank, int nrank);
void zeroNegSurf(Data **data, Bath *bath, Config *setting, int tt);
void limiterCFL(Data **data, Bath *bath, Maps *map, Config *setting);
void explicitVelocity(Data **data, Maps *map, Config *setting);
void advectionTerm(Data **data, Maps *map, Config *setting);
void diffusionTerm(Data **data, Maps *map, Config *setting);
void windTerm(Data **data, BC *bc, Config *setting, int tt);
void dragTerm(Data **data, Maps *map, Config *setting);
void thinLayerDrag(Data **data, Config *setting);
void dragInversion(Data **data, Maps *map, Config *setting);
void rainfall(Data **data, BC *bc, int tt, Config *setting);
void matrixSourceTerm(Data **data, Maps *map, BC *bc, Config *setting, int tt, int irank);
void findInflowLocation(Data **data, Maps *map, Config *setting);
void matrixCoeff(Data **data, Maps *map, Config *setting);
void setupMatrix(Data *data, Maps *map, Config *setting, QMatrix A);
void adjustMatrixBoundary(Data **data, BC *bc, Maps *map, Config *setting, int irank, int nrank);
void adjustTidalBoundary(Data **data, BC *bc, Bath *bath, Config *setting, int tt, int irank, int nrank);
void solveMatrix(Data *data, Config *setting, QMatrix A, Vector x, Vector z);
void getFreeSurface(Data **data, Maps *map, Config *setting, Vector x);
void adjustTidalVelocity(Data **data, Maps *map, Config *setting, int irank, int nrank);
void detectWaterfallLocation(Data **data, Bath *bath, Maps *map, Config *setting);
void waterfallCorrection(Data **data, Bath *bath, Maps *map, Config *setting);
void updateCD(Data **data, Config *setting);
void monitorCFL(Data **data, Bath *bath, int irank, Config *setting, int tt, int root);
void computeEvapRain(Data **data, Bath *bath, Maps *map, BC *bc, Config *setting, int tt);
void infiltration(Data **data, Bath *bath, Maps *map, Config *setting);
