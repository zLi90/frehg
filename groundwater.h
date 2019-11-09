#include "configuration.h"
#include "bathymetry.h"
#include "map.h"

#include "laspack/vector.h"
#include "laspack/qmatrix.h"

void groundwaterExchange(Data **data, Ground **ground, Maps *map, Gmaps *gmap, Config *setting, int irank, int nrank);
void computeConductance(Ground **ground, Gmaps *gmap, Config *setting, int irank, int nrank);
void groundMatrixCoeff(Ground **ground, Data **data, Gmaps *gmap, Config *setting);
void setupGroundMatrix(Ground *ground, Gmaps *gmap, Config *setting, QMatrix A);
void solveGroundMatrix(Ground *ground, Gmaps *gmap, Config *setting, QMatrix A, Vector x, Vector z);
void getHead(Ground **ground, Maps *map, Config *setting, Vector x);
void enforceSidewallBC(Ground **ground, Gmaps *gmap, Config *setting);
void computeSeepage(Data **data, Ground *ground, Maps *map, Gmaps *gmap, Config *setting);
void computeFlowRate(Ground **ground, Data *data, Gmaps *gmap, Config *setting);
void groundScalarAdvection(Ground **ground, Data *data, Gmaps *gmap, Config *setting);
void groundScalarDiffusion(Ground **ground, Data *data, Gmaps *gmap, Config *setting);
void groundScalarTransport(Ground **ground, Data **data, Gmaps *gmap, Maps *map, Config *setting, int irank, int nrank);
double updateSaturation(double h, Config *setting);
double updateDSDH(double h, Config *setting);
double relativePerm(double h, double hP, Config *setting);
