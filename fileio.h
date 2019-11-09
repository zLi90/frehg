

#include "configuration.h"
#include "initialize.h"

void ReadFile(double *arr, char filename[], int N);
void ReadFileInt(int *arr, char filename[], int N);
void DataOutput(Data **data, Ground **ground, Bath *bath, Config *setting, int tt, int root, int irank);
void SubOutput(Sub **sub, Config *setting, int tt, int root, int irank);
void WriteFile(double *ally, char *filename, Config *setting, int tt, int N);
void combineSerial(double *ally, double *y, Config *setting, int N);
void writeRestartFile(Data **data, Bath *bath, Config *setting, int tt, int root, int irank);
