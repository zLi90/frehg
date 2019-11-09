// Utility functions for SubFREHD-C, ZhiLi 20170502


double getMax(double *arr, int N);
double getMin(double *arr, int N);
void dataSplit(double *t, double *value, double *arr, int N);
void dataInterp(double *y, double *t, double *value, int ind, int Nt, double dt, double tNStart);
time_t dateNum(char *t);
int getTime(char *str, int N, int n);
struct tm initTime(struct tm st);
int searchInd(double *arr, double targ, double itval, int N);
void writeText(char text[100], int irank);
double interpSurf(double x, double y, double x1, double y1, double x2, double y2);
void interpRestart(double *x, double *y, int ind1, int ind2, int r, int NX);
