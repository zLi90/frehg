// Utility functions for SubFREHD-C, ZhiLi 20170502
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

// -----------------------------------------------------------------------------
// This file contains auxiliary functions that are not specifically designed
// for FREHD. They are mostly used for data processing purposes.
// - Zhi Li 2017-05-04 -
// -----------------------------------------------------------------------------

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

// =============== Get the maximum and minimum of an array ===============
double getMax(double *arr, int N)
{
  int ii;
  double maxvalue = arr[0];
  for (ii = 1; ii < N; ii++)
  {
    if (arr[ii] > maxvalue)
    {maxvalue = arr[ii];}
  }
  return maxvalue;
}
double getMin(double *arr, int N)
{
  int ii;
  double minvalue = arr[0];
  for (ii = 1; ii < N; ii++)
  {
    if (arr[ii] < minvalue)
    {minvalue = arr[ii];}
  }
  return minvalue;
}

// =============== Split the BCs into time and value columns ===============
void dataSplit(double *t, double *value, double *arr, int N)
{
  int ii;
  for (ii = 0; ii < N; ii++)
  {
    t[ii] = arr[2*ii];
    value[ii] = arr[2*ii+1];
  }
}

// ===== Interpolate 1d array into an array with required length =====
// x = original array, with lenth Nx
// y = output array, with length Ny
// note: For BC, usually Ny > Nx
void dataInterp(double *y, double *t, double *value, int ind, int Nt, double dt, double tNStart)
{
  int ii, kk = ind;
  double slope, intercept, t0 = tNStart, tc;
  for (ii = 0; ii < Nt; ii++)
  {
    // compute the date number of each step
    tc = t0 + ii*dt;
    if (ii > 0 & tc > t[kk+1])
    {
      // marching to the next bc interval if necessary
      kk = kk + 1;
    }
    // create the linear relation within a bc interval
    slope = (value[kk+1]-value[kk]) / (t[kk+1]-t[kk]);
    intercept = value[kk] - slope * t[kk];
    // Interpolate
    y[ii] = slope * tc + intercept;
  }
//4 - 3 5 7 9
//9 - 3 3.75 4.5 5.25 6 6.75 7.5 8.25 9
}

// =============== Initialize the time struct ===============
struct tm initTime(struct tm st)
{
  st.tm_sec = 0; st.tm_min = 0; st.tm_hour = 0;
  st.tm_mday = 0; st.tm_mon = 0; st.tm_year = 0; st.tm_isdst = 0;
  return(st);
}

// =============== Get date number from date string ===============
time_t dateNum(char *t)
{
  struct tm dateStruct;
  dateStruct = initTime(dateStruct);
  int year = getTime(t, 4, 0);
  int mon = getTime(t, 2, 5);
  int day = getTime(t, 2, 8);
  dateStruct.tm_year = year - 1900;
  dateStruct.tm_mon = mon-1;
  dateStruct.tm_mday = day;
  time_t tNum = mktime(&dateStruct);
  return(tNum);
}

// =============== Get field value from date string ==============
// N = length of the string, N=4 for year, N=2 for month and day
// n = starting index to read in the string
int getTime(char *str, int N, int n)
{
  int y = 0, ii;
  char val[N];
  for(ii = 0; ii < N; ii++)
   { val[ii] = str[ii+n]; }
  y = (int) atof(val);
  memset(val, 0, N);
  return(y);
}

// ====== Search the index in an array closest to the target value ======
// arr = an array whose elements will be searched
// targ = the target value to be compared
// itval = interval between elements in arr
int searchInd(double *arr, double targ, double itval, int N)
{
  int count = 0, ind = 0;
  while(1)
  {
    //printf("count,arr, targ, itval = %d,%lf,%lf,%lf\n",count,arr[count],targ,itval);
    if (fabs(arr[count]-targ) <= itval/2)
    {
      ind = count;
      break;
    }
    count += 1;
    if (count > N)
    {
      printf("WARNING: No time index was found for BC...\n");
      break;
    }
  }
  return(ind);
}

// ========== write text on screen ==========
void writeText(char text[100], int irank)
{if (irank == 0) {printf("%s \n",text);} }

// ========== interpolate the surface elevation ==========
double interpSurf(double x, double y, double x1, double y1, double x2, double y2)
{
  if (x2 == x1)
  {y = y1;}
  else
  {y = (y2-y1)*(x-x2)/(x2-x1) + y2;}
  return y;
}

// ========== interpolate coarse grid data onto fine grid ==========
void interpRestart(double *x, double *y, int ind1, int ind2, int r, int NX)
{
  int ii, jj, kk, ll = 0, col, row, base;
  for (ii = ind1; ii < ind2; ii++)
  {
    col = floor(ll / (NX / r));
    row = ll - col*(NX / r);
    base = col * NX * r;
    for (jj = 0; jj < r; jj++)
    {
      for (kk = 0; kk < r; kk++)
      {y[base + NX*kk + row*r + jj] = x[ii];}
    }
    ll++;
  }
}
