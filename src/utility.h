// Head file for utility.c
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"configuration.h"
#include"initialize.h"
#include"map.h"

double compute_wch(Data *data, int ii, Config *param);
double compute_hwc(Data *data, int ii, Config *param);
double compute_ch(Data *data, int ii, Config *param);
double compute_K(Data *data, double *Ksat, int ii, Config *param);
double compute_dKdwc(Data *data, double *Ksat, int ii, Config *param);
double compute_dKdh(Data *data, double *Ksat, int ii, Config *param);
double compute_dwcdh(Data *data, int ii, Config *param);
double tvd_superbee(double sp, double sc, double sm, double u, double delta, Config *param);
char* read_one_input(char field[], char fname[]);
double read_one_input_double(char field[], char fname[]);
int * read_one_input_array(char field[], char fname[], int n);
double * read_one_input_array_double(char field[], char fname[], int n);
void load_data(double *arr, char filename[], int n);
void load_bc(double *value, double *tVec, char filename[], int n);
int get_index_sort(double *arr, double targ, double itval, int n);
void reorder_surf(double *out, double *root, Config *param);
void reorder_subsurf(double *out, double *root, Map *gmap, Config *param);
void write_output(Data **data, Map *gmap, Config *param, int tt, int root, int irank);
void write_one_file(double *ally, char *filename, Config *param, int tt, int n);
void append_to_file(char *filename, double val, Config *param);
void root_to_rank(double *root_array, double *rank_array, Config *param, int irank, int nrank, int offind, double offset);
double interp_bc(double *tVec, double *value, double t_current, int n_dat);
double interp_sub(double *layers, double **subvar, double eta, int ii, int kk);
void mpi_print(char pstr[], int irank);
double getMin(double *arr, int n);
double getMax(double *arr, int n);
int exist(char *fname);
