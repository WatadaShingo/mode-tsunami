#ifndef mallocutil_h
#define mallocutil_h

#include "clib.h"

short int *malloc_vshortint();
short int **malloc_mshortint();
int *malloc_vint();
int **malloc_mint();
float *malloc_vfloat();
float **malloc_mfloat();
double *malloc_vdouble();
double **malloc_mdouble();
double ***malloc_m3double();
complex *malloc_vcomplex();
complex **malloc_mcomplex();
dcomplex *malloc_vdcomplex();
dcomplex **malloc_mdcomplex();
dcomplex ***malloc_m3dcomplex();

void free_vshortint();
void free_vint();
void free_mshortint();
void free_mint();
void free_vfloat();
void free_mfloat();
void free_vdouble();
void free_mdouble();
void free_m3double();
void free_vcomplex();
void free_mcomplex();
void free_vdcomplex();
void free_mdcomplex();
void free_m3dcomplex();

#endif /* mallocutil_h */
