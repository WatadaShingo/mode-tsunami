#ifndef minmax_h
#define minmax_h

#define MAX(A,B)	((A) > (B) ? (A) : (B) )
#define MIN(A,B)	((A) < (B) ? (A) : (B) )
/* ABS doesnt' work sometimes why? 
#define ABS(A)    ((A) < 0.0 ? (-(A)) : (A)) _ABS() macro is defied in <math.h> */
#define IABS(A)   ((A) < 0 ? (-(A)) : (A))
#define DELTA(A,B) ((A) == (B) ? 1.0 : 0.0 )
#define IDELTA(A,B) ((A) == (B) ? 1 : 0 )
#define OUT(x,a,b) ((x) < (a) || (x) > (b) )

int set_imax();
int set_imin();
float set_max();
float set_min();
double set_dmax();
double set_dmin();
void swapi();
void swapf();
void sort_f();
void sort_d();
int search_nearf();
int search_neard();
float get_fmean();
double get_dmean();

#endif

