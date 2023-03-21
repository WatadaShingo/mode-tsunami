#include <math.h>
#ifndef clib_h
#define clib_h

typedef struct {
	float re;
	float im;
} complex;

typedef struct {
	double re;
	double im;
} dcomplex;

typedef struct {
	float r;	/* modulus */
	float t;	/* argument in radian */
} polar;

typedef struct {
	double r;
	double t;
} dpolar;

complex comp(),ccnj(),cmk(),csum(),csub(),cmul(),cscl(),cdiv(),cexp(),cpow(),clog();
polar cpol();
dcomplex dcomp(),dccnj(),dcmk(),dcsum(),dcsub(),dcmul(),dcscl(),dcdiv(),dcexp(),dcpow(),dclog();
dpolar dcpol();

#endif /* clib_h */
