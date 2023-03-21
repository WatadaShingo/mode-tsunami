#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "calc_a.h"
#include "m_atmos.h"
#include "read_atmos.h"
#include "m_solid.h"
#include "read_solid.h"

#define FR	fprintf(stderr
/*
	compute coefficents of linear differentioal equation
	of free oscillation of the earth.
	dz(r)/dr= A(z,earth_model, r)*z(r)
	z: eigen functions
	r: radius from the center
*/

void mult_m_v(zout,a,zin,n)
/*
	for each i, (0<=i,=n-1)

	zout[i] = sum_over_j(j=0;j<=n-1) a[i][j]*zin[j]	
*/
double *zout,(*a)[6],*zin;
int n;
{
	int i,j;
	register double dum;

	for(i=0;i<n;i++)
	{
		for(dum=0.0,j=0;j<n;j++)
			dum += a[i][j]*zin[j];
		zout[i]=dum;
	}
}
void mult_m_v6(zout,a,zin)
/*
	for each i, (0<=i,=n-1)

	zout[i] = sum_over_j(j=0;j<=n-1) a[i][j]*zin[j]	
*/
double *zout,a[][6],*zin;
{
	int i,j;
	register double dum;

	for(i=0;i<6;i++)
	{
		for(dum=0.0,j=0;j<6;j++)
			dum += a[i][j]*zin[j];
		zout[i]=dum;
	}
}
void mult_m_v4(zout,a,zin)
double *zout,a[][4],*zin;
{
	int i,j;
	register double dum;

	for(i=0;i<4;i++)
	{
		for(dum=0.0,j=0;j<4;j++)
			dum += a[i][j]*zin[j];
		zout[i]=dum;
	}
}
void mult_m_v2(zout,a,zin)
double *zout,a[][2],*zin;
{
	int i,j;
	register double dum;

	zout[0] = a[0][0]*zin[0]+a[0][1]*zin[1];
	zout[1] = a[1][0]*zin[0]+a[1][1]*zin[1];
/*	for(i=0;i<2;i++)
	{
		for(dum=0.0,j=0;j<2;j++)
			dum += a[i][j]*zin[j];
		zout[i]=dum;
	}
*/
}
void dzdr_atmos2by2(dzdr,z,r,l,omg)
	double *dzdr,/* dz(r)/dr */
	*z;		/* eigen function z(r) */
	int l;	/* angular order of mode */
	double r,	/* radius from center */
	omg;		/* angular frequency */
/*
	for given (l,omg,r,z) compute 1D array of dz(r)/dr
*/
{
	double a[2][2];

	calc_atmos2by2(a,r,l,omg);
	mult_m_v2(dzdr,a,z);
}
void calc_atmos2by2(a,r,l,omg2)
	double a[][2];
	int l;	/* angular order of mode */
	double r,	/* radius from center */
	omg2;		/* (angular frequency )^2*/
/*
	for given (l,omg,r) compute 2D array of A
	fluid medium without gravity potential perturbation
	(same as Cawling approximation)
*/
{
	double ll,r2;
	struct model_funcatmos m;	/* atmosphere model */


	find_modelatmos(&m,r);

	ll = l*(l+1.0);
	r2 = r*r;
	omg2/=ll;
	a[0][0] = -1.0/r + m.g/(omg2*r2);
	a[0][1] = 1.0/m.kappa - 1.0/(omg2*r2*m.rho);
	a[1][0] = m.rho*(-omg2*ll + m.g*m.g/(omg2*r2) - 4.0*m.g/r);
	a[1][1] = -a[0][0];
}
void y_initial_atmos2by2(y0,r,type)
	double *y0;
	double r;
	int type;
{
	/* fixed boundary condition */
	switch(type){
	case FIXED_BC:
		y0[0] = 0.0;
		y0[1] = 1.0;
	break;
	case FREE_BC:
		y0[0] = 1.0;
		y0[1] = 0.0;
	break;
	}
}

void get_eigenfunc_atmos2by2_from_x(u,dudr,x,dxdr,n_grid,r_grid,l,rho,drhodr,kappa,g,omg2)
/*
	construct eigenfunctions u,v,du/dr,dv/dr
	from normalized eigen functions x, dx/dr
	u[0] = U
	u[1] = V
	dudr[0] = dU/dr
	dudr[1] = dV/dr
*/
double **u,**dudr,**x,**dxdr,*r_grid,omg2;
/*float *rho,*drhodr,*kappa,*g;*/
double *rho,*drhodr,*kappa,*g;
int n_grid,l;
{
	int i;
	double r;
	double x3,sqrt_ll,ll;

	ll = l*(1.0+l);
	for(i=0;i<n_grid;i++)
	if( (r=r_grid[i])<=0.0 ) 
	{
		u[0][i]=u[1][i]=dudr[0][i]=dudr[1][i]=0.0;
	}
	else
	{
		u[0][i] = x[0][i]/r;
/*		dudr[0][i] = (dxdr[0][i] - u[0][i])/r; */
		/* from definition of X2 */
/*		u[1][i] = (x[1][i]/kappa[i] - r*dudr[0][i] - 2*u[0][i])/(-l*(l+1)); 
	## NOTE ##
	Since dx/dr is not computed in this method (see NOTE in solid_rk.c)
	use of dx/dr should be avoided.
	## NOTE END ## AUg. 26, 1994. see consrtut_eigen_uvf() in solid_sub.c
*/
/* this is test 
	from eqution of motion in fluid 
		sqrt_ll=sqrt(l*(l+1.0));
		x3 = sqrt_ll*(x[0][i]*g[1]-x[1][i]/rho[i])/(omg2*r);
		printf("x3/(V*sqrt(l(l+1))*r) = %e r=%e\n",
			x3/(u[1][i]*sqrt_ll*r),r);
	if solution is correct, the ratio should be 1.000000 
*/
		u[1][i] = (x[0][i]*g[1]-x[1][i]/rho[i])/(omg2*r*r);
		dudr[0][i] = (x[1][i]/kappa[i]-2.0*u[0][i]+ll*u[1][i])/r;
		dudr[1][i] = (dxdr[0][i]*g[i]-dxdr[1][i]/rho[i]+(rho[i]-2*g[i]/r)*x[0][i]
		+drhodr[i]*x[1][i]/(rho[i]*rho[i]))/(omg2*r*r) - 2*u[1][i]/r;
	}
	/* ## NOTE ##
	functin dudr[i] (== V) may not be accureate if
	drhodr[] is computed correctrly.
	For example,
	function rho[] is computed from linear interporatino of coarse
	earth model grid points.
	drhodr[] is also computed from linear interporation of coarse model
	grid points.
	But the coarse grid point of drhodr[] is computed from rho[] by
	the spline() routine.
	So near the coarse model grid points, the function drhodr[] differs
	significantly from the one computed directly from the function rho[]
	## NOTE END ## Aug 13, 1994
	*/
}
void get_x3_atmos2by2_from_x(x3,x,n_grid,r_grid,l,omg2,rho,g)
double *x3,**x,*r_grid,omg2;
double *rho,*g;
int n_grid,l;
{
	int i;
	double sqrt_ll;
	sqrt_ll = sqrt(l*(l+1.0));
		/* from eqution of motion in fluid */
	for(i=0;i<n_grid;i++)
		x3[i] = sqrt_ll*(x[0][i]*g[i]-x[1][i]/rho[i])/(omg2*r_grid[i]); 
}

void dzdr_solid6by6(dzdr,z,r,l,omg)
	double *dzdr,/* dz(r)/dr */
	*z;		/* eigen function z(r) */
	int l;	/* angular order of mode */
	double r,	/* radius from center */
	omg;		/* angular frequency */
/*
	for given (l,omg,r,z) compute 1D array of dz(r)/dr
*/
{
	double a[6][6];

	calc_solid6by6(a,r,l,omg);
	mult_m_v6(dzdr,a,z);
}
void dzdr_solid4by4cowl(dzdr,z,r,l,omg)
	double *dzdr,/* dz(r)/dr */
	*z;		/* eigen function z(r) */
	int l;	/* angular order of mode */
	double r,	/* radius from center */
	omg;		/* angular frequency */
/*
	for given (l,omg,r,z) compute 1D array of dz(r)/dr
*/
{
	double a[4][4];

	calc_solid4by4cowl(a,r,l,omg);
	mult_m_v4(dzdr,a,z);
}
void calc_solid6by6(a,r,l,omg2)
	double a[][6];
	int l;	/* angular order of mode */
	double r,	/* radius from center */
	omg2;		/* (angular frequency )^2*/
/*
	for given (l,omg,r) compute 6D array of A
	solid medium
*/
{
	double ll,r2,sqrt_ll,fc,acfn,rhogr;
	struct model_funcsolid m;	/* solid earth model */
	extern int region;

	find_modelsolid(&m,r,region);

	ll = l*(l+1.0);
	sqrt_ll = sqrt(ll);
	r2 = r*r;
	fc = m.f/m.c;
	acfn = m.a-m.f*fc-m.n;
	rhogr = m.rho*m.g*r;
	a[0][0] = (1.0-2.0*fc)/r;
	a[0][1] = 1.0/m.c;
	a[0][2] = sqrt_ll*fc/r;
	a[1][0] = -m.rho*omg2 + 4.0*(acfn-rhogr)/r2;
	a[1][1] = -a[0][0];
	a[1][2] = -sqrt_ll*(2.0*acfn-rhogr)/r2;
	a[1][3] = sqrt_ll/r;
	a[1][4] = -(l+1.0)*m.rho/r;
	a[1][5] = m.rho;
	a[2][0] = -a[1][3];
	a[2][2] =  2.0/r;
	a[2][3] = 1/m.l;
	a[3][0] = a[1][2];
	a[3][1] = -a[0][2];
	a[3][2] = -m.rho*omg2 + (ll*(acfn+m.n)-2.0*m.n)/r2;
	a[3][3] = -a[2][2];
	a[3][4] = sqrt_ll*m.rho/r;
	a[4][0] = -a[1][5];
	a[4][4] = -l/r;
	a[4][5] = 1.0;
	a[5][0] = a[1][4];
	a[5][2] = a[3][4];
	a[5][5] = -a[4][4];
	a[0][3] = a[0][4] = a[0][5] = a[2][1] = a[2][4] = a[2][5] = a[3][5]=
	a[4][1] = a[4][2] = a[4][3] = a[5][1] = a[5][3] = a[5][4] = 0.0; 
}
void calc_solid4by4cowl(a,r,l,omg2)
	double a[][4];
	int l;	/* angular order of mode */
	double r,	/* radius from center */
	omg2;		/* (angular frequency )^2*/
/*
	for given (l,omg,r) compute 4D array of A
	solid medium, with the Colwing approximation.
*/
{
	double ll,r2,sqrt_ll,fc,acfn,rhogr;
	struct model_funcsolid m;	/* solid earth model */
	extern int region;

	find_modelsolid(&m,r,region);

	ll = l*(l+1.0);
	sqrt_ll = sqrt(ll);
	r2 = r*r;
	fc = m.f/m.c;
	acfn = m.a-m.f*fc-m.n;
	rhogr = m.rho*m.g*r;
	a[0][0] = (1.0-2.0*fc)/r;
	a[0][1] = 1.0/m.c;
	a[0][2] = sqrt_ll*fc/r;
	a[1][0] = -m.rho*omg2 + 4.0*(acfn-rhogr)/r2;
	a[1][1] = -a[0][0];
	a[1][2] = -sqrt_ll*(2.0*acfn-rhogr)/r2;
	a[1][3] = sqrt_ll/r;
	a[2][0] = -a[1][3];
	a[2][2] =  2.0/r;
	a[2][3] = 1/m.l;
	a[3][0] = a[1][2];
	a[3][1] = -a[0][2];
	a[3][2] = -m.rho*omg2 + (ll*(acfn+m.n)-2.0*m.n)/r2;
	a[3][3] = -a[2][2];
	a[0][3] = a[2][1] = 0.0;
}
void dzdr_fluid4by4(dzdr,z,r,l,omg)
	double *dzdr,/* dz(r)/dr */
	*z;		/* eigen function z(r) */
	int l;	/* angular order of mode */
	double r,	/* radius from center */
	omg;		/* angular frequency */
/*
	for given (l,omg,r,z) compute 1D array of dz(r)/dr
*/
{
	double a[4][4];

	calc_fluid4by4(a,r,l,omg);
	mult_m_v4(dzdr,a,z);
}
void dzdr_fluid2by2cowl(dzdr,z,r,l,omg)
	double *dzdr,/* dz(r)/dr */
	*z;		/* eigen function z(r) */
	int l;	/* angular order of mode */
	double r,	/* radius from center */
	omg;		/* angular frequency */
/*
	for given (l,omg,r,z) compute 1D array of dz(r)/dr
*/
{
	double a[2][2];

	calc_fluid2by2cowl(a,r,l,omg);
	mult_m_v2(dzdr,a,z);
}
void calc_fluid4by4(a,r,l,omg2)
	double a[][4];
	int l;	/* angular order of mode */
	double r,	/* radius from center */
	omg2;		/* (angular frequency )^2*/
/*
	for given (l,omg,r) compute 4D array of A
	fluid medium
*/
{
	double ll,r2,lomgr;
	struct model_funcsolid m;	/* solid earth model */
	extern int region;

	find_modelsolid(&m,r,region);

	ll = l*(l+1.0);
	r2 = r*r;
	lomgr = ll/(omg2*r2); /* == l(l+1)/(omg^2*r^2) */
	a[0][0] = -1.0/r + m.g*lomgr;
	a[0][1] = 1.0/m.a - lomgr/m.rho;
	a[0][2] = lomgr;
	a[1][0] = m.rho*(-omg2 + m.g*m.g*lomgr - 4.0*m.g/r);
	a[1][1] = -a[0][0];
	a[1][2] = m.rho*(m.g*lomgr - (1.0+l)/r);
	a[1][3] = m.rho;
	a[2][0] = -a[1][3];
	a[2][2] = -l/r;
	a[2][3] = 1.0;
	a[3][0] = a[1][2];
	a[3][1] = -a[0][2];
	a[3][2] = lomgr*m.rho;
	a[3][3] = -a[2][2];
	a[0][3] = a[2][1] = 0.0;
}

void calc_fluid2by2cowl(a,r,l,omg2)
	double a[][2];
	int l;	/* angular order of mode */
	double r,	/* radius from center */
	omg2;		/* (angular frequency )^2*/
/*
	for given (l,omg,r) compute 2D array of A
	fluid medium
*/
{
	double ll,r2,lomgr;
	struct model_funcsolid m;	/* solid earth model */
	extern int region;

	find_modelsolid(&m,r,region);

	ll = l*(l+1.0);
	r2 = r*r;
	lomgr = ll/(omg2*r2); /* == l(l+1)/(omg^2*r^2) */
	a[0][0] = -1.0/r + m.g*lomgr;
	a[0][1] = 1.0/m.a - lomgr/m.rho;
	a[1][0] = m.rho*(-omg2 + m.g*m.g*lomgr - 4.0*m.g/r);
	a[1][1] = -a[0][0];
}
void y_initial_solid6by6(z1,z2,z3,l,r,omg2)
	double *z1,*z2,*z3;
	int l;	/* angular order of mode */
	double r,	/* radius from center */
	omg2;		/* (angular frequency )^2*/
/*
 * all earth model parameters, r and omg2 are nomralized.
 * see note for details
 */
{
	double ll,sqrt_ll,r2,a,f,pr2,qr2;
	double dum1,dum2,zz;
	struct model_funcsolid m;	/* solidearth model */
	double sbess_ratio();
	int region;

	set_regionsolid(&region,r,LOWER_REGION);
	find_modelsolid(&m,r,region);

	ll = l*(l+1.0);
	sqrt_ll = sqrt(ll);
	r2 = r*r;
	a = m.rho/3.0;
	/* normalized by 4*pi*UG because a has demension 1/time^2*/
	dum1 = omg2/m.l - (omg2 + 4*a)/(m.a);
	dum2 = sqrt(dum1*dum1+4.0*a*a*ll/(m.n*m.a));
	dum1 = omg2/m.l + (omg2 + 4*a)/(m.a);

	pr2 = 0.5*(dum1 - dum2)*m.rho;
	f = (pr2*m.l/m.rho - omg2)/a; /* dimension less */
	pr2 *= r2;

	zz = sbess_ratio((double)l,pr2);
	z1[0] = -f*zz;
	z1[1] = (-f*m.a + 2.0*m.l*(ll+2.0*f)*zz)/r;
	z1[2] = sqrt_ll*zz;
	z1[3] = sqrt_ll*m.l*(1.0 - 2.0*(f+1.0)*zz)/r;
	z1[4] = (-3.0*a*f - (omg2-a*l)*(f-l-1.0))*r/pr2;
	z1[5] = (2.0*l+1.0)*z1[4]/r;

	qr2 = 0.5*(dum1 + dum2)*m.rho;
	f = (qr2*m.l/m.rho - omg2)/a; /* dimension less */
	qr2 *= r2;

	zz = sbess_ratio((double)l,qr2);
	z2[0] = -f*zz;
	z2[1] = (-f*m.a + 2.0*m.l*(ll+2.0*f)*zz)/r;
	z2[2] = sqrt_ll*zz;
	z2[3] = sqrt_ll*m.l*(1.0 - 2.0*(f+1.0)*zz)/r;
	z2[4] = (-3.0*a*f - (omg2-a*l)*(f-l-1.0))*r/qr2;
	z2[5] = (2.0*l+1.0)*z2[4]/r;

	if( z3 == (double *)NULL) return;
	z3[0] = l;
	z3[1] = 2.0*m.l*l*(l-1.0)/r;
	z3[2] = sqrt_ll;
	z3[3] = 2.0*m.l*(l-1.0)*sqrt_ll/r;
	z3[4] = r*(omg2-a*l);
	z3[5] = ((2.0*l+1.0)*(omg2-a*l)+3.0*a*l);
}
void y_initial_fluid4by4(z1,z2,l,r,omg2)
	double *z1,*z2;
	int l;	/* angular order of mode */
	double r,	/* radius from center */
	omg2;		/* (angular frequency )^2*/
/*
 * all earth model parameters, r and omg2 are nomralized.
 * see note for details
 */
{
	double ll,r2,a,f,pr2;
	double zz;
	struct model_funcsolid m;	/* solidearth model */
	double sbess_ratio();
	int region;

	set_regionsolid(&region,r,LOWER_REGION);
	find_modelsolid(&m,r,region);

	ll = l*(l+1.0);
	r2 = r*r;
	a = m.rho/3.0;
	/* normalized by 4*pi*UG because a has demension 1/time^2*/

	pr2 = (omg2+4.0*a-a*a*ll/omg2)*m.rho/m.a*r2;
	f =  - omg2/a; /* dimension less */

	zz = sbess_ratio((double)l,pr2);
	z1[0] = -f*zz;
	z1[1] = (-f*m.a)/r;
	z1[2] = (-3.0*a*f - (omg2-a*l)*(f-l-1.0))*r/pr2;
	z1[3] = (2.0*l+1.0)*z1[2]/r;

	if( z2 == (double *)NULL) return;
	z2[0] = l;
	z2[1] = 0.0;
	z2[2] = r*(omg2-a*l);
	z2[3] = ((2.0*l+1.0)*(omg2-a*l)+3.0*a*l);
}
void continue_eigenfunc(xs,xf,type)
double xs[][6],xf[][4];
int type;
/*
 * continuation of eigenfunctions
	xs[i][j] i-th solution of eigenfunction j in solid i=1,2,3 j=1-6
	xf[i][j] i-th solution of eigenfunction j in fluid i=1,2 j=1-4
see Takeuchi & Saito (1972) p.254-255.
 */
{
	int i,j;
	double x4_13,x4_23;

	switch(type){
	case SOLID_TO_FLUID:
		x4_13 = xs[0][3]/xs[2][3];
		x4_23 = xs[1][3]/xs[2][3];
	i=0;	xf[0][i] = xs[0][i] - x4_13*xs[2][i];
		xf[1][i] = xs[1][i] - x4_23*xs[2][i];
	i=1;	xf[0][i] = xs[0][i] - x4_13*xs[2][i];
		xf[1][i] = xs[1][i] - x4_23*xs[2][i];
	i=2;	xf[0][i] = xs[0][4] - x4_13*xs[2][4];
		xf[1][i] = xs[1][4] - x4_23*xs[2][4];
	i=3;	xf[0][i] = xs[0][5] - x4_13*xs[2][5];
		xf[1][i] = xs[1][5] - x4_23*xs[2][5];
	break;
	case FLUID_TO_SOLID:
/* first and second solutions */
		j=0;
	i=0;	xs[i][j]=xf[i][j];
	i=1;	xs[i][j]=xf[i][j];
		j=1;
	i=0;	xs[i][j]=xf[i][j];
	i=1;	xs[i][j]=xf[i][j];
		j=2;
	i=0;	xs[i][4]=xf[i][j];
	i=1;	xs[i][4]=xf[i][j];
		j=3;
	i=0;	xs[i][5]=xf[i][j];
	i=1;	xs[i][5]=xf[i][j];
		xs[0][2] = xs[1][2] = xs[0][3] = xs[1][3] = 0.0;
/* third solutions */
		xs[2][2] = 1.0;
		xs[2][0] = xs[2][1] = xs[2][3] = xs[2][4] = xs[2][5]=0.0;
	break;
	case SOLID_TO_FLUID_COMPLEX:
	case FLUID_TO_SOLID_COMPLEX:
	default:
		FR,"continue_eigenfunc: unknown type=%d.\n",type);
		exit(-1);
	break;
	}
}
void continue_eigenfunccowl(xs,xf,type)
double xs[][6],xf[][4];
/* ### keep this dimension to be the interface consistent ### */
int type;
/*
 * continuation of eigenfunctions
	xs[i][j] i-th solution of eigenfunction j in solid i=1,2 j=1-4
	xf[i][j] i-th solution of eigenfunction j in fluid i=1 j=1-2
see Takeuchi & Saito (1972) p.254-255.
 */
{
	int i,j;
	double x4_12,x43_abs,ab1,ab2;

	switch(type){
	case SOLID_TO_FLUID:
		x4_12 = xs[0][3]/xs[1][3];
	i=0;	xf[0][i] = xs[0][i] - x4_12*xs[1][i];
	i=1;	xf[0][i] = xs[0][i] - x4_12*xs[1][i];
	break;
	case FLUID_TO_SOLID:
/* first solutions */
		xs[0][0]=xf[0][0];
		xs[0][1]=xf[0][1];
		xs[0][2] = xs[0][3] = 0.0;
/* second solutions */
		xs[1][2] = 1.0;
		xs[1][0] = xs[1][1] = xs[1][3] = 0.0;
	break;
	case SOLID_TO_FLUID_COMPLEX:
	x43_abs = (xs[2][3]*xs[2][3] + xs[2][7]*xs[2][7]);
	i=0;/* eigenfunction 1 cowling approx has only one solution*/
	j=i+4;
	ab1 = xs[0][3]*xs[2][i] - xs[0][7]*xs[2][j];/* ar*br-ai*bi */
	ab2 = xs[0][3]*xs[2][j] + xs[0][7]*xs[2][i];/* ar*bi+ai*br */
	xf[0][i] = xs[0][i] - ( ab1*xs[2][3] + ab2*xs[2][7])/x43_abs;
	xf[0][j] = xs[0][j] - (-ab1*xs[2][7] + ab2*xs[2][3])/x43_abs;

	i=1;/* eigenfunction 2 cowling approx has only one solution*/
	j=i+4;
	ab1 = xs[0][3]*xs[2][i] - xs[0][7]*xs[2][j];/* ar*br-ai*bi */
	ab2 = xs[0][3]*xs[2][j] + xs[0][7]*xs[2][i];/* ar*bi+ai*br */
	xf[0][i] = xs[0][i] - ( ab1*xs[2][3] + ab2*xs[2][7])/x43_abs;
	xf[0][j] = xs[0][j] - (-ab1*xs[2][7] + ab2*xs[2][3])/x43_abs;
	break;
	case FLUID_TO_SOLID_COMPLEX:
	break;
	default:
		FR,"continue_eigenfunccowl: unknown type=%d.\n",type);
		exit(-1);
	break;
	}
}

#define OUT(x,a,b)	(x<a || x>b)
double determinant3by3(i,j,k,y)
double y[][6]; /* y[i][j] i=1-3, j=1-6 */
int i,j,k;
/* i,j,k are 1-6 
|  y[0][i] y[1][i] y[2][i] |
|  y[0][j] y[1][j] y[2][j] |
|  y[0][k] y[1][k] y[2][k] |
*/
{
	double d;
	if( OUT(i,1,5) || OUT(j,1,5)|| OUT(k,1,5) ){
		FR,"determinant3by3: input error i,j,k=(%d,%d,%d).\n",i,j,k);
		exit(-1); }
	d = y[0][i]*(y[1][j]*y[2][k] - y[2][j]*y[1][k])
	  + y[1][i]*(y[2][j]*y[0][k] - y[0][j]*y[2][k])
	  + y[2][i]*(y[0][j]*y[1][k] - y[1][j]*y[0][k]);
	return d;
}
double determinant2by2(i,j,y)
double y[][4]; /* y[i][j] i=1-2, j=1-4 */
/* i,j,k are 1-4
|  y[0][i] y[1][i] |
|  y[0][j] y[1][j] |
*/
int i,j;
{
	double d;
	if( OUT(i,1,3) || OUT(j,1,3) ){
		FR,"determinant3by3: input error i,j=(%d,%d).\n",i,j);
		exit(-1);}
	d = y[0][i]*y[1][j] - y[1][i]*y[0][j];
	return d;
}
#undef OUT
double vector_cos(x,y,n)
double *x, *y;
int n;
{
	double xx,yy,a;
	int i;
	xx=yy=a=0.0;
	for(i=0;i<n;i++)
	{
		xx +=x[i]*x[i];
		yy +=y[i]*y[i];
		a +=x[i]*y[i];
	}
	return (a/sqrt(xx*yy));
}

void calc_fluid2by2_character(gamma,r,l,omg2,k,m)
	double *gamma;
	int l;	/* angular order of mode */
	double r,	/* radius from center */
	omg2;		/* (angular frequency )^2*/
	double k;	/* 1.0/(local scale height of density) */
	struct model_funcsolid m;
/*
	for given (l,omg,r) compute gamma, (growth rate of energy density)
	fluid medium. see note.
	if k < 0.0 then use the real local scale hight in the model
*/
{
	double omg2n, vgn,dum;

	non_dimensional_numbers(&omg2n, &vgn, r,&m);
	omg2/=omg2n;
	if( k < 0.0 ) k=- m.drhodr/m.rho;

	k *=  r;
	dum = k - 2.0*vgn + 4.0;
	*gamma = dum*dum + 4.0*(l*(l+1.0)/omg2-vgn)*(omg2-k+vgn);
}
void calc_fluid_omg(omg_buoy2, omg_lamb2,r,l,m)
/*
	for a given (l,r,m) compute buoyancy frequency (Brunt Vaisala freq)
	and Lamb frquency. in (angular grequency)
	l,    angular order of mode
	r,    radius
	m,    model
OUTPUT
	omg_buoy2, (buoyancy frequency in angular frequency)^2
	omg_lamb2, (Lamb frequency in angular frequency)^2
*/
	double *omg_buoy2,*omg_lamb2,r;
	int l;
	struct model_funcsolid m;
{
	double k,vgn,omg2n;

	k = -m.drhodr/m.rho;
	non_dimensional_numbers(&omg2n, &vgn, r,&m);

	k *=r;
	*omg_buoy2 = (k-vgn)*m.g/r;
	*omg_lamb2 = l*(l+1.0)/(r*r)*m.a/m.rho;
}

void boundary_coeff_fluid2by2(gamma,lam1,lam2,coeff1,coeff2,l,
				omg2_r,omg2_i,k,m,omg2n,vgn)
/* INPUT
	l,		angular order of mode
	omg2,		(angular frequency)^2
	k, 		1.0/(local scale height of density)
	m,		model
	omg2n,	frequency normalization factor GM/R^3 (g_bar/a)
	vgn,		Vg = rho g/c^2
   OUTPUT
	gamma,	[0] for real [1] for imaginary
	lam1,		0.5*(2+r*k + sqrt(gamma)) 
	lam2,		0.5*(2+r*k - sqrt(gamma))
			if gamma <0 or complex [0] is for real, [1] is for imaginary
			Definition of sqrt(gamma): Re(sqrt(gamma)) > 0
	*coeff boundary condition is given by
	sum_over_i (coeff1[i]*z[i])=0,
	complex case:
	i=0-1 for real i=2-3 for imaginary, lam1
	 coeff2[] is for lam2 solution
	real case:
	coeff1[0]=0 for lam1, coeff1[1] for lam2;
*/
int l;
double *gamma,*lam1,*lam2,*coeff1,*coeff2,omg2_r,omg2_i,k,omg2n,vgn;
struct model_funcsolid *m;
{
	double dum,dumabs,ll,rhog,a_dum,ivomg2_r,ivomg2_i;
	double r,thet;

	omg2_r/=omg2n;
	omg2_i/=omg2n;
	ll = l*(l+1.0);
	if( k < 0.0 ) k = -m->drhodr/m->rho;

	rhog = m->rho*m->g;
	k *= m->r;
	dum = k - 2.0*vgn + 4;
	a_dum = k-vgn;
	gamma[0] = dum*dum + 4.0*(ll/omg2_r-vgn)*(omg2_r-a_dum);
	if ( gamma[0] >=0.0 )
	{
		gamma[1] = 0.0;
		lam1[0] = 0.5*(2+k-sqrt(gamma[0]));
		lam1[1] = lam2[0] = lam2[1] = 0.0;

		coeff1[0] = rhog*(-1.0 + ll/omg2_r - lam1[0]);
		coeff1[1] = vgn - ll/omg2_r;
		coeff1[2] = coeff1[3] = 0.0;
		coeff2[0] = coeff2[1] = coeff2[2] = coeff2[3] = 0.0;
		dumabs = sqrt(coeff1[0]*coeff1[0] + coeff1[1]*coeff1[1]);
		coeff1[0] /= dumabs;
		coeff1[1] /= dumabs;
	}else{
		dumabs = omg2_r*omg2_r + omg2_i*omg2_i;
		ivomg2_r =  omg2_r/dumabs;
		ivomg2_i = -omg2_i/dumabs;
		gamma[0] = dum*dum + 4.0*(ll+(vgn-ll*ivomg2_r)*a_dum-vgn*omg2_r);
		gamma[1] = 4.0*(-vgn*omg2_i - ll*a_dum*ivomg2_i);
		r        = sqrt(gamma[0]*gamma[0] + gamma[1]*gamma[1]);
		thet     = atan2(gamma[1],gamma[0])/2.0;
		lam1[0] = r*cos(thet);
		lam1[1] = r*sin(thet);
		if( lam1[1] < 0.0 ){lam1[0] = -lam1[0];lam1[1] = - lam1[1];}
		lam2[0] = -lam1[0];
		lam2[1] = -lam1[1];

		lam1[0] = 0.5*(2+k+lam1[0]);
		lam1[1] = 0.5*lam1[1];
		lam2[0] = 0.5*(2+k+lam2[0]);
		lam2[1] = 0.5*lam2[1];

		coeff1[0] = rhog*(-1.0 + ll*ivomg2_r - lam1[0]) ;
		coeff2[0] = rhog*(-1.0 + ll*ivomg2_r - lam2[0]) ;
		coeff1[2] = rhog*(ll*ivomg2_i - lam1[1] );
		coeff2[2] = rhog*(ll*ivomg2_i - lam2[1] );
		coeff1[1] = coeff2[1] = vgn - ll*ivomg2_r;
		coeff1[3] = coeff2[3] = -ll*ivomg2_i;
		dumabs = sqrt(coeff1[0]*coeff1[0]+coeff1[1]*coeff1[1]+
		            coeff1[2]*coeff1[2]+coeff1[3]*coeff1[3]);
		coeff1[0] /=dumabs;
		coeff1[1] /=dumabs;
		coeff1[2] /=dumabs;
		coeff1[3] /=dumabs;
		dumabs = sqrt(coeff2[0]*coeff2[0]+coeff2[1]*coeff2[1]+
		            coeff2[2]*coeff2[2]+coeff2[3]*coeff2[3]);
		coeff2[0] /=dumabs;
		coeff2[1] /=dumabs;
		coeff2[2] /=dumabs;
		coeff2[3] /=dumabs;
	}
}
