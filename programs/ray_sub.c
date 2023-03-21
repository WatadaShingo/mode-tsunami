#include <stdio.h>
#include <math.h>
#include "mallocutil.h"
#include "ray_sub.h"

#define FR	fprintf(stderr

void evanescent_level(level,turning_level,r_grid,v_grid,n_grid,omg,el,eps)
/*
for a given omg, el (angular order) and a velocity profile, evanescent
level where amplitude of a mode of becomes eps*(amplitude of mode at
turning point) v_grid[] can contain discontinuities where
consecutive r_grid[] have the same value.
<< OUTPUT >>
level: mode becomes negligible at this level
turning_level: equivalent ray turning level (fist level where
vertical wave number becomes imaginary from the center )
<< INPUT >>
r_grid[]: distance from the center of the earth. increasing order. in meter
v_grid[]: seismic velocity (in fluid P, in solid S wave velocity) in meter/sec
omg:	     angular frequency of mode ==2*PI/(period in sec)
l:        angular order of mode
eps: tolerance
*/
double *r_grid,*v_grid,omg,eps;
int *level,*turning_level,el,n_grid;
{
	double *kr2,tol,sum;
	int i,j;

	/* test input */
	if( eps <=0.0 )
	{FR,"evanescent_level(%s) %d: eps = %e less than 0.0.\n",
		__FILE__,__LINE__,eps);exit(-1);}
	if( el <= 0 ) 
	{FR,"evanescent_level(%s) %d: l = %d less than 0.\n",
		__FILE__,__LINE__,el);exit(-1);}
	if( n_grid <= 0 ) 
	{FR,"evanescent_level(%s) %d: n_grid = %d less than 0.\n",
		__FILE__,__LINE__,n_grid);exit(-1);}
	/* end test input */

	kr2=malloc_vdouble(0,n_grid-1,"evanescent_level kr");
	radial_wn(kr2,r_grid,v_grid,n_grid,omg,el);

	for(i=0;i<n_grid;i++)
	{
		if (r_grid[i] == 0.0) continue;
		if (kr2[i] > 0.0 ) break;
		kr2[i] = sqrt(-kr2[i]);
	}
	*turning_level = i-1;
	/* compute amplitude decay below the turning point */
	tol = -log(eps);
	sum=0.0;
	for(i=*turning_level;i>0;i--)
	{
		j=i-1;
		if( r_grid[i] == r_grid[j] ) continue;
		sum += 0.5*(kr2[i] + kr2[j])*(r_grid[i]-r_grid[j]);
		if( sum > tol ) break;
	}
	*level = j;
	free_vdouble(kr2,0);
}

void radial_wn(kr2_grid,r_grid,v_grid,n_grid,omg,el)
/*
 for a given omg, el (angular order) and a velocity profile, compute the radial
wave number.
v_grid[] can contain discontinuities where consecutive r_grid[] have the
same value.
kr(r)^2 = omg^2/v(r)^2 - l(l+1)/r^2 

kr(r)^2 >0 above turning point <0 below turning point.
<< OUTPUT >>
 kr[]:     square of radial wave number 1/meter
<< INPUT >>
 r_grid[]: distance from the center of the earth. increasing order. in meter
 v_grid[]: seismic velocity (in fluid P, in solid S wave velocity) in meter/sec
 omg:	     angular frequency of mode == 2*PI/(period in sec)
 el:        angular order of mode
 */
double *kr2_grid,*r_grid,*v_grid,omg;
int n_grid, el;
{
	int i;
	double v,r;

	if(kr2_grid == (double *)NULL){
		FR,"radial_wn(%s) %d: pointer kr2_grid is null\n",
		__FILE__,__LINE__);exit(-1);}
	if(v_grid == (double *)NULL){
		FR,"radial_wn(%s) %d: pointer v_grid is null\n",
		__FILE__,__LINE__);exit(-1);}
	if(r_grid == (double *)NULL){
		FR,"radial_wn(%s) %d: pointer r_grid is null\n",
		__FILE__,__LINE__);exit(-1);}

	for(i=0;i<n_grid;i++)
	{
		if( (v=v_grid[i]) <=0.0 ){
		FR,"radial_wn(%s) %d: negative velocity at r=%f v=%f grid=%d.\n",
		r_grid[i],v,i); exit(-1);}

		if( (r=r_grid[i]) <0.0 ){
		FR,"radial_wn(%s) %d: negative radius at r=%f grid=%d.\n",
		r,i); exit(-1);}
		if (r==0.0) {kr2_grid[i] = 0.0;continue;}
		kr2_grid[i] = omg*omg/(v*v) - el*(el+1.0)/(r*r);
	}
}
