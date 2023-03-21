#include <stdio.h>
#include <math.h>
#include "m_solid.h"
#include "solid_sub.h"

void construct_eigen_x0(yout,yin)
/*
	## RECOMMENDATION ##
	DO NOT USE THIS FUNCTION. USE construct_eigen_x() INSTEAD.
	both functions have the same functionality but construct_eigen_x0()
	requires external struct model_solid.
	## RECOMMENDATION END ##
*/
double **yout,***yin;
/*
Construct eigeunfunction which satisfy free surface boundary condition
from three (in fluid two) independent eigenfunctions.
<< output >>
yout[6][]   eigenfunction which satisfy free surface
<< input >>
yin[3][6][] three independent solutions
*/
{
	int i,j,k,n,upper_layer,lower_layer,old_upper_layer;
	double q[3];
	extern struct model_solid msolid;

	if( msolid.value[N_SOLID][msolid.nmodel-1] == 0.0 ) /* fluid top */
	{
		upper_layer = msolid.region[msolid.nregion-1][UPPER_BOUND_SOLID];
		lower_layer = msolid.region[msolid.nregion-1][LOWER_BOUND_SOLID];
		q[0] = 1.0;
		q[1] = -yin[0][1][upper_layer]/yin[0][1][upper_layer]*q[0];
		for(j=lower_layer;j<=upper_layer;j++)
		{
			yout[0][j] = q[0]*yin[0][0][j];
			yout[1][j] = q[0]*yin[0][1][j];
			yout[4][j] = q[0]*yin[0][4][j];
			yout[5][j] = q[0]*yin[0][5][j];
			for(k=1;k<2;k++)
			{
				yout[0][j] += q[k]*yin[k][0][j];
				yout[1][j] += q[k]*yin[k][1][j];
				yout[4][j] += q[k]*yin[k][4][j];
				yout[5][j] += q[k]*yin[k][5][j];
			}
		}
	for(i=msolid.nregion-1;i>0;i--)
	{
		i--;
		old_upper_layer = upper_layer;
		upper_layer = msolid.region[i][UPPER_BOUND_SOLID];
		lower_layer = msolid.region[i][LOWER_BOUND_SOLID];
		if( msolid.value[N_SOLID][upper_layer] == 0.0 )
		/* in fluid region */
		for(j=lower_layer;j<=upper_layer;j++)
		{
			yout[0][j] = q[0]*yin[0][0][j];
			yout[1][j] = q[0]*yin[0][1][j];
			yout[4][j] = q[0]*yin[0][4][j];
			yout[5][j] = q[0]*yin[0][5][j];
			for(k=1;k<2;k++)
			{
				yout[0][j] += q[k]*yin[k][0][j];
				yout[1][j] += q[k]*yin[k][1][j];
				yout[4][j] += q[k]*yin[k][4][j];
				yout[5][j] += q[k]*yin[k][5][j];
			}
		}
		else 
		{
			/* in solid region just beneath of fluid region */
			if( msolid.value[N_SOLID][old_upper_layer] == 0.0 )
			q[2] = -q[0]*yin[0][3][upper_layer]/yin[2][3][upper_layer] 
				- q[1]*yin[1][3][upper_layer]/yin[2][3][upper_layer];
		for(j=lower_layer;j<=upper_layer;j++)
		for(n=0;n<6;n++)
		{
			yout[n][j] = q[0]*yin[0][n][j];
			for(k=1;k<3;k++)
			yout[n][j] += q[k]*yin[k][n][j];
		}
		} /* end of if else */
	}/* end of region for */
	}else{ /* solid top */
		upper_layer = msolid.region[msolid.nregion-1][UPPER_BOUND_SOLID];
		q[0] = yin[1][1][upper_layer]*yin[2][3][upper_layer]
		     - yin[2][1][upper_layer]*yin[1][3][upper_layer];
		q[1] = yin[2][1][upper_layer]*yin[0][3][upper_layer]
		     - yin[0][1][upper_layer]*yin[2][3][upper_layer];
		q[2] = yin[0][1][upper_layer]*yin[1][3][upper_layer]
		     - yin[1][1][upper_layer]*yin[0][3][upper_layer];
		q[1] /=q[0];
		q[2] /=q[0];
		q[0] = 1.0;
	for(i=msolid.nregion;i>0;i--)
	{
		i--;
		if( msolid.value[N_SOLID][upper_layer] == 0.0 )
		/* in fluid region */
		for(j=lower_layer;j<=upper_layer;j++)
		{
			yout[0][j] = q[0]*yin[0][0][j];
			yout[1][j] = q[0]*yin[0][1][j];
			yout[4][j] = q[0]*yin[0][4][j];
			yout[5][j] = q[0]*yin[0][5][j];
			for(k=1;k<2;k++)
			{
				yout[0][j] += q[k]*yin[k][0][j];
				yout[1][j] += q[k]*yin[k][1][j];
				yout[4][j] += q[k]*yin[k][4][j];
				yout[5][j] += q[k]*yin[k][5][j];
			}
		}
		else 
		/* in solid layer */
		for(j=lower_layer;j<=upper_layer;j++)
		for(n=0;n<6;n++)
		{
			yout[n][j] = q[0]*yin[0][n][j];
			for(k=1;k<3;k++)
			yout[n][j] += q[k]*yin[k][n][j];
		}
	}/* end of for */
	}/* end of if else */
}
void construct_eigen_x(yout,yin,el,grid_start,grid_end)
double **yout,***yin;
/*float *el;*/
double *el;
int grid_start,grid_end;
/*
Construct eigeunfunction which satisfy free surface boundary condition
from three (in fluid two) independent eigenfunctions.
<< output >>
yout[6][]   eigenfunction which satisfy free surface
yout[j][k]  eigenfunction j at grid k
<< input >>
yin[3][6][] three independent solutions
yin[i][j][k]  i-th solution of eigenfunction j at grid k
el[],	elastic constant L (in fluid zero )
eigenfunction j (0<=j<=5) corresponds z_j in note p.17.
*/
{
	int j,k,n;
	double q[3];

	if( el[grid_end] == 0.0 ) /* fluid top *//* free sufrace is assumed */
	{
		q[0] = 1.0;
		q[1] = -q[0]*yin[0][1][grid_end]/yin[1][1][grid_end];
		yout[0][grid_end] = q[0]*yin[0][0][grid_end];
		yout[1][grid_end] = q[0]*yin[0][1][grid_end];
		yout[4][grid_end] = q[0]*yin[0][4][grid_end];
		yout[5][grid_end] = q[0]*yin[0][5][grid_end];
		for(k=1;k<2;k++)
		{
			yout[0][grid_end] += q[k]*yin[k][0][grid_end];
			yout[1][grid_end] += q[k]*yin[k][1][grid_end];
			yout[4][grid_end] += q[k]*yin[k][4][grid_end];
			yout[5][grid_end] += q[k]*yin[k][5][grid_end];
		}
		grid_end--;
		for(j=grid_end;j>=grid_start;j--)
		{
			/* fluid over solid region boundary */
			if( el[j+1] == 0.0 && el[j] !=0.0 )
			q[2] = -q[0]*yin[0][3][j]/yin[2][3][j] 
				- q[1]*yin[1][3][j]/yin[2][3][j];
			if( el[j] == 0.0 )
			/* in fluid */
			{
				yout[0][j] = q[0]*yin[0][0][j];
				yout[1][j] = q[0]*yin[0][1][j];
				yout[4][j] = q[0]*yin[0][4][j];
				yout[5][j] = q[0]*yin[0][5][j];
				for(k=1;k<2;k++)
				{
					yout[0][j] += q[k]*yin[k][0][j];
					yout[1][j] += q[k]*yin[k][1][j];
					yout[4][j] += q[k]*yin[k][4][j];
					yout[5][j] += q[k]*yin[k][5][j];
				}
			}else{
			/* in solid */
				for(n=0;n<6;n++)
				{
					yout[n][j] = q[0]*yin[0][n][j];
					for(k=1;k<3;k++)
					yout[n][j] += q[k]*yin[k][n][j];
				}
			}
		}
	}else{ /* solid top *//* free surface is assumed */
		q[0] = yin[1][1][grid_end]*yin[2][3][grid_end]
		     - yin[2][1][grid_end]*yin[1][3][grid_end];
		q[1] = yin[2][1][grid_end]*yin[0][3][grid_end]
		     - yin[0][1][grid_end]*yin[2][3][grid_end];
		q[2] = yin[0][1][grid_end]*yin[1][3][grid_end]
		     - yin[1][1][grid_end]*yin[0][3][grid_end];
		q[1] /=q[0];
		q[2] /=q[0];
		q[0] = 1.0;
		for(n=0;n<6;n++)
		{
			yout[n][grid_end] = q[0]*yin[0][n][grid_end];
			for(k=1;k<3;k++)
			yout[n][grid_end] += q[k]*yin[k][n][grid_end];
		}
		grid_end --;
		for(j=grid_end;j>=grid_start;j--)
		{
			/* fluid over solid region boundary */
			if( el[j+1] == 0.0 && el[j] !=0.0 )
			q[2] = -q[0]*yin[0][3][j]/yin[2][3][j] 
				- q[1]*yin[1][3][j]/yin[2][3][j];
			if( el[j] == 0.0 )
			/* in fluid */
			{
				yout[0][j] = q[0]*yin[0][0][j];
				yout[1][j] = q[0]*yin[0][1][j];
				yout[4][j] = q[0]*yin[0][4][j];
				yout[5][j] = q[0]*yin[0][5][j];
				for(k=1;k<2;k++)
				{
					yout[0][j] += q[k]*yin[k][0][j];
					yout[1][j] += q[k]*yin[k][1][j];
					yout[4][j] += q[k]*yin[k][4][j];
					yout[5][j] += q[k]*yin[k][5][j];
				}
			}else{
			/* in solid */
				for(n=0;n<6;n++)
				{
					yout[n][j] = q[0]*yin[0][n][j];
					for(k=1;k<3;k++)
					yout[n][j] += q[k]*yin[k][n][j];
				}
			}
		}/* end of for */
	}/* end of if else */
}
void construct_eigen_xcowl(yout,yin,el,grid_start,grid_end)
double **yout,***yin;
double *el;
int grid_start,grid_end;
/*
Construct eigeunfunction which satisfy free surface boundary condition
from three (in fluid two) independent eigenfunctions.
<< output >>
yout[6][]   eigenfunction which satisfy free surface
yout[j][k]  eigenfunction j at grid k
<< input >>
yin[3][6][] three independent solutions
yin[i][j][k]  i-th solution of eigenfunction j at grid k
el[],	elastic constant L (in fluid zero )
eigenfunction j (0<=j<=5) corresponds Z_j in note p.17.
*/
/* ### keep this dimension to be the interface consistent ### */
{
	int j,k,n;
	double q[2];

	if( el[grid_end] == 0.0 ) /* fluid top *//* not necessarily free surface */
	{
		q[0] = 1.0;
		yout[0][grid_end] = yin[0][0][grid_end];
		yout[1][grid_end] = yin[0][1][grid_end];
		grid_end--;
		for(j=grid_end;j>=grid_start;j--)
		{
			/* fluid over solid region boundary */
			if( el[j+1] == 0.0 && el[j] !=0.0 )
			q[1] = -q[0]*yin[0][3][j]/yin[1][3][j] ;
			if( el[j] == 0.0 )
			/* in fluid */
			{
				yout[0][j] = yin[0][0][j];
				yout[1][j] = yin[0][1][j];
			}else{
			/* in solid */
				for(n=0;n<4;n++)
				{
					yout[n][j] = q[0]*yin[0][n][j];
					for(k=1;k<2;k++)
					yout[n][j] += q[k]*yin[k][n][j];
				}
			}
		}
	}else{ /* solid top *//* free surface is assumed */
		q[0] = 1.0;
		q[1] = -q[0]*yin[0][3][grid_end]/yin[1][3][grid_end];
		for(n=0;n<4;n++)
		{
			yout[n][grid_end] = q[0]*yin[0][n][grid_end];
			for(k=1;k<2;k++)
			yout[n][grid_end] += q[k]*yin[k][n][grid_end];
		}
		grid_end --;
		for(j=grid_end;j>=grid_start;j--)
		{
			/* fluid over solid region boundary */
			if( el[j+1] == 0.0 && el[j] !=0.0 )
			q[1] = -q[0]*yin[0][3][j]/yin[1][3][j];
			if( el[j] == 0.0 )
			/* in fluid */
			{
				yout[0][j] = yin[0][0][j];
				yout[1][j] = yin[0][1][j];
			}else{
			/* in solid */
				for(n=0;n<4;n++)
				{
					yout[n][j] = q[0]*yin[0][n][j];
					for(k=1;k<2;k++)
					yout[n][j] += q[k]*yin[k][n][j];
				}
			}
		}/* end of for */
	}/* end of if else */
}

void construct_eigen_xradial(yout,yin,rho,r_grid,grid_start,grid_end,grid_surface)
double **yout,**yin;
double *rho;
double *r_grid;
int grid_start,grid_end,grid_surface;
/*
Construct eigeunfunction which satisfy free surface boundary condition
<< output >>
yout[6][]   eigenfunction which satisfy free surface
yout[i][k]  eigenfunction j at grid k
<< input >>
yin[6][] independent solutions
yin[j][k]  solution of eigenfunction j at grid k
rho[] el[], density, elastic constant L (in fluid zero )
grid_start<=k<=grid_end		eigenfunction at k is computed 
grid_surace,	grid number of surface.
r_grid[], radius at grid points.

eigenfunction j (0<=j<=5) corresponds z_j in note p.17.
j=0,1 are given in yin[][]. consturct j=4,5 in this function.
*/
/* ### keep this dimension to be the interface consistent ### */
{
	int i,j,k;
	double val1,val2,r1,r2,total;

	for(i=grid_start;i<=grid_surface;i++)
	{
		yout[0][i]=yin[0][i];
		yout[1][i]=yin[1][i];
		yout[2][i]=(r_grid[i] <= 0.0 ? 0.0 : yin[0][i]/r_grid[i]);
		/* yout[2][] is not in use. Use it as temporal storage */ 
	}
/*                               |r'=a
   potential perturbation(r) =   |  rho(r') U(r') dr'
				         |r'=r
*/
	yout[5][grid_surface]=0.0;
	total = 0.0;
	for(i=grid_end;i<grid_surface;i++)
	{
		j=i+1;
		r1 = r_grid[i];
		r2 = r_grid[j];
		if( r1==r2 ) continue;
		val1 = rho[i]*yout[2][i];
		val2 = rho[j]*yout[2][j];
		total += (val1+val2)*0.5*(r2-r1);
	}
	yout[5][grid_end] = total; 

	for(i=grid_end-1;i>=grid_start;i--)
	{
		j=i+1;
		r1 = r_grid[i];
		r2 = r_grid[j];
		if( r1==r2 ){yout[5][i]=yout[5][j]; continue;}
		val1 = rho[i]*yout[2][i];
		val2 = rho[j]*yout[2][j];
		yout[5][i] = yout[5][j] + (val1+val2)*0.5*(r2-r1);
	}

/*                               |r'=r
   potential perturbation(r) = - |  rho(r') U(r') dr'
				         |r'=0
  ## NOTE ## physically both methods are equivalent. both satisfies
  d potential perturbation(r)
  -------------------------- = 4*PI*G*rho(r)*U(r)
  dr
  but this one does NOTnot satisfy the energy integral I1*omg^2 = I2 
  ## END NOTE ##

	yout[5][grid_start]=0.0;
	for(i=grid_start;i<grid_end;i++)
	{
		j=i+1;
		r1 = r_grid[i];
		r2 = r_grid[j];
		if( r1==r2 ){yout[5][j]=yout[5][i]; continue;}
		val1 = rho[i]*yout[2][i];
		val2 = rho[j]*yout[2][j];
		yout[5][j] = yout[5][i] - (val1+val2)*0.5*(r2-r1);
	}
*/


	for(i=grid_start;i<=grid_end;i++)
	{
		yout[4][i] = yout[5][i]*r_grid[i];
		yout[2][i] = yout[3][i] = 0.0;
	}
}

void construct_eigen_uvfradial(u,du,v,dv,fi,dfi,x,r_grid,c,el,f,rho,grid_start,grid_end)
double *u,*du,*v,*dv,*fi,*dfi,*r_grid,**x;
double *c,*el,*f,*rho;
int grid_start,grid_end;
/*
	convert eigenfunctions from x1,x2,x3
	to u,v,fi,du/dr,dv/dr,dfi/dr
	x[j][k] is eigenfunction j at grid k.
	eigenfunction j (0<=j<=5) corresponds z_j in note p.17.
	## NOTE ##
	in fluid region dv/dr is not set.
	dv/dr will be set by 1) analytical form or 2) interpolation of V
	## NOTE END ##
*/
{
	int i;
	double r;
	for(i=grid_start;i<=grid_end;i++)
	{
		r = r_grid[i];
		if( r==0.0 )
			u[i]=du[i]=fi[i]=dfi[i] = 0.0;
		else if( el[i] != 0.0 ) /* in solid */
		{
			u[i] = x[0][i]/r;
			du[i] = (x[1][i]-2.0*f[i]*u[i])/(c[i]*r);
			fi[i] = x[5][i];
			dfi[i] =  -rho[i]*u[i];
		}else{
			u[i] = x[0][i]/r;
			du[i] = (x[1][i]/c[i]-2.0*u[i])/r;
			fi[i] = x[5][i];
			dfi[i] =  -rho[i]*u[i];
		}
		v[i]=dv[i]=0.0;
	}
}

void construct_eigen_uvf(u,du,v,dv,fi,dfi,x,dx,r_grid,l,omg2,c,el,f,rho,g,grid_start,grid_end)
double *u,*du,*v,*dv,*fi,*dfi,omg2,*r_grid,**x,**dx;
double *c,*el,*f,*rho,*g;
int l,grid_start,grid_end;
/*
	convert eigenfunctions from x1,x2,x3,dx1/dr,dx2/dr and dx3/dr 
	to u,v,fi,du/dr,dv/dr,dfi/dr
	x[j][k] is eigenfunction j at grid k.
	eigenfunction j (0<=j<=5) corresponds z_j in note p.17.
	## NOTE ##
	in fluid region dv/dr is not set.
	dv/dr will be set by 1) analytical form or 2) interpolation of V
	## NOTE END ##
*/
{
	int i;
	double r,ll,sqrt_ll;

	ll = l*(l+1.0);
	sqrt_ll=sqrt(ll);
	for(i=grid_start;i<=grid_end;i++)
	{
		r=r_grid[i];
		if( r == 0.0 )
			u[i]=v[i]=du[i]=dv[i]=fi[i]=dfi[i] = 0.0;
		else if( el[i] != 0.0 ) /* in solid */
		{
			u[i]  = x[0][i]/r;
			v[i]  = x[2][i]/(r*sqrt_ll);
			du[i] = (x[1][i]-f[i]*(2.0*u[i]-ll*v[i]))/(c[i]*r);
			dv[i] = (x[3][i]/(sqrt_ll*el[i]) - (u[i]-v[i]))/r;
			fi[i] = x[4][i]/r;
			dfi[i]= (x[5][i]-rho[i]*r*u[i]-(l+1.0)*fi[i])/r;
		}else{ /* in fluid */
			u[i]  = x[0][i]/r;
/*			du[i] = (dx[0][i]-u[i])/r; 
			v[i]  = (x[1][i]/c[i]-r*du[i]-2.0*u[i])/(-ll);
	## NOTE ##
	Since dx/dr is not computed in this method (see NOTE in solid_rk.c) 
	use of dx/dr should be avoided.
	## NOTE END ## AUg. 26, 1994.
*/
			v[i] =(g[i]*x[0][i]-x[1][i]/rho[i]+x[4][i])/(omg2*r*r);
			du[i] = (x[1][i]/c[i]-2.0*u[i]+ll*v[i])/r;
			fi[i] = x[4][i]/r;
			dfi[i]= (x[5][i]-rho[i]*r*u[i]-(l+1.0)*fi[i])/r;
		}
	}
}
void construct_eigen_uvcowl(u,du,v,dv,x,dx,r_grid,l,omg2,c,el,f,rho,g,grid_start,grid_end)
double *u,*du,*v,*dv,omg2,*r_grid,**x,**dx;
double *c,*el,*f,*rho,*g;
int l,grid_start,grid_end;
/*
	convert eigenfunctions from x1,x2,x3,dx1/dr,dx2/dr and dx3/dr 
	to u,v,fi,du/dr,dv/dr,dfi/dr
	x[j][k] is eigenfunction j at grid k.
	eigenfunction j (0<=j<=5) corresponds z_j in note p.17.
	## NOTE ##
	in fluid region dv/dr is not set.
	dv/dr will be set by 1) analytical form or 2) interpolation of V
	## NOTE END ##
*/
{
	int i;
	double r,ll,sqrt_ll;

	ll = l*(l+1.0);
	sqrt_ll=sqrt(ll);
	for(i=grid_start;i<=grid_end;i++)
	{
		r=r_grid[i];
		if( r == 0.0 )
			u[i]=v[i]=du[i]=dv[i] = 0.0;
		else if( el[i] != 0.0 ) /* in solid */
		{
			u[i]  = x[0][i]/r;
			v[i]  = x[2][i]/(r*sqrt_ll);
			du[i] = (x[1][i]-f[i]*(2.0*u[i]-ll*v[i]))/(c[i]*r);
			dv[i] = (x[3][i]/(sqrt_ll*el[i]) - (u[i]-v[i]))/r;
		}else{ /* in fluid */
			u[i]  = x[0][i]/r;
/*			du[i] = (dx[0][i]-u[i])/r; 
			v[i]  = (x[1][i]/c[i]-r*du[i]-2.0*u[i])/(-ll);
	## NOTE ##
	Since dx/dr is not computed in this method (see NOTE in solid_rk.c) 
	use of dx/dr should be avoided.
	## NOTE END ## AUg. 26, 1994.
*/
			v[i] =(g[i]*x[0][i]-x[1][i]/rho[i])/(omg2*r*r);
			du[i] = (x[1][i]/c[i]-2.0*u[i]+ll*v[i])/r;
		}
	}
}
void construct_eigen_x3(x3,x4,x1,x2,x5,r_grid,l,omg2,el,rho,g,grid_start,grid_end)
double *x3,*x4,*x1,*x2,*x5,*r_grid,omg2;
double *rho,*g,*el;
int l,grid_start,grid_end;
/*
	set x3,x4 variables in fluid region
	Saito (1988) eq.47, Woodhouse (1988) eq.5.38 or see note p.24
*/
{
	int i;
	double sqrt_ll;

	sqrt_ll = sqrt(l*(l+1.0));
	/* from equation of motion in fluid */
	for(i=grid_start;i<=grid_end;i++)
	if( el[i] == 0.0 )
	{
		x3[i] = sqrt_ll*(g[i]*x1[i]-x2[i]/rho[i]+x5[i])/(omg2*r_grid[i]);
		x4[i] = 0.0;
	}
}
void construct_eigen_x3cowl(x3,x4,x1,x2,r_grid,l,omg2,el,rho,g,grid_start,grid_end)
double *x3,*x4,*x1,*x2,*r_grid,omg2;
double *rho,*g,*el;
int l,grid_start,grid_end;
/*
	set x3,x4 variables in fluid region
	Saito (1988) eq.47, Woodhouse (1988) eq.5.38 or see note p.24
*/
{
	int i;
	double sqrt_ll;

	sqrt_ll = sqrt(l*(l+1.0));
	/* from equation of motion in fluid */
	for(i=grid_start;i<=grid_end;i++)
	if( el[i] == 0.0 )
	{
		x3[i] = sqrt_ll*(g[i]*x1[i]-x2[i]/rho[i])/(omg2*r_grid[i]);
		x4[i] = 0.0;
	}
}

void construct_energy(energy_comp,energy_shear,u,du,rho,a,c,el,en,f,l,r_grid,grid_start,grid_end)
/* construct energy eigenfunction from eigenfunctions and
earth model parameters 
energy_comp,	compression energy density
energy_shear,	shear energy density
u[0][]	U
u[1][]	V
u[2][]	PHI
du[0][]	dU/dr
du[1][]	dV/dr
du[2][]	dPHI/dr
rho[][]	density
eigenfunction in normalized according to
|r=a
|rho*(U^2+l(l+1.0)*V^2) = 1.0
|r=0
*/
double *energy_comp,*energy_shear;
double *r_grid,**u,**du,*rho;
double *a,*c,*el,*en,*f;
int l,grid_start,grid_end;
{
	int i;
	double ll,r,kappa,mue,c1,c2,u0,v0,du0,dv0;
	double du0_r,dv0_r,ullv;

	ll = l*(1.0+l);
	for(i=grid_start;i<=grid_end;i++)
	{
		r = r_grid[i];
		u0 = u[0][i];
		v0 = u[1][i];
		du0_r = du[0][i]*r;
		dv0_r = du[1][i]*r;
		ullv = 2.0*u0-ll*v0;
		kappa = (4.0*(a[i]+f[i]-en[i])+c[i])/9.0;
		mue = (a[i]+c[i]-2.0*f[i]+5.0*en[i]+6.0*el[i])/15.0;
		c1 = kappa/(kappa+4.0/3.0*mue);
		c2 = kappa/(kappa-2.0/3.0*mue);
		energy_comp[i] = c1*c[i]*du0_r*du0_r
				   + 2.0*c2*f[i]*du0_r*ullv
				   + c1*a[i]*ullv*ullv;
		energy_shear[i] = (1.0-c1)*c[i]*du0_r*du0_r
				   + 2.0*(1.0-c2)*f[i]*du0_r*ullv
				   + (1.0-c1)*a[i]*ullv*ullv
				   + el[i]*ll*(u0-v0+dv0_r)*(u0-v0+dv0_r)
				   + en[i]*(ll*(l+2.0)*(l-1.0)*v0*v0-ullv*ullv);
	}
}

void integral_q(q,ratio_comp,ratio_shear,e_comp,e_shear,qk,qu,r_grid,grid_start,grid_end)
/* compute 1/Q for a given modal compressional and shear elastic
density and attenuation structure 1/Qkappa, 1/Qmue*/
double *q,*ratio_comp,*ratio_shear,*e_comp,*e_shear,*qk,*qu,*r_grid;
int grid_start,grid_end;
{
	double r1,r2,val1,val2;
	int i,j;

	*q=*ratio_comp=*ratio_shear=0.0;
	for(i=grid_start;i<grid_end;i++)
	{
		j = i+1;
		r1 = r_grid[i];
		r2 = r_grid[j];
		if( r1==r2 ) continue;
		val1 = qk[i]*e_comp[i]+qu[i]*e_shear[i];
		val2 = qk[j]*e_comp[j]+qu[j]*e_shear[j];
		*ratio_comp += (e_comp[i]+e_comp[j])*0.5*(r2-r1);
		*ratio_shear += (e_shear[i]+e_shear[j])*0.5*(r2-r1);
		*q += (val1+val2)*0.5*(r2-r1);
	}
}
