#include <stdio.h>
#include <math.h>
#include "earthconst.h"
#include "model.h"

#define PI	3.14159265358979324

double mass_r(rho,r_grid,grid_start,grid_end)
/* compute total mass of a spherically symmetric spghere
from grid radius data
                | r_grid[grid_end]
total_mass =    | rho(r)*r*r dr
                | r_grid[start]

to get real total mass scale_rho*scale_rad^3 should be
multiplied.
r_grid[] should be increasing order
if( grid_start >= grid_end
*/
double rho[],r_grid[];
int grid_start,grid_end;
{
	double total_mass,r1,r2,rho1,rho2;
	int i;

	total_mass=0.0;
	for(i=grid_start;i<grid_end;i++)
	{
		r1=r_grid[i];
		r2=r_grid[i+1];
		if( r1==r2) continue;
/*
		total_mass += (rho[i]*r1*r1+rho[i+1]*r2*r2)*0.5*(r2-r1);
*/
		rho1=rho[i];
		rho2=rho[i+1];
		total_mass +=(r2*rho1-r1*rho2)*(r2*r2+r1*r1+r1*r2)/3.0
				+(rho2-rho1)*(r2*r2*r2+r2*r2*r1+r2*r1*r1+r1*r1*r1)/4.0;
	}
	return 4*PI*total_mass;
}

double grav_r(g0,rho,r_grid,grid_start,grid_end)
/*
compute gravity at r_grid[grid_end].
density are given at rho[] and gravity at r_grid[grid_start]
is g0.

g0 should be normalized value= (g0 real_value)/(UG*scale_rho*scale_rad)

to get real gravity, UG*scale_rad*scale_rho should be
multiplied.
*/
double g0;
double rho[],r_grid[];
int grid_start,grid_end;
{
	double tmp,r1,r2;

	r1 = r_grid[grid_start];
	r2 = r_grid[grid_end];

	tmp = 
	tmp=g0*r1*r1 + mass_r(rho,r_grid,grid_start,grid_end);

	if( r2 <=0.0 ) tmp = 0.0;
	else tmp/=(r2*r2);

	return tmp;
}

double find_normalization_from_uv(u,v,l,rho,r_grid,n_grid)
double *u,*v,*r_grid;
double *rho;
int l,n_grid;
{
	double total,r1,r2,rho1,rho2,ll;
	int i,j;

	total=0.0;
	ll = l*(l+1.0);
	for(i=0;i<n_grid-1;i++)
	{
		j = i+1;
		r1=r_grid[i];
		r2=r_grid[j];
		if( r1==r2) continue;
/*
		total_mass += (rho[i]*r1*r1+rho[i+1]*r2*r2)*0.5*(r2-r1);
*/
/*
		rho1=rho[i]*(u[i]*u[i] + ll*v[i]*v[i]);
		rho2=rho[j]*(u[j]*u[j] + ll*v[j]*v[j]);
		total +=(r2*rho1-r1*rho2)*(r2*r2+r1*r1+r1*r2)/3.0
				+(rho2-rho1)*(r2*r2*r2+r2*r2*r1+r2*r1*r1+r1*r1*r1)/4.0;
*/
		rho1=rho[i]*(u[i]*u[i] + ll*v[i]*v[i])*r1*r1;
		rho2=rho[j]*(u[j]*u[j] + ll*v[j]*v[j])*r2*r2;
		total += (rho1+rho2)*0.5*(r2-r1);
	}
	return total;
}

double find_normalization_from_x(x1,x3,rho,r_grid,n_grid)
double *x1,*x3,*r_grid;
double *rho;
int n_grid;
{
	double total,r1,r2,rho1,rho2;
	int i,j;

	total=0.0;
	for(i=0;i<n_grid-1;i++)
	{
		j = i+1;
		r1=r_grid[i];
		r2=r_grid[j];

		if( r1==r2) continue;

		rho1=rho[i]*(x1[i]*x1[i] + x3[i]*x3[i]);
		rho2=rho[j]*(x1[j]*x1[j] + x3[j]*x3[j]);

		total += (rho1+rho2)*0.5*(r2-r1);
	}
	return total;
}

double integral_I1(x1,x3,rho,r_grid,grid_start,grid_end)
double *x1,*x3,*r_grid;
double *rho;
int grid_start,grid_end;
{
	int i,j;
	double r1,r2,total,rho1,rho2;

	total=0.0;
	for(i=grid_start;i<grid_end;i++)
	{
		j = i+1;
		r1 = r_grid[i];
		r2 = r_grid[j];
		if( r1==r2 ) continue;
		rho1=rho[i]*(x1[i]*x1[i] + x3[i]*x3[i]);
		rho2=rho[j]*(x1[j]*x1[j] + x3[j]*x3[j]);
		total += (rho1+rho2)*0.5*(r2-r1);
	}
	return total;
}

double integral_I1_fluid(x1,rho,r_grid,grid_start,grid_end)
double *x1,*r_grid;
double *rho;
int grid_start,grid_end;
{
	int i,j;
	double r1,r2,total,rho1,rho2;

	total=0.0;
	for(i=grid_start;i<grid_end;i++)
	{
		j = i+1;
		r1 = r_grid[i];
		r2 = r_grid[j];
		if( r1==r2 ) continue;
		rho1=rho[i]*(x1[i]*x1[i]);
		rho2=rho[j]*(x1[j]*x1[j]);
		total += (rho1+rho2)*0.5*(r2-r1);
	}
	return total;
}

double integral_I2(l,x1,x2,x3,x4,x5,x6,rho,g,a,c,el,en,f,r_grid,grid_start,grid_end)
double *x1,*x2,*x3,*x4,*x5,*x6,*r_grid;
double *rho,*g,*a,*c,*el,*en,*f;
int l,grid_start,grid_end;
/*
## NOTE ##
x3[] and x4[] are should be set correctly.
## NOTE END ##
*/
{
	int i,j;
	double r1,r2,total,rho1,rho2,ll,zz1,zz2,sqrt_ll;

	total=0.0;

	ll = l*(l+1.0);
	sqrt_ll = sqrt(ll);

	if( x5 == (double *)NULL || x6 == (double *)NULL )
	/* no potential eigenfunction */
	{
		for(i=grid_start;i<grid_end;i++)
		{
			j = i+1;
			r1 = r_grid[i];
			r2 = r_grid[j];
			if( r1==r2 ) continue;
			zz1 = 2.0*x1[i]-sqrt_ll*x3[i];
			zz2 = 2.0*x1[j]-sqrt_ll*x3[j];
			rho1= r1==0.0 ? 0.0 :
			     x2[i]*x2[i]/c[i] -
			     2.0*rho[i]*g[i]*x1[i]*zz1/r1;
			rho2=x2[j]*x2[j]/c[j] - 
			     2.0*rho[j]*g[j]*x1[j]*zz2/r2;
		
			total += (rho1+rho2)*0.5*(r2-r1);
			/* NOT fluid region */
			if( x4 != NULL )
			if( el[i] !=0.0 )
			{
			rho1= r1==0.0 ? 0.0 :
			     x4[i]*x4[i]/el[i] +
			     (a[i]-en[i]-f[i]*f[i]/c[i])*zz1*zz1/(r1*r1) +
			     (l-1.0)*(l+2.0)*en[i]*x3[i]*x3[i]/(r1*r1);
			rho2=x4[j]*x4[j]/el[j] +
			     (a[j]-en[j]-f[j]*f[j]/c[j])*zz2*zz2/(r2*r2) + 
			     (l-1.0)*(l+2.0)*en[j]*x3[j]*x3[j]/(r2*r2); 

			total += (rho1+rho2)*0.5*(r2-r1);
			}
		}
		/* NOT fluid region */
/*
		if( x4 != (double *)NULL )
		for(i=grid_start;i<grid_end;i++)
		{
			j = i+1;
			r1 = r_grid[i];
			r2 = r_grid[j];
			if( r1==r2 ) continue;
			zz1 = 2.0*x1[i]-sqrt_ll*x3[i];
			zz2 = 2.0*x1[j]-sqrt_ll*x3[j];
			rho1= r1==0.0 ? 0.0 :
			     x4[i]*x4[i]/el[i] +
			     (a[i]-en[i]-f[i]*f[i]/c[i])*zz1*zz1/(r1*r1) +
			     (l-1.0)*(l+2.0)*en[i]*x3[i]*x3[i]/(r1*r1);
			rho2=x4[j]*x4[j]/el[j] +
			     (a[j]-en[j]-f[j]*f[j]/c[j])*zz2*zz2/(r2*r2) + 
			     (l-1.0)*(l+2.0)*en[j]*x3[j]*x3[j]/(r2*r2); 

			total += (rho1+rho2)*0.5*(r2-r1);
		}
*/
	}
	else/* with potential perturbation eigenfunctions */
	{
		/* x3 is computed from x1,x2 and x5 */
		for(i=grid_start;i<grid_end;i++)
		{
			j = i+1;
			r1 = r_grid[i];
			r2 = r_grid[j];
			if( r1==r2 ) continue;
			zz1 = 2.0*x1[i]-sqrt_ll*x3[i];
			zz2 = 2.0*x1[j]-sqrt_ll*x3[j];
			rho1= r1==0.0 ? 0.0 :
			     x2[i]*x2[i]/c[i] + 
			     x6[i]*x6[i] -
			     2.0*rho[i]*g[i]*x1[i]*zz1/r1 -
			     2.0*(l+1.0)*rho[i]*x1[i]*x5[i]/r1 +
			     2.0*sqrt_ll*rho[i]*x3[i]*x5[i]/r1;
			rho2=x2[j]*x2[j]/c[j] + 
			     x6[j]*x6[j] -
			     2.0*rho[j]*g[j]*x1[j]*zz2/r2 -
			     2.0*(l+1.0)*rho[j]*x1[j]*x5[j]/r2 +
			     2.0*sqrt_ll*rho[j]*x3[j]*x5[j]/r2;
		
			total += (rho1+rho2)*0.5*(r2-r1);
			if( x4 != NULL )
			if( el[i] != 0.0 )/* NOT fluid region */
			{
			rho1= r1==0.0 ? 0.0 :
			     x4[i]*x4[i]/el[i] +
			     (a[i]-en[i]-f[i]*f[i]/c[i])*zz1*zz1/(r1*r1) +
			     (l-1.0)*(l+2.0)*en[i]*x3[i]*x3[i]/(r1*r1);
			rho2=x4[j]*x4[j]/el[j] +
			     (a[j]-en[j]-f[j]*f[j]/c[j])*zz2*zz2/(r2*r2) + 
			     (l-1.0)*(l+2.0)*en[j]*x3[j]*x3[j]/(r2*r2); 

			total += (rho1+rho2)*0.5*(r2-r1);
			}
		}
		/* NOT fluid region */
/*
		if( x4 != NULL )
		for(i=grid_start;i<grid_end;i++)
		{
			j = i+1;
			r1 = r_grid[i];
			r2 = r_grid[j];
			if( r1==r2 ) continue;
			zz1 = 2.0*x1[i]-sqrt_ll*x3[i];
			zz2 = 2.0*x1[j]-sqrt_ll*x3[j];
			rho1= r1==0.0 ? 0.0 :
			     x4[i]*x4[i]/el[i] +
			     (a[i]-en[i]-f[i]*f[i]/c[i])*zz1*zz1/(r1*r1) +
			     (l-1.0)*(l+2.0)*en[i]*x3[i]*x3[i]/(r1*r1);
			rho2=x4[j]*x4[j]/el[j] +
			     (a[j]-en[j]-f[j]*f[j]/c[j])*zz2*zz2/(r2*r2) + 
			     (l-1.0)*(l+2.0)*en[j]*x3[j]*x3[j]/(r2*r2); 

			total += (rho1+rho2)*0.5*(r2-r1);
		}
*/
	}
	return total;
}

double integral_I2_fluid(l,omg2,x1,x2,x5,x6,rho,g,kappa,r_grid,grid_start,grid_end)
double *x1,*x2,*x5,*x6,*r_grid,omg2;
double *rho,*g,*kappa;
int l,grid_start,grid_end;
{
	int i,j;
	double r1,r2,total,rho1,rho2,ll,ll_r2omg2_1,ll_r2omg2_2;

	total=0.0;

	ll = l*(l+1.0);
	if( x5 == NULL || x6 == NULL) /* no potential eigenfunction */
	for(i=grid_start;i<grid_end;i++)
	{
		j = i+1;
		r1 = r_grid[i];
		r2 = r_grid[j];
		ll_r2omg2_1 = r1==0.0 ? 0.0 : ll/(r1*r1*omg2);
		ll_r2omg2_2 = ll/(r2*r2*omg2);
		if( r1==r2 ) continue;
		rho1= r1==0.0 ? 0.0 :
		     (1.0/kappa[i]-ll_r2omg2_1/rho[i])*x2[i]*x2[i]+
		     rho[i]*g[i]*(ll_r2omg2_1*g[i]-4.0/r1)*x1[i]*x1[i];
		rho2=(1.0/kappa[j]-ll_r2omg2_2/rho[j])*x2[j]*x2[j]+
		     rho[j]*g[j]*(ll_r2omg2_2*g[j]-4.0/r2)*x1[j]*x1[j];
	
		total += (rho1+rho2)*0.5*(r2-r1);
	}
	else
	for(i=grid_start;i<grid_end;i++)
	{
		j = i+1;
		r1 = r_grid[i];
		r2 = r_grid[j];
		ll_r2omg2_1 = ll/(r1*r1*omg2);
		ll_r2omg2_2 = ll/(r2*r2*omg2);
		if( r1==r2 ) continue;
		rho1= r1==0.0 ? 0.0 : 
		     (1.0/kappa[i]-ll_r2omg2_1/rho[i])*x2[i]*x2[i]+
		     x6[i]*x6[i]+
		     rho[i]*g[i]*(ll_r2omg2_1*g[i]-4.0/r1)*x1[i]*x1[i]+
		     ll_r2omg2_1*rho[i]*x5[i]*x5[i]+
		     2.0*(ll_r2omg2_1*g[i]-(l+1.0)/r1)*rho[i]*x1[i]*x5[i];
		rho2=(1.0/kappa[j]-ll_r2omg2_2*rho[j])*x2[j]*x2[j]+
		     x6[j]*x6[j]+
		     rho[j]*g[j]*(ll_r2omg2_2*g[j]-4.0/r2)*x1[j]*x1[j]+
		     ll_r2omg2_2*rho[j]*x5[j]*x5[j]+
		     2.0*(ll_r2omg2_2*g[j]-(l+1.0)/r2)*rho[j]*x1[j]*x5[j];
	
		total += (rho1+rho2)*0.5*(r2-r1);
	}
	return total;
}

double integral_I2radial(x1,x2,x5,x6,rho,g,a,c,en,f,r_grid,grid_start,grid_end)
double *x1,*x2,*x5,*x6,*r_grid;
double *rho,*g,*a,*c,*en,*f;
int grid_start,grid_end;
{
	int i,j;
	double r1,r2,total,rho1,rho2;

	total=0.0;

	if( x5 == NULL || x6 == NULL) /* no potential eigenfunction */
	for(i=grid_start;i<grid_end;i++)
	{
		j = i+1;
		r1 = r_grid[i];
		r2 = r_grid[j];
		if( r1==r2 ) continue;
		rho1= r1==0.0 ? 0.0 :
		     x2[i]*x2[i]/c[i]+
		     4.0*((a[i]-f[i]*f[i]/c[i]-en[i])/r1-rho[i]*g[i])*x1[i]*x1[i]/r1;
		rho2=x2[j]*x2[j]/c[j]+
		     4.0*((a[j]-f[j]*f[j]/c[j]-en[j])/r2-rho[j]*g[j])*x1[j]*x1[j]/r2;
	
		total += (rho1+rho2)*0.5*(r2-r1);
	}
	else
	for(i=grid_start;i<grid_end;i++)
	{
		j = i+1;
		r1 = r_grid[i];
		r2 = r_grid[j];
		if( r1==r2 ) continue;
		rho1= r1==0.0 ? 0.0 :
		     x2[i]*x2[i]/c[i]+
		     4.0*((a[i]-f[i]*f[i]/c[i]-en[i])/r1-rho[i]*g[i])*x1[i]*x1[i]/r1+
		     x6[i]*x6[i] - 2.0*rho[i]*x1[i]*x5[i]/r1;

		rho2=x2[j]*x2[j]/c[j]+
		     4.0*((a[j]-f[j]*f[j]/c[j]-en[j])/r2-rho[j]*g[j])*x1[j]*x1[j]/r2+
		     x6[j]*x6[j] - 2.0*rho[j]*x1[j]*x5[j]/r2;

	
		total += (rho1+rho2)*0.5*(r2-r1);
	}
	return total;
}
double integral_I3(l,x1,x2,x3,x4,x5,x6,rho,g,a,c,el,en,f,r_grid,grid_start,grid_end)
double *x1,*x2,*x3,*x4,*x5,*x6,*r_grid;
double *rho,*g,*a,*c,*el,*en,*f;
int l,grid_start,grid_end;
/*
## NOTE ##
x3[] and x4[] are should be set correctly.
## NOTE END ##
*/
{
	int i,j;
	double r1,r2,total,rho1,rho2,ll,sqrt_ll,af2c1,af2c2;

	total=0.0;

	ll = l*(l+1.0);
	sqrt_ll = sqrt(ll);
	if( x5 == NULL || x6 == NULL)
	/* no potential eigenfunction */
	{
		for(i=grid_start;i<grid_end;i++)
		{
			j = i+1;
			r1 = r_grid[i];
			r2 = r_grid[j];
			if( r1==r2 ) continue;
			rho1 = r1==0.0 ? 0.0 :
				( -f[i]*x2[i]/c[i]
				  +rho[i]*g[i]*x1[i])*x3[i]/(sqrt_ll*r1);
			rho2 =( -f[j]*x2[j]/c[j]
				  +rho[j]*g[j]*x1[j])*x3[j]/(sqrt_ll*r2);
			total += (rho1+rho2)*0.5*(r2-r1);
			/* NOT fluid region */
			if( x4 != NULL )
			if( el[i] !=0.0 )
			{
			af2c1=a[i]-f[i]*f[i]/c[i];
			af2c2=a[j]-f[j]*f[j]/c[j];
			rho1 = r1==0.0 ? 0.0 :
				 ( af2c1*x3[i]*x3[i]/r1
				  +(x4[i]
				   -2.0*(af2c1-en[i])*x3[i]/r1)*x1[i]/sqrt_ll)/r1;
			rho2 = ( af2c2*x3[j]*x3[j]/r2
				  +(x4[j]
				   -2.0*(af2c2-en[j])*x3[j]/r2)*x1[j]/sqrt_ll)/r2;
			total += (rho1+rho2)*0.5*(r2-r1);
			}
		}
		/* NOT fluid region */
/*
		if( x4 != (double *)NULL )
		for(i=grid_start;i<grid_end;i++)
		{
			j = i+1;
			r1 = r_grid[i];
			r2 = r_grid[j];
			if( r1==r2 ) continue;
			af2c1=a[i]-f[i]*f[i]/c[i];
			af2c2=a[j]-f[j]*f[j]/c[j];
			rho1 = r1==0.0 ? 0.0 : 
			       ( af2c1*x3[i]*x3[i]/r1
				  +(x4[i]
				   -2.0*(af2c1-en[i])*x3[i]/r1)*x1[i]/sqrt_ll)/r1;
			rho2 = ( af2c2*x3[j]*x3[j]/r2
				  +(x4[j]
				   -2.0*(af2c2-en[j])*x3[j]/r2)*x1[j]/sqrt_ll)/r2;
			total += (rho1+rho2)*0.5*(r2-r1);
		}
*/
	}
	else/* with potential perturbation eigenfunctions */
	{
		/* x3 is computed from x1,x2 and x5 */
		for(i=grid_start;i<grid_end;i++)
		{
			j = i+1;
			r1 = r_grid[i];
			r2 = r_grid[j];
			if( r1==r2 ) continue;
			rho1 = r1==0.0 ? 0.0 : 
			       x5[i]*(x6[i]-rho[i]*x1[i])/(r1*(l+0.5))
				+(-f[i]/c[i]*x2[i]+rho[i]*(g[i]*x1[i]+x5[i]))*x3[i]
				/(r1*sqrt_ll);
			rho2 = x5[j]*(x6[j]-rho[j]*x1[j])/(r2*(l+0.5))
				+(-f[j]/c[j]*x2[j]+rho[j]*(g[j]*x1[j]+x5[j]))*x3[j]
				/(r2*sqrt_ll);
		
			total += (rho1+rho2)*0.5*(r2-r1);
			if( x4 != NULL )
			if( el[i] != 0.0 ) /* NOT fluid region */
			{
			af2c1=a[i]-f[i]*f[i]/c[i];
			af2c2=a[j]-f[j]*f[j]/c[j];
			rho1 = r1==0.0 ? 0.0 :
				 ( af2c1*x3[i]*x3[i]/r1
				  +(x4[i]
				   -2.0*(af2c1-en[i])*x3[i]/r1)*x1[i]/sqrt_ll)/r1;
			rho2 = ( af2c2*x3[j]*x3[j]/r2
				  +(x4[j]
				   -2.0*(af2c2-en[j])*x3[j]/r2)*x1[j]/sqrt_ll)/r2;
			total += (rho1+rho2)*0.5*(r2-r1);
			}
		}
		/* NOT fluid region */
/*
		if( x4 != (double *)NULL )
		for(i=grid_start;i<grid_end;i++)
		{
			j = i+1;
			r1 = r_grid[i];
			r2 = r_grid[j];
			if( r1==r2 ) continue;
			af2c1=a[i]-f[i]*f[i]/c[i];
			af2c2=a[j]-f[j]*f[j]/c[j];
			rho1 = r1==0.0 ? 0.0 :
				 ( af2c1*x3[i]*x3[i]/r1
				  +(x4[i]
				   -2.0*(af2c1-en[i])*x3[i]/r1)*x1[i]/sqrt_ll)/r1;
			rho2 = ( af2c2*x3[j]*x3[j]/r2
				  +(x4[j]
				   -2.0*(af2c2-en[j])*x3[j]/r2)*x1[j]/sqrt_ll)/r2;
			total += (rho1+rho2)*0.5*(r2-r1);
		}
*/
	}
	return total;
}
