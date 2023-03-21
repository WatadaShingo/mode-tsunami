#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include "m_solid.h" /* UPPER_REGION and LOWER_REGION are defined */
#include "eigen.h"

#define FR	fprintf(stderr
#define YES		1
#define NO		0

void normalize_eigen(u,du,r_grid,n_grid,a0,as,rho0,rhos,g0,i1,is_cowl)
double **u,**du,*r_grid,i1;
double a0,as,rho0,rhos,g0;
int n_grid,is_cowl;
/*
 * normalize eigenfunction to store
<< OUTPUT >>
u[3][],	scaled eigenfunc U,V,FI
du[3][],	and their derivative
<< INPUT >>
u[3][],	computed eigenfunctions with normalization a0,rho,g0
du[3][],	and their derivative
r_grid[],	radius of grid point of u[][] and du[][], from center in unit of a0.
n_grid,	number of grid points
a0,		radius normalization used in the computation of eigenfunctions
as,		radius normalization used to store the eigenfunctions
rho0,		density normalization used in the computation of eigenfunctions
rhos,		density normalization used to store the eigenfunctions
g0,		gravity normalization used in the computation of eigenfunctions
i1,		integral of I1=rho*(U^2+(l+1)*l*V^2)*r^2 using computed
		eigenfunctions
a0,as,rho0 ans rhos should have the same unit system. (eg. MKS, cgs...)
is_cowl,    if Cowling approximation is used ==YES, if not == NO
 */
{
	int i;
	double norm;

	norm = sqrt(i1*rho0*a0*a0*a0/(rhos*as*as*as));
	for(i=0;i<n_grid;i++)
	{
		u[0][i] /= norm; /* U */
		u[1][i] /= norm; /* V */
		du[0][i] *=(as/(a0*norm)); /* dU/dr */
		du[1][i] *=(as/(a0*norm)); /* dV/dr */
		if( is_cowl == NO )
		{	 u[2][i] /= norm; /* PHI */
			du[2][i] *=(as/(a0*norm)); /* dPHI/dr */

	/* ## NOTE ##
	Up to here
	Laplacian( potential perturbation ) = - div (rho displacement) 
	should hold when scaling of earth model parameters are derived
	from a0, rho0 and (Gravitational constant).
	If so, following scaling to the potential perturbation will result in
	the same convention used in the B.sprem eigenfunction file.
	See note for details.
	## NOTE END ##  Aug. 26, 1994 */
		u[2][i] *=g0/as;
		du[2][i] *=g0/as;
		}
	}
}
void re_normalize_eigen(u,du,r_grid,n_grid,a0,as,rho0,rhos,g0,is_cowl)
double **u,**du,*r_grid;
double a0,as,rho0,rhos,g0;
int n_grid,is_cowl;
/*
 * re-normalize eigenfunction for computation 
<< OUTPUT >>
u[3][],	comuted eigenfunc U,V,FI with normalization a0,rho0,g0
du[3][],	and their derivative
## MEMO ##
eigenfunctions are normalized so that
|r=surface
|rho(r)*(U(r)^2+l(l+1)*V^2) r^2 dr = 1.0,
|r=0
where rho,r are normalize by a0 and rho0
## MEMO END ##
<< INPUT >>
u[3][],	stored eigenfunctions with normalization as,rhs,gs
du[3][],	and their derivative
r_grid[],	radius of grid point of u[][] and du[][], from center in unit of a0.
n_grid,	number of grid points
a0,		radius normalization used in the computation of eigenfunctions
as,		radius normalization used to store the eigenfunctions
rho0,		density normalization used in the computation of eigenfunctions
rhos,		density normalization used to store the eigenfunctions
g0,		gravity normalization used in the computation of eigenfunctions
a0,as,rho0 ans rhos should have the same unit system. (eg. MKS, cgs...)
is_cowl,    if Cowling approximation is used ==YES, if not == NO
*/
{
	int i;
	double norm;

	norm = sqrt(rho0*a0*a0*a0/(rhos*as*as*as));
	for(i=0;i<n_grid;i++)
	{
		u[0][i] *= norm; /* U */
		u[1][i] *= norm; /* V */
		du[0][i] *=(a0*norm/as); /* dU/dr */
		du[1][i] *=(a0*norm/as); /* dV/dr */
		if( is_cowl == NO )
		{	 u[2][i] *= norm; /* PHI */
			du[2][i] *=(a0*norm/as); /* dPHI/dr */

	/* ## NOTE ##
	Up to here
	Laplacian( potential perturbation ) = - div (rho displacement) 
	should hold when scaling of earth model parameters are derived
	from a0, rho0 and (Gravitational constant).
	If so, following scaling to the potential perturbation will result in
	the same convention used in the B.sprem eigenfunction file.
	See note for details.
	## NOTE END ##  Aug. 26, 1994 */
		u[2][i] *=as/g0;
		du[2][i] *=as/g0;
		}
	}
}

void print_eigen(file_prefix,t_grid,n_grid,u,du,rho,a0,l,n_solution,freq,is_cowl)
char *file_prefix;
double *t_grid,**u,**du,*rho;
double a0,freq;
int n_grid,l,n_solution,is_cowl;
{
	int i;
	char filename[100];
	double ll,r;
	FILE *fp;

	ll = l*(1.0+l);
	sprintf(filename,"%s.%ds%d",file_prefix,n_solution,l);
	if((fp=fopen(filename,"w"))==NULL){FR,"print_eigen:(%s) %d: cannot create %s.\n",__FILE__,__LINE__,filename);exit(-1);}

	fprintf(fp,"%d s %d %.2f\n",n_solution,l,1.0/freq);
	fprintf(fp,"%d 7\n",n_grid);
	fprintf(fp,"njunk\n");

	if( is_cowl == YES )
	{
		if( rho == (double *)NULL )
		for(i=0;i<n_grid;i++)
		{
		r = t_grid[i];
		fprintf(fp,"%+9.3e %+9.3e %+9.3e %+9.3e %+9.3e\n",
		r*a0/1e3,u[0][i],du[0][i],u[1][i],du[1][i]);
		}
		else /* if *rho is NOT null */
		for(i=0;i<n_grid;i++)
		{
		r = t_grid[i];
		fprintf(fp,"%+9.3e %+9.3e %+9.3e %+9.3e %+9.3e %+9.3e\n",
		r*a0/1e3,u[0][i],du[0][i],u[1][i],du[1][i],
		rho[i]*(u[0][i]*u[0][i]+ll*u[1][i]*u[1][i])*r*r);
		}
	}
	else
	{
		if( rho == (double *)NULL )
		for(i=0;i<n_grid;i++)
		{
		r = t_grid[i];
		fprintf(fp,"%+9.3e %+9.3e %+9.3e %+9.3e %+9.3e %+9.3e %+9.3e\n",
		r*a0/1e3,u[0][i],du[0][i],u[1][i],du[1][i],u[2][i],du[2][i]);
		}
		else /* if *rho is NOT null */
		for(i=0;i<n_grid;i++)
		{
		r = t_grid[i];
		fprintf(fp,"%+9.3e %+9.3e %+9.3e %+9.3e %+9.3e %+9.3e %+9.3e %9.3e\n",
		r*a0/1e3,u[0][i],du[0][i],u[1][i],du[1][i],u[2][i],du[2][i],
		rho[i]*(u[0][i]*u[0][i]+ll*u[1][i]*u[1][i])*r*r);
		}
	}

	fclose(fp);
}

void write_eigen(fd,omg,l,n,u,du,n_grid)
int fd,l,n,n_grid;
double omg,**u,**du;
{
	int i;
	struct eigen_mode mode;

	if(fd<0){ FR,"write_eigen(%s) %d: input error.\n",
		__FILE__,__LINE__);exit(-1);}

	mode.omg = (float)omg;
	mode.l = l;
	mode.n = n;

	for(i=0;i<n_grid;i++)
	{
		mode.u[i] =(float) u[0][i];
		mode.v[i] =(float) u[1][i];
		mode.p[i] =(float) u[2][i];
		mode.du[i] =(float) du[0][i];
		mode.dv[i] =(float) du[1][i];
		mode.dp[i] =(float) du[2][i];
	}
		write(fd,&mode,sizeof(struct eigen_mode));
}
void write_eigencowl(fd,omg,l,n,u,du,n_grid)
int fd,l,n,n_grid;
double omg,**u,**du;
{
	int i;
	struct eigen_modecowl mode;

	if(fd<0){ FR,"write_eigencowl(%s) %d: input error.\n",
		__FILE__,__LINE__);exit(-1);}

	mode.omg = (float)omg;
	mode.l = l;
	mode.n = n;

	for(i=0;i<n_grid;i++)
	{
		mode.u[i] =(float) u[0][i];
		mode.v[i] =(float) u[1][i];
		mode.du[i] =(float) du[0][i];
		mode.dv[i] =(float) du[1][i];
	}
		write(fd,&mode,sizeof(struct eigen_modecowl));
}

int read_eigen(fd,mode)
int fd;
struct eigen_mode *mode;
{
	if(fd<0 || mode == NULL)
	{ FR,"read_eigen(%s) %d: input error.\n",
		__FILE__,__LINE__);exit(-1);}

		return (read(fd, mode,sizeof(struct eigen_mode)));
}

int read_eigencowl(fd,mode)
int fd;
struct eigen_modecowl *mode;
{
	if(fd<0 || mode == NULL)
	{ FR,"read_eigencowl(%s) %d: input error.\n",
		__FILE__,__LINE__);exit(-1);}

		return (read(fd, mode,sizeof(struct eigen_modecowl)));
}

void get_layer(layer,r,r_grid,n_grid,type)
/*
 * set layer which contains r
<< INPUT >>
r: radius
r_grid[]:  grid_point values of r
n_grid: number of grid points of r_grid[]
type: UPPER_REGION or LOWER_REGION
<< OUTPUT >>
layer: which contains (bottom of layer) <=r< (top of layer)

grid # R
6     r[6]
5     r[5]
3,4   r[3]==r[4]
2     r[2]
1     r[1]
0     r[0]
when r=r[3]=r[4]
if type == UPPER_REGION, layer = 4
if type == LOWER_REGION, layer = 2

*/
int *layer;
double r,*r_grid;
int n_grid,type;
{
	int i;
	if(r < r_grid[0] || r >r_grid[n_grid-1])
	{FR,"get_layer(%s) %d: input error r=%e.\n",
	__FILE__,__LINE__,r);exit(-1);}

	for(i=0;r>=r_grid[i] && i<n_grid;i++);

	*layer = i-1;

	/* if layer is near the boundary */
	if (r_grid[i] == r_grid[i-1] )
	{
	switch( type ){
	case UPPER_REGION:
		FR,"get_layer(%s) %d : cannot find upper layer of r=%e.\n",
		__FILE__,__LINE__,r);exit(-1);
		break;
	case LOWER_REGION:
		*layer = i-2;
		break;
		}
	}

	if( r_grid[i-1] == r_grid[i-2] )
	{
	switch( type ){
	case UPPER_REGION:
		*layer = i-1; 
		break;
	case LOWER_REGION:
		if( i <=2 ){
		FR,"get_layer(%s) %d: cannot find lower layer of r=%e.\n",
		__FILE__,__LINE__,r);exit(-1);}
		*layer = i-3;
		break;
		}
	}
}

void get_eigen_at_r(r,s0,s1,s)
/*
 * get eigenfunction 's' from s0,s1
 */
struct eigen_func *s,s0,s1;
double r;
{
	
	float fdum;
	s->i = s0.i;
	s->r = r;
	quadra_interp(s0.u,
			s0.du,
			s1.u,
			s1.du,
			(float)(s1.r-s0.r),
			(float)(r-s0.r),
			&(s->u),
			&(s->du),
			&fdum
			);
	quadra_interp(s0.v,
			s0.dv,
			s1.v,
			s1.dv,
			(float)(s1.r-s0.r),
			(float)(r-s0.r),
			&(s->v),
			&(s->dv),
			&fdum
			);
	quadra_interp(s0.p,
			s0.dp,
			s1.p,
			s1.dp,
			(float)(s1.r-s0.r),
			(float)(r-s0.r),
			&(s->p),
			&(s->dp),
			&fdum
			);
}
void get_eigen_at_rcowl(r,s0,s1,s)
/*
 * get eigenfunction 's' from s0,s1
 */
struct eigen_func *s,s0,s1;
double r;
{
	
	float fdum;
	s->i = s0.i;
	s->r = r;
	quadra_interp(s0.u,
			s0.du,
			s1.u,
			s1.du,
			(float)(s1.r-s0.r),
			(float)(r-s0.r),
			&(s->u),
			&(s->du),
			&fdum
			);
	quadra_interp(s0.v,
			s0.dv,
			s1.v,
			s1.dv,
			(float)(s1.r-s0.r),
			(float)(r-s0.r),
			&(s->v),
			&(s->dv),
			&fdum
			);
	s->p = s->dp = 0.0;
}
