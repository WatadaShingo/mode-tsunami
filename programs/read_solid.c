#include <stdio.h>
#include <string.h>
#include <math.h>
#include "earthconst.h"
#include "m_solid.h"
#include "read_solid.h"
#include "model.h"
#define FR		fprintf(stderr
#define PI		3.14159265358979324
#define ABS(x)		((x)>0.0 ? (x) : -(x))

int layer;

struct model_solid msolid;

#define LINEAR1(param,value) \
	value = (double)(param[layer+1] - param[layer])*hx + param[layer]

#define NATURAL         1.0e30
int read_solid(modelfilename)
	char *modelfilename;
{
	int i,j,region,start_layer;
	FILE *fp;
	char line[200];
	struct model_funcsolid model_unit;
	double tmp_array[SAMPLE_SOLID],dum,dum1;

	strcpy(msolid.modelfilename,modelfilename);
	if((fp = fopen( msolid.modelfilename,"r"))==NULL)
	{
		FR,"read_solid(%s) %d: cannot find file %s.\n",
		__FILE__,__LINE__,msolid.modelfilename);exit(-1);
	}

	fgets(line,200,fp);
	fscanf(fp,"%s %d",msolid.modelname,&(msolid.nmodel));
	fscanf(fp,"%*lf %*lf %lf %lf %lf %lf %lf %lf %*lf %*lf %*lf %lf %lf %lf %lf %lf",
			&(model_unit.r),
			&(model_unit.rho),
			&(model_unit.vpv),
			&(model_unit.vph),
			&(model_unit.vsv),
			&(model_unit.vsh),
			&(model_unit.a),
			&(model_unit.c),
			&(model_unit.l),
			&(model_unit.n),
			&(model_unit.f));
	for(i=0;i<msolid.nmodel;i++)
	{
	fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		&(msolid.value[REGOIN_SOLID][i]),
		&(msolid.value[LAYER_SOLID][i]),
		&(msolid.value[RAD_SOLID][i]),
		&(msolid.value[RHO_SOLID][i]),
		&(msolid.value[VPV_SOLID][i]),
		&(msolid.value[VPH_SOLID][i]),
		&(msolid.value[VSV_SOLID][i]),
		&(msolid.value[VSH_SOLID][i]),
		&(msolid.value[ETA_SOLID][i]),
		&(msolid.value[QU_SOLID][i]),
		&(msolid.value[QK_SOLID][i]),
		&(msolid.value[A_SOLID][i]),
		&(msolid.value[C_SOLID][i]),
		&(msolid.value[L_SOLID][i]),
		&(msolid.value[N_SOLID][i]),
		&(msolid.value[F_SOLID][i]));

		msolid.value[RAD_SOLID][i]*=model_unit.r;
		msolid.value[RHO_SOLID][i]*=model_unit.rho;
		msolid.value[VPV_SOLID][i]*=model_unit.vpv;
		msolid.value[VPH_SOLID][i]*=model_unit.vph;
		msolid.value[VSV_SOLID][i]*=model_unit.vsv;
		msolid.value[VSH_SOLID][i]*=model_unit.vsh;
		msolid.value[A_SOLID][i]*=model_unit.a;
		msolid.value[C_SOLID][i]*=model_unit.c;
		msolid.value[L_SOLID][i]*=model_unit.l;
		msolid.value[N_SOLID][i]*=model_unit.n;
		msolid.value[F_SOLID][i]*=model_unit.f;
		msolid.value[QK_SOLID][i] = 1.0/msolid.value[QK_SOLID][i];
		msolid.value[QU_SOLID][i] = ( msolid.value[QU_SOLID][i] < 0.0 ?
		0.0 : 1.0/msolid.value[QU_SOLID][i]);
	}
	fclose(fp);
/*
 * set gravity
 */
	for(i=0;i<msolid.nmodel;i++)
	{
	msolid.value[G_SOLID][i] = grav_r(0.0,msolid.value[RHO_SOLID],
			msolid.value[RAD_SOLID],0,i);
	msolid.value[G_SOLID][i] *= UG;
	}

/*
 * set up region boundary layer
 */
	region = 0;
	msolid.region[region][LOWER_BOUND_SOLID]=0;
	for(i=0;i<msolid.nmodel-1;i++)
	if( msolid.value[RAD_SOLID][i] == msolid.value[RAD_SOLID][i+1] )
	{
		msolid.region[region][UPPER_BOUND_SOLID]=i;
		region++;
		msolid.region[region][LOWER_BOUND_SOLID]=i+1;
	}
	msolid.region[region][UPPER_BOUND_SOLID]=msolid.nmodel-1;
	msolid.nregion = region+1;

	for(i=0;i<msolid.nregion;i++)
		msolid.region[i][NUMBER_OF_GRIDS] = 
			msolid.region[i][UPPER_BOUND_SOLID]-
			msolid.region[i][LOWER_BOUND_SOLID]+1;
/*
 * set up dho/dr
 */
	for(i=0;i<msolid.nregion;i++)
	{
	start_layer = msolid.region[i][LOWER_BOUND_SOLID];
	splined(&(msolid.value[RAD_SOLID][start_layer]),
	        &(msolid.value[RHO_SOLID][start_layer]),
		  msolid.region[i][NUMBER_OF_GRIDS],
		NATURAL,NATURAL,tmp_array);
	for(j=0;j<msolid.region[i][NUMBER_OF_GRIDS];j++)
		splintd(&(msolid.value[RAD_SOLID][start_layer]),
			  &(msolid.value[RHO_SOLID][start_layer]),
		tmp_array,msolid.region[i][NUMBER_OF_GRIDS],
		msolid.value[RAD_SOLID][j+start_layer],
		&dum,&(msolid.value[DRHODR_SOLID][j+start_layer]),&dum1);
	}
	return 0;
}

void find_modelsolid(model,r1,region)
	double r1;
	struct model_funcsolid *model;
	int region;
/* ## NOTE ##
r1 may not be exactly within the region range.
e.g.
r1-(bottom of region) < 0.0 but abs(r1-(bottom of region)) is very small
r1-(top of region) > 0.0 but abs(r1-(top of region)) is very small
*/
{
	int upper_bound_layer,lower_bound_layer,number_of_grids;
	double *r,*rho,*a,*c,*l,*n,*f,*g,*drhodr;
	double hx;

	upper_bound_layer = msolid.region[region][UPPER_BOUND_SOLID];
	lower_bound_layer = msolid.region[region][LOWER_BOUND_SOLID];
	number_of_grids   = msolid.region[region][NUMBER_OF_GRIDS];
	model->r = r1;

	r= &(msolid.value[RAD_SOLID][lower_bound_layer]);
	if( ABS(r1-r[0]) < 1e-10*(r[1]-r[0]) )
	{
		layer = 0;
		model->rho = msolid.value[RHO_SOLID][lower_bound_layer];
		model->a = msolid.value[A_SOLID][lower_bound_layer];
		model->c = msolid.value[C_SOLID][lower_bound_layer];
		model->l = msolid.value[L_SOLID][lower_bound_layer];
		model->n = msolid.value[N_SOLID][lower_bound_layer];
		model->f = msolid.value[F_SOLID][lower_bound_layer];
		model->g = msolid.value[G_SOLID][lower_bound_layer];
		model->drhodr = msolid.value[DRHODR_SOLID][lower_bound_layer];
	}
	else if( ABS(r[number_of_grids-1]-r1) < 1e-10*(r[number_of_grids-1]-r[number_of_grids-2]))
	{
		layer = number_of_grids-1;
		model->rho = msolid.value[RHO_SOLID][upper_bound_layer];
		model->a = msolid.value[A_SOLID][upper_bound_layer];
		model->c = msolid.value[C_SOLID][upper_bound_layer];
		model->l = msolid.value[L_SOLID][upper_bound_layer];
		model->n = msolid.value[N_SOLID][upper_bound_layer];
		model->f = msolid.value[F_SOLID][upper_bound_layer];
		model->g = msolid.value[G_SOLID][upper_bound_layer];
		model->drhodr = msolid.value[DRHODR_SOLID][upper_bound_layer];
	}else{

	set_layersolid(r1,region);
	if( layer < 0 || layer > number_of_grids-2 )
	{FR,"find_modelsolid (%s) %d: layer is out of range r=%e\n",
		__FILE__,__LINE__,r1);exit(-1);}

	rho = &(msolid.value[RHO_SOLID][lower_bound_layer]);
	a = &(msolid.value[A_SOLID][lower_bound_layer]);
	c = &(msolid.value[C_SOLID][lower_bound_layer]);
	l = &(msolid.value[L_SOLID][lower_bound_layer]);
	n = &(msolid.value[N_SOLID][lower_bound_layer]);
	f = &(msolid.value[F_SOLID][lower_bound_layer]);
	g = &(msolid.value[G_SOLID][lower_bound_layer]);
	drhodr = &(msolid.value[DRHODR_SOLID][lower_bound_layer]);

	hx = (r1 - r[layer]) /( r[layer+1] -r[layer]);
	if( hx < 0.0 || hx > 1.0 )
		printf("hx=%f,r=%f\n",hx,r1);

	LINEAR1(rho,model->rho);
	LINEAR1(a,model->a);
	LINEAR1(c,model->c);
	LINEAR1(l,model->l);
	LINEAR1(n,model->n);
	LINEAR1(f,model->f);
	LINEAR1(g,model->g);
	LINEAR1(drhodr,model->drhodr);
	}
}

void find_modelsolidQ(model,r1,region)
	double r1;
	struct model_funcsolid *model;
	int region;
{
	int upper_bound_layer,lower_bound_layer,number_of_grids;
	double *r,*qk,*qu;
	double hx;

	upper_bound_layer = msolid.region[region][UPPER_BOUND_SOLID];
	lower_bound_layer = msolid.region[region][LOWER_BOUND_SOLID];
	number_of_grids   = msolid.region[region][NUMBER_OF_GRIDS];
	model->r = r1;

	r= &(msolid.value[RAD_SOLID][lower_bound_layer]);
	if( ABS(r1-r[0]) < 1e-10*(r[1]-r[0]) )
	{
		layer = 0;
		model->qk = msolid.value[QK_SOLID][lower_bound_layer];
		model->qu = msolid.value[QU_SOLID][lower_bound_layer];
	}
	else if( ABS(r[number_of_grids-1]-r1) < 1e-10*(r[number_of_grids-1]-r[number_of_grids-2]))
	{
		layer = number_of_grids-1;
		model->qk = msolid.value[QK_SOLID][upper_bound_layer];
		model->qu = msolid.value[QU_SOLID][upper_bound_layer];
	}else{

	set_layersolid(r1,region);
	if( layer < 0 || layer > number_of_grids-2 )
	{FR,"find_modelsolid (%s) %d: layer is out of range r=%e\n",
		__FILE__,__LINE__,r1);exit(-1);}

	qk = &(msolid.value[QK_SOLID][lower_bound_layer]);
	qu = &(msolid.value[QU_SOLID][lower_bound_layer]);

	hx = (r1 - r[layer]) /( r[layer+1] -r[layer]);
	if( hx < 0.0 || hx > 1.0 )
		printf("hx=%f,r=%f\n",hx,r1);

	LINEAR1(qk,model->qk);
	LINEAR1(qu,model->qu);
	}
}
void set_modelsolid_array(a,c,el,en,f,rho,g,qk,qu,r,region,layer_start,layer_end)
double *a,*c,*el,*en,*f,*rho,*g,*qk,*qu,*r;
int region,layer_start,layer_end;
{
	int j;
	struct model_funcsolid model;
      for(j=layer_start;j<=layer_end;j++)
	{
		find_modelsolid(&model,r[j],region);
		find_modelsolidQ(&model,r[j],region);
		rho[j] = model.rho;
		a[j]=model.a;
		c[j]=model.c;
		el[j]=model.l;
		en[j]=model.n;
		f[j]=model.f;
		g[j]=model.g;
		qk[j]=model.qk;
		qu[j]=model.qu;
	}
}

void normalize_modelsolid(a0,rho0,p0,vp0,g0,omg20)
/* normalize
struct model_solid msolid;
input a0    normalization of radius, in meter
     rho0   normalization of density, in kg/m^3
output
     p0    normalizatin of pressure & elastic constants, in Pa == 1e8 kbar
     vp0   normalization of seismic velocity, in m/sec
     g0    normalization of gravity, in m/sec^2
     omg20 normalization of omg^2 in sec^2

Here we assume three quantities,
universal graviry constant, rho and density, are given and compute
normalization factor of other parameters.
See note for details
*/
/* float a0,rho0,*p0,*vp0,*g0,*omg20; */
double a0,rho0,*p0,*vp0,*g0,*omg20;
{
	int i,j;

	*g0 = 4.0*PI*UG*rho0*a0;
	*p0 = *g0*rho0*a0;
	*vp0 = sqrt((double)*g0*a0);
	*omg20 = 4.0*PI*UG*rho0;
	for(i=0;i< msolid.nmodel;i++)
	{
		msolid.value[RAD_SOLID][i]	/=a0;
		msolid.value[RHO_SOLID][i]	/=rho0;
		msolid.value[VPV_SOLID][i]    /=*vp0;
		msolid.value[VPH_SOLID][i]    /=*vp0;
		msolid.value[VSV_SOLID][i]    /=*vp0;
		msolid.value[VSH_SOLID][i]    /=*vp0;
		msolid.value[A_SOLID][i] 	/=*p0;
		msolid.value[C_SOLID][i] 	/=*p0;
		msolid.value[L_SOLID][i] 	/=*p0;
		msolid.value[N_SOLID][i] 	/=*p0;
		msolid.value[F_SOLID][i] 	/=*p0;
		msolid.value[G_SOLID][i]	/=*g0;
		msolid.value[DRHODR_SOLID][i]	/=rho0/a0;
	 }
	 for(j=0;j<PARAMETERS_SOLID;j++)
	 for(i=0;i<msolid.nmodel;i++)
	 msolid.value0[j][i]=msolid.value[j][i];
/* test the effect of gravity 
	 for(i=0;i<msolid.nmodel;i++)
	 msolid.value0[G_SOLID][i]=msolid.value[G_SOLID][i]=0.0;
*/
}
void set_dispersion(a,c,el,en,f,a0,c0,el0,en0,f0,qk,qu,period,n_start,n_end)
/*
	correct elastic constants for physical dispersion for a 
	transversely isotropic body
	reference period is 1 sec
	See note for details
*/
double *a,*c,*el,*en,*f,*a0,*c0,*el0,*en0,*f0,*qk,*qu,period;
int n_start,n_end;
{
	int i;
	double lnt,qa,ql,qf,kappa,mue,kq,uq;
	lnt = -2.0*log(period)/PI;
	for(i=n_start;i<=n_end;i++)
	{
		/* Voigt shear and bulk moduli */
		kappa = (4.0*a0[i]+c0[i]+4.0*f0[i]-4.0*en0[i])/9.0;
		mue = (a0[i]+c0[i]-2.0*f0[i]+5.0*en0[i]+6.0*el0[i])/15.0;
		kq = kappa*qk[i];
		uq = 2.0*mue/3.0;
		qa = (kq+2.0*uq*qu[i])/(kappa+2.0*uq)*lnt+1.0;
		ql = qu[i]*lnt+1.0;
		qf = (kq-uq*qu[i])/(kappa-uq)*lnt+1.0;
		a[i]  =qa*a0[i];
		c[i]  =qa*c0[i];
		el[i] =ql*el0[i];
		en[i] =ql*en0[i];
		f[i]  =qf*f0[i];
	}
}

void set_layersolid(r,region)
double r;
int region;
{
	void huntd();
	huntd(&(msolid.value[RAD_SOLID][msolid.region[region][LOWER_BOUND_SOLID]])
	,msolid.region[region][NUMBER_OF_GRIDS],r,&layer);
}

void set_regionsolid(region,r,type)
int *region,type;
double r;
/*
 * find region from r
 */
{
	int i;

	if( r < msolid.value[RAD_SOLID][msolid.region[0][LOWER_BOUND_SOLID]] ||
	    r > msolid.value[RAD_SOLID][msolid.region[msolid.nregion-1][UPPER_BOUND_SOLID]] )
	{
	FR,"set_regionsolid: r=%f is out of range\n",r);exit(-1);
	}

	switch ( type ){
	case UPPER_REGION: /* if r is just on the discontinuity take upper region */
	if( r ==msolid.value[RAD_SOLID][msolid.region[msolid.nregion-1][UPPER_BOUND_SOLID]] )
		*region = msolid.nregion-1; /* in case r == very top */
	else
	for(i=0;i<msolid.nregion;i++)
	{
		if( msolid.value[RAD_SOLID][msolid.region[i][LOWER_BOUND_SOLID]]<=r
		&&r<msolid.value[RAD_SOLID][msolid.region[i][UPPER_BOUND_SOLID]])
		{ *region=i; break; }
	}
	break;
	case LOWER_REGION: /* if r is just on the discontinuity take lower region */
	if( r==msolid.value[RAD_SOLID][msolid.region[0][LOWER_BOUND_SOLID]] )
		*region = 0; /* in case r == very bottom */
	else
	for(i=0;i<msolid.nregion;i++)
	{
		if(  msolid.value[RAD_SOLID][msolid.region[i][LOWER_BOUND_SOLID]]<r
		&&r<=msolid.value[RAD_SOLID][msolid.region[i][UPPER_BOUND_SOLID]])
		{ *region=i; break;}
	}
	break;
	default:
		FR,"set_regionsolid: unknown type.\n");
		exit(-1);
		break;
	}
}
void non_dimensional_numbers(omg2n, vgn, r1, model)
/*
 * find nondimensional numbers at a give radius and earthmodel
 * take limit r1 -> R. see note p.12-2 for details
 */
	double *omg2n,*vgn; 
	double r1;
	struct model_funcsolid *model;
{
/* model is normalized by normalize_modelsolid() */

	*omg2n = model->g/r1;
	*vgn =model->g*r1*model->rho/model->a;
}
