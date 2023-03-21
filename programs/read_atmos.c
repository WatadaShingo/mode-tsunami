#include <stdio.h>
#include <string.h>
#include <math.h>
#include "earthconst.h"
#include "m_atmos.h"
#include "read_atmos.h"
#define FR		fprintf(stderr
#define PI		3.141592653589793

int layer;

struct model_atmos matmos;
double linear_interp();

#define LINEAR(param,r1,value) \
	value = linear_interp(param[layer],\
				param[layer+1],\
				r[layer+1] -r[layer],\
				r1 - r[layer])
#define LINEAR1(param,value) \
	value = (double)(param[layer+1] - param[layer])*hx + param[layer]
#define NATURAL		1.0e30
int read_atmos(modelfilename)
	char *modelfilename;
{
	FILE *fp;
	char line[100];
	struct model_funcatmos model_unit;
	int i;
	double tmp_array[100],dum,dum1;
	extern void spline();
	extern void splint();

	strcpy(matmos.modelfilename,modelfilename);
	if((fp = fopen(matmos.modelfilename,"r"))==NULL)
	{
		FR,"read_atmos(%s) %d: cannot find file %s.\n",__FILE__,__LINE__,
		matmos.modelfilename);exit(-1);
	}

	fgets(line,100,fp);
	fscanf(fp,"%s %d",matmos.modelname,&(matmos.nmodel));
	fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&(model_unit.r),
				&(model_unit.t),
				&(model_unit.p),
				&(model_unit.rho),
				&(model_unit.vp),
				&(model_unit.kappa),
				&(model_unit.g));
	for(i=0;i<matmos.nmodel;i++)
	{
		fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",
		&(matmos.value[RAD_ATMOS][i]),
		&(matmos.value[T_ATMOS][i]),
		&(matmos.value[P_ATMOS][i]),
		&(matmos.value[RHO_ATMOS][i]),
		&(matmos.value[VP_ATMOS][i]),
		&(matmos.value[KAPPA_ATMOS][i]),
		&(matmos.value[G_ATMOS][i])
		);
/* convert all variables to MKSA unit system */
		matmos.value[RAD_ATMOS][i] += 6371;
		matmos.value[RAD_ATMOS][i] *= model_unit.r;
		matmos.value[T_ATMOS][i] *= model_unit.t;
		matmos.value[P_ATMOS][i] *= model_unit.p;
		matmos.value[RHO_ATMOS][i] *= model_unit.rho;
		matmos.value[VP_ATMOS][i] *= model_unit.vp;
		matmos.value[KAPPA_ATMOS][i] *= model_unit.kappa;
		matmos.value[G_ATMOS][i] *= model_unit.g;
	}
	fclose(fp);

	splined(matmos.value[RAD_ATMOS],matmos.value[RHO_ATMOS],matmos.nmodel,
	NATURAL,NATURAL,tmp_array);
	for(i=0;i<matmos.nmodel;i++)
		splintd(matmos.value[RAD_ATMOS],matmos.value[RHO_ATMOS],
		tmp_array,matmos.nmodel,matmos.value[RAD_ATMOS][i],
		&dum,&(matmos.value[DRHODR_ATMOS][i]),&dum1);

	return 0;
}

void find_modelatmos(model,r1)
	double r1;
	struct model_funcatmos *model;
{
	double *r,*rho,*kappa,*g,*drhodr;
	double hx;

	model->r = r1;
	r = &(matmos.value[RAD_ATMOS][0]);

	/* special treatment for r=min, max. because hunt() routine
	is not reliable(?) when r=min or max */
	if (r1-r[0] < 1e-10*(r[1]-r[0]))
	{
		layer = 0;
		model->rho = matmos.value[RHO_ATMOS][0];
		model->kappa = matmos.value[KAPPA_ATMOS][0];
		model->g = matmos.value[G_ATMOS][0];
		model->drhodr = matmos.value[DRHODR_ATMOS][0];
	}
	else if( r[matmos.nmodel-1]-r1 < 1e-10*(r[matmos.nmodel-1]-r[matmos.nmodel-2]))
	{
		layer = matmos.nmodel-1;
		model->rho = matmos.value[RHO_ATMOS][matmos.nmodel-1];
		model->kappa = matmos.value[KAPPA_ATMOS][matmos.nmodel-1];
		model->g = matmos.value[G_ATMOS][matmos.nmodel-1];
		model->drhodr = matmos.value[DRHODR_ATMOS][matmos.nmodel-1];
	}else{

	set_layeratmos(r1);
	if( layer < 0 || layer > matmos.nmodel-2 )
	{FR,"find_modelatmos (%s) %d: layer is out of range r=%e\n",
	__FILE__,__LINE__,r1);exit(-1);}

	rho = &(matmos.value[RHO_ATMOS][0]);
	kappa = &(matmos.value[KAPPA_ATMOS][0]);
	g = &(matmos.value[G_ATMOS][0]);
	drhodr = &(matmos.value[DRHODR_ATMOS][0]);

	hx = (r1 - r[layer]) /( r[layer+1] -r[layer]);
	if( hx < 0.0 || hx > 1.0 ) printf("hx=%f,r=%f\n",hx,r1);
	/* pointer r is used in LINEAR macro */
/*
	LINEAR(matmos.value[RHO_ATMOS],r1,model->rho);
	LINEAR(matmos.value[KAPPA_ATMOS],r1,model->kappa);
	LINEAR(matmos.value[G_ATMOS],r1,model->g);
	LINEAR(matmos.value[DRHODR_ATMOS],r1,model->drhodr);
*/
	LINEAR1(rho,model->rho);
	LINEAR1(kappa,model->kappa);
	LINEAR1(g,model->g);
	LINEAR1(drhodr,model->drhodr);

	}
}
void normalize_modelatmos(a0,rho0,p0,vp0,g0,omg20)
/* normalize 
struct model_atmos matmos;
input a0	normalization of radius, in meter
     rho0	normalization of density, in kg/m^3
output 
	p0	normalizatin of pressure & elastic constants, in Pa == 1e8 kbar
	vp0	normalization of seismic velocity, in m/sec
	g0	normalization of gravity, in m/sec^2
	omg20	normalization of omg^2 in sec^2

Here we assume three quantities,
universal graviry constant, rho and density, are given and compute
normalization factor of other parameters.
See note for details
*/
double a0,rho0,*p0,*vp0,*g0,*omg20;
{
	int i;

	*g0 = 4.0*PI*UG*rho0*a0;
	*p0 = *g0*rho0*a0;
	*vp0 =sqrt(*g0*a0);
	*omg20 = 4.0*PI*UG*rho0;
	for(i=0;i< matmos.nmodel;i++)
	{
		matmos.value[RAD_ATMOS][i]   /=a0;
		matmos.value[RHO_ATMOS][i]   /=rho0;
		matmos.value[P_ATMOS][i]     /=*p0;
		matmos.value[VP_ATMOS][i]    /=*vp0;
		matmos.value[KAPPA_ATMOS][i] /=*p0;
		matmos.value[G_ATMOS][i]     /=*g0;
		matmos.value[DRHODR_ATMOS][i] /=(rho0/a0);
	}
}

void set_layeratmos(r)
/*
----	r=a   node 0----
	layer 0
----	r=r1	node n----
	layer n
-----	r=r2	node n+1--
	layer n+1

a,r2,r3 are radius from the center of the earth
laynumber starts with 0, the inner most layer

return the minimum node number which satisfies
	r1 <=r <=r2
i.e. if r=r2, return n, not n+1.
*/
double r;
{
	void huntd();
	huntd(matmos.value[RAD_ATMOS],matmos.nmodel, r, &layer);
}
void non_dimensional_numbers_atmos(omg2n,vgn,r1,model)
/*
 * find nondimensional numbers at a give radius and earthmodel
 * take limit r1 -> R. see note p.12-2 for details
 */
	double *omg2n,*vgn;
	double r1;
	struct model_funcatmos *model;
{
/* model is normalized by normalize_modelatmos() */

	*omg2n = model->g/r1;
	*vgn = model->g*r1*model->rho/model->kappa;
}
