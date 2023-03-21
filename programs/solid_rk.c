/*
 test programs of runge-kutta integration
 integration through solid earth (including ocean layer)
 */
#include <stdio.h>
#include <math.h>
#include <sys/file.h>
#include <stdlib.h>
#include "minmax.h"
#include "rk.h"
#include "m_solid.h"
#include "read_solid.h"
#include "model.h"
#include "calc_a.h"
#include "mallocutil.h"
#include "solid_sub.h"
#include "eigen.h"
#define TINY	1e-30
#define YES		1
#define NO		0
#define GRID_NUMBER	5000

#define FLUID_GRID	1
#define SOLID_GRID	2

#define GROUP_V		0x0001
#define PHASE_V		0x0002
#define EIGEN_FUNC	0x0004
#define EIGEN_STORE	0x0008

#define FREE_SURFACE		1
#define RIGID_SURFACE		2
#define NON_REFLECTION		3


#define	FR	fprintf(stderr
#define	IF(x)	((flag & x ) == x)	
#define PERMS		0644	/* o+rwx */	

#define PI	3.141592653589793

int region; /* referred in  calc_a.c */

main()
{
	double i1,i2,i3,gv,pv,eps;
	double old_omg,omg_start,omg_end,domg,omg,tol,best_omg;
	double best_omg2_norm,omg2_norm,old_omg2_norm;
	void solid_driver();
	void zbrent_solid();
	void eigen_process();
	void layer_divider();
	void test_vectors();
	double boundary_condition(),bc_value,old_bc_value,best_bc_value;
	double a0,rho0,p0,vp0,g0,omg20;
	int flag;
	double freq;
	double as,rhos;

	int i,j,k,l,tot_n_solution,n_solution,tot_n_grid,eig_n_grid;
	int layer_start,layer_end,eig_layer_start,eig_layer_end,
	region_start,region_end,start_n_grid,start_eig_n_grid;
	double ys[3][6],yf[2][4],**eig_u,**eig_du,**u,**du;
	double *t_grid,eig_t_grid[GRID_NUMBER],
	***y_grid,***dydt_grid,**dum_y_grid,**dum_dydt_grid;
	extern struct model_solid msolid; /* defined in read_solid.c */
	double *a,*c,*el,*en,*f,*rho,*g;
	double *qk,*qu;
	struct model_funcsolid model;
	int key,fd_eigenset,print_counter,is_disper,is_cowl,bc_type;
	char modelfilename[100],
	ascii_eigenfile[100],binary_eigensetname[100],
	binary_eigeninfofile[100];
	int n_divide_info,idum,divide_grid;
	double ddum,bc_coeff[6];
	double eig_divide_info[MAX_N_REGIONS];/* layer is divided into this thick */
	int eig_region_info[MAX_N_REGIONS][3];/* same to msolid.region */
	char eig_layer_sfflag[GRID_NUMBER];/* if layer is sold 1,fluid 0 */
	char layer_sfflag[GRID_NUMBER];
	int layer_table[GRID_NUMBER];
	/* table j=table[i], where t_grid[i] == eig_t_grid[j] */
	FILE *fp;
	char ckey;

	eps = 1e-7;
	tol = 1e-6;
	a0  = 6.371e6;
	rho0 = 5.5e3;
	flag = 0;
	flag = GROUP_V | PHASE_V;
	as = 1000.0; /* 1000 meter */
	rhos=1000.0; /* 1000 Kg/m^3 */
	print_counter = 0;
	
	printf("model filename >>");
	scanf("%s",modelfilename);
	read_solid(modelfilename);
/*
 * set computation flags
 */
	do{
	printf("include physical dispersion (y/n)");
	scanf("%c",&ckey);printf("\n");
	}
	while(!(ckey =='y'|| ckey=='n'));
	is_disper = (ckey == 'y' ? YES : NO );

	do{
	printf("Cowling approximation (y/n)");
	scanf("%c",&ckey);printf("\n");
	}
	while(!(ckey =='y'|| ckey=='n'));
	is_cowl = (ckey == 'y' ? YES : NO );
	do{
	printf("Top boundary condition stress free surface    1:\n");
	printf("                       rigid surface          2:\n");
	printf("                       non-reflection surface 3:\n");
	scanf("%d",&bc_type);
	}while(bc_type <1 || bc_type >3);
	do {
	printf("top boundary grid points 1-%d -- (0 to set default %d)>>",
		msolid.nmodel,msolid.nmodel);
	scanf("%d",&tot_n_grid);
	}while(tot_n_grid< 0 || tot_n_grid > msolid.nmodel);
	if( tot_n_grid ==0 ) tot_n_grid = msolid.nmodel;
	normalize_modelsolid(a0,rho0,&p0,&vp0,&g0,&omg20);
/*
 *set eig_t_grid[],eig_divide_info[],eig_region_info[][],eig_layer_sfflag[]
 */
	for(i=0;i<MAX_N_REGIONS;i++)eig_divide_info[i]=0.0;
	printf("number of regions to be divided >>");
	scanf("%d",&n_divide_info); printf("\n");
	if( n_divide_info > MAX_N_REGIONS )
	{FR,"input error: too many region divide info=%d\n",n_divide_info);
	exit(-1);}
	for(i=0;i<n_divide_info;i++)
	{
		printf("i-th region and maximum grid interval in meter >>");
		scanf("%d %lf",&idum,&ddum); printf("\n");
		if( idum >= msolid.nregion || idum < 1 || ddum <0.0){
		FR,"input error for layer divider: layer=%d, thickness=%lf\n",
		idum,ddum);exit(-1);}
		eig_divide_info[idum] = ddum/a0;
	}
	eig_n_grid=0;
	for(i=0;i<msolid.nregion;i++)
	{
		layer_start = msolid.region[i][LOWER_BOUND_SOLID];
		eig_region_info[i][LOWER_BOUND_SOLID]=eig_n_grid;
		layer_divider(&(eig_t_grid[eig_n_grid]),
			&divide_grid,&(layer_table[layer_start]),
			eig_divide_info[i],
			&(msolid.value[RAD_SOLID][layer_start]),
			msolid.region[i][NUMBER_OF_GRIDS]);
		for(j=0;j<msolid.region[i][NUMBER_OF_GRIDS];j++)
			layer_table[layer_start+j] += eig_n_grid;
		for(j=0;j<divide_grid;j++)
			eig_layer_sfflag[j+eig_n_grid]=
			( msolid.value[N_SOLID][msolid.region[i][LOWER_BOUND_SOLID]]==0.0 ? 
			FLUID_GRID : SOLID_GRID);
		eig_n_grid += divide_grid ;
		eig_region_info[i][UPPER_BOUND_SOLID]=eig_n_grid-1;
		eig_region_info[i][NUMBER_OF_GRIDS] = divide_grid;
		printf("region[%d] is divided into from %d to %d layer(s)\n",
			i,msolid.region[i][NUMBER_OF_GRIDS],divide_grid);
	}
	for(i=0;i<msolid.nmodel;i++)
		layer_sfflag[i]=(msolid.value[N_SOLID][i]==0.0?FLUID_GRID:SOLID_GRID);
/*
 * test
	for(i=0;i<msolid.nregion;i++)
		printf("eig_region_info[%2d] %d %d %d\n",i,
		eig_region_info[i][LOWER_BOUND_SOLID],
		eig_region_info[i][UPPER_BOUND_SOLID],
		eig_region_info[i][NUMBER_OF_GRIDS]);
	for(i=0;i<msolid.nmodel;i++)
		printf("eig_t_grid[%3d]=%7.4f t_grid[%3d]=%7.4f\n",
		i,eig_t_grid[i],i,msolid.value[RAD_SOLID][i]);
	for(i=msolid.nmodel;i<eig_n_grid;i++)
		printf("eig_t_grid[%2d]=%7.4f\n", i,eig_t_grid[i]);
	for(i=0;i<msolid.nmodel;i++)
		printf("layer_table[%3d]=%3d\n",i,layer_table[i]);
 */
/*
 * initialization
 */
	msolid.nmodel = tot_n_grid;;
	msolid.region[msolid.nregion-1][UPPER_BOUND_SOLID] = tot_n_grid-1;
	eig_region_info[msolid.nregion-1][UPPER_BOUND_SOLID] = 
	layer_table[msolid.nmodel-1];
	eig_n_grid = layer_table[msolid.nmodel-1]+1;
	tot_n_grid= msolid.nmodel;
	tot_n_solution = n_solution = 0;
	layer_start=1;

	y_grid = malloc_m3double(0,2,0,5,0,eig_n_grid-1,"main y_grid");
	dydt_grid = malloc_m3double(0,2,0,5,0,eig_n_grid-1,"main dydt_grid");
	dum_y_grid = malloc_mdouble(0,6-1,0,eig_n_grid-1,"main dum_y_grid");
	dum_dydt_grid = malloc_mdouble(0,6-1,0,eig_n_grid-1,"main dum_dydt_grid");
	a = malloc_vdouble(0,eig_n_grid-1,"main a");
	c = malloc_vdouble(0,eig_n_grid-1,"main c");
	el = malloc_vdouble(0,eig_n_grid-1,"main el");
	en = malloc_vdouble(0,eig_n_grid-1,"main en");
	f = malloc_vdouble(0,eig_n_grid-1,"main f");
	rho = malloc_vdouble(0,eig_n_grid-1,"main rho");
	g = malloc_vdouble(0,eig_n_grid-1,"main g");
	qk = malloc_vdouble(0,eig_n_grid-1,"main qk");
	qu = malloc_vdouble(0,eig_n_grid-1,"main qu");
	u = malloc_mdouble(0,2,0,tot_n_grid-1,"main u");
	du = malloc_mdouble(0,2,0,tot_n_grid-1,"main du");
	eig_u = malloc_mdouble(0,2,0,eig_n_grid-1,"main eig_u");
	eig_du = malloc_mdouble(0,2,0,eig_n_grid-1,"main eig_du");

/*	for(i=0;i<tot_n_grid;i++) t_grid[i] = (double)msolid.value[RAD_SOLID][i];*/
	t_grid = msolid.value[RAD_SOLID];
	/*
	 * if no dispersion, set the earthmodel at eig_t_grid[]
	 * if dispersion, set after period is fixed and earth model is set
	 * by set_dispersion().
	 */
	if( is_disper != YES )
	for(i=0;i<msolid.nregion;i++)
	{
/*
		layer_start = msolid.region[i][LOWER_BOUND_SOLID];
		layer_end = msolid.region[i][UPPER_BOUND_SOLID];
*/
		layer_start = eig_region_info[i][LOWER_BOUND_SOLID];
		layer_end = eig_region_info[i][UPPER_BOUND_SOLID];
		for(j=layer_start;j<=layer_end;j++)
		{
/*			find_modelsolid(&model,t_grid[j],i);*/
			find_modelsolid(&model,eig_t_grid[j],i);
			rho[j] = model.rho;
			a[j]=model.a;
			c[j]=model.c;
			el[j]=model.l;
			en[j]=model.n;
			f[j]=model.f;
			g[j]=model.g;
			/* a, c, el, en, f are used in eigenprocess() */ 
		}
	}

	for(j=0;j<5;j++)
	for(k=0;k<eig_n_grid;k++)
	{
		dum_y_grid[j][k]= dum_dydt_grid[j][k]= 0.0;
		for(i=0;i<3;i++)
		y_grid[i][j][k]= dydt_grid[i][j][k]= 0.0;
	}
/*
 * print info
 */
	printf("<< compute eigenfunction >>>\n");
	printf("model name      : %s [in file %s]\n",
		msolid.modelname,msolid.modelfilename);
	printf("number of grid  : %d\n",msolid.nmodel);
	printf("number of region: %d\n",msolid.nregion);

	puts("print out ascii eigenfunctions 1:");
	puts("store binary eigenfunctions    2:");
	puts("both                           3:");
	puts("do nothing                     0:");
	do{
		scanf("%d",&key);
		printf("\n");
		if( key== 1 || key == 3 )
		{ 	flag |= EIGEN_FUNC;
			puts("file name prefix for ascii eigenfunction >>"); 
			scanf("%s",ascii_eigenfile);
			printf("eigenfunction is stored in %s.?s'angular order'\n",
			ascii_eigenfile);
			printf("\n");
		}
		if( key== 2 || key == 3 ) 
		{	flag |= EIGEN_STORE;
			puts("file name of a set of eigenfunctions >>"); 
			scanf("%s",binary_eigensetname);
			printf("eigenfunctions are stored in %s\n",
			binary_eigensetname);
			printf("\n");
			if( (fd_eigenset=creat(binary_eigensetname,PERMS)) == -1)
			{
				FR,"main: cannot create file %s\n",
				binary_eigensetname); exit(-1);
			}
			sprintf(binary_eigeninfofile,"%s.info",binary_eigensetname);
			if((fp=fopen(binary_eigeninfofile,"w")) == NULL)
			{
				FR,"main: cannot creat file %s\n",
				binary_eigeninfofile);exit(-1);
			}
			fprintf(fp,"model : %s\n",msolid.modelfilename);
			if( is_disper == YES )
				fprintf(fp,"physical dispersion of solid Earth is included\n");
			else
				fprintf(fp,"physical dispersion of solid Earth is NOT included\n");
			if( is_cowl == YES )
				fprintf(fp,"Cowling approximation is used\n");
			else
				fprintf(fp,"Cowling approximation is NOT used\n");
			fprintf(fp,"grid  : 1-%d ( 0.0-%.fkm)\n",
				tot_n_grid,msolid.value[RAD_SOLID][tot_n_grid-1]*a0/1e3);
			fprintf(fp,"===================================================\n");
			fprintf(fp,"angular order and angular frequency \n");
		}
	}while( key < 0 || key > 3);
/*
 * Big loop start
 */
/*
 * determine top is fluid or not
 */
	printf("input l, starting omg, ending omg, domg, start level to stop l=-9999>>");
	printf("\n");
	while( (key=scanf("%d %lf %lf %lf %d",&l,&omg_start,&omg_end,&domg,&layer_start))!=EOF)
	{
	if( l < -1000 )break;
	if( IF(EIGEN_STORE)) fprintf(fp,"%5d ",l);
	n_solution = 0;
	print_counter = 0;
	if ( key == 4 ) continue;
	else if( key != 5 ){FR,"input error!\n");exit(0);}
	if( IF(EIGEN_STORE))
	printf("angular_order=%d omg_start=%lf omg_end=%lf domg=%lf start grid=%d\n",
	l,omg_start,omg_end,domg,layer_start);

/* initial condition */
	omg = omg_start;
	old_omg = omg_start; 

	omg2_norm= omg*omg/omg20;
	old_omg2_norm = old_omg*old_omg/omg20;

	set_regionsolid(&region_start, msolid.value[RAD_SOLID][layer_start],LOWER_REGION);
	region_end = msolid.nregion-1;
	layer_end = msolid.region[region_end][UPPER_BOUND_SOLID];
	/* check layer_start region_start */
	if( layer_start > msolid.region[region_start][UPPER_BOUND_SOLID])
		region_start ++;
	if( layer_start < msolid.region[region_start][LOWER_BOUND_SOLID])
		layer_start++;

	if( layer_start > msolid.region[region_start][UPPER_BOUND_SOLID] || 
	    layer_start < msolid.region[region_start][LOWER_BOUND_SOLID])
	{FR,"layer_start=%d is not in region=%d.\n",layer_start,region_start);
	exit(-1);}

	eig_layer_start = layer_table[layer_start];
	eig_layer_end   = layer_table[layer_end];
/*
 * compute eigen frequency
 */
	start_n_grid = msolid.region[region_start][UPPER_BOUND_SOLID]-layer_start+1;
	if( start_n_grid == 1 ) { layer_start -- ; start_n_grid++; }
	start_eig_n_grid = eig_region_info[region_start][UPPER_BOUND_SOLID]-eig_layer_start+1;
	if( start_n_grid == 1 ) { eig_layer_start -- ; start_eig_n_grid++; }

/* first try */
	if( is_disper == YES )
	set_dispersion(msolid.value[A_SOLID],msolid.value[C_SOLID],
	msolid.value[L_SOLID],msolid.value[N_SOLID],msolid.value[F_SOLID],
	msolid.value0[A_SOLID],msolid.value0[C_SOLID],
	msolid.value0[L_SOLID],msolid.value0[N_SOLID],msolid.value0[F_SOLID],
	msolid.value0[QK_SOLID],msolid.value0[QU_SOLID],
	2.0*PI/sqrt(omg2_norm*omg20),0,221);/* <<====THIS IS AD HOC */
	if( layer_sfflag[layer_start]== FLUID_GRID )
		y_initial_fluid4by4(yf[0],yf[1],l,t_grid[layer_start],omg2_norm);
	else 
		y_initial_solid6by6(ys[0],ys[1],ys[2],l,t_grid[layer_start],omg2_norm);
	solid_driver(y_grid,dydt_grid,ys,yf,&layer_end,layer_sfflag,t_grid,start_n_grid,eps,l,omg2_norm,omg20,
		NO,NO,is_disper,is_cowl,layer_start,msolid.region,region_start,region_end,dum_y_grid,dum_dydt_grid);

	old_bc_value = boundary_condition(ys,yf,l,omg2_norm,
		layer_sfflag[layer_end],is_cowl,bc_type,bc_coeff);
	printf("omg=%11.8f bc_value=%+11.4e\n", omg,old_bc_value);

	if( old_bc_value == 0.0)
	{
		best_omg = sqrt(omg2_norm*omg20);
		freq = best_omg/2/PI;
		printf("found[%d] l=%4d omg = %15.8e period=%.2f bc= %15.8e\n",
		n_solution,l,best_omg,1.0/freq,old_bc_value);
		if( is_disper == YES )
		for(i=0;i<msolid.nregion;i++)
			set_modelsolid_array(a,c,el,en,f,rho,g,qk,qu,eig_t_grid,i,
			eig_region_info[i][LOWER_BOUND_SOLID],eig_region_info[i][UPPER_BOUND_SOLID]);

		eigen_process(&i1,&i2,&i3,u,du,y_grid,dydt_grid,t_grid,start_n_grid,
		layer_start,msolid.region,omg2_norm,l,flag,dum_y_grid,dum_dydt_grid,
		a,c,el,en,f,rho,g,is_cowl);
		pv=sqrt(omg2_norm)*6371e3/a0/(l+0.5);
		gv=i3/i1*pv*(l+0.5)*(l+0.5)/omg2_norm;
		/* group & phase velocity depends on the radius we refer */
		printf("I2/I1 : omg2 = %12.5e Vc=%12.5e Vg=%12.5e\n",
			i2/(i1*omg2_norm),pv*vp0,gv*vp0);

		n_solution ++;
	}

	while( omg  < omg_end)
	{
		omg = omg + domg;
		omg2_norm= omg*omg/omg20;
		old_omg2_norm = old_omg*old_omg/omg20;

	if( is_disper == YES ){
	set_dispersion(msolid.value[A_SOLID],msolid.value[C_SOLID],
	msolid.value[L_SOLID],msolid.value[N_SOLID],msolid.value[F_SOLID],
	msolid.value0[A_SOLID],msolid.value0[C_SOLID],
	msolid.value0[L_SOLID],msolid.value0[N_SOLID],msolid.value0[F_SOLID],
	msolid.value0[QK_SOLID],msolid.value0[QU_SOLID],
	2.0*PI/sqrt(omg2_norm*omg20),0,221);/* <<====THIS IS AD HOC */
	if( layer_sfflag[layer_start]== FLUID_GRID )
		y_initial_fluid4by4(yf[0],yf[1],l,t_grid[layer_start],omg2_norm);
	else 
		y_initial_solid6by6(ys[0],ys[1],ys[2],l,t_grid[layer_start],omg2_norm);
	}
		solid_driver(y_grid,dydt_grid,ys,yf,&layer_end,layer_sfflag,t_grid,start_n_grid,eps,l,omg2_norm,omg20,
		NO,NO,is_disper,is_cowl,layer_start,msolid.region,region_start,region_end,dum_y_grid,dum_dydt_grid);
		bc_value = boundary_condition(ys,yf,l,omg2_norm,
			layer_sfflag[layer_end],is_cowl,bc_type,bc_coeff);
		if( bc_value == 0.0)
		{
			best_omg = sqrt(omg2_norm*omg20);
			freq = best_omg/2/PI;
			printf("found[%d] l=%4d omg = %15.8e period=%.2f bc= %15.8e\n",
			n_solution,l,best_omg,1.0/freq,bc_value);
			if( is_disper == YES )
			for(i=0;i<msolid.nregion;i++)
				set_modelsolid_array(a,c,el,en,f,rho,g,qk,qu,eig_t_grid,i,
				eig_region_info[i][LOWER_BOUND_SOLID],eig_region_info[i][UPPER_BOUND_SOLID]);

			eigen_process(&i1,&i2,&i3,u,du,y_grid,dydt_grid,eig_t_grid,eig_n_grid,
			eig_layer_start,eig_region_info,omg2_norm,l,flag,dum_y_grid,dum_dydt_grid,
			a,c,el,en,f,rho,g,is_cowl);
			pv=sqrt(omg2_norm)*6371e3/a0/(l+0.5);
			gv=i3/i1*pv*(l+0.5)*(l+0.5)/omg2_norm;
		/* group & phase velocity depends on the radius we refer */
			printf("I2/I1 : omg2 = %12.5e Vc=%12.5e Vg=%12.5e\n",
				i2/i1/omg2_norm,pv*vp0,gv*vp0);

			n_solution ++;
		}
		else if( bc_value * old_bc_value < 0.0 )
		/* solution is bounded between old_omg and omg */
		{
			zbrent_solid(&best_omg2_norm,old_omg2_norm,omg2_norm,omg20,
			is_disper,is_cowl,old_bc_value,
			bc_value,y_grid,dydt_grid,dum_y_grid,dum_dydt_grid,
			tol,ys,yf,layer_sfflag,t_grid,start_n_grid,eps,l,
			layer_start,msolid.region,region_start,region_end,
			bc_type,bc_coeff);

/* gready eigen frequency search using fine grid can be avoided here 
			zbrent_solid(&best_omg2_norm,old_omg2_norm,omg2_norm,omg20,
			is_disper,is_cowl,old_bc_value,
			bc_value,y_grid,dydt_grid,dum_y_grid,dum_dydt_grid,
			tol,ys,yf,eig_layer_sfflag,eig_t_grid,start_eig_n_grid,eps,l,
			eig_layer_start,eig_region_info,region_start,region_end,
			bc_type,bc_coeff);
*/
		/*
		 * eigenfrequency is found.
		 * compute eigen functions at eig_t_grid[].
		 */
			for(j=0;j<5;j++)
/*			for(k=layer_start;k<eig_n_grid;k++)*/
			for(k=eig_layer_start;k<eig_n_grid;k++)
			{
				dum_y_grid[j][k]= dum_dydt_grid[j][k]= 0.0;
				for(i=0;i<3;i++)
				y_grid[i][j][k]= dydt_grid[i][j][k]= 0.0;
			}
	if( is_disper == YES ){
	set_dispersion(msolid.value[A_SOLID],msolid.value[C_SOLID],
	msolid.value[L_SOLID],msolid.value[N_SOLID],msolid.value[F_SOLID],
	msolid.value0[A_SOLID],msolid.value0[C_SOLID],
	msolid.value0[L_SOLID],msolid.value0[N_SOLID],msolid.value0[F_SOLID],
	msolid.value0[QK_SOLID],msolid.value0[QU_SOLID],
	2.0*PI/sqrt(best_omg2_norm*omg20),0,221);/* <<====THIS IS AD HOC */
	if( layer_sfflag[layer_start]== FLUID_GRID )
		y_initial_fluid4by4(yf[0],yf[1],l,t_grid[layer_start],best_omg2_norm);
	else 
		y_initial_solid6by6(ys[0],ys[1],ys[2],l,t_grid[layer_start],best_omg2_norm);
	}
/*
			solid_driver(y_grid,dydt_grid,ys,yf,&layer_end,layer_sfflag,t_grid,start_n_grid,eps,l,best_omg2_norm,omg20,
			NO,YES,is_disper,is_cowl,layer_start,msolid.region,region_start,region_end,dum_y_grid,dum_dydt_grid);
*/
			solid_driver(y_grid,dydt_grid,ys,yf,&eig_layer_end,eig_layer_sfflag,eig_t_grid,start_eig_n_grid,eps,l,best_omg2_norm,omg20,
			NO,YES,is_disper,is_cowl,eig_layer_start,eig_region_info,region_start,region_end,dum_y_grid,dum_dydt_grid);
/*			best_bc_value =
			boundary_condition(ys,yf,l,best_omg2_norm,layer_sfflag[layer_end],is_cowl,bc_type,bc_coeff); 
*/
			best_bc_value =
			boundary_condition(ys,yf,l,best_omg2_norm,eig_layer_sfflag[eig_layer_end],is_cowl,bc_type,bc_coeff); 
/*			test_vectors(ys,yf,eig_layer_sfflags[eig_layer_end],is_cowl);*/
			best_omg = sqrt(best_omg2_norm*omg20);
			freq = best_omg/2/PI;
			printf("found[%d] l=%4d omg = %15.8e period=%.2f bc= %15.8e\n",
			n_solution,l,best_omg,1.0/freq,best_bc_value);
			if( is_disper == YES )
			for(i=0;i<msolid.nregion;i++)
				set_modelsolid_array(a,c,el,en,f,rho,g,qk,qu,eig_t_grid,i,
				eig_region_info[i][LOWER_BOUND_SOLID],eig_region_info[i][UPPER_BOUND_SOLID]);
/*
			eigen_process(&i1,&i2,&i3,u,du,y_grid,dydt_grid,t_grid,tot_n_grid,
			layer_start,msolid.region,best_omg2_norm,l,flag,dum_y_grid,dum_dydt_grid,
			a,c,el,en,f,rho,g,is_cowl);
*/
			eigen_process(&i1,&i2,&i3,eig_u,eig_du,y_grid,dydt_grid,eig_t_grid,eig_n_grid,
			eig_layer_start,eig_region_info,best_omg2_norm,l,flag,dum_y_grid,dum_dydt_grid,
			a,c,el,en,f,rho,g,is_cowl);

			pv=sqrt(best_omg2_norm)*6371.0e3/a0/(l+0.5);
			gv=i3/i1*pv*(l+0.5)*(l+0.5)/best_omg2_norm;
	if(IF(EIGEN_FUNC)  || IF(EIGEN_STORE))
	{
		for(i=0;i<tot_n_grid;i++)
		for(j=0;j<3;j++)
		{
			u[j][i] = eig_u[j][layer_table[i]];
			du[j][i] = eig_du[j][layer_table[i]];
		}
			normalize_eigen(u,du,t_grid,tot_n_grid,a0,as,
			rho0,rhos,g0,i1,is_cowl);
			if( IF(EIGEN_FUNC) )
			print_eigen(ascii_eigenfile,
			t_grid,tot_n_grid,u,du,&(msolid.value[RHO_SOLID][0]),
			a0,l,n_solution,freq,is_cowl);
			if( IF(EIGEN_STORE) )
			{
			if( is_cowl == YES)
			write_eigencowl(fd_eigenset,best_omg,l,n_solution,u,du,tot_n_grid);
			else
			write_eigen(fd_eigenset,best_omg,l,n_solution,u,du,tot_n_grid);
			if( print_counter++==5)
				{fprintf(fp,"\n      ");print_counter=0;}
			fprintf(fp,"%10.4e/%10.4e ",best_omg,gv*vp0);
			}
	}
		/* group & phase velocity depends on the radius we refer */
			printf("I2/I1 : omg2 = %15.8e Vc=%15.8e Vg=%15.8e\n",
				i2/i1/best_omg2_norm,pv*vp0,gv*vp0);

			n_solution ++;
		}
		old_omg = omg;
		old_bc_value = bc_value;

#ifdef DBG
/*
 * print out 
 */
		printf("omg=%11.8f bc_value=%+11.4e\n",
		omg,bc_value);
#endif DBG
		
	}
/*
 * Big loop end
 */
	printf("input l, satarting omg, ending omg, domg, start level to stop l=-9999>>");
	printf("\n");
	if( IF(EIGEN_STORE)) fprintf(fp,"\n");
	tot_n_solution += n_solution;
	}
	printf("in total %d eigenfunctions are found\n",tot_n_solution);
	free_m3double(y_grid,0,2,0,5,0);
	free_m3double(dydt_grid,0,2,0,5,0);
	free_mdouble(dum_y_grid,0,5,0);
	free_mdouble(dum_dydt_grid,0,5,0);
	free_vdouble(a,0);
	free_vdouble(c,0);
	free_vdouble(el,0);
	free_vdouble(en,0);
	free_vdouble(f,0);
	free_vdouble(rho,0);
	free_vdouble(g,0);
	free_vdouble(qk,0);
	free_vdouble(qu,0);
	free_mdouble(eig_u,0,2,0);
	free_mdouble(eig_du,0,2,0);
	free_mdouble(u,0,2,0);
	free_mdouble(du,0,2,0);
	if( IF(EIGEN_STORE)){ close(fd_eigenset); fclose(fp); }

	return (0);
}

#define	NATURAL	1.0e30
void eigen_process(i1,i2,i3,u,du,x,dx,r_grid,n_grid,layer_start,region_info,omg2,l,flag,x0,dx0,
a,c,el,en,f,rho,g,is_cowl) 
double *i1,*i2,*i3,**u,**du,***x,***dx,*r_grid,omg2,**x0,**dx0; 
int n_grid,layer_start,l,flag; 
double *a,*c,*el,*en,*f,*rho,*g;
int region_info[MAX_N_REGIONS][3];/* same to msolid.region */
int is_cowl;
/*
	compute egenfcuntions from x[i] i=1-6, normalize and group and phase
	velocities
<< output >>
i1,	integral I1
i2,	integral I2
i3,	integral I3
u[3][] eigenfunctions U, V and PHI
du[3][] its derivative
<< input >>
x[3][6][]	three independent solutions
dx[3][6][]	their derivatives
n_grid	number of grid points of x[][][j],dx[][][j],u[][j] and du[][j]
layer_start	eigenfunction is computed at and above this grid
region_info[][3] region information
r_grid[]	radius of grid points
omg2		eigenfrequency of this mode
l		angular order
flag		GROUP_V, PHASE_V, EIGEN_FUNC
x0[6][]	working area
dx0[6][]	working area
a c el en f rho g  Earth model parameters
is_cowl	if YES Cowling approximation is used
*/
{ 
	int i,j,upper_layer,lower_layer,n_layer;
	double *tmp_array,dum0,dum2,*dnull;
	extern struct model_solid msolid; /* defined in read_solid.c */
	extern void splined();
	extern void splintd();
	dnull = (double *)NULL;

	tmp_array = malloc_vdouble(0,n_grid-1,"eigen_process tmp_array");
	for(i=0;i<3;i++)
	for(j=0;j<n_grid;j++)
	u[i][j]=du[i][j] = 0.0;

	if(is_cowl == YES) construct_eigen_xcowl(x0,x,el,layer_start,n_grid-1);
	else               construct_eigen_x(x0,x,el,layer_start,n_grid-1);
/*	construct_eigen_x(dx0,dx,el,layer_start,n_grid-1);
	## NOTE ##
	this routine tries to construct a dx/dr from the solution of
	three independent solutions of dx/dr. Since dx/dr does not satisfy
	the same boundary consditions x does, this routine does not
	proudce correct dx/dr. dx/dr should not be used.
	## NOTE END ## Aug. 26, 1994 
*/
	if( is_cowl == YES )
	construct_eigen_x3cowl(x0[2],x0[3],x0[0],x0[1],r_grid,l,omg2,el,rho,g,layer_start,n_grid-1);
	else
	construct_eigen_x3(x0[2],x0[3],x0[0],x0[1],x0[4],r_grid,l,omg2,el,rho,g,layer_start,n_grid-1);

	if( (flag & PHASE_V) == PHASE_V  ||
	    (flag & GROUP_V) == GROUP_V )
	{
		*i1 = integral_I1(x0[0],x0[2],rho,r_grid,layer_start,n_grid-1);
		if( is_cowl == YES){
		*i2 = integral_I2(l,x0[0],x0[1],x0[2],x0[3],dnull,dnull,rho,g,a,c,el,en,f,
			r_grid,layer_start,n_grid-1);
		*i3 = integral_I3(l,x0[0],x0[1],x0[2],x0[3],dnull,dnull,rho,g,a,c,el,en,f,
			r_grid,layer_start,n_grid-1);
		}else{
		*i2 = integral_I2(l,x0[0],x0[1],x0[2],x0[3],x0[4],x0[5],rho,g,a,c,el,en,f,
			r_grid,layer_start,n_grid-1);
		*i3 = integral_I3(l,x0[0],x0[1],x0[2],x0[3],x0[4],x0[5],rho,g,a,c,el,en,f,
			r_grid,layer_start,n_grid-1);
		}
	}

	if(IF(EIGEN_FUNC)  || IF(EIGEN_STORE))
	{
		if( is_cowl == YES )
		construct_eigen_uvcowl(u[0],du[0],u[1],du[1],x0,dx0,
		r_grid,l,omg2,c,el,f,rho,g,layer_start,n_grid-1);
		else
		construct_eigen_uvf(u[0],du[0],u[1],du[1],u[2],du[2],x0,dx0,
		r_grid,l,omg2,c,el,f,rho,g,layer_start,n_grid-1);
	/*
	 * for fluid region dV/dr =du[1] is set by 
	 * the interpolation of V=u[1] 
	 */
		for(j=0;j<msolid.nregion;j++)
		if( el[region_info[j][UPPER_BOUND_SOLID]] == 0.0 )
		{
			upper_layer = region_info[j][UPPER_BOUND_SOLID];
			lower_layer = region_info[j][LOWER_BOUND_SOLID];
			if( upper_layer <= layer_start ) continue;
			if( lower_layer < layer_start )
				lower_layer = layer_start;

			n_layer = upper_layer - lower_layer +1;

			splined(&(r_grid[lower_layer]),&(u[1][lower_layer]),n_layer,
			NATURAL,NATURAL,tmp_array);
			for(i=0;i<n_layer;i++)
			splintd(&(r_grid[lower_layer]),&(u[1][lower_layer]),
				tmp_array,n_layer, r_grid[lower_layer+i],
				&dum0, &(du[1][lower_layer+i]), &dum2);
		}
	}
	free_vdouble(tmp_array,0);
}
#undef NATURAL

double boundary_condition(ys,yf,l,omg2,solid_fluid,is_cowl,bc_type,coeff)
	double ys[][6],yf[][4],omg2,coeff[6];
	int l,solid_fluid,is_cowl,bc_type;
/*
bc_type is either FREE_SURFACE, RIGID_SURFACE or NON_REFLECTION
if bc_type==NON_REFLECTION then sum_over_i (coeff[i]*y[i])= 0 is the B.C.
*/
{
	double value;
	if( is_cowl == YES )
	{
		switch(solid_fluid ){
		case FLUID_GRID:
			switch(bc_type){
			case FREE_SURFACE: value = yf[0][1];break;/* normal stress */
			case RIGID_SURFACE: value = yf[0][0];break;
			case NON_REFLECTION: value=coeff[0]*yf[0][0]+coeff[0]*yf[0][1];
			break;
			}
			break;
		case SOLID_GRID:
			switch(bc_type){
			case FREE_SURFACE:
			value = ys[0][1]*ys[1][3] - ys[1][1]*ys[0][3];break;
			case RIGID_SURFACE:
			value = ys[0][0]*ys[1][2] - ys[1][0]*ys[0][2];break;
			case NON_REFLECTION: goto error; break;
			}
			break;
		default:
			FR,"boundary_condition(%s) %d: unknown type=%d.\n",
			__FILE__,__LINE__,solid_fluid);
			exit(-1);
			break;
		}
	}else{
	switch(solid_fluid){
		case FLUID_GRID:
			switch(bc_type){
			case FREE_SURFACE: value = determinant2by2(1,3,yf);break;
			case RIGID_SURFACE: value = determinant2by2(0,2,yf);break;
			case NON_REFLECTION: goto error; break;
			}
			break;
		case SOLID_GRID:
			switch(bc_type){
			case FREE_SURFACE: value = determinant3by3(1,3,5,ys); break;
			case RIGID_SURFACE: value = determinant3by3(0,2,5,ys); break;
			case NON_REFLECTION: goto error; break;
			}
			break;
		default:
			FR,"boundary_condition(%s) %d: unknown type=%d.\n",
			__FILE__,__LINE__,solid_fluid);
			exit(-1);
			break;
		}
	}
	return value;
	error:
	FR,"boundary_condition(%s) %d: this boundary case is not supported.\n",
	__FILE__,__LINE__);
	FR,"is_cowl     =%d\n",is_cowl);
	FR,"solid_folid =%d\n",solid_fluid);
	FR,"bc_type     =%d\n",bc_type);
	exit(-1);
}

void solid_driver(y_grid,dydt_grid,ys,yf,layer1,sfflag,t_grid,n_grid,eps,l,omg2,omg20,
is_scale,is_func,is_disper,is_cowl,layer0,region_info,region_start,region_end,dum_y_grid,dum_dydt_grid)
	double ***y_grid, ***dydt_grid,ys[3][6],yf[2][4],*t_grid,eps,omg2,omg20,
	**dum_y_grid,**dum_dydt_grid;
	int l,n_grid,layer0,*layer1;
	int region_info[MAX_N_REGIONS][3],region_start,region_end;
	int is_scale,is_func,is_disper,is_cowl;
	char *sfflag;
/*
	integrate eq. of motion from 'layer' in 'region_start' with initial value
	'ys' in solid or 'yf' in fluid, depeinding on 'sf' flag through the
	earth model to the 'region_end'( region_end is integrated ).

<< OUTPUT >>
	y_grid[3][6][n_grid],	 eigenfunction solution
	dydt_grid[3][6][n_grid], eigenfunction solution derivative
	ys[3][6],	value of y_grid[][][] at top of region_end when sf == SOLID_GRID
	yf[2][4],	value of y_grid[][][] at top of region_end when sf == FLUID_GRID
	layer1,	integration stops at this layer
<< INPUT >>
	sfflag[],	SOLID_GRID or FLUID_GRID flags of t_grid[]
	t_grid[],	grid of earth model radius t_grid[]
	n_grid,	# of grid points of t_grid[] within this region from layer
	eps,		required accuracy
	l,		angular order
	omg2,		angular frequency (normalized)
	omg20,	normalization factor of (angular freq)^2
	is_scale,	do scaling == YES, if not== NO
	is_funcs,   store all eigenfunctions == YES if not == NO
	is_disper,	if physical dispersion included ==YES, if not == NO
	is_cowl,	if Cowling approximation is used ==YES, if not == NO
	layer0,	integration starts at this layer
	region_info[][] region infomation
	region_start, layer is in this region;
	region_end,	integration stops after integrated through reagion_end
	dum_y_grid[6][n_grid], working area
	dum_dydt_grid[6][n_grid], working area
*/
{
	int i,j,k,n,n_cross,sf;
	extern struct model_solid msolid;	/* defined in read_solid.c */

	region = region_start;
	if( (sf=sfflag[layer0]) == FLUID_GRID )
	y_initial_fluid4by4(yf[0],yf[1],l,t_grid[layer0],omg2);
	else 
	y_initial_solid6by6(ys[0],ys[1],ys[2],l,t_grid[layer0],omg2);

/*
 * Cowling approximation
 */
	if( is_cowl == YES )
	for(j=region_start;j<=region_end;j++)
	{
		region=j; /* used in calc_a */
	if( sf == FLUID_GRID ){
		runge_driver(dum_y_grid,dum_dydt_grid,yf[0],&(t_grid[layer0]),n_grid,
		dzdr_fluid2by2cowl,2,eps,l,omg2,is_scale,is_func,&n_cross);
		for(i=0;i<2;i++) yf[0][i]=dum_y_grid[i][n_grid-1];

		if( is_func==YES )
		{
		for(n=0;n<n_grid;n++) y_grid[0][0][n+layer0] = dum_y_grid[0][n];
		for(n=0;n<n_grid;n++) y_grid[0][1][n+layer0] = dum_y_grid[1][n];
		for(n=0;n<n_grid;n++) dydt_grid[0][0][n+layer0] = dum_dydt_grid[0][n];
		for(n=0;n<n_grid;n++) dydt_grid[0][1][n+layer0] = dum_dydt_grid[1][n];
		}
		if(j+1 == msolid.nregion ) break;
		layer0 = region_info[j+1][LOWER_BOUND_SOLID];
		if( (sf=sfflag[layer0]) == SOLID_GRID )
			continue_eigenfunccowl(ys,yf,FLUID_TO_SOLID);
	}else{
		for(k=0;k<2;k++)
		{
			runge_driver(dum_y_grid,dum_dydt_grid,ys[k],&(t_grid[layer0]),n_grid,
			dzdr_solid4by4cowl,4,eps,l,omg2,is_scale,is_func,&n_cross);
			for(i=0;i<4;i++) ys[k][i]=dum_y_grid[i][n_grid-1];

			if( is_func==YES )
			for(i=0;i<4;i++)
			for(n=0;n<n_grid;n++){
				y_grid[k][i][n+layer0] = dum_y_grid[i][n];
				dydt_grid[k][i][n+layer0] = dum_dydt_grid[i][n];
			}
		}
		if(j+1 == msolid.nregion ) break;
		layer0 = region_info[j+1][LOWER_BOUND_SOLID];
		if( (sf=sfflag[layer0]) == FLUID_GRID )
			continue_eigenfunccowl(ys,yf,SOLID_TO_FLUID);
	}
	n_grid = region_info[j+1][UPPER_BOUND_SOLID]
		 - region_info[j+1][LOWER_BOUND_SOLID]+1;
	}
/*
 * No Cowling approximation
 */
	else
	for(j=region_start;j<=region_end;j++)
	{
		region=j; /* used in calc_a */
	if( sf == FLUID_GRID ){
		for(k=0;k<2;k++)
		{
			runge_driver(dum_y_grid,dum_dydt_grid,yf[k],&(t_grid[layer0]),n_grid,
			dzdr_fluid4by4,4,eps,l,omg2,is_scale,is_func,&n_cross);
			for(i=0;i<4;i++) yf[k][i]=dum_y_grid[i][n_grid-1];

			if( is_func==YES )
			{
			for(n=0;n<n_grid;n++) y_grid[k][0][n+layer0] = dum_y_grid[0][n];
			for(n=0;n<n_grid;n++) y_grid[k][1][n+layer0] = dum_y_grid[1][n];
			for(n=0;n<n_grid;n++) y_grid[k][4][n+layer0] = dum_y_grid[2][n];
			for(n=0;n<n_grid;n++) y_grid[k][5][n+layer0] = dum_y_grid[3][n];
			for(n=0;n<n_grid;n++) dydt_grid[k][0][n+layer0] = dum_dydt_grid[0][n];
			for(n=0;n<n_grid;n++) dydt_grid[k][1][n+layer0] = dum_dydt_grid[1][n];
			for(n=0;n<n_grid;n++) dydt_grid[k][4][n+layer0] = dum_dydt_grid[2][n];
			for(n=0;n<n_grid;n++) dydt_grid[k][5][n+layer0] = dum_dydt_grid[3][n];
			}
		}
		if(j+1 == msolid.nregion ) break;
		layer0 = region_info[j+1][LOWER_BOUND_SOLID];
		if( (sf=sfflag[layer0]) == SOLID_GRID )
			continue_eigenfunc(ys,yf,FLUID_TO_SOLID);
	}else{
		for(k=0;k<3;k++)
		{
			runge_driver(dum_y_grid,dum_dydt_grid,ys[k],&(t_grid[layer0]),n_grid,
			dzdr_solid6by6,6,eps,l,omg2,is_scale,is_func,&n_cross);
			for(i=0;i<6;i++) ys[k][i]=dum_y_grid[i][n_grid-1];

			if( is_func==YES )
			for(i=0;i<6;i++)
			for(n=0;n<n_grid;n++){
				y_grid[k][i][n+layer0] = dum_y_grid[i][n];
				dydt_grid[k][i][n+layer0] = dum_dydt_grid[i][n];
			}
		}
		if(j+1 == msolid.nregion ) break;
		layer0 = region_info[j+1][LOWER_BOUND_SOLID];
		if( (sf=sfflag[layer0]) == FLUID_GRID )
			continue_eigenfunc(ys,yf,SOLID_TO_FLUID);
	}
	n_grid = region_info[j+1][UPPER_BOUND_SOLID]
		 - region_info[j+1][LOWER_BOUND_SOLID]+1;
	}
	*layer1 = region_info[region_end][UPPER_BOUND_SOLID];
}
#define SMALL     1e-8
#define MAXSTEPS	100
#define     ABS(a)      ((a) > 0.0 ? (a) : -(a))
void zbrent_solid(best_omg,omg,omg1,omg20,is_disper,is_cowl,bc_value,bc_value1,y_grid,dydt_grid,
dum_y_grid,dum_dydt_grid,tol,ys0,yf0,sfflag,t_grid,n_grid,eps,l,layer0,
region_info,region0, region_end,bc_type,bc_coeff)
/* input
	omg,		angular frequency lower bound
	omg1,		angular frequency upper bound
	omg20,	normalization factor of (angular freq)^2
	is_disper,	if physical dispersion included ==YES, if not == NO
	is_cowl,	if Cowling approximation is used ==YES, if not == NO
	bc_value,	boundary_condition() at omg
	bc_value1,	boundary_condition() at omg1
	***y_grid,	dummy array used in driver_solid()
	***dydt_grid,dummy array used in driver_solid()
	**dum_y_grid,dummy array used in driver_solid()
	**dum_dydt_grid,dummy array used in driver_solid()
	tol,		required accuracy for zbrent() function
	ys0[3][6],	initial value of ys at t_grid[layer0]
	yf0[2][4],	initial value of yf at t_grid[layer0]
	sfflag,	FLIUD_GRID or SOLID_GRID of t_grid[]
	t_grid[],	grid values of t axis where y and dydt are evaluated
	n_grid,	number of grio points 0<=i < n_grid
	eps,		required accuracy for runge-kuttta integration
	l,		angular order
	layer0,	integration start from this layer
	region0,	layer0 is in this region
	region_end,	integrate up to this region.
	region_info[][], region information
	bc_type,	boudnary condition type
	bc_coeff,	boundary condition coefficient
  output
	best_omg,	estimated eigenfrequency
-------------------------------------------------------------------

	this routeine based on Numerical Recipes p.251
	chapter 9.3 Van Wijngaarden-dekker-brent Method
*/
	double *best_omg,omg,omg1,omg20,tol,***y_grid,***dydt_grid,t_grid[],eps;
	double **dum_y_grid,**dum_dydt_grid,ys0[3][6],yf0[2][4];
	double bc_value,bc_value1,bc_coeff[6];
	int l,n_grid,layer0,region_info[MAX_N_REGIONS][3],region0,region_end;
	int is_disper,is_cowl,bc_type;
	char *sfflag;
{
	double a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm,min1,min2;
	double ys[3][6],yf[2][4];
	int i,j,layer_end;
	extern struct model_solid msolid; /* defined in read_solid.c */
	void solid_driver();
	double boundary_condition();

	a = omg;
	b = omg1;

	fa = bc_value;
	fb = bc_value1;
/*
	solid_driver(y_grid,dydt_grid,ys,yf,&layer_end,sfflag,t_grid,n_grid,eps,l,a,omg20,
			 NO,NO,is_disper,is_cowl,layer0,region_info,region0,region_end,dum_y_grid,dum_dydt_grid);
	fa = boundary_condition(ys,yf,l,a,sfflag[layer_end],is_cowl,bc_type,bc_coeff);

	solid_driver(y_grid,dydt_grid,ys,yf,&layer_end,sfflag,t_grid,n_grid,eps,l,a,omg20,
			 NO,NO,is_disper,is_cowl,layer0,region_info,region0,region_end,dum_y_grid,dum_dydt_grid);
	fb = boundary_condition(ys,yf,l,b,sfflag[layer_end],is_cowl,bc_type,bc_coeff);
*/

	if( fb*fa > 0.0 )
	{
		FR,"zbrent(%s) %d: Root must be bracketed for zbrent.\n", __FILE__,__LINE__);
		exit(-1);
	}
	fc = fb;
	for(i=0;i<MAXSTEPS;i++)
	{
		if(fb*fc > 0.0)
		{
			c = a;
			fc = fa;
			e = d = b - a;
		}
		if( ABS(fc) < ABS(fb))
		{
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}
		tol1 = 2.0*SMALL*fabs(b)+0.5*tol;
		xm = 0.5*(c-b);
		if( ABS(xm) < tol1 || fb == 0.0 )
		{
			*best_omg = b;
			return;
		}
		if( ABS(e) >= tol1 && ABS(fa) > ABS(fb))
		{
			s = fb/fa;
			if( a == c)
			{
				p = 2.0*xm*s;
				q = 1.0 - s;
			}
			else
			{
				q = fa/fc;
				r = fb/fc;
				p = s*(2.0*xm*q*(q-r) - (b-a)*(r-1.0));
				q = (q-1.0)*(r-1.0)*(s-1.0);
			}
			if( p > 0.0 ) q = -q;
			p = ABS(p);
			min1 = 3.0*xm*q-fabs(tol1*q);
			min2 = fabs(e*q);
			if( 2.0*p < (min1 < min2 ? min1 : min2 ))
			{
				e =d;
				d=p/q;
			}
			else
			{
				d = xm;
				e = d;
			}
		}
		else
		{
			d = xm;
			e = d;
		}
		a = b;
		fa = fb;
		if( ABS(d) > tol1 ) 
			b += d;
		else
			b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));

		if( is_disper == YES ){
		set_dispersion(msolid.value[A_SOLID],msolid.value[C_SOLID],
		msolid.value[L_SOLID],msolid.value[N_SOLID],msolid.value[F_SOLID],
		msolid.value0[A_SOLID],msolid.value0[C_SOLID],
		msolid.value0[L_SOLID],msolid.value0[N_SOLID],msolid.value0[F_SOLID],
		msolid.value0[QK_SOLID],msolid.value0[QU_SOLID],
		2.0*PI/sqrt(b*omg20),0,221);/* <<====THIS IS AD HOC */
		if( sfflag[layer0]== FLUID_GRID )
			y_initial_fluid4by4(yf[0],yf[1],l,t_grid[layer0],b);
		else 
			y_initial_solid6by6(ys[0],ys[1],ys[2],l,t_grid[layer0],b);
		}else{
		if(sfflag[layer0] == FLUID_GRID)
			for(i=0;i<2;i++) for(j=0;j<4;j++) yf[i][j]=yf0[i][j];
		else
			for(i=0;i<3;i++) for(j=0;j<6;j++) ys[i][j]=ys0[i][j];
		}

		solid_driver(y_grid,dydt_grid,ys,yf,&layer_end,sfflag,t_grid,n_grid,eps,l,b,omg20,
			 NO,NO,is_disper,is_cowl,layer0,region_info,region0,region_end,dum_y_grid,dum_dydt_grid);
		fb = boundary_condition(ys,yf,l,b,sfflag[layer_end],is_cowl,bc_type,bc_coeff);
	}
	FR,"max step reached in zbrent_solid.\n");
}

void layer_divider(fine_y,fine_n_grid,table,max_dy,coarse_y,coarse_n_grid)
double *fine_y,max_dy,*coarse_y;
int *fine_n_grid,*table,coarse_n_grid;
/*
<<OUTPUT>>
fine_y[]: divided grid
fine_n_grid: # of grid points of fine_y[]
table[]: j=table[i] i=0,fin_n_grid,coarse_y[i] == fine_y[j]
<<INPUT>>
max_dy: maximum inteval of fine_y[] if max_dy <= 0.0 do not divide
coarse_y[]:	This grid is divided into fine_y[]
coarse_n_grid: # of grid points of coarse_y[]

if coarse_y[i+1]-coarse_y[i] > 1.5*max_dy, divide coarse_y[] into
small layers so that max(coarse_y[i+1]-coarse_y[i]) <1.5*max_dy
*/
{
	int i,counter;
	static int limit=0;

	if( max_dy > 0.0 )
	{
		counter=0;
		for(i=0;i<coarse_n_grid-1;i++)
		{
			table[i] = counter;
			fine_y[counter++] = coarse_y[i];limit++;
			while( coarse_y[i+1]-fine_y[counter-1] > 1.5*max_dy )
			{
				fine_y[counter] = max_dy+fine_y[counter-1];
				counter++;limit++;
			if( limit > GRID_NUMBER -1){
				FR,"# of divided layers (=%d) exceeds the limit.\n",
				limit);exit(-1);}
			}
		}
		fine_y[counter] = coarse_y[coarse_n_grid-1];
		table[coarse_n_grid-1] = counter;
		*fine_n_grid = ++counter;
	}else{
		for(i=0;i<coarse_n_grid;i++)
		{
			fine_y[i] = coarse_y[i];
			table[i] = i;
		}
		*fine_n_grid = coarse_n_grid;
		limit += coarse_n_grid;
	}
}
void test_vectors(ys,yf,sfflag,is_cowl)
double ys[6][3],yf[4][2];
int sfflag;
{
	double x1,x2,x3;

	if( is_cowl == YES)
	switch( sfflag ){
	case SOLID_GRID:
		x1 = vector_cos(ys[0],ys[1],4);
		printf("orthogonality of two solution vectors\n %le\n",x1);
		break;
	}
	else
	switch( sfflag ){
	case SOLID_GRID:
		x1 = vector_cos(ys[0],ys[1],6);
		x2 = vector_cos(ys[1],ys[2],6);
		x3 = vector_cos(ys[2],ys[0],6);
		printf("orthogonality of three solution vectors\n %le %le %le\n",
		x1,x2,x3);
		break;
	case FLUID_GRID:
		x1 = vector_cos(yf[0],yf[1],4);
		printf("orthogonality of two solution vectors\n %le\n",x1);
		break;
	}
}
