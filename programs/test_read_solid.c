#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "m_solid.h"
#include "read_solid.h"
#include "mallocutil.h"
#include "ray_sub.h"
#include "calc_a.h"
#include "search.h"
#define 	FR	fprintf(stderr

extern struct model_solid msolid;

void local_analysis();
void modal_analysis();

int region; /* referred in calc_a.c */

main()
{
	int i,nregion,key;
	int level,turning_level,el;
	double r,max_r,min_r;
	double *r_grid,*v_grid,omg,eps;
	struct model_funcsolid model;
	char modelfilename[100];

	do{
	puts(" print out earth model          --->1");
	puts(" find ray turning point of mode --->2");
	puts(" local analysis of mode         --->3");
	puts(" examine modal behavior         --->4");
	scanf("%d",&key);
	}
	while( key < 1 || key > 4);

      printf("model filename >>");
	scanf("%s",modelfilename);
	read_solid(modelfilename);

	if( key == 1)
	{
	for(nregion=0;nregion<msolid.nregion;nregion++)
	{
		min_r = msolid.value[RAD_SOLID][msolid.region[nregion][LOWER_BOUND_SOLID]];
		max_r = msolid.value[RAD_SOLID][msolid.region[nregion][UPPER_BOUND_SOLID]];
		for(i=0,r=min_r;r<=max_r;r=min_r+10.0e3*i)
		{
			find_modelsolid(&model,r,nregion);
			printf("%7.2f %9.4e %9.4e %9.4e %9.4e %9.4e %9.4e %9.4e %9.4e\n",
			model.r/1e3, model.rho,model.a,model.c,model.l,model.n,model.f,model.g,model.drhodr/model.rho);
			i++;
		}
		find_modelsolid(&model,max_r,nregion);
		printf("%7.2f %9.4e %9.4e %9.4e %9.4e %9.4e %9.4e %9.4e %9.4e\n",
		model.r/1e3, model.rho,model.a,model.c,model.l,model.n,model.f,model.g,model.drhodr/model.rho);
	}
	}
	else if (key == 2 ) 
	{
		while(1)
		{
		printf("angular order of mode l(to end -999) = ");
		scanf("%d",&el); printf("\n");
		if( el < 0 ) exit(0);
		printf("period of mode in sec = ");
		scanf("%lf",&omg); printf("\n");
		printf("tolerance  level eps  = ");
		scanf("%lf",&eps); printf("\n");

		omg = 2*3.141592653589793/omg;

		r_grid = malloc_vdouble(0,msolid.nmodel-1,"main r_grid");
		v_grid = malloc_vdouble(0,msolid.nmodel-1,"main v_grid");

		for(i=0;i<msolid.nmodel;i++)
		{
			r_grid[i] = (double)msolid.value[RAD_SOLID][i];
			v_grid[i] = 
			(msolid.value[L_SOLID][i]==0.0 ? msolid.value[VPV_SOLID][i]:
								   msolid.value[VSV_SOLID][i]);
		}
		evanescent_level(&level,&turning_level,r_grid,v_grid,msolid.nmodel,omg,el,eps);
		printf("turning level:r[%d]=%9.4e, evanescent level:r[%d]=%9.4e\n",
		turning_level,r_grid[turning_level],level,r_grid[level]);
		}
	}
	else if (key == 3 )
		local_analysis(modelfilename);
	else if ( key == 4 );
		modal_analysis(modelfilename);
}
void local_analysis(modelfilename)
char *modelfilename;
{
	double best_omg,best_omg2,omg,old_omg,omg2,old_omg2,omg_start,omg_end,domg;
	double r; 
	double best_gamma, gamma, old_gamma,tol,k;
	int l,key,layer_top;
	struct model_funcsolid m;
	FILE *fp;
	char filename[124];
	int find_flag;
	double omg_acou,omg_grav,omg_dum;
	
	fp = (FILE *)NULL;
	printf("tolerance  level eps  (recommend 1e-15) = ");
	scanf("%lf",&tol); printf("\n");
	printf("input layer(0-%d) of top boundary (0 to set default %d)>>",msolid.nmodel-1,msolid.nmodel-1);
	printf("\n");
	scanf("%d",&layer_top); if(layer_top == 0 ) layer_top = msolid.nmodel-1;
	printf("1/(top bondary density scale height) k \n");
	printf("  if k<0 isothermal continuation with scalehight at the top.\n");
	printf("e.g.  homogeneous continuation k = 0.0\n");
	scanf("%lf",&k); printf("\n");
	r = msolid.value[RAD_SOLID][layer_top];
	set_regionsolid(&region,r,UPPER_REGION);
	find_modelsolid(&m,r,region);

	puts("print out boundary of gamma 1:");
	puts("do nothing                  0:");
	do{
		scanf("%d",&key);
		printf("\n");
		if( key ==1)
		{
			printf("file name of gamma info output >>");
			scanf("%s",filename);
			printf("\n");
			if((fp=fopen(filename,"w")) == NULL)
			{
				FR,"local_analysis: cannot creat file %s\n",
				filename);exit(-1);
			}
			fprintf(fp,"model : %s\n",modelfilename);
			fprintf(fp,"local analysis at r[%d]=%10.4e\n",layer_top,r);
			fprintf(fp,"omg<omg1 or omg>om2 evanescent mode, omg1<omg<omg2 propagating mode\n");
			fprintf(fp,"   l    omg1   <  omg2\n");
		}
	}while(key <0 || key >2);
	
	printf("input l, starting omg, ending omg, domg  to stop l=-9999>>");
	printf("\n");


	while( (key=scanf("%d %lf %lf %lf",&l,&omg_start,&omg_end,&domg))!=EOF)
	{
		if( l <-1000 )break;
		if( key !=4 ){FR,"input error!\n");exit(0);}
	find_flag = 0;omg_acou=omg_grav=0.0;
	find_omg_bound(&find_flag,&omg_acou,&omg_grav,omg_start,omg_end,r,l,&m,tol);
	printf("acoustic boundary at %10.4e\n",omg_acou);
	printf("gravity  boundary at %10.4e\n",omg_grav);
	if(fp){ fprintf(fp,"%4d %10.4e %10.4e\n",l,omg_grav,omg_acou);}
/*
 * this part is omitted by using find_omg_bound() function.
 * 
		if( fp ){ 
			fprintf(fp,"\n");
			fprintf(fp,"%4d",l);
		}
		omg = old_omg = omg_start;
		old_omg2 = old_omg*old_omg;
		calc_fluid2by2_character(&old_gamma,r,l,old_omg2,k,m);

		printf("l=%d omg=%9.7lf gamma=%10.4e\n",l,omg,old_gamma);
		while(omg < omg_end)
		{
			omg += domg;
			omg2 = omg*omg;
			calc_fluid2by2_character(&gamma,r,l,omg2,k,m);
			printf("l=%d omg=%9.7lf gamma=%10.4e\n",l,omg,gamma);

			if( gamma *old_gamma  <0.0)
			{
				search_omg(&best_omg2,omg2,old_omg2,gamma,old_gamma,r,l,k,&m,tol);
				calc_fluid2by2_character(&best_gamma,r,l,best_omg2,k,m);
				best_omg = sqrt(best_omg2);
			printf("found   omg=%9.7lf gamma=%10.4e\n",best_omg,best_gamma);
			if(fp) fprintf(fp," %10.4e",best_omg);
			}
			old_omg2 = omg2;
			old_gamma = gamma;
		}
*/
	}
	if(fp) fprintf(fp,"\n");
}
void modal_analysis(modelfilename)
char *modelfilename;
{
	int i,key,l,start_grid,counter;
	struct model_funcsolid model;
	double gamma[SAMPLE_SOLID],r[SAMPLE_SOLID],omg,omg2;
	FILE *fp;
	char filename[124],file_prefix[124];

	puts("print out results of gamma 1:");
	puts("do nothing                 0:");
	do{
		scanf("%d",&key);
		printf("\n");
		if( key ==1)
		{
			printf("file name prefix of gamma info output >>");
			scanf("%s",file_prefix);
			printf("\n");
		}
		else fp = stderr;
	}while(key <0 || key >2);
	
	counter =0;
	while(1){
		printf("mode anguler order (to quit -1) and angular freq =");
		scanf("%d %lf",&l,&omg);
		if( l < 0 ) break;
		omg2 = omg*omg;

		sprintf(filename,"%s.%ds%d",file_prefix,counter,l);
		if((fp=fopen(filename,"w")) == NULL)
		{
			FR,"local_analysis: cannot creat file %s\n",
			filename);exit(-1);
		}
		fprintf(fp,"%d s %d %.2f %e\n",counter,l,2*3.1415926535/omg,omg);
		fprintf(fp,"%d 1\n",msolid.region[13][NUMBER_OF_GRIDS]);
		fprintf(fp,"njunk\n");
		start_grid = msolid.region[13][LOWER_BOUND_SOLID];
		for(i=0;i<msolid.region[13][NUMBER_OF_GRIDS];i++)
		{
			r[i] = msolid.value[RAD_SOLID][start_grid+i];
			set_regionsolid(&region,r[i],UPPER_REGION);
			find_modelsolid(&model,r[i],region);
			calc_fluid2by2_character(&(gamma[i]),r[i],l,omg2,-1.0,model);
		}
		for(i=0;i<msolid.region[13][NUMBER_OF_GRIDS];i++)
		{
			fprintf(fp,"%10.4e %+10.4e\n",r[i]/1e3,gamma[i]);
		}
		counter ++;
		fclose(fp);
	}
}
