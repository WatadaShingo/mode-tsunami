#include "clear.h"
#include "clib.h"

void clear_vshortint(p,i_begin,i_end)
	short int *p;
	int i_begin,i_end;
{
	while(i_end >= i_begin ) p[i_begin++]=0;
}
void clear_mshortint(p,i_begin,i_end,j_begin,j_end)
	short int **p;
	int i_begin,i_end,j_begin,j_end;
{
	int i,j;
	for(i=i_begin;i<=i_end;i++)
	for(j=j_begin;j<=j_end;j++)
		p[i][j] = 0; 
}
void clear_vint(p,i_begin,i_end)
	int *p;
	int i_begin,i_end;
{
	while(i_end >= i_begin ) p[i_begin++]=0;
}
void clear_mint(p,i_begin,i_end,j_begin,j_end)
	int **p;
	int i_begin,i_end,j_begin,j_end;
{
	int i,j;
	for(i=i_begin;i<=i_end;i++)
	for(j=j_begin;j<=j_end;j++)
		p[i][j] = 0; 
}
void clear_vfloat(p,i_begin,i_end)
	float *p;
	int i_begin,i_end;
{
	while(i_end >= i_begin ) p[i_begin++]=0e0;
}
void clear_mfloat(p,i_begin,i_end,j_begin,j_end)
	float **p;
	int i_begin,i_end,j_begin,j_end;
{
	int i,j;
	for(i=i_begin;i<=i_end;i++)
	for(j=j_begin;j<=j_end;j++)
		p[i][j] = 0e0; 
}
void clear_vdouble(p,i_begin,i_end)
	double *p;
	int i_begin,i_end;
{
	while(i_end >= i_begin ) p[i_begin++]=0e0;
}
void clear_mdouble(p,i_begin,i_end,j_begin,j_end)
	double **p;
	int i_begin,i_end,j_begin,j_end;
{
	int i,j;
	for(i=i_begin;i<=i_end;i++)
	for(j=j_begin;j<=j_end;j++)
		p[i][j] = 0e0; 
}
void clear_vcomplex(p,i_begin,i_end)
	complex *p;
	int i_begin,i_end;
{
	while(i_end >= i_begin ) 
		p[i_begin].re=
		p[i_begin++].im=0e0;
	/* ### NOTE ### do not do 
		p[i_begin].re=
		p[i_begin++].re=0e0;
		because some compilers initialize p[i_begin].re with
		the minumum value( i.e. .xxxxxE-310) greater than zero.
		Sept 5 1992
	### END NOTE ### */
}
void clear_mcomplex(p,i_begin,i_end,j_begin,j_end)
	complex **p;
	int i_begin,i_end,j_begin,j_end;
{
	int i,j;
	for(i=i_begin;i<=i_end;i++)
	for(j=j_begin;j<=j_end;j++)
		p[i][j].re=
		p[i][j].im=0e0;
}
void clear_vdcomplex(p,i_begin,i_end)
	dcomplex *p;
	int i_begin,i_end;
{
	while(i_end >= i_begin )
		p[i_begin].re=
		p[i_begin++].im=0e0;
}
void clear_mdcomplex(p,i_begin,i_end,j_begin,j_end)
	dcomplex **p;
	int i_begin,i_end,j_begin,j_end;
{
	int i,j;
	for(i=i_begin;i<=i_end;i++)
	for(j=j_begin;j<=j_end;j++)
		p[i][j].re=
		p[i][j].im=0e0;
}
void clear_m3dcomplex(p,i_begin,i_end,j_begin,j_end,k_begin,k_end)
	dcomplex ***p;
	int i_begin,i_end,j_begin,j_end,k_begin,k_end;
{
	int i,j,k;
	for(i=i_begin;i<=i_end;i++)
	for(j=j_begin;j<=j_end;j++)
	for(k=k_begin;k<=k_end;k++)
		p[i][j][k].re=
		p[i][j][k].im=0e0;
}
