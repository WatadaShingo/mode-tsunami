/*
	miscellaneous functions to allocate and free vector and matrix variables.

	Feb. 15, 192
	Shingo Watada seismo Lab Caltech
	## NOTE ##
	I do not recommend to use clear() functions within mallocutil.c 
	From some experiments, it turned out that clear() takes huge memory space
	when used with mallocutil.c. I don't know why.
	## END NOTE Oct. 2, 1992 ##
*/
#include <stdio.h>
#include <stdlib.h>
#include "mallocutil.h"
#include "clear.h"
#define FR	fprintf(stderr

short int *malloc_vshortint(i_begin,i_end,mesg)
	int i_end,i_begin;
	char *mesg;
{
	short int *p;
	if( i_end <i_begin )
	{
		FR,"%s(malloc_vshort int): inputs error[%d %d].\n",mesg,i_begin,i_end);
		exit(-1);
	}
	if((p=(short int *)malloc((i_end-i_begin+1)*sizeof(short int))) == NULL )
	{
		FR,"%s(malloc_vshort int): cannot allocate.\n",mesg);
		exit(-1);
	}
/*	clear_vshortint(p,0,i_end - i_begin);*/
	return p - i_begin;
}
int *malloc_vint(i_begin,i_end,mesg)
	int i_end,i_begin;
	char *mesg;
{
	int *p;
	if( i_end <i_begin )
	{
		FR,"%s(malloc_vint): inputs error[%d %d].\n",mesg,i_begin,i_end);
		exit(-1);
	}
	if((p=(int *)malloc((i_end-i_begin+1)*sizeof(int))) == NULL )
	{
		FR,"%s(malloc_vint): cannot allocate.\n",mesg);
		exit(-1);
	}
/*	clear_vint(p,0,i_end - i_begin);*/
	return p - i_begin;
}
float *malloc_vfloat(i_begin,i_end,mesg)
	int i_end,i_begin;
	char *mesg;
{
	float *p;
	if( i_end <i_begin )
	{
		FR,"%s(malloc_vfloat): inputs error[%d %d].\n",mesg,i_begin,i_end);
		exit(-1);
	}
	if((p=(float *)malloc((i_end-i_begin+1)*sizeof(float))) == NULL )
	{
		FR,"%s(malloc_vfloat): cannot allocate.\n",mesg);
		exit(-1);
	}
/*	clear_vfloat(p,0,i_end - i_begin);*/
	return p - i_begin;
}
double *malloc_vdouble(i_begin,i_end,mesg)
	int i_end,i_begin;
	char *mesg;
{
	double *p;
	if( i_end <i_begin )
	{
		FR,"%s(malloc_vdouble): inputs error[%d %d].\n",mesg,i_begin,i_end);
		exit(-1);
	}
	if((p=(double *)malloc((i_end-i_begin+1)*sizeof(double))) == NULL )
	{
		FR,"%s(malloc_vdouble): cannot allocate.\n",mesg);
		exit(-1);
	}
/*	clear_vdouble(p,0,i_end - i_begin);*/
	return p - i_begin;
}
complex *malloc_vcomplex(i_begin,i_end,mesg)
	int i_end,i_begin;
	char *mesg;
{
	complex *p;
	if( i_end <i_begin )
	{
		FR,"%s(malloc_vcomplex): inputs error[%d %d].\n",mesg,i_begin,i_end);
		exit(-1);
	}
	if((p=(complex *)malloc((i_end-i_begin+1)*sizeof(complex))) == NULL )
	{
		FR,"%s(malloc_vcomplex): cannot allocate.\n",mesg);
		exit(-1);
	}
/*	clear_vcomplex(p,0,i_end - i_begin);*/
	return p - i_begin;
}
dcomplex *malloc_vdcomplex(i_begin,i_end,mesg)
	int i_end,i_begin;
	char *mesg;
{
	dcomplex *p;
	if( i_end <i_begin )
	{
		FR,"%s(malloc_vdcomplex): inputs error[%d %d].\n",mesg,i_begin,i_end);
		exit(-1);
	}
	if((p=(dcomplex *)malloc((i_end-i_begin+1)*sizeof(dcomplex))) == NULL )
	{
		FR,"%s(malloc_vdcomplex): cannot allocate.\n",mesg);
		exit(-1);
	}
/*	clear_vdcomplex(p,0,i_end - i_begin);*/
	return p - i_begin;
}
short int **malloc_mshortint(i_begin,i_end,j_begin,j_end,mesg)
	int i_end,i_begin,j_end,j_begin;
	char *mesg;
{
	int i;
	short int **p;
	if( i_end < i_begin || j_end < j_begin )
	{
		FR,"%s(malloc_mshortint): inputs error[%d %d][%d %d].\n",
		mesg,i_begin,i_end,j_begin,j_end);
		exit(-1);
	}
	if((p=(short int**)malloc((i_end-i_begin+1)*sizeof(short int *))) == NULL )
	{
		FR,"%s(malloc_mshortint 0):cannot callocate.\n",mesg);
		exit(-1);
	}
	p -= i_begin;
	for(i=i_begin;i<=i_end;i++)
	{
		if((p[i]=(short int*)malloc((j_end-j_begin+1)*sizeof(short int))) == NULL )
		{
			FR,"%s(malloc_mshortint 1):cannot callocate.\n",mesg);
			exit(-1);
		}
		p[i] -= j_begin;
	}
/*	clear_mshortint(p,0,i_end - i_begin,0,j_end - j_begin);*/
	return p;
}

int **malloc_mint(i_begin,i_end,j_begin,j_end,mesg)
	int i_end,i_begin,j_end,j_begin;
	char *mesg;
{
	int i;
	int **p;
	if( i_end < i_begin || j_end < j_begin )
	{
		FR,"%s(malloc_mint): inputs error[%d %d][%d %d].\n",
		mesg,i_begin,i_end,j_begin,j_end);
		exit(-1);
	}
	if((p=(int**)malloc((i_end-i_begin+1)*sizeof(int *))) == NULL )
	{
		FR,"%s(malloc_mint 0):cannot callocate.\n",mesg);
		exit(-1);
	}
	p -= i_begin;
	for(i=i_begin;i<=i_end;i++)
	{
		if((p[i]=(int*)malloc((j_end-j_begin+1)*sizeof(int))) == NULL )
		{
			FR,"%s(malloc_mint 1):cannot callocate.\n",mesg);
			exit(-1);
		}
		p[i] -= j_begin;
	}
/*	clear_mint(p,0,i_end - i_begin,0,j_end - j_begin);*/
	return p;
}
float **malloc_mfloat(i_begin,i_end,j_begin,j_end,mesg)
	int i_end,i_begin,j_end,j_begin;
	char *mesg;
{
	int i;
	float **p;
	if( i_end < i_begin || j_end < j_begin )
	{
		FR,"%s(malloc_mfloat): inputs error[%d %d][%d %d].\n",
		mesg,i_begin,i_end,j_begin,j_end);
		exit(-1);
	}
	if((p=(float**)malloc((i_end-i_begin+1)*sizeof(float *))) == NULL )
	{
		FR,"%s(malloc_mfloat 0):cannot callocate.\n",mesg);
		exit(-1);
	}
	p -= i_begin;
	for(i=i_begin;i<=i_end;i++)
	{
		if((p[i]=(float*)malloc((j_end-j_begin+1)*sizeof(float))) == NULL )
		{
			FR,"%s(malloc_mfloat 1):cannot callocate.\n",mesg);
			exit(-1);
		}
		p[i] -= j_begin;
	}
/*	clear_mfloat(p,0,i_end - i_begin,0,j_end - j_begin);*/
	return p;
}
double **malloc_mdouble(i_begin,i_end,j_begin,j_end,mesg)
	int i_end,i_begin,j_end,j_begin;
	char *mesg;
{
	int i;
	double **p;
	if( i_end < i_begin || j_end < j_begin )
	{
		FR,"%s(malloc_mdouble): inputs error[%d %d][%d %d].\n",
		mesg,i_begin,i_end,j_begin,j_end);
		exit(-1);
	}
	if((p=(double**)malloc((i_end-i_begin+1)*sizeof(double *))) == NULL )
	{
		FR,"%s(malloc_mdouble 0):cannot callocate.\n",mesg);
		exit(-1);
	}
	p -= i_begin;
	for(i=i_begin;i<=i_end;i++)
	{
		if((p[i]=(double*)malloc((j_end-j_begin+1)*sizeof(double))) == NULL )
		{
			FR,"%s(malloc_mdouble 1):cannot callocate.\n",mesg);
			exit(-1);
		}
		p[i] -= j_begin;
	}
/*	clear_mdouble(p,0,i_end - i_begin,0,j_end - j_begin);*/
	return p;
}
complex **malloc_mcomplex(i_begin,i_end,j_begin,j_end,mesg)
	int i_end,i_begin,j_end,j_begin;
	char *mesg;
{
	int i;
	complex **p;
	if( i_end < i_begin || j_end < j_begin )
	{
		FR,"%s(malloc_mcomplex): inputs error[%d %d][%d %d].\n",
		mesg,i_begin,i_end,j_begin,j_end);
		exit(-1);
	}
	if((p=(complex**)malloc((i_end-i_begin+1)*sizeof(complex *))) == NULL )
	{
		FR,"%s(malloc_mcomplex 0):cannot callocate.\n",mesg);
		exit(-1);
	}
	p -= i_begin;
	for(i=i_begin;i<=i_end;i++)
	{
		if((p[i]=(complex*)malloc((j_end-j_begin+1)*sizeof(complex))) == NULL )
		{
			FR,"%s(malloc_mcomplex 1):cannot callocate.\n",mesg);
			exit(-1);
		}
		p[i] -= j_begin;
	}
/*	clear_mcomplex(p,0,i_end - i_begin,0,j_end - j_begin);*/
	return p;
}
dcomplex **malloc_mdcomplex(i_begin,i_end,j_begin,j_end,mesg)
	int i_end,i_begin,j_end,j_begin;
	char *mesg;
{
	int i;
	dcomplex **p;
	if( i_end < i_begin || j_end < j_begin )
	{
		FR,"%s(malloc_mdcomplex): inputs error[%d %d][%d %d].\n",
		mesg,i_begin,i_end,j_begin,j_end);
		exit(-1);
	}
	if((p=(dcomplex**)malloc((i_end-i_begin+1)*sizeof(dcomplex *))) == NULL )
	{
		FR,"%s(malloc_mdcomplex 0):cannot callocate.\n",mesg);
		exit(-1);
	}
	p -= i_begin;
	for(i=i_begin;i<=i_end;i++)
	{
		if((p[i]=(dcomplex*)malloc((j_end-j_begin+1)*sizeof(dcomplex))) == NULL )
		{
			FR,"%s(malloc_mdcomplex 1):cannot callocate.\n",mesg);
			exit(-1);
		}
		p[i] -= j_begin;
	}
	clear_mdcomplex(p,i_begin,i_end,j_begin,j_end);
	return p;
}
double ***malloc_m3double(i_begin,i_end,j_begin,j_end,k_begin,k_end,mesg)
	int i_end,i_begin,j_end,j_begin,k_end,k_begin;
	char *mesg;
{
	int i,j;
	double ***p;
	if( i_end < i_begin || j_end < j_begin || k_end < k_begin)
	{
		FR,"%s(malloc_m3double): inputs error[%d %d][%d %d][%d %d].\n",
		mesg,i_begin,i_end,j_begin,j_end,k_begin,k_end);
		exit(-1);
	}
	if((p=(double***)malloc((i_end-i_begin+1)*sizeof(double **))) == NULL )
	{
		FR,"%s(malloc_m3double 0):cannot callocate.\n",mesg);
		exit(-1);
	}
	p -= i_begin;
	for(i=i_begin;i<=i_end;i++)
	{
		if((p[i]=(double**)malloc((j_end-j_begin+1)*sizeof(double *))) == NULL )
		{
			FR,"%s(malloc_m3double 1):cannot callocate.\n",mesg);
			exit(-1);
		}
		p[i] -= j_begin;
	}
	for(i=i_begin;i<=i_end;i++)
	{
		for(j=j_begin;j<=j_end;j++)
		{
			if((p[i][j]=(double*)malloc((k_end-k_begin+1)*sizeof(double ))) == NULL )
			{
				FR,"%s(malloc_m3double 2):cannot callocate.\n",mesg);
				exit(-1);
			}
			p[i][j] -= k_begin;
		}
	}
	return p;
}
dcomplex ***malloc_m3dcomplex(i_begin,i_end,j_begin,j_end,k_begin,k_end,mesg)
	int i_end,i_begin,j_end,j_begin,k_end,k_begin;
	char *mesg;
{
	int i,j;
	dcomplex ***p;
	if( i_end < i_begin || j_end < j_begin || k_end < k_begin)
	{
		FR,"%s(malloc_m3dcomplex): inputs error[%d %d][%d %d][%d %d].\n",
		mesg,i_begin,i_end,j_begin,j_end,k_begin,k_end);
		exit(-1);
	}
	if((p=(dcomplex***)malloc((i_end-i_begin+1)*sizeof(dcomplex **))) == NULL )
	{
		FR,"%s(malloc_m3dcomplex 0):cannot callocate.\n",mesg);
		exit(-1);
	}
	p -= i_begin;
	for(i=i_begin;i<=i_end;i++)
	{
		if((p[i]=(dcomplex**)malloc((j_end-j_begin+1)*sizeof(dcomplex *))) == NULL )
		{
			FR,"%s(malloc_m3dcomplex 1):cannot callocate.\n",mesg);
			exit(-1);
		}
		p[i] -= j_begin;
	}
	for(i=i_begin;i<=i_end;i++)
	{
		for(j=j_begin;j<=j_end;j++)
		{
			if((p[i][j]=(dcomplex*)malloc((k_end-k_begin+1)*sizeof(dcomplex ))) == NULL )
			{
				FR,"%s(malloc_m3dcomplex 2):cannot callocate.\n",mesg);
				exit(-1);
			}
			p[i][j] -= k_begin;
		}
	}
/*	clear_m3dcomplex(p,i_begin,i_end,j_begin,j_end,k_begin,k_end); */
	return p;
}
void free_vshortint(p,i_begin)
	short int *p;
	int i_begin;
{
	free(p+i_begin);
}
void free_vint(p,i_begin)
	int *p;
	int i_begin;
{
	free(p+i_begin);
}
void free_vfloat(p,i_begin)
	float *p;
	int i_begin;
{
	free(p+i_begin);
}
void free_vdouble(p,i_begin)
	double *p;
	int i_begin;
{
	free(p+i_begin);
}
void free_vcomplex(p,i_begin)
	complex *p;
	int i_begin;
{
	free(p+i_begin);
}
void free_vdcomplex(p,i_begin)
	dcomplex *p;
	int i_begin;
{
	free(p+i_begin);
}
void free_mshortint(p,i_begin,i_end,j_begin)
	short int **p;
	int i_begin,i_end,j_begin;
{
	int i;
	for(i=i_begin;i<=i_end;i++)
		free(p[i]+j_begin);
	free(p+i_begin);
}
void free_mint(p,i_begin,i_end,j_begin)
	int **p;
	int i_begin,i_end,j_begin;
{
	int i;
	for(i=i_begin;i<=i_end;i++)
		free(p[i]+j_begin);
	free(p+i_begin);
}
void free_mfloat(p,i_begin,i_end,j_begin)
	float **p;
	int i_begin,i_end,j_begin;
{
	int i;
	for(i=i_begin;i<=i_end;i++)
		free(p[i]+j_begin);
	free(p+i_begin);
}
void free_mdouble(p,i_begin,i_end,j_begin)
	double **p;
	int i_begin,i_end,j_begin;
{
	int i;
	for(i=i_begin;i<=i_end;i++)
		free(p[i]+j_begin);
	free(p+i_begin);
}
void free_mcomplex(p,i_begin,i_end,j_begin)
	complex **p;
	int i_begin,i_end,j_begin;
{
	int i;
	for(i=i_begin;i<=i_end;i++)
		free(p[i]+j_begin);
	free(p+i_begin);
}
void free_mdcomplex(p,i_begin,i_end,j_begin)
	dcomplex **p;
	int i_begin,i_end,j_begin;
{
	int i;
	for(i=i_begin;i<=i_end;i++)
		free(p[i]+j_begin);
	free(p+i_begin);
}
void free_m3double(p,i_begin,i_end,j_begin,j_end,k_begin)
	double ***p;
	int i_begin,i_end,j_begin,j_end,k_begin;
{
	int i,j;
	for(i=i_begin;i<=i_end;i++)
	{
		for(j=j_begin;j<=j_end;j++)
			free(p[i][j]+k_begin);
		free(p[i]+j_begin);
	}
	free(p+i_begin);
}
void free_m3dcomplex(p,i_begin,i_end,j_begin,j_end,k_begin)
	dcomplex ***p;
	int i_begin,i_end,j_begin,j_end,k_begin;
{
	int i,j;
	for(i=i_begin;i<=i_end;i++)
	{
		for(j=j_begin;j<=j_end;j++)
			free(p[i][j]+k_begin);
		free(p[i]+j_begin);
	}
	free(p+i_begin);
}
