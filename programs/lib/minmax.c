#include "minmax.h"

int set_imax(data,ndata)
	int *data;
	int ndata;
{
	int i=0;
	int dum=*data;
	for(;i<ndata;data++,i++) dum = MAX(dum,*data);
	data-=ndata;
	return(dum);
}
int set_imin(data,ndata)
	int *data;
	int ndata;
{
	int i=0;
	int dum=*data;
	for(;i<ndata;data++,i++) dum = MIN(dum,*data);
	data-=ndata;
	return(dum);
}
float set_max(data,ndata)
	float *data;
	int ndata;
{
	int i=0;
	float dum=*data;
	for(;i<ndata;data++,i++) dum = MAX(dum,*data);
	data-=ndata;
	return(dum);
}
float set_min(data,ndata)
	float *data;
	int ndata;
{
	int i=0;
	float dum=*data;
	for(;i<ndata;data++,i++) dum = MIN(dum,*data);
	data-=ndata;
	return(dum);
}

double set_dmax(data,ndata)
	double *data;
	int ndata;
{
	int i=0;
	double dum=*data;
	for(;i<ndata;data++,i++) dum = MAX(dum,*data);
	data-=ndata;
	return(dum);
}
double set_dmin(data,ndata)
	double *data;
	int ndata;
{
	int i=0;
	double dum=*data;
	for(;i<ndata;data++,i++) dum = MIN(dum,*data);
	data-=ndata;
	return(dum);
}
float get_fmean(data,ndata)
	float *data;
	int ndata;
{
	int i;
	double sum;
	for(sum=0.0,i=0;i<ndata;i++) sum+=data[i];
	return ((float) sum/(double)ndata);
}
double get_dmean(data,ndata)
	double *data;
	int ndata;
{
	int i;
	double sum;
	for(sum=0.0,i=0;i<ndata;i++) sum+=data[i];
	return (sum/(double)ndata);
}

void swapi(i1,i2)
	int *i1,*i2;
{
	int dum;
	dum = *i1;
	*i1 = *i2;
	*i2 = dum;
}
void swapf(f1,f2)
	float *f1,*f2;
{
	float dum;
	dum = *f1;
	*f1 = *f2;
	*f2 = dum;
}
void sort_f(data,ndata)
	float *data;
	int ndata;
{
	int i,flag;
	float dum;

	if ( ndata == 1 ) return;
	flag = 1;
	while(flag == 1)
	{
	flag = 0;
	for(i=0;i<ndata -1;i++)
	{
		if( data[i] > data[i+1] )
		{
			flag = 1;
			dum = data[i];data[i]=data[i+1];data[i+1]=dum;
		}
	}	
	}
}
void sort_d(data,ndata)
	double *data;
	int ndata;
{
	int i,flag;
	double dum;

	if ( ndata == 1 ) return;
	flag = 1;
	while(flag == 1)
	{
	flag = 0;
	for(i=0;i<ndata -1;i++)
	{
		if( data[i] > data[i+1] )
		{
			flag = 1;
			dum = data[i];data[i]=data[i+1];data[i+1]=dum;
		}
	}	
	}
}
int search_nearf(data,ndata,target,value)
	float *data,target,*value;
	int ndata;
/*
 * return the nearest value of data to the target and # of the data 
 * assume that data is sorted by sort_f()
 */
{
	int i;

	if( target < data[0]) return 0;
	if( target > data[ndata-1]) return ndata-1;
	for(i=0;i<ndata-1;i++)
		if( !OUT(target,data[i],data[i+1]) ) break;

	if( data[i+1]-target < target - data[i] ) i++;

	*value = data[i];
	return i; 
}
int search_neard(data,ndata,target,value)
	double *data,target,*value;
	int ndata;
/*
 * return the nearest value of data to the target and # of the data 
 * assume that data is sorted by sort_d()
 */
{
	int i;

	if( target < data[0]) return 0;
	if( target > data[ndata-1]) return ndata-1;
	for(i=0;i<ndata-1;i++)
		if( !OUT(target,data[i],data[i+1]) ) break;

	if( data[i+1]-target < target - data[i] ) i++;

	*value = data[i];
	return i; 
}
