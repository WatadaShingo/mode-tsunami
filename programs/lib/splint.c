#include <stdio.h>
void splint(xa,ya,y2a,n,x,y,dy,ddy)
float xa[],ya[],y2a[],x,*y,*dy,*ddy;
int n;
{
	int klo,khi;
	static int k=0;
	float h,b,a;


	klo=0;
	khi=n-1;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1 ;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) 
	{
		fprintf(stderr,"Bad XA input to routine SPLINT\n");
		exit(-1);
	}
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
	*dy=(ya[khi]-ya[klo])/h-((3.0*a*a-1.0)*y2a[klo]-(3.0*b*b-1.0)*y2a[khi])*h/6.0;
	*ddy=a*y2a[klo]+b*y2a[khi];
}
void splintd(xa,ya,y2a,n,x,y,dy,ddy)
double xa[],ya[],y2a[],x,*y,*dy,*ddy;
int n;
{
	int klo,khi;
	static int k=0;
	double h,b,a;


	klo=0;
	khi=n-1;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1 ;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) 
	{
		fprintf(stderr,"Bad XA input to routine SPLINT\n");
		exit(-1);
	}
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
	*dy=(ya[khi]-ya[klo])/h-((3.0*a*a-1.0)*y2a[klo]-(3.0*b*b-1.0)*y2a[khi])*h/6.0;
	*ddy=a*y2a[klo]+b*y2a[khi];
}
