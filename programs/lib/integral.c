/*
 * integral by simpson rule
 * qsimp() and trapzd() from Numerical Recipes C
 *
	Shingo Watada Jan. 4, 1991 Seismo Lab, Caltech
 */
#include <stdio.h>
#include <math.h>
#include "clib.h"


#define EPS 1.0e-6
/* #define JMAX 20 original */
#define JMAX 15 

double qsimp(func,a,b)
float a,b;
double (*func)();
{
	int j;
	double s,st,ost,os,trapzd();

	ost = os = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		st=trapzd(func,(double)a,(double)b,j);
		s=(4.0*st-ost)/3.0;
		if (fabs(s-os) <= EPS*fabs(os) || fabs(s)<1e-35 ) return s;
		os=s;
		ost=st;
	}
	fprintf(stderr,"qsimp: Too many steps in routine QSIMP val=%15.8e\n",s);
	return s;
}

dcomplex dcqsimp(func,a,b)
float a,b;
dcomplex (*func)();
{
	int j;
	dcomplex s,st,ost,os,dctrapzd();

	ost.re = ost.im = os.re = os.im =  -1.0e30;
	 
	for (j=1;j<=JMAX;j++) {
		st=dctrapzd(func,(double)a,(double)b,j);
		s.re=(4.0*st.re-ost.re)/3.0;
		s.im=(4.0*st.im-ost.im)/3.0;
		if( (fabs(s.re-os.re) <= EPS*fabs(os.re) || fabs(s.re)<1e-35 )&&
		    (fabs(s.im-os.im) <= EPS*fabs(os.im) || fabs(s.im)<1e-35 ) )
		return s;
		os=s;
		ost=st;
	}
	fprintf(stderr,"dcqsimp: Too many steps in routine DCQSIMP val=(%13.6e,%13.6e)\n",s.re,s.im);
	return s;
}

double qtrap(func,a,b,n)
float a,b;
double (*func)();
{
	int i;
	double s,dx;

	dx = (b-a)/n;
	s = 0.5*(func(a)+func(b));
	for(i=1;i<n;i++)
		s += func(a+i*dx);

	s *= dx;
	return s;
}
dcomplex dcqtrap(func,a,b,n)
float a,b;
dcomplex (*func)();
{	
	int i;
	dcomplex s,dum;
	double dx;

	dx = (b-a)/n;
	dum = func(a);
	s.re = dum.re;
	s.im = dum.im;
	dum = func(b);
	s.re += dum.re;
	s.im += dum.im;
	s.re *=0.5;
	s.im *=0.5;
	for(i=1;i<n;i++)
	{
		dum = func(a+i*dx);
		s.re += dum.re;
		s.im += dum.im;
	}

	s.re *= dx;
	s.im *= dx;
	return s;
}
#undef EPS
#undef JMAX
#define  FUNC(x)	((*func)(x))

double trapzd(func,a,b,n)
double a,b;
double (*func)();	/* ANSI: float (*func)(float); */
int n;
{
	double x,tnm,sum,del;
	static double s;
	static int it;
	int j;

	if (n == 1) {
		it=1;
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=0;j<it;j++,x+=del) sum += FUNC(x);
		it *= 2;
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}

dcomplex dctrapzd(func,a,b,n)
double a,b;
dcomplex (*func)();	/* ANSI: dcomplex (*func)(float); */
int n;
{
	double x,tnm,del;
	dcomplex sum;
	static dcomplex s;
	static int it;
	dcomplex dctmp1,dctmp2;
	double tmp;
	int j;

	if (n == 1) {
		it=1;
/*		return (s=0.5*(b-a)*((*func)(a)+(*func)(b)); */
		tmp = 0.5*(b-a);
		dctmp1 = (*func)(a);
		dctmp2 =(*func)(b);
		s.re = tmp*(dctmp1.re+dctmp2.re);
		s.im = tmp*(dctmp1.im+dctmp2.im);
		return s;
	} else {
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
/*		for (sum=0.0,j=0;j<it;j++,x+=del) sum += FUNC(x); */
		for (sum.re=sum.im=0.0,j=0;j<it;j++,x+=del) 
		{
			dctmp1 = (*func)(x);
			sum.re += dctmp1.re;
			sum.im += dctmp1.im;
		}
		it *= 2;
/*		s=0.5*(s+(b-a)*sum/tnm); */
		tmp = (b-a)/tnm;
		s.re = 0.5*(s.re+tmp*sum.re);
		s.im = 0.5*(s.im+tmp*sum.im);
		return s;
	}
}


void quadra_interp(f0,df0,f1,df1,h,x,y,dy,ddy)
/*
interpolation of a function from the functinal vaules and first derivatives
at the end points.

f0:  value at x=0
df0: first derivative at x=0
f1:  value at x=h
df1: fist derivative at x=h
h:   distance between end points
x:   estimate at x
y:   interpolated value 
dy:  interpolated value of first derivative
ddy: interpolated value of second derivative
return interpolated value at x
*/
	float f0,df0,f1,df1,h,x;
	float *y,*dy,*ddy;
{
	double a,b;

	a=((df0+df1)+2.0*(f0-f1)/h)/h/h;
	b=(-(df1+2*df0)+3.0*(f1-f0)/h)/h;

	*y  = ((a*x+b)*x+df0)*x+f0;
	*dy = (3.0*a*x+2.0*b)*x+df0;
	*ddy = 6.0*a*x+2.0*b;
}

double linear_interp(f0,f1,h,x)
/*
linear interpolation of two end points.
f0: value at x=0
f1: value at x=h
h:  distance between two end points
x:  estimate at x

return interpolated value at x
*/
	float f0,f1,h,x;
{
	return ((double)f1-f0)/h*x+f0;
}
