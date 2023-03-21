#include "bess.h"
#define FR		fprintf(stderr
#define MAX(a,b)	((a)>(b) ? (a):(b))

double sbess_ratio(rn, xx)
/* compute
	j   (x) / (x*j   (x))
	 rn+1         rn
   where
   x = sqrt(xx)
   j  (x) is spherical Bessel function of first kind

copied from  DISPER80 Package DSPZN()
*/
	double rn; /* order can be non-integer */
	double xx; /* xx may be less than 0 */
{
	double s,fm,z1,z2;
	int i,l,m,n;
	if (rn < 0.0)
	{FR,"sbess_ratio(%s): input error n=%f\n",__FILE__,rn);
	exit(-1);}

	n=(int)rn;
	s=sqrt(fabs(xx));
	if( s > 100.)        l = 30 + (int)(s*0.05);
	else if ( s > 10.)   l = 20 + (int)(s*0.15);
	else if ( s > 1. )   l =  9 + (int)s;
	else                 l = 10;

	m = (int)s;
	l = l + MAX(n,m) - n;
	fm = 2.0*(rn+l)+1;
	z1 = 0;
/*
 * continued fraction method
 */
	for (i=0;i<l;i++)
	{
		z2 = z1;
		z1 =  1/(fm - xx*z2);
		fm = fm - 2.0;
	}
	return z1;
}

double sbess_ratio1(rn, xx)
/* compute
	j   (x) / (x*j   (x))
	 rn+1         rn
   where
   x = sqrt(x)
   j  (x) is spherical Bessel function of first kind

copied from  DISPER80 Package DSPZN()
*/
	double rn; /* order can be non-integer */
	double xx; /* xx may be less than 0 */
{
	double s,fm,z1,z2;
	int i,l,m,n;
	if (rn < 0.0)
	{FR,"sbess_ratio(%s): input error n=%f\n",__FILE__,rn);
	exit(-1);}

	n=(int)rn;
	s=sqrt(fabs(xx));
	if( s > 100.)        l = 60 + (int)(s*0.05);
	else if ( s > 10.)   l = 40 + (int)(s*0.15);
	else if ( s > 1. )   l = 18 + (int)s;
	else                 l = 20;

	m = (int)s;
	l = l + MAX(n,m) - n;
	fm = 2.0*(rn+l)+1;
	z1 = 0;
/*
 * continued fraction method
 */
	for (i=0;i<l;i++)
	{
		z2 = z1;
		z1 =  1/(fm - xx*z2);
		fm = fm - 2.0;
	}
	return z1;
}
