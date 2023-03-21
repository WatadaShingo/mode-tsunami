#include <stdio.h>
#include <math.h>
#include "clib.h"

#define FR	fprintf(stderr
complex comp(a)
	dcomplex a;
{
	complex b;
	b.re = (float)a.re;
	b.im = (float)a.im;
	return b;
}
dcomplex dcomp(a)
	complex a;
{
	dcomplex b;
	b.re = (double)a.re;
	b.im = (double)a.im;
	return b;
}
complex ccnj(a)
	complex a;
{
	a.im = -a.im;
	return a;
}
dcomplex dccnj(a)
	dcomplex a;
{
	a.im = -a.im;
	return a;
}
complex cmk(a,b)
	float a,b;
{
	complex c;
	c.re = a; c.im = b;
	return c;
}
dcomplex dcmk(a,b)
	double a,b;
{
	dcomplex c;
	c.re = a; c.im = b;
	return c;
}
complex cpow(a,b)
	complex a,b;
{
	return ( cexp(cmul(b,clog(a))) );

}
dcomplex dcpow(a,b)
	dcomplex a,b;
{
	return ( dcexp(dcmul(b,dclog(a))) );

}
complex clog(a)
	complex a;
{
	complex c;
	polar b;
	b = cpol(a);
	if( b.r <= 0.0 )
	{
	   FR,"error in complex func log(a) a.re<=0.\n");
	   exit(-1);
	} 
	c.re = log(b.r);
	c.im = b.t;
	return c;
}
dcomplex dclog(a)
	dcomplex a;
{
	dcomplex c;
	dpolar b;
	b = dcpol(a);
	if( b.r <= 0.0 )
	{
	   FR,"error in complex func log(a) a.re<=0.\n");
	   exit(-1);
	} 
	c.re = log(b.r);
	c.im = b.t;
	return c;
}
complex cexp(a)
	complex a;
{
	complex b;
	double c,s;
	sincos(a.im,&s,&c);
	b.re = exp(a.re)*c;
	b.im = exp(a.re)*s;
/* dbg 
	FR," cexp %f %f\n",b.re,b.im);
 dbg end */
	return b;
}
dcomplex dcexp(a)
	dcomplex a;
{
	dcomplex b;
	double c,s;
	sincos(a.im,&s,&c);
	b.re = exp(a.re)*c;
	b.im = exp(a.re)*s;
/* dbg 
	FR," cexp %f %f\n",b.re,b.im);
 dbg end */
	return b;
}
complex csum(a,b)
	complex a,b;
{
	a.re +=b.re; a.im+=b.im;
	return a;
}
dcomplex dcsum(a,b)
	dcomplex a,b;
{
	a.re +=b.re; a.im+=b.im;
	return a;
}
complex csub(a,b)
	complex a,b;
{
	a.re -=b.re; a.im-=b.im;
/* dbg
	FR," csub %f %f\n",a.re,a.im);
dbg end */
	return a;
}
dcomplex dcsub(a,b)
	dcomplex a,b;
{
	a.re -=b.re; a.im-=b.im;
/* dbg
	FR," csub %f %f\n",a.re,a.im);
dbg end */
	return a;
}
complex cscl(s,a)
	float s;
	complex a;
{
	complex c;

	c.re = s*a.re;
	c.im = s*a.im;
	return c;
}
dcomplex dcscl(s,a)
	double s;
	dcomplex a;
{
	dcomplex c;

	c.re = s*a.re;
	c.im = s*a.im;
	return c;
}
complex cmul(a,b)
	complex a,b;
{
	complex c;
	c.re = a.re*b.re - a.im*b.im;
	c.im = a.re*b.im + a.im*b.re;
/* dbg 
	FR," cmul %f %f\n",c.re,c.im);
 dbg end */
	return c;
}
dcomplex dcmul(a,b)
	dcomplex a,b;
{
	dcomplex c;
	c.re = a.re*b.re - a.im*b.im;
	c.im = a.re*b.im + a.im*b.re;
/* dbg 
	FR," cmul %f %f\n",c.re,c.im);
 dbg end */
	return c;
}
complex cdiv(a,b)
	complex a,b;
{
	complex c;
	double norm;
	if((norm=b.re*b.re+b.im*b.im) == 0.0 )
	{
		FR,"\tdivided by zero in cdiv\n");
		exit(-1);
	}
	else if( norm > 1.0e100)
	{
		FR,"\tdivided by infinit in cdiv\n");
		exit(-1);
	}
	c.re = (a.re*b.re + a.im*b.im)/norm;
	c.im = (a.im*b.re - a.re*b.im)/norm;
/* dbg 
	FR," cdiv %f %f\n",c.re,c.im);
  dbg end */
	return c;
}
dcomplex dcdiv(a,b)
	dcomplex a,b;
{
	dcomplex c;
	double norm;
	if((norm=b.re*b.re+b.im*b.im) == 0.0 )
	{
		FR,"\tdivided by zero in cdiv\n");
		exit(-1);
	}
	else if( norm > 1.0e400)
	{
		FR,"\tdivided by infinit in dcdiv\n");
		exit(-1);
	}
	c.re = (a.re*b.re + a.im*b.im)/norm;
	c.im = (a.im*b.re - a.re*b.im)/norm;
/* dbg 
	FR," cdiv %f %f\n",c.re,c.im);
  dbg end */
	return c;
}
polar cpol(a)
	complex a;
{
	polar b;
	b.r = sqrt(a.re*a.re + a.im*a.im );
	b.t = (a.re == 0.0 && a.im == 0.0 ? 0.0 : atan2( a.im, a.re )); 
	return b;
}
dpolar dcpol(a)
	dcomplex a;
{
	dpolar b;
	b.r = sqrt(a.re*a.re + a.im*a.im );
	b.t = (a.re == 0.0 && a.im == 0.0 ? 0.0 : atan2( a.im, a.re )); 
	return b;
}
