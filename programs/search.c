#include <stdio.h>
#include "search.h"
#include "math.h"
#include "m_solid.h"
#include "calc_a.earth.h"

#define FR		fprintf(stderr

void cross_2d(type,x_guess,y_guess,f1,f2,x,y)
/* determine zero cross pattern of 2 functions
input
double f1[], values of function f1 at four corners at (0-1-2-3)
double f2[], values of function f2 at four corners at (0-1-2-3)
output
int type
*/
int *type;
double *x_guess,*y_guess,f1[4],f2[4],x[4],y[4];
{
	int i,j,k;
	int sf1[4],sf2[4];
	char f1_flag,f2_flag;
	double xc[2],yc[2],a[2],b[2],c[2],d;
	
	f1_flag=f2_flag=0x0;
	for(i=0;i<4;i++)
	{
		sf1[i]=(f1[i] > 0.0 ? 1 : -1 );
		sf2[i]=(f2[i] > 0.0 ? 1 : -1 );
		f1_flag|=(sf1[i] == 1 ? 0x01 : 0x00 )<<i;
		f2_flag|=(sf2[i] == 1 ? 0x01 : 0x00 )<<i;
	}

/* no crossing */
	if (f1_flag ==0x00 || f1_flag ==0x0f ||
	    f2_flag ==0x00 || f2_flag ==0x0f )
	{ *type = NO_CROSSING; return; }
/* find crossing coordinate (x,y) points */
	j=0;
	/* f1 crosing */
	for(i=0;i<4;i++)
	{
		k = (i+1)%4;
		if ( f1[i] == 0.0 || f1[k] == 0.0)
		{ 
			FR,"cross_2d: change grid points f1[%d] = 0.0.\n",i);
			exit(-1);
		}
		else if( sf1[i]*sf1[k] == -1)
		{
			if( j++ > 1){
				FR,"cross_2d(%d) %d: Warning: zero-crossing %d times\n",
				__FILE__,__LINE__,j);
				FR,"regard as no crossing.\n");
				*type = NO_CROSSING; return;
			}
			cross_1d(&(xc[j-1]),x[i],x[k],f1[i],f1[k]);
			cross_1d(&(yc[j-1]),y[i],y[k],f1[i],f1[k]);
		} 
	}
	if( j != 2 )
	{FR,"cross_2d(%s) %d: f1 zero-crossing only %d times\n",
	__FILE__,__LINE__,j);exit(-1);}
	d = xc[0]*yc[1]-xc[1]*yc[0];
	if( d == 0.0 ){
		a[0] = -yc[0];
		b[0] =  xc[0];
		c[0] =  0.0;
	}else{
		a[0] = (yc[1]-yc[0])/d;
		b[0] = (xc[0]-xc[1])/d;
		c[0] = 1.0;
	}

	j=0;
	/* f2 crossing */
	for(i=0;i<4;i++)
	{
		k = (i+1)%4;
		if ( f2[i] == 0 || f2[k] == 0.0 )
		{ 
			FR,"cross_2d: change grid points f1[%d] = 0.0.\n",i);
			exit(-1);
		}
		else if( sf2[i]*sf2[k] == -1)
		{
			cross_1d(&(xc[j]),x[i],x[k],f2[i],f2[k]);
			cross_1d(&(yc[j]),y[i],y[k],f2[i],f2[k]);
			j ++;
		} 
	}
	if( j != 2 )
	{FR,"cross_2d(%s) %d: f2 zero-crossing only %d times\n",
	__FILE__,__LINE__,j);exit(-1);}
	d = xc[0]*yc[1]-xc[1]*yc[0];
	if( d == 0.0 ){
		a[1] = -yc[0];
		b[1] =  xc[0];
		c[1] =  0.0;
	}else{
		a[1] = (yc[1]-yc[0])/d;
		b[1] = (xc[0]-xc[1])/d;
		c[1] = 1.0;
	}

	d = a[0]*b[1]-a[1]*b[0];
	*x_guess = ( b[1]*c[0] - b[0]*c[1])/d;
	*y_guess = (-a[1]*c[0] + a[0]*c[1])/d;

	*type = CROSSING;
}

void cross_1d(xc,x1,x2,f1,f2)
/* find crossing point of xc of a linear function which goes through
(x1,f1) and (x2,f2) f1*f2 < 0 
*/
double *xc,x1,x2,f1,f2;
{
	if( f1*f2 >= 0.0 ){FR,"cross_1d input error f1=%e,f1=%e\n",f1,f2);
	exit(-1);}
	if( x1 == x2 ) {*xc = x1; return;}
	if( f1 > 0.0 )
		*xc = x1 + (x2-x1)*f1/(f1-f2);
	else
		*xc = x1 + (x2-x1)*(-f1)/(-f1+f2);
}
#define OUT(x,a,b) ((x) < (a) || (b) < (x))

void inside_outside(type,x1,y1,x2,y2,xp,yp,x_margine,y_margine)
/*
	determine (xp,yp) is inside of box (x1,y1)-(x2,y2) or outside
	or withtin margine
	x1<x2, y1<y2
	dx = x2-x1
	dy = y2-y1
	inside of margine ==
	box (x1-dx*x_margine,y1-dy*y_margine)-(x2+dx*x_margine,y2+dy*y_margine)
*/
int *type;
double x1,y1,x2,y2,xp,yp,x_margine,y_margine;
{
	double dx,dy;
	if( x1 >=x2 || y1 >=y2 ){FR,"inside_outside: input error.\n");
	FR,"(x1,y1)=(%lf,%lf)\n",x1,y1);
	FR,"(x2,y2)=(%lf,%lf)\n",x2,y2);
	}
	if( !OUT(xp,x1,x2) && !OUT(yp,y1,y2) )
		{*type = INSIDE_XY; return;}
	else 
	{
		dx = (x2-x1)*x_margine;
		dy = (y2-y1)*y_margine;
		if( !OUT(xp,x1-dx,x2+dx) && !OUT(yp,y1-dy,y2+dy) )
		{*type = INSIDE_MARGIN; return;}
	}
	*type = OUTSIDE_MARGIN;
}

#undef OUT
void improve_dx(dx1,dx2,f1,f2,delta_x1)
double *dx1,*dx2;
double f1[2],f2[2],delta_x1;
/* compute the improvement (dx1,dx2)
of a trial solution (x1,x2) of a analytic function
f(x1+i*x2)=f1(x1+i*x2) + i*f2(x1+i*x2)
for given values
f1[0] = f1(x1+delta_x1,x2)
f1[1] = f1(x1-delta_x1,x2)
f2[0] = f2(x1+delta_x1,x2)
f2[1] = f2(x1-delta_x1,x2)
*/
{
	double d1,d2,s1,s2,ratio,factor;
	d1 = f1[0]-f1[1];
	d2 = f2[0]-f2[1];
	s1 = f1[0]+f1[1];
	s2 = f2[0]+f2[1];
/*
	ratio = d2/d1;
	factor = delta_x1/(d1+ratio*d2);
	*dx1 = -factor*(s1+ratio*s2);
	*dx2 = -factor*(s2-ratio*s1);
*/
	if( d1 == 0.0 && d2 == 0.0 )
	{FR,"improve input error d1=d2=delta_x1=0.0.\n"); exit(-1);}
	factor = -delta_x1/(d1*d1+d2*d2);
	*dx1 = (d1*s1+d2*s2)*factor;
	*dx2 = (-d2*s1+d1*s2)*factor;
}

#define SMALL     1e-10 /* tol is about 1e-15 is balanced value */
#define MAXSTEPS  100
#define     ABS(a)      ((a) > 0.0 ? (a) : -(a))
void search_omg(best_omg,omg,omg1,gamma,gamma1,r1,l,k,m,tol)
double *best_omg,omg,omg1,gamma,gamma1;
double r1,k,tol;
int l;
struct model_funcsolid *m;
/* best_omg, omg, omg1 are square of omg */
{
	double a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm,min1,min2;
	int i;

	a = omg;
	b = omg1;
	fa = gamma;
	fb = gamma1;

	if(fa*fb > 0.0 )
	{
		FR,"search_omg(%s) %d: Root must be brecketed.\n",__FILE__,__LINE__);
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

		calc_fluid2by2_character(&fb,r1,l,b,k,*m);
		
	}
	FR,"max step reached in search_omg.\n");
}
#define	MAX_LOOP	1024
void find_omg_bound(find_flag,omg_acou,omg_grav,omg_start,omg_end,r,l,m,tol)
/* find for a given fluid earth mode and radius and modal angular orger,
search the frequency boundary of acoustic mode and graviry mode
INPUT
	omg_start,	starting value of angular frequency
	omg_end,	ending value of angular frequency
	r,		radius from center normalized as with m 
	m,		earth model structure at r
	tol,		tolerance level of searching
	## NOTE ##
	omg_start, omg_end,r and m must be normalized consistently.
	## NOTE END ##
OUTPUT
	find_flag	FOUND_OMG_[NOTHING_EVA|NOTHING_PRO|ACOUSTIC|GRAVITY|BOTH]
	omg_acou,	angular frequency boundary between acoustic mode
			region and evanescent region
	omg_grav,	angular frequency boundary between evanescent mode
			and gravity mode region
*/
int *find_flag;
double *omg_acou,*omg_grav,omg_start,omg_end,r,tol;
int l;
struct model_funcsolid *m;
{
	int ndiv,i;
	double omg2,old_omg2,gamma,old_gamma,gamma_start,gamma_end;
	double omg2_start,omg2_end,best_omg2,dum;
	double dlogomg, logomg;
	double omg_buoy2,omg_lamb2,omg2_mean;

	if( omg_start > omg_end)
	{dum=omg_start;omg_start=omg_end;omg_end=dum;}
	if(omg_start <=0.0 )
	{FR,"find_omg_bound(%s) input error: omg_start=%f omg_end=%f\n",
	__FILE__,omg_start,omg_end);exit(-1);}

	omg2_start = omg_start * omg_start;
	calc_fluid2by2_character(&gamma_start,r,l,omg2_start,-1.0,*m);
	omg2_end = omg_end * omg_end;
	calc_fluid2by2_character(&gamma_end,r,l,omg2_end,-1.0,*m);

	*find_flag = -999;
	if(gamma_start > 0.0  && gamma_end > 0.0 )
	{
	/* no solution */
		*omg_acou = *omg_grav =0.0;
		*find_flag = FOUND_OMG_NOTHING_EVA;
	}else if(gamma_start < 0.0 && gamma_end  > 0.0
	 ||      gamma_start > 0.0 && gamma_end < 0.0 )
	{
		search_omg(&best_omg2,omg2_start,omg2_end,gamma_start,gamma_end,r,l,-1.0,m,tol);
	/* gravity region boundary */
		if( gamma_start < 0.0 )
		{
			*find_flag = FOUND_OMG_GRAVITY;
			*omg_grav = sqrt(best_omg2);
		}
	/* acoustic region boundary */
		else{
			*find_flag = FOUND_OMG_ACOUSTIC;
			*omg_acou = sqrt(best_omg2);
		}

	}else{
	/* both gamma's are negative */
		ndiv = 1; /* this division algorithm is similar to the
		one in trapzd() in Numerical Recipe */
		while(ndiv < MAX_LOOP )
		{
		dlogomg = (log(omg_end) - log(omg_start))/ndiv;
		logomg = log(omg_start) + 0.5*dlogomg;
		for(i=0;i<ndiv;i++)
		{
			dum = exp(logomg);
			omg2 = dum*dum;
			calc_fluid2by2_character(&gamma, r,l,omg2,-1.0,*m);
			if( gamma > 0.0 )
			{
				dum = exp(logomg - 0.5*dlogomg);
				old_omg2 = dum*dum;
				calc_fluid2by2_character(&old_gamma,r,l,old_omg2,-1.0,*m);
				search_omg(&best_omg2,omg2,old_omg2,gamma,old_gamma,r,l,-1.0,m,tol);
				*omg_grav = sqrt(best_omg2);

				old_gamma = gamma;
				old_omg2 = omg2;
				dum = exp(logomg + 0.5*dlogomg);
				omg2 = dum*dum;
				calc_fluid2by2_character(&gamma,r,l,omg2,-1.0,*m);
				search_omg(&best_omg2,omg2,old_omg2,gamma,old_gamma,r,l,-1.0,m,tol);
				*omg_acou = sqrt(best_omg2);
				*find_flag = FOUND_OMG_BOTH;
				break;
			}
			logomg += dlogomg;
		}
		if( *find_flag == FOUND_OMG_BOTH ) break;
		ndiv *=2;
		}
		if( ndiv == MAX_LOOP)
		{
			calc_fluid_omg(&omg_buoy2, &omg_lamb2,r,l,*m);
			omg2_mean = sqrt(omg_buoy2*omg_lamb2);
			if (omg2_start > omg2_mean && omg2_end > omg2_mean )
				*find_flag = FOUND_OMG_NOTHING_ACO;
			else if (omg2_start < omg2_mean && omg2_end < omg2_mean )
				*find_flag = FOUND_OMG_NOTHING_GRA;
			else
			{
				*find_flag = FOUND_OMG_NOTHING_PRO;
				FR,"find_omg_boound(%s) %d: Warning.\n");
				FR,"omg_start = %10.4e omg_end = %10.4e is either in acoustic or gravity region. cannot be determined.\n",omg_start,omg_end);
			}
		}
	}
	if( *find_flag == -999)
	{FR,"find_omg_bound(%s) %d: coding error.\n");exit(-1);}
}

#define OUT(x,a,b)	(x<a || x>b)
void determine_omg_domain(domain,omg,flag,omg_acou,omg_grav)
/* 
   from flag,omg_acou,omg_grav of find_omg_bound() determine
   domain of omg
INPUT
	omg,		angular frequency to be examined
	flag,		FOUND_OMG_BOTH,FOUND_OMG_ACOUSTIC,FOUND_OMG_GRAVITY,
			or FOUND_NOTHING
	omg_acou,	angular frequency boundary between acoustic mode
			region and evanescent region
	omg_grav,	angular frequency boundary between evanescent mode
			and gravity mode region
OUTPUT
	domain, 	EVANESCENT_DOMAIN,ACOUSTIC_DOMAIN,GRAVITY_DOMAIN
			or PROPAGATING_DOMAIN
*/
int *domain,flag;
double omg,omg_acou,omg_grav;
{
	switch( flag ){
	case FOUND_OMG_BOTH:
		if( omg_acou <= omg_grav || omg_acou <=0.0 || omg_grav <= 0.0)
			{FR,"input error.\n");goto err;}
		if(!OUT(omg,omg_acou,omg_grav)) *domain = EVANESCENT_DOMAIN;
		else if( omg > omg_acou) *domain = ACOUSTIC_DOMAIN;
		else if( omg < omg_grav) *domain = GRAVITY_DOMAIN;
		else goto err;
		break;
	case FOUND_OMG_ACOUSTIC:
		if( omg_acou <=0.0 )
			{FR,"input error.\n");goto err;}
		if( omg > omg_acou ) *domain = ACOUSTIC_DOMAIN;
		else *domain = EVANESCENT_DOMAIN;
		break;
	case FOUND_OMG_GRAVITY:
		if( omg_grav <= 0.0 )
			{FR,"input error.\n");goto err;}
		if( omg < omg_grav ) *domain = GRAVITY_DOMAIN;
		else *domain = EVANESCENT_DOMAIN;
		break;
	case FOUND_OMG_NOTHING_EVA:
		*domain = EVANESCENT_DOMAIN;
		break;
	case FOUND_OMG_NOTHING_PRO:
		*domain = PROPAGATING_DOMAIN;
		break;
	case FOUND_OMG_NOTHING_GRA:
		*domain = GRAVITY_DOMAIN;
		break;
	case FOUND_OMG_NOTHING_ACO:
		*domain = ACOUSTIC_DOMAIN;
		break;
	default:
		FR,"determine_omg_domain(%s) %d: unknown domain type=%d.\n",
		__FILE__,__LINE__,flag);exit(-1);
	}
	return;

	err:
	FR,"determine_omg_domain(%s) %d: flag=%d omg=%10.4e,omg_acou=%10.4e,omg_grav=%10.4e\n");
	exit(-1);
}
#undef OUT
