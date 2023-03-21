#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rk.h"

#define MAX_FUNCTION_NUMBER	10
#define FR		fprintf(stderr

#define 	MAX(a,b)	((a)>(b) ? (a) : (b))
#define 	MIN(a,b)	((a)<(b) ? (a) : (b))
#define	ABS(a)	((a) > 0.0 ? (a) : -(a))

int rungev4(y,dydt,yout,t,dt,fdydt,n,l,omg)
/*
 * input y[],       initial value y at t		not modified upon exit
	   dydy[],	initial value dy/dy at t	not modofied upon exit
	   t,       parameter variable
         dt,      incremental interval
         fdydt    function to give dy/dt(t,y) 
	   n,       dimension of input y, output yout
	   l,		angular order
	   omg,	angular frequency
   output yout[],   estimated value at t+dt
*/
	double y[],dydt[],yout[],t,dt;
	void (*fdydt)();
	int n,l;
	double omg;
{
	int i;	/* i-th equation */
	double k1[MAX_FUNCTION_NUMBER],
	       k2[MAX_FUNCTION_NUMBER],
		 k3[MAX_FUNCTION_NUMBER],
		 k4[MAX_FUNCTION_NUMBER],t1;
	double ytemp[MAX_FUNCTION_NUMBER],
		ddum[MAX_FUNCTION_NUMBER];

/* fourth order Runge Kutta. See Numrical recipes Chapter 15 */

	if( MAX_FUNCTION_NUMBER < n )
	{
		fprintf(stderr,"runge4v: too many array dimention\n");
		exit(-1);
	}
	t1 = 0.5*dt + t;

	for(i=0;i<n;i++) k1[i]    = dt*dydt[i];
	for(i=0;i<n;i++) ytemp[i] = y[i] + 0.5*k1[i];

	(* fdydt)(ddum,ytemp,t1,l,omg);
	for(i=0;i<n;i++) k2[i]    = dt*ddum[i];
	for(i=0;i<n;i++) ytemp[i] = y[i] + 0.5*k2[i];

	(* fdydt)(ddum,ytemp,t1,l,omg);
	for(i=0;i<n;i++) k3[i]    = dt*ddum[i];
	for(i=0;i<n;i++) ytemp[i] = y[i] + k3[i];
	(* fdydt)(ddum,ytemp,t+dt,l,omg);
	for(i=0;i<n;i++) k4[i]    = dt*ddum[i];

	for(i=0;i<n;i++)
		yout[i] = y[i] + (k1[i]+k4[i])/6.0+(k2[i]+k3[i])/3.0;
}

int rungevck(y,dydt,yout,yerr,t,dt,fdydt,n,l,omg)
/*
	Cash-Karp type embedded runge-kutta step.
	equivalent to the rkrc() routine in Numerical Recipes 2nd ed.

 * input y[],       initial value y at t		not modified upon exit
	   dydy[],	initial value dy/dy at t	not modofied upon exit
	   t,       parameter variable
         dt,      incremental interval
         fdydt    function to give dy/dt(t,y) 
	   n,       dimension of input y, output yout
	   l,		angular order
	   omg,	angular frequency
   output yout[],   estimated value at t+dt
          yerr[],   estimate of error of yout[]
*/
	double y[],dydt[],yout[],yerr[],t,dt;
	void (*fdydt)();
	int n,l;
	double omg;
{
	int i;	/* i-th equation */
	double k1[MAX_FUNCTION_NUMBER],
	       k2[MAX_FUNCTION_NUMBER],
		 k3[MAX_FUNCTION_NUMBER],
		 k4[MAX_FUNCTION_NUMBER],
		 k5[MAX_FUNCTION_NUMBER],
		 k6[MAX_FUNCTION_NUMBER];
	double ytemp[MAX_FUNCTION_NUMBER],
		ddum[MAX_FUNCTION_NUMBER];
	static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
		b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
		b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
		b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
		b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
		c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
		dc5 = -277.0/14336.0;
	double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
		dc4=c4-13525.0/55296.0,dc6=c6-0.25;

/* fourth order Runge Kutta. See Numrical recipes Chapter 15 */

	if( MAX_FUNCTION_NUMBER < n )
	{
		fprintf(stderr,"runge4v: too many array dimention\n");
		exit(-1);
	}

	for(i=0;i<n;i++) k1[i]    = dt*dydt[i];
	for(i=0;i<n;i++) ytemp[i] = y[i] + b21*k1[i];

	(* fdydt)(ddum,ytemp,t+a2*dt,l,omg);
	for(i=0;i<n;i++) k2[i]    = dt*ddum[i];
	for(i=0;i<n;i++) ytemp[i] = y[i] + b31*k1[i] + b32*k2[i];

	(* fdydt)(ddum,ytemp,t+a3*dt,l,omg);
	for(i=0;i<n;i++) k3[i]    = dt*ddum[i];
	for(i=0;i<n;i++) ytemp[i] = y[i] + b41*k1[i] + b42*k2[i] + b43*k3[i];

	(* fdydt)(ddum,ytemp,t+a4*dt,l,omg);
	for(i=0;i<n;i++) k4[i]    = dt*ddum[i];
	for(i=0;i<n;i++) ytemp[i] = y[i] + b51*k1[i] + b52*k2[i] + b53*k3[i]
					  + b54*k4[i];

	(* fdydt)(ddum,ytemp,t+a5*dt,l,omg);
	for(i=0;i<n;i++) k5[i]    = dt*ddum[i];
	for(i=0;i<n;i++) ytemp[i] = y[i] + b61*k1[i] + b62*k2[i] + b63*k3[i]
					  + b64*k4[i] + b65*k5[i];

	(* fdydt)(ddum,ytemp,t+a6*dt,l,omg);
	for(i=0;i<n;i++) k6[i]    = dt*ddum[i];

	for(i=0;i<n;i++)
		yout[i] = y[i] + c1*k1[i] + c3*k3[i] + c4*k4[i] + c6*k6[i];
	for(i=0;i<n;i++)
		yerr[i] = (dc1*k1[i] + dc3*k3[i] + dc4*k4[i] + dc5*k5[i] 
			  + dc6*k6[i])*dt;
}

#define SAFTY	0.9
#define PGROW	0.2
#define PSHRNK	0.25
#define FCOR	15.0
#define ERRCON    6e-4
void stepper(y,yout,dydt,t,dt,fdydt,n,yscale,eps,dtdid,dtnext,l,omg)
/*
 * input y[],       inital value y at t		not modified	upon exit
	   dydy[],	initial value dy/dy at t	not modified	upon exit
	   t,       parameter variable
         dt,      incremental interval
         fdydt()    function to give dy/dt(t,y) 
	   n,       dimension of input y, output yout
	   yscale[],  present scale of y		not modified upon exit
	   eps,     required accuracy
	   l,		angular order
	   omg,	angular frequency
   output yout[],   estimated value at t+dt
          dtdid   increased y--(dt)-->yout
	    dtnext  dt for next calculation
*/
	double y[],yout[],dydt[],yscale[],t,dt,*dtdid,*dtnext,eps;
	void (*fdydt)();
	int n,l;
	double omg;
{
	int i;
	double dthalf,maxerr,tsav;
	double ytemp[ MAX_FUNCTION_NUMBER ];
	double dydttemp[ MAX_FUNCTION_NUMBER ];

	tsav = t;
	while(1)
	{
		maxerr = 0.0;
		dthalf = dt/2.0;
	      rungev4(y,dydt,ytemp,tsav, dthalf,fdydt,n,l,omg);
		t = tsav+dthalf;
		(*fdydt)(dydttemp,ytemp,t,l,omg);
     		rungev4(ytemp,dydttemp,  yout, t,    dthalf,fdydt,n,l,omg);
		t = tsav + dt;
		if( t == tsav )
		{
			fprintf(stderr,"stepper: too small time interval\n");
			exit(-1);
		}
     		rungev4(y,dydt,ytemp,tsav,    dt,    fdydt,n,l,omg);
		for(i=0;i<n;i++)
		{
			ytemp[i]=yout[i]-ytemp[i];
			maxerr = MAX(maxerr,fabs(ytemp[i]/yscale[i]));
		}
		maxerr = eps/maxerr;
		if ( maxerr >= 1.0 ) break;
		dt = SAFTY*dt*pow(maxerr,PSHRNK);
	}
	*dtdid = dt;
	if( maxerr < 1.0/ERRCON )
		*dtnext = SAFTY*dt*pow(maxerr,PGROW);
	else
		*dtnext = 4.0 * dt;

	for(i=0;i<n;i++)
		yout[i] = yout[i]+ytemp[i]/FCOR;
}
#undef SAFTY
#undef PGROW
#undef PSHRNK
#undef FCOR
#undef ERRCON

#define SAFETY	0.9
#define PGROW	-0.2
#define PSHRNK	-0.25
#define ERRCON    1.89e-4

void stepper_ck(y,yout,dydt,t,dt,fdydt,n,yscale,eps,dtdid,dtnext,l,omg)
/*
	modified version of stepper().
	Using Cash-Karp type embedded runge-kutta, a variation of fehlberg type.
	this routine call rungevck() insted rungev4();
	equivalent to the rkqs() routine in Numerical Recipes 2nd ed.

 * input y[],       inital value y at t		not modified	upon exit
	   dydy[],	initial value dy/dy at t	not modified	upon exit
	   t,       parameter variable
         dt,      incremental interval
         fdydt()    function to give dy/dt(t,y) 
	   n,       dimension of input y, output yout
	   yscale[],  present scale of y		not modified upon exit
	   eps,     required accuracy
	   l,		angular order
	   omg,	angular frequency
   output yout[],   estimated value at t+dt
          dtdid   increased y--(dt)-->yout
	    dtnext  dt for next calculation
*/
	double y[],yout[],dydt[],yscale[],t,dt,*dtdid,*dtnext,eps;
	void (*fdydt)();
	int n,l;
	double omg;
{
	int i;
	double maxerr,tnew,dttemp;
	double yerr[ MAX_FUNCTION_NUMBER ];

	while(1)
	{
		maxerr = 0.0;
		rungevck(y,dydt,yout,yerr,t,dt,fdydt,n,l,omg);
		for(i=0;i<n;i++)
			maxerr = MAX(maxerr,fabs(yerr[i]/yscale[i]));
		maxerr /= eps;
		if ( maxerr > 1.0 ){
			dttemp = SAFETY*dt*pow(maxerr,PSHRNK);
			dt=(dt >=0.0 ? MAX(dttemp,0.1*dt) : MIN(dttemp,0.1*dt));
			tnew = t+dt;
			if( tnew == t )
			{
				fprintf(stderr,"stepper_ck: too small time interval\n");
				exit(-1);
			}
			continue;
		}else{
			if( maxerr > ERRCON ) *dtnext = SAFETY*dt*pow(maxerr,PGROW);
			else *dtnext = 5.0*dt;
			*dtdid = dt;
			break;
		}
	}
}
#define NO			0
#define YES			1
void runge_driver_match(y_over,y_grid,dydt_grid,y0_botm,y0_top,t_grid,
	layerup_botm,layerup_top,layerdown_botm,layerdown_top,
	fdydt,n_order,eps,l,omg,is_funcs)
/* input
	y0_top[n],	initial value of y at t_grid[layerdown_top]
	y0_botm[n],	initial value of y at t_grid[layerup_botm]
	t_grid[i],	grid values of t axis where y and dydt are evaluated
	layerup_botm,	integrate upward from here 
	layerup_top,	integrate upward to here
	layerdown_botm,	integrate downward to here
	layerdown_top,	integrate downward from here here
		t_grid[layerup_top] == t_grid[layerdown_botm]
	fdydt,	function to give dy/dt(t,y)
	n_odrder
	eps,		required accuracy
	l,          angular order
	omg,        angular frequency
	is_funcs,   store all eigenfunctions == YES if not == NO
  output
	y_over[n],	integrated value at grid=layerdown_botm
	y_grid[n][i]	estimated value of y at t_grid[i],
		layerup_botm<=i<=layerdown_top.
		if i=layerup_top==layerdown_botm
		the value at layer=i is from the integration downward.
	dydt_grid[n][i]	dummy
*/
double y_over[],**y_grid,**dydt_grid,y0_top[],y0_botm[],t_grid[],eps,omg;
int layerup_botm,layerup_top,layerdown_botm,layerdown_top,n_order,l,is_funcs;
void (*fdydt)();
{
	int i,n,n_grid,n_cross;

	if( layerup_top != layerdown_botm )
	{
		FR,"runge_driver_match(%s) %d input error:\n",__FILE__,__LINE__);
		FR,"upward integration stops at layer=%d\n",layerup_top);
		FR,"downward integration stops at layer=%d\n",layerdown_botm);
		exit(-1);
	}
	/* integrate downward */
	n_grid = layerdown_top-layerdown_botm+1;
	runge_driver_downto(y_grid,dydt_grid,y0_top,&(t_grid[layerdown_botm]),
	n_grid,fdydt,n_order,eps,l,omg,NO,is_funcs,&n_cross);

	for(i=0;i<n_order;i++) y_over[i] = y_grid[i][0];
	if( is_funcs == YES )
	for(i=0;i<n_order;i++)
	for(n=n_grid;n-->0;)
		y_grid[i][n+layerdown_botm] = y_grid[i][n];

	/* integrate upward */
	n_grid = layerup_top-layerup_botm+1;
	runge_driver(y_grid,dydt_grid,y0_botm,&(t_grid[layerup_botm]),
	n_grid,fdydt,n_order,eps,l,omg,NO,is_funcs,&n_cross);

	if( is_funcs == YES )
	for(i=0;i<n_order;i++)
	for(n=n_grid;n-->0;)
		y_grid[i][n+layerup_botm] = y_grid[i][n];
}

#define TINY		1e-30
#define BIG_SCALE		1e+10
#define SMALL_SCALE	1e-10
#define MAX_STEPS		20000
void runge_driver(y_grid,dydt_grid,y0,t_grid,n_grid,fdydt,n_order,eps,l,omg,is_scale,is_funcs,n_cross)
/* input
	y0[],		initial value of y at t_grid[0]	not modified upon exit
	t_grid[],	grid values of t axis where y and dydt are evaluated
								not modified upon exit
	n_grid,	number of grid points 0<=i < n_grid
	fdydt,	function to give dy/dt(t,y)
	n_order,	dimension of linear differentail equation
	eps,		required accuracy
	l,		angular order
	omg,		angular frequency
	is_scale,	do scaling == YES, if not== NO
	is_funcs,	store all eigenfunctions == YES if not == NO
  output
	y_grid[n][i],	estimated value of y at t_grid[i], 0<=n< n_order
	dydt_grid[n][i],	estimated value of dy/dt at t_grid[i], 0<=n< n_order
	n_corss,		number of zero crossing of y[0]
*/
double **y_grid,**dydt_grid,y0[],t_grid[],eps,omg;
int n_grid,n_order,l,is_scale,is_funcs,*n_cross;
void (*fdydt)();
{
	double y[MAX_FUNCTION_NUMBER],dydt[MAX_FUNCTION_NUMBER];
	double yout[MAX_FUNCTION_NUMBER],yscale[MAX_FUNCTION_NUMBER];
	double t,dt,dtdid,dtnext;
	int i,j,t_index,overrun,counter,scale,n_grid_1;

	t_index=0;
	*n_cross=0;
	n_grid_1 = n_grid - 1;
	if( n_grid_1 == 0 )
	{
		for(i=0;i<n_order;i++) y_grid[i][0] = y0[i];
		return;
	}else if ( n_grid < 0 )
	{
		FR,"runge_driver(%s) input error n_grid=%d.\n",
		__FILE__,n_grid);
		exit(-1);
	}
	for(i=0;i<n_order;i++)	y[i] = y0[i];
	t=t_grid[0];
	dt = 0.5*(t_grid[1]-t_grid[0]);
	overrun = 1;
	counter=MAX_STEPS;

	while(t_index<n_grid)
	{
		(*fdydt)(dydt,y,t,l,omg);
		/*
		 * store the eigenfunc
		 */
		if( overrun )
		{
			if( is_funcs == YES || t_index == n_grid_1 )
			for(i=0;i<n_order;i++)
			{
				y_grid[i][t_index] = y[i];
				dydt_grid[i][t_index] = dydt[i];
			}
			overrun=0;
			if( t_index == n_grid_1 ) break; /* reached to the end of grid*/
		}
		/* if may over run try with shorter dt */
		if( t_grid[t_index+1] < t + dt )
		{
			dt = t_grid[t_index+1]-t;
			overrun = 1;
		}

		for(i=0;i<n_order;i++)
			yscale[i] = fabs(y[i]) + fabs(dt*dydt[i]) + TINY;
		stepper_ck(y,yout,dydt,t,dt,fdydt,n_order,yscale,eps,&dtdid,&dtnext,l,omg);
		/* confirm dtdid is dt at a grid point */
		if( overrun==1 && dt == dtdid ) t_index ++;
		else overrun = 0;
		counter--;
/*
 * test scale of solutions
 */
		if( is_scale == YES )
		{
			scale = 0;
			for(i=0;i<n_order;i++)
			if( yout[i] > BIG_SCALE ) scale = 1;
			else if( yout[i] < SMALL_SCALE) scale =-1;
		}

		if( is_scale == YES && scale == 1  )
		{
			FR,"runge_driver: eigenfunc are scaled to shrink.\n");
			for(i=0;i<n_order;i++)
			{
				y[i] *= SMALL_SCALE;
				yout[i] *= SMALL_SCALE;
				for(j=0;j<=t_index;j++)
					y_grid[i][j]*=SMALL_SCALE;
			}
		}
/*
		else if(is_scale == YES &&  scale == -1 )
		{
			FR,"runge_driver: eigenfunc are scaled to enlarge.\n");
			for(i=0;i<n_order;i++)
			{
				y[i] *= BIG_SCALE;
				yout[i] *= BIG_SCALE;
				for(j=0;j<=t_index;j++)
					y_grid[i][j]*=BIG_SCALE;
			}
		}
*/
		if( y[0]*yout[0] < 0 ) (*n_cross)++;
		for(i=0;i<n_order;i++)
			y[i] = yout[i];
		dt = dtnext;
		t += dtdid;
		if( !counter) 
		{
			FR,"runge_driver(%s) %d: too many integration steps.\n",
			__FILE__,__LINE__); exit(-1);
		}
	}
}

void runge_driver_downto(y_grid,dydt_grid,y0,t_grid,n_grid,fdydt,n_order,eps,l,omg,is_scale,is_funcs,n_cross)
/* input
	y0[],		initial value of y at t_grid[n_grid-1] not modified upon exit
	t_grid[],	grid values of t axis where y and dydt are evaluated
								not modified upon exit
	n_grid,	number of grid points 0<=i < n_grid
	fdydt,	function to give dy/dt(t,y)
	n_order,	dimension of linear differentail equation
	eps,		required accuracy
	l,		angular order
	omg,		angular frequency
	is_scale,	do scaling == YES, if not== NO
	is_funcs,	store all eigenfunctions == YES if not == NO
  output
	y_grid[n][i],	estimated value of y at t_grid[i], 0<=n< n_order
	dydt_grid[n][i],	estimated value of dy/dt at t_grid[i], 0<=n< n_order
	n_corss,		number of zero crossing of y[0]
*/
double **y_grid,**dydt_grid,y0[],t_grid[],eps,omg;
int n_grid,n_order,l,is_scale,is_funcs,*n_cross;
void (*fdydt)();
{
	double y[MAX_FUNCTION_NUMBER],dydt[MAX_FUNCTION_NUMBER];
	double yout[MAX_FUNCTION_NUMBER],yscale[MAX_FUNCTION_NUMBER];
	double t,dt,dtdid,dtnext;
	int i,t_index,overrun,counter,n_grid_1;

	*n_cross=0;
	n_grid_1 = n_grid - 1;
	if( n_grid_1 == 0 )
	{
		for(i=0;i<n_order;i++) y_grid[i][0] = y0[i];
		return;
	}else if ( n_grid < 0 )
	{
		FR,"runge_driver_down(%s) input error n_grid=%d.\n",
		__FILE__,n_grid);
		exit(-1);
	}
	t_index=n_grid_1;
	for(i=0;i<n_order;i++)	y[i] = y0[i];
	t=t_grid[n_grid_1];
	dt = 0.5*(t_grid[n_grid_1 -1]-t_grid[n_grid_1]);
	/* if t_grid[] is increasing order dt < 0.0 */
	overrun = 1;
	counter=MAX_STEPS;

	while(t_index>-1)
	{
		(*fdydt)(dydt,y,t,l,omg);
		/*
		 * store the eigenfunc
		 */
		if( overrun )
		{
			if( is_funcs == YES || t_index == 0 )
			for(i=0;i<n_order;i++)
			{
				y_grid[i][t_index] = y[i];
				dydt_grid[i][t_index] = dydt[i];
			}
			overrun=0;
			if( t_index == 0 ) break; /* reached to the bottom of grid*/
		}
		/* if may over run try with shorter dt */
		if( t_grid[t_index-1] > t + dt )
		{
			dt = t_grid[t_index-1]-t;
			overrun = 1;
		}

		for(i=0;i<n_order;i++)
			yscale[i] = fabs(y[i]) + fabs(dt*dydt[i]) + TINY;
		stepper_ck(y,yout,dydt,t,dt,fdydt,n_order,yscale,eps,&dtdid,&dtnext,l,omg);
		/* confirm dtdid is dt at a grid point */
		if( overrun==1 && dt == dtdid ) t_index --;
		else overrun = 0;
		counter--;

		if( y[0]*yout[0] < 0 ) (*n_cross)++;
		for(i=0;i<n_order;i++)
			y[i] = yout[i];
		dt = dtnext;
		t += dtdid;
		if( !counter) 
		{
			FR,"runge_driver_downto(%s) %d: too many integration steps.\n",
			__FILE__,__LINE__); exit(-1);
		}
	}
}
void runge_driver_sol(y_grid,dydt_grid,y0,t_grid,n_grid,fdydt,n_order,eps,l,omg,is_scale,is_funcs,n_cross)
/* a short version of runge_driver_sol().
	## NOTE ## Aug. 6 1994
	This function does not
	compute the eigenfunction at the grid points of an earth model.
	This routine would save CPU time, however, runs slower than
	runge_driver_sol() for a real earth model in which model parameters
	do not change smoothly. For example
	model a_model.1
	l=1000, omg1=0.01 omg2=0.02 domg=0.0004
	with runge_drivder() comp_rk.c takes 13 sec CPU time.
	with runge_drivder_sol() comp_rk.c takes 25 sec CPU time.
	## END NOTE ##
  input
	y0[],		initial value of y at t_grid[0]	not modified upon exit
	t_grid[],	grid values of t axis where y and dydt are evaluated
								not modified upon exit
	n_grid,	number of grid points 0<=i < n_grid
	fdydt,	function to give dy/dt(t,y)
	n_order,	dimension of linear differentail equation
	eps,		required accuracy
	l,		angular order
	omg,		angular frequency
	is_scale,	do scaling == YES, if not== NO
	is_funcs,	store all eigenfunctions == YES if not == NO
  output
	y_grid[n][i],	estimated value of y at t_grid[i], 0<=n< n_order
	dydt_grid[n][i],	estimated value of dy/dt at t_grid[i], 0<=n< n_order
	n_corss,		number of zero crossing of y[0]
*/
double **y_grid,**dydt_grid,y0[],t_grid[],eps,omg;
int n_grid,n_order,l,is_scale,is_funcs,*n_cross;
void (*fdydt)();
{
	double y[MAX_FUNCTION_NUMBER],dydt[MAX_FUNCTION_NUMBER];
	double yout[MAX_FUNCTION_NUMBER],yscale[MAX_FUNCTION_NUMBER];
	double t,dt,dtdid,dtnext;
	int i,overrun,counter,scale;

	*n_cross=0;
	for(i=0;i<n_order;i++)	y[i] = y0[i];
	t=t_grid[0];
	dt = 0.5*(t_grid[1]-t_grid[0]);
	overrun = 0;
	counter=MAX_STEPS;

	while(1)
	{
		(*fdydt)(dydt,y,t,l,omg);
		/*
		 * store the eigenfunc
		 */
		if( overrun )
		{
			for(i=0;i<n_order;i++)
			{
				y_grid[i][n_grid-1] = y[i];
				dydt_grid[i][n_grid-1] = dydt[i];
			}
			break;
		}

		if( t_grid[n_grid-1] < t + dt )
		{
			dt = t_grid[n_grid-1]-t;
			overrun = 1;
		}
		for(i=0;i<n_order;i++)
			yscale[i] = fabs(y[i]) + fabs(dt*dydt[i]) + TINY;
		stepper_ck(y,yout,dydt,t,dt,fdydt,n_order,yscale,eps,&dtdid,&dtnext,l,omg);
		if( overrun==1 && dt != dtdid ) overrun=0;
		counter--;
/*
 * test scale of solutions
 */
		if( is_scale == YES )
		{
			scale = 0;
			for(i=0;i<n_order;i++)
			if( yout[i] > BIG_SCALE ) scale = 1;
			else if( yout[i] < SMALL_SCALE) scale =-1;
		}

		if( is_scale == YES && scale == 1  )
		{
			FR,"runge_driver: eigenfunc are scaled to shrink.\n");
			for(i=0;i<n_order;i++)
			{
				y[i] *= SMALL_SCALE;
				yout[i] *= SMALL_SCALE;
			}
		}
/*
		else if(is_scale == YES &&  scale == -1 )
		{
			FR,"runge_driver: eigenfunc are scaled to enlarge.\n");
			for(i=0;i<n_order;i++)
			{
				y[i] *= BIG_SCALE;
				yout[i] *= BIG_SCALE;
			}
		}
*/
		if( y[0]*yout[0] < 0 ) (*n_cross)++;
		for(i=0;i<n_order;i++)
			y[i] = yout[i];
		dt = dtnext;
		t += dtdid;
		if( !counter) 
		{
			FR,"runge_driver(%s) %d: too many integration steps.\n",
			__FILE__,__LINE__); exit(-1);
		}
	}
}
#undef TINY
#undef MAXSTEPS

#define SMALL 	1e-8
#define MAXSTEPS 	100
int zbrent(best_omg,driver,boundary_condition,bc,omg,omg1,bc_value,bc_value1,y_grid,dydt_grid,tol,y0,t_grid,n_grid,fdydt,n_order,eps,l)
/* input
	driver,	runge kutta driver
	boundary,	function of boundary condition
	bc,		boundary condition
	y0[],		initial value of y at t_grid[0]
	**y_grid,	dummy array used in (*driver)();
	**dydt_grid,dummy array used in (*driver)();
	t_grid[],	grid values of t axis where y and dydt are evaluated
	n_grid,	number of grid points 0<=i < n_grid
	fdydt,	function to give dy/dt(t,y)
	n_order,	dimension of linear differentail equation
	eps,		required accuracy for runge-kuttta integration
	tol,		required accuracy for zbrent() function
	l,		angular order
	omg,		angular frequency lower bound
	omg1,		angular frequency upper bound
	bc_value,	boundary() at omg
	bc_value1,	boundary() at omg1
  output
	best_omg,	estimated eigenfrequency
-------------------------------------------------------------------
	find root which satisfies (*boundary)()=0.0;

	this routeine based on Numerical Recipes p.251
	chapter 9.3 Van Wijngaarden-dekker-brent Method
*/
	double *best_omg,omg,omg1,tol,**y_grid,**dydt_grid,y0[],t_grid[],eps;
	double (*boundary_condition)(),bc_value,bc_value1;
	void (*driver)();
	void (*fdydt)();
	int l,n_grid,n_order,bc;
{
	double a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm,min1,min2;
	double y_bound[MAX_FUNCTION_NUMBER];
	double dydt_bound[MAX_FUNCTION_NUMBER];
	int i,n_cross,bound_grid;

	a = omg;
	b = omg1;

	fa = bc_value;
	fb = bc_value1;
/*
	(*driver)(y_grid,dydt_grid,y0,t_grid,n_grid,fdydt,n_order,eps,l,a,NO,NO,&n_cross);
	for(i=0;i<n_order;i++) y_bound[i] = y_grid[i][n_grid-1];
	for(i=0;i<n_order;i++) dydt_bound[i] = dydt_grid[i][n_grid-1];
	fa = (*boundary_condition)(y_bound,dydt_bound,n_order,l,omg,bc);

	(*driver)(y_grid,dydt_grid,y0,t_grid,n_grid,fdydt,n_order,eps,l,b,NO,NO,&n_cross);
	for(i=0;i<n_order;i++) y_bound[i] = y_grid[i][n_grid-1];
	for(i=0;i<n_order;i++) dydt_bound[i] = dydt_grid[i][n_grid-1];
	fb = (*boundary_condition)(y_bound,dydt_bound,n_order,l,omg1,bc);
*/
	if( driver == runge_driver || driver == runge_driver_sol ){
		bound_grid = n_grid-1;
	}else if (driver == runge_driver_downto ) {
		bound_grid = 0;
	}else{
		FR,"zbrent(%s) %d: unknown driver function.\n",__FILE__,__LINE__);
		exit(-1);
	}

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
		(*driver)(y_grid,dydt_grid,y0,t_grid,n_grid,fdydt,n_order,eps,l,b,NO,NO,&n_cross);
		for(i=0;i<n_order;i++) y_bound[i] = y_grid[i][bound_grid];
		for(i=0;i<n_order;i++) dydt_bound[i] = dydt_grid[i][bound_grid];
		fb = (*boundary_condition)(y_bound,dydt_bound,n_order,l,b,bc);
	}
}
int zbrent_match(best_omg,driver,
boundary_condition,bc,omg,omg1,
bc_value,bc_value1,y_grid,dydt_grid,tol,y0_botm,y0_top,
t_grid,n_grid,layer_match,fdydt,n_order,eps,l)
/* input
	driver,	runge kutta driver
	boundary,	function of boundary condition
	bc,		boundary condition
	y0_botm[],		initial value of y at t_grid[0]
	y0_top[],		initial value of y at t_grid[n_grid-1]
	**y_grid,	dummy array used in (*driver)();
	**dydt_grid,dummy array used in (*driver)();
	t_grid[],	grid values of t axis where y and dydt are evaluated
	n_grid,	number of grid points 0<=i < n_grid
	layer_match,matching condtion is evaluated at this layer
	fdydt,	function to give dy/dt(t,y)
	n_order,	dimension of linear differentail equation
	eps,		required accuracy for runge-kuttta integration
	tol,		required accuracy for zbrent() function
	l,		angular order
	omg,		angular frequency lower bound
	omg1,		angular frequency upper bound
	bc_value,	boundary() at omg
	bc_value1,	boundary() at omg1
  output
	best_omg,	estimated eigenfrequency
-------------------------------------------------------------------
	find root which satisfies (*boundary)()=0.0;

	this routeine based on Numerical Recipes p.251
	chapter 9.3 Van Wijngaarden-dekker-brent Method
*/
	double *best_omg,omg,omg1,tol,**y_grid,**dydt_grid,
		y0_botm[],y0_top[],t_grid[],eps;
	double (*boundary_condition)(),bc_value,bc_value1;
	void (*driver)();
	void (*fdydt)();
	int l,n_grid,layer_match,n_order,bc;
{
	double a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm,min1,min2;
	double y_bound[MAX_FUNCTION_NUMBER];
	double dydt_bound[MAX_FUNCTION_NUMBER];
	int i;
	double y_over[MAX_FUNCTION_NUMBER];

	a = omg;
	b = omg1;

	fa = bc_value;
	fb = bc_value1;
/*
	(*driver)(y_grid,dydt_grid,y0,t_grid,n_grid,fdydt,n_order,eps,l,a,NO,NO,&n_cross);
	for(i=0;i<n_order;i++) y_bound[i] = y_grid[i][n_grid-1];
	for(i=0;i<n_order;i++) dydt_bound[i] = dydt_grid[i][n_grid-1];
	fa = (*boundary_condition)(y_bound,dydt_bound,n_order,l,omg,bc);

	(*driver)(y_grid,dydt_grid,y0,t_grid,n_grid,fdydt,n_order,eps,l,b,NO,NO,&n_cross);
	for(i=0;i<n_order;i++) y_bound[i] = y_grid[i][n_grid-1];
	for(i=0;i<n_order;i++) dydt_bound[i] = dydt_grid[i][n_grid-1];
	fb = (*boundary_condition)(y_bound,dydt_bound,n_order,l,omg1,bc);
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
		(*driver)(y_over,y_grid,dydt_grid,y0_botm,y0_top,t_grid,
		0,layer_match,layer_match,n_grid-1,
		fdydt,n_order,eps,l,b,NO);

		for(i=0;i<n_order;i++) y_bound[i] = y_grid[i][layer_match];
		for(i=0;i<n_order;i++) dydt_bound[i] = dydt_grid[i][layer_match];
		fb = (*boundary_condition)(y_bound,y_over,n_order,l,b,bc);
	}
}
