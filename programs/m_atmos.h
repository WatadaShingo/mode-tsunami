#ifndef m_atmos_h
#define m_atmos_h

#define	SAMPLE_ATMOS	57
#define	PARAMETERS_ATMOS	8

#define	RAD_ATMOS		0	/* radius from the center of the earth */
#define	RHO_ATMOS		1	/* density */
#define	T_ATMOS		2	/* temparature in Kelvin*/
#define	P_ATMOS		3	/* pressure in Pa=N/m^2 */
#define	VP_ATMOS		4	/* sound velocity in km/sec */
#define	KAPPA_ATMOS		5	/* bulk modulous Pa=N/m^2 */
#define	G_ATMOS		6	/* gracity km/sec^2 */
#define	DRHODR_ATMOS	7	/* density gradient */

struct model_atmos{
	int nmodel; /* number of radial samples in the model */
/*	float value[PARAMETERS_ATMOS][SAMPLE_ATMOS]; */
	double value[PARAMETERS_ATMOS][SAMPLE_ATMOS];
	char modelname[20];
	char modelfilename[100];
	};
struct eigen_atmos{
	int l,n;
/*	float omg,vg,u[SAMPLE_ATMOS],du[SAMPLE_ATMOS],*/
	double omg,vg,u[SAMPLE_ATMOS],du[SAMPLE_ATMOS],
	p[SAMPLE_ATMOS],dp[SAMPLE_ATMOS];
	};
struct model_funcatmos{
/*	float r,rho,t,p,vp,kappa,g,drhodr;*/
	double r,rho,t,p,vp,kappa,g,drhodr;
	};

#endif
