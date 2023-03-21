#ifndef m_solid_h
#define m_solid_h

/*#define 	SAMPLE_SOLID	222	*//* for PREM */
#define 	SAMPLE_SOLID	279   /* for PREM + atmosphere */
#define 	PARAMETERS_SOLID  18
#define	MAX_N_REGIONS	14

#define	REGOIN_SOLID	0
#define	LAYER_SOLID		1
#define 	RAD_SOLID		2 /* radius from the center of the earth */
#define 	RHO_SOLID		3 /* density */
#define	VPV_SOLID		4 
#define	VPH_SOLID		5 
#define	VSV_SOLID		6 
#define	VSH_SOLID		7 
#define	ETA_SOLID		8
#define	QU_SOLID		9
#define	QK_SOLID		10
#define	A_SOLID		11
#define	C_SOLID		12
#define	L_SOLID		13
#define	N_SOLID		14
#define	F_SOLID		15
#define	G_SOLID		16
#define	DRHODR_SOLID	17

#define	UPPER_BOUND_SOLID		0
#define	LOWER_BOUND_SOLID		1
#define	NUMBER_OF_GRIDS		2

#define	UPPER_REGION		0
#define	LOWER_REGION		1

struct model_solid{
	int nmodel; /* number of radial samples in the model */
	int nregion;
	int region[MAX_N_REGIONS][3];
/*	float value[PARAMETERS_SOLID][SAMPLE_SOLID];*/
	double value[PARAMETERS_SOLID][SAMPLE_SOLID];
	double value0[PARAMETERS_SOLID][SAMPLE_SOLID];
	/* reference model at T=1 sec */
	char modelname[20];
	char modelfilename[100];
	};
struct model_funcsolid{
/*	float r,rho,vpv,vph,vsv,vsh,eta,qu,qk,a,c,l,n,f,g;*/
	double r,rho,vpv,vph,vsv,vsh,eta,qu,qk,a,c,l,n,f,g,drhodr;
	};

#endif
