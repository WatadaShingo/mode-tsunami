#ifndef eigen_h
#define eigen_h

#define MAX_EIGEN_GRID_POINT	300
#define UPPER_REGION_EIGEN	0
#define LOWER_REGION_EIGEN	1

struct eigen_func{		/* eigen function at depth r */
	int i;	/* knot number */
	int n,l;	/* radial and angular order */
	float u,du,v,dv,p,dp;
	double r;
	};

struct eigen_mode{	/* eigen function from r=0 to r=ra */
	int l,n;
	float omg,dum[9];
	float u[MAX_EIGEN_GRID_POINT],du[MAX_EIGEN_GRID_POINT],
		v[MAX_EIGEN_GRID_POINT],dv[MAX_EIGEN_GRID_POINT],
		p[MAX_EIGEN_GRID_POINT],dp[MAX_EIGEN_GRID_POINT];
	};
struct eigen_modecowl{	/* eigen function from r=0 to r=ra */
	int l,n;
	float omg,dum[9];
	float u[MAX_EIGEN_GRID_POINT],du[MAX_EIGEN_GRID_POINT],
		v[MAX_EIGEN_GRID_POINT],dv[MAX_EIGEN_GRID_POINT];
	};
	
void normalize_eigen();
void re_normalize_eigen();
void print_eigen();
void write_eigen();
void write_eigencowl();
int read_eigen();
int read_eigencowl();
void get_layer();
void get_eigen_at_r();
void get_eigen_at_rcowl();

#endif
