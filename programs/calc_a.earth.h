#ifndef calc_a_earth_h
#define calc_a_earth_h

void mult_m_v();
void mult_m_v2();
void mult_m_v4();
void mult_m_v6();
void dzdr_atmos2by2();
void dzdr_solid6by6();
void dzdr_solid4by4cowl();
void dzdr_solid4by4cowlcomplex();
void dzdr_solid2by2radial();
void dzdr_solid2by2radialcomplex();
void dzdr_fluid4by4();
void dzdr_fluid2by2cowl();
void dzdr_fluid2by2cowlcomplex();
void dzdr_fluid2by2radial();
void dzdr_fluid2by2radialcomplex();
void calc_atmos2by2();
void calc_solid6by6();
void calc_solid4by4cowl();
void calc_solid4by4cowlcomplex();
void calc_solid2by2radial();
void calc_solid2by2radialcomplex();
void calc_fluid4by4();
void calc_fluid2by2cowl();
void calc_fluid2by2cowlcomplex();
void calc_fluid2by2_character();
void calc_fluid2by2radial();
void calc_fluid2by2radialcomplex();
void calc_fluid_omg();
void y_initial_atmos2by2();
void y_initial_solid6by6();
void y_initial_solid4by4cowl();
void y_initial_fluid4by4();
void y_initial_fluid2by2cowl();
void get_eigenfunc_atmos2by2_from_x();
void get_x3_atmos2by2_from_x();

void continue_eigenfunc();
void continue_eigenfunccowl();
double determinant3by3complex();
double determinant2by2complex();
double vector_cos();
void boundary_coeff_fluid2by2();

#define FIXED_BC	1
#define FREE_BC	2

#define SOLID_TO_FLUID	1
#define FLUID_TO_SOLID	2
#define SOLID_TO_FLUID_COMPLEX 3
#define FLUID_TO_SOLID_COMPLEX 4
#endif
