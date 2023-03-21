#ifndef search_h
#define search_h

#define NO_CROSSING	0
#define CROSSING		1

#define INSIDE_XY		1
#define INSIDE_MARGIN	2
#define OUTSIDE_MARGIN	3

#define EVANESCENT_DOMAIN     0x1 /* used in domain */
#define ACOUSTIC_DOMAIN       0x2
#define GRAVITY_DOMAIN        0x4
#define PROPAGATING_DOMAIN	(ACOUSTIC_DOMAIN | GRAVITY_DOMAIN)

#define FOUND_OMG_NOTHING_EVA	0
#define FOUND_OMG_NOTHING_GRA	1
#define FOUND_OMG_NOTHING_ACO	2
#define FOUND_OMG_NOTHING_PRO	3
#define FOUND_OMG_ACOUSTIC	4
#define FOUND_OMG_GRAVITY	5
#define FOUND_OMG_BOTH		6

void cross_2d();
void cross_1d();
void inside_outside();
void improve_dx();
void search_omg();
void find_omg_bound();
void determine_omg_domain();

#endif
