#ifndef rk_h
#define rk_h
int rungev4();
int rungevck();
void stepper();
void stepper_ck();
void runge_driver();
void runge_driver_match();
void runge_driver_downto();
void runge_driver_sol();
int zbrent();
int zbrent_match();
#endif
