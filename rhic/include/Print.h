#ifndef PRINT_H_
#define PRINT_H_

void print_hydro_center(double t, double e, int s);

void print_parameters(int nx, int ny, int nz, double dt, double dx, double dy, double dz, double t0, double T_switch, double etabar, int adaptive_time_step);


#endif