#ifndef VISCOUSGUBSER_H_
#define VISCOUSGUBSER_H_

#include "Parameters.h"

double run_semi_analytic_viscous_gubser(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro);

void set_viscous_gubser_initial_condition(double T0_hat, int nx, int ny, int nz, double dt, double dx, double dy, double dz, hydro_parameters hydro, initial_condition_parameters initial);

#endif