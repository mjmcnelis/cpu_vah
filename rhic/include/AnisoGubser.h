
#ifndef ANISOGUBSER_H_
#define ANISOGUBSER_H_

#include "Parameters.h"

double rho_function(double t, double r, double q);
void gubser_rho_evolution(double * e_hat, double * pl_hat, double * rho_array, int rho_pts, double drho, hydro_parameters hydro);
void run_semi_analytic_aniso_gubser(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro);

#endif