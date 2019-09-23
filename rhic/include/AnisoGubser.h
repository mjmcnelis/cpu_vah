
#ifndef ANISOGUBSER_H_
#define ANISOGUBSER_H_

#include "Parameters.h"

const double eps = 1.e-5;			// prevent to cubic spline interpolation errors (if rho outside range)
const double drho_default = 1.e-4;	// warning: if too small, get a segfault (e_hat, pl_hat array size is too large)
const double T0_hat_error = 1.e-8;	// tolerance error for T0_hat search

double rho_function(double t, double r, double q);
void gubser_rho_evolution(double * e_hat, double * pl_hat, double * rho_array, int rho_pts, double drho, double t, hydro_parameters hydro);
double run_semi_analytic_aniso_gubser(lattice_parameters lattice, initial_condition_parameters initial, hydro_parameters hydro);
void set_aniso_gubser_energy_density_and_flow_profile(double T0_hat, int nx, int ny, int nz, double dt, double dx, double dy, double dz, hydro_parameters hydro, initial_condition_parameters initial);

#endif