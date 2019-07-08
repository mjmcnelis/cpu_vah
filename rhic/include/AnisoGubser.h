
#ifndef ANISOGUBSER_H_
#define ANISOGUBSER_H_


double rho_function(double t, double r, double q);
void gubser_rho_evolution(double * e_hat, double * pl_hat, double * rho_array, int rho_pts, double drho, double etas);

#endif