
#include <stdlib.h>

#ifndef GAUSS_INTEGRATION_H

#define GAUSS_INTEGRATION_H


double Gauss1D(double thermal_1D_integrand(double pbar, double mbar), double * pbar_root, double * pbar_weight, int pbar_pts, double mbar);


#endif
