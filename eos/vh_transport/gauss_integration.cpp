
#include <stdlib.h>
#include "gauss_integration.hpp"




// for thermal equilibrium functions (1D integral over pbar coordinate from 0 to infinity)
double Gauss1D(double thermal_1D_integrand(double pbar, double mbar), double * pbar_root, double * pbar_weight, int pbar_pts, double mbar)
{
	double sum = 0.0;
	for(int k = 0; k < pbar_pts; k++) sum += pbar_weight[k] * thermal_1D_integrand(pbar_root[k], mbar);
	return sum;
}

