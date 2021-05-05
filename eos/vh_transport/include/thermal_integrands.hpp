
#include <stdlib.h>

#ifndef THERMAL_INTEGRANDS_H

#define THERMAL_INTEGRANDS_H


double I00_integrand(double pbar, double mbar);
double I01_integrand(double pbar, double mbar);

double I11_integrand(double pbar, double mbar);


double I21_integrand(double pbar, double mbar);
double I22_integrand(double pbar, double mbar);

double I32_integrand(double pbar, double mbar);

double I40_integrand(double pbar, double mbar);
double I41_integrand(double pbar, double mbar);
double I42_integrand(double pbar, double mbar);

#endif