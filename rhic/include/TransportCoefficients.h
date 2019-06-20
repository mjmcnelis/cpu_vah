#ifndef TRANSPORTCOEFFICIENTS_H_
#define TRANSPORTCOEFFICIENTS_H_

#include "DynamicalVariables.h"

double transversePressureHat(double e, double p, double pl);
double Rbar0_fun(double a);
double Rbar0P_fun(double a);

void secondOrderTransportCoefficientsZ(double e, double p, double pl, double cs2, double T,
double *beta_lPi, double *delta_lPi, double *lambda_piPi, double *beta_PiPi, double *delta_PiPi, double *lambda_Pipi);


// perhaps I should make a struct containing all the transport coefficients
// it would be faster to calculate them

#endif
