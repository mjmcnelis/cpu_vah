#include <stdlib.h>
#include <math.h>
#include "thermal_integrands.hpp"

double I32_integrand(double pbar, double mbar, int sign)
{
	// gauss laguerre (a = 3) boltzman
	double pbar2 = pbar * pbar;
	double Ebar = sqrt(pbar2 + mbar*mbar);

	return pbar*pbar2 / (Ebar*Ebar) * exp(pbar-Ebar);
}

double I11_integrand(double pbar, double mbar, int sign)
{
	// gauss laguerre (a = 1) boltzmann
	double pbar2 = pbar * pbar;
	double Ebar = sqrt(pbar2 + mbar*mbar);

	return pbar*pbar2 / (Ebar*Ebar) * exp(pbar-Ebar);
}