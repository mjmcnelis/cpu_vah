#include <stdlib.h>
#include <math.h>
#include "thermal_integrands.hpp"


double I00_integrand(double pbar, double mbar)
{
	// gauss laguerre (a = 0) boltzmann
	double pbar2 = pbar * pbar;
	double Ebar = sqrt(pbar2  +  mbar * mbar);

	return pbar2 / Ebar * exp(pbar - Ebar);
}

double I01_integrand(double pbar, double mbar)
{
	// gauss laguerre (a = 0) boltzmann
	double pbar2 = pbar * pbar;
	double Ebar = sqrt(pbar2  +  mbar * mbar);

	return pbar2 * pbar2 / (Ebar * Ebar * Ebar) * exp(pbar - Ebar);
}

double I11_integrand(double pbar, double mbar)
{
	// gauss laguerre (a = 1) boltzmann
	double pbar2 = pbar * pbar;
	double Ebar = sqrt(pbar2  +  mbar * mbar);

	return pbar * pbar2 / (Ebar * Ebar) * exp(pbar - Ebar);
}

double I21_integrand(double pbar, double mbar)
{
	// gauss laguerre (a = 2) boltzmann
	double pbar2 = pbar * pbar;
	double Ebar = sqrt(pbar2  +  mbar * mbar);

	return pbar2 / Ebar * exp(pbar - Ebar);
}

double I22_integrand(double pbar, double mbar)
{
	// gauss laguerre (a = 2) boltzmann
	double pbar2 = pbar * pbar;
	double Ebar = sqrt(pbar2  +  mbar * mbar);

	return pbar2 * pbar2 / (Ebar * Ebar * Ebar) * exp(pbar - Ebar);
}

double I32_integrand(double pbar, double mbar)
{
	// gauss laguerre (a = 3) boltzman
	double pbar2 = pbar * pbar;
	double Ebar = sqrt(pbar2  +  mbar * mbar);

	return pbar * pbar2 / (Ebar * Ebar) * exp(pbar - Ebar);
}

double I40_integrand(double pbar, double mbar)
{
	// gauss laguerre (a = 4) boltzman
	double pbar2 = pbar * pbar;
	double Ebar = sqrt(pbar2  +  mbar * mbar);

	return Ebar * Ebar * Ebar / pbar2 * exp(pbar - Ebar);
}

double I41_integrand(double pbar, double mbar)
{
	// gauss laguerre (a = 4) boltzman
	double pbar2 = pbar * pbar;
	double Ebar = sqrt(pbar2  +  mbar * mbar);

	return Ebar * exp(pbar - Ebar);
}

double I42_integrand(double pbar, double mbar)
{
	// gauss laguerre (a = 4) boltzman
	double pbar2 = pbar * pbar;
	double Ebar = sqrt(pbar2  +  mbar * mbar);

	return pbar2 / Ebar * exp(pbar - Ebar);
}










