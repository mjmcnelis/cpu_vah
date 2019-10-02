#include <stdlib.h>
#include <math.h>
#include <cmath>
#include "../include/Hydrodynamics.h"
#include "../include/Precision.h"
#include "../include/Parameters.h"


precision eta_over_s(precision T, hydro_parameters hydro)
{
	if(hydro.temperature_etas)								// temperature dependent piecewise parameterization
	{
		precision etas_min = hydro.etas_min;				// etas value for T < Tc
		precision etas_slope = hydro.etas_slope * hbarc;	// etas slope (converted to fm)
		precision Tc = 0.154 / hbarc;						// psuedocritical temperature (fm^-1)

		return fmax(etas_min, etas_min  +  etas_slope * (T - Tc));
	}
	return hydro.constant_etas;
}


precision zeta_over_s(precision T, hydro_parameters hydro)
{
	precision Tpeak = hydro.zetas_peak_temperature_GeV / hbarc;
	precision x = T / Tpeak;

	precision zetas;

	if(x > 1.05)
	{
		zetas = 0.9 * exp(-(x - 1.) / 0.025)  +  0.25 * exp(-(x - 1.) / 0.13)  +  0.001;
	}
	else if(x < 0.995)
	{
		zetas = 0.9 * exp((x - 1.) / 0.0025)  +  0.22 * exp((x - 1.) / 0.022)  +  0.03;
	}
	else
	{
		zetas = - 13.77 * x * x  +  27.55 * x  -  13.45;
	}
	return zetas * hydro.zetas_normalization_factor;
}




