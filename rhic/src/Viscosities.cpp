#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include "../include/Hydrodynamics.h"
#include "../include/Precision.h"
#include "../include/Parameters.h"

inline precision sign(precision x)
{
	if(x > 0.) return 1.;
	else if(x < 0.) return -1.;
	else return 0.;
}

inline precision Theta(precision x)
{
	if(x > 0.)
	{
		return 1.;
	}

	return 0.;
}

precision eta_over_s(precision T, hydro_parameters hydro)
{
	if(hydro.temperature_etas == 1)								// temperature dependent piecewise parameterization
	{
		// precision etas_min = hydro.etas_min;				// etas value for T < Tc
		// precision etas_slope = hydro.etas_slope * hbarc;	// etas slope (converted to fm)
		// precision Tc = 0.154 / hbarc;						// psuedocritical temperature (fm^-1)

		// return fmax(etas_min, etas_min  +  etas_slope * (T - Tc));

		// precision aL = -0.77 * hbarc;
		// precision aH = 0.21 * hbarc;
		// precision Tk = 0.22 / hbarc;
		// precision etask = 0.093;

		precision aL = hydro.etas_aL * hbarc;			// left slope [fm]
		precision aH = hydro.etas_aH * hbarc;			// right slope [fm]
		precision Tk = hydro.etas_Tk_GeV / hbarc;		// kink temperature fm^-1]
		precision etask = hydro.etas_etask;				// kink etas value

		return etask  +  (T - Tk) * (aL * Theta(Tk - T)  +  aH * Theta(T - Tk));
	}

	return hydro.constant_etas;
}


precision zeta_over_s(precision T, hydro_parameters hydro)
{
	//precision Tpeak = hydro.zetas_peak_temperature_GeV / hbarc;
	// precision x = T / Tpeak;

	// precision zetas;

	// if(x > 1.05)
	// {
	// 	zetas = 0.9 * exp(-(x - 1.) / 0.025)  +  0.25 * exp(-(x - 1.) / 0.13)  +  0.001;
	// }
	// else if(x < 0.995)
	// {
	// 	zetas = 0.9 * exp((x - 1.) / 0.0025)  +  0.22 * exp((x - 1.) / 0.022)  +  0.03;
	// }
	// else
	// {
	// 	zetas = - 13.77 * x * x  +  27.55 * x  -  13.45;
	// }
	// return zetas * hydro.zetas_normalization_factor;


	// precision norm = hydro.zetas_normalization_factor;				// normalization factor
	// precision Tpeak = hydro.zetas_peak_temperature_GeV / hbarc;		// peak temperature [fm^-1]
	// precision width = 0.089 / hbarc;								// width [fm^-1]
	// precision skew = -0.15;											// skew

	precision norm = hydro.zetas_normalization_factor;				// normalization factor
	precision Tpeak = hydro.zetas_peak_temperature_GeV / hbarc;		// peak temperature [fm^-1]
	precision width = hydro.zetas_width_GeV / hbarc;				// width [fm^-1]
	precision skew = hydro.zetas_skew;								// skew

	precision Lambda = width * (1.  +  skew * sign(T - Tpeak));

	return norm * Lambda * Lambda / (Lambda * Lambda  +  (T - Tpeak) * (T - Tpeak));
}




